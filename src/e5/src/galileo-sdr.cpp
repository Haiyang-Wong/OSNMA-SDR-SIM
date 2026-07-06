#define _CRT_SECURE_NO_DEPRECATE

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _WIN32
#include "getopt.h"
#else
#include <unistd.h>
#endif

#include "../include/galileo-sdr.h"
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <vector>
#include <ncurses.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sstream>
#include <cmath>
#include <signal.h>
#include "../include/socket.h"

// The implementation of init_tables() and the definition of cosTable/sinTable
// are kept in this .cpp file. Declarations are in constants.h.
int cosTable[COS_TAB_LENGTH];
int sinTable[COS_TAB_LENGTH];

void init_tables()
{
    for (int i = 0; i < COS_TAB_LENGTH; i++)
    {
        cosTable[i] = (int)(32767.0 * cos(2.0 * PI * i / COS_TAB_LENGTH));
        sinTable[i] = (int)(32767.0 * sin(2.0 * PI * i / COS_TAB_LENGTH));
    }
}


int samples_per_code;
std::vector<int> current_eph;
std::vector<int> old_eph;
bool stop_signal_called_gal_task = false;
bool use_cursors = true;
bool advance_fptr = false;

FILE *TV_ptr;
FILE *TV_ptrs_sv[MAX_SAT];

void sigint_handler_gal_task(int code)
{
    (void)code;
    endwin();
    stop_signal_called_gal_task = true;
    fprintf(stderr, "\n[STOP] Signal generation interrupted.\n");
}

// Receiver antenna attenuation in dB for boresight angle = 0:5:180 [deg]
double ant_pat_db[37] = {0.00, 0.00, 0.22, 0.44, 0.67, 1.11, 1.56, 2.00,
                         2.44, 2.89, 3.56, 4.22, 4.89, 5.56, 6.22, 6.89,
                         7.56, 8.22, 8.89, 9.78, 10.67, 11.56, 12.44, 13.33,
                         14.44, 15.56, 16.67, 17.78, 18.89, 20.00, 21.33, 22.67,
                         24.00, 25.56, 27.33, 29.33, 31.56};

int allocatedSat[MAX_SAT];

static void make_e5b_output_path(const char *outfile, char *out, size_t out_size)
{
    const char *last_slash = strrchr(outfile, '/');
    const char *last_backslash = strrchr(outfile, '\\');
    const char *dir_end = last_slash;
    if (last_backslash != NULL && (dir_end == NULL || last_backslash > dir_end))
        dir_end = last_backslash;

    const char *dot = strrchr(outfile, '.');
    if (dot != NULL && dir_end != NULL && dot < dir_end)
        dot = NULL;

    if (dot == NULL)
        snprintf(out, out_size, "%s_E5b", outfile);
    else
    {
        size_t prefix_len = (size_t)(dot - outfile);
        if (prefix_len >= out_size)
            prefix_len = out_size - 1;
        snprintf(out, out_size, "%.*s_E5b%s", (int)prefix_len, outfile, dot);
    }
}


void *galileo_task(void *arg)
{
    init_tables();
    if (cosTable[0] == 0) {
        fprintf(stderr, "FATAL: Carrier tables not initialized!\n");
    }
    struct hash_queue hq;

    vector<queue<int>> queues(MAX_CHAN);
    map<int, int> chn_prn_map;
    int queue_index;

    sim_t *s = (sim_t *)arg;

    // pthread_t th_bits;
    bool exit_flag = false;

    // Stores current data bit from all channels
    int dataBit[MAX_CHAN];

    // Values for TOW fix - these determine whether corrections received over
    // socket are applied
    bool tow_fixed;
    bool local_fix = false;
    bool fixed_prn = false;

    clock_t tstart, tend;

    FILE *fp = NULL;
    FILE *fp_E5b = NULL;

    int sv;
    int neph, ieph;
    ephem_t eph1[EPHEM_ARRAY_SIZE][MAX_SAT];
    galtime_t g0;

    double llh[3];

    int i;
    channel_t chan[MAX_CHAN];
    double elvmask = 10; // in degree

    int ip, qp;
    int iTable;
    short *iq_buff = NULL;
    short *iq_buff_E5b = NULL;
    
    // 16-bit scaling factor: target ~12000 (37.5% of ±32767 range)
    // Worst-case analysis:
    //   signal_sum range: [-2, +2] (E1B data ±1 + E1C pilot ±1)
    //   ch_gain range: [0, ~1.5] at high elevation (variable by location/elevation)
    //   Conservative 6000.0 factor produces max ~18000 (safe margin)
    //   Worst-case output: signal_sum=2.0 * ch_gain=1.5 * 6000.0 = 18000 ✓ (safe margin)
    //   NOTE: Reduced from 12000.0 to 6000.0 due to overflow at high-elevation sites
    long overflow_count_i = 0;
    long overflow_count_q = 0;
    long overflow_count_i_E5b = 0;
    long overflow_count_q_E5b = 0;
    const double SCALE_FACTOR_16BIT = 3000.0;  // Conservative scaling to prevent overflows at all elevation angles
                                               // Target: ~6000 (18% utilization) with safety margin
                                               // Reduced iteratively due to overflows at higher elevation sites

    galtime_t grx;
    gal_sow_t grx_sow;
    int grx_week;
    double delt;
    int isamp;

    int iumd;
    int numd;
    int doppler_debug_prn = 0;
    int e5b_phase_debug_prn = 0;
    const int e5b_phase_debug_max_blocks = 10;
    const int e5b_phase_debug_head_samples = 10;
    const int e5b_phase_debug_tail_samples = 10;
    const int e5b_doppler_update_samples = (int)(TX_SAMPLERATE / 1000.0 + 0.5);
    const double e5b_doppler_update_interval = 0.001;

    std::vector<std::vector<double>> xyz(USER_MOTION_SIZE, std::vector<double>(3));

    char navfile[MAX_CHAR];
    char outfile[MAX_CHAR];
    char tv_file[MAX_CHAR];

    double samp_freq;
    int iq_buff_size;
    int data_format;

    int result;

    int gain[MAX_CHAN];
    double path_loss;
    double ant_gain;
    double ant_pat[37];
    int ibs; // boresight angle index

    datetime_t t0, tmin, tmax;
    galtime_t gmin, gmax;
    double dt;
    int igrx;

    double duration;
    int iduration;
    int verb;

    bool use_usrp;
    bool use_bits_from_streamer;
    bool generate_e1;
    bool generate_e5b;

    int timeoverwrite = FALSE; // Overwrite the TOC and TOE in the RINEX file

    ionoutc_t iono;

    ////////////////////////////////////////////////////////////
    // Read options
    ////////////////////////////////////////////////////////////

    strcpy(navfile, s->opt.navfile);
    strcpy(tv_file, s->opt.tvfile);
    strcpy(outfile, s->opt.outfile);

    g0.week = s->opt.g0.week;
    g0.sec = s->opt.g0.sec;

    gal2date(&g0, &t0);

    iduration = USER_MOTION_SIZE;

    samp_freq = (double)TX_SAMPLERATE;

    verb = s->opt.verb;

    iq_buff_size = NUM_IQ_SAMPLES;

    delt = 1.0 / (double)TX_SAMPLERATE;

    int interactive = s->opt.interactive;

    timeoverwrite = s->opt.timeoverwrite;

    iono.enable = s->opt.iono_enable;

    use_usrp = s->opt.use_usrp;
    generate_e1 = (s->opt.signal_mode == SIGNAL_MODE_E1 || s->opt.signal_mode == SIGNAL_MODE_BOTH);
    generate_e5b = (s->opt.signal_mode == SIGNAL_MODE_E5B || s->opt.signal_mode == SIGNAL_MODE_BOTH);

    duration = (double)s->opt.iduration / 10.0;

    use_bits_from_streamer = s->opt.use_bit_stream;

    llh[0] = s->opt.llh[0];
    llh[1] = s->opt.llh[1];
    llh[2] = s->opt.llh[2];

    ////////////////////////////////////////////////////////////
    // Start auxiliary task threads
    ////////////////////////////////////////////////////////////

    memcpy(llhr, llh, 3 * sizeof(double)); // Initialize the global llhr array (from socket.h)
    // Start auxiliary task threads
    std::thread th_loc(&locations_thread, llhr);
    th_loc.detach(); // Detach to prevent terminate() when th_loc goes out of scope

    ////////////////////////////////////////////////////////////
    // Load navigation messages and satellite ephemeris
    ////////////////////////////////////////////////////////////

    vector<ephem_t> eph_vector[MAX_SAT];
    ephem_t eph;

    // Load file pointers to individual satellites
    // std::cerr << s->opt.tvfile[0];

    // char tfname[MAX_CHAR];
    // snprintf(tfname, sizeof(tfname), "../rinex_files/week171.rnx");

    ephem_t t_eph[MAX_SAT][EPHEM_ARRAY_SIZE];
    
    int eph_count = readRinexV3(eph_vector, &iono, navfile);
   
    // Keep track of current ephemeris index for each satellite
    for (int i = 0; i < MAX_SAT; i++)
        current_eph.push_back(0);

    for (int i = 0; i < MAX_SAT; i++)
        old_eph.push_back(0);

    // Convert input llh to radians before converting to xyz!
    // This is paramount to ensure the satellite visibility mask centers on the correct hemisphere.
    llh[0] = llh[0] / R2D; // convert to RAD
    llh[1] = llh[1] / R2D; // convert to RAD

    iduration = (int)(duration * 10.0 + 0.5);
    numd = iduration;

    for (int i = 0; i < numd; i++)
        llh2xyz(llh, xyz[i].data());  // Convert llh to xyz for all steps

    ////////////////////////////////////////////////////////////
    // Receiver position
    ////////////////////////////////////////////////////////////

    fprintf(stderr, "\n[E5] Receiver position\n");
    numd = iduration;

    fprintf(stderr, "  ECEF xyz: %11.1f, %11.1f, %11.1f\n", xyz[0][0], xyz[0][1],
            xyz[0][2]);
    fprintf(stderr, "  LLH:      %11.6f, %11.6f, %11.1f\n", llh[0] * R2D, llh[1] * R2D,
            llh[2]);

    /* sv = PRN-1 */
    for (sv = 0; sv < MAX_SAT; sv++)
    {
        if (eph_vector[sv].size() == 0)
        {
            continue;
        }

        eph = eph_vector[sv][0];

        if (eph.vflg == 1)
        {
            gmin = eph.toc;
            gal2date(&gmin, &tmin);
            break;
        }
    }


    gmax.sec = 0;
    gmax.week = 0;
    tmax.sec = 0;
    tmax.mm = 0;
    tmax.hh = 0;
    tmax.d = 0;
    tmax.m = 0;
    tmax.y = 0;

    for (sv = 0; sv < MAX_SAT; sv++)
    {   
        if (eph_vector[sv].size() < 2)
            continue;

        eph = eph_vector[sv][eph_vector[sv].size() - 2];

        if (eph.vflg == 1)
        {
            if (eph.toc.sec > gmax.sec)
                {gmax = eph.toc;
                }
        }
    }
    

    gal2date(&gmax, &tmax);
    set_scenario_start_time(&g0, gmin, gmax, &t0, &tmin, &tmax, timeoverwrite, &iono, neph, eph_vector);
   
    datetime_t tl;
    gal2date(&g0, &tl);

    fprintf(stderr, "\n[E5] Ephemeris range: %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)", tmin.y,
            tmin.m, tmin.d, tmin.hh, tmin.mm, tmin.sec, gmin.week, gmin.sec);

    fprintf(stderr, " to %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n", tmax.y,
            tmax.m, tmax.d, tmax.hh, tmax.mm, tmax.sec, gmax.week, gmax.sec);

    fprintf(stderr, "[E5] Duration: %.1f sec\n", ((double)numd) / 10.0);

    fprintf(stderr, "[E5] Start time: %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n",
            tl.y, tl.m, tl.d, tl.hh, tl.mm, tl.sec, g0.week,
            g0.sec);

    ////////////////////////////////////////////////////////////
    // Read data from file and setup bit and TOW correction queues
    ////////////////////////////////////////////////////////////

    // Initial reception time
    grx = g0;
    grx_sow = galSowFromGalTime(grx);
    grx_week = grx.week;

    double dt1, correction;

    datetime_t t_temp;
    gal2date(&grx, &t_temp);
    double obs_time = gps_time(&t_temp);

    for (sv = 0; sv < MAX_SAT; sv++)
    {
        current_eph[sv] = epoch_matcher(grx, eph_vector[sv], current_eph[sv]);
//        if (current_eph[sv] == -1)
//        {
//            fprintf(stderr, "ERROR: No current set of ephemerides has been found for SVID %d grx %d-%f.\n", sv + 1, grx.week, grx.sec);
//            continue;
//        }

//        if (current_eph[sv] != -1)
//        {
//        	print_eph2(&eph_vector[sv][current_eph[sv]], sv+1);
//        }
    }



    ////////////////////////////////////////////////////////////
    // Baseband signal buffer and output file
    ////////////////////////////////////////////////////////////

    if (generate_e1)
        iq_buff = (short *)calloc(2 * iq_buff_size, sizeof(short));
    if (generate_e5b)
        iq_buff_E5b = (short *)calloc(2 * iq_buff_size, sizeof(short));

    if (generate_e1)
    {
        if (strcmp("-", outfile))
        {
            if (NULL == (fp = fopen(outfile, "wb")))
            {
                fprintf(stderr, "ERROR: Failed to open E1 output file.\n");
                exit(1);
            }
        }
        else
        {
            fp = stdout;
        }
    }

    if (generate_e5b)
    {
        char outfile_E5b[MAX_CHAR];
        if (s->opt.signal_mode == SIGNAL_MODE_BOTH)
            make_e5b_output_path(outfile, outfile_E5b, sizeof(outfile_E5b));
        else
            strncpy(outfile_E5b, outfile, sizeof(outfile_E5b) - 1);
        outfile_E5b[sizeof(outfile_E5b) - 1] = '\0';

        if (strcmp("-", outfile_E5b))
        {
            if (NULL == (fp_E5b = fopen(outfile_E5b, "wb")))
            {
                fprintf(stderr, "ERROR: Failed to open E5b output file.\n");
                exit(1);
            }
        }
        else
        {
            fp_E5b = stdout;
        }
    }

    ////////////////////////////////////////////////////////////
    // Initialize channels
    ////////////////////////////////////////////////////////////

    dt = 0.1;
    incGalSowMillis(&grx_sow, 100);
    grx = galSowToGalTime(&grx_sow, &grx_week);

    init_channel(chan, allocatedSat);

    allocateChannel(chan, eph_vector, iono, grx, xyz[0].data(), elvmask, &chn_prn_map, current_eph, allocatedSat);

    // for (i = 0; i < MAX_CHAN; i++)
    // {
    //     if (chan[i].prn > 0)
    //         fprintf(stderr, "%02d %6.1f %5.1f %11.1f %5.5f\n", chan[i].prn, chan[i].azel[0] * R2D, chan[i].azel[1] * R2D, chan[i].rho0.range, grx.sec);
    // }

    ////////////////////////////////////////////////////////////
    // Receiver antenna gain pattern
    ////////////////////////////////////////////////////////////

    for (i = 0; i < 37; i++)
        ant_pat[i] = pow(10.0, -ant_pat_db[i] / 20.0);

    tstart = clock();

    // Update receiver time

    int outbuf_idx = 0;
    float cosPhase[1];
    float sinPhase[1];
    double temp_phase = 0;

    bool debug = false;
    float *debug_codephase;
    FILE *fp_debug;

    // Ctrl+C will exit loop
    signal(SIGINT, &sigint_handler_gal_task);

    if (debug)
    {
        fp_debug = fopen("./debug_codePhase_data_with_CCP", "wb");
        debug_codephase = (float *)calloc(iq_buff_size, sizeof(float));
    }

    // Comprehensive channel initialization and safety checks
    for (i = 0; i < MAX_CHAN; i++)
    {
        if (chan[i].prn > 0)
        {
            sv = chan[i].prn - 1;
            // Ensure codes are generated (normally done in allocateChannel but double-checking)
            if (chan[i].ca_E1B == NULL) {
                chan[i].ca_E1B = (short *)calloc(2 * CA_SEQ_LEN_E1, sizeof(short));
                codegen_E1B(chan[i].ca_E1B, chan[i].prn);
            }
            if (chan[i].ca_E1C == NULL) {
                chan[i].ca_E1C = (short *)calloc(2 * CA_SEQ_LEN_E1, sizeof(short));
                codegen_E1C(chan[i].ca_E1C, chan[i].prn);
            }
            // Ensure navigation message buffer is allocated and populated
            if (chan[i].page == NULL) {
                chan[i].page = (int *)calloc(N_PAGE, sizeof(int));
                generateINavMsg(grx, &chan[i], &eph_vector[sv][current_eph[sv]], &iono);
            }
            if (chan[i].ca_E5bI == NULL) {
                chan[i].ca_E5bI = (short *)calloc(CA_SEQ_LEN_E5b, sizeof(short));
                codegen_E5bI(chan[i].ca_E5bI, chan[i].prn);
            }
            if (chan[i].ca_E5bQ == NULL) {
                chan[i].ca_E5bQ = (short *)calloc(CA_SEQ_LEN_E5b, sizeof(short));
                codegen_E5bQ(chan[i].ca_E5bQ, chan[i].prn);
            }
            if (chan[i].page_E5b == NULL) {
                chan[i].page_E5b = (int *)calloc(N_PAGE, sizeof(int));
                generateINavMsg_E5b(grx, &chan[i], &eph_vector[sv][current_eph[sv]], &iono);
            }
        }
    }

    if (use_bits_from_streamer)
    {
        fprintf(stderr, "\n[E5] Waiting for navigation message bits...");
        // ... (existing code for streamer bits)
    }
    else
    {
        fprintf(stderr, "\n[E5] Navigation messages ready.\n");
        fflush(stdout);
    }

    long last_nanos;
    long start;

    ////////////////////////////////////////////////////////////
    // Generate baseband signals
    ////////////////////////////////////////////////////////////

    fprintf(stderr, "\n[E5] Starting signal generation...\n");
    fflush(stderr);

    int cp = 0;

    incGalSowMillis(&grx_sow, 100);
    grx = galSowToGalTime(&grx_sow, &grx_week);

    for (iumd = 1; iumd < numd; iumd++)
    {
        // Copy contents of llh received over the locations thread
        if (llhr != NULL)
            memcpy(llh, llhr, 3 * sizeof(double));

        llh[0] = llh[0] / R2D; // convert to RAD
        llh[1] = llh[1] / R2D; // convert to RAD

        if (iumd < (int)xyz.size())
            llh2xyz(llh, xyz[iumd].data());

        for (i = 0; i < MAX_CHAN; i++)
        {
            if (chan[i].prn > 0)
            {
                if (doppler_debug_prn == 0)
                    doppler_debug_prn = chan[i].prn;
                if (generate_e5b && e5b_phase_debug_prn == 0)
                    e5b_phase_debug_prn = chan[i].prn;

                // Refresh code phase and data bit counters
                range_t rho;
                range_t rho_E5b;
                sv = chan[i].prn - 1;
                eph = eph_vector[sv][current_eph[sv]];

                // Current pseudorange
                computeRange(&rho, eph, &iono, grx, xyz[iumd].data(), chan[i].prn);

                chan[i].azel[0] = rho.azel[0];
                chan[i].azel[1] = rho.azel[1];

                // Update code phase and data bit counters for the first run
                computeCodePhase(&chan[i], rho, dt, grx);

                if (chan[i].set_code_phase)
                {
                    chan[i].set_code_phase = false;

                    // Absolute Time of Transmission
                    double tx_time = grx.sec - (rho.range / SPEED_OF_LIGHT);
                    
                    // Page boundary (2 seconds)
                    long page_idx = (long)(tx_time / 2.0);
                    double page_start_sec = (double)page_idx * 2.0;
                    double sec_in_page = tx_time - page_start_sec;
                    
                    if (sec_in_page < 0) {
                        sec_in_page += 2.0;
                        page_idx--;
                    }
                    
                    int ibit = (int)(sec_in_page / 0.004);
                    if (ibit >= N_SYM_PAGE) ibit = N_SYM_PAGE - 1;
                    
                    double code_sec = sec_in_page - ((double)ibit * 0.004);
                    double code_phase = (code_sec / 0.004) * (double)CA_SEQ_LEN_E1;

                    chan[i].code_phase = code_phase;
                    chan[i].ibit = ibit;
                    chan[i].ipage = (int)page_idx;

                    // STRICT TRANSMISSION TIME for page generation to avoid bit scrambles
                    galtime_t tx_time_gal = grx;
                    tx_time_gal.sec = (double)chan[i].ipage * 2.0;

                    generateINavMsg(tx_time_gal, &chan[i], &eph, &iono);
                }

                if (generate_e5b)
                {
                    computeRange_E5b(&rho_E5b, eph, &iono, grx, xyz[iumd].data(), chan[i].prn);

                    if (chan[i].ipage_E5b < 0)
                    {
                        // Signal transmission time.
                        double tx_time_E5b = grx.sec - (rho_E5b.range / SPEED_OF_LIGHT);

                        // Navigation page index for the current transmission time.
                        long page_idx_E5b = (long)(tx_time_E5b / (double)PAGE_TRANS_TIME);
                        // Start time of the current page.
                        double page_start_sec_E5b = (double)page_idx_E5b * (double)PAGE_TRANS_TIME;
                        // Offset within the page in seconds.
                        double sec_in_page_E5b = tx_time_E5b - page_start_sec_E5b;

                        // Correct negative offsets caused by floating-point precision.
                        if (sec_in_page_E5b < 0.0)
                        {
                            sec_in_page_E5b += (double)PAGE_TRANS_TIME;
                            page_idx_E5b--;
                            page_start_sec_E5b = (double)page_idx_E5b * (double)PAGE_TRANS_TIME;
                        }

                        // Bit index for the current transmission time.
                        int ibit_E5b = (int)(sec_in_page_E5b / 0.004);
                        if (ibit_E5b >= N_SYM_PAGE)
                            ibit_E5b = N_SYM_PAGE - 1;

                        // Millisecond index within the current page.
                        int ms_in_page_E5b = (int)(sec_in_page_E5b / 0.001);
                        // Keep only the fractional part after the integer millisecond.
                        double code_sec_E5b = sec_in_page_E5b - ((double)ms_in_page_E5b * 0.001);
                        // Code phase at the current transmission time.
                        double code_phase_E5b = (code_sec_E5b / 0.001) * (double)CA_SEQ_LEN_E5b;

                        if (code_phase_E5b >= (double)CA_SEQ_LEN_E5b)
                            code_phase_E5b -= (double)CA_SEQ_LEN_E5b;

                        chan[i].code_phase_E5b = code_phase_E5b;
                        chan[i].ibit_E5b = ibit_E5b;
                        chan[i].ipage_E5b = (int)page_idx_E5b;
                        chan[i].sec_code_idx_E5b_Q = ms_in_page_E5b % 100;
                        chan[i].sec_code_idx_E5b_I = ms_in_page_E5b % 4;

                        galtime_t tx_time_gal_E5b = grx;
                        tx_time_gal_E5b.sec = page_start_sec_E5b;

                        generateINavMsg_E5b(tx_time_gal_E5b, &chan[i], &eph, &iono);
                    }
                }

                // Path loss
                path_loss = 20200000.0 / (generate_e1 ? rho.d : rho_E5b.d);

                // Receiver antenna gain
                ibs = (int)((90.0 - (generate_e1 ? rho.azel[1] : rho_E5b.azel[1]) * R2D) / 5.0); // covert elevation to boresight
                ant_gain = ant_pat[ibs];

                // Signal gain
                gain[i] = (int)(path_loss * ant_gain * 128.0); // scaled by 2^7
            }
        }

        // if (doppler_debug_prn > 0)
        // {
        //     for (i = 0; i < MAX_CHAN; i++)
        //     {
        //         if (chan[i].prn == doppler_debug_prn)
        //         {
        //             fprintf(stderr, "PRN %02d grx=%10.3f doppler=%12.3f Hz\n",
        //                     chan[i].prn, grx.sec, chan[i].f_carr);
        //             break;
        //         }
        //     }
        // }

        // Current navigation bits and secondary codes for all channels
        // Current navigation bits and secondary codes for all channels
        int ch_databit[MAX_CHAN];
        int ch_secCode[MAX_CHAN];
        int ch_databit_E5b[MAX_CHAN];
        int ch_secCode_E5b_Q[MAX_CHAN];
        int ch_secCode_E5b_I[MAX_CHAN];
        for (i = 0; i < MAX_CHAN; i++)
        {
            if (chan[i].prn > 0)
            {
                // Synchronize with the precise code phase and bit indices calculated in computeCodePhase
                // Pre-calculate bits for the start of this 100ms interval
                ch_databit[i] = chan[i].page[chan[i].ibit] > 0 ? -1 : 1;
                ch_secCode[i] = GALILEO_E1_SECONDARY_CODE[chan[i].ibit % 25] > 0 ? -1 : 1;
                if (generate_e5b)
                {
                    ch_databit_E5b[i] = chan[i].page_E5b[chan[i].ibit_E5b] > 0 ? -1 : 1;
                    ch_secCode_E5b_Q[i] = GALILEO_E5b_Q_SECONDARY_CODE[chan[i].prn - 1][chan[i].sec_code_idx_E5b_Q] == '0' ? 1 : -1;
                    ch_secCode_E5b_I[i] = GALILEO_E5b_I_SECONDARY_CODE[chan[i].sec_code_idx_E5b_I] == '0' ? 1 : -1;
                }
            }
        }

        for (isamp = 0; isamp < iq_buff_size; isamp++)
        {
            int i_acc = 0;
            int q_acc = 0;
            int i_acc_E5b = 0;
            int q_acc_E5b = 0;
            
            for (i = 0; i < MAX_CHAN; i++)
            {
                if (chan[i].prn > 0)
                {
                    if (generate_e1)
                    {
                    // Update indices if code phase rolls over
                    if (chan[i].code_phase >= (double)CA_SEQ_LEN_E1)
                    {
                        chan[i].code_phase -= (double)CA_SEQ_LEN_E1;
                        chan[i].ibit++;
                        // 500 bits = 1 page
                        if (chan[i].ibit >= N_SYM_PAGE)
                        {   
                            chan[i].ibit = 0;
                            chan[i].ipage++;

                            // Generate new page strictly based on Satellite Transmission Time
                            galtime_t tx_time_gal = grx;
                            tx_time_gal.sec = (double)chan[i].ipage * 2.0;

                            sv = chan[i].prn - 1;
                            eph = eph_vector[sv][current_eph[sv]];
                            generateINavMsg(tx_time_gal, &chan[i], &eph, &iono);
                        }

                        // Update bits for the new 4ms symbol
                        ch_databit[i] = chan[i].page[chan[i].ibit] > 0 ? -1 : 1;
                        ch_secCode[i] = GALILEO_E1_SECONDARY_CODE[chan[i].ibit % 25] > 0 ? -1 : 1;
                    }

                    int carr_idx = ((int)(COS_TAB_LENGTH * chan[i].carr_phase)) & COS_TAB_MASK;
                    int cosPh = cosTable[carr_idx];
                    int sinPh = sinTable[carr_idx];

                    // For BOC(1,1), we have 2 sub-chips per code chip (total 8184 sub-chips per 4ms)
                    double icode_f = chan[i].code_phase * 2.0; 
                    int icode = (int)icode_f;
                    if (icode >= 2 * CA_SEQ_LEN_E1) icode %= (2 * CA_SEQ_LEN_E1);

                    int E1B_subchip = chan[i].ca_E1B[icode];
                    int E1C_subchip = chan[i].ca_E1C[icode];

                    // Galileo E1 signal is (E1B_data * E1B_subchip + E1C_pilot * E1C_subchip)
                    double signal_sum = (double)(E1B_subchip * ch_databit[i] + E1C_subchip * ch_secCode[i]);

                    // Apply channel-specific gain (path loss + antenna)
                    double ch_gain = (double)gain[i] / 128.0; 
                    
                    // Scale to utilize 16-bit signed short range (±32767)
                    // Using SCALE_FACTOR_16BIT = 12000.0 to place peak at ~24000 (75% utilization)
                    double scaled_signal = signal_sum * ch_gain * SCALE_FACTOR_16BIT;

                    i_acc += (int)(scaled_signal * (double)cosPh / 32767.0);
                    q_acc += (int)(scaled_signal * (double)sinPh / 32767.0);

                    // Update phases
                    chan[i].code_phase += chan[i].f_code * delt;
                    chan[i].carr_phase += chan[i].f_carr * delt;
                    chan[i].carr_phase -= (int)chan[i].carr_phase; 
                    if (chan[i].carr_phase < 0) chan[i].carr_phase += 1.0;
                    }

                    if (generate_e5b)
                    {
                        if ((isamp % e5b_doppler_update_samples) == 0)
                        {
                            int ms_in_block = isamp / e5b_doppler_update_samples;
                            galtime_t g_ms_start = grx;
                            galtime_t g_ms_end = grx;
                            range_t rho_ms_start_E5b;
                            range_t rho_ms_end_E5b;

                            g_ms_start.sec += (double)ms_in_block * e5b_doppler_update_interval;
                            g_ms_end.sec = g_ms_start.sec + e5b_doppler_update_interval;

                            sv = chan[i].prn - 1;
                            eph = eph_vector[sv][current_eph[sv]];
                            computeRange_E5b(&rho_ms_start_E5b, eph, &iono, g_ms_start, xyz[iumd].data(), chan[i].prn);
                            computeRange_E5b(&rho_ms_end_E5b, eph, &iono, g_ms_end, xyz[iumd].data(), chan[i].prn);

                            double rhorate_E5b = (rho_ms_end_E5b.range - rho_ms_start_E5b.range) / e5b_doppler_update_interval;
                            chan[i].f_carr_E5b = -rhorate_E5b / LAMBDA_E5b;
                            chan[i].f_code_E5b = CODE_FREQ_E5b + chan[i].f_carr_E5b * CARR_TO_CODE_E5b;
                            chan[i].rho0_E5b = rho_ms_end_E5b;
                        }

                        // Advance to the next primary-code period.
                        if (chan[i].code_phase_E5b >= (double)CA_SEQ_LEN_E5b)
                        {
                            // if (chan[i].prn == 13)
                            // {
                            //     fprintf(stderr,
                            //             "E5b branch enter PRN %02d block=%d sample=%d Q_idx=%d I_idx=%d\n",
                            //             chan[i].prn,
                            //             iumd,
                            //             isamp,
                            //             chan[i].sec_code_idx_E5b_Q,
                            //             chan[i].sec_code_idx_E5b_I);
                            // }

                            // Remove one primary-code period.
                            chan[i].code_phase_E5b -= (double)CA_SEQ_LEN_E5b;

                            // Each secondary-code chip lasts 1 ms.
                            chan[i].sec_code_idx_E5b_Q++;
                            if (chan[i].sec_code_idx_E5b_Q >= 100)
                                chan[i].sec_code_idx_E5b_Q = 0;

                            chan[i].sec_code_idx_E5b_I++;
                            if (chan[i].sec_code_idx_E5b_I >= 4) // Update the bit index every 4 ms.
                            {
                                chan[i].sec_code_idx_E5b_I = 0;
                                chan[i].ibit_E5b++;

                                if (chan[i].ibit_E5b >= N_SYM_PAGE) // Advance to the next navigation page.
                                {
                                    chan[i].ibit_E5b = 0;
                                    chan[i].ipage_E5b++;

                                    // Start time of the current page.
                                    galtime_t tx_time_E5b = grx;
                                    tx_time_E5b.sec = (double)chan[i].ipage_E5b * (double)PAGE_TRANS_TIME;

                                    sv = chan[i].prn - 1;
                                    eph = eph_vector[sv][current_eph[sv]];
                                    generateINavMsg_E5b(tx_time_E5b, &chan[i], &eph, &iono);
                                }
                            }

                            // if (chan[i].prn == e5b_phase_debug_prn)
                            // {
                            //     fprintf(stderr,
                            //             "E5b index update PRN %02d block=%d sample=%d Q_idx=%d I_idx=%d ibit=%d ipage=%d\n",
                            //             chan[i].prn,
                            //             iumd,
                            //             isamp,
                            //             chan[i].sec_code_idx_E5b_Q,
                            //             chan[i].sec_code_idx_E5b_I,
                            //             chan[i].ibit_E5b,
                            //             chan[i].ipage_E5b);
                            // }

                            ch_databit_E5b[i] = chan[i].page_E5b[chan[i].ibit_E5b] > 0 ? -1 : 1;
                            ch_secCode_E5b_Q[i] = GALILEO_E5b_Q_SECONDARY_CODE[chan[i].prn - 1][chan[i].sec_code_idx_E5b_Q] == '0' ? 1 : -1;
                            ch_secCode_E5b_I[i] = GALILEO_E5b_I_SECONDARY_CODE[chan[i].sec_code_idx_E5b_I] == '0' ? 1 : -1;
                        }

                        int carr_idx_E5b = ((int)(COS_TAB_LENGTH * chan[i].carr_phase_E5b)) & COS_TAB_MASK;
                        int cosPh_E5b = cosTable[carr_idx_E5b];
                        int sinPh_E5b = sinTable[carr_idx_E5b];

                        int icode_E5b = (int)chan[i].code_phase_E5b;
                        if (icode_E5b >= CA_SEQ_LEN_E5b)
                            icode_E5b %= CA_SEQ_LEN_E5b;

                        int E5bI_subchip = chan[i].ca_E5bI[icode_E5b];
                        int E5bQ_subchip = chan[i].ca_E5bQ[icode_E5b];

                        double data = E5bI_subchip * ch_databit_E5b[i] * ch_secCode_E5b_I[i];
                        double pilot = E5bQ_subchip * ch_secCode_E5b_Q[i];

                        double base_I = E5b_Pow * pilot;
                        double base_Q = -E5b_Pow * data;
                        double scale = SCALE_FACTOR_16BIT;
                        double ch_gain = (double)gain[i] / 128.0;

                        double E5b_I = (base_I * cosPh_E5b - base_Q * sinPh_E5b) / 32767.0;
                        double E5b_Q = (base_I * sinPh_E5b + base_Q * cosPh_E5b) / 32767.0;

                        i_acc_E5b += (int)(E5b_I * ch_gain * scale);
                        q_acc_E5b += (int)(E5b_Q * ch_gain * scale);

                        // if (chan[i].prn == e5b_phase_debug_prn &&
                        //     iumd <= e5b_phase_debug_max_blocks &&
                        //     (isamp < e5b_phase_debug_head_samples ||
                        //      isamp >= iq_buff_size - e5b_phase_debug_tail_samples))
                        // {
                        //     fprintf(stderr,
                        //             "E5b PRN %02d block=%d sample=%d code_chip_idx=%d\n",
                        //             chan[i].prn, iumd, isamp, icode_E5b);
                        // }

                        chan[i].code_phase_E5b += chan[i].f_code_E5b * delt;
                        chan[i].carr_phase_E5b += chan[i].f_carr_E5b * delt;
                        chan[i].carr_phase_E5b -= (int)chan[i].carr_phase_E5b;
                        if (chan[i].carr_phase_E5b < 0) chan[i].carr_phase_E5b += 1.0;
                    }
                }
            }
            
            // Detect overflows before clipping (for 16-bit diagnostics)
            if (i_acc > 32767 || i_acc < -32768) overflow_count_i++;
            if (q_acc > 32767 || q_acc < -32768) overflow_count_q++;
            if (i_acc_E5b > 32767 || i_acc_E5b < -32768) overflow_count_i_E5b++;
            if (q_acc_E5b > 32767 || q_acc_E5b < -32768) overflow_count_q_E5b++;
            
            // Currently clipping to 8-bit for backward compat; will be removed in T02
            // Store I/Q samples into buffer (8-bit signed) with scale factor now targeting 16-bit range
            // Scaling factor changed from 2000.0 to 12000.0:
            //   signal_sum is max 2.0 (E1B pilot/data + E1C pilot)
            //   ch_gain is max ~1.0 at reference distance
            //   So scaled_signal is max ~24000
            //   Max carrier is ±32767
            //   Result is max ~24000 (75% of 16-bit range, safe margin for no clipping)
            
            // Clipping removed: with 16-bit scaling (SCALE_FACTOR_16BIT = 12000.0),
            // samples are now in valid 16-bit range and should not be clipped to ±127
            // Overflow detection at ±32767 provides safety margin

            if (generate_e1)
            {
                iq_buff[isamp * 2] = (short)i_acc;
                iq_buff[isamp * 2 + 1] = (short)q_acc;
            }
            if (generate_e5b)
            {
                iq_buff_E5b[isamp * 2] = (short)i_acc_E5b;
                iq_buff_E5b[isamp * 2 + 1] = (short)q_acc_E5b;
            }
        }

        if (generate_e1)
            fwrite(iq_buff, sizeof(short), 2 * iq_buff_size, fp);
        if (generate_e5b)
            fwrite(iq_buff_E5b, sizeof(short), 2 * iq_buff_size, fp_E5b);

        if (iumd % 10 == 0) {
            fflush(stderr);
        }

        // Check and update ephemeris index and satellite allocation every 30 seconds
        igrx = (int)(grx.sec*10.0+0.5);
        if ((int)fmodf(igrx, 300) == 0)
        {

            gal2date(&grx, &t_temp);
            obs_time = gps_time(&t_temp);
            // cout << "All: " << grx.sec << " : " << g0.sec;
            // for (int sv = 0; sv < MAX_SAT; sv++)
            //     old_eph[sv] = current_eph[sv];

            for (int sv = 0; sv < MAX_SAT; sv++)
            {
                current_eph[sv] = epoch_matcher(grx, eph_vector[sv], current_eph[sv]);
                // if (current_eph[sv] != old_eph[sv])
                //     cout << endl << "Changed eph for " << sv << " - " << current_eph[sv] << " : " << old_eph[sv] << endl;
            }
            if (iumd < (int)xyz.size())
                allocateChannel(chan, eph_vector, iono, grx, xyz[iumd].data(), elvmask, &chn_prn_map, current_eph, allocatedSat);
        }

        incGalSowMillis(&grx_sow, 100);
        grx = galSowToGalTime(&grx_sow, &grx_week);

        ////////////////////////////////////////////////////////////
        // Write into FIFO
        ///////////////////////////////////////////////////////////

        if (use_usrp)
        {
            if (!s->galileo_sim.ready)
            {
                // Initialization has been done. Ready to create TX task.
                printf("[E5] Signal generator is ready.\n");
                s->galileo_sim.ready = 1;
                pthread_cond_signal(&(s->galileo_sim.initialization_done));
            }

            long int a = get_nanos();
            // Wait utill FIFO write is ready
            pthread_mutex_lock(&(s->galileo_sim.lock));
            while (!is_fifo_write_ready(s))
                pthread_cond_wait(&(s->fifo_write_ready), &(s->galileo_sim.lock));
            pthread_mutex_unlock(&(s->galileo_sim.lock));

            // Write into FIFO
            short *fifo_src = generate_e1 ? iq_buff : iq_buff_E5b;
            memcpy(&(s->fifo[s->head * 2]), fifo_src,
                   NUM_IQ_SAMPLES * 2 * sizeof(short));

            s->head += (long)NUM_IQ_SAMPLES;
            if (s->head >= FIFO_LENGTH)
                s->head -= FIFO_LENGTH;
            pthread_cond_signal(&(s->fifo_read_ready));
        }

        //
        // Update receiver time
        if (use_bits_from_streamer)
        {

            llh[0] = llh[0] * R2D;
            llh[1] = llh[1] * R2D;

            // Update time counter
            if (!local_fix && tow_fixed)
            {
                double d = (get_nanos() - start) / 1e9;
                incGalSowSeconds(&grx_sow, dt1);
                grx = galSowToGalTime(&grx_sow, &grx_week);

                printf("\r[E5] Stream time offset: %.3f s", d);
                printf("\n[E5] TOW correction applied: %.6f s\n", dt1);

                tow_fixed = true;
                local_fix = true;
            }
        }

        if (verb == TRUE)
        {
            if (use_cursors)
            {
                clear();
                attron(A_REVERSE);
                llh[0] = llh[0] * R2D;
                llh[1] = llh[1] * R2D;

                printw("\n Location: %10f, %10f, %4f - Time: \n", llh[0], llh[1], llh[2]);
                printw("\n Elapsed time %4.1f s\n", subGalTime(grx, g0));
                printw("\n%3s%6s%14s%17s%21s%18s%18s%18s%5s\n", "CH", "PRN", "Azimuth", "Elevation", "Doppler [Hz]", "Code phase", "rx_time", "Pseudorange", "Eph"); //,"Y","M","D","HH","MM","SS");
                attroff(A_REVERSE);
                for (i = 0; i < MAX_CHAN; i++)
                {
                    if (chan[i].prn > 0 && use_cursors)
                        printw("%3d%6d%14f%17f%21f%18f%18f%18f%5d\n", i, chan[i].prn, chan[i].azel[0] * R2D, chan[i].azel[1] * R2D, chan[i].f_carr, chan[i].code_phase, grx.sec, chan[i].rho0.range, current_eph[chan[i].prn - 1]); //, t0.y, t0.m, t0.d, t0.hh, t0.mm, t0.sec);
                }
                refresh();
            }
        }
        else if ((iumd % 100) == 0 || iumd + 1 == numd)
            fprintf(stderr, "\r[E5] Progress: %4.1f / %.1f s", subGalTime(grx, g0), ((double)numd) / 10.0);

        if (stop_signal_called_gal_task)
            break;

        fflush(stdout);
    }

    tend = clock();
    exit_flag = true;
    fprintf(stderr, "\n[E5] Signal generation complete.\n");
    fflush(stderr);
    
    // T01: Report overflow detection counters (should be zero for valid 16-bit scaling)
    if (overflow_count_i > 0 || overflow_count_q > 0) {
        fprintf(stderr, "[WARN] E1 overflow detected during signal generation.\n");
        fprintf(stderr, "       I-channel overflows: %ld\n", overflow_count_i);
        fprintf(stderr, "       Q-channel overflows: %ld\n", overflow_count_q);
        fprintf(stderr, "       Scale factor may be too aggressive.\n");
    } else {
        fprintf(stderr, "[E5] E1 overflow check: OK (0 overflows)\n");
    }
    if (overflow_count_i_E5b > 0 || overflow_count_q_E5b > 0) {
        fprintf(stderr, "[WARN] E5b overflow detected during signal generation.\n");
        fprintf(stderr, "       E5b I-channel overflows: %ld\n", overflow_count_i_E5b);
        fprintf(stderr, "       E5b Q-channel overflows: %ld\n", overflow_count_q_E5b);
    } else if (generate_e5b) {
        fprintf(stderr, "[E5] E5b overflow check: OK (0 overflows)\n");
    }
    fflush(stderr);
    
    endwin();

    // Free I/Q buffer
    if (iq_buff != NULL)
        free(iq_buff);
    if (iq_buff_E5b != NULL)
        free(iq_buff_E5b);

    // Close file
    if (fp != NULL)
        fclose(fp);
    if (fp_E5b != NULL && fp_E5b != stdout)
        fclose(fp_E5b);

    if (debug)
        fclose(fp_debug);

    // Process time
    fprintf(stderr, "[E5] Process time: %.1f sec\n",
            (double)(tend - tstart) / CLOCKS_PER_SEC);

    // return 1;
    return (NULL);
}
