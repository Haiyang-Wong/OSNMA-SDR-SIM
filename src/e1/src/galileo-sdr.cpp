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


int samples_per_code;
std::vector<int> current_eph;
std::vector<int> old_eph;
vector<int> candidate_prn_list;
vector<int> initial_visible_prn_list;
vector<int> all_visible_prn_list;
vector<int> future_added_prn_list;
vector<vector<int>> nearest_prn_by_visible;
map<int, int> prn_first_seen_tow;
map<int, int> prn_last_seen_tow;
vector<pair<int, vector<int>>> tow_visible_prn_sets;
int visibility_scan_step_seconds = 30;
int initial_visible_prn_flag[MAX_SAT + 1] = {0};

bool stop_signal_called_gal_task = false;
bool location_thread_started = false;
const bool enable_runtime_location_updates = false;

bool use_cursors = true;
bool advance_fptr = false;
vector<page_struct> nav_page;
vector<page_struct> nav_page_2;

FILE *TV_ptr;
FILE *TV_ptrs_sv[MAX_SAT];
int prnTable[40]={0};

static void start_location_thread_once(double llh[3])
{
    if (!enable_runtime_location_updates || location_thread_started)
        return;

    std::thread th_loc(locations_thread, llh);
    th_loc.detach();
    location_thread_started = true;
}

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
int getPrnIndex(int a)
{
    int n[]={1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,20 ,21 ,22 ,23 ,24 ,25 ,26 ,27 ,28 ,29 ,30 ,31 ,32 ,33 ,34 ,35,36,37};

    for(int i=0;i<40;i++)
    {
        if(n[i]==a)
        {
            return i;
        }
    }
    return -1;
}

static bool contains_prn(const vector<int> &prns, int prn)
{
    return find(prns.begin(), prns.end(), prn) != prns.end();
}

static void scan_visible_prns_for_duration(vector<ephem_t> eph_vector[MAX_SAT],
                                           const vector<int> &base_current_eph,
                                           galtime_t g0,
                                           double duration,
                                           double *xyz,
                                           double elvmask)
{
    map<int, int> first_seen;
    map<int, int> last_seen;
    vector<pair<int, vector<int>>> visible_sets;
    vector<pair<int, int>> future_seen_order;

    int duration_seconds = static_cast<int>(ceil(duration));
    vector<int> scan_offsets;
    for (int offset = 0; offset <= duration_seconds; offset += visibility_scan_step_seconds)
        scan_offsets.push_back(offset);
    if (scan_offsets.empty() || scan_offsets.back() != duration_seconds)
        scan_offsets.push_back(duration_seconds);

    for (size_t si = 0; si < scan_offsets.size(); si++)
    {
        galtime_t gt = incGalTime(g0, static_cast<double>(scan_offsets[si]));
        vector<int> scan_current_eph = base_current_eph;
        vector<int> visible_now;

        for (int sv = 0; sv < MAX_SAT; sv++)
        {
            if (!eph_vector[sv].size())
                continue;

            scan_current_eph[sv] = epoch_matcher(gt, eph_vector[sv], scan_current_eph[sv]);
            if (scan_current_eph[sv] < 0)
                continue;

            ephem_t eph = eph_vector[sv][scan_current_eph[sv]];
            double azel[2] = {0.0, 0.0};
            if (checkSatVisibility(eph, gt, xyz, elvmask, azel, sv + 1) == 1)
            {
                int prn = sv + 1;
                visible_now.push_back(prn);
                int tow = static_cast<int>(gt.sec);
                if (first_seen.find(prn) == first_seen.end())
                {
                    first_seen[prn] = tow;
                    if (!contains_prn(initial_visible_prn_list, prn))
                        future_seen_order.push_back(make_pair(tow, prn));
                }
                last_seen[prn] = tow;
            }
        }

        sort(visible_now.begin(), visible_now.end());
        visible_sets.push_back(make_pair(static_cast<int>(gt.sec), visible_now));
    }

    all_visible_prn_list = initial_visible_prn_list;
    int fallback_start_tow = static_cast<int>(g0.sec);
    for (size_t i = 0; i < initial_visible_prn_list.size(); i++)
    {
        int prn = initial_visible_prn_list[i];
        if (first_seen.find(prn) == first_seen.end())
            first_seen[prn] = fallback_start_tow;
        if (last_seen.find(prn) == last_seen.end())
            last_seen[prn] = fallback_start_tow;
    }
    sort(future_seen_order.begin(), future_seen_order.end());
    future_added_prn_list.clear();
    for (size_t i = 0; i < future_seen_order.size(); i++)
    {
        int prn = future_seen_order[i].second;
        if (!contains_prn(all_visible_prn_list, prn))
        {
            all_visible_prn_list.push_back(prn);
            future_added_prn_list.push_back(prn);
        }
    }

    prn_first_seen_tow = first_seen;
    prn_last_seen_tow = last_seen;
    tow_visible_prn_sets = visible_sets;
}
string findPage(double tow,int prn)
{
    for(int i=0;i<nav_page.size();i++)
    {
        int tow1=static_cast<int>(nav_page[i].tow);
        int tow2=static_cast<int>(tow);
        if(tow1==tow2&&prn==nav_page[i].prn)
        {
            return nav_page[i].navMsg;
        }
    }
    return "";
}
void strToInt(string s, int index,int* n)
{

    for(int i = 0;i<120;i++)
    {
        n[i]=s[index+i]-'0';
    }
}
void galileo_task1(void *arg)
{
    struct hash_queue hq;

    vector<queue<int>> queues;
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

    FILE *fp;

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
    signed char *iq8_buff = NULL;

    galtime_t grx;
    double delt;
    int isamp;

    int iumd;
    int numd;

    vector<array<double, 3>> xyz;

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

    duration = s->opt.iduration / 10.0;
    iduration = (int)(duration * 10.0 + 0.5);
    numd = iduration;
    xyz.resize(max(numd, 1));

    use_bits_from_streamer = s->opt.use_bit_stream;

    llh[0] = s->opt.llh[0];
    llh[1] = s->opt.llh[1];
    llh[2] = s->opt.llh[2];

    ////////////////////////////////////////////////////////////
    // Start auxiliary task threads
    ////////////////////////////////////////////////////////////

    // Update location via UDP Socket
    start_location_thread_once(llh);

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

    llh[0] = llh[0] / R2D; // convert to RAD
    llh[1] = llh[1] / R2D; // convert to RAD
    llh2xyz(llh, xyz[0].data());  // Convert llh to xyz

    ////////////////////////////////////////////////////////////
    // Receiver position
    ////////////////////////////////////////////////////////////

    fprintf(stderr, "\n[E1] Receiver position\n");
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
    fprintf(stderr, "\n[E1] Ephemeris range: %4d/%02d/%02d,%02d:%02d:%02.0f",
            tmin.y, tmin.m, tmin.d, tmin.hh, tmin.mm, tmin.sec);

    fprintf(stderr, " to %4d/%02d/%02d,%02d:%02d:%02.0f\n", tmax.y,
            tmax.m, tmax.d, tmax.hh, tmax.mm, tmax.sec);
    ////////////////////////////////////////////////////////////
    // Read data from file and setup bit and TOW correction queues
    ////////////////////////////////////////////////////////////

    // Initial reception time
    grx = g0;

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

    // Allocate I/Q buffer
    iq_buff = (short *)calloc(2 * iq_buff_size, sizeof(short));

    // Open output file
    // "-" can be used as name for stdout
    if (strcmp("-", outfile))
    {
        if (NULL == (fp = fopen(outfile, "wb")))
        {
            fprintf(stderr, "ERROR: Failed to open output file.\n");
            exit(1);
        }
    }
    else
    {
        fp = stdout;
    }

    ////////////////////////////////////////////////////////////
    // Initialize channels
    ////////////////////////////////////////////////////////////

    dt = 0.10000002314200000;
    grx = incGalTime(grx, dt);

    init_channel(chan, allocatedSat);

    allocateChannel(chan, eph_vector, iono, grx, xyz[0].data(), elvmask, &chn_prn_map, current_eph, allocatedSat);

    candidate_prn_list.clear();
    initial_visible_prn_list.clear();
    nearest_prn_by_visible.clear();
    for (int flag_i = 0; flag_i <= MAX_SAT; flag_i++)
        initial_visible_prn_flag[flag_i] = 0;

    vector<array<double, 3>> candidate_pos;
    for (sv = 0; sv < MAX_SAT; sv++)
    {
        if (!eph_vector[sv].size())
            continue;
        if (current_eph[sv] < 0)
            continue;

        double pos[3], vel[3], clk[3];
        eph = eph_vector[sv][current_eph[sv]];
        satpos(eph, grx, pos, vel, clk);
        array<double, 3> sat_pos = {{pos[0], pos[1], pos[2]}};
        candidate_prn_list.push_back(eph.svid);
        candidate_pos.push_back(sat_pos);
    }

    for (i = 0; i < MAX_CHAN; i++)
    {
        if (chan[i].prn > 0)
        {
            initial_visible_prn_list.push_back(chan[i].prn);
            if (chan[i].prn <= MAX_SAT)
                initial_visible_prn_flag[chan[i].prn] = 1;
        }
    }

    scan_visible_prns_for_duration(eph_vector, current_eph, g0, duration, xyz[0].data(), elvmask);

    for (size_t v = 0; v < all_visible_prn_list.size(); v++)
    {
        int prn = all_visible_prn_list[v];
        if (contains_prn(candidate_prn_list, prn))
            continue;

        int sv_idx = prn - 1;
        if (sv_idx < 0 || sv_idx >= MAX_SAT || !eph_vector[sv_idx].size())
            continue;

        int first_seen_tow = static_cast<int>(g0.sec);
        map<int, int>::const_iterator seen_it = prn_first_seen_tow.find(prn);
        if (seen_it != prn_first_seen_tow.end())
            first_seen_tow = seen_it->second;

        galtime_t first_seen_time = incGalTime(g0, static_cast<double>(first_seen_tow) - g0.sec);
        int eph_idx = current_eph[sv_idx];
        eph_idx = epoch_matcher(first_seen_time, eph_vector[sv_idx], eph_idx);
        if (eph_idx < 0)
            continue;

        double pos[3], vel[3], clk[3];
        eph = eph_vector[sv_idx][eph_idx];
        satpos(eph, first_seen_time, pos, vel, clk);
        array<double, 3> sat_pos = {{pos[0], pos[1], pos[2]}};
        candidate_prn_list.push_back(prn);
        candidate_pos.push_back(sat_pos);
    }

    for (size_t v = 0; v < all_visible_prn_list.size(); v++)
    {
        int base_prn = all_visible_prn_list[v];
        int base_idx = -1;
        for (size_t c = 0; c < candidate_prn_list.size(); c++)
        {
            if (candidate_prn_list[c] == base_prn)
            {
                base_idx = static_cast<int>(c);
                break;
            }
        }

        vector<pair<double, int>> distances;
        if (base_idx >= 0)
        {
            for (size_t c = 0; c < candidate_prn_list.size(); c++)
            {
                if (candidate_prn_list[c] == base_prn)
                    continue;

                double dx = candidate_pos[base_idx][0] - candidate_pos[c][0];
                double dy = candidate_pos[base_idx][1] - candidate_pos[c][1];
                double dz = candidate_pos[base_idx][2] - candidate_pos[c][2];
                distances.push_back(make_pair(dx * dx + dy * dy + dz * dz, candidate_prn_list[c]));
            }
        }

        sort(distances.begin(), distances.end());
        vector<int> nearest_row;
        for (size_t n = 0; n < distances.size(); n++)
            nearest_row.push_back(distances[n].second);
        nearest_prn_by_visible.push_back(nearest_row);
    }

    cout << "[E1] Candidate PRNs: [";
    for (size_t c = 0; c < candidate_prn_list.size(); c++)
    {
        cout << candidate_prn_list[c];
        if (c + 1 < candidate_prn_list.size())
            cout << ",";
    }
    cout << "]" << endl;

    cout << "[E1] Initially visible PRNs: [";
    for (size_t v = 0; v < initial_visible_prn_list.size(); v++)
    {
        cout << initial_visible_prn_list[v];
        if (v + 1 < initial_visible_prn_list.size())
            cout << ",";
    }
    cout << "]" << endl;

    cout << "[E1] Visible PRNs during simulation: [";
    for (size_t v = 0; v < all_visible_prn_list.size(); v++)
    {
        cout << all_visible_prn_list[v];
        if (v + 1 < all_visible_prn_list.size())
            cout << ",";
    }
    cout << "]" << endl;

    cout << "[E1] Future-added PRNs: [";
    for (size_t v = 0; v < future_added_prn_list.size(); v++)
    {
        cout << future_added_prn_list[v];
        if (v + 1 < future_added_prn_list.size())
            cout << ",";
    }
    cout << "]" << endl;

    cout << "[E1] Cross-auth neighbor table prepared." << endl;

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

    if (use_bits_from_streamer)
    {
        fprintf(stderr, "\n[E1] Waiting for navigation message bits");

        while (1)
        {
            bool flag = false;
            for (int i = 0; i < MAX_CHAN; i++)
            {
                if (!queues[i].empty())
                {
                    flag = true;
                    break;
                }
            }
            if (flag)
            {
                break;
            }
            else
            {
                sleep(1);
                fprintf(stderr, ".");
            }
        }

        
        bool reset_index = false;
        fflush(stdout);
    }
    else
    {
       
        fflush(stdout);
    }

    long last_nanos;
    long start;

    ////////////////////////////////////////////////////////////
    // Generate baseband signals
    ////////////////////////////////////////////////////////////

    initscr();

    int cp = 0;

    grx = incGalTime(grx, dt);

    for (iumd = 1; iumd < numd; iumd++)
    {
        start = get_nanos();

        // Copy contents of llh received over the locations thread
        memcpy(llh, llhr, 3 * sizeof(double));

        llh[0] = llh[0] / R2D; // convert to RAD
        llh[1] = llh[1] / R2D; // convert to RAD

        llh2xyz(llh, xyz[iumd].data());

	        for (i = 0; i < MAX_CHAN; i++)
	        {
	            if (chan[i].prn > 0)
	            {

	                // Refresh code phase and data bit counters
	                range_t rho;
	                sv = chan[i].prn - 1;
	                if (sv < 0 || sv >= MAX_SAT || current_eph[sv] < 0 ||
	                    current_eph[sv] >= static_cast<int>(eph_vector[sv].size()))
	                {
	                    double elapsed = subGalTime(grx, g0);
	                    if (elapsed >= 6500.0 && elapsed <= 6600.0)
	                    {
	                        fprintf(stderr,
	                                "\n[SignalGen] skip invalid eph: elapsed=%.1f tow=%.1f ch=%d prn=%d eph_idx=%d eph_count=%zu\n",
	                                elapsed, grx.sec, i, chan[i].prn,
	                                (sv >= 0 && sv < MAX_SAT) ? current_eph[sv] : -999,
	                                (sv >= 0 && sv < MAX_SAT) ? eph_vector[sv].size() : 0);
	                    }
	                    chan[i].prn = 0;
	                    continue;
	                }
	                eph = eph_vector[sv][current_eph[sv]];

	                // Current pseudorange
	                computeRange(&rho, eph, &iono, grx, xyz[iumd].data(), chan[i].prn);

                chan[i].azel[0] = rho.azel[0];
                chan[i].azel[1] = rho.azel[1];

                // Update code phase and data bit counters for the first run
                computeCodePhase(&chan[i], rho, dt, grx);

                // Path loss
                path_loss = 20200000.0 / rho.d;

	                // Receiver antenna gain
	                ibs = (int)((90.0 - rho.azel[1] * R2D) / 5.0); // covert elevation to boresight
	                if (ibs < 0)
	                    ibs = 0;
	                if (ibs > 36)
	                    ibs = 36;
	                ant_gain = ant_pat[ibs];

                // Signal gain
                gain[i] = (int)(path_loss * ant_gain * 128.0); // scaled by 2^7
            }
        }

	
            // fprintf(stderr, "\n");
    if(static_cast<int>(grx.sec)%2!=0)
    {
        for (sv = 0; sv < MAX_SAT; sv++)
        {
            if (!eph_vector[sv].size())
                continue;

            if (current_eph[sv] < 0)
                continue;

            eph = eph_vector[sv][current_eph[sv]];
            int index=getPrnIndex(eph.svid);

            if(index >= 0 && prnTable[index]<static_cast<int>(grx.sec))//prnTable[index]<static_cast<int>(g.sec)
            {
                channel_t nav_chan = {};
                nav_chan.prn = eph.svid;
                generateINavMsg(grx, &nav_chan, &eph, &iono);
                prnTable[index]=static_cast<int>(grx.sec);
            }
        }
    }

        // Check and update ephemeris index and satellite allocation every 30 seconds
        igrx = (int)(grx.sec*10.0+0.5);
        if ((int)fmodf(igrx, 300) == 0)
        {

            gal2date(&grx, &t_temp);
            obs_time = gps_time(&t_temp);
            for (int sv = 0; sv < MAX_SAT; sv++)
            {
                current_eph[sv] = epoch_matcher(grx, eph_vector[sv], current_eph[sv]);

            }
            allocateChannel(chan, eph_vector, iono, grx, xyz[iumd].data(), elvmask, &chn_prn_map, current_eph, allocatedSat);
        }

        grx = incGalTime(grx, dt);

        // Update receiver time
        if (use_bits_from_streamer)
        {

            llh[0] = llh[0] * R2D;
            llh[1] = llh[1] * R2D;

            // Update time counter
            if (!local_fix && tow_fixed)
            {
                double d = (get_nanos() - start) / 1e9;
                grx = incGalTime(grx, dt1);

                printf("\r[E1] Stream time offset: %.3f s", d);
                printf("\n[E1] TOW correction applied: %.6f s\n", dt1);

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

                
                attroff(A_REVERSE);
                for (i = 0; i < MAX_CHAN; i++)
                {
                    if (chan[i].prn > 0 && use_cursors)
                        printw("%3d%6d%14f%17f%21f%18f%18f%18f%5d\n", i, chan[i].prn, chan[i].azel[0] * R2D, chan[i].azel[1] * R2D, chan[i].f_carr, chan[i].code_phase, grx.sec, chan[i].rho0.range, current_eph[chan[i].prn - 1]); //, t0.y, t0.m, t0.d, t0.hh, t0.mm, t0.sec);
                }
                refresh();
            }
        }
        else ;
         //   fprintf(stderr, "\rTime into run = %4.1f - %4.1ld", subGalTime(grx, g0), clock() - tstart);

        if (stop_signal_called_gal_task)
            break;

        fflush(stdout);
    }

    tend = clock();
    exit_flag = true;

    endwin();

    // Free I/Q buffer
    free(iq_buff);

    // Close file
    fclose(fp);

    if (debug)
        fclose(fp_debug);

    return;
}
void *galileo_task(void *arg)
{
    struct hash_queue hq;

    vector<queue<int>> queues;
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

    FILE *fp;

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
    signed char *iq8_buff = NULL;

    galtime_t grx;
    double delt;
    int isamp;

    int iumd;
    int numd;

    vector<array<double, 3>> xyz;

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

    duration = s->opt.iduration / 10.0;
    iduration = (int)(duration * 10.0 + 0.5);
    numd = iduration;
    xyz.resize(max(numd, 1));

    use_bits_from_streamer = s->opt.use_bit_stream;

    llh[0] = s->opt.llh[0];
    llh[1] = s->opt.llh[1];
    llh[2] = s->opt.llh[2];

    ////////////////////////////////////////////////////////////
    // Start auxiliary task threads
    ////////////////////////////////////////////////////////////

    // Update location via UDP Socket
    start_location_thread_once(llh);
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

    llh[0] = llh[0] / R2D; // convert to RAD
    llh[1] = llh[1] / R2D; // convert to RAD
    llh2xyz(llh, xyz[0].data());  // Convert llh to xyz

    ////////////////////////////////////////////////////////////
    // Receiver position
    ////////////////////////////////////////////////////////////

    
numd = iduration;
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

    

    fprintf(stderr, "\n[E1] Duration: %.1f sec\n", ((double)numd) / 10.0);

 

    ////////////////////////////////////////////////////////////
    // Read data from file and setup bit and TOW correction queues
    ////////////////////////////////////////////////////////////

    // Initial reception time
    grx = g0;

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

    // Allocate I/Q buffer
    iq_buff = (short *)calloc(2 * iq_buff_size, sizeof(short));

    // Open output file
    // "-" can be used as name for stdout
    if (strcmp("-", outfile))
    {
        if (NULL == (fp = fopen(outfile, "wb")))
        {
            fprintf(stderr, "ERROR: Failed to open output file.\n");
            exit(1);
        }
    }
    else
    {
        fp = stdout;
    }

    ////////////////////////////////////////////////////////////
    // Initialize channels
    ////////////////////////////////////////////////////////////

    dt = 0.10000002314200000;
    grx = incGalTime(grx, dt);

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

    if (use_bits_from_streamer)
    {
        fprintf(stderr, "\n[E1] Waiting for navigation message bits");

        while (1)
        {
            bool flag = false;
            for (int i = 0; i < MAX_CHAN; i++)
            {
                if (!queues[i].empty())
                {
                    flag = true;
                    break;
                }
            }
            if (flag)
            {
                break;
            }
            else
            {
                sleep(1);
                fprintf(stderr, ".");
            }
        }

        fprintf(stderr, "\n[E1] Navigation bits received; starting generator.\n");
        bool reset_index = false;
        fflush(stdout);
    }
    else
    {
        fprintf(stderr, "\n[E1] Generating Galileo navigation signals.\n");
        fflush(stdout);
    }

    long last_nanos;
    long start;

    ////////////////////////////////////////////////////////////
    // Generate baseband signals
    ////////////////////////////////////////////////////////////

    initscr();

    int cp = 0;

    grx = incGalTime(grx, dt);

    for (iumd = 1; iumd < numd; iumd++)
    {
        start = get_nanos();

        // Copy contents of llh received over the locations thread
        memcpy(llh, llhr, 3 * sizeof(double));

        llh[0] = llh[0] / R2D; // convert to RAD
        llh[1] = llh[1] / R2D; // convert to RAD

        llh2xyz(llh, xyz[iumd].data());

        for (i = 0; i < MAX_CHAN; i++)
        {
            if (chan[i].prn > 0)
            {

                // Refresh code phase and data bit counters
                range_t rho;
                sv = chan[i].prn - 1;
                eph = eph_vector[sv][current_eph[sv]];

                // Current pseudorange
                computeRange(&rho, eph, &iono, grx, xyz[iumd].data(), chan[i].prn);

                chan[i].azel[0] = rho.azel[0];
                chan[i].azel[1] = rho.azel[1];

                // Update code phase and data bit counters for the first run
                computeCodePhase(&chan[i], rho, dt, grx);

                // Path loss
                path_loss = 20200000.0 / rho.d;

                // Receiver antenna gain
                ibs = (int)((90.0 - rho.azel[1] * R2D) / 5.0); // covert elevation to boresight
                ant_gain = ant_pat[ibs];

                // Signal gain
                gain[i] = (int)(path_loss * ant_gain * 128.0); // scaled by 2^7
            }
        }

        for (isamp = 0; isamp < iq_buff_size; isamp++)
        {
            int i_acc = 0;
            int q_acc = 0;
            double a = get_nanos();
            // fprintf(stderr, "\n");
            for (i = 0; i < MAX_CHAN; i++)
            {
                if (chan[i].prn > 0)
                {
                    if (chan[i].code_phase >= CA_SEQ_LEN_E1)
                    {
                        chan[i].code_phase -= CA_SEQ_LEN_E1;
                        chan[i].ibit++;

                        // 500 bits = 1 page
                        if (chan[i].ibit >= N_SYM_PAGE)
                        {   
                            chan[i].ibit = 0;
                            chan[i].ipage++;

	                            // Generate new page
	                            sv = chan[i].prn - 1;
	                            if (sv < 0 || sv >= MAX_SAT || current_eph[sv] < 0 ||
	                                current_eph[sv] >= static_cast<int>(eph_vector[sv].size()))
	                            {
	                                double elapsed = subGalTime(grx, g0);
	                                if (elapsed >= 6500.0 && elapsed <= 6600.0)
	                                {
	                                    fprintf(stderr,
	                                            "\n[SignalGen] skip page update invalid eph: elapsed=%.1f tow=%.1f ch=%d prn=%d eph_idx=%d eph_count=%zu\n",
	                                            elapsed, grx.sec, i, chan[i].prn,
	                                            (sv >= 0 && sv < MAX_SAT) ? current_eph[sv] : -999,
	                                            (sv >= 0 && sv < MAX_SAT) ? eph_vector[sv].size() : 0);
	                                }
	                                chan[i].prn = 0;
	                                continue;
	                            }
	                            eph = eph_vector[sv][current_eph[sv]];
	                            // the key is chan[i].page=page
	                            string ss= findPage(grx.sec,eph.svid);
	                            if (ss.size() < 240)
	                            {
	                                double elapsed = subGalTime(grx, g0);
	                                if (elapsed >= 6500.0 && elapsed <= 6600.0)
	                                {
	                                    fprintf(stderr,
	                                            "\n[SignalGen] missing nav page: elapsed=%.1f tow=%.1f ch=%d prn=%d eph_svid=%d page_len=%zu\n",
	                                            elapsed, grx.sec, i, chan[i].prn, eph.svid, ss.size());
	                                }
	                                chan[i].prn = 0;
	                                continue;
	                            }
	                            //cout<<"\n\r"<<grx.sec<<","<<eph.svid<<","<<ss<<"\n\r"<<endl;
	    			     int * od=(int*)malloc(sizeof(int)*120);
	   		   	     int * ev=(int*)malloc(sizeof(int)*120);
   			     strToInt(ss,0,ev);
      			     strToInt(ss,120,od);
      			     int *page = (int *)malloc(PAGE_SIZE * sizeof(int));
      			     generateFrame(ev, &page[0]);
          	 	     generateFrame(od, &page[250]);
                            chan[i].page=page;
                        }
                    }

                    int cosPh = cosTable512[((int)(511 * chan[i].carr_phase)) & 511];
                    int sinPh = sinTable512[((int)(511 * chan[i].carr_phase)) & 511];

                    int icode = (int)(chan[i].code_phase * 2);

                    int E1B_chip = chan[i].ca_E1B[icode];
                    int E1C_chip = chan[i].ca_E1C[icode];

                    int databit = chan[i].page[chan[i].ibit] > 0 ? -1 : 1;
                    int secCode = GALILEO_E1_SECONDARY_CODE[chan[i].ibit % 25] > 0 ? -1 : 1;

                    ip = (E1B_chip * databit - E1C_chip * secCode) * cosPh; // * gain[i];
                    qp = (E1B_chip * databit - E1C_chip * secCode) * sinPh; // * gain[i];

                    // Accumulate for all visible satellites
                    i_acc += ip;
                    q_acc += qp;

                    // Update code phase
                    chan[i].code_phase += chan[i].f_code * delt;

                    // Update carrier phase
                    chan[i].carr_phase += (chan[i].f_carr) * delt;
                    chan[i].carr_phase -= (long)chan[i].carr_phase; // (carr_phase %1)
                }
            }
            // Store I/Q samples into buffer
            iq_buff[isamp * 2] = (short)i_acc;
            iq_buff[isamp * 2 + 1] = (short)q_acc;
            // advance_fptr = true;
        }


        fwrite(iq_buff, sizeof(short), 2 * iq_buff_size, fp);

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
            allocateChannel(chan, eph_vector, iono, grx, xyz[iumd].data(), elvmask, &chn_prn_map, current_eph, allocatedSat);
        }

        grx = incGalTime(grx, dt);

        ////////////////////////////////////////////////////////////
        // Write into FIFO
        ///////////////////////////////////////////////////////////

        if (use_usrp)
        {
            if (!s->galileo_sim.ready)
            {
                // Initialization has been done. Ready to create TX task.
                printf("[E1] Signal generator is ready.\n");
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
            memcpy(&(s->fifo[s->head * 2]), iq_buff,
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
                grx = incGalTime(grx, dt1);

                printf("\r[E1] Stream time offset: %.3f s", d);
                printf("\n[E1] TOW correction applied: %.6f s\n", dt1);

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
            fprintf(stderr, "\r[E1] Progress: %4.1f / %.1f s", subGalTime(grx, g0), ((double)numd) / 10.0);

        if (stop_signal_called_gal_task)
            break;

        fflush(stdout);
    }

    tend = clock();
    exit_flag = true;
    fprintf(stderr, "\n[E1] Signal generation complete.\n");
    endwin();

    // Free I/Q buffer
    free(iq_buff);

    // Close file
    fclose(fp);

    if (debug)
        fclose(fp_debug);

    // Process time
    fprintf(stderr, "[E1] Process time: %.1f sec\n",
            (double)(tend - tstart) / CLOCKS_PER_SEC);

    // return 1;
    return (NULL);
}
