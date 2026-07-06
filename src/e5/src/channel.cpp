#include "../include/galileo-sdr.h"

void init_channel(channel_t *chan, int *allocatedSat)
{
    // Clear all channels
    memset(chan, 0, MAX_CHAN * sizeof(channel_t));
    for (int i = 0; i < MAX_CHAN; i++)
    {
        chan[i].set_code_phase = true;
        chan[i].prn = 0;
        chan[i].ca_E1B = NULL;
        chan[i].ca_E1C = NULL;
        chan[i].ca_E5bI = NULL;
        chan[i].ca_E5bQ = NULL;
        chan[i].page = NULL;
        chan[i].page_E5b = NULL;
        chan[i].ipage_E5b = -1;
    }

    // Clear satellite allocation flag
    for (int sv = 0; sv < MAX_SAT; sv++)
        allocatedSat[sv] = -1;
}

int allocateChannel(channel_t *chan,
                    // map<int, vector<Rinex3Nav::DataGAL>> *navGAL,
                    vector<ephem_t> *eph_vector,
                    ionoutc_t ionoutc, galtime_t grx, double *xyz,
                    double elvMask, map<int, int> *sm,
                    vector<int> current_eph,
                    int *allocatedSat)
{
    int nsat = 0;
    int i, sv;
    double azel[2];

    range_t rho;
    range_t rho_E5b;
    double ref[3] = {0.0};
    double r_ref, r_xyz;
    double phase_ini;
    double r_ref_E5b, r_xyz_E5b;
    double phase_ini_E5b;

    ephem_t eph;

    for (sv = 0; sv < MAX_SAT; sv++)
    {
        if (!eph_vector[sv].size())
            continue;

        if (!eph_vector[sv][0].vflg)
            continue;
        
        datetime_t t;
        gal2date(&grx, &t);
        double obs_time = gps_time(&t);
        
        current_eph[sv] = epoch_matcher(grx, eph_vector[sv], current_eph[sv]);
        // cout << "here : " << current_eph[sv] << " : " << sv<<  endl;
        if (current_eph[sv] < 0)
            continue;

        eph = eph_vector[sv][current_eph[sv]];
        //cout << "GRX: " << grx.sec << " : " << current_eph[sv] << endl;

        if (checkSatVisibility(eph, grx, xyz, 10, azel, sv + 1) == 1)
        {
            nsat++; // Number of visible satellites

            if (allocatedSat[sv] == -1) // Visible but not allocated
            {
                // Allocated new satellite
                for (i = 0; i < MAX_CHAN; i++)
                {
                    if (chan[i].prn == 0)
                    {
                        // Initialize channel
                        chan[i].prn = sv + 1;
                        chan[i].azel[0] = azel[0];
                        chan[i].azel[1] = azel[1];
                        chan[i].g0 = grx;

                        // Insert latest channel assignment to the map
                        // sm->insert({chan[i].prn, i});

                        // C/A code generation
                        if (chan[i].ca_E1B == NULL)
                            chan[i].ca_E1B = (short *)calloc(2 * CA_SEQ_LEN_E1, sizeof(short));
                        codegen_E1B(chan[i].ca_E1B, chan[i].prn);

                        if (chan[i].ca_E1C == NULL)
                            chan[i].ca_E1C = (short *)calloc(2 * CA_SEQ_LEN_E1, sizeof(short));
                        codegen_E1C(chan[i].ca_E1C, chan[i].prn);

                        // Generate the primary codes for both components.
                        if (chan[i].ca_E5bI == NULL)
                            chan[i].ca_E5bI = (short *)calloc(CA_SEQ_LEN_E5b, sizeof(short));
                        codegen_E5bI(chan[i].ca_E5bI, chan[i].prn);

                        if (chan[i].ca_E5bQ == NULL)
                            chan[i].ca_E5bQ = (short *)calloc(CA_SEQ_LEN_E5b, sizeof(short));
                        codegen_E5bQ(chan[i].ca_E5bQ, chan[i].prn);

                        // Generate navigation messages for the reception time.
                        // Generate navigation message
                        //generateNavMsg(grx, &chan[i], advance_fptr);
                        if (chan[i].page == NULL)
                            chan[i].page = (int *)calloc(N_PAGE, sizeof(int));
                        generateINavMsg(grx, &chan[i], &eph, &ionoutc);

                        if (chan[i].page_E5b == NULL)
                            chan[i].page_E5b = (int *)calloc(N_PAGE, sizeof(int));
                        generateINavMsg_E5b(grx, &chan[i], &eph, &ionoutc);

                        // Initialize pseudorange
                        computeRange(&rho, eph, &ionoutc, grx, xyz, chan[i].prn);
                        chan[i].rho0 = rho;

                        // Initialize carrier phase
                        r_xyz = rho.range;

                        // Compute the initial carrier phase.
                        computeRange(&rho, eph, &ionoutc, grx, ref, chan[i].prn);
                        r_ref = rho.range;

                        phase_ini = (2.0 * r_ref - r_xyz) / LAMBDA_L1;
                        chan[i].carr_phase = phase_ini - floor(phase_ini);

                        computeRange_E5b(&rho_E5b, eph, &ionoutc, grx, xyz, chan[i].prn);
                        chan[i].rho0_E5b = rho_E5b;
                        r_xyz_E5b = rho_E5b.range;

                        computeRange_E5b(&rho_E5b, eph, &ionoutc, grx, ref, chan[i].prn);
                        r_ref_E5b = rho_E5b.range;

                        phase_ini_E5b = (2.0 * r_ref_E5b - r_xyz_E5b) / LAMBDA_E5b;
                        chan[i].carr_phase_E5b = phase_ini_E5b - floor(phase_ini_E5b);
                        chan[i].ipage_E5b = -1;
                        chan[i].ibit_E5b = 0;
                        chan[i].sec_code_idx_E5b_Q = 0;
                        chan[i].sec_code_idx_E5b_I = 0;
                        chan[i].code_phase_E5b = 0.0;

                        //print_eph2(&eph, sv+1);
                        break;
                    }
                }

                // Set satellite allocation channel
                if (i < MAX_CHAN)
                    allocatedSat[sv] = i;
            }
        }
        else if (allocatedSat[sv] >= 0) // Not visible but allocated
        {
            // Clear channel
            chan[allocatedSat[sv]].prn = 0;

            // Clear satellite allocation flag
            allocatedSat[sv] = -1;
        }
    }
    // advance_fptr = true;
    return (nsat);
}
