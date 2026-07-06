#define _CRT_SECURE_NO_WARNINGS

#include "../include/galileo-sdr.h"
#include <math.h>
#include <unistd.h>
#include <vector>
#include<malloc.h>
#include<iostream>
#include<string>
#include <unordered_set>
#include<queue>
#include <cstring>
#include <cstdint>
#include<fstream>
#include <iomanip>
#include <sstream>
#include <bitset>
#include <unordered_map>
#include <algorithm>
#include <bitset>
#include<stdio.h>
#include <cerrno>
#include <climits>
#include <cstdlib>
#include <cctype>
#include <sys/wait.h>

#define BIT_ISSET(var, n) !!((long)(var) & (1 << (n)))

#define EXECUTE_OR_GOTO(label, ...) \
    if (__VA_ARGS__)                \
    {                               \
        return_code = EXIT_FAILURE; \
        goto label;                 \
    }

using namespace std;

int durring_time=0;
int numSub = 0;
string PATH="";
int neb_num = 3;  // test
int cross_auth_mode = 1;
vector<int> no_osnma_prns;
vector<int> selected_focus_prns;
int subNub=10;
const string mapping = "0123456789ABCDEF";
double start_time = 1000000000.0;
vector<page_struct> nav_page_all;
int cur_wn = 0;
bool stop_signal_called = false;
    int HF=0;
    int MF=0;
    int dic=0;
vector<vector<string>> mltSequence;
int ks_num[9]={96,104,112,120,128,160,192,224,256};
int ts_num[5]={20,24,28,32,40};
int KS;
int TS;
int KL;
int TL;
void sigint_handler(int code)
{
    (void)code;
    stop_signal_called = true;
    fprintf(stderr, "\n[STOP] Interrupted by user.\n");
}

long int samples_consumed = 0l;

void init_sim(sim_t *s)
{
   

    pthread_mutex_init(&(s->tx.lock), NULL);
    // s->tx.error = 0;

    pthread_mutex_init(&(s->galileo_sim.lock), NULL);
    // s->galileo_sim.error = 0;
    s->galileo_sim.ready = 0;

    pthread_cond_init(&(s->galileo_sim.initialization_done), NULL);

    s->status = 0;
    s->head = 0;
    s->tail = 0;
    s->sample_length = 0;

    pthread_cond_init(&(s->fifo_write_ready), NULL);

    pthread_cond_init(&(s->fifo_read_ready), NULL);

    s->time = 0.0;
}

void *tx_task(void *arg)
{
    sim_t *s = (sim_t *)arg;
    size_t samples_populated;
    size_t num_samps_sent = 0;
    size_t samples_per_buffer = SAMPLES_PER_BUFFER;
    fprintf(stderr, "\n[TX] Transmission task started.\n");
    // sleep(1);
    int k = 0;

    while (1)
    {
        if (stop_signal_called)
            goto out;

        short *tx_buffer_current = s->tx.buffer;
        int buffer_samples_remaining = SAMPLES_PER_BUFFER;
        const void **buffs_ptr = NULL;

        while (buffer_samples_remaining > 0)
        {
            pthread_mutex_lock(&(s->galileo_sim.lock));

            while (get_sample_length(s) == 0)
            {
                pthread_cond_wait(&(s->fifo_read_ready), &(s->galileo_sim.lock));
            }

            samples_populated = fifo_read(tx_buffer_current, buffer_samples_remaining, s);
            pthread_mutex_unlock(&(s->galileo_sim.lock));

            pthread_cond_signal(&(s->fifo_write_ready));

            // Advance the buffer pointer.
            buffer_samples_remaining -= (unsigned int)samples_populated;

            tx_buffer_current += (2 * samples_populated);
        }

        k++;

        buffs_ptr = (const void **)tx_buffer_current;

        // vector<double> tx_buffer_vector; //(tx_buffer_current, tx_buffer_current + SAMPLES_PER_BUFFER);
        //  vector<complex<float>> fc_buffer;
        //  for(int i=0; i < SAMPLES_PER_BUFFER; i=i+2)
        //  {
        //      fc_buffer.push_back(complex<float>(s->tx.buffer[i], s->tx.buffer[i+1]));// * 2], tx_buffer_current[i * 2 + 1]) );
        //      //fprintf(stderr, "\nF - %f/ %d", s->tx.buffer[i], i);
        //  }

        // uhd_tx_streamer_send(s->tx.stream, s->tx.buffer_ptr, SAMPLES_PER_BUFFER, &s->tx.md, 1000, &samples_populated);
        // printf("\nHere %ld\n", tx_buffer_vector.size());

        size_t num_tx_samps = s->tx.stream->send(s->tx.buffer, SAMPLES_PER_BUFFER, s->tx.md, 1000);
        samples_consumed = samples_consumed + num_tx_samps;
        // fprintf(stderr, "\nSent %d", k);
        // tx_buffer_vector.clear();
	//usleep(10000);
        if (is_fifo_write_ready(s))
        {
            fprintf(stderr, "\r[TX] elapsed %.1f s, buffers sent %d", s->time, k);
            s->time += 0.1;
            fflush(stdout);
        }
        else if (is_finished_generation(s))
        {
            goto out;
        }
    }
out:
    return NULL;
}

int start_tx_task(sim_t *s)
{
    int status;

    status = pthread_create(&(s->tx.thread), NULL, tx_task, s);

    return (status);
}

int start_galileo_task(sim_t *s)
{
    int status;

    status = pthread_create(&(s->galileo_sim.thread), NULL, galileo_task, s);

    return (status);
}

void usage(char *progname)
{
    printf(
"Usage: %s -f <config.ini>\n"
        "Options:\n"
        "  -f <config.ini>  Run configuration file\n",
        progname);

    return;
}

string trimConfigValue(const string &value)
{
    size_t first = 0;
    while (first < value.size() && isspace(static_cast<unsigned char>(value[first]))) {
        first++;
    }

    size_t last = value.size();
    while (last > first && isspace(static_cast<unsigned char>(value[last - 1]))) {
        last--;
    }

    return value.substr(first, last - first);
}

bool parseConfigInt(const map<string, string> &config, const string &key, int &value)
{
    map<string, string>::const_iterator it = config.find(key);
    if (it == config.end()) {
        cerr << "ERROR: Missing required config field: " << key << endl;
        return false;
    }

    char *endptr = NULL;
    errno = 0;
    long parsed = strtol(it->second.c_str(), &endptr, 10);
    if (errno != 0 || endptr == it->second.c_str() || *endptr != '\0' ||
        parsed < INT_MIN || parsed > INT_MAX) {
        cerr << "ERROR: Invalid integer for config field " << key << ": " << it->second << endl;
        return false;
    }

    value = static_cast<int>(parsed);
    return true;
}

bool parseOptionalConfigBool(const map<string, string> &config,
                             const string &key,
                             bool defaultValue,
                             bool &value)
{
    map<string, string>::const_iterator it = config.find(key);
    if (it == config.end()) {
        value = defaultValue;
        return true;
    }

    string parsed = trimConfigValue(it->second);
    transform(parsed.begin(), parsed.end(), parsed.begin(),
              [](unsigned char c) { return static_cast<char>(tolower(c)); });

    if (parsed == "true") {
        value = true;
        return true;
    }
    if (parsed == "false") {
        value = false;
        return true;
    }

    cerr << "ERROR: Invalid boolean for config field " << key
         << ": " << it->second << ". Expected true or false." << endl;
    return false;
}

bool requireConfigString(const map<string, string> &config, const string &key, string &value)
{
    map<string, string>::const_iterator it = config.find(key);
    if (it == config.end() || it->second.empty()) {
        cerr << "ERROR: Missing required config field: " << key << endl;
        return false;
    }

    value = it->second;
    return true;
}

bool parsePrnListText(const string &text, vector<int> &prns)
{
    string normalized = text;
    for (size_t i = 0; i < normalized.size(); i++) {
        if (normalized[i] == '[' || normalized[i] == ']' || normalized[i] == ',') {
            normalized[i] = ' ';
        }
    }

    stringstream ss(normalized);
    int prn;
    vector<int> parsed;
    while (ss >> prn) {
        parsed.push_back(prn);
    }
    if (!ss.eof()) {
        return false;
    }

    sort(parsed.begin(), parsed.end());
    parsed.erase(unique(parsed.begin(), parsed.end()), parsed.end());
    prns = parsed;
    return true;
}

bool parseOptionalConfigPrnList(const map<string, string> &config, const string &key, vector<int> &prns)
{
    map<string, string>::const_iterator it = config.find(key);
    if (it == config.end()) {
        prns.clear();
        return true;
    }

    string value = trimConfigValue(it->second);
    if (value.empty() || value == "[]") {
        prns.clear();
        return true;
    }

    if (!parsePrnListText(value, prns)) {
        cerr << "ERROR: Invalid PRN list for config field " << key << ": " << it->second << endl;
        return false;
    }
    return true;
}

string formatPrnVector(const vector<int> &prns)
{
    ostringstream out;
    for (size_t i = 0; i < prns.size(); i++) {
        if (i > 0) {
            out << ",";
        }
        out << prns[i];
    }
    return out.str();
}

bool containsPrn(const vector<int> &prns, int prn)
{
    return find(prns.begin(), prns.end(), prn) != prns.end();
}

bool isNoOsnmaPrn(int prn)
{
    return containsPrn(no_osnma_prns, prn) || containsPrn(selected_focus_prns, prn);
}

bool isOsnmaEnabledPrn(int prn)
{
    if (cross_auth_mode != 1 &&
        !all_visible_prn_list.empty() &&
        !containsPrn(all_visible_prn_list, prn)) {
        return false;
    }
    return !isNoOsnmaPrn(prn);
}

bool canBeCrossAuthTargetPrna(int prn, int mode)
{
    if (mode == 1) {
        return true;
    }
    return isOsnmaEnabledPrn(prn);
}

bool canBeCrossAuthSourcePrnd(int prna, int prnd, int mode)
{
    if (mode == 1) {
        return true;
    }
    if (prna == prnd) {
        return false;
    }
    if (mode == 2) {
        return true;
    }
    if (mode == 3) {
        return !isOsnmaEnabledPrn(prnd);
    }
    return false;
}

vector<int> defaultFocusPrnsFromSelected(const vector<int> &selectedPrns)
{
    vector<int> sortedPrns = selectedPrns;
    sort(sortedPrns.begin(), sortedPrns.end());
    vector<int> result;
    for (size_t i = 1; i < sortedPrns.size(); i += 2) {
        result.push_back(sortedPrns[i]);
    }
    return result;
}

bool validateSelectedPrns(const vector<int> &inputPrns, const vector<int> &selectedPrns)
{
    if (inputPrns.empty()) {
        cout << "[ERROR] At least one PRN must be selected." << endl;
        return false;
    }
    for (size_t i = 0; i < inputPrns.size(); i++) {
        if (!containsPrn(selectedPrns, inputPrns[i])) {
            cout << "[ERROR] PRN " << inputPrns[i] << " is not available in the selected PRN list." << endl;
            return false;
        }
    }
    return true;
}

void promptFocusPrns(const vector<int> &selectedPrns)
{
    if (cross_auth_mode == 1) {
        selected_focus_prns.clear();
        no_osnma_prns.clear();
        return;
    }

    if (no_osnma_prns.empty()) {
        selected_focus_prns = defaultFocusPrnsFromSelected(selectedPrns);
        cout << "[OSNMA] no_osnma_prns not set; using default disabled PRNd candidates: ["
             << formatPrnVector(selected_focus_prns) << "]" << endl;
        return;
    }

    selected_focus_prns = no_osnma_prns;
    cout << "[OSNMA] Disabled PRNd candidates: ["
         << formatPrnVector(selected_focus_prns) << "]" << endl;
}

bool validateNoOsnmaPrnsAgainstSelected(const vector<int> &selectedPrns)
{
    if (cross_auth_mode == 1) {
        return true;
    }

    for (size_t i = 0; i < no_osnma_prns.size(); i++) {
        if (!containsPrn(selectedPrns, no_osnma_prns[i]) &&
            !containsPrn(candidate_prn_list, no_osnma_prns[i])) {
            cerr << "ERROR: no_osnma_prns contains PRN " << no_osnma_prns[i]
                 << ", but it is not in selected or candidate PRN list." << endl;
            return false;
        }
    }
    return true;
}

void printCrossAuthModeSummary(const vector<int> &selectedPrns)
{
    vector<int> enabledPrns;
    vector<int> disabledPrns;
    for (size_t i = 0; i < selectedPrns.size(); i++) {
        int prn = selectedPrns[i];
        if (cross_auth_mode != 1 && isNoOsnmaPrn(prn)) {
            disabledPrns.push_back(prn);
        } else {
            enabledPrns.push_back(prn);
        }
    }

    cout << "\n== Cross Authentication ==" << endl;
    cout << "Mode: " << cross_auth_mode << endl;
    cout << "OSNMA enabled PRNs:  " << formatPrnVector(enabledPrns) << endl;
    cout << "OSNMA disabled PRNs: " << formatPrnVector(disabledPrns) << endl;
}

bool readRunConfig(const string &path, map<string, string> &config)
{
    ifstream in(path.c_str());
    if (!in.is_open()) {
        cerr << "ERROR: Failed to open config file: " << path << endl;
        return false;
    }

    string line;
    int lineNumber = 0;
    while (getline(in, line)) {
        lineNumber++;
        string trimmed = trimConfigValue(line);
        if (trimmed.empty() || trimmed[0] == '#' || trimmed[0] == ';') {
            continue;
        }

        size_t equals = trimmed.find('=');
        if (equals == string::npos) {
            cerr << "ERROR: Invalid config line " << lineNumber << ": " << line << endl;
            return false;
        }

        string key = trimConfigValue(trimmed.substr(0, equals));
        string value = trimConfigValue(trimmed.substr(equals + 1));
        if (key.empty() || value.empty()) {
            cerr << "ERROR: Invalid config line " << lineNumber << ": " << line << endl;
            return false;
        }

        config[key] = value;
    }

    return true;
}

bool applyRunConfig(const string &configPath, sim_t &s, double &duration, bool &use_usrp, bool &enable_e5)
{
    map<string, string> config;
    if (!readRunConfig(configPath, config)) {
        return false;
    }

    string location;
    string startTime;
    string osnmaFile;
    string rinexFile;
    string outputFile;
    if (!requireConfigString(config, "location", location) ||
        !requireConfigString(config, "start_time", startTime) ||
        !requireConfigString(config, "osnma_file", osnmaFile) ||
        !requireConfigString(config, "rinex_file", rinexFile) ||
        !requireConfigString(config, "output_file", outputFile)) {
        return false;
    }

    double lat;
    double lon;
    double hgt;
    char extra;
    if (sscanf(location.c_str(), "%lf,%lf,%lf%c", &lat, &lon, &hgt, &extra) != 3) {
        cerr << "ERROR: Invalid location. Expected location=lat,lon,hgt" << endl;
        return false;
    }
    s.opt.llh[0] = lat;
    s.opt.llh[1] = lon;
    s.opt.llh[2] = hgt;

    datetime_t t0;
    if (sscanf(startTime.c_str(), "%d/%d/%d,%d:%d:%lf%c",
               &t0.y, &t0.m, &t0.d, &t0.hh, &t0.mm, &t0.sec, &extra) != 6) {
        cerr << "ERROR: Invalid start_time. Expected YYYY/MM/DD,hh:mm:ss" << endl;
        return false;
    }
    if (t0.y <= 1980 || t0.m < 1 || t0.m > 12 || t0.d < 1 || t0.d > 31 ||
        t0.hh < 0 || t0.hh > 23 || t0.mm < 0 || t0.mm > 59 || t0.sec < 0.0 ||
        t0.sec >= 60.0) {
        cerr << "ERROR: Invalid start_time value." << endl;
        return false;
    }
    t0.sec = floor(t0.sec);
    date2gal(&t0, &s.opt.g0);

    if (rinexFile.size() >= sizeof(s.opt.navfile) || outputFile.size() >= sizeof(s.opt.outfile)) {
        cerr << "ERROR: rinex_file or output_file path is too long." << endl;
        return false;
    }
    strcpy(s.opt.navfile, rinexFile.c_str());
    strcpy(s.opt.outfile, outputFile.c_str());
    PATH = osnmaFile;

    int usrpConfig = 0;
    int bitConfig = 0;
    int durationConfig = 0;
    int neighborConfig = 0;
    if (!parseConfigInt(config, "usrp", usrpConfig) ||
        !parseConfigInt(config, "bit", bitConfig) ||
        !parseConfigInt(config, "duration", durationConfig) ||
        !parseConfigInt(config, "cross_auth_neighbor_count", neighborConfig)) {
        return false;
    }

    if (usrpConfig != 0 && usrpConfig != 1) {
        cerr << "ERROR: Invalid usrp. Expected 0 or 1." << endl;
        return false;
    }
    if (bitConfig != 0 && bitConfig != 1) {
        cerr << "ERROR: Invalid bit. Expected 0 or 1." << endl;
        return false;
    }
    if (durationConfig <= 0) {
        cerr << "ERROR: Invalid duration. It must be a positive integer." << endl;
        return false;
    }
    if (neighborConfig <= 0) {
        cerr << "ERROR: Invalid cross_auth_neighbor_count. It must be a positive integer." << endl;
        return false;
    }

    map<string, string>::const_iterator modeIt = config.find("cross_auth_mode");
    if (modeIt != config.end()) {
        int modeConfig = 1;
        map<string, string> modeOnly;
        modeOnly["cross_auth_mode"] = modeIt->second;
        if (!parseConfigInt(modeOnly, "cross_auth_mode", modeConfig)) {
            return false;
        }
        if (modeConfig < 1 || modeConfig > 3) {
            cerr << "ERROR: Invalid cross_auth_mode. Expected 1, 2, or 3." << endl;
            return false;
        }
        cross_auth_mode = modeConfig;
    }

    if (!parseOptionalConfigPrnList(config, "no_osnma_prns", no_osnma_prns)) {
        return false;
    }
    if (!parseOptionalConfigBool(config, "enable_e5", true, enable_e5)) {
        return false;
    }

    s.opt.use_usrp = (usrpConfig == 0);
    use_usrp = s.opt.use_usrp;
    s.opt.use_bit_stream = (bitConfig != 0);
    duration = static_cast<double>(durationConfig);
    numSub = durationConfig;
    s.opt.iduration = static_cast<int>(duration * 10.0 + 0.5);
    neb_num = neighborConfig;

    return true;
}

string deriveE5OutputPath(const string &e1OutputPath)
{
    size_t slash = e1OutputPath.find_last_of("/\\");
    size_t dot = e1OutputPath.find_last_of('.');

    if (dot == string::npos || (slash != string::npos && dot < slash)) {
        return e1OutputPath + "_e5";
    }

    return e1OutputPath.substr(0, dot) + "_e5" + e1OutputPath.substr(dot);
}

string executableDir(const char *argv0)
{
    if (argv0 == NULL) {
        return ".";
    }

    string path(argv0);
    size_t slash = path.find_last_of("/\\");
    if (slash == string::npos) {
        return ".";
    }
    if (slash == 0) {
        return "/";
    }
    return path.substr(0, slash);
}

bool isExecutableFile(const string &path)
{
    return access(path.c_str(), X_OK) == 0;
}

string findE5GeneratorExecutable(const char *argv0)
{
    vector<string> candidates;
    string dir = executableDir(argv0);
    candidates.push_back(dir + "/osnma-sdr-sim-e5");
    candidates.push_back("./osnma-sdr-sim-e5");
    candidates.push_back("../galileo-sdr-sim/build/usrp_galileo");
    candidates.push_back("galileo-sdr-sim/build/usrp_galileo");

    for (size_t i = 0; i < candidates.size(); i++) {
        if (isExecutableFile(candidates[i])) {
            return candidates[i];
        }
    }

    return "";
}

string formatLocation(const double llh[3])
{
    ostringstream oss;
    oss << setprecision(12) << llh[0] << "," << llh[1] << "," << llh[2];
    return oss.str();
}

pid_t startE5Generation(const char *argv0,
                        const string &configPath,
                        const sim_t &s,
                        double duration)
{
    map<string, string> config;
    if (!readRunConfig(configPath, config)) {
        return -1;
    }

    string startTime;
    if (!requireConfigString(config, "start_time", startTime)) {
        return -1;
    }

    string e5OutputPath;
    map<string, string>::const_iterator outE5 = config.find("output_file_e5");
    if (outE5 != config.end() && !trimConfigValue(outE5->second).empty()) {
        e5OutputPath = trimConfigValue(outE5->second);
    } else {
        e5OutputPath = deriveE5OutputPath(s.opt.outfile);
    }

    string e5Executable = findE5GeneratorExecutable(argv0);
    if (e5Executable.empty()) {
        cerr << "ERROR: Failed to find E5 generator executable. "
             << "Expected osnma-sdr-sim-e5 next to osnma-sdr-sim." << endl;
        return -1;
    }

    int durationSeconds = static_cast<int>(duration + 0.5);
    if (durationSeconds <= 0) {
        cerr << "ERROR: Invalid E5 duration: " << duration << endl;
        return -1;
    }

    vector<string> args;
    args.push_back(e5Executable);
    args.push_back("-e");
    args.push_back(s.opt.navfile);
    args.push_back("-o");
    args.push_back(e5OutputPath);
    args.push_back("-l");
    args.push_back(formatLocation(s.opt.llh));
    args.push_back("-t");
    args.push_back(startTime);
    args.push_back("-d");
    args.push_back(to_string(durationSeconds));
    args.push_back("-U");
    args.push_back("-b");
    args.push_back("-S");
    args.push_back("e5b");

    cout << "\n== E5 Signal Generation ==" << endl;
    cout << "Output file: " << e5OutputPath << endl;

    vector<char *> execArgs;
    for (size_t i = 0; i < args.size(); i++) {
        execArgs.push_back(const_cast<char *>(args[i].c_str()));
    }
    execArgs.push_back(NULL);

    pid_t pid = fork();
    if (pid < 0) {
        perror("ERROR: Failed to fork E5 generator");
        return -1;
    }

    if (pid == 0) {
        execv(e5Executable.c_str(), execArgs.data());
        perror("ERROR: Failed to execute E5 generator");
        _exit(127);
    }

    return pid;
}

int waitE5Generation(pid_t pid)
{
    if (pid <= 0) {
        return -1;
    }

    int status = 0;
    if (waitpid(pid, &status, 0) < 0) {
        perror("ERROR: Failed to wait for E5 generator");
        return -1;
    }

    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
        cerr << "ERROR: E5 generator failed with status " << status << endl;
        return -1;
    }

    return 0;
}


// typedef struct GST gst;

const int osnma_start=138;

const int tagSize=40;

//int getNumberPrn(vector<int> &cur_prn);

//----------------------hex and bin-----------------------------
char binary4ToHex(const string& bin4) {
    static const unordered_map<string, char> binToHex = {
        {"0000", '0'}, {"0001", '1'}, {"0010", '2'}, {"0011", '3'},
        {"0100", '4'}, {"0101", '5'}, {"0110", '6'}, {"0111", '7'},
        {"1000", '8'}, {"1001", '9'}, {"1010", 'A'}, {"1011", 'B'},
        {"1100", 'C'}, {"1101", 'D'}, {"1110", 'E'}, {"1111", 'F'}
    };
    auto it = binToHex.find(bin4);
    return (it != binToHex.end()) ? it->second : '0'; 
}

string hexToBinary(const string& hexStr) {
    const int hexLength = 60;

    string binaryStr;
    binaryStr.reserve(hexStr.size() * 4); 
    
    
    
    for (char c : hexStr) {
        int value;
        if (isdigit(c)) {
            value = c - '0';
           
        } else if (isxdigit(c)) {
           
            value = 10 + (toupper(c) - 'A');
           
        } else {
       
            throw invalid_argument("Invalid hexadecimal character: " + string(1, c));
        }

  
        bitset<4> bits(value);
        binaryStr += bits.to_string();
    }
    return binaryStr;
}
string binaryToHexBitset(const string& binaryStr) {
    string hexStr;
    string paddedBinary = binaryStr;
    
 
    int padding = (4 - paddedBinary.length() % 4) % 4;
   // cout<<"pad :"<<padding<<endl;
    paddedBinary = string(padding, '0') + paddedBinary;
    
  
    for (size_t i = 0; i < paddedBinary.length(); i += 4) {
        string group = paddedBinary.substr(i, 4);
        hexStr += binary4ToHex(group);
    }
    return hexStr;
}
//----------------------hex and bin-----------------------------

//----------------------hash 256--------------------------------
#define ROTRIGHT(a,b) (((a) >> (b)) | ((a) << (32-(b))))

const uint32_t k[64] = {
  0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,
  0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
  0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,
  0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
  0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,
  0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
  0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,
  0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
  0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,
  0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
  0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,
  0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
  0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,
  0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
  0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,
  0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
};

class SHA256 {
public:
    SHA256() {
        reset();
    }

    void update(const uint8_t *data, size_t len) {
        for (size_t i = 0; i < len; ++i) {
            dataBuffer[dataLen++] = data[i];
            if (dataLen == 64) {
                transform(dataBuffer);
                bitLen += 512;
                dataLen = 0;
            }
        }
    }

    void final(uint8_t *hash) {
        size_t i = dataLen;

        // Pad whatever data is left
        if (dataLen < 56) {
            dataBuffer[i++] = 0x80;
            while (i < 56)
                dataBuffer[i++] = 0x00;
        } else {
            dataBuffer[i++] = 0x80;
            while (i < 64)
                dataBuffer[i++] = 0x00;
            transform(dataBuffer);
            memset(dataBuffer, 0, 56);
        }

        bitLen += dataLen * 8;
        dataBuffer[63] = bitLen;
        dataBuffer[62] = bitLen >> 8;
        dataBuffer[61] = bitLen >> 16;
        dataBuffer[60] = bitLen >> 24;
        dataBuffer[59] = bitLen >> 32;
        dataBuffer[58] = bitLen >> 40;
        dataBuffer[57] = bitLen >> 48;
        dataBuffer[56] = bitLen >> 56;
        transform(dataBuffer);

        for (i = 0; i < 4; ++i) {
            for (int j = 0; j < 8; ++j) {
                hash[i + j * 4] = (state[j] >> (24 - i * 8)) & 0x000000ff;
            }
        }
    }

private:
    uint8_t dataBuffer[64];
    uint32_t dataLen = 0;
    uint64_t bitLen = 0;
    uint32_t state[8];

    void reset() {
        dataLen = 0;
        bitLen = 0;
        state[0] = 0x6a09e667;
        state[1] = 0xbb67ae85;
        state[2] = 0x3c6ef372;
        state[3] = 0xa54ff53a;
        state[4] = 0x510e527f;
        state[5] = 0x9b05688c;
        state[6] = 0x1f83d9ab;
        state[7] = 0x5be0cd19;
    }

    void transform(const uint8_t *data) {
        uint32_t m[64];
        uint32_t a, b, c, d, e, f, g, h;

        for (uint32_t i = 0, j = 0; i < 16; ++i, j += 4)
            m[i] = (data[j] << 24) | (data[j + 1] << 16) |
                   (data[j + 2] << 8) | (data[j + 3]);
        for (uint32_t i = 16; i < 64; ++i) {
            uint32_t s0 = ROTRIGHT(m[i - 15], 7) ^ ROTRIGHT(m[i - 15], 18) ^ (m[i - 15] >> 3);
            uint32_t s1 = ROTRIGHT(m[i - 2], 17) ^ ROTRIGHT(m[i - 2], 19) ^ (m[i - 2] >> 10);
            m[i] = m[i - 16] + s0 + m[i - 7] + s1;
        }

        a = state[0];
        b = state[1];
        c = state[2];
        d = state[3];
        e = state[4];
        f = state[5];
        g = state[6];
        h = state[7];

        for (uint32_t i = 0; i < 64; ++i) {
            uint32_t S1 = ROTRIGHT(e, 6) ^ ROTRIGHT(e, 11) ^ ROTRIGHT(e, 25);
            uint32_t ch = (e & f) ^ ((~e) & g);
            uint32_t temp1 = h + S1 + ch + k[i] + m[i];
            uint32_t S0 = ROTRIGHT(a, 2) ^ ROTRIGHT(a, 13) ^ ROTRIGHT(a, 22);
            uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
            uint32_t temp2 = S0 + maj;

            h = g;
            g = f;
            f = e;
            e = d + temp1;
            d = c;
            c = b;
            b = a;
            a = temp1 + temp2;
        }

        state[0] += a;
        state[1] += b;
        state[2] += c;
        state[3] += d;
        state[4] += e;
        state[5] += f;
        state[6] += g;
        state[7] += h;
    }
};

// -------------------- HMAC-SHA256 --------------------




//-------------------------hash  256-----------------------------

void getNavMsg()
{
    string line;
    ifstream inf;
    inf.open("E://Vs_project//test.csv");
    
    while(getline(inf,line))
    {   
        std::vector<std::string> columns;
        std::stringstream ss(line);
        std::string item;
        while(std::getline(ss, item, ','))
        {
            columns.push_back(item);
        }
        int temp_num=stoi(columns[0]);
        page_struct a;
        a.tow=stod(columns[0]);
        a.wn=cur_wn;
        a.prn=stoi(columns[2]);
        a.navMsg=binaryToHexBitset(columns[3]);
        nav_page.push_back(a);
        if(a.prn==2)
        {
           // cout<<a.tow<<"---------------"<<endl;
            nav_page_2.push_back(a);
        }
    }
}
void show(vector<page_struct> navPage)
{
    for(int i=0;i<navPage.size();i++)
    {
        cout<<navPage[i].tow<<","<<cur_wn<<","<<navPage[i].prn<<","<<binaryToHexBitset(navPage[i].navMsg)<<endl;
    }
}

void saveNavPageTxt(const vector<page_struct> &navPage, const string &outputPath)
{
    ofstream out(outputPath);
    if(!out.is_open())
    {
        cerr<<"[WARN] Failed to write navigation data to "<<outputPath<<endl;
        return;
    }

    for(size_t i=0;i<navPage.size();i++)
    {
        out<<navPage[i].tow<<","<<cur_wn<<","<<navPage[i].prn<<","<<binaryToHexBitset(navPage[i].navMsg)<<endl;
    }
}

static string joinPrnList(const vector<int> &prns)
{
    ostringstream oss;
    for (size_t i = 0; i < prns.size(); i++)
    {
        if (i > 0)
            oss << ", ";
        oss << prns[i];
    }
    return oss.str();
}

static size_t navMsgBitLength(const string &navMsg)
{
    if (navMsg.empty())
        return 0;
    bool binary = true;
    for (size_t i = 0; i < navMsg.size(); i++)
    {
        if (navMsg[i] != '0' && navMsg[i] != '1')
        {
            binary = false;
            break;
        }
    }
    return binary ? navMsg.size() : hexToBinary(navMsg).size();
}

static map<int, int> countNavPagesByPrn(const vector<page_struct> &navPage)
{
    map<int, int> counts;
    for (size_t i = 0; i < navPage.size(); i++)
        counts[navPage[i].prn]++;
    return counts;
}

static bool validateNavPageCoverage(const vector<page_struct> &navPage,
                                    const vector<int> &requiredPrns)
{
    bool ok = true;
    map<int, int> counts;
    map<int, int> firstTow;
    map<int, int> lastTow;

    for (size_t i = 0; i < navPage.size(); i++)
    {
        int prn = navPage[i].prn;
        int tow = static_cast<int>(navPage[i].tow);
        counts[prn]++;
        if (firstTow.find(prn) == firstTow.end() || tow < firstTow[prn])
            firstTow[prn] = tow;
        if (lastTow.find(prn) == lastTow.end() || tow > lastTow[prn])
            lastTow[prn] = tow;

        size_t bitLength = navMsgBitLength(navPage[i].navMsg);
        if (bitLength < 240)
        {
            cerr << "[ERROR] Missing or short nav page: PRN " << prn
                 << " TOW " << tow << " bit_length=" << bitLength << endl;
            ok = false;
        }
    }

    for (size_t i = 0; i < requiredPrns.size(); i++)
    {
        int prn = requiredPrns[i];
        if (counts[prn] == 0)
        {
            cerr << "[ERROR] all_visible_prn_list contains PRN " << prn
                 << " but nav_page count is 0" << endl;
            ok = false;
        }
        else
        {
            cout << "[COVERAGE] PRN " << prn << ": " << counts[prn]
                 << " pages, tow " << firstTow[prn] << ".." << lastTow[prn] << endl;
        }
    }

    return ok;
}

static bool isPageAlignedToSimulationGrid(int tow)
{
    int startTow = static_cast<int>(start_time);
    return tow >= startTow && ((tow - startTow) % 2) == 0;
}

static bool isPrnVisibleAtTow(int prn, int tow)
{
    if (tow_visible_prn_sets.empty())
        return true;

    for (size_t i = 0; i < tow_visible_prn_sets.size(); i++)
    {
        int rangeStart = tow_visible_prn_sets[i].first;
        int rangeEnd = (i + 1 < tow_visible_prn_sets.size())
                           ? tow_visible_prn_sets[i + 1].first
                           : INT_MAX;
        if (tow >= rangeStart && tow < rangeEnd)
        {
            const vector<int> &visiblePrns = tow_visible_prn_sets[i].second;
            return find(visiblePrns.begin(), visiblePrns.end(), prn) != visiblePrns.end();
        }
    }

    return false;
}

static void writeVisibilityReport(const string &outputPath,
                                  double startTow,
                                  double duration,
                                  const vector<page_struct> &navPage)
{
    ofstream out(outputPath);
    if (!out.is_open())
    {
        cerr << "[WARN] Failed to write visibility report to " << outputPath << endl;
        return;
    }

    map<int, int> counts = countNavPagesByPrn(navPage);
    map<int, int> navFirstTow;
    map<int, int> navLastTow;
    for (size_t i = 0; i < navPage.size(); i++)
    {
        int prn = navPage[i].prn;
        int tow = static_cast<int>(navPage[i].tow);
        if (navFirstTow.find(prn) == navFirstTow.end() || tow < navFirstTow[prn])
            navFirstTow[prn] = tow;
        if (navLastTow.find(prn) == navLastTow.end() || tow > navLastTow[prn])
            navLastTow[prn] = tow;
    }

    out << "Start TOW: " << static_cast<int>(startTow) << "\n";
    out << "Duration: " << static_cast<int>(duration) << "\n";
    out << "Scan step: " << visibility_scan_step_seconds << " seconds\n\n";

    out << "Initial visible PRNs:\n" << joinPrnList(initial_visible_prn_list) << "\n\n";
    out << "All visible PRNs during simulation:\n" << joinPrnList(all_visible_prn_list) << "\n\n";
    out << "Future added PRNs:\n" << joinPrnList(future_added_prn_list) << "\n\n";

    out << "PRN first/last seen:\n";
    for (size_t i = 0; i < all_visible_prn_list.size(); i++)
    {
        int prn = all_visible_prn_list[i];
        out << "PRN " << prn << ": first=" << prn_first_seen_tow[prn]
            << ", last=" << prn_last_seen_tow[prn] << "\n";
    }

    out << "\nTOW ranges:\n";
    for (size_t i = 0; i < tow_visible_prn_sets.size();)
    {
        size_t j = i + 1;
        while (j < tow_visible_prn_sets.size() &&
               tow_visible_prn_sets[j].second == tow_visible_prn_sets[i].second)
        {
            j++;
        }
        int rangeStart = tow_visible_prn_sets[i].first;
        int rangeEnd = (j < tow_visible_prn_sets.size())
                           ? tow_visible_prn_sets[j].first
                           : static_cast<int>(startTow + duration);
        out << "[" << rangeStart << ", " << rangeEnd << "): "
            << joinPrnList(tow_visible_prn_sets[i].second) << "\n";
        i = j;
    }

    out << "\nNav page count:\n";
    for (size_t i = 0; i < all_visible_prn_list.size(); i++)
    {
        int prn = all_visible_prn_list[i];
        out << "PRN " << prn << ": " << counts[prn] << " pages";
        if (counts[prn] > 0)
            out << " (tow " << navFirstTow[prn] << ".." << navLastTow[prn] << ")";
        out << "\n";
    }
}

void getWeekTow()
{
    string line1;
    ifstream inf;
    //cout<<"asdasd"<<PATH<<endl;
    inf.open(PATH); //Read the original OSNMA file

    vector<page_struct> osnma_r;
    
    while(getline(inf,line1))
    {   
   
        std::vector<std::string> columns;
        std::stringstream ss(line1);
        std::string item;
        while(std::getline(ss, item, ','))
        {
            columns.push_back(item);
        }
        if (columns.size() < 2)
        {
            continue;
        }
        if (columns.size() == 2)
        {
            inf.close();
            return;
        }
        
            page_struct a;

            cur_wn = stoi(columns[1]);
            
            //cout<<"-----++---------"<<cur_wn<<start_time<<endl;
            inf.close();
	     return;
            
           // cout<<columns[3]<<endl;
        
    }
    
}
// ------------------------------ADKD4------------------------

// ------------------------------ADKD12------------------------

void show2(vector<page_struct> nav_ss)
{
    for(int i=0;i<nav_ss.size();i++)
    {
        cout<<nav_ss[i].tow<<","<<cur_wn<<","<<nav_ss[i].prn<<","<<nav_ss[i].navMsg<<endl;
    }
    cout<<"----------------2----------------"<<endl;
}

vector<KeyItem> getPageKeyItems(vector<page_struct> &cur_page)
{
    vector<KeyItem> keyItems;
    string mack;
    for(size_t i = 0; i < cur_page.size(); i++)
    {
        int tow = static_cast<int>(cur_page[i].tow);
        if(getPageNumInSub(tow) >= 10 && getPageNumInSub(tow) <= 14)
        {
            mack.append(getMackInPage(cur_page[i].navMsg));
        }
        if(getPageNumInSub(tow) == 14)
        {
            if(mack.size() >= 36)
            {
                KeyItem item;
                item.key = mack.substr(4, 32);
                item.firstTow = getFirstTow(tow);
                keyItems.push_back(item);
            }
            mack.clear();
        }
    }
    return keyItems;
}

vector<string> alignKeysToGst(const vector<KeyItem> &keyItems, const vector<gst> &gstList, size_t authCount, int keyDelayTow)
{
    vector<string> alignedKeys;
    for(size_t i = 0; i < authCount && i < gstList.size(); i++)
    {
        int targetFirstTow = stoi(gstList[i].tow) + keyDelayTow;
        string key;
        for(size_t j = 0; j < keyItems.size(); j++)
        {
            if(keyItems[j].firstTow == targetFirstTow)
            {
                key = keyItems[j].key;
                break;
            }
        }
        alignedKeys.push_back(key);
    }
    return alignedKeys;
}

vector<string> alignKeysToAuthDataGst(const vector<KeyItem> &keyItems, const vector<string> &authData, int keyDelayTow)
{
    vector<string> alignedKeys;
    for(size_t i = 0; i < authData.size(); i++)
    {
        string key;
        if(authData[i].size() >= 48)
        {
            int tow = static_cast<int>(strtol(authData[i].substr(28, 20).c_str(), NULL, 2));
            int targetFirstTow = tow + keyDelayTow;
            for(size_t j = 0; j < keyItems.size(); j++)
            {
                if(keyItems[j].firstTow == targetFirstTow)
                {
                    key = keyItems[j].key;
                    break;
                }
            }
        }
        alignedKeys.push_back(key);
    }
    return alignedKeys;
}

vector<string> expandAdkd4TagsToFillSlots(vector<page_struct> &cur_page, const vector<string> &authData, const vector<string> &tags)
{
    map<int, string> tagByTow;
    for(size_t i = 0; i < authData.size() && i < tags.size(); i++)
    {
        if(authData[i].size() >= 48)
        {
            int tow = static_cast<int>(strtol(authData[i].substr(28, 20).c_str(), NULL, 2));
            tagByTow[tow] = tags[i];
        }
    }

    vector<string> expanded;
    for(size_t i = 0; i < cur_page.size(); i++)
    {
        int tow = static_cast<int>(cur_page[i].tow);
        int firstTow = getFirstTow(tow);
        if(getPageNumInSub(tow) == 3 && is60page(firstTow) && firstTow > static_cast<int>(start_time))
        {
            map<int, string>::const_iterator it = tagByTow.find(firstTow - 1);
            expanded.push_back(it == tagByTow.end() ? "0000000000" : it->second);
        }
    }
    return expanded;
}

void fillmecseq(vector<string> macseq, vector<page_struct> &nav_page) {
    if (cross_auth_mode == 3) {
        return;
    }
    for (size_t i = 0; i < nav_page.size(); i++) {
        int tow = static_cast<int>(nav_page[i].tow);
        if (getPageNumInSub(tow) == 0 && tow > start_time) // Found the start of a subframe.
        {
            i++;
            string navBin= hexToBinary(nav_page[i].navMsg);
            string macBin= hexToBinary(macseq[0]);

            //cout << "navMsg before filling: " << nav_page[i].navMsg << endl;

            //cout<< "macBin to fill: " << macseq[0] << endl;

            int start = 154; // Start position, equivalent to 146 + 8 in Python.
            navBin.replace(start, macBin.size(), macBin);
            string navHex= binaryToHexBitset(navBin);
		    nav_page[i].navMsg=navHex;
	            
            macseq.erase(macseq.begin());
            if(macseq.empty()) return;
        }
        
    }
}
void printSequence(const vector<vector<string>>& seq) {
	
    cout<<"MACLT sequence:"<<endl;
    if (seq.empty()) {
        cout << "[WARN] MACLT sequence not found for the selected ID." << endl;
        return;
    }

    for (const auto& row : seq) {
        for (const auto& s : row) {
            cout << s << " ";
        }
        cout << endl;
    }
}



int main(int argc, char *argv[])
{
    if (argc != 3 || strcmp(argv[1], "-f") != 0)
    {
        usage(argv[0]);
        exit(1);
    }

    // Set default values
    sim_t s;

    s.finished = false;
    s.opt.navfile[0] = 0;
    s.opt.tvfile[0] = 0;
    s.opt.outfile[0] = 0;

    s.opt.g0.week = -1;
    s.opt.g0.sec = 0.0;
    s.opt.iduration = USER_MOTION_SIZE;
    s.opt.verb = FALSE;
    s.opt.llh[0] = 42.3601;
    s.opt.llh[1] = -71.0589;
    s.opt.llh[2] = 2;
    s.opt.interactive = FALSE;
    s.opt.timeoverwrite = FALSE;
    s.opt.iono_enable = TRUE;
    s.udp_port = 5671;
    s.opt.use_usrp = true;
    s.opt.use_bit_stream = true;

    // Options
    double duration = s.opt.iduration / 10.0;
    double gain = 0;
    std::string device_args = "";
    size_t channel = 1;
    bool verbose = false;
    int return_code = EXIT_SUCCESS;
    char error_string[512];
    double rate = TX_SAMPLERATE;
    double freq = TX_FREQUENCY;
    bool use_usrp = true;
    bool enable_e5 = true;
    // Buffer sizes
    size_t samps_per_buff = SAMPLES_PER_BUFFER;
    float *buffer = NULL;
    const void **buffer_ptr = NULL;

    if (!applyRunConfig(argv[2], s, duration, use_usrp, enable_e5))
    {
        exit(1);
    }

    if (s.opt.navfile[0] == 0 && s.opt.tvfile[0] == 0)
    {
        printf("ERROR: Galileo ephemeris/nav_msg file is not specified.\n");
        exit(1);
    }

    if (s.opt.outfile[0] == 0)
    {
        printf("[CONFIG] output_file not set; using galileosim.bin\n");
        snprintf(s.opt.outfile, sizeof(s.opt.outfile), "galileosim.ishort");
    }
	printf("\n== Galileo E1/OSNMA Generation ==\n");
    // Initialize simulator
    numSub=numSub/30;
    	//cout<<"numSub:"<<numSub<<endl;
    init_sim(&s);
    sim_t s2=s;
    cur_wn = s.opt.g0.week - 1024;
    //cout<<"--------------"<<cur_wn<<endl;
    getWeekTow();
    galileo_task1(&s2);

    //getNavMsg(); //getpage
   //cout<<PATH<<endl;
    //cout<<"----asdasdas-------"<<cur_wn<<endl;
    vector<int> cur_prn;

     for(int i=0;i<initial_visible_prn_list.size();i++)
    {
       //cout<<cur_prn[i]<<endl;
       nav_page.erase(nav_page.begin());
    }
    nav_page_all = nav_page;

	// Keep the generation input on the simulation TOW grid. Do not trim by
	// visibility here, because OSNMA tag/key indexes are subframe-relative.
	vector<int> cross_auth_prn_list = all_visible_prn_list;
	for (size_t i = 0; i < nearest_prn_by_visible.size(); i++)
	{
	    for (size_t j = 0; j < nearest_prn_by_visible[i].size(); j++)
	    {
		int neighborPrn = nearest_prn_by_visible[i][j];
		if (!containsPrn(cross_auth_prn_list, neighborPrn))
		{
		    cross_auth_prn_list.push_back(neighborPrn);
		}
	    }
	}
	cout << "[CrossAuth] Prepared " << cross_auth_prn_list.size()
	     << " PRN sources for authentication data." << endl;

	vector<page_struct> filtered_nav_page;

		for (size_t i = 0; i < nav_page_all.size(); i++)
		{
		    int pagePrn = nav_page_all[i].prn;
		    int pageTow = static_cast<int>(nav_page_all[i].tow);
		    if (!isPageAlignedToSimulationGrid(pageTow))
		    {
			continue;
		    }
		    for (size_t j = 0; j < cross_auth_prn_list.size(); j++)
		    {
			if (pagePrn == cross_auth_prn_list[j])
			{
			    filtered_nav_page.push_back(nav_page_all[i]);
			    break;
			}
		    }
		}

	// nav_page_all is the read-only source for neighbor auth data. It must
	// match the aligned generation grid, but final output is filtered later.
    nav_page_all = filtered_nav_page;
    nav_page = filtered_nav_page;
    //show2(nav_page);
    //show2(nav_page_all);
    cur_prn = all_visible_prn_list;
    int prnNum=static_cast<int>(cur_prn.size());
    
    
    size_t pos = PATH.find_last_of('/');

    // 2. Extract the filename.
    string filename = PATH.substr(pos + 1);

    // 3. Keep the first two characters, for example "07".
    if (filename.size() >= 2) {
        filename =  filename.substr(0, 2);
    }
    //cout <<"filename :" <<filename << endl;
   int prrn;
   
   if(filename == "07"){
   	prrn=8;
   }
   else if(filename == "16"){
   	prrn=2;
   }
   else{
   	prrn=30;
   }
    
    
    vector<KeyItem> key= getKeyItem(prrn,PATH);
    
   for (size_t i = 0; i < key.size(); i++)
   {
    /* code */
      //cout<<"key: "<<key[i].key<<"  Tow: "<<key[i].firstTow<<endl;
   }
    

    
    cout << "Selected PRNs: [";
	for (int i = 0; i < prnNum; i++)
	{
	    cout << cur_prn[i];
	    if (i != prnNum - 1)
		cout << ",";
	}
	cout << "]" << endl;
	    promptFocusPrns(cur_prn);
	    if (!validateNoOsnmaPrnsAgainstSelected(cur_prn))
	    {
	        exit(1);
	    }
	    printCrossAuthModeSummary(cur_prn);
	            //cout<<"---------------------------------------------"<<endl;
	    
	        cout<<"\n== OSNMA Configuration =="<<endl;
    
    //cout<<nav_page.size()<<endl;

    int fl=0;

    

    for(int i=0;i<prnNum;i++)
    {
    	
    	int prn = cur_prn[i];
    	
    	//cout<<"Generate navigation data for prn "<<prn<<endl;
	    	vector<page_struct> cur_page;
	    
	    	fixNavPage(prn,cur_page);
	    	if(cur_page.empty())
	    	{
	    	    cout<<"[WARN] E"<<prn<<": no visible navigation pages; skipped."<<endl;
	    	    continue;
	    	}
	    	if(cross_auth_mode != 1 && !isOsnmaEnabledPrn(prn))
	    	{
	    	    cout<<"[E"<<prn<<"] OSNMA disabled by cross_auth_mode "<<cross_auth_mode
	    	        <<"; generating navigation data only."<<endl;
	    	    addWord10(cur_page);
	    	    fillCrc(cur_page);
	    	    singleIn(prn,cur_page);
	    	    cout<<"[E"<<prn<<"] Done."<<endl;
	    	    continue;
	    	}
		
		    	fillOSNMA(prrn,cur_page,PATH);
	    	if(fl==0){
		getDic(cur_page);
		if(dic!=33 && dic!=34 && dic!=35 && dic!=29)
		{
		    dic=34;
		    HF=0;
		    MF=0;
		    KS=4;
		    TS=9;
		}
		fl=1;
		
		findSeqByID(dic);
	printSequence(mltSequence);
	KL=ks_num[KS];
	TL=ts_num[TS-5];
	TL=TL;
	cout<<"\nMessage authentication:"<<endl;
	cout<<"Hash function: SHA-256"<<endl;
	cout<<"Key size:      "<<KL<<" bits"<<endl;
	cout<<"Tag size:      "<<TL<<" bits"<<endl;
	//cout<<"HF: "<<HF<<"MF: "<<MF<<"dic: "<<dic<<" KL:"<<KL<<" TL:"<<TS<<endl;
		    cout<<endl;
	    cout<<"== Galileo I/NAV Message Generation =="<<endl;
	    cout<<"Generating navigation pages and OSNMA fields..."<<endl;

	}
	
	cleanTag(cur_page);
	fillTagInfo(prn,cur_page);
    	addWord10(cur_page);
    	cout<<"[E"<<prn<<"] Loading HKROOT/TESLA data and generating tags ";
   	if(duration>361)
   	{
   		cout<<"(ADKD 0, 4, 12)..."<<endl;
   	}
   	else{
   	cout<<"(ADKD 0, 4)..."<<endl;
   	
   	}
	   	
	        vector<string> teslaKey=getKey(cur_page);
	        vector<KeyItem> pageKeyItems = getPageKeyItems(cur_page);
	        //----------------------mac-------------------------
        vector<string> teslaKeyformac=teslaKey;
        
	        if(!teslaKeyformac.empty())
	        {
	            teslaKeyformac.erase(teslaKeyformac.begin());
	        }
        
        vector<string> macauthData=getmacauth(prn,cur_page);
        
        vector<string> macList=commac(teslaKeyformac,macauthData);

		if(!macList.empty())
		{
		    macList.erase(macList.begin());
		}
		fillmecseq(macList,cur_page);
		
	        //----------------------mac-------------------------
	        if(!teslaKey.empty())
	        {
	            teslaKey.erase(teslaKey.begin());
	        }
	        if(!teslaKey.empty())
	        {
	            teslaKey.erase(teslaKey.begin());
	        }

        vector<string> authDataADKD4 = getAuthDataForADKD4(prn,cur_page);
        vector<string> teslaKeyforADKD4=teslaKey;
        
	        vector<string> KeyforADKD4=alignKeysToAuthDataGst(pageKeyItems, authDataADKD4, 31);
        
	
		if(!teslaKeyforADKD4.empty())
		{
			teslaKeyforADKD4.erase(teslaKeyforADKD4.begin());
		}
	
	
	//cout<<teslaKeyforADKD4.size()<<"++++++"<<KeyforADKD4.size()<<endl;
	for(int i = 0; i < teslaKeyforADKD4.size(); i++) {
	  // cout<<"1:"<<teslaKeyforADKD4[i]<<endl;
	   	   //cout<<"2:"<<KeyforADKD4[i]<<endl;
	}
	
	
	for(int i = static_cast<int>(teslaKeyforADKD4.size()) - 1; i >= 0; i--) {
	    if(i % 2 == 1) { // Odd index.
		teslaKeyforADKD4.erase(teslaKeyforADKD4.begin() + i);
	    }
	}
	size_t adkd4CommonSize = min(authDataADKD4.size(), min(teslaKeyforADKD4.size(), KeyforADKD4.size()));
	authDataADKD4.resize(adkd4CommonSize);
	teslaKeyforADKD4.resize(adkd4CommonSize);
	KeyforADKD4.resize(adkd4CommonSize);
	//cout<<"authchangdu:"<<authDataADKD4.size()<<","<<KeyforADKD4.size()<<endl;
	
	vector<string> tagListADKD4=computeTagADKD4(KeyforADKD4,authDataADKD4);
	tagListADKD4 = expandAdkd4TagsToFillSlots(cur_page, authDataADKD4, tagListADKD4);
	if(dic==33)
	{
		tagListADKD4.insert(tagListADKD4.begin(), "bfab8f967d");
	}
	//-----------------------ADKD12-----------------------
	if(duration>361)
	{
			vector<string> authDataADKD12=getAuthDatADKD12(cur_page,prn);
			vector<gst> adkd12Gst = getGstADKD12(cur_page);
			vector<string> teslaKeyforADKD12 = alignKeysToGst(pageKeyItems, adkd12Gst, authDataADKD12.size(), 331);
		if(!teslaKeyforADKD12.empty())
		{
			vector<string> tagListADKD12=computeTagADKD12(teslaKeyforADKD12,authDataADKD12);
			vector<string> authDataADKD12Counter6=getAuthDataADKD12ForCounter(cur_page,prn,6);
			vector<string> tagListADKD12Counter6=computeTagADKD12(teslaKeyforADKD12,authDataADKD12Counter6);
			
			for(int i = tagListADKD12.size() - 1; i >= 0; i--) {
			    if(i % 2 == 1) { // Odd index.
				//tagListADKD12.erase(tagListADKD12.begin() + i);
			    }
			}
			if(!tagListADKD12.empty())
			{
				fillTagADKD12(cur_page,tagListADKD12);
			}
			if(!tagListADKD12Counter6.empty())
			{
				fillTagADKD12Counter(cur_page,tagListADKD12Counter6,6);
			}
		}
		else
		{
			cout<<"[WARN] E"<<prn<<": ADKD12 skipped because the delayed TESLA key is missing."<<endl;
		}
	}
	
	//-----------------------ADKD12-----------------------


	        vector<string> authData=getAuthData(prn ,cur_page);
		
		        vector<gst> adkd0Gst = getGst(cur_page);
	        vector<string> teslaKeyForADKD0 = alignKeysToGst(pageKeyItems, adkd0Gst, authData.size(), 31);
	        vector<string> tagList=computeTag(teslaKeyForADKD0,authData);
	        
		fillTagADKD4(cur_page,tagListADKD4);
	
	 
	        fillTag(tagList,cur_page);
	        fillTagADKD0ExtraSlots(teslaKeyForADKD0, authData, cur_page);
	        //cout<<"Test:"<<tagList.size()<<endl;
        //show2(cur_page);
        fillCrc(cur_page);
        
cout<<"[E"<<prn<<"] Done."<<endl;
        singleIn(prn,cur_page);
        
    }
    
	    // Cross-authentication.
	    for(int i=0;i<prnNum;i++) // Iterate over visible satellites.
		{
		    int targetPrn = cur_prn[i];
		    if(!canBeCrossAuthTargetPrna(targetPrn, cross_auth_mode))
		    {
		        // Debug-only detail: target PRN has no OSNMA field in this mode.
		        continue;
		    }
		    runCrossAuthenticationForVisiblePrn(i, targetPrn);
		}

    
    

    cout<<"[I/NAV] Filtering generated pages to visible PRNs."<<endl;
    for(int i=0;i<prnNum;i++)
    {
    	//cout<<"Generate osnma data for prn "<<cur_prn[i]<<endl;
    }
        vector<page_struct> visible_nav_page;
        for (size_t i = 0; i < nav_page.size(); i++)
        {
            int pagePrn = nav_page[i].prn;
            int pageTow = static_cast<int>(nav_page[i].tow);
            if (isPrnVisibleAtTow(pagePrn, pageTow))
            {
                visible_nav_page.push_back(nav_page[i]);
            }
        }
        nav_page = visible_nav_page;

	    bool coverage_ok = validateNavPageCoverage(nav_page, all_visible_prn_list);
	    writeVisibilityReport("../tow_satellite_visibility_report.txt",
	                          s.opt.g0.sec,
	                          duration,
	                          nav_page);
	    if (!coverage_ok)
	    {
	        cerr << "[ERROR] Navigation page coverage check failed before signal generation." << endl;
	        exit(1);
	    }

	    saveNavPageTxt(nav_page, "../nav_page.txt");
	 //   show(nav_page);
	    writeGeneratedOsnmaTagLog("../osnma_tag_log.csv");
	    writeGeneratedOsnmaTagStatisticsTxt("../osnma_tag_statistics.txt");
	
	    flag_write_nav=1;
    //exit(1);
    

    // Allocate FIFOs to hold 0.1 seconds of I/Q samples each.
    s.fifo = (short *)malloc(FIFO_LENGTH * sizeof(short) * 2); // for 32-bit I and Q samples

    if (s.fifo == NULL)
    {
        printf("ERROR: Failed to allocate I/Q sample buffer.\n");
        // goto out;
    }

    if (use_usrp)
    {   
        usrp_conf_t usrp_conf;  
        usrp_conf.carr_freq = GALILEO_E1_FREQ_HZ;
        usrp_conf.gain = gain;
        usrp_conf.device_args = device_args;
        usrp_conf.samp_rate = TX_SAMPLERATE;

        init_usrp(usrp_conf, &s);
    }

    // Start Galileo task.
    pid_t e5_pid = -1;
    if (enable_e5)
    {
        e5_pid = startE5Generation(argv[0], argv[2], s, duration);
        if (e5_pid <= 0)
        {
            fprintf(stderr, "\n[ERROR] E5 generation failed to start.\n");
            exit(1);
        }
    }
    else
    {
        cout << "\n== E5 Signal Generation ==" << endl;
        cout << "Status: disabled by enable_e5=false" << endl;
    }

    s.status = start_galileo_task(&s);
    if (s.status < 0)
    {
        fprintf(stderr, "[ERROR] Failed to start Galileo task.\n");
    }
    else
        //printf("\nCreating Galileo task...\n");

    if (use_usrp)
    { // Wait until Galileo task is initialized
        pthread_mutex_lock(&(s.tx.lock));
        while (!s.galileo_sim.ready)
            pthread_cond_wait(&(s.galileo_sim.initialization_done), &(s.tx.lock));

        pthread_mutex_unlock(&(s.tx.lock));

        // Fillfull the FIFO.
        if (is_fifo_write_ready(&s))
            pthread_cond_signal(&(s.fifo_write_ready));

        // // Start TX task
        s.status = start_tx_task(&s);

        if (s.status < 0)
        {
            fprintf(stderr, "[ERROR] Failed to start TX task.\n");
            // goto out;
        }
        else
            printf("[TX] Transmission task created.\n");

        // Running...
        printf("[RUN] Signal generation is running.\n"
               "Press Ctrl+C to abort.\n");

        // Wainting for TX task to complete.
        pthread_join(s.tx.thread, NULL);
    }
    else
        pthread_join(s.galileo_sim.thread, NULL);

    fprintf(stderr, "\n[TX] Total samples consumed: %ld", samples_consumed);

    if (enable_e5 && waitE5Generation(e5_pid) != 0)
    {
        fprintf(stderr, "\n[ERROR] E5 generation failed.\n");
        exit(1);
    }

    printf("\n[DONE] Signal generation completed.\n");
}
