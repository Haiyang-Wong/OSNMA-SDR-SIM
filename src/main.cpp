#define _CRT_SECURE_NO_WARNINGS

#include "../include/galileo-sdr.h"
#include <math.h>
#include <unistd.h>
#include <vector>
#include<malloc.h>
#include<iostream>
#include<string>
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
int subNub=10;
const string mapping = "0123456789ABCDEF";
double start_time = 1000000000.0;

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
    fprintf(stderr, "Done\n");
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
    fprintf(stderr, "\nTx task");
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
            fprintf(stderr, "\rTime = %4.1f - %d", s->time, k);
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
"Usage: %s [options]\n"
        "Options:\n"
        "  -e <Ephemeris>   RINEX navigation file for Galileo ephemerides (required)\n"
        "  -o <File sink>   File to store IQ samples\n"
        "  -l <location>    Lat,Lon,Hgt (static mode) e.g. 35.274,137.014,100\n"
        "  -t <date,time>   Scenario start time YYYY/MM/DD,hh:mm:ss\n"
        "  -T <date,time>   Overwrite TOC and TOE to scenario start time\n"
        "  -d <duration>    Duration [sec] (max: %.0f)\n"
        "  -a <rf_gain>     Absolute RF gain in [0 ... 60] (default: 30)\n"
        "  -s <PATH>        Custom path for config/resource/anything you want\n"
        "  -U               Disable USRP (-U 1)\n"
        "  -b               Disable Bit stream (-b 1)\n",
        progname,

        ((double)USER_MOTION_SIZE) / 10.0);

    return;
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
       
            throw invalid_argument("无效的十六进制字符: " + string(1, c));
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

void fillmecseq(vector<string> macseq, vector<page_struct> &nav_page) {
    for (size_t i = 0; i < nav_page.size(); i++) {
        int tow = static_cast<int>(nav_page[i].tow);
        if (getPageNumInSub(tow) == 0 && tow > start_time) // 成功，这是一个开头
        {
            i++;
            string navBin= hexToBinary(nav_page[i].navMsg);
            string macBin= hexToBinary(macseq[0]);

            //cout << "填充前的navMsg: " << nav_page[i].navMsg << endl;

            //cout<< "填充的macBin: " << macseq[0] << endl;

            int start = 154; // 你的开始位置，比如 Python 里的 146+8
            navBin.replace(start, macBin.size(), macBin);
            string navHex= binaryToHexBitset(navBin);
	    nav_page[i].navMsg=navHex;
            
            macseq.erase(macseq.begin());
            if(macseq.empty()) return;
        }
        
    }
}
void printSequence(const vector<vector<string>>& seq) {
	
    cout<<"MAC Lookup Table (MACLT)::"<<endl;
    if (seq.empty()) {
        cout << "❌ 未找到对应ID的数据！" << endl;
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
    if (argc < 3)
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
    int result;
    double duration;
    datetime_t t0;
    double gain = 0;
    std::string device_args = "";
    size_t channel = 1;
    bool verbose = false;
    int return_code = EXIT_SUCCESS;
    char error_string[512];
    double rate = TX_SAMPLERATE;
    double freq = TX_FREQUENCY;
    bool use_usrp = true;
    // Buffer sizes
    size_t samps_per_buff = SAMPLES_PER_BUFFER;
    float *buffer = NULL;
    const void **buffer_ptr = NULL;

    while ((result = getopt(argc, argv, "e:n:o:u:g:l:T:t:s:d:G:a:p:iI:U:b:v")) != -1)
    {
        switch (result)
        {
        case 'e':
            strcpy(s.opt.navfile, optarg);
            break;

        case 'n':
            strcpy(s.opt.tvfile, optarg);
            break;

        case 'o':
            strcpy(s.opt.outfile, optarg);
            break;

        case 'l':
            sscanf(optarg, "%lf,%lf,%lf", &s.opt.llh[0], &s.opt.llh[1],
                   &s.opt.llh[2]);
            break;

        case 'T':
            s.opt.timeoverwrite = TRUE;
            if (strncmp(optarg, "now", 3) == 0)
            {
                time_t timer;
                struct tm *gmt;

                time(&timer);
                gmt = gmtime(&timer);

                t0.y = gmt->tm_year + 1900;
                t0.m = gmt->tm_mon + 1;
                t0.d = gmt->tm_mday;
                t0.hh = gmt->tm_hour;
                t0.mm = gmt->tm_min;
                t0.sec = (double)gmt->tm_sec;

                date2gal(&t0, &s.opt.g0);

                break;
            }

        case 't':
        
            sscanf(optarg, "%d/%d/%d,%d:%d:%lf", &t0.y, &t0.m, &t0.d, &t0.hh, &t0.mm,
                   &t0.sec);
            if (t0.y <= 1980 || t0.m < 1 || t0.m > 12 || t0.d < 1 || t0.d > 31 ||
                t0.hh < 0 || t0.hh > 23 || t0.mm < 0 || t0.mm > 59 || t0.sec < 0.0 ||
                t0.sec >= 60.0)
            {
                printf("ERROR: Invalid date and time.\n");
                exit(1);
            }
            
            t0.sec = floor(t0.sec);
            date2gal(&t0, &s.opt.g0);
            break;

        case 'd':
            duration = atof(optarg);

            numSub = static_cast<int>(duration);

            s.opt.iduration = (int)(duration * 10.0 + 0.5);
            break;

        case 'G':
            gain = atof(optarg);
            if (gain < 0)
                gain = 0;
            if (gain > 60)
                gain = 60;
            break;

        case 'a':
            device_args = strdup(optarg);
            break;

        case 'p':
            // printf("%s\n", optarg);
            s.udp_port = atoi(optarg);

            break;

        case 'i':
            s.opt.interactive = TRUE;
            break;

        case 'I':
            s.opt.iono_enable = FALSE; // Disable ionospheric correction
            break;

        case 'U':
            s.opt.use_usrp = false;
            use_usrp = false;
            break;

        case 'b':
            s.opt.use_bit_stream = false;
            break;

        case 'v':
            s.opt.verb = true;
            break;
            
        case 's':
            PATH=std::string(optarg);
            break;

        case ':':

        case '?':
            usage(argv[0]);
            exit(1);

        default:
            break;
        }
    }

    if (s.opt.navfile[0] == 0 && s.opt.tvfile[0] == 0)
    {
        std::cout << s.opt.navfile;
        std::cout << s.opt.tvfile;
        printf("ERROR: Galileo ephemeris/nav_msg file is not specified.\n");
        exit(1);
    }

    if (s.opt.outfile[0] == 0)
    {
        printf("[+] File sink not specified. Using galileosim.bin\n");
        snprintf(s.opt.outfile, sizeof(s.opt.outfile), "galileosim.ishort");
    }
	printf("\nCreating Galileo task...\n");
    // Initialize simulator
    numSub=numSub/30;
    	//cout<<"numSub:"<<numSub<<endl;
    init_sim(&s);
    sim_t s2=s;
    //cout<<"--------------"<<cur_wn<<endl;
    getWeekTow();
    galileo_task1(&s2);

    //getNavMsg(); //getpage
   //cout<<PATH<<endl;
    //cout<<"----asdasdas-------"<<cur_wn<<endl;
    vector<int> cur_prn;
    //cout<<nav_page.size()<<endl;
   //show2(nav_page);
    int prnNum=getNumberPrn(cur_prn);
    cout<<"prnNum: "<<prnNum<<endl;
    
    
        size_t pos = PATH.find_last_of('/');

    // 2. 提取文件名
    string filename = PATH.substr(pos + 1);

    // 3. 去掉前两个字符（"07"）
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
    
    for(int i=0;i<prnNum;i++)
    {
       //cout<<cur_prn[i]<<endl;
       
    }
    cout << "selected PRN :[";
	for (int i = 0; i < prnNum; i++)
	{
	    cout << cur_prn[i];
	    nav_page.erase(nav_page.begin());
	    if (i != prnNum - 1)
		cout << ",";
	}
	cout << "]" << endl;
            //cout<<"---------------------------------------------"<<endl;
    
        cout<<"================= OSNMA Configuration ================="<<endl;

    int fl=0;

    

    for(int i=0;i<prnNum;i++)
    {
    	
    	int prn = cur_prn[i];
    	
    	//cout<<"Generate navigation data for prn "<<prn<<endl;
    	vector<page_struct> cur_page;
    
    	fixNavPage(prn,cur_page);
	
    	fillOSNMA(prrn,cur_page,PATH);
    	if(fl==0){
	getDic(cur_page);
	fl=1;
	
	findSeqByID(dic);
	printSequence(mltSequence);
	KL=ks_num[KS];
	TL=ts_num[TS-5];
	TL=TL;
	cout<<"\nMessage Authentication:"<<endl;
	cout<<"Hash Function :SHA-256"<<endl;
	cout<<"Key_size      :"<<KL<<" bits"<<endl;
	cout<<"Tag_size      :"<<TL<<" bits"<<endl;
	//cout<<"HF: "<<HF<<"MF: "<<MF<<"dic: "<<dic<<" KL:"<<KL<<" TL:"<<TS<<endl;
		    cout<<"======================================================\n\n"<<endl;
	    cout<<"=========== Galileo I/NAV Message Generation =========="<<endl;
	    cout<<"[INFO] Initializing navigation message generation...\n"<<endl;

	}
	
	cleanTag(cur_page);
	if(prrn!=30){
		fillTagInfo(prn,cur_page);
	}
    	addWord10(cur_page);
    	cout<<"---------------------------------------------"<<endl;
    	cout<<"[PRN "<<prn<<"]"<<endl;
    	cout<<"- Loading HKROOT and TESLA Key from auxiliary OSNMA data"<<endl;
   	if(duration>361)
   	{
   		cout<<"- Generating authentication tags (ADKD: 0, 4, 12)"<<endl;
   	}
   	else{
   	cout<<"- Generating authentication tags (ADKD: 0, 4)"<<endl;
   	
   	}
   	
        vector<string> teslaKey=getKey(cur_page);
        //----------------------mac-------------------------
        vector<string> teslaKeyformac=teslaKey;
        
        teslaKeyformac.erase(teslaKeyformac.begin());
        
        vector<string> macauthData=getmacauth(prn,cur_page);
        
        vector<string> macList=commac(teslaKeyformac,macauthData);

	macList.erase(macList.begin());
	fillmecseq(macList,cur_page);
	
        //----------------------mac-------------------------
        teslaKey.erase(teslaKey.begin());
        teslaKey.erase(teslaKey.begin());

        vector<string> authDataADKD4 = getAuthDataForADKD4(prn,cur_page);
        vector<string> teslaKeyforADKD4=teslaKey;
        
        vector<string> KeyforADKD4=keyforADKD4(key,start_time);
        
	
		teslaKeyforADKD4.erase(teslaKeyforADKD4.begin());
	
	
	//cout<<teslaKeyforADKD4.size()<<"++++++"<<KeyforADKD4.size()<<endl;
	for(int i = 0; i < teslaKeyforADKD4.size(); i++) {
	  // cout<<"1:"<<teslaKeyforADKD4[i]<<endl;
	   	   //cout<<"2:"<<KeyforADKD4[i]<<endl;
	}
	
	
	for(int i = authDataADKD4.size() - 1; i >= 0; i--) {
	    if(i % 2 == 1) { // 奇数下标
		teslaKeyforADKD4.erase(teslaKeyforADKD4.begin() + i);
		authDataADKD4.erase(authDataADKD4.begin() + i);
		KeyforADKD4.erase(KeyforADKD4.begin() + i);
	    }
	}
	//cout<<"authchangdu:"<<authDataADKD4.size()<<","<<KeyforADKD4.size()<<endl;
	
	vector<string> tagListADKD4=computeTagADKD4(KeyforADKD4,authDataADKD4);
	if(dic==33)
	{
		tagListADKD4.insert(tagListADKD4.begin(), "bfab8f967d");
	}
	//-----------------------ADKD12-----------------------

        
	vector<string> teslaKeyforADKD12=keyforADKD12(key,start_time);

	    for (size_t i = 0; i < teslaKeyforADKD12.size(); i++)
	   {
	    /* code */
		//cout<<"key: "<<key12[i]<<endl;
	   }
	   
        
	vector<string> authDataADKD12=getAuthDatADKD12(cur_page,prn);
        
	
	
	vector<string> tagListADKD12=computeTagADKD12(teslaKeyforADKD12,authDataADKD12);
	
     
    	for(int i = tagListADKD12.size() - 1; i >= 0; i--) {
	    if(i % 2 == 1) { // 奇数下标
		//tagListADKD12.erase(tagListADKD12.begin() + i);
	    }
	}
	fillTagADKD12(cur_page,tagListADKD12);
	
	//-----------------------ADKD12-----------------------


        vector<string> authData=getAuthData(prn ,cur_page);
	
        vector<string> tagList=computeTag(teslaKey,authData);
        
	fillTagADKD4(cur_page,tagListADKD4);
	
	 
        fillTag(tagList,cur_page);
        //cout<<"Test:"<<tagList.size()<<endl;
        //show2(cur_page);
        fillCrc(cur_page);
        
cout<<" ✔ Completed"<<endl;
        singleIn(prn,cur_page);
        
    }
    cout<<"---------------------------------------------"<<endl;
    for(int i=0;i<prnNum;i++)
    {
    	//cout<<"Generate osnma data for prn "<<cur_prn[i]<<endl;
    }
    //cout<<start_time<<"++++++++++++++++"<<endl;
    
    
    //show(nav_page);

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
    s.status = start_galileo_task(&s);
    if (s.status < 0)
    {
        fprintf(stderr, "Failed to start Galileo task.\n");
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
            fprintf(stderr, "Failed to start TX task.\n");
            // goto out;
        }
        else
            printf("Creating TX task...\n");

        // Running...
        printf("Running...\n"
               "Press Ctrl+C to abort.\n");

        // Wainting for TX task to complete.
        pthread_join(s.tx.thread, NULL);
    }
    else
        pthread_join(s.galileo_sim.thread, NULL);

    fprintf(stderr, "\nTotal samples consumed: %ld", samples_consumed);
    printf("\nDone!\n");
}
