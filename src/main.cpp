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
string PATH="";
int subNub=11;
double start_time = 1000000000.0;
int cur_wn = 0;
bool stop_signal_called = false;

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

struct GST{
    string tow;
    string wn;
};
typedef struct GST gst;

const int osnma_start=138;
const int mack_start=146;
const int tagSize=40;
int len_tow=20;
const int len_gst=32;

vector<string> getKey(vector<page_struct> &cur_page);
int getPageNumInSub(int);
string getMackInPage(string);
vector<string> getAuthData(int prn, vector<page_struct> &cur_page);
vector<gst> getGst(vector<page_struct> &cur_page);
vector<string> getNavData(vector<page_struct> &cur_page);
vector<string> computeTag(vector<string>,vector<string>);
void fillTag(vector<string>,vector<page_struct> &cur_page);
int getNumberPrn(vector<int> &cur_prn);

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
string hexToBinary(const string& hexStr);

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
std::vector<uint8_t> hmac_sha256(const std::vector<uint8_t> &key, const std::vector<uint8_t> &message) {
    const size_t BLOCK_SIZE = 64;
    std::vector<uint8_t> k_pad(BLOCK_SIZE, 0);

    if (key.size() > BLOCK_SIZE) {
        SHA256 sha;
        uint8_t hashedKey[32];
        sha.update(key.data(), key.size());
        sha.final(hashedKey);
        std::copy(hashedKey, hashedKey + 32, k_pad.begin());
    } else {
        std::copy(key.begin(), key.end(), k_pad.begin());
    }

    std::vector<uint8_t> o_key_pad = k_pad;
    std::vector<uint8_t> i_key_pad = k_pad;
    for (size_t i = 0; i < BLOCK_SIZE; ++i) {
        o_key_pad[i] ^= 0x5c;
        i_key_pad[i] ^= 0x36;
    }

    SHA256 shaInner;
    shaInner.update(i_key_pad.data(), i_key_pad.size());
    shaInner.update(message.data(), message.size());
    uint8_t innerHash[32];
    shaInner.final(innerHash);

    SHA256 shaOuter;
    shaOuter.update(o_key_pad.data(), o_key_pad.size());
    shaOuter.update(innerHash, 32);
    uint8_t finalHash[32];
    shaOuter.final(finalHash);

    return std::vector<uint8_t>(finalHash, finalHash + 32);
}
std::vector<uint8_t> hexstr_to_bytes(const std::string& hex) {
    std::vector<uint8_t> bytes;
    if (hex.size() % 2 != 0) return bytes;  

    for (size_t i = 0; i < hex.size(); i += 2) {
        uint8_t byte = std::stoi(hex.substr(i, 2), nullptr, 16);
        bytes.push_back(byte);
    }
    return bytes;
}
std::vector<uint8_t> bitstring_to_bytes(const std::string& bits) {
    std::vector<uint8_t> bytes;
    size_t len = bits.size();

    for (size_t i = 0; i < len; i += 8) {
        std::string byte_str = bits.substr(i, std::min((size_t)8, len - i));
        while (byte_str.size() < 8) byte_str += '0';  
        uint8_t byte = static_cast<uint8_t>(std::stoi(byte_str, nullptr, 2));
        bytes.push_back(byte);
    }

    return bytes;
}
std::string macToString(const std::vector<uint8_t>& mac) {
    std::ostringstream oss;
    for (size_t i = 0; i < mac.size(); ++i) {
        if (i > 0) oss;  
        oss << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(mac[i]);
    }
    return oss.str();
}
//-------------------------hash  256-----------------------------
void fixNavPage(int a,vector<page_struct> &cur_prn);
int getNumberPrn(vector<int> &cur_prn)
{
    int prn[50]={0};
    for(int i = 0;i<100;i++)
    {
    	int s=nav_page[i].prn;
    	prn[s]++;
    }
    int cont=0;
    for(int i = 0;i<50;i++)
    {
    	if(prn[i]!=0)
    	{
    	    cont++;
    	    cur_prn.push_back(i);
    	}
    }
    return cont;
}
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
void  fillOSNMA(int prn,vector<page_struct> &cur_page)
{
    string line1;
    ifstream inf;
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
        if(stoi(columns[0])<static_cast<int>(start_time))
        {
        	//cout<<columns[0]<<","<<start_time<<endl;
        	continue;
        }
        
        if (stoi(columns[2])!= prn)
            continue;
        
        else // Extract the page corresponding to the satellite
        {
            page_struct a;
            a.tow = stod(columns[0]);
            a.wn = stoi(columns[1]);
            a.prn = stoi(columns[2]);
            columns[3].pop_back();
            a.navMsg = columns[3];
            
            osnma_r.push_back(a);
            
           
        }
    }
    inf.close();
    // show(osnma_r);

    //get the osnma data from osnam_r and fill it in osnma_c interesting
    for(int i=0;i<cur_page.size();i++) 
    {
   // cout<<osnma_r[i].navMsg<<"fillOSNMA"<<endl;
        //cout<<osnma_r[i].navMsg<<endl;
        string nav_mas_bing_r=hexToBinary(osnma_r[i].navMsg);
        
        string nav_mas_bing_c=hexToBinary(cur_page[i].navMsg);
        string osnma_bin=nav_mas_bing_r.substr(138,40);  

        nav_mas_bing_c.replace(138,40,osnma_bin);
        
        string nav_mes_hex=binaryToHexBitset(nav_mas_bing_c);
        //exit(1);
        cur_page[i].navMsg=nav_mes_hex;
        
        
    }
}
void fixNavPage(int a,vector<page_struct> &cur_prn)
{
    for(int i =0 ;i<nav_page.size();i++)
    {
    	if(nav_page[i].prn==a)
    	{
    		cur_prn.push_back(nav_page[i]);
    		//cout<<nav_page[i].tow<<"++++++++++"<<endl;
    	}
    }
}
int getPageNumInSub(int tow)
{
    return ((tow-static_cast<int>(start_time))/2)%15;
}
string getMackInPage(string nav_msg_hex)  // in nav_msg
{
    string nav_msg_bin=hexToBinary(nav_msg_hex);
    string mack_bin=nav_msg_bin.substr(mack_start,32);
    string mack_hex=binaryToHexBitset(mack_bin);

    return mack_hex;
}
vector<string> getKey(vector<page_struct> &cur_page)
{
    //In each subframe, the position of the key is (10,15), (16,144)
    vector<string> teslaKey;
    string mack;
  
    for(int i=0;i<cur_page.size();i++)
    {
        int tow=cur_page[i].tow;
        if(getPageNumInSub(tow)>=10&&getPageNumInSub(tow)<=14) //Starting from here, extract tags for the next five pages
        {
            string temp=getMackInPage(cur_page[i].navMsg);
            
            mack.append(temp);
        }
        if(getPageNumInSub(tow)==14)
        {   
            teslaKey.push_back(mack);
            mack.clear();
        }
    }

    vector<string> teslaKey_bit;
    // teslaKey.erase(teslaKey.begin());
    // teslaKey.erase(teslaKey.begin());
    while(!teslaKey.empty())
    {
        string key=teslaKey.front();
        key=key.substr(4,32);
        //cout<<"key:"<<key<<endl;
        teslaKey_bit.push_back(key);
        teslaKey.erase(teslaKey.begin());
    }
    return teslaKey_bit;
}
vector<string> getNavData(vector<page_struct> &cur_page)
{   

    vector<string> words(5);

    string mack;


    for (size_t i = 0; i < cur_page.size(); i++)
    {
        int tow=static_cast<int>(cur_page[i].tow);
        if(getPageNumInSub(tow)==10)
        {
            string navBin=hexToBinary(cur_page[i].navMsg);
            string s1=navBin.substr(8,106)+navBin.substr(122,14);
            //cout<<s1<<endl;
            words[0]=s1;
        }
        if(getPageNumInSub(tow)==0)
        {
            string navBin=hexToBinary(cur_page[i].navMsg);
            string s1=navBin.substr(8,106)+navBin.substr(122,14);
            words[1]=s1;
        }
        if (getPageNumInSub(tow)==11)
        {
            string navBin=hexToBinary(cur_page[i].navMsg);
            string s1=navBin.substr(8,106)+navBin.substr(122,16);
            words[2]=s1;
        }
        if (getPageNumInSub(tow)==1)
        {
            string navBin=hexToBinary(cur_page[i].navMsg);
            string s1=navBin.substr(8,106)+navBin.substr(122,14);
            words[3]=s1;
        }
        if (getPageNumInSub(tow)==12)
        {
            string navBin=hexToBinary(cur_page[i].navMsg);
            string s1=navBin.substr(8,67);
            words[4]=s1;
        }
        if (words[4]!="")
        {
            return words;
        }
    }

    return words;
}
vector<gst> getGst(vector<page_struct> &cur_page)
{
    vector<gst> GST;   
    for (size_t i = 0; i < cur_page.size(); i++)
    {
        int tow = static_cast<int>(cur_page[i].tow);
        if((tow>static_cast<int>(start_time)+28)&&(tow-(static_cast<int>(start_time)+30))%30==0)
        {
            gst temp;
            temp.tow=to_string(tow-1);
            //cout<<columns[0]<<endl;
            temp.wn=to_string(cur_page[i].wn);
            GST.push_back(temp);
        }
    }
    return GST;
}
vector<string> getAuthData(int prn, vector<page_struct> &cur_page) //Should return a binary string
{
    //What Tag0 needs，self.prn_a + self.gst_subframe.bitarray + BitArray(uint=self.ctr,length=8) + self.nma_status + self.nav_data.nav_data_stream
    vector<gst> GST=getGst(cur_page);
    //  for(int i=0;i<GST.size();i++)
    // {
    //     cout<<GST[i].tow<<","<<GST[i].wn<<endl;
    // }
    vector<string> navData = getNavData(cur_page);

    vector<string> authData;
    string navDataBin;
    for(int i=0;i<navData.size();i++)
    {   
        navDataBin+=navData[i];
    }
    //cout<<navDataBin<<endl;
    for(int i=0;i<subNub;i++)
    {
        string temp;
        string sss="";
        
        if(prn<10)
        {
        	sss="0"+to_string(prn);
        }
        else if(prn<30)
        {
        	sss=to_string(prn-6);
        }
        else
        {
        sss=to_string(prn-12);}
        
        string prn_a=hexToBinary(sss);
        string bitAr=hexToBinary("01");
        string nav_status="01";
        bitset<len_gst> bitArray((stoi(GST[i].wn) << len_tow | stoi(GST[i].tow)));
        string gst_bin=bitArray.to_string();
        // cout<<gst_bin<<endl;
        temp=prn_a+gst_bin+bitAr+nav_status+navDataBin;
        authData.push_back(temp);
        //cout<<temp<<endl;
    }

    return authData;
}
vector<string> computeTag(vector<string> key,vector<string> authData) //
{
    vector<string> tagList;
    for(int i=0;i<subNub;i++)
    {
        vector<uint8_t>k1=hexstr_to_bytes(key[i]);

        vector<uint8_t>d1=bitstring_to_bytes(authData[i]);

        auto mac= hmac_sha256(k1,d1);

        string temp=macToString(mac);
        temp=temp.substr(0,10);
       // cout<<temp<<endl;
        tagList.push_back(temp);
    }
    return tagList;
}
void fillTag(vector<string> tagList,vector<page_struct> &cur_page)
{
    for (size_t i = 0; i < cur_page.size(); i++)
    {
        int tow=static_cast<int>(cur_page[i].tow);
        // cout<<tow<<endl;
        //cout<<"cur_page[i].navMsg:1"<<cur_page[i].navMsg<<endl;
        if(getPageNumInSub(tow)==0&&tow>start_time) // 
        {
            
            string navBin=hexToBinary(cur_page[i].navMsg);
            string osnma_bin = navBin.substr(138, 40);


int zero_count = count(osnma_bin.begin(), osnma_bin.end(), '0');
if (zero_count > 30) {
    continue; 
}
            string temp=navBin.substr(0,146);
            string tag=hexToBinary(tagList[0]);
            string tag1=tag.substr(0,32);
            string tag2=tag.substr(32);

            temp+=tag1;
            temp+=navBin.substr(178);
            // cout<<temp<<endl;

            cur_page[i].navMsg=binaryToHexBitset(temp); //The current index i has been filled in, now fill in i+1
            if(i+1>=cur_page.size())
            {
                return;
            }
            
            tow=static_cast<int>(cur_page[i+1].tow);

            navBin=hexToBinary(cur_page[i+1].navMsg);
            temp=navBin.substr(0,146);
            temp+=tag2;
            temp+=navBin.substr(154);
            cur_page[i+1].navMsg=binaryToHexBitset(temp); 
            
            
            if(!tagList.empty())
            {
                tagList.erase(tagList.begin());
            }
            else{
                return;
            }
        }
        //cout<<"cur_page[i].navMsg:2"<<cur_page[i].navMsg<<endl;
    }
}
void binStringToBinInt(const std::string &binaryStr, int *bitArray) {
    if (binaryStr.length() != 240) {
        std::cerr << "Error: Input string must be exactly 240 bits long." << std::endl;
        return;
    }

    for (int i = 0; i < 240; ++i) {
        if (binaryStr[i] == '0') {
            bitArray[i] = 0;
        } else if (binaryStr[i] == '1') {
            bitArray[i] = 1;
        } else {
            std::cerr << "Error: Invalid character '" << binaryStr[i] << "' at position " << i << std::endl;
            return;
        }
    }
}
void trsData(int *page,int *data)
{
    

    for(int i=0;i<114;i++)
    {
        data[i]=page[i];
    }
    for (int i=114;i<196;i++)
    {
        data[i]=page[i+6];
    }
    
}

void encode_int_to_bits(int *page, unsigned int value, int num_bits) {
    // Make sure the value fits within the specified number of bits
    value &= (1 << num_bits) - 1;

    // for (int i = 0; i < num_bits; ++i) {
    //     page[(num_bits - i - 1)] = ((((value >> i) & 1) == 1) ? 1:0);
    // }
    int index = 202;

    for (int j = num_bits-1; j >= 0; j--)
    {
        page[index] = BIT_ISSET(value, j);
        index++;
    }
}
unsigned int Crc24qEncode1(int *BitStream, int Length)
{
    unsigned int crc_result = 0;
    unsigned char dataByte = 0;

    int j;
    for (j = 0; j < Length; j++)
    {
        if (j > 0 && j % 8 == 0)
        {
            crc_result = ((crc_result | dataByte) << 8) ^ Crc24q[(crc_result >> 24) & 0xFF];
            dataByte = 0;
        }

        dataByte = (dataByte << 1) | (char)BitStream[j];
    }

    unsigned char trail_bits_count = j % 8;
    dataByte <<= (8U - trail_bits_count);
    crc_result = ((crc_result | dataByte) << trail_bits_count) ^ Crc24q[(crc_result >> (32 - trail_bits_count)) & 0xFF];

    for (int i = 0; i < 3; i++)
    {
        crc_result = (crc_result << 8) ^ Crc24q[(crc_result >> 24) & 0xFF];
    }
    return (crc_result >> 8);
}
void fillCrc(vector<page_struct> &cur_page)
{
    for (size_t i = 0; i < cur_page.size(); i++)
    {
        //cout<<"cur_page[i].navMsg:"<<cur_page[i].navMsg<<endl;
        string navBin=hexToBinary(cur_page[i].navMsg);
        int binBit[240];
        int data[196]={};
        binStringToBinInt(navBin,binBit);
        trsData(binBit,data);
        unsigned int crcq = Crc24qEncode1(data, 196);
        encode_int_to_bits(binBit,crcq,24);
        navBin.clear();
        for(int i=0;i<240;i++)
        {
            //cout<<binBit[i];
            if(binBit[i]==0)
            {
                navBin.push_back('0');
            }
            else
            {
                navBin.push_back('1');
            }
        }
        string navdata=binaryToHexBitset(navBin);
        cur_page[i].navMsg=navdata;
    }
}
void singleIn(int a,vector<page_struct> cur_page)
{
    for (size_t i = 0; i < nav_page.size(); i++)
    {
        if (nav_page[i].prn==a&&static_cast<int>(nav_page[i].tow)>start_time-1)
        {
           
            
            for (size_t j = 0; j < cur_page.size(); j++)
            {
                int tow1=static_cast<int>(nav_page[i].tow);
                int tow2=static_cast<int>(cur_page[j].tow);
                
                if(tow1==tow2)
                {
                      nav_page[i].navMsg=hexToBinary(cur_page[j].navMsg);
                }

            }
            
        }

        
    }
}
void show2(vector<page_struct> nav_ss)
{
    for(int i=0;i<nav_ss.size();i++)
    {
        cout<<nav_ss[i].tow<<","<<cur_wn<<","<<nav_ss[i].prn<<","<<nav_ss[i].navMsg<<endl;
    }
    cout<<"----------------2----------------"<<endl;
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
    init_sim(&s);
    sim_t s2=s;
    //cout<<"--------------"<<cur_wn<<endl;
    getWeekTow();
    galileo_task1(&s2);

    //getNavMsg(); //getpage
   // cout<<PATH<<endl;
    
    vector<int> cur_prn;
    //cout<<nav_page.size()<<endl;
   // show2(nav_page);
    int prnNum=getNumberPrn(cur_prn);
    
    
    for(int i=0;i<prnNum;i++)
    {
       //cout<<cur_prn[i]<<endl;
       nav_page.erase(nav_page.begin());
    }
    cout<<"Generate Galileo I/NAV message ... "<<endl;
    for(int i=0;i<prnNum;i++)
    {
    	
    	int prn = cur_prn[i];
    	
    	//cout<<"Generate navigation data for prn "<<prn<<endl;
    	vector<page_struct> cur_page;
    	//cout<<"ssssssssssssssssssssss"<<endl;
    	fixNavPage(prn,cur_page);

    	fillOSNMA(prn,cur_page);
    	cout<<"Generate navigation message for prn :"<<prn<<endl;
   	//if(prn==8)
   	// show(cur_page);
        vector<string> teslaKey=getKey(cur_page);
        teslaKey.erase(teslaKey.begin());
        teslaKey.erase(teslaKey.begin());

        vector<string> authData=getAuthData(prn ,cur_page);

        vector<string> tagList=computeTag(teslaKey,authData);

        fillTag(tagList,cur_page);
        
        fillCrc(cur_page);
	
        singleIn(prn,cur_page);
        
    }
    cout<<"---------------------------------------------"<<endl;
    for(int i=0;i<prnNum;i++)
    {
    	cout<<"Generate osnma data for prn "<<cur_prn[i]<<endl;
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
