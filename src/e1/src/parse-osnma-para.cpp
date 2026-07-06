#include "../include/galileo-sdr.h"
#include <fstream>   

const int mack_start=146;
int std_time = 277201;
using json = nlohmann::json;
int getFirstTow(int a)
{
    int numPage=getPageNumInSub(a);
    return a-numPage*2;
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
int is60page(int tow)
{
    return (tow-(std_time))%60==0;
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
vector<string> getmacauth(int prn,vector<page_struct> cur_page)
{
    vector<string> authData;
    for(int i=0;i<cur_page.size();i++)
    {
        int tow=static_cast<int>(cur_page[i].tow);
        
        if((tow-277201)%30==0)
        {
            //cout<<"getmacauth:"<<tow<<endl;
            string wn_bin=toBinary12(cur_page[i].wn);
            string tow_bin=toBinary20(tow-1);

            string prn_bin=toBinary8(prn);
            
            string temp=prn_bin+wn_bin+tow_bin;
            if(dic!=33){
            	if(is60page(getFirstTow(tow)))
		    {
		    	temp+="0000000000000000";
		    	temp+="0000000000000000";
		    }
		    else{
		    	temp+="0000000000000000";	
		    }
            }
            authData.push_back(temp);   
            //cout<<temp<<endl;
        }
    }
    return authData;
}
std::string toBinary8(int num) {
    return std::bitset<8>(num).to_string();
}

std::string toBinary12(int num) {
    // bitset<12> converts automatically and pads with zeros to 12 bits.
    return std::bitset<12>(num).to_string();
}

// Convert to a 20-bit binary string.
std::string toBinary20(int num) {
    return std::bitset<20>(num).to_string();
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
void getDic(vector<page_struct> cur_page)
{   
    int  count = 0;
    string dics="";
    dsm_id dsm;
    const size_t page_limit = min<size_t>(120, cur_page.size());
    for (size_t i = 0; i < page_limit; i++)
    {

        string nav_msg_hex=hexToBinary(cur_page[i].navMsg);
        string s1=binaryToHexBitset(nav_msg_hex.substr(138,8));
        dics+=s1;
        
        if ((i+1)%15==0)
        {
            dics=dics.substr(4);
            //cout<<"  dic:"<<dic<<endl;
            if(dics[3]=='0')
            {
                //cout<<"dics:"<<bitset<8>(hexToBinary(dics.substr(6,2))).to_ulong()<<"  all:"<<dics<<endl;
               // cout<<"HF:"<<hexToBinary(string(1, dics[3])).substr(0,2)<<endl;
                //cout<<"MF:"<<hexToBinary(string(1, dics[3])).substr(2,2)<<endl;
                
                HF=bitset<8>(hexToBinary(string(1, dics[3])).substr(0,2)).to_ulong();
                MF=bitset<8>(hexToBinary(string(1, dics[3])).substr(2,2)).to_ulong();
                dic=bitset<8>(hexToBinary(dics.substr(6,2))).to_ulong();
                KS=bitset<8>(hexToBinary(string(1, dics[4]))).to_ulong();
                TS=bitset<8>(hexToBinary(string(1, dics[5]))).to_ulong();
                //cout<<"KS:"<<KS<<"TS:"<<TS<<endl;
                if(dic==33||dic==34||dic==35||dic==29) return;
            }
            if (i!=0)
            {
                dics="";
            }
            
        }

    }
    
}
void findSeqByID(int target_id){
        ifstream file("../include/mlt.json");  

        // Parse JSON directly from the file stream using nlohmann/json.
        json json_array = json::parse(file);
        file.close(); // Close the file after parsing.

    for (const auto& item : json_array) {
        int current_id = item["ID"];
        if (current_id == target_id) {
           mltSequence = item["sequence"];
        }
    }
    return; 
}

