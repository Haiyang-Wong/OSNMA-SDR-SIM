#include "../include/galileo-sdr.h"

static bool readMackBits(const vector<page_struct> &cur_page, size_t subframeStart, int bitStart, int bitLen, string &out)
{
    out.clear();
    for (int bit = 0; bit < bitLen; bit++)
    {
        int globalBit = bitStart + bit;
        size_t pageIndex = subframeStart + static_cast<size_t>(globalBit / 32);
        int pageBit = 146 + (globalBit % 32);
        if (pageIndex >= cur_page.size())
            return false;

        string navBin = hexToBinary(cur_page[pageIndex].navMsg);
        if (pageBit < 0 || pageBit >= static_cast<int>(navBin.size()))
            return false;
        out.push_back(navBin[pageBit]);
    }
    return true;
}

static bool writeMackBits(vector<page_struct> &cur_page, size_t subframeStart, int bitStart, const string &bits)
{
    for (size_t bit = 0; bit < bits.size(); bit++)
    {
        int globalBit = bitStart + static_cast<int>(bit);
        size_t pageIndex = subframeStart + static_cast<size_t>(globalBit / 32);
        int pageBit = 146 + (globalBit % 32);
        if (pageIndex >= cur_page.size())
            return false;

        string navBin = hexToBinary(cur_page[pageIndex].navMsg);
        if (pageBit < 0 || pageBit >= static_cast<int>(navBin.size()))
            return false;

        navBin[pageBit] = bits[bit];
        cur_page[pageIndex].navMsg = binaryToHexBitset(navBin);
    }
    return true;
}

static string computeNormalAdkd0Tag(const string &key, const string &tag0AuthData, int counter)
{
    if (key.empty() || tag0AuthData.size() < 50)
        return "";

    string prnBits = tag0AuthData.substr(0, 8);
    string authData = prnBits + tag0AuthData.substr(0, 40) + bitset<8>(counter).to_string() + tag0AuthData.substr(48);
    vector<uint8_t> keyBytes = hexstr_to_bytes(key);
    vector<uint8_t> authBytes = bitstring_to_bytes(authData);
    if (keyBytes.empty() || authBytes.empty())
        return "";

    vector<uint8_t> mac = hmac_sha256(keyBytes, authBytes);
    return macToString(mac).substr(0, 10);
}

void fillTagADKD0ExtraSlots(const vector<string> &keyList,
                            const vector<string> &authData,
                            vector<page_struct> &cur_page)
{
    size_t tagIndex = 0;
    int filledCount = 0;

    for (size_t i = 0; i < cur_page.size(); i++)
    {
        int tow = static_cast<int>(cur_page[i].tow);
        if (getPageNumInSub(tow) != 0 || tow <= start_time)
            continue;

        if (tagIndex >= keyList.size() || tagIndex >= authData.size())
            return;

        int prn = cur_page[i].prn;
        string selfAdkd0Info = bitset<8>(prn).to_string() + "0000" + "0001";
        for (int slot = 1; slot <= 5; slot++)
        {
            string tagInfo;
            int infoStart = slot * 56 + 40;
            if (!readMackBits(cur_page, i, infoStart, 16, tagInfo))
                continue;
            if (tagInfo != selfAdkd0Info)
                continue;

            string tagHex = computeNormalAdkd0Tag(keyList[tagIndex], authData[tagIndex], slot + 1);
            if (tagHex.empty())
                continue;

            string tagBits = hexToBinary(tagHex);
            if (writeMackBits(cur_page, i, slot * 56, tagBits.substr(0, 40)))
            {
                recordGeneratedOsnmaTag(static_cast<int>(cur_page[i].tow),
                                        static_cast<int>(cur_page[i].tow),
                                        prn,
                                        prn,
                                        0,
                                        false,
                                        "self_adkd0_slot",
                                        slot);
                filledCount++;
            }
        }

        tagIndex++;
    }

    if (filledCount > 0)
    {
        cout << "[OSNMA] filled " << filledCount
             << " extra self ADKD0 MACK slots on PRN " << cur_page.front().prn << endl;
    }
}

void fillTag(vector<string> tagList,vector<page_struct> &cur_page)
{
	    int tag0CopSetCount = 0;
	    for (size_t i = 0; i < cur_page.size()&&!tagList.empty(); i++)
	    {
        int tow=static_cast<int>(cur_page[i].tow);
        // cout<<tow<<endl;
        //cout<<"cur_page[i].navMsg:1"<<cur_page[i].navMsg<<endl;
        if(getPageNumInSub(tow)==0&&tow>start_time) // 
        {
            
            string navBin=hexToBinary(cur_page[i].navMsg);
            string temp=navBin.substr(0,146);
            string tag=hexToBinary(tagList[0]);
            string tag1=tag.substr(0,32);
            string tag2=tag.substr(32);

            temp+=tag1;
            temp+=navBin.substr(178);
            //cout<<"tag:"<<binaryToHexBitset(tag)<<endl;

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
	            temp.replace(166, 4, "0001");
	            if (temp.substr(166, 4) != "0001") {
	                cerr << "[OSNMA] tag0 COP write failed: PRN " << cur_page[i].prn
	                     << " TOW " << cur_page[i].tow << endl;
	            } else {
	                tag0CopSetCount++;
	            }
		            cur_page[i+1].navMsg=binaryToHexBitset(temp); 
		            recordGeneratedOsnmaTag(static_cast<int>(cur_page[i].tow),
		                                    static_cast<int>(cur_page[i].tow),
		                                    cur_page[i].prn,
		                                    cur_page[i].prn,
		                                    0,
		                                    false,
		                                    "self_adkd0",
		                                    0);
	            
	            
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
