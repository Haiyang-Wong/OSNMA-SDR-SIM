#include "../include/galileo-sdr.h"
void fillTag(vector<string> tagList,vector<page_struct> &cur_page)
{
    for (size_t i = 0; i < cur_page.size()&&!tagList.empty(); i++)
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
