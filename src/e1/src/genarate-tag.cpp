#include "../include/galileo-sdr.h"

#define ROTRIGHT(a,b) (((a) >> (b)) | ((a) << (32-(b))))
const string mapping = "0123456789ABCDEF";
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
vector<string> commac(vector<string> teslaKey, vector<string> authData) {
    vector<string> macseq;
    for(int i=0;i<authData.size()&&i<teslaKey.size();i++)
    {
        std::vector<uint8_t> key = hexstr_to_bytes(teslaKey[i]);
        std::vector<uint8_t> msg = bitstring_to_bytes(authData[i]);
        auto mac = hmac_sha256(key, msg);
        std:: string Mac2=macToString(mac);
        macseq.push_back(Mac2.substr(0,3));
        //std::cout << "HMAC-SHA256: "<<Mac2<<endl;
        //if(i==10) break;
    }
    return macseq;
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
std::string macToString(const std::vector<uint8_t>& mac) {
    std::ostringstream oss;
    for (size_t i = 0; i < mac.size(); ++i) {
        if (i > 0) oss;  
        oss << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(mac[i]);
    }
    return oss.str();
}

vector<string> computeTag(vector<string> key,vector<string> authData) //
{
    vector<string> tagList;
    size_t tagCount = min(key.size(), authData.size());
    for(size_t i=0;i<tagCount;i++)
    {	
        if(authData[i].empty() || key[i].empty())
        {
            tagList.push_back("0000000000");
            continue;
        }
    	//cout<<"key:"<<key[i]<<" -- adkd0---  auth_data:"<<authData[i]<<endl;
        vector<uint8_t>k1=hexstr_to_bytes(key[i]);
	
        vector<uint8_t>d1=bitstring_to_bytes(authData[i]);
	std::vector<uint8_t> mac;
        if (MF==0){
        	mac= hmac_sha256(k1,d1);
        }

        string temp=macToString(mac);
        //cout<<temp<<endl;
        temp=temp.substr(0,10);
       //cout<<"temp_Tag:"<<temp<<endl;
        tagList.push_back(temp);
        
    }
    return tagList;
}
void cleanTag(vector<page_struct>& cur_page)
{
    for (size_t i = 0; i < cur_page.size(); i++)
    {
        int tow = static_cast<int>(cur_page[i].tow);
        int pageNum = getPageNumInSub(tow);

        if (is60page(getFirstTow(tow)))
        {
            if (pageNum==1)
            {
                string navBin=hexToBinary(cur_page[i].navMsg);
                for(int j=170;j<178;j++)
                {
                    navBin[j]='0';
                }   
                cur_page[i].navMsg=binaryToHexBitset(navBin);    
                i++;
                navBin=hexToBinary(cur_page[i].navMsg);
                for(int j=146;j<178;j++)
                {
                    navBin[j]='0';
                }
                cur_page[i].navMsg=binaryToHexBitset(navBin);
                i++;
                navBin=hexToBinary(cur_page[i].navMsg);
                for(int j=146;j<162;j++)
                {
                    navBin[j]='0';
                }
                cur_page[i].navMsg=binaryToHexBitset(navBin);
            }
            if (pageNum==5)
            {
                string navBin=hexToBinary(cur_page[i].navMsg);
                for(int j=154;j<178;j++)
                {
                    navBin[j]='0';
                }
                cur_page[i].navMsg=binaryToHexBitset(navBin);
                i++;
                navBin=hexToBinary(cur_page[i].navMsg);
                for(int j=146;j<178;j++)
                {
                    navBin[j]='0';
                }
                cur_page[i].navMsg=binaryToHexBitset(navBin);
            }
            if (pageNum==8)
            {
                string navBin=hexToBinary(cur_page[i].navMsg);
                for(int j=170;j<178;j++)
                {
                    navBin[j]='0';
                }
                cur_page[i].navMsg=binaryToHexBitset(navBin);
                i++;
                navBin=hexToBinary(cur_page[i].navMsg);
                for(int j=146;j<178;j++)
                {
                    navBin[j]='0';
                }
                cur_page[i].navMsg=binaryToHexBitset(navBin);
                i++;
                navBin=hexToBinary(cur_page[i].navMsg);
                for(int j=146;j<162;j++)
                {
                    navBin[j]='0';
                }
                cur_page[i].navMsg=binaryToHexBitset(navBin);
            }
            
            
        }
        else
        {
            if (pageNum==1)
            {
                string navBin=hexToBinary(cur_page[i].navMsg);
                for(int j=170;j<178;j++)
                {
                    navBin[j]='0';
                }   
                cur_page[i].navMsg=binaryToHexBitset(navBin);    
                i++;
                navBin=hexToBinary(cur_page[i].navMsg);
                for(int j=146;j<178;j++)
                {
                    navBin[j]='0';
                }
                cur_page[i].navMsg=binaryToHexBitset(navBin);
                i++;
                navBin=hexToBinary(cur_page[i].navMsg);
                for(int j=146;j<178;j++)
                {
                    navBin[j]='0';
                }
                cur_page[i].navMsg=binaryToHexBitset(navBin);
                i++;
                navBin=hexToBinary(cur_page[i].navMsg);
                for(int j=146;j<178;j++)
                {
                    navBin[j]='0';
                }
                cur_page[i].navMsg=binaryToHexBitset(navBin);
                i++;
                navBin=hexToBinary(cur_page[i].navMsg);
                for(int j=146;j<154;j++)
                {
                    navBin[j]='0';
                }
                cur_page[i].navMsg=binaryToHexBitset(navBin);
            }
            if (pageNum==7)
            {
                string navBin=hexToBinary(cur_page[i].navMsg);
                for(int j=170;j<178;j++)
                {
                    navBin[j]='0';
                }
                cur_page[i].navMsg=binaryToHexBitset(navBin);
                i++;
                navBin=hexToBinary(cur_page[i].navMsg);
                for(int j=146;j<178;j++)
                {
                    navBin[j]='0';
                }
                cur_page[i].navMsg=binaryToHexBitset(navBin);
                i++;
                navBin=hexToBinary(cur_page[i].navMsg);
                for(int j=146;j<178;j++)
                {
                    navBin[j]='0';
                }
                cur_page[i].navMsg=binaryToHexBitset(navBin);
                i++;
                navBin=hexToBinary(cur_page[i].navMsg);
                for(int j=146;j<162;j++)
                {
                    navBin[j]='0';
                }
                cur_page[i].navMsg=binaryToHexBitset(navBin);
            }
            
        }
        
         
    }
    
}
void fillTagInfo(int prn ,vector<page_struct>& cur_page)
{
    const string selfAdkd0 = bitset<8>(prn).to_string() + "0000" + "0001";
    const string selfAdkd4 = bitset<8>(prn).to_string() + "0100" + "0001";
    const string selfAdkd12 = bitset<8>(prn).to_string() + "1100" + "0001";

    for (size_t i = 0; i < cur_page.size(); i++)
    {
        int tow= static_cast<int>(cur_page[i].tow);
        int pageNum = getPageNumInSub(tow);
        string navBin=hexToBinary(cur_page[i].navMsg);

        if (is60page(getFirstTow(tow)))
        {
            if (pageNum==3)
            {
                navBin.replace(146,16,selfAdkd0);
            }
            if (pageNum==4)
            {
                navBin.replace(170,8,selfAdkd4.substr(0, 8));
            }
            if (pageNum==5)
            {
                navBin.replace(146,8,selfAdkd4.substr(8, 8));
            }
            if (pageNum==6)
            {
                navBin.replace(162,16,selfAdkd0);
            }
            if (pageNum==8)
            {
                navBin.replace(154,16,selfAdkd12);
            }
            if (pageNum==10)
            {
                navBin.replace(146,16,selfAdkd0);
            }
        }
        else{
            if (pageNum==3)
            {
                navBin.replace(146,16,selfAdkd0);
            }
            if (pageNum==4)
            {
                navBin.replace(170,8,selfAdkd0.substr(0, 8));
            }
            if (pageNum==5)
            {
                navBin.replace(146,8,selfAdkd0.substr(8, 8));
            }
            if (pageNum==6)
            {
                navBin.replace(162,16,selfAdkd12);
            }
            if (pageNum==8)
            {
                navBin.replace(154,16,selfAdkd0);
            }
            if (pageNum==10)
            {
                navBin.replace(146,16,selfAdkd12);
            }
        }

        cur_page[i].navMsg=binaryToHexBitset(navBin);
    }
}
