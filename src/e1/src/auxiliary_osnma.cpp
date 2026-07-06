#include "../include/galileo-sdr.h"
#include<string>
#include<map>

namespace {

const string kZeroOsnma40(40, '0');

struct OsnmaRecord
{
    int tow;
    string osnma40;
};

string trimCsvField(string value)
{
    while (!value.empty() && (value.back() == '\r' || value.back() == '\n' || value.back() == ' ' || value.back() == '\t'))
    {
        value.pop_back();
    }
    size_t first = 0;
    while (first < value.size() && (value[first] == ' ' || value[first] == '\t'))
    {
        first++;
    }
    return value.substr(first);
}

bool isBinaryString(const string &value, size_t expectedSize)
{
    if (value.size() != expectedSize)
    {
        return false;
    }
    for (size_t i = 0; i < value.size(); i++)
    {
        if (value[i] != '0' && value[i] != '1')
        {
            return false;
        }
    }
    return true;
}

vector<OsnmaRecord> readOsnmaKeychain(string PATH)
{
    string line1;
    ifstream inf;
    inf.open(PATH);

    vector<OsnmaRecord> records;

    while(getline(inf,line1))
    {
        if (line1.empty())
        {
            continue;
        }

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

        int tow;
        try
        {
            tow = stoi(trimCsvField(columns[0]));
        }
        catch (...)
        {
            continue;
        }

        if(tow < static_cast<int>(start_time))
        {
            continue;
        }

        string osnma40;
        if (columns.size() == 2 && isBinaryString(trimCsvField(columns[1]), 40))
        {
            osnma40 = trimCsvField(columns[1]);
        }
        else if (columns.size() >= 4)
        {
            string navMsg = trimCsvField(columns[3]);
            if (!navMsg.empty() && (navMsg.back() == '\r' || navMsg.back() == '\n'))
            {
                navMsg.pop_back();
            }
            string navBits = hexToBinary(navMsg);
            if (navBits.size() < 178)
            {
                continue;
            }
            osnma40 = navBits.substr(138, 40);
        }
        else
        {
            continue;
        }

        if (!isBinaryString(osnma40, 40))
        {
            continue;
        }

        OsnmaRecord record;
        record.tow = tow;
        record.osnma40 = osnma40;
        records.push_back(record);
    }

    inf.close();
    return records;
}

map<int, string> makeOsnmaByTow(const vector<OsnmaRecord> &records)
{
    map<int, string> byTow;
    for (size_t i = 0; i < records.size(); i++)
    {
        byTow[records[i].tow] = records[i].osnma40;
    }
    return byTow;
}

} // namespace

void show1(vector<page_struct> navPage)
{
    for(int i=0;i<navPage.size();i++)
    {
        cout<<navPage[i].tow<<","<<cur_wn<<","<<navPage[i].prn<<","<<navPage[i].navMsg<<endl;
    }
}
vector<KeyItem> getKeyItem(int prn,string PATH)
{
    (void)prn;
    vector<OsnmaRecord> osnma_r = readOsnmaKeychain(PATH);
    
    
    vector<KeyItem> teslaKey;
    string mack;
  
    for(int i=0;i+30<static_cast<int>(osnma_r.size());i++)
    {
    	
        int tow=osnma_r[i].tow;
        int firstTow=getFirstTow(tow);
        // if(is60page(firstTow))
        // {   
        //     // cout<<"Not the first tow in subframe: "<<tow<<endl;
        //    continue;
        // }
        if(getPageNumInSub(tow)>=10&&getPageNumInSub(tow)<=14) //Starting from here, extract tags for the next five pages
        {
            string temp=binaryToHexBitset(osnma_r[i].osnma40.substr(8,32));
            
            mack.append(temp);
        }
        if(getPageNumInSub(tow)==14)
        {   
            KeyItem temp;
            temp.key=mack.substr(4,32);
            //cout<<"key "<<temp.key<<endl;
            temp.firstTow=firstTow;
            teslaKey.push_back(temp);
            mack.clear();
            if(firstTow>573781) break;
        }
    }
      
    return teslaKey;
}
void  fillOSNMA(int prn,vector<page_struct> &cur_page,string PATH)
{
    (void)prn;
    vector<OsnmaRecord> osnma_r = readOsnmaKeychain(PATH);
    map<int, string> osnmaByTow = makeOsnmaByTow(osnma_r);
    //cout<<"osnma_r----------"<<endl;
    //show1(osnma_r);
    //cout<<"cur_page----------"<<endl;
   // show1(cur_page);

    //get the osnma data from osnam_r and fill it in osnma_c interesting
    int matchedRows = 0;
    int zeroFilledRows = 0;
    for(int i=0;i<static_cast<int>(cur_page.size());i++) 
    {
        int tow = static_cast<int>(cur_page[i].tow);
        string nav_mas_bing_c=hexToBinary(cur_page[i].navMsg);
        if(nav_mas_bing_c.size() < 178)
        {
            continue;
        }

        string osnma_bin = kZeroOsnma40;
        map<int, string>::const_iterator osnmaIt = osnmaByTow.find(tow);
        if(osnmaIt != osnmaByTow.end())
        {
            osnma_bin = osnmaIt->second;
            matchedRows++;
        }
        else
        {
            zeroFilledRows++;
        }

        nav_mas_bing_c.replace(138,40,osnma_bin);
        
        string nav_mes_hex=binaryToHexBitset(nav_mas_bing_c);
        //exit(1);
        cur_page[i].navMsg=nav_mes_hex;
        
        
    }
    if(zeroFilledRows > 0)
    {
        cout << "[WARN] E" << (cur_page.empty() ? -1 : cur_page[0].prn)
             << ": " << zeroFilledRows
             << " OSNMA pages had no matching TOW and were zero-filled"
             << " (matched " << matchedRows << " rows)." << endl;
    }
}
