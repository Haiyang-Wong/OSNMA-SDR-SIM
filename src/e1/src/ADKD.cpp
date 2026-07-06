#include "../include/galileo-sdr.h"
int len_tow=20;
const string mapping = "0123456789ABCDEF";
const int len_gst=32;

static bool isAuthPageVisible(int prn, int tow)
{
    if (tow_visible_prn_sets.empty())
        return true;

    for (size_t i = 0; i < tow_visible_prn_sets.size(); i++)
    {
        int rangeStart = tow_visible_prn_sets[i].first;
        int rangeEnd = (i + 1 < tow_visible_prn_sets.size()) ? tow_visible_prn_sets[i + 1].first : INT_MAX;
        if (tow >= rangeStart && tow < rangeEnd)
        {
            const vector<int> &visiblePrns = tow_visible_prn_sets[i].second;
            return find(visiblePrns.begin(), visiblePrns.end(), prn) != visiblePrns.end();
        }
    }

    return false;
}

void addWord10(vector<page_struct> &cur_page)
{
    for (size_t i = 0; i < cur_page.size(); i++)
    {
        int cur_tow=static_cast<int>(cur_page[i].tow);
        
        if(getPageNumInSub(cur_tow)==4&&!is60page(getFirstTow(cur_tow)))
        {
            string msg=hexToBinary(cur_page[i].navMsg);
            //cout<<cur_tow<<msg<<endl;
            msg.replace(2,6,"001010");
            //cout<<msg<<endl;

            msg=binaryToHexBitset(msg);
            cur_page[i].navMsg=msg;

        }
    }
}
/*--------------------ADKD 0-----------------*/
vector<vector<string>> getNavData(int prn, vector<page_struct> &cur_page)
{   

    vector<vector<string>> words(max(numSub, 1), vector<string>(5));

    string mack;

    int j=0;
    for (size_t i = 0; i < cur_page.size(); i++)
    {
        int tow=static_cast<int>(cur_page[i].tow);
        bool pageVisible = isAuthPageVisible(prn, tow);
        //cout<<"Tow"<<cur_page[i].tow<<endl;
        if(getPageNumInSub(tow)==10)
        {
            if(j >= static_cast<int>(words.size())) break;
            if(pageVisible)
            {
                string navBin=hexToBinary(cur_page[i].navMsg);
                string s1=navBin.substr(8,106)+navBin.substr(122,14);
                words[j][0]=s1;
            }
        }
        if(getPageNumInSub(tow)==0)
        {
            if(j >= static_cast<int>(words.size())) break;
            if(pageVisible)
            {
                string navBin=hexToBinary(cur_page[i].navMsg);
                string s1=navBin.substr(8,106)+navBin.substr(122,14);
                words[j][1]=s1;
            }
        }
        if (getPageNumInSub(tow)==11)
        {
            if(j >= static_cast<int>(words.size())) break;
            if(pageVisible)
            {
                string navBin=hexToBinary(cur_page[i].navMsg);
                string s1=navBin.substr(8,106)+navBin.substr(122,16);
                words[j][2]=s1;
            }
        }
        if (getPageNumInSub(tow)==1)
        {
            if(j >= static_cast<int>(words.size())) break;
            if(pageVisible)
            {
                string navBin=hexToBinary(cur_page[i].navMsg);
                string s1=navBin.substr(8,106)+navBin.substr(122,14);
                words[j][3]=s1;
            }
        }
        if (getPageNumInSub(tow)==12)
        {
            if(j >= static_cast<int>(words.size())) break;
            if(pageVisible)
            {
                string navBin=hexToBinary(cur_page[i].navMsg);
                string s1=navBin.substr(8,67);
                words[j][4]=s1;
            }
            j++;
        }

    }

    return words;
}

vector<string> getAuthData(int prn, vector<page_struct> &cur_page) //Should return a binary string
{
    //What Tag0 needs，self.prn_a + self.gst_subframe.bitarray + BitArray(uint=self.ctr,length=8) + self.nma_status + self.nav_data.nav_data_stream
    vector<gst> GST=getGst(cur_page);
    //  for(int i=0;i<GST.size();i++)
    // {
    //     cout<<GST[i].tow<<","<<GST[i].wn<<endl;
    // }
    vector<vector<string>> navData = getNavData(prn, cur_page);
    size_t usableSubframes = min(static_cast<size_t>(max(numSub - 1, 0)), min(navData.size(), GST.size()));
    for(size_t i=0;i<usableSubframes;i++){
    string se="";
    	for(int j=0;j<5;j++)
    	{
    	    se+=navData[i][j];
    	}
    	//cout<<"xianzai:"<<i<<" :"<<se<<endl;
    }
    vector<string> authData;
    vector<string> navDataBin(usableSubframes);

    for(size_t i=0;i<usableSubframes;i++){
    navDataBin[i] ="";
    	for(int j=0;j<5;j++)
    	{
    	    navDataBin[i]+=navData[i][j];
    	}
    	
    }
    //cout<<navDataBin<<endl;
    for(size_t i=0;i<usableSubframes;i++)
    {
        if(navDataBin[i].size() != 549)
        {
            authData.push_back("");
            continue;
        }
        string temp;
        string prn_a=bitset<8>(prn).to_string();
        string bitAr=hexToBinary("01");
        string nav_status;
        if(dic==33){
        	nav_status="01";
        }
        else{
        	nav_status="10";
        }
        bitset<len_gst> bitArray((stoi(GST[i].wn) << 20 | stoi(GST[i].tow)));
        string gst_bin=bitArray.to_string();
        // cout<<gst_bin<<endl;
        temp=prn_a+gst_bin+bitAr+nav_status+navDataBin[i];
        authData.push_back(temp);
        //cout<<temp<<endl;
    }

    return authData;
}
/*--------------------ADKD 4-----------------*/
vector<string> keyforADKD4(vector<KeyItem> key,int start)
{
    vector<string> sss;
    //cout<<"key.size :"<<key.size()<<key[1].key<<endl;
    //cout<<"numSub :"<<numSub<<endl;
    //KeyItem t={"adwawd",123};
    //if(!is60page(getFirstTow(key[0].firstTow)))
   // {
	//    key.insert(key.begin(),t);
   // }

    if(!is60page(getFirstTow(start)))
    {
    	//start+=30;
    }
    
    for (size_t i = 0; i < key.size(); i++)
    {
        if(key[i].firstTow==start)
        {
            i+=3;
            if(!is60page(getFirstTow(start)))
	    {
	    	i-=1;
	    }
    
            vector<string> sss;
            for(int count=0;count<numSub-1&&i < key.size();count++)
            {
            	//cout<<"key :"<<key[i].key<<endl;
                sss.push_back(key[i].key);
                i+=1;
            }
            return sss;
        }
    }
    //scout<<"zhe li :"<<key.size()<<endl;
    return sss;
    
}
vector<vector<string>> getNavDataForADKD4(int prn, vector<page_struct> &cur_page) 
{
    // word 6:page 2 ,bit: 8~107. word 7:page 4,bit:88~114 122~138
    vector<vector<string>> navData; // store the nav data for each page
    vector<string> words(3);
    // cout<<"ADKD :"<<endl;
    for(int i=0;i<cur_page.size();i++)
    {
        
        int tow=static_cast<int>(cur_page[i].tow);
        int firstTow=getFirstTow(tow);
        bool pageVisible = isAuthPageVisible(prn, tow);
        // cout<<"tow: "<<tow<<endl;   
        // cout<<"First tow: "<<firstTow<<endl;
        if(is60page(firstTow))
        {  
            // cout<<"Not the first tow in subframe: "<<tow<<endl;
           continue;
        }
        if(getPageNumInSub(tow)==2)
        {   
            if(pageVisible)
            {
                string navBin=hexToBinary(cur_page[i].navMsg);
                string temp=navBin.substr(8,99);
                words[0]=temp;
            }
            // cout<<"wO:"<<words[0]<<endl;
        }
        if(getPageNumInSub(tow)==4)
        {   
            if(pageVisible)
            {
                string navBin=hexToBinary(cur_page[i].navMsg);
                string s1=navBin.substr(88,26);
                string s2=navBin.substr(122,16);
                words[1]=s1+s2;
            }
            // cout<<"w1:"<<words[1]<<" sss"<<firstTow<<endl;
            
            words[2]=to_string(firstTow);
            // cout<<"w2:"<<words[2]<<endl;
            navData.push_back(words);
            
            words = vector<string>(3);
            // cout<<"++++++++++++++++++++++++++++++++++"<<endl;
        }
    }
    
    return navData;
}
vector<string> getAuthDataForADKD4(int prn, vector<page_struct> &cur_page)
{
    // Extract directly from nav_page.
    vector<string> authData;
    vector<gst> GST=getGst(cur_page);
    vector<vector<string>> navData = getNavDataForADKD4(prn, cur_page);
    // cout<<"Nav data for ADKD4: "<<navData.size()<<endl;
    // cout<<"GST size: "<<GST.size()<<endl;
    if(navData.empty() || GST.empty())
    {
        return authData;
    }

    // Each ADKD4 self-auth tag is placed in a 60s subframe and authenticates
    // the WT6/WT10 data from the preceding non-60s subframe.
    for(size_t i=0;i<navData.size();i++)
    {   
        string temp;
        string prn_a=bitset<8>(prn).to_string();
        string prn_d=bitset<8>(prn).to_string();
        string tagcont=hexToBinary("03");
        string nav_status;
        if(dic==33){
        	nav_status="01";
        }
        else{
        	nav_status="10";
        }
        if(navData[i][0].size() != 99 || navData[i][1].size() != 42)
        {
            authData.push_back(string(181, '0'));
            continue;
        }
        int tagTow = stoi(navData[i][2]) + 29;
        bitset<len_gst> bitArrays((cur_page[0].wn << len_tow | tagTow));
        string gst_bin=bitArrays.to_string();
        // cout<<"GST for ADKD4: "<<GST[i].tow<<""<<endl;
        temp =prn_a+prn_d+gst_bin+tagcont+nav_status+navData[i][0]+navData[i][1];
        //cout<<"Auth data for ADKD4: "<<temp<<endl;
        
        authData.push_back(temp);
    }
    return authData;
}
vector<string> computeTagADKD4(vector<string> key,vector<string> authData) //
{
    vector<string> tagList;
    for(int i=0;i<authData.size()&&i<key.size();i++)
    {	
        if(authData[i].empty() || key[i].empty())
        {
            tagList.push_back("0000000000");
            continue;
        }
    	//cout<<"key:"<<key[i]<<" -----  auth_data:"<<authData[i]<<endl;
        vector<uint8_t>k1=hexstr_to_bytes(key[i]);

        vector<uint8_t>d1=bitstring_to_bytes(authData[i]);

        auto mac= hmac_sha256(k1,d1);

        string temp=macToString(mac);
        //cout<<"key: "<<key[i]<<"  temp:"<<temp<<endl;
        temp=temp.substr(0,10);
        
        tagList.push_back(temp);
    }
    //cout<<"huillm"<<endl;
    return tagList;
}
void fillTagADKD4(vector<page_struct> &nav_page,vector<string> &tagList)
{
    vector<string> pageKeys = getKey(nav_page);
    for(int i=0;i+1<static_cast<int>(nav_page.size());i++)
    {
        
        int tow=static_cast<int>(nav_page[i].tow);
        int firstTow=getFirstTow(tow);
        int conts=0;
        if(getPageNumInSub(tow)==3&&is60page(firstTow)&&firstTow>static_cast<int>(start_time))
        {   
            string tag;
            bool dummyTag = false;
            if(tagList.empty())
            {
                int subframeIndex = (firstTow - static_cast<int>(start_time)) / 30;
                int keyIndex = subframeIndex + 1;
                if(keyIndex < 0 || keyIndex >= static_cast<int>(pageKeys.size()))
                {
                    return;
                }

                int prn = nav_page[i].prn;
                string prnBits = bitset<8>(prn).to_string();
                bitset<len_gst> gstBits((nav_page[i].wn << len_tow) | (firstTow - 1));
                string navStatus = dic==33 ? "01" : "10";
                string authData = prnBits + prnBits + gstBits.to_string() + hexToBinary("03") + navStatus + string(141, '0');
                vector<uint8_t> k1=hexstr_to_bytes(pageKeys[keyIndex]);
                vector<uint8_t> d1=bitstring_to_bytes(authData);
                tag=macToString(hmac_sha256(k1,d1)).substr(0,10);
                dummyTag = true;
            }
            else
            {
                tag=tagList.front();
                tagList.erase(tagList.begin());
            }
            tag=hexToBinary(tag.substr(0,10)); 
            string navBin1=hexToBinary(nav_page[i].navMsg);
            string navBin2=hexToBinary(nav_page[i+1].navMsg);
            string tag1=tag.substr(0,16);
            string tag2=tag.substr(16);
            // cout<<tag1<<","<<tag2<<endl;

            navBin1.replace(146+16,16,tag1);
            navBin2.replace(146,24,tag2);
	            nav_page[i].navMsg=binaryToHexBitset(navBin1);
	            nav_page[i+1].navMsg=binaryToHexBitset(navBin2);
		            recordGeneratedOsnmaTag(static_cast<int>(nav_page[i].tow),
		                                    static_cast<int>(nav_page[i].tow),
		                                    nav_page[i].prn,
		                                    nav_page[i].prn,
		                                    4,
		                                    false,
		                                    "self_adkd4",
		                                    1);
	            if(dummyTag && i+2 < static_cast<int>(nav_page.size()))
	            {
	                int prn = nav_page[i].prn;
	                string infoBits = bitset<8>(prn).to_string() + "0100" + "0000";
                string navBinInfo1=hexToBinary(nav_page[i+1].navMsg);
                string navBinInfo2=hexToBinary(nav_page[i+2].navMsg);
                navBinInfo1.replace(170,8,infoBits.substr(0,8));
                navBinInfo2.replace(146,8,infoBits.substr(8,8));
                nav_page[i+1].navMsg=binaryToHexBitset(navBinInfo1);
                nav_page[i+2].navMsg=binaryToHexBitset(navBinInfo2);
            }
            
            i++;
        }
        else
        {
            
        }
    }
    // cout<<"sssssssssssssssss"<<endl;
    return;
}
/*--------------------ADKD 12-----------------*/
vector<gst> getGstADKD12(const vector<page_struct> &cur_page)
{
  
    vector<gst> GST;
    for (size_t i = 0; i < cur_page.size(); i++)
    {
        int tow = static_cast<int>(cur_page[i].tow);
        if((tow>static_cast<int>(start_time)+28)&&(tow-(static_cast<int>(start_time)+30))%30==0)
        {
            
            if(1)
            {
                // cout<<"Tow of GST: "<<tow<<endl;
                gst temp;
                // cout<<tow<<endl;
                temp.tow=to_string(tow-1);
                //cout<<columns[0]<<endl;
                temp.wn=to_string(cur_page[i].wn);
                GST.push_back(temp);
                bitset<len_gst> bitArray((stoi(temp.wn) << len_tow | stoi(temp.tow)));
                // cout<<bitArray.to_string()<<endl;
            }
        }
    }
    return GST;

}
vector<string> getNavDataADKD12(const vector<page_struct> &cur_page)
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

vector<string> getNavDataADKD12BySubframe(const vector<page_struct> &cur_page, int prn)
{
    vector<string> navDataList;
    vector<string> words(5);

    for (size_t i = 0; i < cur_page.size(); i++)
    {
        int tow=static_cast<int>(cur_page[i].tow);
        bool pageVisible = isAuthPageVisible(prn, tow);
        if(getPageNumInSub(tow)==10)
        {
            if(pageVisible)
            {
                string navBin=hexToBinary(cur_page[i].navMsg);
                words[0]=navBin.substr(8,106)+navBin.substr(122,14);
            }
        }
        if(getPageNumInSub(tow)==0)
        {
            if(pageVisible)
            {
                string navBin=hexToBinary(cur_page[i].navMsg);
                words[1]=navBin.substr(8,106)+navBin.substr(122,14);
            }
        }
        if (getPageNumInSub(tow)==11)
        {
            if(pageVisible)
            {
                string navBin=hexToBinary(cur_page[i].navMsg);
                words[2]=navBin.substr(8,106)+navBin.substr(122,16);
            }
        }
        if (getPageNumInSub(tow)==1)
        {
            if(pageVisible)
            {
                string navBin=hexToBinary(cur_page[i].navMsg);
                words[3]=navBin.substr(8,106)+navBin.substr(122,14);
            }
        }
        if (getPageNumInSub(tow)==12)
        {
            if(pageVisible)
            {
                string navBin=hexToBinary(cur_page[i].navMsg);
                words[4]=navBin.substr(8,67);
            }
            string navDataBin;
            for(size_t j=0;j<words.size();j++)
            {
                navDataBin+=words[j];
            }
            navDataList.push_back(navDataBin);
            words = vector<string>(5);
        }
    }

    return navDataList;
}

vector<string> getAuthDatADKD12(const vector<page_struct> &cur_page,int prn) // Return a binary string.
{
    // Tag0 fields: self.prn_a + self.gst_subframe.bitarray + BitArray(uint=self.ctr,length=8) + self.nma_status + self.nav_data.nav_data_stream.
    vector<gst> GST=getGstADKD12(cur_page);
    vector<string> navDataBySubframe = getNavDataADKD12BySubframe(cur_page, prn);
    vector<string> authData;
    if(navDataBySubframe.empty() || GST.empty())
    {
        return authData;
    }
    //cout<<navDataBin<<endl;
    for(size_t i=0;i<GST.size();i++)
    {
        size_t navDataIndex = min(i, navDataBySubframe.size() - 1);
        string navDataBin = navDataBySubframe[navDataIndex];
        if(navDataBin.size() != 549)
        {
            authData.push_back("");
            continue;
        }
        string temp;
        string prn_a=bitset<8>(prn).to_string();
        string prn_d=bitset<8>(prn).to_string();
        string bitAr="";
        if((stoi(GST[i].tow)-277200)%60==0)
        {
            bitAr=hexToBinary("05");
        }
        else
        {
            bitAr=hexToBinary("04");
        }

        // cout<<"Nav data length: "<<GST[i].tow<<endl;
        string nav_status;
        if(dic==33){
        	nav_status="01";
        }
        else{
        	nav_status="10";
        }
        bitset<len_gst> bitArray((stoi(GST[i].wn) << len_tow | stoi(GST[i].tow)));
        // cout<<"GST for ADKD12: "<<bitArray<<""<<endl;
        string gst_bin=bitArray.to_string();

        temp=prn_a+prn_d+gst_bin+bitAr+nav_status+navDataBin;
        authData.push_back(temp);
        // cout<<temp<<endl;
    }

    return authData;
}

vector<string> getAuthDataADKD12ForCounter(const vector<page_struct> &cur_page, int prn, int counter)
{
    vector<gst> GST=getGstADKD12(cur_page);
    vector<string> navDataBySubframe = getNavDataADKD12BySubframe(cur_page, prn);
    vector<string> authData;
    if(navDataBySubframe.empty() || GST.empty() || counter < 0 || counter > 255)
    {
        return authData;
    }

    for(size_t i=0;i<GST.size();i++)
    {
        size_t navDataIndex = min(i, navDataBySubframe.size() - 1);
        string navDataBin = navDataBySubframe[navDataIndex];
        if(navDataBin.size() != 549)
        {
            authData.push_back("");
            continue;
        }
        string prn_a=bitset<8>(prn).to_string();
        string prn_d=bitset<8>(prn).to_string();
        string nav_status = dic==33 ? "01" : "10";
        bitset<len_gst> bitArray((stoi(GST[i].wn) << len_tow | stoi(GST[i].tow)));
        authData.push_back(prn_a + prn_d + bitArray.to_string() + bitset<8>(counter).to_string() + nav_status + navDataBin);
    }

    return authData;
}

vector<string> computeTagADKD12(vector<string> key,vector<string> authData) //
{
    vector<string> tagList;
    for(size_t i=0;i<key.size() && i<authData.size();i++)
    {
        if(authData[i].empty())
        {
            tagList.push_back("0000000000");
            continue;
        }
        vector<uint8_t>k1=hexstr_to_bytes(key[i]);
        // 
        vector<uint8_t>d1=bitstring_to_bytes(authData[i]);

        auto mac= hmac_sha256(k1,d1);

        string temp=macToString(mac);
        //cout<<"Key for ADKD12: "<<key[i]<<"  auth:"<<authData[i]<<endl;
        // temp=temp.substr(0,10);
        //cout<<temp<<endl;
        tagList.push_back(temp);
    }
    return tagList;
}
void fillTagADKD12(vector<page_struct> &nav_page,vector<string> tagList)
{
    if(tagList.empty())
    {
        return;
    }

    for(int i=0;i<nav_page.size();i++)
    {
        
        int tow=static_cast<int>(nav_page[i].tow);
        int firstTow=getFirstTow(tow);
        if(is60page(firstTow))
        {   
          
           if(getPageNumInSub(tow)==7&&firstTow>static_cast<int>(start_time))
            {   
                //cout<<"zhe li shi:is60page"<<endl;
                if(tagList.empty() || i+1 >= static_cast<int>(nav_page.size()))
                {
                    return;
                }
                string tag=tagList.front();
                //cout<<"cur tag:"<<tag<<"   cur tow"<<tow<<endl;
                tag=hexToBinary(tag.substr(0,10)); 
                
                tagList.erase(tagList.begin());
                
                string navBin1=hexToBinary(nav_page[i].navMsg);
                string navBin2=hexToBinary(nav_page[i+1].navMsg);
                string tag1=tag.substr(0,32);
                string tag2=tag.substr(32);
                // cout<<tag1<<","<<tag2<<endl;

	                navBin1.replace(146,32,tag1);
	                navBin2.replace(146,8,tag2);
	                nav_page[i].navMsg=binaryToHexBitset(navBin1);
	                nav_page[i+1].navMsg=binaryToHexBitset(navBin2);
		                recordGeneratedOsnmaTag(static_cast<int>(nav_page[i].tow),
		                                        static_cast<int>(nav_page[i].tow),
		                                        nav_page[i].prn,
		                                        nav_page[i].prn,
		                                        12,
		                                        false,
		                                        "self_adkd12",
		                                        0);
	                 if(tagList.empty()) return;
	                i++;
 
            }
        }
        else
        {
        
	   if(getPageNumInSub(tow)==5&&firstTow>static_cast<int>(start_time))
            {   
            //cout<<"zhe li not shi:is60page"<<endl;
                if(tagList.empty() || i+1 >= static_cast<int>(nav_page.size()))
                {
                    return;
                }
                string tag=tagList.front();
                //cout<<"cur tag:"<<tag<<"   cur tow"<<tow<<endl;
                tag=hexToBinary(tag.substr(0,10)); 
                
                tagList.erase(tagList.begin());
                
                string navBin1=hexToBinary(nav_page[i].navMsg);
                string navBin2=hexToBinary(nav_page[i+1].navMsg);
                string tag1=tag.substr(0,24);
                string tag2=tag.substr(24);
                // cout<<tag1<<","<<tag2<<endl;

	                navBin1.replace(146+8,24,tag1);
	                navBin2.replace(146,16,tag2);
	                nav_page[i].navMsg=binaryToHexBitset(navBin1);
	                nav_page[i+1].navMsg=binaryToHexBitset(navBin2);
		                recordGeneratedOsnmaTag(static_cast<int>(nav_page[i].tow),
		                                        static_cast<int>(nav_page[i].tow),
		                                        nav_page[i].prn,
		                                        nav_page[i].prn,
		                                        12,
		                                        false,
		                                        "self_adkd12",
		                                        0);
	                if(tagList.empty()) return;
	                i++;
 
            }

        }
        
    }
}

static bool readMackBits(const vector<page_struct> &nav_page, int subframeStart, int bitOffset, int bitLen, string &out)
{
    out.clear();
    int remaining = bitLen;
    int curOffset = bitOffset;
    while(remaining > 0)
    {
        int pageOffset = curOffset / 32;
        int bitInPage = curOffset % 32;
        int pageIndex = subframeStart + pageOffset;
        if(pageIndex < 0 || pageIndex >= static_cast<int>(nav_page.size()))
        {
            return false;
        }

        string navBin = hexToBinary(nav_page[pageIndex].navMsg);
        int navBit = 146 + bitInPage;
        int chunk = min(remaining, 32 - bitInPage);
        if(navBit < 0 || navBit + chunk > static_cast<int>(navBin.size()))
        {
            return false;
        }

        out += navBin.substr(navBit, chunk);
        curOffset += chunk;
        remaining -= chunk;
    }

    return true;
}

static bool writeMackBits(vector<page_struct> &nav_page, int subframeStart, int bitOffset, const string &bits)
{
    int remaining = static_cast<int>(bits.size());
    int curOffset = bitOffset;
    int written = 0;
    while(remaining > 0)
    {
        int pageOffset = curOffset / 32;
        int bitInPage = curOffset % 32;
        int pageIndex = subframeStart + pageOffset;
        if(pageIndex < 0 || pageIndex >= static_cast<int>(nav_page.size()))
        {
            return false;
        }

        string navBin = hexToBinary(nav_page[pageIndex].navMsg);
        int navBit = 146 + bitInPage;
        int chunk = min(remaining, 32 - bitInPage);
        if(navBit < 0 || navBit + chunk > static_cast<int>(navBin.size()))
        {
            return false;
        }

        navBin.replace(navBit, chunk, bits.substr(written, chunk));
        nav_page[pageIndex].navMsg=binaryToHexBitset(navBin);
        curOffset += chunk;
        written += chunk;
        remaining -= chunk;
    }

    return true;
}

void fillTagADKD12Counter(vector<page_struct> &nav_page, vector<string> tagList, int counter)
{
    if(tagList.empty() || counter < 2)
    {
        return;
    }

    int tagIndex = 0;
    int tagOffset = 56 + (counter - 2) * 56;
    int infoOffset = tagOffset + 40;

    for(int i=0;i<static_cast<int>(nav_page.size());i++)
    {
        int tow=static_cast<int>(nav_page[i].tow);
        int firstTow=getFirstTow(tow);
        if(getPageNumInSub(tow)!=0 || firstTow<=static_cast<int>(start_time))
        {
            continue;
        }
        if(tagIndex >= static_cast<int>(tagList.size()))
        {
            return;
        }

        string infoBits;
        if(!readMackBits(nav_page, i, infoOffset, 16, infoBits))
        {
            tagIndex++;
            continue;
        }

        string expectedInfo = bitset<8>(nav_page[i].prn).to_string() + "1100" + "0001";
        if(infoBits != expectedInfo)
        {
            tagIndex++;
            continue;
        }

        string tag=hexToBinary(tagList[tagIndex].substr(0,10));
        if(writeMackBits(nav_page, i, tagOffset, tag))
        {
            int firstTagPage = i + tagOffset / 32;
            int writeTow = firstTagPage < static_cast<int>(nav_page.size()) ? static_cast<int>(nav_page[firstTagPage].tow) : tow;
            recordGeneratedOsnmaTag(writeTow,
                                    writeTow,
                                    nav_page[i].prn,
                                    nav_page[i].prn,
                                    12,
                                    false,
                                    "self_adkd12_ctr" + to_string(counter),
                                    counter);
        }
        tagIndex++;
    }
}
vector<string> keyforADKD12(vector<KeyItem> key,int start)
{
	vector<string> sss;
    for (size_t i = 0; i < key.size(); i++)
    {
        if(key[i].firstTow==start)
        {
            i+=12;
            vector<string> sss;
            for(int count=0;count<10&&i < key.size();count++)
            {
                sss.push_back(key[i].key);
                i+=1;
            }
            return sss;
        }
    }
    return sss;
    
}
