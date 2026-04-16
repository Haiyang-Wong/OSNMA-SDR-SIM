#include "../include/galileo-sdr.h"
int len_tow=20;
const string mapping = "0123456789ABCDEF";
const int len_gst=32;
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
vector<vector<string>> getNavData(vector<page_struct> &cur_page)
{   

    vector<vector<string>> words(30, vector<string>(5));

    string mack;

    int j=0;
    for (size_t i = 0; i < cur_page.size(); i++)
    {
        int tow=static_cast<int>(cur_page[i].tow);
        //cout<<"Tow"<<cur_page[i].tow<<endl;
        if(getPageNumInSub(tow)==10)
        {
            string navBin=hexToBinary(cur_page[i].navMsg);
            string s1=navBin.substr(8,106)+navBin.substr(122,14);
            //cout<<s1<<endl;
            words[j][0]=s1;
        }
        if(getPageNumInSub(tow)==0)
        {
            string navBin=hexToBinary(cur_page[i].navMsg);
            string s1=navBin.substr(8,106)+navBin.substr(122,14);
            //cout<<"Tow"<<cur_page[i].tow<<"j:"<<binaryToHexBitset(s1)<<endl;
            words[j][1]=s1;
        }
        if (getPageNumInSub(tow)==11)
        {
            string navBin=hexToBinary(cur_page[i].navMsg);
            string s1=navBin.substr(8,106)+navBin.substr(122,16);
            words[j][2]=s1;
        }
        if (getPageNumInSub(tow)==1)
        {
            string navBin=hexToBinary(cur_page[i].navMsg);
            string s1=navBin.substr(8,106)+navBin.substr(122,14);
            words[j][3]=s1;
        }
        if (getPageNumInSub(tow)==12)
        {
            string navBin=hexToBinary(cur_page[i].navMsg);
            string s1=navBin.substr(8,67);
            words[j][4]=s1;
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
    vector<vector<string>> navData = getNavData(cur_page);
    for(int i=0;i<numSub-1;i++){
    string se="";
    	for(int j=0;j<5;j++)
    	{
    	    se+=navData[i][j];
    	}
    	//cout<<"xianzai:"<<i<<" :"<<se<<endl;
    }
    vector<string> authData;
    vector<string> navDataBin(numSub);

    for(int i=0;i<numSub-1;i++){
    navDataBin[i] ="";
    	for(int j=0;j<5;j++)
    	{
    	    navDataBin[i]+=navData[i][j];
    	}
    	
    }
    //cout<<navDataBin<<endl;
    for(int i=0;i<numSub-1;i++)
    {
        string temp;
        string sss="";
        
        /*if(prn<10)
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
        */
        if(prn<16)
        {
        	//cout<<mapping[prn];
        	sss=string("0")+mapping[prn];
        	
        	//cout<<sss<<endl;
        }
        else if(prn<30){
        	sss=to_string(prn-6);
        }
        else
        {
        	sss=to_string(prn-12);
        }
        //cout<<"ssss:------getAuthDatagetAuthData---------"<<sss<<endl;
        
        string prn_a=hexToBinary(sss);
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
vector<vector<string>> getNavDataForADKD4(vector<page_struct> &cur_page) 
{
    // word 6:page 2 ,bit: 8~107. word 7:page 4,bit:88~114 122~138
    vector<vector<string>> navData; // store the nav data for each page
    vector<string> words(3);
    // cout<<"ADKD :"<<endl;
    for(int i=0;i<cur_page.size();i++)
    {
        
        int tow=static_cast<int>(cur_page[i].tow);
        int firstTow=getFirstTow(tow);
        // cout<<"tow: "<<tow<<endl;   
        // cout<<"First tow: "<<firstTow<<endl;
        if(is60page(firstTow))
        {  
            // cout<<"Not the first tow in subframe: "<<tow<<endl;
           continue;
        }
        if(getPageNumInSub(tow)==2)
        {   
            
            string navBin=hexToBinary(cur_page[i].navMsg);
            string temp=navBin.substr(8,99);
            words[0]=temp;
            // cout<<"wO:"<<words[0]<<endl;
        }
        if(getPageNumInSub(tow)==4)
        {   
            // cout<<"tow 4 :"<<tow<<endl;
            string navBin=hexToBinary(cur_page[i].navMsg);
            string s1=navBin.substr(88,26);
            string s2=navBin.substr(122,16);
            words[1]=s1+s2;
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
    // 这里直接从nav_page中提取
    vector<string> authData;
    vector<gst> GST=getGst(cur_page);
    vector<vector<string>> navData = getNavDataForADKD4(cur_page);
    // cout<<"Nav data for ADKD4: "<<navData.size()<<endl;
    // cout<<"GST size: "<<GST.size()<<endl;
    
    string nav_dates=navData[0][0]+navData[0][1];
    // cout<<"Nav data for ADKD4: "<<endl;
    for(int i=0;i<navData.size()*2-1;i++)
    {   
        string temp;
        string sss="";
        if(prn<16)
        {
        	//cout<<mapping[prn];
        	sss=string("0")+mapping[prn];
        	
        	//cout<<sss<<endl;
        }
        else if(prn<30){
        	sss=to_string(prn-6);
        }
        else
        {
        	sss=to_string(prn-12);
        }
        string prn_a=hexToBinary(sss);
        string prn_d=hexToBinary(sss);
        string tagcont=hexToBinary("03");
        string nav_status;
        if(dic==33){
        	nav_status="01";
        }
        else{
        	nav_status="10";
        }
        bitset<len_gst> bitArrays((stoi(GST[i].wn) << len_tow | stoi(GST[i].tow)));
        string gst_bin=bitArrays.to_string();
        // cout<<"GST for ADKD4: "<<GST[i].tow<<""<<endl;
        temp =prn_a+prn_d+gst_bin+tagcont+nav_status+nav_dates;
        //cout<<"Auth data for ADKD4: "<<temp<<endl;
        
        authData.push_back(temp);
    }
    authData.erase(authData.begin());
    return authData;
}
vector<string> computeTagADKD4(vector<string> key,vector<string> authData) //
{
    vector<string> tagList;
    for(int i=0;i<authData.size()&&i<key.size();i++)
    {	
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
    
    for(int i=0;i+1<numSub*30;i++)
    {
        
        int tow=static_cast<int>(nav_page[i].tow);
        int firstTow=getFirstTow(tow);
        int conts=0;
        if(getPageNumInSub(tow)==3&&is60page(firstTow)&&firstTow!=571501)
        {   
           
            string tag=tagList.front();
            tag=hexToBinary(tag.substr(0,10)); 
            tagList.erase(tagList.begin());
            string navBin1=hexToBinary(nav_page[i].navMsg);
            string navBin2=hexToBinary(nav_page[i+1].navMsg);
            string tag1=tag.substr(0,16);
            string tag2=tag.substr(16);
            // cout<<tag1<<","<<tag2<<endl;

            navBin1.replace(146+16,16,tag1);
            navBin2.replace(146,24,tag2);
            nav_page[i].navMsg=binaryToHexBitset(navBin1);
            nav_page[i+1].navMsg=binaryToHexBitset(navBin2);
            
            i++;
 	    if(tagList.empty()) return;
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
vector<string> getAuthDatADKD12(const vector<page_struct> &cur_page,int prn) //返回一个二进制的字符串
{
    //Tag0所需要的东西，self.prn_a + self.gst_subframe.bitarray + BitArray(uint=self.ctr,length=8) + self.nma_status + self.nav_data.nav_data_stream
    vector<gst> GST=getGstADKD12(cur_page);
    vector<string> navData = getNavDataADKD12(cur_page);
    vector<string> authData;
    string navDataBin;
    for(int i=0;i<navData.size();i++)
    {   
        navDataBin+=navData[i];
    }
    //cout<<navDataBin<<endl;
    for(int i=0;i<10;i++)
    {
        string temp;
        string sss="";
        if(prn<16)
        {
        	//cout<<mapping[prn];
        	sss=string("0")+mapping[prn];
        	
        	//cout<<sss<<endl;
        }
        
        else
        {
        	sss=to_string(prn-12);
        }
        string prn_a=hexToBinary(sss);
        string prn_d=hexToBinary(sss);
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
vector<string> computeTagADKD12(vector<string> key,vector<string> authData) //
{
    vector<string> tagList;
    for(int i=0;i<key.size();i++)
    {
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

    for(int i=0;i<nav_page.size();i++)
    {
        
        int tow=static_cast<int>(nav_page[i].tow);
        int firstTow=getFirstTow(tow);
        if(is60page(firstTow))
        {   
          
           if(getPageNumInSub(tow)==7&&firstTow!=571501)
            {   
                //cout<<"zhe li shi:is60page"<<endl;
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
                 if(tagList.empty()) return;
                i++;
 
            }
        }
        else
        {
        
	   if(getPageNumInSub(tow)==5&&firstTow!=571501)
            {   
            //cout<<"zhe li not shi:is60page"<<endl;
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
                if(tagList.empty()) return;
                i++;
 
            }

        }
        
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
