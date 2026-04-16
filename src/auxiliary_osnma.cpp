#include "../include/galileo-sdr.h"
#include<string>
void show1(vector<page_struct> navPage)
{
    for(int i=0;i<navPage.size();i++)
    {
        cout<<navPage[i].tow<<","<<cur_wn<<","<<navPage[i].prn<<","<<navPage[i].navMsg<<endl;
    }
}
vector<KeyItem> getKeyItem(int prn,string PATH)
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
    
    
    vector<KeyItem> teslaKey;
    string mack;
  
    for(int i=0;i<osnma_r.size()-30;i++)
    {
    	
        int tow=static_cast<int>(osnma_r[i].tow);
        int firstTow=getFirstTow(tow);
        // if(is60page(firstTow))
        // {   
        //     // cout<<"Not the first tow in subframe: "<<tow<<endl;
        //    continue;
        // }
        if(getPageNumInSub(tow)>=10&&getPageNumInSub(tow)<=14) //Starting from here, extract tags for the next five pages
        {
            //cout<<osnma_r[i].navMsg<<endl;
            string temp=getMackInPage(osnma_r[i].navMsg);
            
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
    //cout<<"osnma_r----------"<<endl;
    //show1(osnma_r);
    //cout<<"cur_page----------"<<endl;
   // show1(cur_page);

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
