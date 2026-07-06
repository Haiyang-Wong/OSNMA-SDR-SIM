#include "../include/galileo-sdr.h"
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
void fill_nav_page(galtime_t g, ephem_t *eph,int *odd_page,int *even_page)
{
        page_struct a;
	a.tow=g.sec;
	if(start_time>g.sec)
	{
		start_time=g.sec;
	}
        a.wn=static_cast<int>(cur_wn);
        a.prn=eph->svid;
        //cout<<a.prn<<"++++++++++++++++++++++++"<<endl;
        ostringstream oss;

        for (int i = 0; i < 120; ++i) 
        {
            oss << even_page[i];
        }
        for (int i = 0; i < 120; ++i) 
        {
            oss << odd_page[i];
        }
        string ss = oss.str();
        a.navMsg=binaryToHexBitset(ss);
        nav_page.push_back(a);
	//if(a.prn==2)
        //{
            
           // cout<<a.tow<<"---------------"<<endl;
            //nav_page_2.push_back(a);
       // }
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

