#include "newdy.h"

NewDy::NewDy()
{

}
//
void NewDy::D00(int rr, double probAlpha, double probBeta, int llnn, iMatrix &aa, int trr, int ln,int net)
{
    act=0;
    if(trr<rr)
    {
        double deactP=(probBeta/(probAlpha+probBeta));
        //cout<<deactP<<" "<<probAlpha<<endl;
        //cin.get();
        double sumOFf=0;
        random_device rd;
        mt19937 gen(rd());  // to seed mersenne twister.
        uniform_real_distribution<> dist1(0, 1);
        //uniform_int_distribution<> distall(0,llnn-1);
        ostringstream fn;
        vector<int> dy(llnn);
        vector<int> infects;
        //vector<int> rinfects;
        for(int i=0;i<llnn;++i)
        {
            dy[i]=1;
            infects.push_back(i);

        }
        /*
        for(int r=0;r<llnn;++r)
        {
            uniform_int_distribution<> rinf(0,rinfects.size()-1);
            int cal=rinf(gen);
            int rcal=rinfects[cal];
            infects.push_back(rcal);
            rinfects.shrink_to_fit();
        }
        */
        //cout<<infects.size()<<endl;
        //cin.get();
        //=============================================================
        fn<<"./data/"<<net<<"/"<<probAlpha<<"/infection_fraction_"<<rr<<"_"<<probAlpha<<"_netE"<<net<<"_landaE"<<ln<<".dat";
        ofstream out1(fn.str().c_str(),ios_base::binary);
        //
        int isum=0;
        for(int ii=0;ii<rr;++ii)
        {
            act=ii;
            isum=0;
            random_device rdinf;
            mt19937 geninf(rdinf());
            uniform_int_distribution<> infs(0,(infects.size()-1));
            uniform_real_distribution<> dist1l(0, 1);
            double actt=0;
            int infector=0,iactt=0;
            iactt=infs(geninf);
            infector=infects[iactt];
            actt=dist1(gen);
            //cout<<iactt<<" "<<infector<<endl;
            //cout<<actt<<endl;
            //cin.get();
            if(actt<=deactP)
            {
                infects.erase(infects.begin()+iactt);
                dy[infector]=0;
            }else
            {
                vector<int> suss;
                for(int m=0;m<aa[infector].size();++m)
                {
                    int mm=aa[infector][m];
                    if(dy[mm]==0)
                    {
                        suss.push_back(mm);
                    }
                }
                if(suss.size()>0)
                {
                    for(int j=0;j<suss.size();++j)
                    {
                      int s2i=suss[j];
                      double sp=dist1l(gen);
                      if(sp<=probAlpha)
                      {
                        dy[s2i]=1;
                        infects.push_back(s2i);
                      }
                    }
                }
                infects.erase(infects.begin()+iactt);
                dy[infector]=0;
                suss.clear();
                suss.shrink_to_fit();
            }
            //cout<<"a"<<endl;
            //---------------------------------------------
            //---------------------------------------------
            isum=infects.size();
            double f1=0,f2=0,f3=0;
            f1=isum;
            f2=llnn;
            f3=f1/f2;
            //cout<<f1<<" "<<f2<<" "<<f3<<endl;
            if(f3==0)
            {
                //cin.get();
                avgI=0;
                sumOFf=0;
                ii=rr;
                break;
            }
            int um=rr-trr;
            if(ii >= um)
                sumOFf=sumOFf+f3;
            //
            out1<< ii+1 <<' '<<f3<<endl;
        }
        avgI=sumOFf/trr;

    }else
    {
        cout<<"SIRS stoped!"<<endl<<"avrage loop is biger than loop try."<<endl;
    }
}
//
//
void NewDy::D0(int rr, double probAlpha, double probBeta, int llnn, iMatrix &aa, int trr, int ln, int net)
{
    act=0;
    if(trr<rr)
    {
        double deactP=(probBeta/(probAlpha+probBeta));
        //cout<<deactP<<" "<<probAlpha<<endl;
        //cin.get();
        double sumOFf=0;
        random_device rd;
        mt19937 gen(rd());  // to seed mersenne twister.
        uniform_real_distribution<> dist1(0, 1);
        //uniform_int_distribution<> distall(0,llnn-1);
        ostringstream fn;
        vector<int> dy(llnn);
        vector<int> infects;
        vector<int> rinfects;
        for(int i=0;i<llnn;++i)
        {
            dy[i]=1;
            rinfects.push_back(i);

        }
        for(int r=0;r<llnn;++r)
        {
            uniform_int_distribution<> rinf(0,rinfects.size()-1);
            int cal=rinf(gen);
            int rcal=rinfects[cal];
            infects.push_back(rcal);
            rinfects.shrink_to_fit();
        }
        //cout<<infects.size()<<" of "<<llnn<<endl;
        //cin.get();
        //=============================================================
        fn<<"./data/"<<net<<"/"<<probAlpha<<"/infection_fraction_"<<rr<<"_"<<probAlpha<<"_netE"<<net<<"_landaE"<<ln<<".dat";
        ofstream out1(fn.str().c_str(),ios_base::binary);
        //
        int isum=0;
        for(int ii=0;ii<rr;++ii)
        {
            act=ii;
            isum=0;
            //random_device rdinf;
            //mt19937 geninf(rdinf());
            //int dislimit=infects.size()-1;
            //cout<<dislimit<<endl;
            //cin.get();
            //uniform_int_distribution<> infs(0,dislimit);
            double actt=0;
            int infector=0,iactt=0;
            //iactt=infs(geninf);
            infector=infects[iactt];
            actt=dist1(gen);
            //cout<<iactt<<" "<<infector<<endl;
            //cout<<actt<<endl;
            //cin.get();
            if(actt<=deactP)
            {
                infects.erase(infects.begin()+iactt);
                dy[infector]=0;
            }else
            {
                int ics=aa[infector].size();
                for(int m=0;m<ics;++m)
                {
                    int mm=aa[infector][m];
                    if(dy[mm]==0)
                    {
                        double c2p=0;
                        c2p=dist1(gen);
                        if(c2p<=probAlpha)
                        {
                            dy[mm]=1;
                            infects.push_back(mm);
                        }
                    }
                }
                infects.erase(infects.begin()+iactt);
                dy[infector]=0;
            }
            //cout<<"a"<<endl;
            //---------------------------------------------
            //---------------------------------------------
            isum=infects.size();
            double f1=0,f2=0,f3=0;
            f1=isum;
            f2=llnn;
            f3=f1/f2;
            //cout<<f1<<" "<<f2<<" "<<f3<<endl;
            if(f3==0)
            {
                //cin.get();
                avgI=0;
                sumOFf=0;
                ii=rr;
                break;
            }
            int um=rr-trr;
            if(ii >= um)
                sumOFf=sumOFf+f3;
            //
            out1<< ii+1 <<' '<<f3<<endl;
        }
        avgI=sumOFf/trr;

    }else
    {
        cout<<"SIRS stoped!"<<endl<<"avrage loop is biger than loop try."<<endl;
    }
}
//
//
void NewDy::D01(int rr, double probAlpha, double probBeta, int llnn, iMatrix &aa, int trr, int ln,int net)
{
    //parallel update sis with infection abality for recently deactivated nodes
    act=0;
    if(trr<rr)
    {
        double sumOFf=0;
        random_device rd;
        mt19937 gen(rd());  // to seed mersenne twister.
        uniform_real_distribution<> dist1(0, 1);
        //uniform_int_distribution<> distall(0,llnn-1);
        ostringstream fn;
        vector<int> dy(llnn);
        vector<int> infects;
        vector<int> ninfects;
        vector<int> rinfects;
        //vector<int> resus;
        for(int i=0;i<llnn;++i)
        {
            dy[i]=1;
            rinfects.push_back(i);
        }
        for(int r=0;r<llnn;++r)
        {
            uniform_int_distribution<> rinf(0,rinfects.size()-1);
            int cal=rinf(gen);
            int rcal=rinfects[cal];
            infects.push_back(rcal);
            rinfects.shrink_to_fit();
        }
        //=============================================================
        fn<<"./data/"<<net<<"/"<<probAlpha<<"/infection_fraction_"<<rr<<"_"<<probAlpha<<"_netE"<<net<<"_landaE"<<ln<<".dat";
        ofstream out1(fn.str().c_str(),ios_base::binary);
        //
        cout<<"run..."<<endl;
        int isum=0;
        for(int ii=0;ii<rr;++ii)
        {
            //cout<<ii<<endl;
            act=ii;
            isum=0;
            int isize=infects.size();
            for(int i=0;i<isize;++i)
            {
                int infector=infects[0];
                double actt=dist1(gen);
                if(actt<=(probAlpha/(probAlpha+probBeta)))
                {
                    for(int m=0;m<aa[infector].size();++m)
                    {
                        int mm=aa[infector][m];
                        if(dy[mm]==0)
                        {
                            double c2p=0;
                            c2p=dist1(gen);
                            if(c2p<=probAlpha)
                            {
                                dy[mm]=1;
                                ninfects.push_back(mm);
                            }
                        }
                    }
                }
                infects.erase(infects.begin()+0);
                dy[infector]=0;
            }
            infects.clear();
            infects.shrink_to_fit();
            //
            //---------------------------------------------
            for(int j=0;j<ninfects.size();++j)
            {
                infects.push_back(0);
                infects[j]=ninfects[j];
            }
            ninfects.clear();
            ninfects.shrink_to_fit();
            //---------------------------------------------
            isum=infects.size();
            double f1=0,f2=0,f3=0;
            f1=isum;
            f2=llnn;
            f3=f1/f2;
            if(f3==0)
            {
                avgI=0;
                sumOFf=0;
                ii=rr;
                break;
            }
            int um=rr-trr;
            if(ii >= um)
                sumOFf=sumOFf+f3;
            //
            out1<< ii+1 <<' '<<f3<<endl;
        }
        avgI=sumOFf/trr;
    }else
    {
        cout<<"SIRS stoped!"<<endl<<"avrage loop is biger than loop try."<<endl;
    }
}
//
//
void NewDy::D1(int rr, double probAlpha, double probBeta, int llnn, iMatrix &aa, int trr, int ln,int net)
{
    //parallel update sis with infection abality for recently deactivated nodes
    act=0;
    if(trr<rr)
    {
        double sumOFf=0;
        random_device rd;
        mt19937 gen(rd());  // to seed mersenne twister.
        uniform_real_distribution<> dist1(0, 1);
        //uniform_int_distribution<> distall(0,llnn-1);
        ostringstream fn;
        vector<int> dy(llnn);
        vector<int> infects;
        vector<int> ninfects;
        //vector<int> resus;
        for(int i=0;i<llnn;++i)
        {
            dy[i]=1;
            infects.push_back(i);
        }
        //=============================================================
        fn<<"./data/"<<net<<"/"<<probAlpha<<"/infection_fraction_"<<rr<<"_"<<probAlpha<<"_netE"<<net<<"_landaE"<<ln<<".dat";
        ofstream out1(fn.str().c_str(),ios_base::binary);
        //
        int isum=0;
        for(int ii=0;ii<rr;++ii)
        {
            act=ii;
            isum=0;
            int isize=infects.size();
            for(int i=0;i<isize;++i)
            {
                int infector=infects[0];
                double actt=dist1(gen);
                if(actt<=(probAlpha/(probAlpha+probBeta)))
                {
                    for(int m=0;m<aa[infector].size();++m)
                    {
                        int mm=aa[infector][m];
                        if(dy[mm]==0)
                        {
                            double c2p=0;
                            c2p=dist1(gen);
                            if(c2p<=probAlpha)
                            {
                                dy[mm]=1;
                                ninfects.push_back(mm);
                            }
                        }
                    }
                }
                infects.erase(infects.begin()+0);
                dy[infector]=0;
            }
            infects.clear();
            infects.shrink_to_fit();
            //---------------------------------------------
            for(int j=0;j<ninfects.size();++j)
            {
                infects.push_back(0);
                infects[j]=ninfects[j];
            }
            ninfects.clear();
            ninfects.shrink_to_fit();
            //---------------------------------------------
            isum=infects.size();
            double f1=0,f2=0,f3=0;
            f1=isum;
            f2=llnn;
            f3=f1/f2;
            if(f3==0)
            {
                avgI=0;
                sumOFf=0;
                ii=rr;
                break;
            }
            int um=rr-trr;
            if(ii >= um)
                sumOFf=sumOFf+f3;
            //
            out1<< ii+1 <<' '<<f3<<endl;
        }
        avgI=sumOFf/trr;
    }else
    {
        cout<<"SIRS stoped!"<<endl<<"avrage loop is biger than loop try."<<endl;
    }
}
//
//
void NewDy::D2(int rr, double probAlpha, double probBeta, double teta, int rt, int llnn, iMatrix &aa, int trr, int ln, int net)
{
    act=0;
    if(trr<rr)
    {
        double sumOFf=0;
        random_device rd;
        mt19937 gen(rd());  // to seed mersenne twister.
        uniform_real_distribution<> dist1(0, 1);
        //uniform_int_distribution<> distall(0,llnn-1);
        ostringstream fn;
        vector<int> dy(llnn);
        vector<int> infects;
        vector<int> recovers;
        //vector<int> resus;
        for(int i=0;i<llnn;++i)
        {
            dy[i]=1;
            infects.push_back(i);

        }
        //=============================================================
        fn<<"./data/"<<net<<"/"<<probAlpha<<"/infection_fraction_"<<rr<<"_"<<probAlpha<<"_netE"<<net<<"_landaE"<<ln<<".dat";
        ofstream out1(fn.str().c_str(),ios_base::binary);
        //
        int isum=0;
        if(probAlpha<=1)
        {
            for(int ii=0;ii<rr;++ii)
            {
                act=ii;
                isum=0;
                uniform_int_distribution<> infs(0,(infects.size()-1));
                double infector=0,iactt=0,actt=0;
                iactt=infs(gen);
                infector=infects[iactt];
                actt=dist1(gen);
                if(actt<=(probAlpha/(probAlpha+probBeta)))
                {
                    for(int m=0;m<aa[infector].size();++m)
                    {
                        int mm=aa[infector][m];
                        if(dy[mm]==0)
                        {
                            double c2p=0;
                            c2p=dist1(gen);
                            if(c2p<=probAlpha)
                            {
                                dy[mm]=1;
                                infects.push_back(mm);
                            }
                        }
                    }
                }
                infects.erase(infects.begin()+iactt);
                dy[infector]=2;
                recovers.push_back(infector);
                //cout<<"a"<<endl;
                //---------------------------------------------
                //
                int k=0;
                while(k<recovers.size())
                {
                    int kk=recovers[k];
                    if(dy[kk]>=rt)
                    {
                        dy[kk]=0;
                        recovers.erase(recovers.begin()+k);
                        k--;
                    }else
                        dy[kk]++;
                    k++;
                }
                //---------------------------------------------
                isum=infects.size();
                double f1=0,f2=0,f3=0;
                f1=isum;
                f2=llnn;
                f3=f1/f2;
                if(f3==0)
                {
                    avgI=0;
                    sumOFf=0;
                    ii=rr;
                    break;
                }
                int um=rr-trr;
                if(ii >= um)
                    sumOFf=sumOFf+f3;
                //
                out1<< ii+1 <<' '<<f3<<endl;
            }
            avgI=sumOFf/trr;
        }else
        {
            for(int ii=0;ii<rr;++ii)
            {
                act=ii;
                isum=0;
                uniform_int_distribution<> infs(0,(infects.size()-1));
                double infector=0,iactt=0,actt=0;
                iactt=infs(gen);
                infector=infects[iactt];
                actt=dist1(gen);
                if(actt<=(probAlpha/(probAlpha+probBeta)))
                {
                    for(int m=0;m<aa[infector].size();++m)
                    {
                        int mm=aa[infector][m];
                        if(dy[mm]==0)
                        {
                            dy[mm]=1;
                            infects.push_back(mm);
                        }
                    }
                }
                actt=dist1(gen);
                if(actt<=(teta/(probBeta+teta)))
                {
                    infects.erase(infects.begin()+iactt);
                    dy[infector]=2;
                    recovers.push_back(iactt);
                }else
                {
                    infects.erase(infects.begin()+iactt);
                    dy[infector]=0;
                }
                //---------------------------------------------
                int k=0;
                while(k<recovers.size())
                {
                    int kk=recovers[k];
                    if(dy[kk]==5)
                    {
                        dy[kk]=0;
                        recovers.erase(recovers.begin()+k);
                        k--;
                    }else
                        dy[kk]++;
                    k++;
                }
                /*
                if(recovers.size()!=0)
                {
                    uniform_int_distribution<> recs(0,(recovers.size()-1));
                    double recovector=0,ractt=0,ract=0;
                    ractt=recs(gen);
                    recovector=recovers[ractt];
                    ract=dist1(gen);
                    if(ract<=(teta/(teta+probBeta)))
                    {
                        for(int m=0;m<aa[recovector].size();++m)
                        {
                            int mm=aa[recovector][m];
                            if(dy[mm]==1)
                            {
                                double c2p=0;
                                c2p=dist1(gen);
                                if(c2p<=teta)
                                {
                                    dy[mm]=2;
                                    recovers.push_back(mm);
                                    infects.erase(remove(infects.begin(), infects.end(), mm), infects.end());
                                }
                            }
                        }
                    }
                    recovers.erase(recovers.begin()+ractt);
                    dy[recovector]=0;
                }
                */
                //cout<<"b"<<endl;
                //          }
                //---------------------------------------------
                isum=infects.size();
                double f1=0,f2=0,f3=0;
                f1=isum;
                f2=llnn;
                f3=f1/f2;
                if(f3==0)
                {
                    avgI=0;
                    sumOFf=0;
                    ii=rr;
                    break;
                }
                int um=rr-trr;
                if(ii >= um)
                    sumOFf=sumOFf+f3;
                //
                out1<< ii+1 <<' '<<f3<<endl;
            }
            avgI=sumOFf/trr;
        }
    }else
    {
        cout<<"SIRS stoped!"<<endl<<"avrage loop is biger than loop try."<<endl;
    }
}
//
//
void NewDy::D3(int rr, double probAlpha, double probBeta, double teta, int llnn, iMatrix &aa, int trr, int ln,int net)
{
    act=0;
    if(trr<rr)
    {
        double sumOFf=0;
        random_device rd;
        mt19937 gen(rd());  // to seed mersenne twister.
        uniform_real_distribution<> dist1(0, 1);
        //uniform_int_distribution<> distall(0,llnn-1);
        ostringstream fn;
        vector<int> dy(llnn);
        vector<int> infects;
        vector<int> recovers;
        //vector<int> resus;
        for(int i=0;i<llnn;++i)
        {
            dy[i]=1;
            infects.push_back(i);
        }
        /*
        cout<<"inf:"<<infects.size()<<endl;
        cout<<"recs:"<<recovers.size()<<endl;
        cout<<llnn<<endl;
        */
        //=============================================================
        fn<<"./data/"<<net<<"/"<<probAlpha<<"/infection_fraction_"<<rr<<"_"<<probAlpha<<"_netE"<<net<<"_landaE"<<ln<<".dat";
        ofstream out1(fn.str().c_str(),ios_base::binary);
        //
        int isum=0;
        for(int ii=0;ii<rr;++ii)
        {
            act=ii;
            isum=0;
            uniform_int_distribution<> infs(0,(infects.size()-1));
            double infector=0,iactt=0,actt=0;
            iactt=infs(gen);
            infector=infects[iactt];
            actt=dist1(gen);
            if(actt<=(probAlpha/(probAlpha+probBeta)))
            {
                for(int m=0;m<aa[infector].size();++m)
                {
                    int mm=aa[infector][m];
                    if(dy[mm]==0)
                    {
                        double c2p=0;
                        c2p=dist1(gen);
                        if(c2p<=probAlpha)
                        {
                            dy[mm]=1;
                            infects.push_back(mm);
                        }
                    }
                }
            }
            infects.erase(infects.begin()+iactt);
            dy[infector]=2;
            recovers.push_back(infector);
            //cout<<"a"<<endl;
            //---------------------------------------------
            //
            int k=0;
            while(k<recovers.size())
            {
                int kk=recovers[k];
                if(dy[kk]==3)
                {
                    double r2s=dist1(gen);
                    if(r2s<=teta)
                    {
                        dy[kk]=0;
                        recovers.erase(recovers.begin()+k);
                        k--;
                    }else
                        dy[kk]++;
                }else
                    dy[kk]++;
                k++;
            }
            //
            //---------------------------------------------
            isum=infects.size();
            double f1=0,f2=0,f3=0;
            f1=isum;
            f2=llnn;
            f3=f1/f2;
            if(f3==0)
            {
                avgI=0;
                sumOFf=0;
                ii=rr;
                break;
            }
            int um=rr-trr;
            if(ii >= um)
                sumOFf=sumOFf+f3;
            //
            out1<< ii+1 <<' '<<f3<<endl;
        }
        avgI=sumOFf/trr;
    }else
    {
        cout<<"SIRS stoped!"<<endl<<"avrage loop is biger than loop try."<<endl;
    }
}
//
//
void NewDy::D4_PRE_sis(int rr, double probAlpha, double probBeta, int llnn, iMatrix &aa, int trr, int ln, int net)
{
    act=0;
    if(trr<rr)
    {
        double deltaT=0;
        double sumOFf=0;
        random_device rd;
        mt19937 gen(rd());  // to seed mersenne twister.
        uniform_real_distribution<> dist1(0, 1);
        //uniform_int_distribution<> distall(0,llnn-1);
        ostringstream fn;
        vector<int> dy(llnn);
        vector<int> infects;
        vector<int> recovers;
        //vector<int> resus;
        for(int i=0;i<llnn;++i)
        {
            dy[i]=1;
            infects.push_back(i);
        }
        //=============================================================
        fn<<"./data/"<<probAlpha<<"/infection_fraction_"<<rr<<"_"<<probAlpha<<"_ans"<<ln<<".dat";
        ofstream out1(fn.str().c_str(),ios_base::binary);
        //
        int isum=0;
        for(int ii=0;ii<rr;++ii)
        {
            act=ii;
            isum=0;
            uniform_int_distribution<> infs(0,(infects.size()-1));
            int infector=0,iactt=0;
            double actt=0;
            iactt=infs(gen);
            infector=infects[iactt];
            //
            double Ni=infects.size();
            double Nn=0;
            for(int nj=0;nj<Ni;++nj)
            {
                int c=infects[nj];
                Nn=Nn+aa[c].size();
            }
            //
            actt=dist1(gen);
            if(actt<=((probAlpha*Nn)/(Ni+probAlpha*Nn)))
            {
                int inftry=0,istate=0;
                uniform_int_distribution<> sess(0,(aa[infector].size()-1));
                while(inftry<aa[infector].size() && istate==0)
                {
                    int sinfec=0,sactt=0;
                    sactt=sess(gen);
                    sinfec=aa[infector][sactt];
                    if(dy[sinfec]==0)
                    {
                        dy[sinfec]=1;
                        infects.push_back(sinfec);
                        istate=1;
                    }
                    inftry++;
                }
            }else{
                infects.erase(infects.begin()+iactt);
                dy[infector]=0;
                //recovers.push_back(infector);
            }
            //cout<<"a"<<endl;
            //---------------------------------------------
            /*
            int k=0;
            while(k<recovers.size())
            {
                int kk=recovers[k];
                if(dy[kk]>=3)
                {
                    double r2s=dist1(gen);
                    if(r2s<=probBeta)
                    {
                    dy[kk]=0;
                    recovers.erase(recovers.begin()+k);
                    k--;
                    }else
                        dy[kk]++;
                }else
                    dy[kk]++;
                k++;
            }
            */
            deltaT=deltaT+(1/(Ni+(probAlpha*Nn)));
            //---------------------------------------------
            isum=infects.size();
            double f1=0,f2=0,f3=0;
            f1=isum;
            f2=llnn;
            f3=f1/f2;
            if(f3==0)
            {
                avgI=0;
                sumOFf=0;
                ii=rr;
                break;
            }
            int um=rr-trr;
            if(ii >= um)
                sumOFf=sumOFf+f3;
            //
            out1<< deltaT <<' '<<f3<<endl;
        }
        avgI=sumOFf/trr;
    }else
    {
        cout<<"SIS stoped!"<<endl<<"avrage loop is biger than loop try."<<endl;
    }
}
//
//
void NewDy::D4_PRE_sirs(int rr, double probAlpha, double probBeta, int llnn, iMatrix &aa, int trr, int ln,int net)
{
    act=0;
    if(trr<rr)
    {
        double deltaT=0;
        double sumOFf=0;
        random_device rd;
        mt19937 gen(rd());  // to seed mersenne twister.
        uniform_real_distribution<> dist1(0, 1);
        //uniform_int_distribution<> distall(0,llnn-1);
        ostringstream fn;
        vector<int> dy(llnn);
        vector<int> infects;
        vector<int> recovers;
        //vector<int> resus;
        for(int i=0;i<llnn;++i)
        {
            dy[i]=1;
            infects.push_back(i);

        }
        //=============================================================
        fn<<"./data/"<<net<<"/"<<probAlpha<<"/infection_fraction_"<<rr<<"_"<<probAlpha<<"_netE"<<net<<"_landaE"<<ln<<".dat";
        ofstream out1(fn.str().c_str(),ios_base::binary);
        //
        int isum=0;
        for(int ii=0;ii<rr;++ii)
        {
            act=ii;
            isum=0;
            uniform_int_distribution<> infs(0,(infects.size()-1));
            int infector=0,iactt=0;
            double actt=0;
            iactt=infs(gen);
            infector=infects[iactt];
            //
            double Ni=infects.size();
            double Nn=0;
            for(int nj=0;nj<Ni;++nj)
            {
                int c=infects[nj];
                Nn=Nn+aa[c].size();
            }
            //
            actt=dist1(gen);
            if(actt<=((probAlpha*Nn)/(Ni+probAlpha*Nn)))
            {
                int inftry=0,istate=0;
                uniform_int_distribution<> sess(0,(aa[infector].size()-1));
                while(inftry<aa[infector].size() && istate==0)
                {
                    int sinfec=0,sactt=0;
                    sactt=sess(gen);
                    sinfec=aa[infector][sactt];
                    if(dy[sinfec]==0)
                    {
                        dy[sinfec]=1;
                        infects.push_back(sinfec);
                        istate=1;
                    }
                    inftry++;
                }
            }else{
                infects.erase(infects.begin()+iactt);
                dy[infector]=2;
                recovers.push_back(infector);
            }
            //cout<<"a"<<endl;
            //---------------------------------------------
            //
            int k=0;
            while(k<recovers.size())
            {
                int kk=recovers[k];
                if(dy[kk]>=3)
                {
                    double r2s=dist1(gen);
                    if(r2s<=probBeta)
                    {
                        dy[kk]=0;
                        recovers.erase(recovers.begin()+k);
                        k--;
                    }else
                        dy[kk]++;
                }else
                    dy[kk]++;
                k++;
            }
            deltaT=deltaT+(1/(Ni+(probAlpha*Nn)));
            //---------------------------------------------
            isum=infects.size();
            double f1=0,f2=0,f3=0;
            f1=isum;
            f2=llnn;
            f3=f1/f2;
            if(f3==0)
            {
                avgI=0;
                sumOFf=0;
                ii=rr;
                break;
            }
            int um=rr-trr;
            if(ii >= um)
                sumOFf=sumOFf+f3;
            //
            out1<< deltaT <<' '<<f3<<endl;
        }
        avgI=sumOFf/trr;
    }else
    {
        cout<<"SIRS stoped!"<<endl<<"avrage loop is biger than loop try."<<endl;
    }
}
//
