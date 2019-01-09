#include "nsirs.h"
using namespace std;
NSIRS::NSIRS()
{
    act=0;
    avgI=0;
}
//
//
void NSIRS::ns0(int rr, double probAlpha, double probBeta, int llnn, iMatrix &aa, int trr, int ln)
{
    act=0;
    if(trr<rr)
    {
        double sumOFf=0;
        random_device rd;
        mt19937 gen(rd());  // to seed mersenne twister.
        uniform_real_distribution<> dist1(0, 1);
        ostringstream fn;
        vector<int> dy(llnn);
        vector<int> infects;
        vector<int> ninfects;
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
            //cout<<"hihihi"<<endl;
            act=ii;
            isum=0;
            int isize=infects.size();
            for(int ij=0;ij<isize;++ij)
            {
                uniform_int_distribution<> infs(0,(infects.size()-1));
                int infector=0,iactt=0;
                iactt=infs(gen);
                infector=infects[iactt];
                double pp=dist1(gen);
                if(pp<probAlpha/(probAlpha+probBeta))
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
                infects.erase(infects.begin()+iactt);
                dy[infector]=0;

            }
            infects.shrink_to_fit();
            //cout<<"hihihi"<<endl;
            //--------------------------------------
            //*-*-*-*-*-*-*-*-*-*-*-*-*-*
            int ninf=ninfects.size();
            for(int ni=0;ni<ninf;++ni)
            {
                infects.push_back(0);
                infects[ni]=ninfects[ni];
            }
            //*-*-*-*-*-*-*-*-*-*-*-*-*-*-
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
void NSIRS::ns1(int rr, double probAlpha, double probBeta, double teta, int llnn, iMatrix &aa, int trr, int ln)
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
        vector<int> ninfects;
        for(int i=0;i<llnn/2;++i)
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
            //cout<<"hihihi"<<endl;
            act=ii;
            isum=0;
            int isize=infects.size();
            for(int ij=0;ij<isize;++ij)
            {
                uniform_int_distribution<> infs(0,(infects.size()-1));
                double infector=0,iactt=0;
                iactt=infs(gen);
                infector=infects[iactt];
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
                //cout<<"loob"<<ii<<"--"<<iactt<<endl;
                infects.erase(infects.begin()+iactt);
                dy[infector]=2;
                //cout<<"loob"<<ii<<"--"<<iactt<<endl;
                recovers.push_back(infector);
            }
            //cout<<"hihihi"<<endl;
            //--------------------------------------
            int k=0;
            while(k<recovers.size())
            {
                int kk=recovers[k];
                if(dy[kk]>=3)
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
            //*-*-*-*-*-*-*-*-*-*-*-*-*-*
            int ninf=ninfects.size();
            for(int ni=0;ni<ninf;++ni)
            {
                infects.push_back(0);
                infects[ni]=ninfects[ni];
            }
            //*-*-*-*-*-*-*-*-*-*-*-*-*-*-
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
void NSIRS::ns2(int rr, double probAlpha,int llnn, iMatrix &aa,int trr,int ln)
{
    act=0;
    if(trr<rr)
    {
        double sumOFf=0;
        random_device rd;
        mt19937 gen(rd());  // to seed mersenne twister.
        uniform_real_distribution<> dist1(0, 1);
        uniform_int_distribution<> distall(0,llnn-1);
        ostringstream fns;//,fnss,fnsss;
        vector<int> dy1(llnn);
        vector<int> infects;
        vector<int> ninfects;
        vector<int> recovers;
#pragma omp parallel
        {
#pragma omp for
            for(int i=0;i<llnn;++i)
            {
                dy1[i]=0;
            }
        }
        int c1=0,cr=0;
        do{
            c1=distall(gen);
            if(dy1[c1]==0)
            {
                dy1[c1]=1;
                infects.push_back(c1);
                cr++;
            }
        }while (cr<1);

        //
        fns<<"./data"<<probAlpha<<"/infection_fraction_"<<rr<<"_"<<probAlpha<<"_ans"<<ln<<".dat";
        ofstream out1(fns.str().c_str(),ios_base::binary);
        //
        int isum=0;
        for(int ii=0;ii<rr;++ii)
        {
            act=ii;
            isum=0;
            while(infects.size()>0)
            {
                int m=infects[0];
                for(int jj=0;jj<aa[m].size();++jj)
                {
                    int j=aa[m][jj];
                    if(dy1[j]==0)
                    {
                        double probel=0;
                        probel=dist1(gen);
                        if(probel<probAlpha)
                        {
                            dy1[j]=1;
                            ninfects.push_back(j);
                        }
                    }

                }
                dy1[m]=2;
                infects.erase(infects.begin()+0);
                recovers.push_back(m);
            }
            //cout<<"hioioioi"<<endl;
            infects.clear();
            infects.shrink_to_fit();
            //---------------------------------------------
            int k=0;
            while(k<recovers.size())
            {
                int kk=recovers[k];
                if(dy1[kk]>=4)
                {
                    dy1[kk]=0;
                    recovers.erase(recovers.begin()+k);
                    k--;
                }else
                    dy1[kk]++;
                k++;
            }
            //*-*-*-*-*-*-*-*-*-*-*-*-*-*
            int ninf=ninfects.size();
            for(int ni=0;ni<ninf;++ni)
            {
                infects.push_back(0);
                infects[ni]=ninfects[ni];
            }
            //*-*-*-*-*-*-*-*-*-*-*-*-*-*-
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
        //
        avgI=sumOFf/trr;
        //cout<<avgI<<endl;
    }else
    {
        cout<<"SIRS stoped!"<<endl<<"avrage loop is biger than loop try."<<endl;
    }
}
