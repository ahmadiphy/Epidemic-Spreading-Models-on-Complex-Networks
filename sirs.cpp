#include "sirs.h"
#include "matrixf.h"
using namespace::std;
SIRS::SIRS()
{
    avgI=0;
}
//
//
void SIRS::NSLdynamic2(int rr, double probAlpha,int llnn, iMatrix &aa,int trr,int ln)
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
        vector<int> dy2(llnn);
        vector<double> vdeg(llnn);
        vector<double> xi(llnn);
        Matrixf mat1;
        mat1.vDeg(llnn,vdeg,aa);
#pragma omp parallel
        {
#pragma omp for
            for(int i=0;i<llnn;++i)
            {
                dy1[i]=0;
                dy2[i]=0;
            }
        }
        int c1=0,cr=0;
        do{
            c1=distall(gen);
            if(dy1[c1]==0)
            {
                dy1[c1]=1;
                cr++;
            }
        }while (cr<1);

        //
        fns<<"./data/infection_fraction_"<<rr<<"_"<<probAlpha<<"_ans"<<ln<<".dat";
        ofstream out1(fns.str().c_str(),ios_base::binary);
        int isum=0,isumr=0,isums=0;
        for(int ii=0;ii<rr;++ii)
        {
            act=ii;
            //

#pragma omp parallel
            {
#pragma omp for
                for(int ix=0;ix<llnn;++ix)
                {
                    xi[ix]=0;
                    for(int jx=0;jx<llnn;++jx)
                    {
                        if(aa[ix][jx]==1 && dy1[jx]==1)
                            xi[ix]=xi[ix]+1;
                    }
                    xi[ix]=(xi[ix]/vdeg[ix]);
                }
            }
            double var=0;
            var=mat1.varianceVec(llnn,vdeg);
            isum=0,isumr=0,isums=0;
            vector<double> probel(llnn);
#pragma omp parallel
            {
#pragma omp for
                for(int m=0;m<llnn;++m)
                {
                    if(dy1[m]==0)
                    {

                        probel[m]=0;
                        probel[m]=dist1(gen);
                        double lprob=0,pval=0; //probebality of infection
                        pval=-1*(pow((xi[m]-probAlpha),2))/(2*var);
                        lprob=exp(pval);
                        if(probel[m]<lprob)
                        {
                            dy2[m]=1;
                        }
                    }
                }

            }
            //---------------------------------------------
            for(int kk=0;kk<llnn;++kk)
            {
                if(dy1[kk]>0)
                    dy1[kk]++;
                if(dy1[kk]==0)
                    dy1[kk]=dy1[kk]+dy2[kk];
                dy2[kk]=0;
                if(dy1[kk]>=4)
                    dy1[kk]=0;
                if(dy1[kk]==1)
                {
                    isum=isum+1;
                }
            }
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
void SIRS::NSLdynamic1(int rr, double probAlpha,int llnn, iMatrix &aa,int trr,int ln)
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
        vector<int> dy2(llnn);
        vector<double> vdeg(llnn);
        Matrixf mat1;
        mat1.vDeg(llnn,vdeg,aa);
#pragma omp parallel
        {
#pragma omp for
            for(int i=0;i<llnn;++i)
            {
                dy1[i]=0;
                dy2[i]=0;
            }
        }
        int c1=0,cr=0;
        do{
            c1=distall(gen);
            if(dy1[c1]==0)
            {
                dy1[c1]=1;
                cr++;
            }
        }while (cr<1);

        //
        fns<<"./data/infection_fraction_"<<rr<<"_"<<probAlpha<<"_ans"<<ln<<".dat";
        ofstream out1(fns.str().c_str(),ios_base::binary);
        int isum=0,isumr=0,isums=0;
        for(int ii=0;ii<rr;++ii)
        {
            act=ii;
            isum=0,isumr=0,isums=0;
            vector<double> probel(llnn); //---------------------------------------------
#pragma omp parallel
            {
#pragma omp for
                for(int m=0;m<llnn;++m)
                {
                    if(dy1[m]==0)
                    {
                        int deg_a=0;
                        for(int j=0;j<llnn;j++)
                        {
                            if(dy1[j]==1 && aa[m][j]==1)
                            {
                                deg_a++;
                            }
                        }
                        double d3=0,d1=deg_a,d2=vdeg[m];
                        d3=d1/d2;
                        probel[m]=0;
                        probel[m]=dist1(gen);
                        if(probel[m]<(d3/2+probAlpha/2))
                        {
                            dy2[m]=1;
                        }
                    }
                }

            }
            //---------------------------------------------
            for(int kk=0;kk<llnn;++kk)
            {
                if(dy1[kk]>0)
                    dy1[kk]++;
                if(dy1[kk]==0)
                    dy1[kk]=dy1[kk]+dy2[kk];
                dy2[kk]=0;
                if(dy1[kk]>=4)
                    dy1[kk]=0;
                if(dy1[kk]==1)
                {
                    isum=isum+1;
                }
            }
            //
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
void SIRS::Sdynamic(int rr, double probAlpha,int llnn, iMatrix &aa)
{
    double usep=1-probAlpha;
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_real_distribution<> dist1(0, 1);
    uniform_int_distribution<> distall(0,llnn-1);
    ostringstream fn;
    iMatrix dy(llnn, iRow(2));
    vector<int> infects;
    vector<int> recovers;
    vector<int> resus;
    for(int i=0;i<llnn;++i)
    {
        dy[i][0]=0;
        dy[i][1]=0;
    }
    int ksum=0;
    int c1=0;
    do{
        c1=0;
        c1=distall(gen);
        if(dy[c1][0]==0)
        {
            dy[c1][0]=1;
            infects.push_back(c1);
            ksum++;
        }
    }while(ksum<10);
    //
    fn<<"./infection_fraction"<<rr<<"_"<<probAlpha<<".dat";
    ofstream out1(fn.str().c_str(),ios_base::binary);
    //
    int isum=0;
    for(int ii=0;ii<rr;++ii)
    {
        cout<<ii<<endl;
        isum=0;

        //---------------------------------------------
        for(int m=0;m<infects.size();++m)
        {
            int mm=infects[m];
            for(int j=0;j<llnn;j++)
            {
                if(aa[mm][j]==1 && dy[j][0]==0)
                {
                    double c2p=0;
                    c2p=dist1(gen);
                    if(c2p>usep)
                    {
                        dy[j][0]=1;
                        infects.push_back(j);
                    }
                }
            }
            dy[mm][1]++;
        }
        //---------------------------------------------
        for(int l=0;l<recovers.size();++l)
        {
            int ll=recovers[l];
            dy[ll][1]++;

        }
        //---------------------------------------------
        for(int m=0;m<infects.size();++m)
        {
            int mm=infects[m];
            if(dy[mm][1]==4)
            {
                double itr=0;
                itr=dist1(gen);
                if(itr>0)
                {
                    dy[mm][0]=2;
                    dy[mm][1]=0;
                    infects.erase (infects.begin()+m);
                    recovers.push_back(mm);
                }
            }
        }
        //---------------------------------------------
        for(int l=0;l<recovers.size();++l)
        {
            int ll=recovers[l];
            if(dy[ll][1]==9)
            {
                double rts=0;
                rts=dist1(gen);
                if(rts>0)
                {
                    dy[ll][0]=0;
                    dy[ll][1]=0;
                    recovers.erase (recovers.begin()+l);
                }
            }
        }
        //---------------------------------------------
        isum=infects.size();
        double f1=0,f2=0,f3=0;
        f1=isum;
        f2=llnn;
        f3=f1/f2;
        out1<< ii+1 <<' '<<f3<<endl;
    }
}
//
void SIRS::ENSdynamic(int rr, double probAlpha,int llnn, iMatrix &aa)
{
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_real_distribution<> dist1(0, 1);
    uniform_int_distribution<> distall(0,llnn-1);
    ostringstream fn;
    vector<int> dy1(llnn);
    vector<int> dy2(llnn);
#pragma omp parallel
    {
#pragma omp for
        for(int i=0;i<llnn;++i)
        {
            dy1[i]=0;
            dy2[i]=0;
        }
    }
    int c1=0,cr=0;
    do{
        c1=distall(gen);
        if(dy1[c1]==0)
        {
            dy1[c1]=1;
            cr++;
        }
    }while (cr<1);

    //
    fn<<"./infection_fraction"<<rr<<"_"<<probAlpha<<".dat";
    ofstream out1(fn.str().c_str(),ios_base::binary);
    //
    int isum=0;
    for(int ii=0;ii<rr;++ii)
    {
        isum=0;
        //---------------------------------------------
#pragma omp parallel
        {
#pragma omp for
            for(int m=0;m<llnn;++m)
            {
                if(dy1[m]==1 || dy1[m]==2)
                {
                    for(int j=0;j<aa[m].size();++j)
                    {
                        int am=aa[m][j];
                        if(dy1[am]==0)
                        {
                            double c2p=0;
                            c2p=dist1(gen);
                            if(c2p<probAlpha)
                            {
                                dy2[am]=1;
                            }
                        }
                    }
                }
            }
        }
        //---------------------------------------------
        for(int kk=0;kk<llnn;++kk)
        {
            if(dy1[kk]>0)
                dy1[kk]++;
            dy1[kk]=dy1[kk]+dy2[kk];
            dy2[kk]=0;
            if(dy1[kk]==4)
                dy1[kk]=0;
            if(dy1[kk]==1)
                isum=isum+1;
        }
        double f1=0,f2=0,f3=0;
        f1=isum;
        f2=llnn;
        f3=f1/f2;
        out1<< ii+1 <<' '<<f3<<endl;
    }
}
//
//
void SIRS::NSdynamic(int rr, double probAlpha,int llnn, iMatrix &aa,int trr,int ln)
{
    act=0;
    if(trr<rr)
    {
        double sumOFf=0;
        random_device rd;
        mt19937 gen(rd());  // to seed mersenne twister.
        uniform_real_distribution<> dist1(0, 1);
        uniform_int_distribution<> distall(0,llnn-1);
        ostringstream fns;
        vector<int> dy1(llnn);
        vector<int> dy2(llnn);
#pragma omp parallel
        {
#pragma omp for
            for(int i=0;i<llnn;++i)
            {
                dy1[i]=0;
                dy2[i]=0;
            }
        }
        int c1=0,cr=0;
        do{
            c1=distall(gen);
            if(dy1[c1]==0)
            {
                dy1[c1]=1;
                cr++;
            }
        }while (cr<1);

        //
        fns<<"./data/infection_fraction_"<<rr<<"_"<<probAlpha<<"_ans"<<ln<<".dat";
        ofstream out1(fns.str().c_str(),ios_base::binary);
        //
        int isum=0,isumr=0,isums=0;
        for(int ii=0;ii<rr;++ii)
        {
            act=ii;
            isum=0,isumr=0,isums=0;
            vector<double> probel(llnn); //---------------------------------------------
#pragma omp parallel
            {
#pragma omp for
                for(int m=0;m<llnn;++m)
                {
                    if(dy1[m]==1)
                    {
                        for(int j=0;j<llnn;j++)
                        {
                            if(dy1[j]==0 && aa[m][j]==1 && dy2[j]==0)
                            {
                                probel[m]=0;
                                probel[m]=dist1(gen);
                                if(probel[m]<probAlpha)
                                {
                                    dy2[j]=1;
                                }
                            }
                        }
                    }
                }

            }
            //---------------------------------------------
            for(int kk=0;kk<llnn;++kk)
            {
                if(dy1[kk]>0)
                    dy1[kk]++;
                if(dy1[kk]==0)
                    dy1[kk]=dy1[kk]+dy2[kk];
                dy2[kk]=0;
                if(dy1[kk]>=4)
                    dy1[kk]=0;
                if(dy1[kk]==1)
                {
                    isum=isum+1;
                }
            }
            double f1=0,f2=0,f3=0;
            f1=isum;
            f2=llnn;
            f3=f1/f2;
            //f44=f4/f2;
            //f55=f5/f2;
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
            //infectionF=f3;
            //cout<<f3<<endl;
            out1<< ii+1 <<' '<<f3<<endl;//' '<<f44<<' '<<f55<<endl;
            //out2<< ii+1 <<' '<<f44<<endl;
            //out3<< ii+1 <<' '<<f55<<endl;
        }
        //for(int dyi=0;dyi<llnn;++dyi)
        //out2<<dy4[dyi]<<' '<<dy3[dyi]<<endl;
        avgI=sumOFf/trr;
    }else
    {
        cout<<"SIRS stoped!"<<endl<<"avrage loop is biger than loop try."<<endl;
    }
}
//
//
void SIRS::NSdynamic2(int rr, double probAlpha,double probBeta,int llnn, iMatrix &aa,int trr,int ln)
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
        vector<int> dy2(llnn);
        vector<int> dy3(llnn);
        //vector<int> dy4(llnn);
#pragma omp parallel
        {
#pragma omp for
            for(int i=0;i<llnn;++i)
            {
                dy1[i]=1;
                dy2[i]=1;
                //dy3[i]=1;
                /*dy4[i]=0;
            for(int dy4i=0;dy4i<llnn;++dy4i)
                dy4[i]=dy4[i]+aa[i][dy4i];
                */
            }
        }
        /*
        int c1=0,cr=0;
        do{
            c1=distall(gen);
            if(dy1[c1]==0)
            {
                dy1[c1]=1;
                //dy3[c1]++;
                cr++;
            }
        }while (cr<1);
        */
        //
        fns<<"./data/"<<probAlpha<<"/infection_fraction_"<<rr<<"_"<<probAlpha<<"_ans"<<ln<<".dat";
        //fnss<<"./data/r_fraction_"<<rr<<"_"<<probAlpha<<".dat";
        //fnsss<<"./data/s_fraction_"<<rr<<"_"<<probAlpha<<".dat";
        ofstream out1(fns.str().c_str(),ios_base::binary);
        //ofstream out2(fnss.str().c_str(),ios_base::binary);
        //ofstream out3(fnsss.str().c_str(),ios_base::binary);
        //
        int isum=0,isumr=0,isums=0;
        for(int ii=0;ii<rr;++ii)
        {
            act=ii;
            // cout<<ii<<endl;
            isum=0,isumr=0,isums=0;
            //
            int state=0;
            int m;
            while(state==0){
                m=distall(gen);
                if(dy1[m]==1)
                {
                    double probel=0;
                    state=1;
                    probel=dist1(gen);
                    if(probel<=(probAlpha/(probAlpha+probBeta)))
                    {
                        for(int j=0;j<llnn;j++)
                        {
                            if(dy1[j]==0 && aa[m][j]==1 && dy2[j]==0)
                            {
                                probel=0;
                                probel=dist1(gen);
                                if(probel<=probAlpha)
                                {
                                    dy2[j]=1;
                                }
                            }
                        }
                    }
                    dy1[m]=2;
                    dy2[m]=0;
                }

            }
            //---------------------------------------------
            for(int kk=0;kk<llnn;++kk)
            {
                if(dy1[kk]>1)
                    dy1[kk]++;
                if(dy1[kk]==0)
                    dy1[kk]=dy1[kk]+dy2[kk];
                dy2[kk]=0;
                if(dy1[kk]>=3)
                    dy1[kk]=0;
                if(dy1[kk]==1)
                {
                    isum=isum+1;
                }
            }
            //
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
            //infectionF=f3;
            //cout<<f3<<endl;
            out1<< ii+1 <<' '<<f3<<endl;//' '<<f44<<' '<<f55<<endl;
            //out2<< ii+1 <<' '<<f44<<endl;
            //out3<< ii+1 <<' '<<f55<<endl;
        }
        //for(int dyi=0;dyi<llnn;++dyi)
        //out2<<dy4[dyi]<<' '<<dy3[dyi]<<endl;
        avgI=sumOFf/trr;
    }else
    {
        cout<<"SIRS stoped!"<<endl<<"avrage loop is biger than loop try."<<endl;
    }
}
//
//
void SIRS::NSdynamic3(int rr, double probAlpha,double probBeta,double teta,int llnn, iMatrix &aa,int trr,int ln)
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
        vector<int> dy2(llnn);
        vector<int> dy3(llnn);
        //vector<int> dy4(llnn);
#pragma omp parallel
        {
#pragma omp for
            for(int i=0;i<llnn;++i)
            {
                dy1[i]=1;
                dy2[i]=0;
                //dy3[i]=1;
                /*dy4[i]=0;
            for(int dy4i=0;dy4i<llnn;++dy4i)
                dy4[i]=dy4[i]+aa[i][dy4i];
                */
            }
        }
        /*
        int c1=0,cr=0;
        do{
            c1=distall(gen);
            if(dy1[c1]==0)
            {
                dy1[c1]=1;
                //dy3[c1]++;
                cr++;
            }
        }while (cr<1);
        */
        //
        fns<<"./data/"<<probAlpha<<"/infection_fraction_"<<rr<<"_"<<probAlpha<<"_ans"<<ln<<".dat";
        //fnss<<"./data/r_fraction_"<<rr<<"_"<<probAlpha<<".dat";
        //fnsss<<"./data/s_fraction_"<<rr<<"_"<<probAlpha<<".dat";
        ofstream out1(fns.str().c_str(),ios_base::binary);
        //ofstream out2(fnss.str().c_str(),ios_base::binary);
        //ofstream out3(fnsss.str().c_str(),ios_base::binary);
        //
        int isum=0,isumr=0,isums=0;
        for(int ii=0;ii<rr;++ii)
        {
            act=ii;
            // cout<<ii<<endl;
            isum=0,isumr=0,isums=0;
            //
            int state=0;
            int m;
            while(state==0){
                m=distall(gen);
                if(dy1[m]==1)
                {
                    double probel=0;
                    state=1;
                    probel=dist1(gen);
                    if(probel<(probAlpha/(probAlpha+probBeta)))
                    {
                        for(int j=0;j<llnn;j++)
                        {
                            if(dy1[j]==0 && aa[m][j]==1 && dy2[j]==0)
                            {
                                probel=0;
                                probel=dist1(gen);
                                if(probel<probAlpha)
                                {
                                    dy2[j]=1;
                                }
                            }
                        }
                    }
                    dy1[m]=2;
                    dy2[m]=0;
                }
            }
            //---------------------------------------------
            for(int kk=0;kk<llnn;++kk)
            {
                if(dy1[kk]>1)
                    dy1[kk]++;
                if(dy1[kk]==0)
                    dy1[kk]=dy1[kk]+dy2[kk];
                dy2[kk]=0;
                if(dy1[kk]>=3)
                {
                    double rRate=0;
                    rRate=dist1(gen);
                    if(rRate<=teta)
                        dy1[kk]=0;
                }
                if(dy1[kk]==1)
                {
                    isum=isum+1;
                }
            }
            //
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
            //infectionF=f3;
            //cout<<f3<<endl;
            out1<< ii+1 <<' '<<f3<<endl;//' '<<f44<<' '<<f55<<endl;
            //out2<< ii+1 <<' '<<f44<<endl;
            //out3<< ii+1 <<' '<<f55<<endl;
        }
        //for(int dyi=0;dyi<llnn;++dyi)
        //out2<<dy4[dyi]<<' '<<dy3[dyi]<<endl;
        avgI=sumOFf/trr;
    }else
    {
        cout<<"SIRS stoped!"<<endl<<"avrage loop is biger than loop try."<<endl;
    }
}
//
//
void SIRS::a_sis(int rr, double probAlpha,double probBeta,int llnn, iMatrix &aa,int trr,int ln)
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
            dy[infector]=0;
            //---------------------------------------------
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
void SIRS::a_sirs(int rr, double probAlpha,double probBeta,double teta,int llnn, iMatrix &aa,int trr,int ln)
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
            //---------------------------------------------
            for(int k=0;k<llnn;++k)
            {
                if(dy[k]==2)
                {
                    double r2s=0;
                    r2s=dist1(gen);
                    if(r2s<teta)
                        dy[k]=0;
                }
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
void SIRS::a_sirsE(int rr, double probAlpha,double probBeta,double teta,int llnn, iMatrix &aa,int trr,int ln)
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
        vector<int> recovers;
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
            //---------------------------------------------
            for(int k=0;k<recovers.size();++k)
            {
                int kk=recovers[k];
                if(dy[kk]==2)
                {
                    double r2s=0;
                    r2s=dist1(gen);
                    if(r2s<teta)
                    {
                        dy[kk]=0;
                        recovers.erase(recovers.begin()+k);
                    }
                }
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
void SIRS::new_sirs(int rr, double probAlpha,double probBeta,double teta,int llnn, iMatrix &aa,int trr,int ln)
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
        vector<int> recovers;
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
    }else
    {
        cout<<"SIRS stoped!"<<endl<<"avrage loop is biger than loop try."<<endl;
    }
}
//
//
void SIRS::new_sirs2(int rr, double probAlpha,double probBeta,double teta,int llnn, iMatrix &aa,int trr,int ln)
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
        vector<int> recovers;
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
            double infector=0,iactt=0,actt=0;
            iactt=infs(gen);
            infector=infects[iactt];
            actt=dist1(gen);
            if(actt<=(probAlpha/(probAlpha+probBeta+teta)))
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
            }else if(actt<=((probAlpha+teta)/(probAlpha+probBeta+teta)))
            {
                dy[infector]=2;
                recovers.push_back(infector);
            }else
            {
                dy[infector]=0;
            }
            infects.erase(infects.begin()+iactt);
            cout<<"cdc"<<endl;
            //---------------------------------------------
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
                    if(dy[mm]==0)
                    {
                        double c2p=0;
                        c2p=dist1(gen);
                        if(c2p<=teta)
                        {
                            dy[mm]=2;
                            recovers.push_back(mm);
                        }
                    }
                }
            }
            recovers.erase(recovers.begin()+ractt);
            dy[infector]=0;
       
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
void SIRS::new_sirs3(int rr, double probAlpha,double probBeta,double teta,int llnn, iMatrix &aa,int trr,int ln)
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
            double cha=0;
            cha=dist1(gen);
            if(cha<=0.5)
            {
                dy[i]=1;
                infects.push_back(i);
            }else if(cha<=0.75)
            {
                dy[i]=2;
                recovers.push_back(i);
            }else
            {
                dy[i]=0;
            }
        }
        /*
    int ksum=0;
    int c1=0;
    do{
        c1=distall(gen);
        if(dy[c1]==0)
        {
            dy[c1]=1;
            infects.push_back(c1);
            ksum++;
        }
    }while(ksum<llnn/2);
    */
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
            double infector=0,iactt=0,actt=0;
            iactt=infs(gen);
            infector=infects[iactt];
            actt=dist1(gen);
            if(actt<=(probAlpha/(probAlpha+teta)))
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
            }//else if(actt<=((probAlpha+teta)/(probAlpha+probBeta+teta)))
            //{
            infects.erase(infects.begin()+iactt);
            dy[infector]=0;
            //---------------------------------------------
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
    }else
    {
        cout<<"SIRS stoped!"<<endl<<"avrage loop is biger than loop try."<<endl;
    }
}

//
//
double SIRS::get_avgI()
{
    return avgI;
}
