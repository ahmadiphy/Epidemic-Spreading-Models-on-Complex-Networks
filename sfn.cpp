#include "sfn.h"
//
using namespace std;
//
void SFn::SFnetwork(int m,double pm,int m0,int nn, iMatrix& aa)
{
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_real_distribution<> dist1(0,1);
    vector<int> vdeg(nn);
    //
    for(int i=0;i<nn;++i)
        vdeg[i]=0;
    //
    int mm=pm*m;
    double sumk=0;
    for(int i=0;i<m;++i)
    {
        for(int j=0;j<mm;++j)
        {
            double pp=0;
            pp=dist1(gen);
            if(pp<pm && i!=j)
            {
                aa[i][j]=1;
            }
        }
    }
    for(int i=0;i<m;++i)
    {
        for(int j=0;j<m;++j)
            vdeg[i]=vdeg[i]+aa[i][j];
        sumk=sumk+vdeg[i];
    }
    for(int i=m;i<nn;++i)
    {
        uniform_int_distribution<> disti(0,i-1);
        int ii=0;
        do{
            double valn=0,pch=0;
            int rch=-1;
            rch=disti(gen);
            valn=vdeg[rch]/sumk;
            pch=dist1(gen);
            if(pch < valn && aa[i][rch]!=1)
            {
                aa[i][rch]=1;
                aa[rch][i]=1;
                vdeg[i]=vdeg[i]+1;
                vdeg[rch]=vdeg[rch]+1;
                sumk=sumk+2;
                ii++;
            }
        }while (ii<m0);
    }
}

