#include "wsn.h"
using namespace std;
void WSn::WSnetwork(double reW, int Z, int nn, iMatrix& aa)
{
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_int_distribution<> distall(0,nn-1);
    uniform_real_distribution<> dist1(0,1);
    uniform_int_distribution<> dist0(0,1);
    int re=(Z/2)+1;
    for(int i=0;i<nn;++i)
    {
        for(int j=1;j<re;++j)
        {
            int ii=i+j;
            if(ii>=nn)
                ii=(ii-nn);
            aa[i][ii]=1;
            aa[ii][i]=1;
        }
    }
    for(int i=0;i<nn;++i)
    {
        for(int j=0;j<nn;++j)
        {
            if(aa[i][j]==1)
            {
                double r1=0;
                r1=dist1(gen);
                if(r1<reW)
                {
                    aa[i][j]=0;
                    aa[j][i]=0;
                    int r2=-1;
                    r2=dist0(gen);
                    if(r2==0)
                    {
                        int ch=0;
                        do{
                            int nj=-1;
                            nj=distall(gen);
                            if(nj!=j && nj!=i)
                            {
                                aa[i][nj]=1;
                                aa[nj][i]=1;
                                ch=3;
                            }
                        }while(ch<2);
                    }else if(r2==1)
                    {
                        int ch=0;
                        do{
                            int ni=-1;
                            ni=distall(gen);
                            if(ni!=i && ni!=j)
                            {
                                aa[j][ni]=1;
                                aa[ni][j]=1;
                                ch=3;
                            }
                        }while(ch<2);
                    }
                }
            }
        }
    }

}

