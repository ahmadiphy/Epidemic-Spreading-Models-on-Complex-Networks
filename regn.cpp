#include "regn.h"
using namespace std;

RegN::RegN()
{

}
//
//
void RegN::Reg0(int nm, int sm, iMatrix &m1)
{
    Hire testh;
    int N=nm*sm,k=sm-1;
    Reg1(N,k,m1);
}

//
//
void RegN::Reg1(int N, int k, iMatrix &m1)
{
    if(N>1 && k>0 && k<N)
    {
        int l=k/2;
        for(int i=0;i<N;++i)
        {
            for(int j=i+1;j<(i+l+1);++j)
            {
                if(j<N)
                {
                    m1[i][j]=1;
                    m1[j][i]=1;
                }else
                {
                    m1[i][j-N]=1;
                    m1[j-N][i]=1;
                }
            }
        }
    }else
    {
        cout<<"check your input !"<<endl;
    }
}
//
//
void RegN::Reg2(int N, int k, iMatrix &m1, int l)
{
    if(l<(N-k) && l>0)
    {
    Reg1(N,k,m1);
    for(int i=0;i<N;i+=2)
    {
       int l=i+(k/2)+1;
        if(l<N)
        {
            m1[i][l]=1;
            m1[l][i]=1;
        }else
        {
            m1[i][l-N]=1;
            m1[l-N][i]=1;
        }
    }
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    //uniform_real_distribution<> dist1(0, 1);
    uniform_int_distribution<> distall(0,N-1);
    int ll=0;
    do{
        int n1=0,n2=0;
        n1=distall(gen);
        n2=distall(gen);
        if(n1!=n2 && m1[n1][n2]==0)
        {
            m1[n1][n2]=1;
            m1[n2][n1]=1;
            ll++;
        }
    }while(ll<l);
    }else
    {
        cout<<"check input!"<<endl;
    }
}
