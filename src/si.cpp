#include "si.h"
using namespace::std;
//
SI::SI()
{
    Check_result=0;
}
void SI::ECheck_SI(int llnn, iMatrix &ac)
{
    Check_result=0;
    int rr=llnn;
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_int_distribution<> distall(0,llnn-1);
    vector<int> dy1(llnn);
    vector<int> infects;
#pragma omp parallel
    {
#pragma omp for
        for(int i=0;i<llnn;++i)
        {
            dy1[i]=0;
        }
    }
    int c1=0;
    c1=distall(gen);
    dy1[c1]=1;
    infects.push_back(c1);
    int iii=0;
    for(int ii=0;ii<rr;++ii)
    {
        int isize=infects.size();
        for(int m=iii;m<isize;++m)
        {
            int im=infects[m];
            for(int j=0;j<ac[im].size();j++)
            {
                int jc=ac[im][j];
                if(dy1[jc]==0)
                {
                    dy1[jc]=1;
                    infects.push_back(jc);
                }
            }
        }
        iii=isize-1;
    }
    cout<<llnn<<endl;
    cout<<"infs:"<<infects.size()<<"**"<<infects.capacity()<<endl;
    if(llnn==infects.size())
        Check_result=1;
    else
        Check_result=0;
    int sum=0;
    for(int i=0;i<llnn;++i)
    {
        if(dy1[i]==1)
            sum++;
    }
    cout<<sum<<endl;
}
//
int SI::get_Check_result()
{
    return Check_result;
}
