#include "si.h"
using namespace::std;
//
SI::SI()
{
    //Check_result=0;
}
//
/*
void SI::Check_SI(int llnn, iMatrix &ac)
{
    Check_result=0;
    int rr=llnn;
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_int_distribution<> distall(0,llnn-1);
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
    int c1=0;
    c1=distall(gen);
    dy1[c1]=1;
    dy2[c1]=1;
    double f0=0,f00=0;
    f00=llnn;
    f0=1/f00;
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
                if(dy2[m]==1)
                {
                    for(int j=0;j<llnn;j++)
                    {
                        if(dy1[j]==0 && ac[m][j]==1)
                        {
                            dy1[j]=1;
                            dy2[j]=1;
                        }
                    }

                }
                dy2[m]=0;
            }

        }
        //---------------------------------------------
        for(int kk=0;kk<llnn;++kk)
        {
            if(dy1[kk]==1)
            {
                isum=isum+1;
            }

        }
        double f1=0,f2=0,f3=0;
        f1=isum;
        f2=llnn;
        f3=f1/f2;
        if(f3==1)
        {
            Check_result=1;
            ii=rr;
            break;
        }else if(f3==f0)
        {
            Check_result=0;
            ii=rr;
            break;
        }
        f0=f3;
    }
}
*/
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
        //isum=0;
        //---------------------------------------------
        int isize=infects.size();
        //cout<<"isize="<<isize<<endl;
        for(int m=iii;m<isize;++m)
        {
            //cout<<m;
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
        //---------------------------------------------
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
/*
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
    double f0=0,f00=0;
    f00=llnn;
    f0=1/f00;
    int isum=1;
    int iii=0,cc=0;
    for(int ii=0;ii<rr;++ii)
    {
        //isum=0;
        //---------------------------------------------
        int isize=infects.size();
        //cout<<"isize="<<isize<<endl;
        for(int m=iii;m<isize;++m)
        {
            //cout<<m;
            int im=infects[m];
            for(int j=0;j<ac[im].size();j++)
            {
                int jc=ac[im][j];
                if(dy1[jc]==0)
                {
                    dy1[jc]=1;
                    infects.push_back(jc);
                    isum++;
                }
            }
        }
        iii=isize-1;
        //---------------------------------------------
        //isum=infects.size();
        double f1=0,f2=0,f3=0;
        f1=isum;
        f2=llnn;
        f3=f1/f2;
        if(isum==llnn)
        {
            Check_result=1;
            ii=rr;
            break;
        }else if(f3==f0)
        {
            Check_result=0;
            ii=rr;
            break;

        }
        f0=f3;

    }
    cout<<llnn<<endl;
    cout<<"infs:"<<isum<<endl;
    if(llnn==isum)
        Check_result=1;
}
*/
