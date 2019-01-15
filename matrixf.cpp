#include "matrixf.h"
//
using namespace::std;

//
void Matrixf::zMatrix(int nn,iMatrix& aa)
{
    for(int i=0;i<nn;i++)
    {
        for(int j=0;j<aa[i].size();j++)
        {
            aa[i][j]=0;
        }
    }
}
//---------------------------------------
void Matrixf::pMatrix(int nn,iMatrix& aa)
{
    for(int i=0;i<nn;i++)
    {
        for(int j=0;j<aa[i].size();j++)
        {
            cout<<aa[i][j]<<" ";
        }
        cout<<endl;
    }
}
//
void Matrixf::CnMatrix(int nn, iMatrix &aa)
{
    double sumk=0;
    for(int i=0;i<nn;++i)
    {
        for(int j=0;j<nn;++j)
            sumk=sumk+aa[i][j];
    }
    cout<<"average degree is: "<<double(sumk/nn)<<endl;
}
//
void Matrixf::dD(int nn, iMatrix &aa)
{
    ostringstream fnk;
    fnk<<"./data/VertexDegree.dat";
    ofstream kout(fnk.str().c_str(),ios_base::binary);
    for(int l=0;l<nn;++l)
    {
        int ssk=0;
        for(int ll=0;ll<nn;++ll)
            ssk=ssk+aa[l][ll];
        kout<<l<<" "<<ssk<<endl;
    }
}
//
void Matrixf::allMatrix(int nn,iMatrix& aa)
{
    for(int i=0;i<nn;i++)
    {
        for(int j=0;j<aa[i].size();j++)
        {
            aa[i][j]=1;
        }
        aa[i][i]=0;
    }
}
//
void Matrixf::vDeg(int nn, dRow &v, iMatrix &aa)
{
    for(int i=0;i<nn;++i)
    {
        v[i]=0;
        for(int j=0;j<nn;++j)
            v[i]=v[i]+aa[i][j];
    }
}
//
double Matrixf::avrageVec(int nn, dRow &v)
{
    double sum=0,avg=0;
    for(int i=0;i<nn;++i)
        sum=sum+v[i];
    double n=nn;
    avg=sum/n;
    return avg;
}
//
double Matrixf::varianceVec(int nn, dRow &v)
{
    double avg=avrageVec(nn,v);
    double sum=0,var=0;
    for(int i=0;i<nn;++i)
    {
        double value=v[i];
        sum=sum+pow((value-avg),2);
    }
    double n=nn;
    var=sum/n;
    return var;
}
//
void Matrixf::conMatrix(int nn, iMatrix &a, iMatrix &b)
{
    for(int i=0;i<nn;++i)
    {
        for(int j=0;j<nn;++j)
        {
            if(a[i][j]==1)
                b[i].push_back(j);
        }
    }
}
//
void Matrixf::DeleteMatElement(int nn, iMatrix &a)
{
    for(int i=0;i<nn;++i)
    {
        a[i].clear();
        a[i].shrink_to_fit();
    }
}
//
void Matrixf::mat_shrink_to_fit(int nn, iMatrix &a)
{
    for(int i=0;i<nn;++i)
    {
        a[i].shrink_to_fit();
    }
}
//
void Matrixf::dD_eff(int nn, iMatrix &aa, int ans)
{
    ostringstream fnk;
    fnk<<"./data/VertexDegree_ans"<<ans<<".dat";
    ofstream kout(fnk.str().c_str(),ios_base::binary);
    for(int l=0;l<nn;++l)
    {
        int ssk=0;
        ssk=aa[l].size();
        kout<<l<<" "<<ssk<<endl;
    }
}
