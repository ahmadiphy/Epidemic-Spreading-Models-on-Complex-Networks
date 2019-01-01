#include "in.h"
#include "matrixf.h"
#include "sirs.h"
#include "sfn.h"
#include "ern.h"
#include "wsn.h"
#include "regn.h"
#include "si.h"
#include "newdy.h"
#include "nsirs.h"
//
using namespace::std;
//
//
int main()
{
    omp_set_num_threads(3);
    //==== core dump error checking =================
    struct rlimit core_limit;
    core_limit.rlim_cur = RLIM_INFINITY;
    core_limit.rlim_max = RLIM_INFINITY;
    if (setrlimit(RLIMIT_CORE, &core_limit) < 0) {
        /* ERROR */
    }
    //================================================
    //|==========================| make directory |============================|
    const int dir_err = mkdir("data", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (-1 == dir_err)
    {
        cout<<"Error creating directory"<<endl;
        cout<<"-------------- No Problem ---------------"<<endl;
    }
    //|========================================================================
    //
    SIRS s1;
    Matrixf mf1;
    RegN rn1;
    SFn sf1;
    SI si1;
    NewDy nd1;
    //
    ostringstream infof;
    infof<<"./data/info.txt";
    ofstream outinf(infof.str().c_str(),ios_base::binary);
    //.........................................................................
    outinf<<"*-*-*-*-*-* Run Information *-*-*-*-*-*"<<endl;
    time_t now = time(0)+12600;
    char* dt = ctime(&now);
    tm *gmtm = gmtime(&now);
    dt = asctime(gmtm);
    outinf<<"The Run Start date and time is:"<< dt << endl;//iran date and time
    //==========================================================================
    //-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    int state=0;
    int l=13;
    int m0=2,ll=pow(2,l);
    int N=m0*ll;
    iMatrix a(N,iRow());
    //-------------------------------------------------------------------------
    //-------------------------------------------------------------------------
    //-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    cout<<"dynamics..."<<endl;
    //-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    //------------------ Ensembles --------------------------------------------
    //-------------------------------------------------------------------------
    int netans=10,Nans=50;
    //-------------------------------------------------------------------------
    double mpv,mpvstep,mpvlimit;
    //
    for(int net=0;net<netans;++net)
    {
        mpv=0.36;
        mpvstep=0.1;
        mpvlimit=0.37;
        ostringstream netdir;
        netdir<<"./data/"<<net;
        mkdir(netdir.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        //
        ostringstream fnm;
        fnm<<"./data/"<<net<<"/phase_diagram.dat";
        ofstream outM(fnm.str().c_str(),ios_base::binary);
        //oooooooooooooooooo Network ooooooooooooooooooooooooooo
	double avgK=0.5;
        rn1.ERnetwork(avgK,N,a);
        mf1.mat_shrink_to_fit(N,a);
        //oooooooooooooooooooooooooooooooooooooooooooooooooooooo
        while(mpv < mpvlimit)
        {
            ostringstream mdn;
            mdn<<"./data/"<<net<<"/"<<mpv;
            mkdir(mdn.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            double suminf=0,avginf=0;
            for(int i=0;i<Nans;++i)
            {
                cout<<"network : "<<net<<"  landa = "<<mpv<<" ans:"<<i<<endl;
                nd1.D00(10000000,mpv,1,N,a,1000,i,net);
                suminf=suminf+nd1.avgI;
            }
            avginf=suminf/Nans;
            outM<< mpv <<' '<<avginf<<endl;
            mpv=mpv+mpvstep;
        }
        outM.close();
        mf1.DeleteMatElement(N,a);
    }
    //-------
    outinf<<"Network: "<<m0<<"_"<<2<<"^"<<l<<"("<<ll<<") - N="<<m0*ll<<endl;
    outinf<<"beg: "<<mpv<<" - step: "<<mpvstep<<" - end: "<<mpvlimit<<endl;
    outinf<<"Number of Ensembles on Network : "<<Nans<<endl;
    outinf<<"Number of Ensembles on Landa   : "<<netans<<endl;
    //--------
    //-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    //+++++++++++++++++++++++++++++++++++++++++++
    now = time(0)+12600;
    dt = ctime(&now);
    gmtm = gmtime(&now);
    dt = asctime(gmtm);
    outinf<<endl;
    outinf<<"The Run End date and time is:"<< dt << endl;//iran date and time
    outinf.close();
}

