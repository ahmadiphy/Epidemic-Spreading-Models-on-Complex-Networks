//A class to simulate new spreading dynamic
#include "newdy.h"

NewDy::NewDy()
{
    avgI = 0;
    act = 0;
}
//
void NewDy::D00(int rr, double probAlpha, double probBeta, int llnn, iMatrix &aa, int trr, int ln, int net)
{
    // a new version of SIRS model
    if(trr < rr)
    {
        double deactP = (probBeta/(probAlpha+probBeta));
        double sumOFf = 0;
        random_device rd;//random devide to generate random numbers
        mt19937 gen(rd());  // to seed mersenne twister.
        uniform_real_distribution<> dist1(0, 1);
        ostringstream fn;
        vector<int> dy(llnn);
        vector<int> infects;
        for(int i=0; i<llnn; ++i)
        {
            dy[i] = 1;
            infects.push_back(i);

        }
        //=============================================================
        fn << "./data/" << net << "/" << probAlpha << "/infection_fraction_" << rr << "_" << probAlpha << "_netE" << net << "_landaE" << ln << ".dat";
        ofstream out1(fn.str().c_str(), ios_base::binary);
        //
        int isum = 0;
        for(int ii=0; ii<rr; ++ii)
        {
            act = ii;
            isum = 0;
            random_device rdinf;
            mt19937 geninf(rdinf());
            uniform_int_distribution<> infs(0, (infects.size()-1));
            uniform_real_distribution<> dist1l(0, 1);
            double actt = 0;
            int infector = 0, iactt = 0;
            iactt = infs(geninf);
            infector = infects[iactt];
            actt = dist1(gen);
            if(actt <= deactP)
            {
                infects.erase(infects.begin() + iactt);
                dy[infector] = 0;
            }else
            {
                vector<int> suss;
                for(int m=0; m<aa[infector].size(); ++m)
                {
                    int mm = aa[infector][m];
                    if(dy[mm] == 0)
                    {
                        suss.push_back(mm);
                    }
                }
                if(suss.size() > 0)
                {
                    for(int j=0; j<suss.size(); ++j)
                    {
                      int s2i = suss[j];
                      double sp = dist1l(gen);
                      if(sp <= probAlpha)
                      {
                        dy[s2i] = 1;
                        infects.push_back(s2i);
                      }
                    }
                }
                infects.erase(infects.begin() + iactt);
                dy[infector] = 0;
                suss.clear();
                suss.shrink_to_fit();
            }
            //---------------------------------------------
            //---------------------------------------------
            isum=infects.size();
            double f1 = 0, f2 = 0, f3 = 0;
            f1 = isum;
            f2 = llnn;
            f3 = f1/f2;
            if(f3 == 0)
            {
                avgI = 0;
                sumOFf = 0;
                ii = rr;
                break;
            }
            int um = rr - trr;
            if(ii >= um)
                sumOFf = sumOFf + f3;
            //
            out1 << ii+1 << ' ' << f3 << endl;
        }
        avgI = sumOFf / trr;

    }else
    {
        cout << "SIRS stoped!" << endl << "avrage loop is biger than loop try." <<endl;
    }
}
//
double NewDy::get_avgI()
{
    return avgI;
}
