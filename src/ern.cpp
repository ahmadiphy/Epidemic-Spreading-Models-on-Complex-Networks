//A class to generate ER random networks
#include "ern.h"
//
using namespace::std;
//
void ERn::ERnetwork(double avgK, int nn, iMatrix &aa)
{
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_int_distribution<> distall(0,nn-1);
    double linkss;
    linkss = double((avgK*nn)/2);
    int links;
    links = int(linkss);
    int i = 0;
    while(i < links){
        int r1 = 0, r2 = 0;
        r1 = distall(gen);
        r2 = distall(gen);
        if(aa[r1][r2] == 0 && r1 != r2)
        {
            aa[r1][r2] = 1;
            aa[r2][r1] = 1;
            ++i;
        }
    }
}
