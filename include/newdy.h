#include "in.h"
#ifndef NEWDY_H
#define NEWDY_H

using namespace std;

class NewDy
{
public:
    NewDy();
    void D00(int rr, double probAlpha, double probBeta, int llnn, iMatrix &aa, int trr, int ln, int net);
    double get_avgI();
private:
    double avgI;
    int act;
};

#endif // NEWDY_H
