#include "in.h"
#ifndef NSIRS_H
#define NSIRS_H

using namespace std;
class NSIRS
{
public:
    NSIRS();
    void ns0(int rr, double probAlpha, double probBeta, int llnn, iMatrix &aa, int trr, int ln);
    void ns1(int rr, double probAlpha, double probBeta, double teta, int llnn, iMatrix &aa, int trr, int ln);
    void ns2(int rr, double probAlpha, int llnn, iMatrix &aa, int trr, int ln);
private:
    double avgI;
    int act;
};

#endif // NSIRS_H
