#include "in.h"
#ifndef NEWDY_H
#define NEWDY_H

using namespace std;

class NewDy
{
public:
    NewDy();
    void D00(int rr, double probAlpha, double probBeta, int llnn, iMatrix &aa, int trr, int ln, int net);
    void D0(int rr, double probAlpha, double probBeta, int llnn, iMatrix &aa, int trr, int ln,int net);
    void D01(int rr, double probAlpha, double probBeta, int llnn, iMatrix &aa, int trr, int ln, int net);
    void D1(int rr, double probAlpha, double probBeta, int llnn, iMatrix &aa, int trr, int ln, int net);
    void D2(int rr, double probAlpha, double probBeta, double teta,int rt, int llnn, iMatrix &aa, int trr, int ln,int net);
    void D3(int rr, double probAlpha, double probBeta, double teta, int llnn, iMatrix &aa, int trr, int ln, int net);
    void D4_PRE_sis(int rr, double probAlpha, double probBeta, int llnn, iMatrix &aa, int trr, int ln,int net);
    void D4_PRE_sirs(int rr, double probAlpha, double probBeta, int llnn, iMatrix &aa, int trr, int ln, int net);
    double avgI;
    int act;
};

#endif // NEWDY_H
