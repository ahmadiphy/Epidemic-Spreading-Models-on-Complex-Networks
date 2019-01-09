#include "in.h"
#ifndef SIRS_H
#define SIRS_H


class SIRS
{
public:
    SIRS();
    void Sdynamic(int rr, double probAlpha, int llnn, iMatrix &aa);
    void ENSdynamic(int rr, double probAlpha, int llnn, iMatrix &aa);
    void NSdynamic(int rr, double probAlpha, int llnn, iMatrix &aa, int trr,int ln);
    void NSLdynamic1(int rr, double probAlpha, int llnn, iMatrix &aa, int trr,int ln);
    void NSLdynamic2(int rr, double probAlpha, int llnn, iMatrix &aa, int trr,int ln);
    void NSdynamic2(int rr, double probAlpha, double probBeta, int llnn, iMatrix &aa, int trr, int ln);
    void NSdynamic3(int rr, double probAlpha, double probBeta,double teta, int llnn, iMatrix &aa, int trr, int ln);
    void article_sis(int rr, double probAlpha, double probBeta, int llnn, iMatrix &aa, int trr, int ln);
    void article_sirs(int rr, double probAlpha, double probBeta, double teta, int llnn, iMatrix &aa, int trr, int ln);
    void article_sirsE(int rr, double probAlpha, double probBeta, double teta, int llnn, iMatrix &aa, int trr, int ln);
    void new_sirs(int rr, double probAlpha, double probBeta, double teta, int llnn, iMatrix &aa, int trr, int ln);
    void new_sirs2(int rr, double probAlpha, double probBeta, double teta, int llnn, iMatrix &aa, int trr, int ln);
    void new_sirs3(int rr, double probAlpha, double probBeta, double teta, int llnn, iMatrix &aa, int trr, int ln);
    double get_avgI();
 private:
    double avgI;
    int act;
};

#endif // SIRS_H
