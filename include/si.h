#include "in.h"
#ifndef SI_H
#define SI_H


class SI
{
public:
    SI();
    void Check_SI(int llnn, iMatrix& ac);
    void ECheck_SI(int llnn, iMatrix& ac);
    int get_Check_result();
private:
    int Check_result;
};

#endif // SI_H
