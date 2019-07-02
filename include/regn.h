#include "in.h"

#ifndef REGN_H
#define REGN_H

using namespace std;

class RegN
{
public:
    RegN();
    void Reg1(int N,int k,iMatrix & m1);//this is a (N node regular looped network) by k link per each node
    void Reg2(int N,int k,iMatrix & m1,int l);//add l random link to reg1 net
};

#endif // REGN_H
