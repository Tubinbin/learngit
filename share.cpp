#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/RR.h>
#include <NTL/mat_RR.h>
#include <iostream>


using namespace std;
using namespace NTL;

///def sharing function
ZZ_p function(mat_ZZ_p coff_fun, int N_fun, ZZ_p x_fun)
{
    ZZ_p y;
y = 0;
    for(int j = N_fun - 1; j >= 0; j--)
    {
        y = y * x_fun + coff_fun[0][j];
    }
    return y;
}


int main()
{
ZZ p;
p=7;
ZZ_p::init(p);
mat_ZZ_p coff;
//cout<<"input coff matrix"<<coff<<"\n";
//cin>>coff;
random(coff, 1, 3);
cout<<"coff="<<coff<<"\n";

int N=3;

ZZ_p x1;
x1=1;
ZZ_p x2;
x2=2;
ZZ_p x3;
x3=3;

mat_ZZ_p coff_fun;
coff_fun=coff;
int N_fun=N;
ZZ_p x_fun;
 
x_fun=x1;
ZZ_p z1;
z1=function(coff_fun, N_fun, x_fun);
cout<<"z1="<<z1<<"\n";

x_fun=x2;
ZZ_p z2;
z2=function(coff_fun, N_fun, x_fun);
cout<<"z2="<<z2<<"\n";

x_fun=x3;
ZZ_p z3;
z3=function(coff_fun, N_fun, x_fun);
cout<<"z3="<<z3<<"\n";

}
