#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/RR.h>
#include <NTL/mat_RR.h>
#include <iostream>


using namespace std;
using namespace NTL;

//input secret matrix
///def combining function

ZZ_p combine(mat_ZZ_p ID_com, mat_ZZ_p TD_com, int T_com)
{
ZZ_p f0;
f0=0;
mat_ZZ_p coff_com;
coff_com.SetDims(1, T_com); 
for (int i=0;i<T_com-1;i++)
{
coff_com[0][i]=1;
for (int j=0;j<T_com-1;j++)
{
if (i==j) {continue;}
coff_com[0][i]=coff_com[0][i]*(-ID_com[0][j])*inv(ID_com[0][i]-ID_com[0][j]);
}
}
for (int i=0; i<T_com-1;i++)
{f0=f0+coff_com[0][i]*TD_com[0][1];}

return f0;
}


int main()
{
ZZ p;
p=7;
ZZ_p::init(p);
//poly0=sk*transpose(coff)；
//输入各个用户的私钥
//cin>>"输入私钥">>sk_ID>>"\n";
//sk_ID_1=[ID_1 y_1]
//根据t个用户身份，计算系数矩阵
//先定义系数矩阵的行列数coff=[0 0 ... 0]
//读入各个用户的私钥信息矩阵sk
mat_ZZ_p ID;
cin>>ID; 
//ID.SetDims(1, T_com);
//ID=[6 1 3];
cout<<"ID="<<ID<<"\n";

mat_ZZ_p TD; 
cin>>TD;
//TD=[3 6 1];
cout<<"TD="<<TD<<"\n";

int T=3;

mat_ZZ_p ID_com;
mat_ZZ_p TD_com;
int T_com;

ID_com=ID;
TD_com=TD;
T_com=T;

ZZ_p td0;
td0=combine(ID_com, TD_com, T_com);
cout<<"td0="<<td0<<"\n";
}




