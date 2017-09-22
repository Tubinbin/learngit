#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/RR.h>
#include <NTL/mat_RR.h>
#include <iostream>
#include <random>

using namespace std;
using namespace NTL;

int main()
{
long n, d, w, m, p; 
ZZ q;//input each parameter
//cout<<"input (long n)"<<"\n";
//cin>>n;
//cout<<"input (long d)"<<"\n";
//cin>>d;
//cout<<"input (long w)"<<"\n";
//cin>>w;
//cout<<"input (ZZ p)"<<"\n";
//cin>>p;
//cout<<"input (ZZ q)"<<"\n";
//cin>>q;
n=100, d=50, w=20, m=20, p=33,q=5971847;
ZZ_p::init(q);

mat_ZZ_p A, Z;//(pk) random A and (sk) random Z 
random(A, n, d);
random(Z, d, w);
cout<<A<<"\n";
cout<<Z<<"\n";

mat_ZZ_p A_Z;//def the product of A and Z
mul(A_Z, A, Z);
cout<< A_Z <<"\n";

//def the message matrix B
long l=5;//n=m*l
//1,def vector b
Mat<RR> b;
b.SetDims(l, 1);
for (int i=0;i<l;i++)
{b[i][0]=pow(2, i);}
cout<<b<<"\n";

mat_RR I;
mat_ZZ_p B;//def message matrix
ident(I,m);
//Mat<RR> B;
B.SetDims(n, m);

int ii=0;
for (int j=0; j<m; j++)
{
for (int i=0; i<m; i++)
{
for (int k=0;k<l;k++)
{
B[ii][j]=conv<ZZ_p>(conv<ZZ>(I[i][j]*b[k][0]*conv<RR>(q)/conv<RR>(p)+0.5));
ii=ii+1;
if (ii==n)
{
ii=0;
}
}
}
}
cout<<B<<"\n";

//def error matrix
double delta;
delta=1.0/(32*33*100);
unsigned seed=time(NULL);
default_random_engine generator (seed);
normal_distribution<double>distribution(0.0,delta);

Mat<ZZ_p> E;
E.SetDims(n, m);
for (int i=0;i<n;i++){
for (int j=0;j<m;j++){
E[i][j]=conv<ZZ_p>(ZZ(5971847*distribution(generator)+0.5));
}
}

cout <<"E="<<E << "\n";

//def ciphertex matrix C (function index injective)
Mat<ZZ_p> C;
C=A_Z+B+E;
cout<<"C="<<C<<"\n";
//*************************************************
//eva y=x*C
//first randomly choose 0/1 strings
ZZ binary;
binary=2;
ZZ_p::init(binary);
mat_ZZ_p x;
random(x, 1, n);
cout<<"x="<<x<<"\n";

ZZ_p::init(q);
mat_ZZ_p y1, y2;
y2=x*C;
cout<<"y2="<<y2<<"\n";
y1=x*A;
cout<<"y1="<<y1<<"\n";
//*************************************************
//def inversion algorithm
mat_RR t;
t.SetDims(1, m);
cout<<"t="<<t<<"\n";
//mat_ZZ_p Z_transpose, Z1, Z2;
mat_ZZ_p Z1, Z2;
//Z_transpose=transpose(Z);
//cout<<"Z_transpose="<<Z_transpose<<"\n";
//Z2=transpose(Z_transpose[0]);
//Z1=y1*Z2;
//cout<<"Z2="<<Z[0]<<"\n";
//Z1=transpose(Z[0]);
//cout<<"Z1="<<Z1<<"\n";
//for (int i=0;i<m;i++)
//{t[0][i]=conv<RR>(C[0][i]-conv<vec>y1*transpose(Z[i]))/conv<RR>(q)}

//t[0][0]=conv<RR>(C[0][0]-y1*Z1)/conv<RR>(q);
//cout<<"t="<<t<<"\n";
//Z1=y1*Z_transpose;
Z1=y1*Z;
cout<<"Z1="<<Z1<<"\n";
Z2=y2-Z1;
cout<<"Z2="<<Z2<<"\n";
for (int i=0;i<m;i++)
{t[0][i]=conv<RR>(conv<ZZ>(Z2[0][i]))/conv<RR>(q);}
cout<<"t="<<t<<"\n";

mat_RR closest, xx;
closest.SetDims(1, p);
xx.SetDims(1, m);


for (int j=0;j<m;j++)
{
for (int i=0;i<p;i++)
{
closest[0][i]=abs(abs(t[0][j]-double(i)/conv<RR>(p))-0.5);
}
//cout<<"closest="<<closest<<"\n";

RR max;
max=0;
int m;
for (int i=0; i<p;i++)
{
if (closest[0][i]>max)
{max=closest[0][i];
m=i;}
}
//cout<<"m="<<m<<"\n";
xx[0][j]=to_RR(m);
}
cout<<"xx="<<xx<<"\n";

//transform into binary
mat_ZZ xxx;
xxx.SetDims(1, n);
ZZ temp;
int iii=0;
for (int j=0; j<m; j++)
{
temp=conv<ZZ>(xx[0][j]);
for (int i=0;i<l;i++)
{
xxx[0][iii]=temp%2;
if (xxx[0][iii]!=0){temp=(temp-1)/2;}
else {temp=temp/2;}
iii++;
}}
cout<<"x="<<x<<"\n";
cout<<"xxx="<<xxx<<"\n";
//check solution
mat_ZZ solution;
solution.SetDims(1,n);
for (int i=0;i<l;i++)
{solution[0][i]=conv<ZZ>(x[0][i])-conv<ZZ>(xxx[0][i]);}
cout<<"solution="<<solution<<"\n";

}
