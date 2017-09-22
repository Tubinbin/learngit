#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/RR.h>
#include <NTL/mat_RR.h>
#include <iostream>
#include <random>

using namespace std;
using namespace NTL;

///def hash function
mat_ZZ_p hashfunction(long ww_hash, long t_hash, long m_hash, mat_ZZ_p a_hash, mat_ZZ_p b_hash, mat_ZZ_p b_abo_hash)
{
mat_ZZ_p yy;//branch matrix
long mm;
mm=m_hash*ww_hash;
yy.SetDims(mm, t_hash);
for(int j=0;j<mm;j++)
	{
		for (int i=0;i<t_hash;i++)
			{
				yy[j][i]=a_hash[0][j]*b_abo_hash[0][i]+b_hash[0][j];
			}
	}
cout<<"yy="<<yy<<"\n";

mat_ZZ_p y;
long w_hash;
w_hash=ww_hash*t_hash;
y.SetDims(m_hash, w_hash);
for (int k=0; k<m_hash; k++)
{
	for (int i=0; i<ww_hash; i++)
	{
		for (int j=0; j<t_hash; j++)
		{
			y[k][j+i*t_hash]=yy[k*ww_hash+i][j];//将yy的m*ww行 t列的矩阵，按行分成ww块，依次链接成m行w列，
		}
	}
}
cout<<"y="<<y<<"\n";
return y;
}

int main()
{

long n, d, w, m; 
ZZ p, q;//input each parameter
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
n=8, d=5, m=4, w=20, p=5 ,q=23981;
ZZ_p::init(q);

mat_ZZ_p A, Z;//(pk) random A and (sk) random Z 
random(A, n, d);
random(Z, d, w);
cout<<A<<"\n";
cout<<Z<<"\n";

mat_ZZ_p A_Z;//def the product of A and Z
mul(A_Z, A, Z);
cout<<"A_Z"<< A_Z <<"\n";


//def message matrix
long t=5;//设计t时，保证是能整除w；
long ww;
ww=w/t;//定义消息矩阵的长度（尽可能保证满秩）w=m+t+t'(而扩展出来的矩阵需要是t的整数倍)
ZZ_p::init(p);

mat_ZZ_p b_abo;
random(b_abo, 1, t);
cout<<"b_abo="<<b_abo<<"\n";
//cin>>b_abo //input lossy branch 

mat_ZZ_p aa,bb;//这里需要的aa和bb个数是m*ww个，因为消息矩阵是m行w列，而w需要保证是t的整数倍（即w=ww*t）
random(aa, 1, ww*m);
random(bb, 1, ww*m);
cout<<"aa="<<aa<<"\n";
cout<<"bb="<<bb<<"\n";

long t_hash, m_hash, ww_hash;
mat_ZZ_p a_hash, b_hash, b_abo_hash;
ww_hash=ww;
m_hash=m;
t_hash=t;
a_hash=aa;
b_hash=bb;
b_abo_hash=b_abo;

mat_ZZ_p hash_b;
hash_b.SetDims(m, w);
cout<<"hash_b="<<hash_b<<"\n";
hash_b=hashfunction(ww_hash, t_hash, m_hash, a_hash, b_hash, b_abo_hash);
cout<<"hash(mat_ZZ_p b_abo)="<<hash_b<<"\n";
//convert into negative 
hash_b=-hash_b;
cout<<"hash_b="<<hash_b<<"\n";



//def the message matrix B
ZZ_p::init(q);

long l=2;//n=m*l
//1,def vector b
Mat<RR> b;
b.SetDims(l, 1);
for (int i=0;i<l;i++)
	{
		b[i][0]=pow(2, i);
	}
cout<<"b="<<b<<"\n";

mat_ZZ_p B;//def message matrix
//Mat<RR> B;
B.SetDims(n, w);

int ii=0;
for (int j=0; j<w; j++)
	{
		for (int i=0; i<m; i++)
			{
				for (int k=0;k<l;k++)
					{
						B[ii][j]=conv<ZZ_p>(conv<ZZ>(conv<RR>(conv<ZZ>(hash_b[i][j]))*b[k][0]*conv<RR>(q)/conv<RR>(p)+0.5));
						ii=ii+1;
						if (ii==n)
							{
								ii=0;
							}
					}
			}
	}
cout<<"B="<<B<<"\n";

//def error matrix
double delta;
delta=1.0/(3200);
unsigned seed=time(NULL);
default_random_engine generator (seed);
normal_distribution<double>distribution(0.0,delta);

Mat<ZZ_p> E;
E.SetDims(n, w);
for (int i=0;i<n;i++)
	{
		for (int j=0;j<w;j++)
			{
				E[i][j]=conv<ZZ_p>(ZZ(23981*distribution(generator)+0.5));
			}
		}

cout <<"E="<<E << "\n";

//def ciphertex matrix C (function index injective)
Mat<ZZ_p> C;
C=A_Z+B+E;
cout<<"C="<<C<<"\n";


//*************************************************
//eva y=x*C
// input a branch b_branch
ZZ_p::init(p);
mat_ZZ_p b_branch;
//b1_branch=b_abo;
random(b_branch, 1, t);
cout<<"b_branch="<<b_branch<<"\n";
//cin>>b_abo //input a branch 

//use hash function to convert the branch into a matrix
mat_ZZ_p b_branch_hash;
b_branch_hash=b_branch;
a_hash=aa;
b_hash=bb;

mat_ZZ_p hash_b_branch;
hash_b_branch.SetDims(m, w);
hash_b_branch=hashfunction(ww_hash, t_hash, m_hash, a_hash, b_hash, b_branch_hash);
cout<<"hash_b_branch="<<hash_b_branch<<"\n";


//def the message matrix B
ZZ_p::init(q);

mat_ZZ_p B_branch;//def message matrix
//Mat<RR> B_branch;
B_branch.SetDims(n, w);

ii=0;
for (int j=0; j<w; j++)
	{
		for (int i=0; i<m; i++)
			{
				for (int k=0;k<l;k++)
					{
						B_branch[ii][j]=conv<ZZ_p>(conv<ZZ>(conv<RR>(conv<ZZ>(hash_b_branch[i][j]))*b[k][0]*conv<RR>(q)/conv<RR>(p)+0.5));
						ii=ii+1;
						if (ii==n)
							{
								ii=0;
							}
					}
			}
	}
cout<<"B_branch"<<B_branch<<"\n";

//renew function index (ciphertext) C
Mat<ZZ_p> C_branch;
C_branch=C+B_branch;
cout<<"C_branch="<<C_branch<<"\n";

//randomly choose 0/1 strings
ZZ binary;
binary=2;
ZZ_p::init(binary);
mat_ZZ_p x;
random(x, 1, n);
cout<<"x="<<x<<"\n";

ZZ_p::init(q);
mat_ZZ_p y1, y2;
y2=x*C_branch;
//cout<<"y2="<<y2<<"\n";
y1=x*A;
//cout<<"y1="<<y1<<"\n";
//*************************************************

//def inversion algorithm
mat_RR tt;
tt.SetDims(1, w);
//cout<<"tt="<<tt<<"\n";

mat_ZZ_p Z1, Z2;

Z1=y1*Z;
//cout<<"Z1="<<Z1<<"\n";
Z2=y2-Z1;
//cout<<"Z2="<<Z2<<"\n";
for (int i=0;i<w;i++)
	{
		tt[0][i]=conv<RR>(conv<ZZ>(Z2[0][i]))/conv<RR>(q);
	}
//cout<<"tt="<<tt<<"\n";
mat_RR closest, xx;
closest.SetDims(1, 5);//def p is ZZ type not int
xx.SetDims(1, w);

for (int j=0;j<w;j++)
	{
		for (int i=0;i<p;i++)
			{
				closest[0][i]=abs(abs(tt[0][j]-double(i)/conv<RR>(p))-0.5);
			}
//cout<<"closest="<<closest<<"\n";

RR max;
max=0;
int m;
for (int i=0; i<p;i++)
	{
		if (closest[0][i]>max)
			{
				max=closest[0][i];
				m=i;
			}
	}
//cout<<"m="<<m<<"\n";
xx[0][j]=to_RR(m);
}
cout<<"xx="<<xx<<"\n";


//def branch matrix
ZZ_p::init(p);
mat_ZZ_p H;
H=hash_b_branch+hash_b;

cout<<"H="<<H<<"\n";

//def matrix v (vH=xx)

//*********************check
//ZZ_p::init(q);
mat_ZZ_p B1_branch;//def message matrix
//Mat<RR> B;
B1_branch.SetDims(n, w);

ii=0;
for (int j=0; j<w; j++)
	{
		for (int i=0; i<m; i++)
			{
				for (int k=0;k<l;k++)
					{
						B1_branch[ii][j]=conv<ZZ_p>(conv<ZZ>(conv<RR>(conv<ZZ>(H[i][j]))*b[k][0]));
						ii=ii+1;
						if (ii==n)
							{
								ii=0;
							}
					}
			}
	}
cout<<"B1_branch="<<B1_branch<<"\n";

ZZ_p::init(p);
mat_ZZ_p xxx;
xxx=x*B1_branch;
cout<<"xxx="<<xxx<<"\n";

mat_ZZ_p H1;
H1.SetDims(m+1, w);
for (int i=0; i<m+1; i++)
	{
		for (int j=0; j<w; j++)
			{
				if (i<4)
					{
						H1[i][j]=H[i][j];
					}
					else
						{
							H1[i][j]=xxx[0][j];
						}
			}
	}
cout<<"H1="<<H1<<"\n";

mat_ZZ_p H2;
H2=transpose(H1);
cout<<"H2="<<H2<<"\n";

long r;
r=gauss(H);
cout<<"r="<<r<<"\n";
cout<<"H="<<H<<"\n";

long r1;
r1=gauss(H1);
cout<<"r1="<<r1<<"\n";
cout<<"H1="<<H1<<"\n";

long r2;
r2=gauss(H2);
cout<<"H2="<<H2<<"\n";
cout<<"r2="<<r2<<"\n";

ZZ aaa;
aaa=-3;
cout<<aaa<<"\n";
ZZ_p aaaa;
aaaa=conv<ZZ_p>(aaa);
cout<<aaaa<<"\n";
ZZ aaaaa;
aaaaa=conv<ZZ>(aaaa);
cout<<aaaaa<<"\n";

mat_ZZ_p v;
v.SetDims(1, m);
for (int i=m-1; i>=0; i--)
	{	v[0][i]=H2[i][m];
		for (int j=m-1; j>i; j--)
			{			
v[0][i]=conv<ZZ_p>(conv<ZZ>(v[0][i])-conv<ZZ>(v[0][j]*H2[i][j]));
//v[0][i]=conv<ZZ_p>(conv<ZZ>(H2[i][m])-conv<ZZ>(v[0][j]*H2[i][j]));
}
		v[0][i]=v[0][i]*inv(H2[i][i]);
	}

cout<<"v="<<v<<"\n";

//transform into binary
mat_ZZ xxxx;
xxxx.SetDims(1, n);
ZZ temp;
int iii=0;
for (int j=0; j<m; j++)
{
temp=conv<ZZ>(v[0][j]);
for (int i=0;i<l;i++)
{
xxxx[0][iii]=temp%2;
if (xxxx[0][iii]!=0){temp=(temp-1)/2;}
else {temp=temp/2;}
iii++;
}}
cout<<"xxxx="<<xxxx<<"\n";

//check solution
mat_ZZ solution;
solution.SetDims(1,n);
for (int i=0;i<n;i++)
	{
		solution[0][i]=conv<ZZ>(x[0][i])-conv<ZZ>(xxxx[0][i]);
	}
cout<<"solution="<<solution<<"\n";
}
