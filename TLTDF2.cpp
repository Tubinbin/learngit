#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <NTL/mat_RR.h>
#include <iostream>
#include <random>
#include <time.h>

using namespace std;
using namespace NTL;

//分享陷门函数,用户输入身份，函数输出用户的私钥分量，最好是输出整个私钥，但是这个和参数w相关（可以写两个函数共同完成该功能）
//需要定义系数矩阵，输入变量（身份）,多项式的最高次幂（T-1）
mat_ZZ_p sharefun1(mat_ZZ_p coff_share0, mat_ZZ_p coff_share1, mat_ZZ_p coff_share2, mat_ZZ_p coff_share3, mat_ZZ_p coff_share4, 
	mat_ZZ_p coff_share5, mat_ZZ_p coff_share6, mat_ZZ_p coff_share7, mat_ZZ_p coff_share8, mat_ZZ_p coff_share9, mat_ZZ_p coff_share10, mat_ZZ_p coff_share11, mat_ZZ_p coff_share12, mat_ZZ_p coff_share13, mat_ZZ_p coff_share14, 
	mat_ZZ_p coff_share15, mat_ZZ_p coff_share16, mat_ZZ_p coff_share17, mat_ZZ_p coff_share18, mat_ZZ_p coff_share19, long T_share, long id_share, 
	long d_share, long w_share)
{
	int i,j;
	mat_ZZ_p id_matrix;//定义身份向量；
	id_matrix.SetDims(1,T_share);
	for(i=0; i<T_share; i++)
		{
			id_matrix[0][i]=power(conv<ZZ_p>(conv<ZZ>(id_share)), i);
			//id_matrix[0][i]=conv<ZZ_p>(conv<ZZ>(power(id_share, i));
		}
	mat_ZZ_p td_share, td_share_transpose;
	td_share_transpose.SetDims(w_share,d_share);
	//vector(1*d)=vector(1*T)*matrix(T*d)
	for (i=0;i<d_share;i++)
		{
			for(j=0;j<T_share;j++)
				{
					td_share_transpose[0][i]=td_share_transpose[0][i]+id_matrix[0][j]*transpose(coff_share0)[j][i];
				}
		}

	for (i=0;i<d_share;i++)
	{
		for(j=0;j<T_share;j++)
			{
				td_share_transpose[1][i]=td_share_transpose[1][i]+id_matrix[0][j]*transpose(coff_share1)[j][i];
			}
	}
	for (i=0;i<d_share;i++)
		{
			for(j=0;j<T_share;j++)
				{
					td_share_transpose[2][i]=td_share_transpose[2][i]+id_matrix[0][j]*transpose(coff_share2)[j][i];
				}
		}
	for (i=0;i<d_share;i++)
		{
			for(j=0;j<T_share;j++)
				{
					td_share_transpose[3][i]=td_share_transpose[3][i]+id_matrix[0][j]*transpose(coff_share3)[j][i];
				}
		}

	for (i=0;i<d_share;i++)
		{
			for(j=0;j<T_share;j++)
				{
					td_share_transpose[4][i]=td_share_transpose[4][i]+id_matrix[0][j]*transpose(coff_share4)[j][i];
				}
		}

	for (i=0;i<d_share;i++)
		{
			for(j=0;j<T_share;j++)
				{
					td_share_transpose[5][i]=td_share_transpose[5][i]+id_matrix[0][j]*transpose(coff_share5)[j][i];
				}
		}

	for (i=0;i<d_share;i++)
		{
			for(j=0;j<T_share;j++)
				{
					td_share_transpose[6][i]=td_share_transpose[6][i]+id_matrix[0][j]*transpose(coff_share6)[j][i];
				}
		}

	for (i=0;i<d_share;i++)
		{
			for(j=0;j<T_share;j++)
				{
					td_share_transpose[7][i]=td_share_transpose[7][i]+id_matrix[0][j]*transpose(coff_share7)[j][i];
				}
		}

	for (i=0;i<d_share;i++)
		{
			for(j=0;j<T_share;j++)
				{
					td_share_transpose[8][i]=td_share_transpose[8][i]+id_matrix[0][j]*transpose(coff_share8)[j][i];
				}
		}
	for (i=0;i<d_share;i++)
		{
			for(j=0;j<T_share;j++)
				{
					td_share_transpose[9][i]=td_share_transpose[9][i]+id_matrix[0][j]*transpose(coff_share9)[j][i];
				}
		}
	for (i=0;i<d_share;i++)
		{
			for(j=0;j<T_share;j++)
				{
					td_share_transpose[10][i]=td_share_transpose[10][i]+id_matrix[0][j]*transpose(coff_share10)[j][i];
				}
		}

	for (i=0;i<d_share;i++)
	{
		for(j=0;j<T_share;j++)
			{
				td_share_transpose[11][i]=td_share_transpose[11][i]+id_matrix[0][j]*transpose(coff_share11)[j][i];
			}
	}
	for (i=0;i<d_share;i++)
		{
			for(j=0;j<T_share;j++)
				{
					td_share_transpose[12][i]=td_share_transpose[12][i]+id_matrix[0][j]*transpose(coff_share12)[j][i];
				}
		}
	for (i=0;i<d_share;i++)
		{
			for(j=0;j<T_share;j++)
				{
					td_share_transpose[13][i]=td_share_transpose[13][i]+id_matrix[0][j]*transpose(coff_share13)[j][i];
				}
		}

	for (i=0;i<d_share;i++)
		{
			for(j=0;j<T_share;j++)
				{
					td_share_transpose[14][i]=td_share_transpose[14][i]+id_matrix[0][j]*transpose(coff_share14)[j][i];
				}
		}

	for (i=0;i<d_share;i++)
		{
			for(j=0;j<T_share;j++)
				{
					td_share_transpose[15][i]=td_share_transpose[15][i]+id_matrix[0][j]*transpose(coff_share15)[j][i];
				}
		}

	for (i=0;i<d_share;i++)
		{
			for(j=0;j<T_share;j++)
				{
					td_share_transpose[16][i]=td_share_transpose[16][i]+id_matrix[0][j]*transpose(coff_share16)[j][i];
				}
		}

	for (i=0;i<d_share;i++)
		{
			for(j=0;j<T_share;j++)
				{
					td_share_transpose[17][i]=td_share_transpose[17][i]+id_matrix[0][j]*transpose(coff_share17)[j][i];
				}
		}

	for (i=0;i<d_share;i++)
		{
			for(j=0;j<T_share;j++)
				{
					td_share_transpose[18][i]=td_share_transpose[18][i]+id_matrix[0][j]*transpose(coff_share18)[j][i];
				}
		}
	for (i=0;i<d_share;i++)
		{
			for(j=0;j<T_share;j++)
				{
					td_share_transpose[19][i]=td_share_transpose[19][i]+id_matrix[0][j]*transpose(coff_share19)[j][i];
				}
		}

	//totel 10 rows, get a 10*d matrix
	//after transposing, get a d*10 matrix
	td_share=transpose(td_share_transpose);
	//cout<<"td_share="<<td_share<<"\n";
	return td_share;
}


//为了方便，直接先生成w个随机的系数矩阵，因为每个随机矩阵的第一列构成了陷门，而且都是Z_q上的随机选择

int main()
{
clock_t tbegin, tend;
int i,j,k; 
long n, d, w, m, p, N;
long T, eta;//n代表加密用的随机串的长度；d代表安全参数，即私钥矩阵的行数；w代表私钥矩阵的列数；
//m代表m*l=n；p消息矩阵的模数；N用户总数；T解密至少需要的用户数； eta代表N^2(error) (N!)^2；
ZZ q;
n=100, d=50, w=20, m=20, p=37,q=5971847, N=3, T=3; 
int eta2=1;
for (i=1;i<=N;i++)
	{
		eta2=eta2*i;
	}

eta=eta2*eta2;

ZZ_p::init(q);

//思考生成随机矩阵的数量随输入w的变化而变化
mat_ZZ_p A, Z, ZZZ, B0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15, B16, B17, B18, B19;
//(pk) random A and (sk) random Z 
tbegin=clock();
random(A, n, d);
Z.SetDims(d, w);
ZZZ=transpose(Z);
random(B0, d, T);
random(B1, d, T);
random(B2, d, T);
random(B3, d, T);
random(B4, d, T);
random(B5, d, T);
random(B6, d, T);
random(B7, d, T);
random(B8, d, T);
random(B9, d, T);
random(B10, d, T);
random(B11, d, T);
random(B12, d, T);
random(B13, d, T);
random(B14, d, T);
random(B15, d, T);
random(B16, d, T);
random(B17, d, T);
random(B18, d, T);
random(B19, d, T);
cout<<"B0="<<transpose(B0)<<"\n";
cout<<"B1="<<transpose(B1)<<"\n";
cout<<"B2="<<transpose(B2)<<"\n";
cout<<"B3="<<transpose(B3)<<"\n";
cout<<"B4="<<transpose(B4)<<"\n";
cout<<"B5="<<transpose(B5)<<"\n";
cout<<"B6="<<transpose(B6)<<"\n";
cout<<"B7="<<transpose(B7)<<"\n";
cout<<"B8="<<transpose(B8)<<"\n";
cout<<"B9="<<transpose(B9)<<"\n";
cout<<"B10="<<transpose(B10)<<"\n";
cout<<"B11="<<transpose(B11)<<"\n";
cout<<"B12="<<transpose(B12)<<"\n";
cout<<"B13="<<transpose(B13)<<"\n";
cout<<"B14="<<transpose(B14)<<"\n";
cout<<"B15="<<transpose(B15)<<"\n";
cout<<"B16="<<transpose(B16)<<"\n";
cout<<"B17="<<transpose(B17)<<"\n";
cout<<"B18="<<transpose(B18)<<"\n";
cout<<"B19="<<transpose(B19)<<"\n";


//提取w个随机矩阵的第一列组成陷门矩阵Z；用到NTL的矩阵转置transpose；因为NTL可以直接表示行；
ZZZ[0]=transpose(B0)[0];
ZZZ[1]=transpose(B1)[0];
ZZZ[2]=transpose(B2)[0];
ZZZ[3]=transpose(B3)[0];
ZZZ[4]=transpose(B4)[0];
ZZZ[5]=transpose(B5)[0];
ZZZ[6]=transpose(B6)[0];
ZZZ[7]=transpose(B7)[0];
ZZZ[8]=transpose(B8)[0];
ZZZ[9]=transpose(B9)[0];
ZZZ[10]=transpose(B10)[0];
ZZZ[11]=transpose(B11)[0];
ZZZ[12]=transpose(B12)[0];
ZZZ[13]=transpose(B13)[0];
ZZZ[14]=transpose(B14)[0];
ZZZ[15]=transpose(B15)[0];
ZZZ[16]=transpose(B16)[0];
ZZZ[17]=transpose(B17)[0];
ZZZ[18]=transpose(B18)[0];
ZZZ[19]=transpose(B19)[0];

Z=transpose(ZZZ);

//定义N个用户的私钥矩阵
mat_ZZ_p ZZ0, ZZ1, ZZ2, ZZ3, ZZ4;

long id0, id1, id2, id3, id4;
id0=1;id1=2;id2=3;id3=4;id4=5;
//根据share函数求分享的矩阵
//分享陷门，计算N个用户各自的私钥ZZ1,ZZ2,ZZ3,...,ZZN；
ZZ0=sharefun1(B0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15, B16, B17, B18, B19,T, id0, d, w);
ZZ1=sharefun1(B0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15, B16, B17, B18, B19,T, id1, d, w);
ZZ2=sharefun1(B0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15, B16, B17, B18, B19,T, id2, d, w);
ZZ3=sharefun1(B0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15, B16, B17, B18, B19,T, id3, d, w);
ZZ4=sharefun1(B0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15, B16, B17, B18, B19,T, id4, d, w);

cout<<"A="<<A<<"\n";
cout<<"Z="<<Z<<"\n";
cout<<"ZZ0="<<ZZ0<<"\n";
cout<<"ZZ1="<<ZZ1<<"\n";
cout<<"ZZ2="<<ZZ2<<"\n";
cout<<"ZZ3="<<ZZ3<<"\n";
cout<<"ZZ4="<<ZZ4<<"\n";
tend=clock();
double t_1;
t_1=tend-tbegin;//the time is cost by generating the pk sk and shared sk  

mat_ZZ_p A_Z;//def the product of A and Z
mul(A_Z, A, Z);
cout<<"A_Z="<< A_Z <<"\n";

//def the message matrix B
tbegin=clock();
long l=5;//n=m*l
//1,def vector b
Mat<RR> b;
b.SetDims(l, 1);
for (i=0;i<l;i++)
	{
		b[i][0]=pow(2, i);
	}
cout<<"b="<<b<<"\n";

mat_RR I;
mat_ZZ_p B;//def message matrix
ident(I,m);
//Mat<RR> B;
B.SetDims(n, m);

int ii=0;
for (j=0; j<m; j++)
	{
		for (i=0; i<m; i++)
			{
				for (k=0;k<l;k++)
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
cout<<"B="<<B<<"\n";

//def error matrix
double delta;
delta=1.0/(32*37*100);
unsigned seed=time(NULL);
default_random_engine generator (seed);
normal_distribution<double>distribution(0.0,delta);

Mat<ZZ_p> E;
E.SetDims(n, m);
for (i=0;i<n;i++)
	{
		for (j=0;j<m;j++)
			{
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
cout<<"y2=x*C"<<y2<<"\n";
y1=x*A;
cout<<"y1=x*A"<<y1<<"\n";
tend=clock();
double t_2;
t_2=tend-tbegin;
//*************************************************
//*************************************************
//*************************************************
//def inversion algorithm

ZZ_p eta1;
eta1=inv(conv<ZZ_p>(conv<ZZ>(eta)));
mat_ZZ_p aa;
aa=eta1*y1;
cout<<"aa=eta1*y1"<<aa<<"\n";
//下面计算各个用户的求逆分享
mat_ZZ_p inv0, inv1, inv2, inv3, inv4;
mat_ZZ_p inv0_transpose, inv1_transpose, inv2_transpose, inv3_transpose, inv4_transpose;
mat_ZZ_p ee0, ee1,ee2,ee3,ee4;

ee0.SetDims(1, w);
for(i=0;i<w;i++)
	{
		ee0[0][i]=conv<ZZ_p>(ZZ(5971847*distribution(generator)+0.5));
	}
cout <<"ee0="<<ee0<<"\n";

cout <<"aa.NumRows()="<<aa.NumRows()<<"\n";
cout <<"aa.NumCols()="<<aa.NumCols()<<"\n";

cout <<"ZZ0.NumRows()="<<ZZ0.NumRows()<<"\n";
cout <<"ZZ0.NumCols()="<<ZZ0.NumCols()<<"\n";

inv0_transpose=aa*ZZ0+ee0;//行向量
inv0=transpose(inv0_transpose);
cout <<"inv0_transpose=aa*ZZ0+ee0="<<inv0_transpose<<"\n";//列向量
//inv0_transpose=aa*ZZ0;
//cout <<"inv0_transpose=aa*ZZ0"<<inv0_transpose<<"\n";//列向量


ee1.SetDims(1, w);
for (i=0;i<w;i++)
	{
		ee1[0][i]=conv<ZZ_p>(ZZ(5971847*distribution(generator)+0.5));
	}
cout <<"ee1="<<ee1<<"\n";

inv1_transpose=aa*ZZ1+ee1;//行向量
inv1=transpose(inv1_transpose);
cout <<"inv1="<<inv1<<"\n";//列向量

ee2.SetDims(1, w);
for (i=0;i<w;i++)
	{
		ee2[0][i]=conv<ZZ_p>(ZZ(5971847*distribution(generator)+0.5));
	}
cout <<"ee2="<<ee2<<"\n";

inv2_transpose=aa*ZZ2+ee2;//行向量
inv2=transpose(inv2_transpose);
cout <<"inv2="<<inv2<<"\n";//列向量

ee3.SetDims(1, w);
for (i=0;i<w;i++)
	{
		ee3[0][i]=conv<ZZ_p>(ZZ(5971847*distribution(generator)+0.5));
	}
cout <<"ee3="<<ee3<<"\n";

inv3_transpose=aa*ZZ3+ee3;//行向量
inv3=transpose(inv3_transpose);
cout <<"inv3="<<inv3<<"\n";//列向量

ee4.SetDims(1, w);
for (i=0;i<w;i++)
	{
		ee4[0][i]=conv<ZZ_p>(ZZ(5971847*distribution(generator)+0.5));
	}
cout <<"ee4="<<ee4<<"\n";

inv4_transpose=aa*ZZ4+ee4;//行向量
inv4=transpose(inv4_transpose);
cout <<"inv4="<<inv4<<"\n";//列向量

//下面设计组合算法combine
//首先输入门限（T）个用户的身份和解密分享（因为T=3,假设输入的是T0=1,T1=2,T2=3；inv0,inv1,inv2）s
//输入用户身份；
tbegin=clock();
long T0=1,T1=2,T2=3;
mat_ZZ_p ID;//定义一个身份矩阵
ID.SetDims(1, T); 
ID[0][0]=T0;
ID[0][1]=T1;
ID[0][2]=T2;

mat_ZZ_p coff_com;//定义一个系数矩阵
coff_com.SetDims(1, T); 

for(i=0;i<T;i++)
	{
		coff_com[0][i]=1;
		for(j=0;j<T;j++)
			{
				if (i==j) 
					{
						continue;
					}
				coff_com[0][i]=coff_com[0][i]*(-ID[0][j])*inv(ID[0][i]-ID[0][j]);
				//此处注意因为ZZ_p类型下ZZ_p a； a=-3; 输出的a=3,并非p-3；一般需要先把a转换成ZZ类型，进行加减，然后？？？？转换成ZZ_p类型因为
				//ZZ a; a=-3; ZZ_p b; b=conv<ZZ_p>(a); 输出b=p-3;？？？？？？似乎没有问题这里，因为在abotdf中最后使用高斯消元法解方程求v时，出现
				//这个问题；
			}
	}
cout<<"coff_com="<<coff_com<<"\n";
cout<<"eta="<<eta<<"\n";
//由于上面用户的身份是1 2 3，下面也用inv0，inv1，inv2；
//将用户的解密分享组合出结果；
mat_ZZ_p yy;//定义新的矩阵用于存放组合后的密文；
yy=(eta*coff_com[0][0])*inv0_transpose+eta*coff_com[0][1]*inv1_transpose+eta*coff_com[0][2]*inv2_transpose;
cout<<"yy=(eta*coff_com[0][0])*inv0_transpose"<<yy<<"\n";


mat_ZZ_p Z2;

Z2=y2-yy;//相当于密文向量减去A和Z的内积以及部分噪声；
cout<<"Z2=y2-yy"<<Z2<<"\n";

mat_RR t;
t.SetDims(1, m);//此处应该是w，但是m=w，也就不区分了；

for (int i=0;i<m;i++)
{t[0][i]=conv<RR>(conv<ZZ>(Z2[0][i]))/conv<RR>(q);}
cout<<"t="<<t<<"\n";

mat_RR closest, xx;
closest.SetDims(1, p);
xx.SetDims(1, m);

//ZZ temp;

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
for (int i=0; i<p-1;i++)
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
cout<<"xxx="<<xxx<<"\n";
tend=clock();
double t_3;
t_3=tend-tbegin;

//check solution
mat_ZZ solution;
solution.SetDims(1,n);
for (int i=0;i<n;i++)
	{
		solution[0][i]=conv<ZZ>(x[0][i])-conv<ZZ>(xxx[0][i]);
	}
cout<<"solution="<<solution<<"\n";
cout<<"t_1="<<t_1<<"\n";
cout<<"t_2="<<t_2<<"\n";
cout<<"t_3="<<t_3<<"\n";

}