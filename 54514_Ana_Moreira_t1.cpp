#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <cmath>
using namespace std;
//1.1.
int main()
{
	int c=16807;
	int ci=3;
	int nmaxi=31;
	int nmax=2147483647;
	int M=3000;
	int X[M];
	int i;
	int Xa[M];
	Xa[0]=1;
	X[0]=1;
	for (i=1; i<M;i++) {
		X[i]=(c*X[i-1])%nmax;
		Xa[i]=(ci*Xa[i-1])%nmaxi;
//		cout<<(double)(X[i])/(double)(nmax)<<endl;	
	}
//Código para os gráficos;
	int j;
	int k=0;
	int l=0;
	for (j=1;j<M;j++) {
		k=j-1;
		l=j-2;
		//cout<<(long double)X[j]/(long double)nmax<<"\t"<<(long double)X[k]/(long double)nmax<<"\t"<<(long double)X[l]/(long double)nmax<<endl;

	}
	int p;
	int Z[M];
	double W[M];
	W[0]=1;
	Z[0]=1;
	for (p=1;p<M;p++){
		Z[p]=rand();
		W[p]=drand48();
		//cout<<(long double)Z[p]/(long double)RAND_MAX<<endl;
		//cout<<(long double)(W[p])<<endl;
	}
	int t;
	int s=0;
	int q=0;
	for (t=1;t<M;t++){
		s=t-1;
                q=t-2;
                //cout<<(long double)Z[t]/(long double)RAND_MAX<<"\t"<<(long double)Z[s]/(long double)RAND_MAX<<"\t"<<(long double)Z[q]/(long double)RAND_MAX<<endl;

	}


//1.2.
	float r,theta,x,y;
	int Y[M];
	int n;
	for (n=1;n<M;n++){
		r=(double)rand()/(double)RAND_MAX;
		theta=(double)(2*M_PI*rand())/(double)RAND_MAX;
		x=(double)sqrt(r)*(double)sin(theta);
		y=(double)sqrt(r)*(double)cos(theta);
//		cout<<x<<"\t"<<y<<endl;
	}



//1.3.
	const int inter=5;
	const int T=1000;
	double prob=1/double(inter);
	int Cont[inter]={0,0,0,0,0};
	double ni=T*prob;
	long double qui=0;
	double Xn;
	double tamanho=0.2;
	int m;
	for (m=1;m<T;m++){
		Xn=abs((long double)(X[m])/(long double)(nmax));
		for (i=1;i<=inter;i++){
			if ((double)Xn<double((i*tamanho)) && double(Xn)>=double((i-1)*tamanho)){ 
				Cont[i-1]=Cont[i-1]+1;    
				//cout<<Cont[0]<<"\t"<<Cont[1]<<"\t"<<Cont[2]<<"\t"<<Cont[3]<<"\t"<<Cont[4]<<endl;
			}
		}
	}
        int Cont1[inter]={0,0,0,0,0};
        long double qui1=0;
        double Xna;
        int mi;
        for (m=1;m<T;m++){
                Xna=abs((long double)(Xa[m])/(long double)(nmaxi));
                for (i=1;i<=inter;i++){
                        if ((double)Xna<double((i*tamanho)) && double(Xna)>=double((i-1)*tamanho)){
                                Cont1[i-1]=Cont1[i-1]+1;
                                //cout<<Cont1[0]<<"\t"<<Cont1[1]<<"\t"<<Cont1[2]<<"\t"<<Cont1[3]<<"\t"<<Cont1[4]<<endl;
                        }
                }
        }
        int Cont2[inter]={0,0,0,0,0};
        long double qui2=0;
        double Xnj;
        int o;
        for (o=1;o<T;o++){
                Xnj=abs(double(rand())/RAND_MAX);
                for (i=1;i<=inter;i++){
                        if ((double)Xnj<double((i*tamanho)) && double(Xnj)>=double((i-1)*tamanho)){
                                Cont2[i-1]=Cont2[i-1]+1;
                               // cout<<Cont2[0]<<"\t"<<Cont2[1]<<"\t"<<Cont2[2]<<"\t"<<Cont2[3]<<"\t"<<Cont2[4]<<endl;
                        }
                }
        }
	for (j=1;j<inter;j++){
		qui=qui+((Cont[j]-ni)*(Cont[j]-ni))/(ni);
		qui2=qui2+((Cont2[j]-ni)*(Cont2[j]-ni))/(ni);
		qui1=qui1+((Cont1[j]-ni)*(Cont1[j]-ni))/(ni);
	}
	cout<<qui<<"\t"<<qui2<<"\t"<<qui1<<endl;
}
