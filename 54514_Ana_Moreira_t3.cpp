#include<iostream>
#include<stdlib.h>
#include<vector>
#include<math.h>
#include<list>
#include<string>
#include<fstream>
#include<cmath>
#include<sstream>
#include<algorithm> 

using namespace std;

int main(){
//3.1
	double mf=3190;
	double mr=6.09;
	double mc=0.00383;
	double m=mf+mr+mc;
	int n=100000;
	int nbins=2000;
	int c=0;
	float dx=pow(10,-8);
	int tamanho=100000;
	double N[tamanho]={0.0};
	int i;
	vector<double>X;
	vector<double>P;
	int M[tamanho]={0};
	for (i=0;i<n;i++){
		long double r=drand48();
                long double x=(-1/m)*(log(1-r));
		long double p=exp(-m*x);
		X.push_back(x);
		P.push_back(p);
	}

	for (c=0;c<n;c++){
		int k=(int)(X[c]/dx);
		if (k<tamanho-1){
			N[k]++;}
		else{
			N[tamanho-1]++;}
	}
	M[0]=n-N[0];
	for (int a=1; a<tamanho;a++){
			M[a]=M[a-1]-N[a];
		}
	for (int w=0;w<tamanho-1;w++){
		//cout<<w*dx<<"\t"<<double(M[w])/double(n)<<endl;
	}

//3.2
	ifstream filename;
	filename.open("data.txt");
        double data1; double data2; double data3;
        vector<double> Vec;
	vector<double> Vec1;
	vector<double>Vec2;
        while (filename>>data1>>data3){
                Vec.push_back(data1);
		Vec2.push_back(data3);}
        filename.close();
	ifstream filename2;
	filename2.open("probabilidades.txt");
	while (filename2>>data2){
		Vec1.push_back(data2);}
	filename2.close();
        long double soma=0;
        for (int i=0;i<Vec1.size();i++){
                soma+=Vec1[i];}

        long double MC=(Vec1.back()-Vec1[0])*soma/Vec1.size();
        for(int i=0;i<Vec.size();i++){
                Vec1[i]=Vec1[i]/MC;
	}
	
        vector<double>CDF1=Vec;
        vector<double>CDF2=Vec1;
	sort(CDF2.begin(),CDF2.end());
	double p=drand48();
	double densidade=2.69890;
	double total=0;
	for(int i=0;i<99;i++){
		if(p<=CDF2[i]+total){break;}
		total+=CDF2[i];
	}
        double hist[10000]={0.0};
	double miut=0;
	for (int i=0;i<Vec2.size();i++){
		miut+=Vec2[i];}
	double miut1=double(miut/Vec2.size()*densidade);
	for(int q=0;q<2000;q++){
		double na=drand48();
		double x=(-1/m)*(log(1-na));
		int k=(int)(x/dx);
                if (k<=10000-1){
                        hist[k]++;}
                else{
                        hist[10000-1]++;}

         	}
	int Ma[10000]={0};
	for(int a=0;a<10000;a++){
		if (a==0){
			Ma[a]=n-hist[a];}
		else{
			Ma[a]=Ma[a-1]-hist[a];}	
                //cout<<a*dx<<"\t"<<log(double(Ma[a])/double(n))<<endl;}

}
	

	
