#include <iostream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <list>
#include <string>
#include <fstream>

using namespace std;
// Todas as constantes e valores obtidos estão em parsec, anos e massas solares.
int exercicio1 () {
	int N = 20000;
	double xA_i, xA, xA_f;
	double yA_i, yA, yA_f;
        double zA_i, zA, zA_f;
        double xB_i, xB, xB_f;
	double yB_i, yB, yB_f;
        double zB_i, zB, zB_f; 
    	double vA_x, vA_y, vA_z;
    	double vB_x, vB_y, vB_z;
    	double aA_x, aA_y, aA_z;
    	double aB_x, aB_y, aB_z;
    	double f_x, f_y, f_z;
    	long double G = 0.000000000000004398499024; 
    	double delta_t = 20000; 
    	double r;
    	double Ep[N] = {0}, Ec[N] = {0};
    	xA_i = -1;
    	yA_i = 0;
    	zA_i = 0;
    	xB_i = 1;
    	yB_i = 0;
    	zB_i = 0;
    	vA_x = 0;
    	vA_y = 2*pow(10,-6)*1.02269032*pow(10,-2); 
    	vA_z = 0;
    	vB_x = 0;
    	vB_y = -2*pow(10,-6)*1.02269032*pow(10,-2); 
    	vB_z = 0;
    	for (int pos=0; pos<N; pos++) {
        	double vA_i = sqrt(pow(vA_x,2) + pow(vA_y,2) + pow(vA_z,2));
        	double vB_i = sqrt(pow(vB_x,2) + pow(vB_y,2) + pow(vB_z,2));
        	double vCM_x = 0.5*(vA_x + vB_x);
        	double vCM_y = 0.5*(vA_y + vB_y);
        	double vCM_z = 0.5*(vA_z + vB_z);
        	double vCM = sqrt(vCM_x*vCM_x + vCM_y*vCM_y + vCM_z*vCM_z);
        	r = sqrt(pow((xA_i-xB_i),2.0) + pow((yA_i-yB_i),2.0) + pow((zA_i-zB_i),2.0));
        	Ep[pos] = -(G/r);  
        	Ec[pos] = 0.5*2*vCM*vCM + 0.5*(pow(vA_i+vB_i,2));
        	f_x = (G*(xA_i - xB_i)) / pow(r,3);
        	f_y = (G*(yA_i - yB_i)) / pow(r,3);
        	f_z = (G*(zA_i - zB_i)) / pow(r,3);
        	aA_x = -f_x;
        	aA_y = -f_y;
        	aA_z = -f_z;
        	aB_x = f_x;
        	aB_y = f_y;
        	aB_z = f_z;
        	xA = xA_i + vA_x * delta_t;
        	yA = yA_i + vA_y * delta_t; 
        	zA = zA_i + vA_z * delta_t; 
        	xB = xB_i + vB_x * delta_t; 
        	yB = yB_i + vB_y * delta_t; 
        	zB = zB_i + vB_z * delta_t; 
        	xA_f = 2*xA - xA_i + aA_x*pow(delta_t,2);
        	yA_f = 2*yA - yA_i + aA_y*pow(delta_t,2);
        	zA_f = 2*zA - zA_i + aA_z*pow(delta_t,2);
        	xB_f = 2*xB - xB_i + aB_x*pow(delta_t,2);
        	yB_f = 2*yB - yB_i + aB_y*pow(delta_t,2);
        	zB_f = 2*zB - zB_i + aB_z*pow(delta_t,2);
        	vA_x = (xA_f - xA_i)/(2.0*delta_t);
        	vA_y = (yA_f - yA_i)/(2.0*delta_t);
        	vA_z = (zA_f - zA_i)/(2.0*delta_t);
        	vB_x = (xB_f - xB_i)/(2.0*delta_t);
        	vB_y = (yB_f - yB_i)/(2.0*delta_t);
        	vB_z = (zB_f - zB_i)/(2.0*delta_t);
        	xA_i = xA;
        	yA_i = yA;
        	zA_i = zA;
        	xB_i = xB;
        	yB_i = yB;
        	zB_i = zB;
        	cout << pos << "\t" << Ec[pos] << "\t" << Ep[pos] << "\t" << Ep[pos]+Ec[pos] << endl;

    }
	return 0;

}
int exercicio2(){
	string myText;
	ifstream datax("100x.txt");//alterar o número de cada ficheiro para obter para 10,100 e 500 estrelas;
	ifstream datay("100y.txt");
	ifstream dataz("100z.txt");
	ifstream datavx("100vx.txt");
	ifstream datavy("100vy.txt");
	ifstream datavz("100vz.txt");
	ifstream datam("100m.txt");
	vector<double> listax;
	vector<double> listay;
	vector<double> listaz;
	vector<double> listavx;
	vector<double> listavy;
	vector<double> listavz;
	vector<double> listam;
	while(getline(datax, myText))
	{
		listax.push_back(stod(myText));
	}
	while(getline(datay, myText))
	{
		listay.push_back(stod(myText));
	}
	while(getline(dataz, myText))
	{
		listaz.push_back(stod(myText));
	}
	while(getline(datavx, myText))
	{
		listavx.push_back(stod(myText));
	}
	while(getline(datavy, myText))
	{
		listavy.push_back(stod(myText));
	}
	while(getline(datavz, myText))
	{
		listavz.push_back(stod(myText));
	}
	while(getline(datam, myText))
	{
		listam.push_back(stod(myText));
	}
	int N=listax.size();
	int t=5*pow(10,4);
	int Nint=200;
	vector<double> X_last = listax;
	vector<double> Y_last = listay;
	vector<double> Z_last = listaz;
	vector<double> X_now (N,0);
	vector<double> Y_now (N,0);
	vector<double> Z_now (N,0);
	for(int i=0; i < N; i++)
	{
		X_now[i] = X_last[i] + listavx[i]*t;
		Y_now[i] = Y_last[i] + listavy[i]*t;
		Z_now[i] = Z_last[i] + listavz[i]*t;
	}
	double rt = 5.0;
	double c1=0;
	double c2=0;
	for (int i=0;i<N;i++){
		double r1=sqrt((X_last[i]*X_last[i])+(Y_last[i]*Y_last[i])+(Z_last[i]*Z_last[i]));
		double r2=sqrt((X_now[i]*X_now[i])+(Y_now[i]*Y_now[i])+(Z_now[i]*Z_now[i]));
		if (r1<rt){
			c1+=1;}
		if (r2<rt){
			c2+=1;}
	}
	//cout<<c1<<endl;
	//cout<<c2<<endl;
	long double G= 4.49390350312651*pow(10,-15);
	double x1; double y2; double z3;
	for(int iteracao=1; iteracao < Nint; iteracao++){
		vector<double> ax (N,0);
		vector<double> ay (N,0);
		vector<double> az (N,0);
		int c3 = 0;
		for(int j =0; j < N; j++){
			for(int k=0; k<N; k++){
				if(k!=j){
					x1=(X_now[j]-X_now[k])*(X_now[j]-X_now[k]);
					y2=(Y_now[j]-Y_now[k])*(Y_now[j]-Y_now[k]);
					z3=(Z_now[j]-Z_now[k])*(Z_now[j]-Z_now[k]);
					double r = sqrt(x1+y2+z3);
					double a = (G*listam[j])/(r*r);
					ax[j] += a*((X_now[k] - X_now[j])/r);
					ay[j] += a*((Y_now[k] - Y_now[j])/r);
					az[j] += a*((Z_now[k] - Z_now[j])/r);
				}			
			}
		}
		for (int b=0; b<N; b++){
			double X_next = -X_last[b]+2*X_now[b]+ax[b]*pow(t,2);
			double Y_next = -Y_last[b]+2*Y_now[b]+ay[b]*pow(t,2);
			double Z_next = -Z_last[b]+2*Z_now[b]+az[b]*pow(t,2);
			X_last[b] = X_now[b];
			Y_last[b] = Y_now[b];
			Z_last[b] = Z_now[b];
			X_now[b] = X_next;
			Y_now[b] = Y_next;
			Z_now[b] = Z_next;
			double raioEstrela=sqrt((X_now[b]*X_now[b])+(Y_now[b]*Y_now[b])+(Z_now[b]*Z_now[b]));
			if (raioEstrela<rt){
				c3 += 1;}

		}
		cout<<(iteracao+1)*t<<"\t"<<c3<<endl;}
	return 0;}
int exercicio3(){
	int N=100;
	int t=5*pow(10,4);
	int Nint=200;
	int Nconfigs=100;
	vector <double> Media(Nint,0);
	double x1; double y2; double z3;
	long double G= 4.49390350312651*pow(10,-15);
	for (int config = 0; config < Nconfigs; config++){
		vector<double> listavx (N,0);
		vector<double> listavy (N,0);
		vector<double> listavz (N,0);
		vector<double> listam (N,0);
		vector<double> X_last (N,0);
		vector<double> Y_last (N,0);
		vector<double> Z_last (N,0);
		vector<double> X_now (N,0);
		vector<double> Y_now (N,0);
		vector<double> Z_now (N,0);
		double r=0;
		for(int i = 0; i < N; i++){
			bool a=false;
			while(a==false){
				double u=drand48();
				double x=drand48();
				if (u<=(((asinh(x) - x/sqrt(26))*0.7508071177*N)/1.33186)){
					a=true;
					r=u;}
			}
			double theta = rand()% 2*M_PI;
			double phi = acos(2*drand48() - 1);
			X_last[i]=r*cos(theta)*sin(phi);
			Y_last[i]=r*sin(theta)*sin(phi);
			Z_last[i]=r*cos(phi);
			double o=drand48();
			double k=0.1965943953; //constante de normalização
			double base;
			double expoente;
			if (o <= 0.0367523392){
				base = ((o/k) + 0.0568725)* (1/1.42857);
				expoente = double(10)/double(7);}
			else if (o > 0.0367523392 && o <= 0.6280158809){
				base = (7.298295 - (o/k))*(1/3.333);
				expoente = double(-10)/double(3);}
			else if (o > 0.6280158809){
				base = (5.088545 - (o/k))*(1/0.769231);
				expoente = double(-10)/double(13);}
			listam[i]=pow(base,expoente);
			if (i%2!=0){
				long double m=1/(3.085678 * pow(10, 16));
				long double s=3.16887646 * pow(10,-8);
				double conv = pow(10,3)*m*pow(s,-1);
				double a1=drand48();
				double a2=drand48();
				double v1 = sqrt(-log(1-a2))*conv*sin(2*M_PI*a1);
				double v2 = sqrt(-log(1-a2))*conv*cos(2*M_PI*a1);
				double v1x = v1*cos(theta)*sin(phi);
				double v1y = v1*sin(theta)*sin(phi);
				double v1z = v1*v1*cos(phi);
				double v2x = v2*cos(theta)*sin(phi);
				double v2y = v2*sin(theta)*sin(phi);
				double v2z = v2*v1*cos(phi);
				listavx[i-1]=v1x*cos(theta)*sin(phi);
				listavx[i]=v2x*cos(theta)*sin(phi);
				listavy[i-1]=v1y*sin(theta)*sin(phi);
				listavy[i]=v2y*sin(theta)*sin(phi);
				listavz[i-1]=v1z*cos(phi);
				listavz[i]=v2z*cos(phi);
			}
		}
		for(int i=0; i < N; i++){
                	X_now[i] = X_last[i] + listavx[i]*t;
                	Y_now[i] = Y_last[i] + listavy[i]*t;
                	Z_now[i] = Z_last[i] + listavz[i]*t;}
        	double rt = 5.0;
        	double c1=0;
        	double c2=0;
        	for (int i=0;i<N;i++){
                	double r1=sqrt((X_last[i]*X_last[i])+(Y_last[i]*Y_last[i])+(Z_last[i]*Z_last[i]));
                	double r2=sqrt((X_now[i]*X_now[i])+(Y_now[i]*Y_now[i])+(Z_now[i]*Z_now[i]));
                	if (r1<rt){
                        	c1+=1;}
			if (r2<rt){
                        	c2+=1;}
        	}
		Media[0]+=c1;
		Media[1]+=c2;
		for(int iteracao=2; iteracao < Nint; iteracao++){
                	vector<double> ax (N,0);
                	vector<double> ay (N,0);
                	vector<double> az (N,0);
                	int c3 = 0;
                	for(int j =0; j < N; j++){
                        	for(int k=0; k<N; k++){
                                	if(k!=j){
                                        	x1=(X_now[j]-X_now[k])*(X_now[j]-X_now[k]);
                                        	y2=(Y_now[j]-Y_now[k])*(Y_now[j]-Y_now[k]);
                                        	z3=(Z_now[j]-Z_now[k])*(Z_now[j]-Z_now[k]);
                                        	double r = sqrt(x1+y2+z3);
                                        	double a = (G*listam[j])/(r*r);
                                        	ax[j] += a*((X_now[k] - X_now[j])/r);
                                        	ay[j] += a*((Y_now[k] - Y_now[j])/r);
                                        	az[j] += a*((Z_now[k] - Z_now[j])/r);
                                }
                        }
                }
		for (int b=0; b<N; b++){
                        double X_next = -X_last[b]+2*X_now[b]+ax[b]*pow(t,2);
                        double Y_next = -Y_last[b]+2*Y_now[b]+ay[b]*pow(t,2);
                        double Z_next = -Z_last[b]+2*Z_now[b]+az[b]*pow(t,2);
                        X_last[b] = X_now[b];
                        Y_last[b] = Y_now[b];
                        Z_last[b] = Z_now[b];
                        X_now[b] = X_next;
                        Y_now[b] = Y_next;
                        Z_now[b] = Z_next;
                        double raioEstrela=sqrt((X_now[b]*X_now[b])+(Y_now[b]*Y_now[b])+(Z_now[b]*Z_now[b]));
                        if (raioEstrela<rt){
                                c3 += 1;}
                }
		Media[iteracao]+=c3;}
	}
	for (int j=0; j<Nint;j++){
		cout<<t*j<<"\t"<<Media[j]/Nconfigs<<endl;}
	return 0;}
int main(){
	//exercicio1();
	//exercicio2();
	//exercicio3();
}
