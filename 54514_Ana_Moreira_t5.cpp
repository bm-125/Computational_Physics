#include <iostream>
#include "sparselib.hh"
#include <cmath>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <fstream>
#include <string.h>
#include <list>

using namespace ::std;
using namespace ::sparselib_load;
using ::sparselib::index_type ;
using ::sparselib_fun::maxval ;
using ::sparselib_fun::minval ;
using ::sparselib_fun::absval ;

long double calor(long double x){
	long double a=(exp(-x)*(-exp(2*x)+exp(20)))/(exp(20)-1);
	return abs(a);}
int main(){
	int L = 10;
	double vals[] = {0.1,0.5,1};
	for (int i=0; i<4;i++){
		Vector<double>x,b;
		CCoorMatrix<double> A;
		double deltax=vals[i];
		int dimensao=L/deltax;
		int nzero=3*dimensao-2;
		x.new_dim(dimensao);
		b.new_dim(dimensao);
		A.new_dim(dimensao, dimensao, nzero);
		b=0;
		b[0]=(deltax*deltax-6)/(6*deltax);
		//cout<<b<<"\t"<<endl;
		for (int i=0; i<dimensao; i++)
         		A.insert(i,i) = - (6 + 2 * deltax * deltax) / (3 * deltax);
		   	A.insert(0,1) = - (deltax * deltax - 6) / (6 * deltax);
		   	A.insert(dimensao-1,dimensao-2) = - ( deltax * deltax - 6) / (6 * deltax);
		for (i=1; i<dimensao-1; i++) {
         		A.insert(i,i-1) = - (deltax * deltax - 6) / (6 * deltax);
		   	A.insert(i,i+1) = - (deltax * deltax - 6) / (6 * deltax);
		}
		//cout<<A<<"\n";
		x=0;
		IdPreco<int> P;
		P.build(A); 
		int iter, max_iter = 10000; 
		double epsi = 0.0001; 
		x = 0;
		double res = bicgstab(A,b,x,P,epsi,max_iter,iter); 
		//cout << x << endl;
		vector <double> dist_sum_vec;
		vector <double> calor_vec;
		long double dist_sum = 0;
		for (int i=0; i<max_iter; i++) {
			long double ai = calor(dist_sum);
			calor_vec.push_back (ai);
			dist_sum_vec.push_back (dist_sum);
			dist_sum += (epsi / 0.1);
		}
		for (int i=0;i<dist_sum_vec.size();i++){
			//cout<<dist_sum_vec[i]<<"\t"<<calor_vec[i]<<endl;
		}
		vector <double> n_dist_vec;
		vector <double> counter_vec;
		n_dist_vec.push_back (1);
		counter_vec.push_back (0);
		if (deltax==0.1){
                        for (int i=0; i<x.size(); i++) {
                                n_dist_vec.push_back (abs (x[i]));
                                double j = (double) i / 10;
                                counter_vec.push_back (j + 0.1);}
                        for (int i=0;i<n_dist_vec.size();i++){
                                //cout<<counter_vec[i]<<"\t"<<n_dist_vec[i]<<endl;
                        }
                }
		if (deltax==0.5){
                        for (int i=0; i<x.size(); i++) {
                                n_dist_vec.push_back (abs (x[i]));
                                double j = (double) i / 2;
                                counter_vec.push_back (j + 0.5);}
                        for (int i=0;i<n_dist_vec.size();i++){
                                //cout<<counter_vec[i]<<"\t"<<n_dist_vec[i]<<endl;
			}
		}
		if (deltax==1){
                        for (int i=0; i<x.size(); i++) {
                                n_dist_vec.push_back (abs (x[i]));
                                counter_vec.push_back (i + 1);}
                        for (int i=0;i<n_dist_vec.size();i++){
                                //cout<<counter_vec[i]<<"\t"<<n_dist_vec[i]<<endl;
			}
		}
		
}
		
}

