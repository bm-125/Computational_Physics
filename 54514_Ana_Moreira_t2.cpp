#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "latticeview.h"
#define ImageWidth 1000
#define ImageHeight 1000

using namespace std;

int main() {
	for(int k=1;k<=500;k++){
        	int lx, ly, i, j, site;	
        	lx = 15;
        	ly = 15;
        	int lat[lx*ly]; 
        	double p, E;
		int passo=0;
		int passop=0;
		int perc=0;
		int topo=0;
		int media=0;
		int num=0;
       		p = 0.002*k; //alterar para obter os diferentes agregados;
		for(int n=1;n<=10;n++){
        // 0- livre; 1- ocupado; 2- a queimar; 3- queimado; 4- queimar na próxima iteração
			for (i=0; i<lx; i++) {
                		for (j=0; j<ly; j++) {
                        		E= drand48();
                        		site= i+j*lx;
                        		i= site%lx;
                        		j= site/lx;
                        		if (E<p){
                                		lat[site]= 1;}
                        		else{
                                		lat[site]= 0;}
                		} 
        		}

        		//Print_lattice (lat, lx, ly, ImageWidth, ImageHeight, "2_1.ppm"); //2.1.

        		int t, q, g;
        		bool queima=false;
			bool prim=true;

        		for (t=0; t<lx*ly; t++) {
				if (drand48()<p)
					lat[t]= 1;
			}

			for (q=0; q<lx; q++) {
				if (lat[q]==1){
					lat[q]= 2;}
					if (lat[q+lx]==1){
						lat[q+lx]= 3;}
					queima= true;
			}

        		//Print_lattice (lat, lx, ly, ImageWidth, ImageHeight, "2_2_1.ppm");
	
		        while (queima) {
				queima= false;		
				for (g=0; g<lx*ly; g++) {
					if (lat[g]==2){
						lat[g]= 4;}
					if (lat[g]==3){
						lat[g]= 2;
						queima= true;}
					if (lat[g]==1){
						if ((lat[g-1]==2) && (g%lx!=0)){
							lat[g]= 3;
							queima= true;}
						if ((lat[g+1]==2) && (g%lx!=lx-1)){
							lat[g]= 3;
							queima= true;}
						if ((lat[g+lx]==2) && (g<(lx*ly)-lx)){
							lat[g]= 3;
							queima= true;}
						if (lat[g-lx]==2){
							lat[g]= 3;
							queima= true;}}
					if ((g>((lx*ly)-lx-1)) && (lat[g]==2) && (prim==true)){
						passop= passo+1;
						prim= false;}
		}
				passo+=1;
				for (i=lx*lx-lx;i<lx*lx;i++){
					if(passop==0)
						if (lat[i]==2)
							{passop=passo;}}
		}
        		if(passop!=0){
				perc=perc+1;
				media+=passop;}	
			cout<<p<<"\t"<<passo<<"\t"<<passop<<"\t"<<perc<<"\t"<<media<<endl;		}
	}
}

