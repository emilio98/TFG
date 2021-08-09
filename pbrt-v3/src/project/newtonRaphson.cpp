#include "project/newtonRaphson.h"
#include "project/ellipticIntegral.h"
#include<cmath>
#include<iostream>

using namespace std;

double omegaDerivative_r2(double Phi, double n, double m, double p, double at_bt){
	double phi = atan(at_bt*tan(Phi));
	double sinphi = sin(phi);
	double sinphi2 = sinphi*sinphi;

	double f = 1/((1-n*sinphi2)*sqrt(1-sinphi2*m));

	double tanPhi = tan(Phi);
	double tanPhi2 = tanPhi*tanPhi;

	double phiDerivative = at_bt * (1+tanPhi2)/(1+at_bt*at_bt*tanPhi2);

	return 1-p*f*phiDerivative;
}

double NR_diskRad(double tol, double n, double m, double p, 
										 double at_bt, double vomegapiover2){
	double xi=M_PI*0.25, f = 20, fd;
	int i = 0, max_it = 10;

	while(abs(f)>tol && i<max_it){
		f=omega_r2(xi,n,m,p,at_bt)-vomegapiover2;
		fd=omegaDerivative_r2(xi,n,m,p,at_bt);
		xi=xi-f/fd;
		i++;
	}	

	xi=0;
	i = 0;
	while(abs(f)>tol && i<max_it){
		f=omega_r2(xi,n,m,p,at_bt)-vomegapiover2;
		fd=omegaDerivative_r2(xi,n,m,p,at_bt);
		xi=xi-f/fd;
		i++;
	}

	xi=M_PI*0.5;
	i = 0;
	while(abs(f)>tol && i<max_it){
		f=omega_r2(xi,n,m,p,at_bt)-vomegapiover2;
		fd=omegaDerivative_r2(xi,n,m,p,at_bt);
		xi=xi-f/fd;
		i++;
	}

	if(abs(f)>tol){
		cout << xi << " xi\n";
		cout << f << " f\n";
		cout << fd << " fd\n";
		cout << vomegapiover2 << " vomegapiover2\n";
	}
	
	return xi;
}

double NR_spherePar(int type, double tol, double vApay, double ax,
										double ay, double xe, double yl){
	double xi=ay*0.5, f = 20, fd;
	int i = 0, max_it = 10;

	if(type == 1){
		while(abs(f) > tol){
			f =Ap1(ax,ay,xi) - vApay;
			fd=Ap1Derivative(ax, ay, xi);
			xi=xi-f/fd;
		}
	}else if (type == 2){
		while(abs(f) > tol){
			f =Ap2(ax,ay,xe,yl,xi) - vApay;
			fd=Ap2Derivative(ax, ay, xe, yl, xi);
			xi=xi-f/fd;
		}
	}else{
		while(abs(f) > tol){
			f =Ap3(ax,ay,xe,yl,xi) - vApay;
			fd=Ap3Derivative(ax, ay, xe, yl, xi);
			xi=xi-f/fd;
		}
	}
	
	return xi;

}