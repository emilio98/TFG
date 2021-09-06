#ifndef NEWTONRAPHSON_H
#define NEWTONRAPHSON_H

#include "project/ellipticIntegral.h"
#include<cmath>

inline double omega_r2(double Phi, double n, double m, double p, double at_bt){
	return Phi-p*elliptic_inc_pim(std::atan(at_bt*std::tan(Phi)),n,m);
}

double omegaDerivative_r2(double Phi, double n, double m, double p, double at_bt);

double NR_diskRad(double tol, double n, double m, double p, 
										 double at_bt, double vomegapiover2, double x0);

inline double omega_pp(double vphi, double ctbt, double n, double m){
	return ctbt*((1-n)*elliptic_inc_pim(vphi,n,m)	- elliptic_inc_fm(vphi, m));
}

double omega_p(double Phi, double ctbt, double n, double m, double bt, double beta);

double omegaDerivative_p(double Phi, double ctbt, double n, double m, double bt);

double NR_diskPar(double tol, double ctbt, double n, double m, 
									double bt, double beta, double vomegabeta, double x0);

inline double I(double u, double w){
	return 0.5*(u*w*std::sqrt(1-u*u) + std::asin(u));
}

double NR_spherePar(int type, double tol, double vApay, double ax,
										double ay, double xe, double yl, double x0);


inline double Ap1(double ax,	double ay, double y){
	return 2*ax*ay*I(y/ay, 1);
}

inline double Ap1Derivative(double ax, double ay, double y){
	return 2*ax*std::sqrt(1-(y*y)/(ay*ay));
}

double Ap3(double ax,	double ay, double xe, 
									double yl, double y);

double Ap3Derivative(double ax, double ay, double xe, double yl, double y);

inline double Ap2(double ax,	double ay, double xe, double yl, double y){
	return Ap1(ax, ay, y) + Ap3(ax, ay, xe, yl, y);
}

inline double Ap2Derivative(double ax, double ay, 
										double xe, double yl, double y){
	return Ap1Derivative(ax, ay, y) + Ap3Derivative(ax,ay,xe,yl,y);
}

double Ar1(double ax,	double ay, double phi);

inline double rmax1(double ax, double cosbeta, double sinphi){
	return ax/sqrt(1-cosbeta*cosbeta*sinphi*sinphi);
}

inline double rmax2(double xe, double sinphi, double cosphi){
	return sqrt(1-xe*xe*sinphi*sinphi)-xe*cosphi;
}

double Ar3(double ax,	double ay, double xe, double phil, double phi);

double Ar3Derivative(double ax, double ay, double xe, 
										double phil, double phi, double beta);

inline double Ar2(double ax,	double ay, double xe, 
									double phil, double phi){
	return Ar1(ax, ay, phi) + Ar3(ax, ay, xe, phil, phi);
}

double Ar2Derivative(double ax, double ay, 
										double xe, double phil, double phi, double beta);

double NR_sphereRad(int type, double tol, double vApay, double ax,
										double ay, double xe, double yl, double beta, double x0);

#endif