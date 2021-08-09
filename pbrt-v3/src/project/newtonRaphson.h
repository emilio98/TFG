#ifndef NEWTONRAPHSON_H
#define NEWTONRAPHSON_H

#include "project/ellipticIntegral.h"
#include<cmath>

inline double omega_r2(double Phi, double n, double m, double p, double at_bt){
	return Phi-p*elliptic_inc_pim(std::atan(at_bt*std::tan(Phi)),n,m);
}

double omegaDerivative_r2(double Phi, double n, double m, double p, double at_bt);

double NR_diskRad(double tol, double n, double m, double p, 
										 double at_bt, double vomegapiover2);

inline double I(double u, double w){
	return 0.5*(u*w*std::sqrt(1-u*u) + std::asin(u));
}

double NR_spherePar(int type, double tol, double vApay, double ax,
										double ay, double xe, double yl);


inline double Ap1(double ax,	double ay, double y){
	return 2*ax*ay*I(y/ay, 1);
}

inline double Ap1Derivative(double ax, double ay, double y){
	return 2*ax*std::sqrt(1-(y*y)/(ay*ay));
}

inline double Ap3(double ax,	double ay, double xe, 
									double yl, double y){
	if(y > yl)
		return I(yl,1) - xe*yl - 0.5*Ap1(ax, ay, yl);
	else
		return I(y,1) - xe*y - 0.5*Ap1(ax, ay, y);
}

inline double Ap3Derivative(double ax, double ay, double xe, 
										double yl, double y){
	if(y > yl)
		return 0;
	else
		return std::sqrt(1-y*y) - xe - 0.5*Ap1Derivative(ax, ay, y);
}

inline double Ap2(double ax,	double ay, double xe, 
									double yl, double y){
	return Ap1(ax, ay, y) + Ap3(ax, ay, xe, yl, y);
}

inline double Ap2Derivative(double ax, double ay, 
										double xe, double yl, double y){
	return Ap1Derivative(ax, ay, y) + Ap3Derivative(ax,ay,xe,yl,y);
}

#endif