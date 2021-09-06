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
										 double at_bt, double vomegapiover2, double x0){
	double xi=x0, f = 20, fd;
	int i = 0, max_it = 10;

	while(abs(f)>tol && i<max_it){
		f=omega_r2(xi,n,m,p,at_bt)-vomegapiover2;
		fd=omegaDerivative_r2(xi,n,m,p,at_bt);
		xi=xi-f/fd;
		i++;
	}	

	if(abs(f)<=tol){
		return xi;
	}

	xi=0;
	i = 0;
	while(abs(f)>tol && i<max_it){
		f=omega_r2(xi,n,m,p,at_bt)-vomegapiover2;
		fd=omegaDerivative_r2(xi,n,m,p,at_bt);
		xi=xi-f/fd;
		i++;
	}

	if(abs(f)<=tol){
		return xi;
	}

	xi=M_PI_2;
	i = 0;
	while(abs(f)>tol && i<max_it){
		f=omega_r2(xi,n,m,p,at_bt)-vomegapiover2;
		fd=omegaDerivative_r2(xi,n,m,p,at_bt);
		xi=xi-f/fd;
		i++;
	}

	if(abs(f)<=tol){
		return xi;
	}

	/*if(abs(f)>tol){
		cout << xi << " xi\n";
		cout << f << " f\n";
		cout << fd << " fd\n";
		cout << vomegapiover2 << " vomegapiover2\n";
	}*/
	
	return x0;
}

double omega_p(double Phi, double ctbt, double n, double m, double bt, double beta){
	if(Phi<0){
		return omega_pp(M_PI_2, ctbt, n, m) - omega_pp(std::asin(std::tan(-Phi)/bt), ctbt, n, m);
	}else{
		return omega_pp(M_PI_2, ctbt, n, m) + omega_pp(std::asin(std::tan(Phi)/bt), ctbt, n, m);
	}
}

double omegaDerivative_p(double Phi, double ctbt, double n, double m, 
												double bt){
	double vphi = std::asin(std::tan(Phi)/bt);
	double sinvphi = std::sin(vphi);
	double sinvphi2 = sinvphi*sinvphi;

	double cosPhi = std::cos(Phi);
	double cosPhi2 = cosPhi*cosPhi;
	double tanPhi = std::tan(Phi);
	double tanPhi2 = tanPhi * tanPhi;

	double f = 1/(std::sqrt(1-m*sinvphi2));
	double pi = 1/((1-n*sinvphi2)*sqrt(1-sinvphi2*m));
	double vphiDerivative = 1/(bt*cosPhi2*std::sqrt(1-tanPhi2/(bt*bt)));

	return ctbt*((1-n)*pi*vphiDerivative - f*vphiDerivative);
}

double NR_diskPar(double tol, double ctbt, double n, double m, 
									double bt, double beta, double vomegabeta, double x0){
	double xi=x0, f = 20, fd;
	int i = 0, max_it = 10;

	while(abs(f)>tol && i<max_it){
		f=omega_p(xi, ctbt, n, m, bt, beta)-vomegabeta;
		fd=omegaDerivative_p(xi, ctbt, n, m, bt);
		xi=xi-f/fd;
		i++;
	}	

	if(abs(f)<=tol){
		return xi;
	}

	xi=beta*0.5;
	i = 0;
	while(abs(f)>tol && i<max_it){
		f=omega_p(xi, ctbt, n, m, bt, beta)-vomegabeta;
		fd=omegaDerivative_p(xi, ctbt, n, m, bt);
		xi=xi-f/fd;
		i++;
	}

	if(abs(f)<=tol){		
		return xi;
	}

	xi=-beta*0.5;
	i = 0;
	while(abs(f)>tol && i<max_it){
		f=omega_p(xi, ctbt, n, m, bt, beta)-vomegabeta;
		fd=omegaDerivative_p(xi, ctbt, n, m, bt);
		xi=xi-f/fd;
		i++;
	}

	if(abs(f)<=tol){		
		return xi;
	}

	/*if(abs(f)>tol){
		cout << xi << " xi\n";
		cout << f << " f\n";
		cout << fd << " fd\n";
		cout << vomegabeta << " vomegabeta\n";
	}*/
	
	return x0;
}

double NR_spherePar(int type, double tol, double vApay, double ax,
										double ay, double xe, double yl, double x0){
	double xi=x0, f = 20, fd;
	int i = 0, max_it = 10;

	if(type == 1){
		while(abs(f) > tol){
			f =Ap1(ax,ay,xi) - vApay;
			fd=Ap1Derivative(ax, ay, xi);
			xi=xi-f/fd;
			i++;
			//Si no converge
			if(isnan(xi) || isinf(xi) || abs(xi)>ay || i>max_it ){
				xi = x0;
				f=0.0;
			}
		}
	}else if (type == 2){
		while(abs(f) > tol){
			f =Ap2(ax,ay,xe,yl,xi) - vApay;
			fd=Ap2Derivative(ax, ay, xe, yl, xi);
			xi=xi-f/fd;
			i++;
			//Si no converge
			if(isnan(xi) || isinf(xi) || abs(xi)>ay || i>max_it ){
				xi = x0;
				f=0;
			}
		}
	}else{
		while(abs(f) > tol){
			f =Ap3(ax,ay,xe,yl,xi) - vApay;
			fd=Ap3Derivative(ax, ay, xe, yl, xi);
			xi=xi-f/fd;
			i++;
			//Si no converge
			if(isnan(xi) || isinf(xi) || abs(xi)>ay || i>max_it ){
				xi = x0;
				f=0;
			}
		}
	}
	
	return xi;
}

double Ap3(double ax,	double ay, double xe, 
									double yl, double y){
	if(y > yl)
		return I(yl,1) - xe*yl - 0.5*Ap1(ax, ay, yl);
	else
		return I(y,1) - xe*y - 0.5*Ap1(ax, ay, y);
}

double Ap3Derivative(double ax, double ay, double xe, 
										double yl, double y){
	if(y > yl)
		return 0;
	else
		return std::sqrt(1-y*y) - xe - 0.5*Ap1Derivative(ax, ay, y);
}

double Ar1(double ax,	double ay, double phi){
	if(phi<=M_PI_2)
		return 0.5*ax*ay*atan((ax/ay)*tan(phi));
	else
		return M_PI + 0.5*ax*ay*atan((ax/ay)*tan(phi));
}

double Ar2Derivative(double ax, double ay, 
										double xe, double phil, double phi, double beta){
	double r;
	if(phi > phil)
		r = rmax1(ax,cos(beta), sin(phi));
	else
		r = rmax2(xe,sin(phi),cos(phi));
	
	return 0.5*r*r;
}

double Ar3(double ax,	double ay, double xe, 
									double phil, double phi){
	double minphi = min(phil,phi);
	double sinphim = sin(minphi);
	return I(sinphim,xe*xe)-I(xe*sinphim,1)-Ar1(ax,ay,minphi);
}

double Ar3Derivative(double ax, double ay, double xe, 
										double phil, double phi, double beta){
	if(phi > phil)
		return 0;
	else{
		double r1 = rmax1(ax,cos(beta), sin(phi)), r2 = rmax2(xe,sin(phi),cos(phi));
		return 0.5*(r2*r2 -r1*r1);
	}
}

double NR_sphereRad(int type, double tol, double vApay, double ax,
										double ay, double xe, double phil, double beta, double x0){
	double xi=x0, f = 20, fd, max_it=15;
	int i = 0;

	if (type == 2){
		while(abs(f) > tol){
			f =Ar2(ax,ay,xe,phil,xi) - vApay;
			fd=Ar2Derivative(ax, ay, xe, phil, xi, beta);
			xi=xi-f/fd;
			i++;
			//Si no converge
			if(isnan(xi) || isinf(xi) || i>max_it ){
				xi = x0;
				f=0;
				//cout << "hola\n";
			}
		}
	}else{
		while(abs(f) > tol){
			f =Ar3(ax,ay,xe,phil,xi) - vApay;
			fd=Ar3Derivative(ax, ay, xe, phil, xi, beta);
			xi=xi-f/fd;
			i++;
			//Si no converge
			if(isnan(xi) || isinf(xi) || i>max_it ){
				xi = x0;
				f=0;
				//cout << "hola\n";

			}
		}
	}
	
	return xi;

}