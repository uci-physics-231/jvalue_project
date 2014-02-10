#ifndef FUN_SPH_STELLAR_H
#define FUN_SPH_STELLAR_H

#include <cmath>
#include <iostream>
#include "GregsMathHdr.h"
#include "AndrewsMathHdr.h"

using namespace std;

class StellarProfile_Sph
{
	public:
		virtual double intensity(double) = 0;// 2d intensity,  spherical symmetry
		virtual void InputParams(double *) = 0;//input parameter, 2-4 depending on the profile 
		virtual double rho_star(double) = 0;// 3d stellar density profile, spherical symmetry
		virtual double rho_star_dr(double)=0;// derivative of 3d stellar density profile wrt to r, spherical symmetry
		virtual double rScale(){return 0.0;}// return the scale radius, for plummer this is 2D rhalf, for king this is the core radius
		virtual double rTidal(){return 0.;}
		
		virtual double mass_total()=0;// integral rho_star dv,  ie the total mass of the system
		virtual int ParameterType(int ){return 0;}
		virtual int NumOfParams() = 0;//number of inputed parameters in the function InputParams(double *).  Needed for pointers when multiple classes are created
		virtual ~StellarProfile_Sph() {}//deconstrustors.  These should be called before when the class is deleted to stop seg faults
};

class Plummer_Sph : public StellarProfile_Sph
{
	//tested 7-29-2013
	//normalization of axisymmetric plummer profile is different
	private:
		double rhalf;// 2D hlaf light radius
	public:
		Plummer_Sph() {}
		~Plummer_Sph() {}
		void InputParams(double *in)
		{
			rhalf = in[0];
		}
		double rho_star(double rin)
		{
			//normalized such that mass  = 1
			double term = pow(1.0+pow(rin/rhalf,2.0), - 2.5) * .75 /PI/pow(rhalf,3.0);
			return term;
		}
		double rho_star_dr(double rin)
		{
			return -5.0*rin/SQR(rhalf) * pow(1.0+pow(rin/rhalf,2.0), -3.5)* .75 /PI/pow(rhalf,3.0);
		}
		double intensity(double r)
		{
			double bot = 1.0+SQR(r/rhalf);
			double tot = 1./SQR(bot) /PI/pow(rhalf,2.0);
			return tot;
		}
		double mass_total()
		{
			return 1.;
		}
		int ParameterType(int ){return 0;}
		double rScale() {return rhalf;}
		double rTidal() {return 0.;}
		int NumOfParams() {return 1;}
};

class Plummer_Sph_Cutoff : public StellarProfile_Sph
{
	//tested 7-29-2013
	//created 1-28-2014,  only added if statements
	//normalization of axisymmetric plummer profile is different
	private:
		double rhalf;// 2D hlaf light radius
		double rtidal;
	public:
		Plummer_Sph_Cutoff() {}
		~Plummer_Sph_Cutoff() {}
		void InputParams(double *in)
		{
			rhalf = in[0];
			rtidal = in[1];
		}
		double rho_star(double rin)
		{
			//normalized such that mass  = 1
			if(rin > rtidal)
				return 0.;
			double term = pow(1.0+pow(rin/rhalf,2.0), - 2.5) * .75 /PI/pow(rhalf,3.0);
			return term;
		}
		double rho_star_dr(double rin)
		{
			if(rin > rtidal)
				return 0.;
			return -5.0*rin/SQR(rhalf) * pow(1.0+pow(rin/rhalf,2.0), -3.5)* .75 /PI/pow(rhalf,3.0);
		}
		double intensity(double r)
		{
			if(r > rtidal)
				return 0.;
			double bot = 1.0+SQR(r/rhalf);
			double tot = 1./SQR(bot) /PI/pow(rhalf,2.0);
			return tot;
		}
		double mass_total()
		{
			return 1.;
		}
		int ParameterType(int ){return 0;}
		double rScale() {return rhalf;}
		double rTidal() {return rtidal;}
		int NumOfParams() {return 2;}
};

class King_Sph : public StellarProfile_Sph
{
	//tested 7-29-2013
	private:
		double rcore;// core radius
		double rtidal; // tidal radius
	public:
		King_Sph() {}
		~King_Sph() {}
		void InputParams(double *in)
		{
			rcore = in[0];
			rtidal = in[1];
		}
		double rho_star(double rin)
		{
			if(rin>rtidal)
				return 0.;
			double z = sqrt((1.+SQR(rin/rcore))/(1.+SQR(rtidal/rcore)));
			double con = PI*rcore*pow(1.+SQR(rtidal/rcore),1.5);
			double term = (acos(z)/z-sqrt(1.-z*z))/z/z;
			return term/con;
		}
		double rho_star_dr(double rin)
		{
			if(rin>rtidal)
				return 0.;
			double z = sqrt((1.+SQR(rin/rcore))/(1.+SQR(rtidal/rcore)));
// 			z = .1;
			double con = -rin/PI/pow(rcore,3.)/pow(1.+rtidal*rtidal/rcore/rcore, 2.5);
			double top = z*z*z - z +3.*sqrt(1.-z*z)*acos(z);
			double bot = pow(z,5.)*sqrt(1.-z*z);
// 			cout <<"king stuff  "<< con << "  " <<  top << "  " << bot << endl;
			return con * top/bot;
		}
		double intensity(double r)
		{
			if(r>rtidal)
				return 0.;
			double term = 1./sqrt(1.+SQR(r/rcore)) - 1./sqrt(1.+SQR(rtidal/rcore));
			return SQR(term);
		}
		double mass_total()
		{
			//added and tested 9-3-2013
			double term = -4. + SQR(rtidal)/(rcore*rcore + SQR(rtidal)) + 4.*rcore /sqrt(SQR(rcore) + SQR(rtidal)) + log(1. + SQR(rtidal/rcore)); 
			return PI * rcore*rcore *term;
		}
		
		double rScale() {return rcore;}
		double rTidal() {return rtidal;}
		int NumOfParams() {return 2;}
};

class Sersic_Sph : public StellarProfile_Sph, GaussLegendre
{
	//dont use this, need better integrator.
	// exp varys too much to use gauss legrende
	private:
		double re;// core radius
		double m; // tidal radius
		
		double integrand(double xx, double bigr)
		{
			double rr = xx /(1.- xx);
			double term = intensity_dr(rr) /sqrt(rr*rr - bigr*bigr);
			return term / SQR(1. - xx);
		}
		
		double integral (double rin)
		{
			double ll = rin/(1. + rin);
			double term = NIntegrate(static_cast <double (GaussianIntegral::*)(double, double)>(&Sersic_Sph::integrand), ll,  1., rin);
			
			return -term /PI ;
		}
		
		double intensity_dr(double rin)
		{
			double term = -exp(- pow(rin/re, 1./m)) * pow(rin/re, 1./m - 1.)/m/re;
			return term;
		}
		
	public:
		Sersic_Sph() : GaussLegendre(50) {}
		~Sersic_Sph() {}
		void InputParams(double *in)
		{
			re = in[0];
			m = in[1];
		}
		double rho_star(double rin)
		{
			double  p = 1. - 0.6097/m + 0.05463/m/m;
			double nuo = Gamma(2. * m) /2./re/Gamma((3.-p)*m);
			double term = pow(rin/re, -p) * exp(-pow(rin/re, 1./m));
			cout << term* nuo << "  " << integral(rin) << endl;
			return term* nuo;
		}
		double rho_star_dr(double rin)
		{
			double  p = 1. - 0.6097/m + 0.05463/m/m;
			double nuo = Gamma(2. * m) /2./re/Gamma((3.-p)*m);
			double term = pow(rin/re, -p) * exp(-pow(rin/re, 1./m));
			double drpart = (m * p + pow(rin/re, 1./m))/m/rin;
			return -term* nuo*drpart;
		}
		double intensity(double r)
		{
			double term =  exp(-pow(r/re,1./m) );
			return term;
		}
		double mass_total()
		{
			//added and tested 9-3-2013
			double term = 2. * PI* m * re*re * Gamma(2. * m);
			return term;
		}
		
		double rScale() {return re;}
		double rTidal() {return 0.;}
		int NumOfParams() {return 2;}
};

class Sersic_Sph_Approximation : public StellarProfile_Sph
{
	//tested 7-29-13
	//not sure if this is the best approximation
	private:
		double re;// core radius
		double m; // tidal radius
	public:
		Sersic_Sph_Approximation() {}
		~Sersic_Sph_Approximation() {}
		void InputParams(double *in)
		{
			re = in[0];
			m = in[1];
		}
		double rho_star(double rin)
		{
			//note that this is an approximation, should be ok ~5%
			double  p = 1. - 0.6097/m + 0.05463/m/m;
			double nuo = Gamma(2. * m) /2./re/Gamma((3.-p)*m);
			double term = pow(rin/re, -p) * exp(-pow(rin/re, 1./m));
			return term* nuo;
		}
		double rho_star_dr(double rin)
		{
			//note that this is an approximation, should be ok ~5%
			double  p = 1. - 0.6097/m + 0.05463/m/m;
			double nuo = Gamma(2. * m) /2./re/Gamma((3.-p)*m);
			double term = pow(rin/re, -p) * exp(-pow(rin/re, 1./m));
			double drpart = (m * p + pow(rin/re, 1./m))/m/rin;
			return -term* nuo*drpart;
		}
		double intensity(double r)
		{
			double term =  exp(-pow(r/re,1./m) );
			return term;
		}
		double mass_total()
		{
			cout << " dont use yet :( king profile spherical" <<endl;
			return 0.;
		}
		
		double rScale() {return re;}
		double rTidal() {return 0.;}
		int NumOfParams() {return 2;}
};

#endif
