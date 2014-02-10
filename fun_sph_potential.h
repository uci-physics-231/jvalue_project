#ifndef FUN_SPH_POTENTIAL_H
#define FUN_SPH_POTENTIAL_H

#include <cmath>
#include <iostream>
#include "GregsMathHdr.h"
#include "AndrewsMathHdr.h"

using namespace std;

class DarkMatterProfile_Sph
{
public:
	virtual void InputParams(double *)=0;
	virtual double phi(double)=0;
	virtual double phi_dr(double) =0;
	virtual double mass(double) =0;
	virtual double rScale(){return 0.0;}
	virtual double RhoS(){return 0.0;}
	virtual int NumOfParams()=0;
	virtual int ParameterType(int ){return 1;}
	virtual ~DarkMatterProfile_Sph() {}
};

class NFW_Sph : public DarkMatterProfile_Sph //spherical nfw potential 
{
//tested 7-27-2013
private:
	double rhoS;
	double rS;
// 	double rMax;
// 	double vMax;
	double Gaa;
	
	inline double dphi_dre(const double ree)
	{
		//double ree = sqrt(R*R+z*z);
		double constants = Gaa * pow(rS, 3.0) * rhoS * 4.0 * PI;
		if(ree < .001 )
		{
			if(ree== 0.0)
				return constants*.5;
			double g = ree / rS;
			double dphidlittler =  constants *(.5 - 2.0 *g/3.0 + .75 *g*g - .8* pow(g,3.0) + 5.0*pow(g,4.0)/6.0)/rS/rS;
			return dphidlittler;
		}
		double dphi_drec = pow(ree, -2.0) * (log(1.0 + ree/rS) - ree/(ree+rS));
		if(dphi_drec*0.0) 
		{
			cout << "NFW_spherical::dphi_dre\t" << ree << "  " << rS << " " << rhoS << endl; 
			return 0.0;
		}
		return constants * dphi_drec;
	}
	
public:
	NFW_Sph() : Gaa(.0000043) {}//Gaa is gravitional constant in useful astro units.  
					// units of [G] = kpc* (solar Mass)^-1 * (km/s)^2
	~NFW_Sph() {}
	void InputParams(double *in)
	{
// 			double xMax = 2.16258;
// 			vMax = pow(10.0, in[0]);
// 			rMax = pow(10.0, in[1]);
// 			rS = rMax/xMax;
// 			rhoS = vMax*vMax*mpcToKm*1.0e12/(4.0*PI*G*rMax*rMax*SM/xMax/(1.0+xMax)/(1.0+xMax));
		rS  = in[0];
		rhoS = in[1];
// 			cout << "Inputparams dark matter " << rS << "\t" << rhoS << endl; getchar();
	}

	double phi(double rin)
	{
		double con =  -4.0*PI*Gaa * rhoS*SQR(rS);
		double rr = rin/rS;
		return con * log(1.0+rr)/rr;
	}
	double phi_dr(double rin)
	{
		double con =  4.0*PI*Gaa * rhoS * rS;
		double xx = rin/rS;
		return con * (log(1.0+xx)/SQR(xx) - 1./xx/(1.+xx));
	}
	double mass (double rin)
	{
		////tested 10/16/2013
		double con = 4. *PI * rhoS * pow(rS, 3.);
		double x = rin/rS;
		double term = log(1. +x ) - x /(1.+x);
		return term * con;
	}
	double OutputParams(int i)
	{
		if(i==0)
			return rS;
		else if(i ==1)
			return rhoS;
		else 
			return 0.;
	}
	double rScale(){return rS;}
	double RhoS(){return rhoS;}

	int NumOfParams(){return 2;}
};

class Burkert_Sph : public DarkMatterProfile_Sph //spherical nfw potential 
{
	//tested 7-27-2013
private:
	double rhos;
	double rs;
	double gravity;
	double con;
	
public:
	Burkert_Sph() : gravity(.0000043) {}//Gaa is gravitional constant in useful astro units.  
					// units of [G] = kpc* (solar Mass)^-1 * (km/s)^2
	~Burkert_Sph() {}
	void InputParams(double * in)
	{
		rs = in[0];
		rhos = in[1];
		con = PI *gravity * rhos* pow(rs, 2.0);
	}
	double phi(double r)
	{
		double x = r/rs;
		double trig_term = 2.*x *acot(x) - 2. * atan(x);
		double log_term = -2. *x *log(x/(1.+x)) - x*log(1.+1./x/x)+2.*log(1.+x)+log(1.+x*x);
		return -con *(trig_term+log_term)/x;
	}
	double phi_dr(double r)
	{
		double x = r/rs;
		double term = -2.* atan(x) +2.*log(1.+x) + log(1.+x*x);
		return term * con/rs /x/x;
	}
	double mass(double r)
	{
		cout << "not made" << endl;
		return 0.;
	}
	double OutputParams(int i)
	{
		if(i==0)
			return rs;
		else if(i ==1)
			return rhos;
		else 
			return 0.;
	}
	double rScale(){return rs;}
	double RhoS(){return rhos;}

	int NumOfParams(){return 2;}
};

class Einasto_Sph : public DarkMatterProfile_Sph //spherical nfw potential 
{
//tested 7-27-2013
private:
	double rhos;
	double rs;
	double Gaa;
	double alpha;
	
public:
	Einasto_Sph() : Gaa(.0000043) {}//Gaa is gravitional constant in useful astro units.  
					// units of [G] = kpc* (solar Mass)^-1 * (km/s)^2
	~Einasto_Sph() {}
	void InputParams(double *in)
	{
		rs  = in[0];
		rhos = in[1];
		alpha = in[2];
	}

	double phi(double rin)
	{
		double con =  -4.0*PI*Gaa * rhos*SQR(rs);
		double xx = rin/rs;
		double term1 = exp(2./alpha)/alpha * pow(2./alpha, -2./alpha) * IncGammaUp(2./alpha, 2./alpha * pow(xx,alpha));
		double term2 = pow(2., -3./alpha) * pow(alpha, 3./alpha -1.)*exp(2./alpha) * IncGamma(3./alpha, 2./alpha*pow(xx,alpha))/xx;
// 			cout << term1 *4.* PI<< "  " << term2 *4.* PI << endl;
		return con * (term1 + term2);
	}
	double phi_dr(double rin)
	{
		double x = rin/rs;
		double con = 4. * PI * exp(2./alpha)/alpha * Gaa * rhos * rs;
		double term1 = 2. * pow(2./alpha, -1.-2./alpha) * x * IncGammaUp(1.+2./alpha, 2./alpha*pow(x,alpha));
		double term2 = -2. *x * pow( 2./alpha, -2./alpha) * IncGammaUp(2./alpha,2./alpha*pow(x,alpha) );
		double term3 = -alpha * pow(x,3.)*exp(-2./alpha * pow(x,alpha));
		double term4 = pow(8.,-1./alpha)*pow(alpha, 3./alpha) * IncGamma(3./alpha, 2./alpha*pow(x,alpha));
// 			cout << term1 << "  " << term2 << "  " << term3 << "  " << term4 << endl;
		return con * (term1+term2+term3 + term4)/x/x;
	}
	double mass(double r)
	{
		cout << "not made" << endl;
		return 0.;
	}
	double OutputParams(int i)
	{
		if(i==0)
			return rs;
		else if(i ==1)
			return rhos;
		else if(i ==2)
			return alpha;
		else 
			return 0.;
	}
	double rScale(){return rs;}
	double RhoS(){return rhos;}

	int NumOfParams(){return 3;}
};

class Plummer_Potential_Sph : public DarkMatterProfile_Sph //spherical plummer potential 
{
//tested  10-16-2013
private:
	double mass_para;
	double rhalf;
	double Gaa;
	
public:
	Plummer_Potential_Sph() : Gaa(.0000043) {}//Gaa is gravitional constant in useful astro units.  
					// units of [G] = kpc* (solar Mass)^-1 * (km/s)^2
	~Plummer_Potential_Sph() {}
	void InputParams(double *in)
	{
		rhalf  = in[0];
		mass_para = in[1];
	}

	double phi(double rin)
	{
		double term =  -mass_para * Gaa/rhalf  /  sqrt(1. + SQR(rin/rhalf));
// 			double rr = rin/rS;
		return term;
	}
	double phi_dr(double rin)
	{
		double con =  Gaa * mass_para/ pow(rhalf, 2.);
		double xx = rin/rhalf;
		return con * xx / pow( 1. + xx*xx, 1.5);
	}
	double mass(double r)
	{//tested 10/16/2013
		double term = mass_para * pow(r, 3.)/pow(r*r + rhalf*rhalf, 1.5);
		return term;
	}
	double OutputParams(int i)
	{
		if(i==0)
			return rhalf;
		else if(i ==1)
			return mass_para;
		else 
			return 0.;
	}
	double rScale(){return rhalf;}
	double RhoS(){return mass_para/pow(rhalf, 3.);}

	int NumOfParams(){return 2;}
};

class Multi_Potential_Sph : public DarkMatterProfile_Sph //spherical nfw potential 
{
//tested  10-16-2013, with plummer and nfw
private:
	DarkMatterProfile_Sph *dark1, *dark2;
	
public:
	Multi_Potential_Sph(DarkMatterProfile_Sph *in1, DarkMatterProfile_Sph *in2) 
	{
		dark1 = in1;
		dark2 = in2;
	}
	~Multi_Potential_Sph() 
	{
		dark1->~DarkMatterProfile_Sph();
		dark2->~DarkMatterProfile_Sph();
		delete dark1;
		delete dark2;
	}
	void InputParams(double *in)
	{
		dark1->InputParams(in);
		dark2->InputParams(in + dark1->NumOfParams());
	}

	double phi(double rin)
	{
		double term = dark1->phi(rin) + dark2->phi(rin);
		return term;
	}
	double phi_dr(double rin)
	{
		double term = dark1->phi_dr(rin) + dark2->phi_dr(rin);
		return term;
	}
	double mass(double r)
	{
		//tested 10/16/2013
		double term = dark1->mass(r) + dark2->mass(r);
		return term;
	}
// 		double OutputParams(int i)
// 		{
// 			if(i<dark1->NumOfParams())
// 				return dark1->OutputParams(i);
// 			else if(i< (dark1->NumOfParams()+dark2->NumOfParams()))
// 				return dark2->OutputParams(i - dark1->NumOfParams());
// 			else 
// 				return 0.;
// 		}
	double rScale(){return dark1->rScale();}
	double RhoS(){return dark1->RhoS();}
	int NumOfParams()
	{
		return dark1->NumOfParams() + dark2->NumOfParams();
	}
};

#endif
