#ifndef LOS_SPH_JEANS_H
#define LOS_SPH_JEANS_H

#include <iostream>
#include <cmath>
#include "GregsMathHdr.h"
#include "AndrewsMathHdr.h"
#include "fun_sph_potential.h"
#include "fun_sph_stellar.h"
#include "fun_beta.h"

using namespace std;

class GeneralLineOfSight_Sph :  public  GaussLegendre
{
//this is stil 5 % off or so
protected:
// 		JeansEqnSph *jeans;
	StellarProfile_Sph *stars;
	double rtidal;
	double los_sph(double rin)
	{
	// 	global_r = rin;
	// 	double rtidal = jeans->rTidal();
		double term;
		if(rtidal == 0.)
		{
			term = NIntegrate(static_cast <double (GaussianIntegral::*)(double, double)>(&GeneralLineOfSight_Sph::los_integrand), rin/(1.+rin), 1.0, rin);
		}
		else 
		{
			if(rin >= rtidal)
				return 0.;
	// integration limits are R -> rt , not R-> sqrt(rt^2 - R^2)  (dz limits are 0 -> sqrt(rt^2 - R^2) )
			term = NIntegrate(static_cast <double (GaussianIntegral::*)(double, double)>(&GeneralLineOfSight_Sph::los_integrand), rin/(1.+rin), rtidal/(1.+ rtidal), rin);
		}

		double in= stars->intensity(rin);
		double temp = 2. * term / in; 
	// 	cout << rin << "  " << temp << "  s" << skip  <<"s  " << term << "  " <<in << '\n';
		return temp;
	}
	double los_integrand(double rin, double bigr)
	{
	// 	double bigr = global_r;   
		double fun_in = rin/(1.-rin);

		double term1 = (1. -beta_profile(fun_in) *SQR(bigr/fun_in) );
		double term2 = vrr(fun_in) * fun_in/sqrt(SQR(fun_in) - SQR(bigr));
		return term1*term2 /SQR(1.-rin);
	}
public:
	GeneralLineOfSight_Sph() : GaussLegendre(10) 
	{
// 			jeans = in_sph;
	}

	~GeneralLineOfSight_Sph()
	{
		stars->~StellarProfile_Sph();
		delete stars;
	}
	void input_stars(StellarProfile_Sph *stars_in)
	{
		stars = stars_in;
	}
	void input_star_params(double *in)
	{
		stars->InputParams(in);
		rtidal = stars->rTidal();
	}

	virtual double vrr(double)=0;
	virtual double beta_profile(double)=0;
	virtual void InputParams(double *) =0;
	virtual int NumOfParams(int i) =0;
	
	virtual void Fill(int *){};//for mcmc 
	
	double get_velocity_los_square(double rin)
	{
// 			double rin = sqrt(SQR(xin) + SQR(yin));
		return los_sph(rin);
	}
};

#endif
