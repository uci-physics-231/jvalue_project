#ifndef SPH_JEANS_NEW_H
#define SPH_JEANS_NEW_H

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include "GregsMathHdr.h"
#include "AndrewsMathHdr.h"
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include "fun_sph_potential.h"
#include "fun_sph_stellar.h"
#include "fun_beta.h"

#include "los_sph_jeans.h"

using namespace std;

class  Sph_Jeans : public GeneralLineOfSight_Sph
{
//tested with plummer and nfw 1-27-2014
//limiting radius untested
protected:
// 	double beta;
	StellarProfile_Sph *stars;
	DarkMatterProfile_Sph *dark;
	BetaFunction *beta;
	
	double integrad(double rin)
	{
		double fun_in = rin/(1.-rin);
		double temp =beta->exp_integrate(fun_in) * stars->rho_star(fun_in) * dark->phi_dr(fun_in)/pow(1. - rin,2.0);
		if(temp*0. ||  temp<0.)
			cout << rin << "  " << fun_in << "  " << temp << "  " << beta->exp_integrate(fun_in) << "  " << stars->rho_star(fun_in) << "  " << dark->phi_dr(fun_in) << endl;
		return temp;
	}
	double vel_r_square(double rin)
	{
		double temp;
		
		if(stars->rTidal() == 0.)
			temp =  NIntegrate(static_cast <double (GaussianIntegral::*)(double)>(&Sph_Jeans::integrad), rin/(1.+rin), 1.0);
		else 
		{
			if(rin > stars->rTidal())
				return 0.;
			temp =  NIntegrate(static_cast <double (GaussianIntegral::*)(double)>(&Sph_Jeans::integrad), rin/(1.+rin), rtidal/(1.+ rtidal));		
		} 
		if(temp*0.0 || temp <0.0) 
		{
			cout << "o noes problems---------->Sph_Jeans_ConstBeta::vel_r_square  " << temp << "  " << rin << "  " << beta->exp_integrate(rin) << endl; 
	// 		return 0.0;
		} 
		double term = beta->exp_integrate(rin);
		return temp/term;
	}

public:
	Sph_Jeans (StellarProfile_Sph *in_star, DarkMatterProfile_Sph *in_dark, BetaFunction * in_beta)
	{
		stars = in_star;
		input_stars(stars);
		dark = in_dark;
		beta = in_beta;
	}
	~Sph_Jeans ()
	{
		stars->~StellarProfile_Sph();
		delete stars;
		dark->~DarkMatterProfile_Sph();
		delete dark;
		beta->~BetaFunction();
		delete beta;
	}
	void InputParams(double *in)
	{
		stars->InputParams(in);
		input_star_params(in);
		dark->InputParams(in+stars->NumOfParams());
// 			cout << in[0] << "  " << in[1] << "  " << in[2] << "  " << in[3] << "  " << in[4] << endl;
// 		double *temp = in + stars->NumOfParams() + dark->NumOfParams();
		beta->InputParams(in+stars->NumOfParams() + dark->NumOfParams());
	}
	void Fill(int * out)
	{
		for (int i=0;i<stars->NumOfParams();i++)
			out[i] = stars->ParameterType(i);
		for (int i=0;i< dark->NumOfParams();i++)
			out[stars->NumOfParams() + i] = dark->ParameterType(i);
		for (int i=0;i< beta->NumOfParams();i++)
			out[stars->NumOfParams()+ dark->NumOfParams() + i] = beta->ParameterType(i);
	}
	double vrr(double r)
	{
		return vel_r_square(r);		
	}
	double beta_profile(double rin)
	{
		return beta->betafun(rin);
	}

	int NumOfParams(int i)
	{
		if(i==0)
			return stars->NumOfParams();
		if(i==1)
			return dark->NumOfParams();
		if(i==2)
			return beta->NumOfParams();
		else 
			return 0;
	}
};

#endif