#ifndef FUN_BETA_H
#define FUN_BETA_H

#include <cmath>
#include <iostream>
#include "GregsMathHdr.h"
#include "AndrewsMathHdr.h"

using namespace std;


class BetaFunction
{
	public:
		virtual double betafun(double, double )=0;
		virtual double betafun_dR(double, double)=0;
		virtual double betafun(double)=0;
		virtual void InputParams(double *)=0;
		virtual int NumOfParams()=0;
		virtual double beta_integrate_factor(double)=0;// integrating factor 
		virtual double exp_integrate(double)=0;// exp of integrating factor
		virtual double exp_integrate(double, double){return 0.;}
		virtual int ParameterType(int ){return 0;}
		virtual ~BetaFunction(){}
};

class ConstantBeta : public BetaFunction
{
	// class tested 7-24-2013
	private:
		double beta; 
		
	public:
		ConstantBeta(){}
		~ConstantBeta(){}
		double betafun_dR(double rin, double zin)
		{
			return 0.0;    
		}
		double betafun(double rin, double zin)
		{
			return beta;    
		}
		double betafun(double ree)
		{
			return beta;
		}
		double exp_integrate(double rin)
		{
			return pow(rin, 2.*beta);
		}
		double exp_integrate(double rin, double zin)
		{//only for axisymmetry and beta phi model
			return pow(rin, beta);
		}
		double beta_integrate_factor(double rin)
		{
			double term = 2. * beta * log(rin);
			return term;
		}
		
		void InputParams(double *ain)
		{
			beta = ain[0]; 
		}
		int ParameterType(int i)
		{
				return 0;
		}
		int NumOfParams(){return 1;}
};

class IsotropicBeta : public BetaFunction
{
	// class tested 7-24-2013
	private:
// 		double beta; 
		
	public:
		double betafun_dR(double rin, double zin)
		{
			return 0.0;    
		}
		double betafun(double rin, double zin)
		{
			return 0.0;    
		}
		double betafun(double rin)
		{
			return 0.0;
		}
		void InputParams(double *ain)
		{
			//beta = ain[0]; 
		}
		double beta_integrate_factor(double rin)
		{
			return 0.;
		}
		double exp_integrate(double rin)
		{
			return 0.;
		}
		int NumOfParams(){return 0;}
};

class MultiBeta_Beta_Phi : public BetaFunction
{
	//tested 1-13-2014
	private:
		double betazero; 
		double betaone;
		double betarad;
// 		double axis;
	public:
		MultiBeta_Beta_Phi() {}
		~MultiBeta_Beta_Phi() {}
		double betafun_dR(double rin, double zin)
		{			
			double bot = SQR(betarad) + SQR(rin);
			
			return 2.* (betaone- betazero)* SQR(betarad/bot) * rin;
		}
		double betafun(double rin, double zin)
		{
// 			double ree = sqrt(SQR(rin) + SQR(zin/axis));
			double term;
			term = (betaone - betazero)*rin*rin/(rin*rin + betarad*betarad) + betazero;

			return term;    
		}
		double betafun(double ree)
		{
			return 0.;    
		}//this class isnt made to use with spherical symmetry
		double beta_integrate_factor(double rin)
		{
			return 0.;
		}
		double exp_integrate(double rin)
		{
			return 0.0;
		}
		double exp_integrate(double rin, double zin)
		{
			return pow(rin, betazero)*pow(rin*rin + betarad*betarad, (betaone- betazero)/2.);
		}
		void InputParams(double *ain)
		{
			betazero = ain[0]; 
			betaone = ain[1];
			betarad = ain[2];
		}
		int ParameterType(int i)
		{
			if(i ==0 || i == 1)
				return 0;
			else if (i ==2 )
				return 1;
			else
				return 0;
		}
		int NumOfParams(){return 3;}
		double Betazero(){return betazero;}
		double Betaone(){return betaone;}
		double Betarad(){return betarad;}
};

class MultiBeta : public BetaFunction
{
	// class tested 7-24-2013
	private:
		double betazero; 
		double betaone;
		double betarad;
		double axis;
	public:
		MultiBeta() {}
		~MultiBeta() {}
		double betafun_dR(double rin, double zin)
		{
			if(betaone==0.0 && betarad==0.0) 
				return 0.0;
			double ree = sqrt(SQR(rin) + SQR(zin/axis));
			
			double bot = SQR(betarad) + SQR(ree);
			
			return 2.* (betaone- betazero)* SQR(betarad/bot) * rin;
		}
		double betafun(double rin, double zin)
		{
			double ree = sqrt(SQR(rin) + SQR(zin/axis));
			double term;
			term = (betaone - betazero)*ree*ree/(ree*ree + betarad*betarad) + betazero;

			return term;    
		}
		double betafun(double ree)
		{
			double term = (betaone - betazero)*ree*ree/(ree*ree + betarad*betarad) + betazero;

			return term;    
		}
		double beta_integrate_factor(double rin)
		{
			return 2. * betazero *log(rin) + (betaone-betazero)*log(SQR(rin)+SQR(betarad));
		}
		double exp_integrate(double rin)
		{
			return pow(rin, 2.* betazero)*pow(rin*rin + betarad*betarad, betaone- betazero);
		}
		void InputParams(double *ain)
		{
			betazero = ain[0]; 
			betaone = ain[1];
			betarad = ain[2];
			axis = ain[3];
		}
		int NumOfParams(){return 4;}
		double Betazero(){return betazero;}
		double Betaone(){return betaone;}
		double Betarad(){return betarad;}
};

class Beta_Exp : public BetaFunction
{
	// class tested 7-24-2013
	private:
		double betazero; 
		double betaone;
		double betarad;
		double axis;
	public:
		Beta_Exp() {}
		~Beta_Exp() {}
		double betafun_dR(double rin, double zin)
		{
			if(betaone==0.0 && betarad==0.0) 
				return 0.0;
			double ree = sqrt(SQR(rin) + SQR(zin/axis));
			
			double term = (betazero - betaone)* rin *exp(2.- ree/betarad)*(ree - 2.*betarad)/4./pow(betarad,3.);
			
			return term;
		}
		double betafun(double rin, double zin)
		{
			double ree = sqrt(SQR(rin) + SQR(zin/axis));
			double term;
			term = (betaone - betazero)*SQR(ree/2./betarad)*exp(2. - ree/betarad) + betazero;

			return term;    
		}
		double betafun(double ree)
		{
			double term = (betaone - betazero)*SQR(ree/2./betarad)*exp(2. - ree/betarad) + betazero;

			return term;    
		}
		double beta_integrate_factor(double rin)
		{
			double term = ((betazero - betaone)*(rin+betarad) * exp(2.)+ 2.*betazero * exp(rin/betarad)*betarad*log(rin))*exp(-rin/betarad)/betarad;
			return term;
		}
		double exp_integrate(double rin)
		{
			double term = pow(rin, 2.* betazero)* exp((betazero - betaone)*(rin+ betarad)/betarad*exp(2.- rin/betarad));
			return term;
		}
		void InputParams(double *ain)
		{
			betazero = ain[0]; 
			betaone = ain[1];
			betarad = ain[2];
			axis = ain[3];
		}
		int NumOfParams(){return 4;}
		double Betazero(){return betazero;}
		double Betaone(){return betaone;}
		double Betarad(){return betarad;}
};

#endif
