#ifndef ALL_SPH_MULTI_H
#define ALL_SPH_MULTI_H

#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>

#include "mcmchdr.h"
#include "GregsMathHdr.h"
#include "AndrewsMathHdr.h"
// #include "sph_jeans.h"
#include "sph_jeans_new.h"
#include "InputStuff.h"

#include "fun_stellar.h"
#include "fun_potential.h"
#include "fun_beta.h"

#include <sstream>
#include <pthread.h>
#include <iomanip>

// pthread_mutex_t mutexs[4] = {PTHREAD_MUTEX_INITIALIZER};

using namespace std;

class Sph_Like_Data : public UCMC, protected obs_parameters
{
protected:
	int numOfData;
	int offset;
	int skip;
	double *ra;
	double *dec;
	double *vel_los;
	double *err_los;
	double *member;
	GeneralLineOfSight_Sph *solver;
	
	int stellar_density_int;//chooses stellar density profile
	int dark_matter_profile_int;// chooses dark matter density profile
	int betafun_int;//chooses beta function
	
	int place_stellar_prior;
// 	int place_mv_v_prior;
	
// 	int vel_perspec_int;
// 	int gauss_order;//order of gauss legendre integration
// 	int likelihood_int;
// 	int lumi_int;
	int distance_int;//kpc or degree input
	
	int params_num;
	
	double fill(double * in, double * output)
	{
		int *array_prior = new int[params_num];
		solver->Fill(array_prior);
		
		int input_place =1; 
		double distance = in[0];
		double conver = distance * PI/ 60./180.;
		if(distance_int == 1)
		{
			conver = 1.;
			input_place = 0;
		}
		double *temp = in + input_place;
		int place =0;
// 		cout << "fill  " << params_num - input_place << endl;
		for(int i =0;i < params_num - input_place;i++)
		{
// 			cout << "I" << in[i] <<"  ";
// 			cout << i << "  "  << array_prior[i] << endl;
			if(array_prior[i] == 0)//normal flat prior
			{
				output[i] = temp[place];
				place++;
			}
			else if(array_prior[i] == 1)
			{
				output[i] = pow(10.,  temp[place]);// log10 prior
				place++;
			}
// 			else if(array_prior[i] == 2)// inclination angle
// 				output[i] = i_angle;
// 			else if(array_prior[i] == 3)// axis ratio
// 				output[i] = intrinsic_stellar_axis;
			else if(array_prior[i] == 4)// convert from arcmin to kpc
			{
				output[i] = temp[place] * conver;
				place++;
			}
// 			cout <<i << "  "  << array_prior[i] << "  " << output[i] << "  " <<endl;
			
		}
// 		getchar();
// 		cout << endl;
// 		cout << "array shit inputend fucj " << endl;
	}

public:
	Sph_Like_Data(const int nin, double * obs_params, int * choice_input) : numOfData(nin),  obs_parameters(obs_params)
	{
		cout << "start  " << numOfData << endl;
		offset = 0;
		skip = 1;

		stellar_density_int = choice_input[0];
		dark_matter_profile_int = choice_input[1];
		betafun_int = choice_input[2];
// 		vel_perspec_int = choice_input[3];
// 		gauss_order = choice_input[4];
// 		likelihood_int = choice_input[5];
// 		lumi_int = choice_input[6];
		distance_int = choice_input[3];
		
// 		cout << "fuck " << endl;
		StellarProfile_Sph * stellar;
		DarkMatterProfile_Sph * dark;
		BetaFunction * betafunction;
// 		cout << "fuck " << endl;
		int input_place = 2;
		if (distance_int ==1 )
			input_place =1;
		place_stellar_prior = input_place;
// 		cout << "fuck " << endl;
		cout << stellar_density_int << "--Stellar Profile--";
		if(stellar_density_int == 0)//interloper_int =0, no interlopers
		{
			cout << "Plummer--No prior||";
			input_place+=1;
// 			place_mv_v_prior = place_stellar_prior + 1;
			stellar = new Plummer_Sph;
		}
		else if(stellar_density_int == 1)//interloper_int =0, no interlopers
		{
			cout << "King--No prior||";
// 			place_stellar_prior = input_place;
// 			place_mv_v_prior = place_stellar_prior + 2;
			input_place+=2;
			stellar = new King_Sph;
		}
		else
		{
			cout << "NONE selected" << endl;
			return;
		}
		
		cout << dark_matter_profile_int << "--Dark Matter Profile--";
		if (dark_matter_profile_int == 0)
		{
			cout << "NFW||";
			input_place +=2;
			dark = new NFW_Sph;
		}
		else if (dark_matter_profile_int == 1)
		{
			cout << "Burkert||";
			input_place +=2;
			dark = new Burkert_Sph;
		}
		else if (dark_matter_profile_int == 2)
		{
			cout << "Einasto||";
			input_place +=3;
			dark = new Einasto_Sph;
		}
		else if (dark_matter_profile_int == 3)
		{
			cout << "NFW and Plummer star potential||";
			input_place+=2;
			dark = new Multi_Potential_Sph(new NFW_Sph, new Plummer_Potential_Sph);
		}
		else
		{
			cout << "NONE selected" << endl;
			return;
		}
		
		cout << betafun_int << "--Beta Function--";
		if (betafun_int == 0)
		{
			cout << "Isotropic||";
			betafunction = new IsotropicBeta;
// 			input_place +=1;
		}
		else if (betafun_int == 1)
		{
			cout << "Constant||";
			betafunction = new ConstantBeta;
			input_place +=1;
		}
		else if(betafun_int == 1)
		{
			cout << "Beta Rad function||";
			betafunction = new MultiBeta;
			input_place += 3;
		}
		else if(betafun_int == 2)
		{
			cout << "Beta Rad function Exp||";
			betafunction = new Beta_Exp;
			input_place += 3;
		}
		else
		{
			cout << "NONE selected" << endl;
			return;
		}
/*		
		cout << lumi_int << "--Lumi stuff--";
		if (lumi_int == 0)
			cout << "no lumi";
		else if(lumi_int == 1)
		{
			cout << "m_v parameter||";
			place_mv_v_prior = input_place;
			input_place += 1;
		}
		else if(lumi_int == 2)
		{
			cout << "m_v parameter & m/l ||";
			place_mv_v_prior = input_place;
			input_place += 2;
		}
		else
		{
			cout << "NONE selected" << endl;
			return;
		}*/
		
		
		cout << endl;
// 		solver = new GeneralLineOfSight_Sph(gauss_order, new Sph_Jeans_Beta(gauss_order,  stellar,  dark, betafunction));
		solver = new Sph_Jeans(stellar,  dark, betafunction);
		params_num = 0;
		for(int i=0;i<5;i++)
		{
// 			if(solver->NumOfParams(i) >= 0)
				params_num += solver->NumOfParams(i);
// 				cout << i << "  " << params_num << "  " << solver->NumOfParams(i) << endl;
		}
// 		cout << "  " << params_num << endl;
	}
	~Sph_Like_Data() 
	{
		solver->~GeneralLineOfSight_Sph();
		delete solver;
	}
	void InputData(double * xin, double * yin, double *vel_in, double * err_in, double * mem_in)
	{
		ra = xin;
		dec = yin;
		vel_los = vel_in;
		err_los = err_in;
		member = mem_in;
	}
	
	double LogLike_standard_sph(double *in)
	{
		int num_params = 0;
// 			cout << "num_params  " << num_params << endl;
		for(int i=0;i<4;i++)
		{
			num_params += solver->NumOfParams(i);
// 				cout << i << "  " << num_params << "  " << solver->NumOfParams(i) << endl;
		}
// 			cout << "num_params  " << num_params << endl;
		double los_params_to_input[num_params];
// 			double params[] = {-290.7, .223, log10(1.5), 7. , -.1};
// 			in = params;
// 			double los_params_to_input[] = ;
		double avr_vel = in[0];//systemic velocity
		double distance = in[1];
		double ra_cen = obs_ra_center, dec_cen = obs_dec_center;
		if (distance_int == 1)
		{
			distance = 180./PI;
			ra_cen = 0.;
			dec_cen = 0.;
		}
// 		double m_v = 
		double top=0.0;
// 		double bot = 0.0;
		fill(in+1, los_params_to_input);
		solver->InputParams(los_params_to_input);

		for (int i = offset; i < numOfData; i+=skip)
		{
			double x = out_x(ra[i], dec[i] , distance, obs_position_angle, ra_cen, dec_cen);
			double y = out_y(ra[i], dec[i] , distance, obs_position_angle, ra_cen, dec_cen);
// 			cout << x << "  " << y << "  " << ra[i] - obs_ra_center << endl;
			double rad = sqrt(x*x + y*y);
			
			double hold2 = solver->get_velocity_los_square(rad) ;
			if(hold2 < 0.0 )
			{
// 					cout << in[0] << "  " << in[1] << "  " <<in[2] << "  " << in[3] << "  " << in[4] << "  " << endl;
				return 1e20;
			}
			double sigma = hold2  + SQR(err_los[i]);
			
// 			bot += log(2.0 * PI * sigma);

			top += member[i] * (pow(vel_los[i] - avr_vel ,2.0)/sigma + log(2.0 * PI * sigma) ); 
// 				if(top*0. || bot*0.)
// 				{
// 					cout << i << "  " << skip << endl;
// 					break;
// 				}
// 				os << i << "  " << (pow(vel_los[i] - avr_vel ,2.0)/sigma + log(2.0 * PI * sigma))/2.0 << "  " << hold2 << "  " << sqrt(hold2) << endl;
// 			if()
		}
// 			os.close();
		double temp =  (top )/2.0;
// 		cout << numOfData << " " << temp << endl;
		if(temp*0.)
			return 1e20;
// 			cout << offset << "  " << temp << "  " << endl;
		return temp;
	}

	double LogPrior_to_thread(double *in)
	{
		double term =0.;
		for(int i= 0;i<ma;i++)
		{
// 				if( i != 2 )
				term += log(fabs(upperLimits[i]-lowerLimits[i]));
		}
		if(stellar_density_int == 0 )
		{
			double rh = in[place_stellar_prior];
// 				cout << rh << "  " <<  rhalf_stellar << "  " <<  rhalf_stellar_error << "  " << place_stellar_prior << endl;
			term -=log(fabs(upperLimits[place_stellar_prior]-lowerLimits[place_stellar_prior]));
			term += LogGaussianPrior(rh, obs_rhalf_2d, obs_rhalf_2d_error, upperLimits[place_stellar_prior], lowerLimits[place_stellar_prior]);
		}
		else if(stellar_density_int == 1)
		{
			double rc = in[place_stellar_prior];
			double rt = in[place_stellar_prior+1];
			
			term -= log(fabs(upperLimits[place_stellar_prior]-lowerLimits[place_stellar_prior]));
			term -= log(fabs(upperLimits[place_stellar_prior+1]-lowerLimits[place_stellar_prior+1]));
			term += LogGaussianPrior(rc, obs_rcore, obs_rcore_error, upperLimits[place_stellar_prior], lowerLimits[place_stellar_prior]);
			term += LogGaussianPrior(rt, obs_rtidal, obs_rtidal_error, upperLimits[place_stellar_prior+1], lowerLimits[place_stellar_prior+1]);
		}
		// if(likelihood_int == 1)
// 		{
// 			double dist = in[1];
// 			double ra_cen = in[2];
// 			double dec_cen = in[3];
// 			for(int i = 0;i<3;i++)
// 				term -= log(fabs(upperLimits[1 + i]-lowerLimits[1+i]));
// 			term += LogGaussianPrior(dist, obs_distance, obs_distance_error, upperLimits[1], lowerLimits[1]);
// 			term += LogGaussianPrior(ra_cen, obs_ra_center, obs_ra_center_error, upperLimits[2], lowerLimits[2]);
// 			term += LogGaussianPrior(dec_cen, obs_dec_center, obs_dec_center_error, upperLimits[3], lowerLimits[3]);
// 		}
// 		else if(likelihood_int == 2)
// 		{
// 			double dist = in[1];
// 			double ra_cen = in[2];
// 			double dec_cen = in[3];
// // 			for(int i = 0;i<3;i++)
// 				term -= log(fabs(upperLimits[1]-lowerLimits[1]));
// 			term += LogGaussianPrior(dist, obs_distance, obs_distance_error, upperLimits[1], lowerLimits[1]);
// // 			term += LogGaussianPrior(ra_cen, obs_ra_center, obs_ra_center_error, upperLimits[2], lowerLimits[2]);
// // 			term += LogGaussianPrior(dec_cen, obs_dec_center, obs_dec_center_error, upperLimits[3], lowerLimits[3]);
// 		}
		
// 		if (lumi_int == 1 || lumi_int == 2)
// 		{
// 			double m_v_para = in[place_mv_v_prior];
// 			term -= log(fabs(upperLimits[place_mv_v_prior]-lowerLimits[place_mv_v_prior]));
// 			term += LogGaussianPrior(m_v_para, obs_m_v, obs_m_v_error, upperLimits[place_mv_v_prior], lowerLimits[place_mv_v_prior]);
// 		}
		
		return term;
	}
	void SetSkip(const int oin, const int sin)
	{
		offset = oin; skip = sin;
	}
};

struct Thing_Sph_Data
{
	//FullLikeBin *louies;
	Sph_Like_Data *louies;
	double *point;
	double like;
// 	double prior;
	
	void Input_Class(int star, double *in_galaxy_params, int * in_choice_params )
	{
	      louies = new Sph_Like_Data(star, in_galaxy_params, in_choice_params );
	}
};

void *f(void *ptr)
{
// 	std::cout << "test" << std::endl;
	Thing_Sph_Data *in = (Thing_Sph_Data *)ptr;
// 	std::cout << "test" << std::endl;LogLike_to_thread
	in->like = in->louies->LogLike_standard_sph(in->point);
	return NULL;
}

class MCMC_Sph_Model : public UCMC, protected inputcrap
{
	private:
 		Thing_Sph_Data *thing;
 		int threadsNum;
		int num;
		
	public:
		MCMC_Sph_Model(int star, double *in_galaxy_params, int * in_choice_params, string name, int * cols, int threads_num_in, double *startpoint, double *upperlimit, double *lowerlimit, int paranumber ) : num(star), threadsNum(threads_num_in), inputcrap(name, cols[0])
		{
			//#include "chi2.h"
//  			threadsNum = threads_num_in;  
			InputPoint(startpoint, upperlimit, lowerlimit, paranumber);
 			thing = new Thing_Sph_Data[threadsNum];
			cout << "first line of data  " << datashit[cols[1]][0] << "  " <<  datashit[cols[2]][0] << "  " <<  datashit[cols[3]][0]<< "  " << datashit[cols[4]][0]  << "threads " << threadsNum << "  " << num << endl;
 			for (int i = 0; i < threadsNum; i++)
 			{
				//Chi2 *temp;// = new Chi2(num);
// 				thing[i].louies = new Chi2(num);
// 				std::cout << "constructor !!!" << std::endl;
				thing[i].Input_Class(num, in_galaxy_params, in_choice_params);
// 				std::cout << "constructor !!" << std::endl;
// 				thing[i].louies->inputdata(mcmcx,mcmcy, mcmcvel, mcmcerr, parameterinputarray);
				thing[i].louies->SetSkip(i, threadsNum);
				if(cols[6] > 0)
				{
					thing[i].louies->InputData(datashit[cols[1]], datashit[cols[2]], datashit[cols[3]], datashit[cols[4]], datashit[cols[5]]);
// 					cout << datashit[cols[2]][0] << "  " <<  datashit[cols[3]][0] << "  " <<  datashit[cols[4]][0]<< "  " << datashit[cols[5]][0] << "  " << datashit[cols[6]][0] << endl;
				}
				else 
				{
					double * mem;
					mem = matrix(star, 1.);
					thing[i].louies->InputData(datashit[cols[1]], datashit[cols[2]], datashit[cols[3]], datashit[cols[4]], mem);
// 					cout << datashit[cols[2]][0] << "  " <<  datashit[cols[3]][0] << "  " <<  datashit[cols[4]][0]<< "  " << datashit[cols[5]][0] << "  " << mem[0] << endl;
				}
				thing[i].louies->InputPoint(startpoint, upperlimit, lowerlimit, paranumber);
// 				std::cout << "constructor ! endl" << std::endl;
// 				thing[i].louies->InputParams(params_in);
 			}
// 			cout << "constructor end" << endl;
		}
		double LogLike(double *a)
		{
// 			std::cout << "test  0" << std::endl;
 			pthread_t *threads = new pthread_t[threadsNum];
// 			std::cout << "test  1" << std::endl;
 			for (int i = 0; i < threadsNum; i++)
 			{
 				thing[i].point = a;
 				pthread_create(&threads[i], NULL, f, (void *)(thing+i));
				
 			}
//  			std::cout << "test  2" << std::endl;
 			for (int i = 0; i < threadsNum; i++)
 			{
 				pthread_join(threads[i], NULL);
 			}
 			double temp = 0;
// 			std::cout << "test  3" << std::endl;
 			for (int i = 0; i < threadsNum; i++)
 				temp += thing[i].like;
// 			std::cout << "test  4" << std::endl;
 			delete[] threads;
// 			cout << "fucking thing  " << temp << endl;
 			return temp;
		}

		double LogPrior(double *a)
		{
// 			cout << "log prior started " << endl;
			double temp = thing[0].louies->LogPrior_to_thread(a);
// 			cout << "log prior worked ?? " << endl;
 			return temp;
		}

		~MCMC_Sph_Model()
		{
// 			cout << " seg " << endl;
// 			for (int i = 0; i < threadsNum; i++)
// 			{
// 				cout << "pre  seg  " << i << endl;
// // 				thing[i].louies->~Likelihood_Prespective_Besancon();
// // 				delete thing[i].louies;
// 				cout << "post  seg  " << i << endl;
// 			}
// 			cout << " seg " << endl;
 			delete[] thing;
// 			cout << " seg " << endl;
		}
};

#endif

