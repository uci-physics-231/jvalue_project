#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>

#include "mcmchdr.h"
#include "GregsMathHdr.h"
#include "AndrewsMathHdr.h"

// #include "all_sph_multi.h"
// #include "null.h"
// #include "all_null_like.h"
// #include "all_multi_pop.h"
#include "all_sph_multi.h"

// #include "combining_Data.h"

#include <sstream>
#include <pthread.h>
#include <iomanip>
#include "/usr/local/boost_1_55_0/boost/program_options.hpp"

// pthread_mutex_t mutexs[4] = {PTHREAD_MUTEX_INITIALIZER};
namespace po = boost::program_options;
using namespace std;

int main(int ac, char **av)  
{
	int samplesize = 2000;//number of parameter points in the nested sampling, 2000 will work for the nested sampling
	int paranumber = 0;
	int numstars = -1;//number of stars in the data set
	string outputname;//name of file that data will go in
	string inputname;//file containing the data
	string parameter_name;// parameter file name
	string model_string;
	
	vector <string> models;
// 
	models.push_back("Spherical_Jeans");

	vector <string> stellar_models;	
	stellar_models.push_back("SPH_Plummer");
	stellar_models.push_back("SPH_King");	

	vector <string> dark_models;
	dark_models.push_back("SPH_NFW");
	dark_models.push_back("SPH_Burkert");
	
	vector <string> beta_models;
	beta_models.push_back("Isotropic");
	beta_models.push_back("Constant");
	beta_models.push_back("Functional");
	beta_models.push_back("Exp_Form");
	
// 	int *like_params_new;
	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("program", po::value<string>(), "set program to run")
		("data", po::value<string>(), "data file--required")
		("output",po::value<string>(),"name of output file--required")
		("parameter", po::value<string>(), "name of parameter file--required")
		("num_star", po::value<int>(), "number of stars--required")
		("num_parameters",po::value<int>(), "number of parameters for mcmc--required")
		("models","List Models")
		("stellar",po::value<int>()->default_value(0), "stellar distribution, -1 for options")
// 		("sph_stellar",po::value<int>()->default_value(0), "spherical stellar distribution, -1 for options")
		("dark", po::value<int>()->default_value(0), "dark matter distribution, -1 for options")
		("beta", po::value<int>()->default_value(0) , "beta function, -1 for options")
		("num_of_threads", po::value<int>()->default_value(1), "Number of threads of multi-threaded programs")
// 		("axi_jeans_model", po::value<int>()->default_value(0), "axi jeans options, beta_phi or beta_z model")
		("distance_prior", po::value<int>()->default_value(0), "jeans options, kpc versus degrees input")
// 		("axis_dm_prior", po::value<int>()->default_value(0), "axi jeans options, axi dm gaussian prior")
		("stellar_prior", po::value<int>()->default_value(0), "jeans options, stellar prior")
// 		("membership", po::value<bool>()->default_value(false), "Is there membership included in the file")
// 		("temp", po::vaule<int>()->default_value(0), "fuck")
	;

	po::variables_map vm;        
	po::store(po::parse_command_line(ac, av, desc), vm);
	po::notify(vm);    
	
	if (vm.count("help"))
	 {
		cout << desc << "\n";
		return 0;
	}
	
	if (vm.count("models"))
	{
		cout << endl;
		for(int i = 0;i<models.size();i++)
		{
			cout << models[i] << "\n";
		}	
		// cout << endl;
		return 0;
	}
	if (vm.count("program")) 
	{
		cout << "Program was set to " << vm["program"].as<string>() << ".\n";
		model_string = vm["program"].as<string>();
	}
	else
	{
		cout << "Program was not set.\n";
		return 0;
	}
	
	if(vm.count("data"))
	{
		inputname = vm["data"].as<string>();
		cout << "Data File Name " <<  inputname <<  ".\n";
	}
	else 
	{
		cout << "No Data file " << endl;
		return 0;
	}
	
	if(vm.count("output"))
	{
		outputname = vm["output"].as<string>();
		cout << "Output File Name " <<  outputname <<  ".\n";
	}
	else 
	{
		cout << "No output name selected --output" << endl;
		return 0;
	}
	
	if(vm.count("parameter"))
	{
		parameter_name = vm["parameter"].as<string>();
		cout << "Parameter File " << parameter_name << ".\n";
	}
	else 
	{
		cout << "No parameter file  --parameter" << endl;
		return 0;
	}
		
	if(vm.count("num_star"))
	{
		numstars = vm["num_star"].as<int>();
		cout << "Number of Stars " <<  numstars <<  ".\n";
	}
	else 
	{
		cout << "Set number of stars --num_star " << endl;
		return 0;
	}
	
	if(vm.count("num_parameters"))
	{
		paranumber = vm["num_parameters"].as<int>();
		cout << "Number of Stars " <<  paranumber <<  ".\n";
	}
	else 
	{
		cout << "Set number of parameters --num_parameters " << endl;
		return 0;
	}
	
	if(vm.count("stellar") && vm["stellar"].as<int>() == -1)
	{
		for (int i=0;i<stellar_models.size();i++)
			cout << i << "  =  " << stellar_models[i] << "\n";
		return 0;
	}
	if(vm.count("dark") && vm["dark"].as<int>() == -1)
	{
		for (int i=0;i<dark_models.size();i++)
			cout << i << "  =  " << dark_models[i] << "\n";
		return 0;
	}
	if(vm.count("beta") && vm["beta"].as<int>() == -1)
	{
		for (int i=0;i<beta_models.size();i++)
			cout << i << "  =  " << beta_models[i] << "\n";
		return 0;
	}
		
	char buffer[250];//needed for inputting data
	vector<string> vhold;
	
	ifstream ifs(parameter_name.c_str());
	
	ifs.getline(buffer, 250); vhold = parse(buffer, 250);// 2D half-light radius
	double rhalf2d = static_cast <double> ( atof( vhold[0].c_str() ) );
	double rhalf2d_error = static_cast <double> ( atof( vhold[1].c_str() ) );
	ifs.getline(buffer, 250); vhold = parse(buffer, 250);// tidal radius
	double rtidal = static_cast <double> ( atof( vhold[0].c_str() ) );
	double rtidal_error = static_cast <double> ( atof( vhold[1].c_str() ) );
	ifs.getline(buffer, 250); vhold = parse(buffer, 250);//core radius
	double rcore = static_cast <double> ( atof( vhold[0].c_str() ) );
	double rcore_error = static_cast <double> ( atof( vhold[1].c_str() ) );
	ifs.getline(buffer, 250); vhold = parse(buffer, 250);// axis ratio
	double obs_axis_ratio = static_cast <double> ( atof( vhold[0].c_str() ) );
	double obs_axis_ratio_error = static_cast <double> ( atof( vhold[1].c_str() ) );
	ifs.getline(buffer, 250); vhold = parse(buffer, 250);// skiping 3d rhalf
	ifs.getline(buffer, 250); vhold = parse(buffer, 250);// skiping postion angle
	double obs_position_angle = static_cast <double> ( atof( vhold[0].c_str() ) );
	double obs_position_angle_error = static_cast <double> ( atof( vhold[1].c_str() ) );
// 	postion_angle = (90. - postion_angle);
	ifs.getline(buffer, 250); vhold = parse(buffer, 250);//limits for velocity
	double avr_vel_lower = static_cast <double> ( atof( vhold[0].c_str() ) );
	double avr_vel_upper = static_cast <double> ( atof( vhold[1].c_str() ) );
	ifs.getline(buffer, 250); vhold = parse(buffer, 250);
	double obs_distance = static_cast <double> ( atof( vhold[0].c_str() ) );
	double obs_distance_error = static_cast <double> ( atof( vhold[1].c_str() ) );
	ifs.getline(buffer, 250); vhold = parse(buffer, 250);
	double obs_ra_center = static_cast <double> ( atof( vhold[0].c_str() ) );
	double obs_ra_center_error = static_cast <double> ( atof( vhold[1].c_str() ) );
	ifs.getline(buffer, 250); vhold = parse(buffer, 250);
	double obs_dec_center = static_cast <double> ( atof( vhold[0].c_str() ) );
	double obs_dec_center_error = static_cast <double> ( atof( vhold[1].c_str() ) );
	ifs.getline(buffer, 250); vhold = parse(buffer, 250);
	double m_v = static_cast <double> ( atof( vhold[0].c_str() ) );
	double m_v_error = static_cast <double> ( atof( vhold[1].c_str() ) );
	ifs.getline(buffer, 250); vhold = parse(buffer, 250);
	double m_l_stellar = static_cast <double> ( atof( vhold[0].c_str() ) );
	ifs.close();  
	
	double galaxy_params[] = {
		obs_distance,
		obs_distance_error,
		rhalf2d,
		rhalf2d_error,
		rcore,
		rcore_error,
		rtidal,
		rtidal_error,
		obs_position_angle,
		obs_position_angle_error,
		obs_axis_ratio,
		obs_axis_ratio_error,
		obs_ra_center,
		obs_ra_center_error,
		obs_dec_center,
		obs_dec_center_error,
		m_v,
		m_v_error,
		m_l_stellar
	};
	
	cout << model_string <<"--model--"<< endl;
	cout << "DID YOU CHANGE THE OUTPUT NAME ?? ?  " << outputname << endl;
	cout << endl;
	
	getchar();
	
	UCMC * mcmc = NULL; 

	if(model_string.compare("Spherical_Jeans") == 0)
	{
		int stellar_profile_int = vm["stellar"].as<int>();
		int dm_profile_int = vm["dark"].as<int>();
		int beta_int = vm["beta"].as<int>();
		
		int distance_int = vm["distance_prior"].as<int>(); //choose degree (ra, dec) vs kpc
// 		int jeans_prior = vm["axi_jeans_model"].as<int>();
		//jeans model is the next entry
// 		double scale_max_temp= rtidal;
// 		if(likelihood_int == 1)
// 			scale_max_temp *= PI/60./180.*obs_distance;
			
		vector <double>  llimit, ulimit;
		llimit.push_back(avr_vel_lower);
		ulimit.push_back(avr_vel_upper);
		if(distance_int == 0)
		{
			llimit.push_back(obs_distance - obs_distance_error * 5.); 
			ulimit.push_back(obs_distance + obs_distance_error * 5.);
		}
		
		cout << stellar_profile_int << "---Stellar Density Profile---";
		if(stellar_profile_int == 0)
		{
			cout << "Plummer Profile" << endl;
			double min = rhalf2d - rhalf2d_error*5.;
			if(min < 0.)
				min = 0.;
			ulimit.push_back(rhalf2d + rhalf2d_error*5.);
			llimit.push_back(min);  
// 			input_place+=1;
		}
		else if (stellar_profile_int == 1)
		{
			cout << "King Profile" << endl;
			double min = rcore - rcore_error*5.;
			double min_tidal = rtidal - rtidal_error*5.;
			if(min < 0.)
				min = 0.;
			if(min_tidal < 0.)
				min_tidal = 0.;
			
			ulimit.push_back(rcore + rcore_error*5.);// rcore
			llimit.push_back(min);  
			ulimit.push_back(rtidal + rtidal_error*5.);// rtidal
			llimit.push_back(min_tidal);  
// 			input_place+=2;
		}
		else
		{
			cout << "WRONG" << endl;
			return -1;
		}
		
		cout << dm_profile_int << "---Dark Matter Profile---";
		if(dm_profile_int == 0)//nfw
		{
			cout << "Spherical NFW" << endl;
			ulimit.push_back(3.);// log10(10. * scale_max_temp * 2.));//rs
			llimit.push_back(-3.);  
			ulimit.push_back(10.);// rhos
			llimit.push_back(0.);  
// 			input_place+=2;
		}
		else if(dm_profile_int == 1)
		{
			cout << "Spherical Burkert" << endl;
			ulimit.push_back(3.); //log10(10. * scale_max_temp * 2.));//rs
			llimit.push_back(-3.);  
			ulimit.push_back(10.);// rhos
			llimit.push_back(0.);  
// 			input_place+=2;
		}
		// 
// 		else if(dm_profile_int == 5 || dm_profile_int == 6)
// 		{
// 			cout << "dm + stellar" << endl;
// 			ulimit.push_back( log10(10. * scale_max_temp));//rs
// 			llimit.push_back(-3.);  
// 			ulimit.push_back(11.);// rhos
// 			llimit.push_back(0.);  
// 			ulimit.push_back(m_v + m_v_error * 5.);
// 			llimit.push_back(m_v - m_v_error * 5.); 
// 			if(ml_prior == 1)
// 			{
// 				ulimit.push_back(2.);// M/L
// 				llimit.push_back(0.1);  
// 			}
// 			input_place+=4;
// 		}
		else
		{
			cout << "WRONG" << endl;
			return -1;
		}

		cout << beta_int << "---Beta Function---";
		if(beta_int ==0)
		{
			cout << "Isotropic, beta = 0" << endl;
		}
		else if (beta_int == 1)
		{
			cout << "beta = constant" << endl;
			ulimit.push_back(1.);// vo
			llimit.push_back(-10.); 
// 			input_place+=1;
		}
		else if (beta_int == 2)
		{
			cout << "functional form" << endl;
			ulimit.push_back(1.);// beta 0
			llimit.push_back(-10.); 
			ulimit.push_back(1.); // beta 1
			llimit.push_back(-10.);
			ulimit.push_back( 3.);// change in beta function slope
			llimit.push_back(-3.);
// 			input_place +=3;
		}
		else if ( beta_int == 3)
		{
			cout << "exp form" << endl;
			ulimit.push_back(1.);// beta 0
			llimit.push_back(-10.); 
			ulimit.push_back(1.); // beta 1
			llimit.push_back(-10.);
			ulimit.push_back(3. );// change in beta function slope
			llimit.push_back(-3.);
// 			input_place +=3;
		}
		else
		{
			cout << "NONE" << endl;
			return -1;
		}
		
		double upperlimit[paranumber];
		double lowerlimit[paranumber];
		for(int i =0;i<paranumber;i++)
		{
			upperlimit[i] = ulimit[i];
			lowerlimit[i] = llimit[i];
			cout << ulimit[i] << "  " << llimit[i] << endl;
		}
		int temp_hold =vm["distance_prior"].as<int>();
		int threads_hold = vm["num_of_threads"].as<int>();
		int fuck[]= {stellar_profile_int, dm_profile_int, beta_int, temp_hold};
		
		int cols[] = {7, 2,3,4,5, 6};

		mcmc = new MCMC_Sph_Model(numstars, galaxy_params, fuck, inputname, cols, threads_hold, lowerlimit, upperlimit, lowerlimit, paranumber);
	}
	mcmc->MonoSample(outputname.c_str(), samplesize);
	cout << mcmc->out_totalpoints() << "  " << mcmc->outlnZ() << "  " << mcmc->outmp() << endl;
 

	cout << outputname << endl;
		
	return 0;
}


