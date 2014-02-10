#include <iostream>
#include <cmath>

#include "/usr/local/boost_1_55_0/boost/program_options.hpp"

namespace po = boost::program_options;
using namespace std;

int main(int ac, char **av)  
{
// 	int samplesize = 2000;//number of parameter points in the nested sampling, 2000 will work for the nested sampling
// 	int paranumber = 0;
// 	int numstars = -1;//number of stars in the data set
// 	string outputname;//name of file that data will go in
// 	string inputname;//file containing the data
// 	string parameter_name;// parameter file name
// 	string model_string;

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
		("axi_jeans_model", po::value<int>()->default_value(0), "axi jeans options, beta_phi or beta_z model")
		("distance_prior", po::value<int>()->default_value(0), "axi jeans options, kpc versus degrees input")
		("axis_dm_prior", po::value<int>()->default_value(0), "axi jeans options, axi dm gaussian prior")
		("stellar_prior", po::value<int>()->default_value(0), "axi jeans options, stellar prior")
		("symmetry", po::value<int>()->default_value(0), "stellar model option, sph/axi symmetry 0/1")
		("rotation", po::value<int>()->default_value(0), "stellar model option, x,y y,z  x,z")
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
	return 0;
}
