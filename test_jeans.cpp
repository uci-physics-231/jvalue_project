#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include "GregsMathHdr.h"
#include "AndrewsMathHdr.h"
// #include "los_axi.h"
// #include "los_sph.h"

#include "sph_jeans_new.h"
#include "axi_jeans_beta_phi.h"

// #include "los_sph.h"

#include <sstream>
#include <pthread.h>
#include <iomanip>

using namespace std;

int main(int argc, char **argv)
{
	if(false)// los axi jeans, isotropic rotator stuff
	{
// 		-178.659                       
		
// 		double conver =  914.36 *PI/180./60.;
		double array[] = {.5, 0.75, PI/4., 1., 1e7,  -.5, .5, .1};
		Axi_Jeans_Beta_Phi * los = new Axi_Jeans_Beta_Phi(new Plummer_axi_oblate, new NFW_spherical_axi, new ConstantBeta);
		los->InputParams(array);
		
		double max = 1.;
		double step = .01;
		ofstream os;
// 		ofstream os("test_jeans_axi");
// 		for(double r = 0.; r<max; r+= step)
// 		{
// 			for(double z = 0.; z<max; z+=step)
// 			{
// 				os << r << "  " << z << "  " << los->vzz_int(r, z) << "  " << los->vrr_int(r, z)<< "  " << los->vphiphi_int(r	, z) << endl;
// 			}
// 		}
		double r =.1, z = .2;
		cout << los->vzz_int(r, z) << "  " << los->vrr_int(r, z)<< "  " << los->vphiphi_int(r, z) << endl;
		cout << "los  1 " << los->get_los_square(.1,.1) << "--los 2--" << los->get_los_square(.5,.5) << "--los3--" << los->get_los_square(.5,.25) << endl;
// 		os.close();	
		os.open("los_test2");
		for(double r = -max; r<max; r+= step)
		{
			for(double z = -max; z<max; z+=step)
			{
				os << "  " << r << "  " << z << "  " << los->get_los_square(r, z ) << endl;	
			}
		}
		
	}
	
	if(false)// los axi jeans, isotropic rotator stuff
	{
// 		-178.659                       
// 		57.2958            -6.26129
// 		double conver =  914.36 *PI/180./60.;
// 		double array[] = {13.7854/60., 0.65, 1.25, pow(10., -0.901346), pow(10., 9.4),  -0.1, .5, .5};
// 		double array[] = {13.1061/60., 0.616264, 1.23513, pow(10., -0.634576), pow(10., 8.89355),  -0.40899, .5, .5};
// 		double array[] = {11.7022/60., 0.61391, 1.21842, pow(10., -0.65523), pow(10., 8.88194),  -1.89881, -0.54987, pow(10., -2.1)};
		double array[] = {0.22552, 0.614218, 1.22719, pow(10., -.6), pow(10., 6.45693), -.5,  -.5 ,pow(10.,  -1.)};
		Axi_Jeans_Beta_Z * los = new Axi_Jeans_Beta_Z(new Plummer_axi_oblate, new NFW_spherical_axi, new ConstantBeta);//MultiBeta_Beta_Phi);
		los->InputParams(array);
		
		double max = 1.;
		double step = .005;
		ofstream os;
// 		ofstream os("test_jeans_axi");
// 		for(double r = 0.; r<max; r+= step)
// 		{
// 			for(double z = 0.; z<max; z+=step)
// 			{
// 				os << r << "  " << z << "  " << los->vzz_int(r, z) << "  " << los->vrr_int(r, z)<< "  " << los->vphiphi_int(r	, z) << endl;
// 			}
// 		}
		double r =.1, z = .2;
		cout << los->vzz_int(r, z) << "  " << los->vrr_int(r, z)<< "  " << los->vphiphi_int(r, z) << endl;
		cout << "los  1 " << los->get_los_square(.1,.1) << "--los 2--" << los->get_los_square(.5,.5) << "--los3--" << los->get_los_square(.5,.25) << endl;
// 		os.close();	
		os.open("los_const_betaz_betasph_neg");
		for(double r = step ; r<max; r+= step)
		{
			for(double z = -max; z<max; z+=step)
			{ 
				double hold = 0.;//los->get_los_square(r, z );
				double vrr = los->vrr_int(r, z);
				double vzz = los->vzz_int(r, z);
				double vphi = los->vphiphi_int(r, z);
				double top = r*r * vrr + z*z * vzz + (r*r + z*z) *vphi;
				double bot = z*z * vrr + r*r * vzz;
// 				cout << r << "  " << z << "  " << vrr << "  " << vzz <<  "  " << 1. - top/bot/2. << endl;
				os  << r << "  " << z << "  " << vrr << "  " << vzz << "  " << vphi  << "  " << 1. - top/bot/2. << endl;	
			}
		}
		
	}
	
	if(false)// test beta z version
	{
// 		-178.659                       
// 		{rh, qnu, qphi, ii, rs, ps, G, beta} = { .5, .75, 1, \[Pi]/4, 1, 1*^7, 4.3*^-6, .5}
// 		double conver =  914.36 *PI/180./60.;
		double array[] = {.5, 0.75, PI/4., 1., 1e7,  .5, .5, .1};
		Axi_Jeans_Beta_Z * los = new Axi_Jeans_Beta_Z(new Plummer_axi_oblate, new NFW_spherical_axi, new ConstantBeta);
		los->InputParams(array);
		
		double max = 1.;
		double step = .01;
		ofstream os;
// 		ofstream os("test_jeans_axi");
// 		for(double r = 0.; r<max; r+= step)
// 		{
// 			for(double z = 0.; z<max; z+=step)
// 			{
// 				os << r << "  " << z << "  " << los->vzz_int(r, z) << "  " << los->vrr_int(r, z)<< "  " << los->vphiphi_int(r	, z) << endl;
// 			}
// 		}
		double r =.4, z = .1;
		cout << los->vzz_int(r, z) << "  " << los->vrr_int(r, z)<< "  " << los->vphiphi_int(r, z) << endl;
// 		cout << "los  1 " << los->get_los_square(.1,.1) << "--los 2--" << los->get_los_square(.5,.5) << "--los3--" << los->get_los_square(.5,.25) << endl;

		
	}
	
	if(true)// los sph jeans testing 
	{
		double array[] = {.5,  1., pow(10., 7.), .25};
		double arrayking[] = {.25, 2.,  1., pow(10., 7.), .25};
		Sph_Jeans * los = new Sph_Jeans(new Plummer_Sph, new NFW_Sph, new ConstantBeta);//MultiBeta_Beta_Phi);
		Sph_Jeans * king = new Sph_Jeans(new King_Sph, new NFW_Sph, new ConstantBeta);
		los->InputParams(array);
		king->InputParams(arrayking);
		
		StellarProfile_Sph * plum = new Plummer_Sph;
		double params [] = {.5, 1e6};
		plum->InputParams(params);
		
		double max = 1.;
		double step = .005;
		ofstream os;
// 		ofstream os("test_jeans_axi");
// 		for(double r = 0.; r<max; r+= step)
// 		{
// 			for(double z = 0.; z<max; z+=step)
// 			{
// 				os << r << "  " << z << "  " << los->vzz_int(r, z) << "  " << los->vrr_int(r, z)<< "  " << los->vphiphi_int(r	, z) << endl;
// 			}
// 		}
		double r =.25, z = .9;
		
		cout << "vrr, second velocity moment--------" << endl;
		cout << "plum  "<< los->vrr(r) << "  " << los->vrr(z)   << endl;
		cout << "mathematica  " << 68.1363 << "  "  << 1.73858 << endl;
		cout << "king  " <<king->vrr(r) << "  " << king->vrr(z)   << endl;
		cout << "mathematica king  " << 25.3907 << "  " << 0.450449 << endl;
		cout << "line-of-sight velocity dispersion-------------" << endl;
		cout <<"plum  " << los->get_velocity_los_square(.1) << "  " << los->get_velocity_los_square(.5) << "  " << los->get_velocity_los_square(1.5) << endl;
		cout << "mathematica  " << 47.3594 << "  " << 25.902 << "  " << 21.8394 << endl;
		cout <<"king  " << king->get_velocity_los_square(.1) << "  " << king->get_velocity_los_square(.5) << "  " <<king->get_velocity_los_square(1.5) << endl;
		cout << "mathematica  " <<35.791 << "  " << 23.4678 << "  " << 7.64636 << endl;
		
		cout << "stellar profile -----------" << endl;
		cout << "c++ "  << plum->rho_star(.1) << "  " << plum->rho_star(.5) << endl;
		cout << "mathematica  " << 1.73148 << "  " << 0.337619 << endl;
		
// 		os.close();	
// 		os.open("los_const_betaphi");
// 		for(double r = -max; r<max; r+= step)
// 		{
// 			for(double z = -max; z<max; z+=step)
// 			{
// 				double hold = los->get_los_square(r, z );
// 				os << "  " << r << "  " << z << "  " << hold << "  " << sqrt(hold) << endl;	
// 			}
// 		}
// 		
	}
// 	cout << "NODE test" << endl;
	
	return 0;
}
