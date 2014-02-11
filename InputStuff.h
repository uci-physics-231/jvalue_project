#ifndef INPUTSTUFF_H
#define INPUTSTUFF_H

// #include <cmath>
// #include "GregsMathHdr.h"
// #include "AndrewsMathHdr.h"
// #include "random.h"

#include <iostream>
#include <fstream>

using namespace std;

class inputcrap
{
	protected:
		double **datashit;
		
		int numOfData;
		int cols;
		
	public:
		inputcrap(string name, int cols_in, int bad_int)
		{
			ifstream in(name.c_str());
			cols = cols_in;
			in >> numOfData;
			datashit = new double*[cols];
			for(int i =0;i<cols;i++)
				datashit[i] = new double[numOfData];
			for (int i = 0; i < numOfData; i++)
			{
				for(int j=0;j<cols;j++)
				{
					double temp;
					string temp2;
					if(j==bad_int)
						in>> temp2;
					else 
						in >> temp;
					datashit[j][i] = temp;
// 					in >>  ;
				}
			}
		}
		inputcrap(string name, int cols_in)
		{
			ifstream in(name.c_str());
			cols = cols_in;
			in >> numOfData;
			datashit = new double*[cols];
			for(int i =0;i<cols;i++)
				datashit[i] = new double[numOfData];
			for (int i = 0; i < numOfData; i++)
			{
				for(int j=0;j<cols;j++)
				{
					double temp;
// 					string temp2;
// 					if(j==0)
// 						in>> temp2;
// 					else 
						in >> temp;
					datashit[j][i] = temp;
// 					in >>  ;
				}
			}
		}
		~inputcrap()
		{
			for(int i =0;i<cols;i++)
				delete[] datashit[i];
			delete [] datashit;
		}
};

class inputcrap_vectors
{
	protected:
		vector < vector <double> >datashit;
		
		int numOfData;
		int cols;
		
	public:
		inputcrap_vectors(string name, int cols_in)
		{
			ifstream in(name.c_str());
			cols = cols_in;
			in >> numOfData;
			cout << "numOfData  " << numOfData << endl;
// 			datashit = new double*[cols]
// 			for(int i =0;i<cols;i++)
// 				datashit[i] = new double[numOfData]
			for (int i = 0; i < numOfData; i++)
			{
				vector <double> TEMP;
				for(int j=0;j<cols;j++)
				{
					double temp;
// 					string temp2;
// 					if(j==0)
// 						in>> temp2;
// 					else 
						in >> temp;
					TEMP.push_back(temp);
// 					cout << i << "  " << j << "   " << temp << endl;
				}
				datashit.push_back(TEMP);
			}
		}
		inputcrap_vectors(string name, int cols_in, int bad_int)
		{
			ifstream in(name.c_str());
			cols = cols_in;
			in >> numOfData;
			cout << "numOfData  " << numOfData << endl;
// 			datashit = new double*[cols]
// 			for(int i =0;i<cols;i++)
// 				datashit[i] = new double[numOfData]
			for (int i = 0; i < numOfData; i++)
			{
				vector <double> TEMP;
				for(int j=0;j<cols;j++)
				{
					double temp;
					string temp2;
					if(j==bad_int)
						in>> temp2;
					else 
						in >> temp;
					TEMP.push_back(temp);
// 					cout << i << "  " << j << "   " << temp << endl;
				}
				datashit.push_back(TEMP);
			}
		}
		~inputcrap_vectors()
		{
// 			for(int i =0;i<cols;i++)
// 				delete[] datashit[i];
// 			delete [] datashit;
		}
};

class obs_parameters
{
protected:
	double obs_distance, obs_distance_error;
	double obs_rhalf_2d, obs_rhalf_2d_error;
	double obs_rcore, obs_rcore_error;
	double obs_rtidal, obs_rtidal_error;
	double obs_position_angle, obs_position_angle_error;
	double obs_axis_ratio, obs_axis_ratio_error;
	double obs_ra_center, obs_ra_center_error;
	double obs_dec_center, obs_dec_center_error;
	double obs_m_v, obs_m_v_error;
	double obs_m_l_stellar;
	
public:
	obs_parameters(string name)
	{
		ifstream ifs(name.c_str());
		ifs >> obs_distance;
		ifs >> obs_distance_error;
		ifs >> obs_rhalf_2d;
		ifs >> obs_rhalf_2d_error;
		ifs >> obs_rcore; 
		ifs >> obs_rcore_error;
		ifs >> obs_rtidal;
		ifs >> obs_rtidal_error;
		ifs >> obs_position_angle;
		ifs >> obs_position_angle_error;
		ifs >> obs_axis_ratio;
		ifs >> obs_axis_ratio_error;
		ifs >> obs_ra_center;
		ifs >> obs_ra_center_error;
		ifs >> obs_dec_center;
		ifs >> obs_dec_center_error;
		ifs >> obs_m_v;
		ifs >> obs_m_v_error;
		ifs >> obs_m_l_stellar;
	}
	obs_parameters(double * name)
	{
		obs_distance = name[0];
		obs_distance_error= name[1];
		obs_rhalf_2d = name[2];
		obs_rhalf_2d_error= name[3];
		obs_rcore= name[4];
		obs_rcore_error= name[5];
		obs_rtidal= name[6];
		obs_rtidal_error= name[7];
		obs_position_angle= name[8];
		obs_position_angle_error= name[9];
		obs_axis_ratio= name[10];
		obs_axis_ratio_error= name[11];
		obs_ra_center= name[12];
		obs_ra_center_error= name[13];
		obs_dec_center= name[14];
		obs_dec_center_error= name[15];
		obs_m_v= name[16];
		obs_m_v_error= name[17];
		obs_m_l_stellar= name[18];
	}
	obs_parameters(vector<double> name)
	{
		obs_distance = name[0];
		obs_distance_error= name[1];
		obs_rhalf_2d = name[2];
		obs_rhalf_2d_error= name[3];
		obs_rcore= name[4];
		obs_rcore_error= name[5];
		obs_rtidal= name[6];
		obs_rtidal_error= name[7];
		obs_position_angle= name[8];
		obs_position_angle_error= name[9];
		obs_axis_ratio= name[10];
		obs_axis_ratio_error= name[11];
		obs_ra_center= name[12];
		obs_ra_center_error= name[13];
		obs_dec_center= name[14];
		obs_dec_center_error= name[15];
		obs_m_v= name[16];
		obs_m_v_error= name[17];
		obs_m_l_stellar= name[18];
	}
	~obs_parameters(){}
};


#endif
