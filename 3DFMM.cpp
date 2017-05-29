#include "3DFMM.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <deque>
#include <vector>
#include <algorithm>
#include <ctime>
#include <queue>
#include <chrono>

#define pi 3.14159265358979323846

int size_lon;
int size_lat;
int size_z;

double dlat;
double dlon;
double dz;

bool geographic = true;

int max_fdc = 2;

std::vector<double> fdc(int accuracy)
{
	std::vector<double> cooef;
	if (accuracy == 1)
	{
		cooef.push_back(1.0);
		cooef.push_back(-1.0);
	}
	else if(accuracy == 2)
	{
		cooef.push_back(3.0/2.0);
		cooef.push_back(-2.0);
		cooef.push_back(1.0/2.0);
	}
	else if(accuracy == 3)
	{
		cooef.push_back(11.0/6.0);
		cooef.push_back(-3.0);
		cooef.push_back(3.0/2.0);
		cooef.push_back(-1.0/3.0);
	}
	else
	{
		cooef.push_back(0.0);
	}
	return cooef;
}

int sub2ind(int y, int x, int z)
{
	return y + (x)*size_lat + (z)*size_lon*size_lat;
}

void ind2sub(int n, int &y, int &x, int &z)
{
	y = (n - 1) % size_lat + 1;
	n = (n - y)/size_lat +1;
	
	x = (n - 1) % size_lon + 1;
	n = (n - x)/size_lon +1;
	
	z = (n - 1) % size_z + 0;
	x--;
}

struct NB
{
	int index;
	float time;
	int prev_index;
	
	NB(int a, float b, int c) : index(a), time(b), prev_index(c) {}
};

struct CompareTime
{
	inline bool operator() (const NB& s1, const NB& s2)
	{
		return(s1.time > s2.time);
	}
};

struct Neighbour
{
	double a;
	double b;
	double c;
	
	Neighbour(double a, double b, double c) : a(a), b(b), c(c) {}
};

double deg2rad(double deg)
{
	return deg * (pi / 180);
}

bool Check(double lat, double lon, double z)
{
	if(lat < 0 || lat >= size_lat)
		return false;
	
	if(lon < 0 || lon >= size_lon)
		return false;
	
	if(z < 0 || z >= size_z)
		return false;
	
	return true;
}

Neighbour GetNeighbours(int lat, int lon, int z, int direction, float *TIME, bool *ACCEPTED, double *LAT, double *H)
{
	int dirlat = 0;
	int dirlon = 0;
	int dirz = 0;
	
	double r = 6371000 + H[z];
	
	double element = 0;
	
	if(geographic)
	{
		if(direction == 1)
		{
			dirlat = 1;
			element = r*dlat;
		}
		else if (direction == 2)
		{
			dirlon = 1;
			element = r*cos(deg2rad(LAT[lat]))*dlon;
		}
		else if ( direction == 3)
		{
			dirz = 1;
			element = dz;
		}
	}
	else
	{
		if(direction == 1)
		{
			dirlat = 1;
			element = dlat;
		}
		else if (direction == 2)
		{
			dirlon = 1;
			element = dlon;
		}
		else if ( direction == 3)
		{
			dirz = 1;
			element = dz;
		}
	}
	
	
	double up = 0;
	double xc = 0;
	
	int forward = 0;
	int backward = 0;
	
	for(int a = 1; a < max_fdc+1; a++)
	{
		//printf("?%d", a);
		if(Check(lat+dirlat*a,lon+dirlon*a,z+dirz*a))
		{
			if(ACCEPTED[sub2ind(lat+dirlat*a,lon+dirlon*a,z+dirz*a)] == 1)
			{
				forward ++;
			}
			else
			{
				break;
			}
		}
		else
		{
			break;
		}
	}
	//printf("\n");
	for(int a = 1; a < max_fdc+1; a++)
	{
		//printf("!%d", a);
		if(Check(lat-dirlat*a,lon-dirlon*a,z-dirz*a))
		{
			if(ACCEPTED[sub2ind(lat-dirlat*a,lon-dirlon*a,z-dirz*a)] == 1)
			{
				backward ++;
			}
			else
			{
				break;
			}
		}
		else
		{
			break;
		}
	}
	//printf("\n");
	
	//printf("C++ FMM f:%d, b:%d\n", forward, backward);
	
	std::vector<double> forward_cooef = fdc(forward);
	std::vector<double> backward_cooef = fdc(backward);
	
	double forward_up = 0;
	double backward_up = 0;
	
	for(int a = 1; a <= forward; a++)
	{
		//printf("f%d", a);
		forward_up += forward_cooef[a]*TIME[sub2ind(lat+dirlat*a,lon+dirlon*a,z+dirz*a)];
	}
	
	//printf("\n");
	for(int a = 1; a <= backward; a++)
	{
		//printf("b%d", a);
		backward_up += backward_cooef[a]*TIME[sub2ind(lat-dirlat*a,lon-dirlon*a,z-dirz*a)];
	}
	//printf("\n");
	
	//printf("fup = %12.8f, bup = %12.8f\n", forward_up, backward_up);
	
	if(forward_up == 0 && backward_up == 0)
	{
		//printf("returned 0, 0, 0\n");
		return Neighbour(0,0,0);
		
	}
	else if(forward_up < backward_up)
	{
		up = forward_up;
		xc = forward_cooef[0];
	}
	else
	{
		up = backward_up;
		xc = backward_cooef[0];
	}
	//printf("up = %12.8f, xc = %12.8f\n", up, xc);
	
	
	/*
	if(ACCEPTED[sub2ind(lat+dirlat,lon+dirlon,z+dirz)] == 1 && ACCEPTED[sub2ind(lat-dirlat,lon-dirlon,z-dirz)] == 1)
	{
		if(TIME[sub2ind(lat+dirlat,lon+dirlon,z+dirz)] <= TIME[sub2ind(lat-dirlat,lon-dirlon,z-dirz)])
		{
			up = TIME[sub2ind(lat+dirlat,lon+dirlon,z+dirz)];
		}
		else
		{
			up = TIME[sub2ind(lat-dirlat,lon-dirlon,z-dirz)];
		}
	}
	else if(ACCEPTED[sub2ind(lat+dirlat,lon+dirlon,z+dirz)] == 1 && ACCEPTED[sub2ind(lat-dirlat,lon-dirlon,z-dirz)] == 0)
	{
		up = TIME[sub2ind(lat+dirlat,lon+dirlon,z+dirz)];
	}
	else if(ACCEPTED[sub2ind(lat+dirlat,lon+dirlon,z+dirz)] == 0 && ACCEPTED[sub2ind(lat-dirlat,lon-dirlon,z-dirz)] == 1)
	{
		up = TIME[sub2ind(lat-dirlat,lon-dirlon,z-dirz)];
	}
	*/
	
	double a = pow(xc,2) / pow(element,2);
	double b = (2*up*xc) / pow(element,2);
	double c = pow(up,2) / pow(element,2);
	//printf("a = %12.8f, b = %12.8f, c = %12.8f\n", a, b, c);
	if(pow(up,2) > 0)
	{
		//std::cout << a << ":" << b << ":" << c<< "\t\t";
		return Neighbour(a,b,c);
		
	}
	else
	{
		//std::cout << 0 << ":" << 0 << ":" << 0 << "\t\t";
		return Neighbour(0,0,0);
		
	}
	
}

void calculate(int lat, int lon, int z, double *LAT, double *LON, double *H, float *TIME, bool *ACCEPTED, float *MODEL, std::priority_queue<NB, std::vector<NB>, CompareTime> &narrow_band, int prev_index)
{
	//printf("C++ FMM: \tcalculated: %d, %d, %d: ", lat, lon, z);
	if(Check(lat, lon, z))
	{
		double u = NAN;
		double u1, u2;
		
		if(ACCEPTED[sub2ind(lat, lon, z)] == 0 && MODEL[sub2ind(lat, lon, z)] > 0) // Nie powtarzamy obilczeń dla komórek już zaakceptowanych
		{
			if(lat >= 0 && lat < size_lat && lon >= 0 && lon < size_lon && z >= 0 && z < size_z ) // Nie liczymy dla komórek na brzegu modelu
			{
				
				Neighbour nlat = GetNeighbours(lat, lon, z, 1, TIME, ACCEPTED, LAT, H);
				Neighbour nlon = GetNeighbours(lat, lon, z, 2, TIME, ACCEPTED, LAT, H);
				Neighbour nz = GetNeighbours(lat, lon, z, 3, TIME, ACCEPTED, LAT, H);
				
				
				double a = nlat.a + nlon.a + nz.a;
				double b = nlat.b + nlon.b + nz.b;
				double c = nlat.c + nlon.c + nz.c - 1/pow(MODEL[sub2ind(lat, lon, z)],2);
				//printf("a = %24.20f, b = %24.20f, c = %24.20f, v = %24.20f\n", a, b, c, MODEL[sub2ind(lat, lon, z)]);
				double delta = sqrt(pow(b,2) - 4*a*c);
				u1 = (-b - delta) / (2*a);
				u2 = (-b + delta) / (2*a);
				u = std::max(u1,u2);
				
				if (TIME[sub2ind(lat, lon, z)] > u || TIME[sub2ind(lat, lon, z)] == 0)
				{
					TIME[sub2ind(lat, lon, z)] = u;
					narrow_band.push(NB(sub2ind(lat, lon, z), u, prev_index));
				}
			}
		}
	}
	//printf("%3.2f \n", u);
	
}


void _FMM3D(float *MODEL, float *TIME, bool *ACCEPTED, double *LAT, double *LON, double *H, int *TRACE, int N)
{
	printf("C++ FMM start\n");
	printf("C++ FMM model size: %d, %d, %d\n", size_lat, size_lon, size_z);
	printf("C++ FMM grid size: %6.4f, %6.4f, %3.2f\n", dlat, dlon, dz);

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point t2;// = std::chrono::high_resolution_clock::now();

    
	
	std::priority_queue<NB, std::vector<NB>, CompareTime> narrow_band;
	int lat, lon, z;
	
	printf("C++ FMM configuring source...\t");
	for(int i = 0; i < size_lat*size_lon*size_z; i++)
	{
		if(TIME[i] > 0)
		{
			ind2sub(i, lat, lon, z);
			narrow_band.push(NB(sub2ind(lat, lon, z), TIME[i], -1));
			//printf("C++ FMM: adding to nb: %d, %d, %d (%d)\n", lat, lon, z, i);
		}
	}
	printf("done.\n");
	
	const int debug_step = 1000000;

	while (N > 0)
	{
		if (N % debug_step == 0)
		{
			t2 = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
			printf("C++ FMM %d (duration: %10.5fs\tgps:%10.5f)\n", N, (double)duration/1000000., ((double)debug_step/(double)duration)*1000000.);
			//printf("C++ FMM %d \n", N);
			t1 = std::chrono::high_resolution_clock::now();
		}
		if(narrow_band.size() == 0)
		{
			break;
		}
		
		while(TIME[narrow_band.top().index] != narrow_band.top().time || ACCEPTED[narrow_band.top().index]==1)
		{
			//printf("C++ FMM: popping:(%d)\n", narrow_band.top().index);
			narrow_band.pop();
			if(narrow_band.size() == 0)
			{
				break;
			}
		}
		
		
		
		ind2sub(narrow_band.top().index, lat, lon, z);
		ACCEPTED[narrow_band.top().index]	=	1;
		TIME[narrow_band.top().index]	=	narrow_band.top().time;
		TRACE[narrow_band.top().index]	=	narrow_band.top().prev_index;
		//printf("C++ FMM: calculating neighbours of: %d, %d, %d (%d, %d) - %d\n", lat, lon, z, narrow_band.top().index, sub2ind(lat, lon, z), N);
		narrow_band.pop();
		
		
		calculate(lat + 1,	lon,		z,		LAT, LON, H, TIME, ACCEPTED, MODEL, narrow_band, sub2ind(lat, lon, z));
		calculate(lat - 1,	lon,		z,		LAT, LON, H, TIME, ACCEPTED, MODEL, narrow_band, sub2ind(lat, lon, z));
		calculate(lat,		lon + 1,	z,		LAT, LON, H, TIME, ACCEPTED, MODEL, narrow_band, sub2ind(lat, lon, z));
		calculate(lat,		lon - 1,	z,		LAT, LON, H, TIME, ACCEPTED, MODEL, narrow_band, sub2ind(lat, lon, z));
		calculate(lat,		lon,		z + 1,	LAT, LON, H, TIME, ACCEPTED, MODEL, narrow_band, sub2ind(lat, lon, z));
		calculate(lat,		lon,		z - 1,	LAT, LON, H, TIME, ACCEPTED, MODEL, narrow_band, sub2ind(lat, lon, z));
		
		
		
		
		N--;
	}
	printf("C++ FMM end\n");
}

void _SetModelSize(int _size_lat, int _size_lon, int _size_z, double _dlat, double _dlon, double _dz)
{
	size_lat = _size_lat;
	size_lon = _size_lon;
	size_z = _size_z;
	
	dlat = _dlat;
	dlon = _dlon;
	dz = _dz;
}