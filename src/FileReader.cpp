/*
* FileReader.cpp
*
* Created on: 23.02.2010
* Author: eckhardw
*/

#include "FileReader.h"
#include "utils/Vector.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>

using namespace std;

log4cxx::LoggerPtr FileReader::logger(log4cxx::Logger::getLogger("FileReader"));

FileReader::FileReader() {
}

FileReader::~FileReader() {
}


void FileReader::readFileList(list<Particle>& particles, char* filename) 
{
	double x[] = {0,0,0};
	double v[] = {1,1,1};
	double m = 1;
	int num_particles = 0;
	
	ifstream input_file(filename);
	string tmp_string;

	if (input_file.is_open()) 
	{
		getline(input_file, tmp_string);
		LOG4CXX_DEBUG(logger, "Read line: " << tmp_string);
		
		while (tmp_string.size() == 0 || tmp_string[0] == '#') {
			getline(input_file, tmp_string);
			LOG4CXX_DEBUG(logger, "Read line: " << tmp_string);
		}
		
		istringstream numstream(tmp_string);
		numstream >> num_particles;
		LOG4CXX_DEBUG(logger, "Reading " << num_particles);
		getline(input_file, tmp_string);
		LOG4CXX_DEBUG(logger, "Read line: " << tmp_string);
		
		for (int i = 0; i < num_particles; i++) 
		{
			istringstream datastream(tmp_string);

			for (int j = 0; j < 3; j++) 
			{
				datastream >> x[j];
			}
			
			for (int j = 0; j < 3; j++) 
			{
				datastream >> v[j];
			}
			
			if (datastream.eof()) 
			{
				LOG4CXX_FATAL(logger, "Error reading file: eof reached unexpectedly reading from line " << i);
				exit(-1);
			}
			
			datastream >> m;
			Particle p(x, v, m);
			particles.push_back(p);

			getline(input_file, tmp_string);
			LOG4CXX_DEBUG(logger, "Read line: " << tmp_string);
		}
	} 
	else 
	{
		LOG4CXX_FATAL(logger, "Error: could not open file " << string(filename));
		exit(-1);
	}
}

void FileReader::readFileCuboid(list<Particle>& particles, char* filename) 
{
	double x[] = {0,0,0};
	double v[] = {1,1,1};
	int N[] = {1,1,1};
	double h = 1;
	double m = 1;
	int num_particles = 0;
	
	ifstream input_file(filename);
	string tmp_string;

	if (input_file.is_open()) 
	{
		getline(input_file, tmp_string);
		LOG4CXX_DEBUG(logger, "Read line: " << tmp_string);
		
		while (tmp_string.size() == 0 || tmp_string[0] == '#') {
			getline(input_file, tmp_string);
			LOG4CXX_DEBUG(logger, "Read line: " << tmp_string);
		}
		
		istringstream numstream(tmp_string);
		numstream >> num_particles;
		LOG4CXX_DEBUG(logger, "Reading " << num_particles);
		getline(input_file, tmp_string);
		LOG4CXX_DEBUG(logger, "Read line: " << tmp_string);
		
		for (int i = 0; i < num_particles; i++) 
		{
			istringstream datastream(tmp_string);

			for (int j = 0; j < 3; j++) 
			{
				datastream >> x[j];
			}
			
			for (int j = 0; j < 3; j++) 
			{
				datastream >> v[j];
			}
			
			for (int j = 0; j < 3; j++) 
			{
				datastream >> N[j];
			}
			
			if (datastream.eof()) 
			{
				LOG4CXX_FATAL(logger, "Error reading file: eof reached unexpectedly reading from line " << i);
				exit(-1);
			}
			
			datastream >> h;
			datastream >> m;
			
			double particlePosition[3];
			particlePosition[0] = x[0];
			
			LOG4CXX_INFO(logger, "Start to generate cuboid of Particles");
			
			for ( int j0 = 0; j0 < N[0]; j0++ )
			{
				particlePosition[0] += h;
				particlePosition[1] = x[1];
				
				for ( int j1 = 0; j1 < N[1]; j1++ )
				{
					particlePosition[1] += h;
					particlePosition[2] = x[2];
					
					for ( int j2 = 0; j2 < N[2]; j2++ )
					{	
						particlePosition[2] += h;
						
						Particle p(particlePosition, v, m);
						particles.push_back(p);
					}
				}
			}
			
			LOG4CXX_INFO(logger, "Finish generate cuboid of Particles");

			getline(input_file, tmp_string);
			LOG4CXX_DEBUG(logger, "Read line: " << tmp_string);
		}
	} 
	else 
	{
		LOG4CXX_FATAL(logger, "Error: could not open file " << string(filename));
		exit(-1);
	}
}