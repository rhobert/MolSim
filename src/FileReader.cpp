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
		cout << "Read line: " << tmp_string << endl;
		
		while (tmp_string.size() == 0 || tmp_string[0] == '#') {
			getline(input_file, tmp_string);
			cout << "Read line: " << tmp_string << endl;
		}
		
		istringstream numstream(tmp_string);
		numstream >> num_particles;
		cout << "Reading " << num_particles << "." << endl;
		getline(input_file, tmp_string);
		cout << "Read line: " << tmp_string << endl;
		
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
				cout << "Error reading file: eof reached unexpectedly reading from line " << i << endl;
				exit(-1);
			}
			
			datastream >> m;
			Particle p(x, v, m);
			particles.push_back(p);

			getline(input_file, tmp_string);
			cout << "Read line: " << tmp_string << endl;
		}
	} 
	else 
	{
		cout << "Error: could not open file " << filename << endl;
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
		cout << "Read line: " << tmp_string << endl;
		
		while (tmp_string.size() == 0 || tmp_string[0] == '#') {
			getline(input_file, tmp_string);
			cout << "Read line: " << tmp_string << endl;
		}
		
		istringstream numstream(tmp_string);
		numstream >> num_particles;
		cout << "Reading " << num_particles << "." << endl;
		getline(input_file, tmp_string);
		cout << "Read line: " << tmp_string << endl;
		
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
				cout << "Error reading file: eof reached unexpectedly reading from line " << i << endl;
				exit(-1);
			}
			
			datastream >> h;
			datastream >> m;
			
			double particlePosition[3];
			particlePosition[0] = x[0];
			
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

			getline(input_file, tmp_string);
			cout << "Read line: " << tmp_string << endl;
		}
	} 
	else 
	{
		cout << "Error: could not open file " << filename << endl;
		exit(-1);
	}
}