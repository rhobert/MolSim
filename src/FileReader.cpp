/*
 * FileReader.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "FileReader.h"
#include "utils/Vector.h"
#include "ParticleGenerator.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>

using namespace std;

FileReader::FileReader() {
}

FileReader::~FileReader() {
}


void FileReader::readFile(std::list<Particle>& particles, char* filename) {
	double x[] = {0,0,0};
	double v[] = {1,1,1};
	double m = 1;
    int num_particles = 0;

	utils::Vector<double, 3> x;
	utils::Vector<double, 3> v;
	utils::Vector<double, 3> n;
	double m;
	double h;
	double meanVelocity;

    std::ifstream input_file(filename);
    string tmp_string;

    if (input_file.is_open()) {

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

    	for (int i = 0; i < num_particles; i++) {
    		istringstream datastream(tmp_string);

    		for (int j = 0; j < 3; j++) {
    			datastream >> x[j];

    		}
    		for (int j = 0; j < 3; j++) {
    			datastream >> v[j];
    		}
    		if (datastream.eof()) {
    			cout << "Error reading file: eof reached unexpectedly reading from line " << i << endl;
    			exit(-1);
    		}
    		datastream >> m;
    		Particle p(x, v, m);
    		particles.push_back(p);

    		getline(input_file, tmp_string);
    		cout << "Read line: " << tmp_string << endl;
    	}
		
		while (!input_file.eof()) {

			istringstream cstream(tmp_string);

			utils::Vector<double, 3> x;
			utils::Vector<double, 3> v;
			utils::Vector<int, 3> n;
			double m;
			double h;
			double meanVelocity;

			for (int j = 0; j < 3; j++) {
    			cstream >> x[j];
			}

			for (int j = 0; j < 3; j++) {
    			cstream >> v[j];
			}

			for (int j = 0; j < 3; j++) {
    			cstream >> n[j];
			}

			cstream >> m;
			cstream >> h;
			cstream >> meanVelocity;

			generateCuboid(particles, x, v, n, h, m, meanVelocity);

			getline(input_file, tmp_string);
		}

    } else {
    	std::cout << "Error: could not open file " << filename << std::endl;
    	exit(-1);
    }

}

