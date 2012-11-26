#ifndef PARTICLEGENERATOR_H_
#define PARTICLEGENERATOR_H_

#include "utils/Vector.h"
#include "Particle.h"
#include <list>


/**
* implements the functionality to create a cuboid where particles are arranged in a 3d rectangular grid
* 
* @param x the coordinate of the lower left front-side corner
* @param v the initial velocity of the particles
* @param N number of particles per dimension
* @param h the distance of the particles (mesh width of the grid)
* @param m the mass of the particles
*/

void generateCuboid(
	std::list<Particle> &particles, 
	utils::Vector<double, 3> x,
	utils::Vector<double, 3> v, 
	utils::Vector<int, 3> N,
	double h, 
	double m
);

/**
* implements the functionality to create a sphere of particles
* 
* @param x the coordinate of the center
* @param v the initial velocity of the particles
* @param N number of particles along the radius
* @param d indicator to which dimension the sphere expands
* @param h the distance of the particles (mesh width of the grid)
* @param m the mass of the particles
*/

void generateSphere(
	std::list<Particle> &particles, 
	utils::Vector<double, 3> x,
	utils::Vector<double, 3> v, 
	int N,
	utils::Vector<int, 3> d, 
	double h, 
	double m
);

#endif /* PARTICLEGENERATOR_H_ */
