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
* @param n number of particles per dimension
* @param h the distance of the particles (mesh width of the grid)
* @param m the mass of the particles
* @param meanVelocity the mean-value of the velocity of the Brownian Motion
*/

	void generateCuboid(std::list<Particle> &particles, utils::Vector<double, 3> x,
							utils::Vector<double, 3> v, utils::Vector<int, 3> n,
								double h, double m, double meanVelocity);


#endif /* PARTICLEGENERATOR_H_ */
