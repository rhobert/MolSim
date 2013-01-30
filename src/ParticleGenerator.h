#ifndef PARTICLEGENERATOR_H_
#define PARTICLEGENERATOR_H_

#include "utils/Vector.h"
#include "Particle.h"
#include "MembraneParticle.h"
#include <list>


/**
* implements the functionality to create a cuboid where particles are arranged in a 3d rectangular grid
* 
* @param particles List of particles to append the generated particles
* @param x the coordinate of the lower left front-side corner
* @param v the initial velocity of the particles
* @param N number of particles per dimension
* @param h the distance of the particles (mesh width of the grid)
* @param m the mass of the particles
* @param sigma Sigma parameter for Lennard-Jones potential of the particles
* @param epsilon Epsilon parameter for Lennard-Jones potential of the particles
* @param type Type of the particles
*/

void generateCuboid(
	std::list<Particle*> &particles, 
	utils::Vector<double, 3> x,
	utils::Vector<double, 3> v, 
	utils::Vector<int, 3> N,
	double h, 
	double m,
	double sigma = PARTCLE_SIGMA,
	double epsilon = PARTCLE_EPSILON,
	int type = PARTCLE_TYPE
	
);

/**
* @brief implements the functionality to create a membran where particles are arranged in a 3d rectangular grid
* 
* @param particles List of particles to append the generated particles
* @param x the coordinate of the lower left front-side corner
* @param v the initial velocity of the particles
* @param N number of particles per dimension
* @param h the distance of the particles (mesh width of the grid)
* @param m the mass of the particles
* @param sigma Sigma parameter for Lennard-Jones potential of the particles
* @param epsilon Epsilon parameter for Lennard-Jones potential of the particles
* @param type Type of the particles
* @param stiffnessConstant the stiffness constant of particle pair
* @param averageBondLength the average bond length of a particle pair
* @param positions list of grid positions on, which particles are additionally stored
* @param positionParticles list of a additionally stored particles which a picked via the position list
*/

void generateCuboidMembrane(
	std::list<Particle*> &particles, 
	utils::Vector<double, 3> x,
	utils::Vector<double, 3> v, 
	utils::Vector<int, 3> N,
	double h, 
	double m,
	double sigma,
	double epsilon,
	int type,
	double stiffnessConstant,
	double averageBondLength,
	std::list<utils::Vector<int,3> >& positions,
	std::list<MembraneParticle*>& positionParticles
);

/**
* @brief implements the functionality to create a sphere of particles
* 
* @param particles List of particles to append the generated particles
* @param x the coordinate of the center
* @param v the initial velocity of the particles
* @param N number of particles along the radius
* @param d dimension count to which the sphere should expand
* @param h the distance of the particles (mesh width of the grid)
* @param m the mass of the particles
* @param sigma Sigma parameter for Lennard-Jones potential of the particles
* @param epsilon Epsilon parameter for Lennard-Jones potential of the particles
* @param type Type of the particles
*/

void generateSphere(
	std::list<Particle*> &particles, 
	utils::Vector<double, 3> x,
	utils::Vector<double, 3> v, 
	int N,
	int d, 
	double h, 
	double m,
	double sigma = PARTCLE_SIGMA,
	double epsilon = PARTCLE_EPSILON,
	int type = PARTCLE_TYPE
);

#endif /* PARTICLEGENERATOR_H_ */
