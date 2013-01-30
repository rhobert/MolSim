
#ifndef MEMBRANEPARTICLE_H_
#define MEMBRANEPARTICLE_H_

#include <vector>

#include "utils/Vector.h"

#include "Particle.h"
#include "Membrane.h"

using namespace std;

/**
 * @brief This class represents a particle of a membrane
 */
class MembraneParticle: public Particle 
{
	/**
	* @brief Type of vector of MembraneParticles
	*/
	typedef vector<MembraneParticle*> NeighbourList;
	
	/**
	* @brief Pointer to the membrane associated with this particle
	*/
	Membrane* membrane;
	
	/**
	* @brief Lists of neighboured particles classified by their type of neighbourhood (direct: 0, diagonal: 1, diagonal over two dimensions: 2)
	*/
	NeighbourList neighbours[3];
	
public:
	
	/**
	* @brief Create an instance of MembraneParticle
	* 
	* @param x_arg Position of the particle
	* @param v_arg Velocity of the particle
	* @param m_arg Mass of the particle
	* @param sigma Sigma parameter for Lennard-Jones potential of the particle
	* @param epsilon Epsilon parameter for Lennard-Jones potential of the particle
	* @param type Type of the particle
	* @param membrane Pointer to the membrane associated with this particle
	*/
	MembraneParticle(
		utils::Vector<double, 3> x_arg,
		utils::Vector<double, 3> v_arg,
		double m_arg,
		double sigma,
		double epsilon,
		int type,
		Membrane* membrane
	);
	
	/**
	 * Destructs this MembraneParticle. This particle will be removed from the neighbourlists of his neighbours
	 */
	virtual ~MembraneParticle();
	
	/**
	 * @brief Returns the stiffness constant of the membrane associated with this particle
	 * 
	 * @return the stiffness constant
	 */
	double getStiffnessConstant();
	
	/**
	 * @brief Returns the average bond length of the membrane associated with this particle
	 * 
	 * @return the average bond length
	 */
	double getAverageBondLength();
	
	/**
	 * @brief Returns the average bond length of the membrane adapted to the type of neighbourhood (direct,diagonal,diagonal over two dimensions)
	 * 
	 * @param type type of neighbourhood (direct: 0, diagonal: 1, diagonal over two dimensions: 2)
	 * 
	 * @return the typed average bond length
	 */
	double getAverageBondLengthTyped( int type );
	
	/**
	 * @brief Add a neighbour to this particle
	 * 
	 * @param type type of neighbourhood (direct: 0, diagonal: 1, diagonal over two dimensions: 2)
	 * @param particle the new neighboured particle
	 */
	void addNeighbour( int type, MembraneParticle* particle );
	
	/**
	 * @brief Removes a particle from the neighbourlist of this type (direct: 0, diagonal: 1, diagonal over two dimensions: 2)
	 * 
	 * @param type type of neighbourhood (direct: 0, diagonal: 1, diagonal over two dimensions: 2)
	 * @param particle the particle to remove from the neighbourlist
	 */
	void removeNeighbour( int type, MembraneParticle* particle );
	
	/**
	 * @brief Iterate over all neighboured particles and applies a function on the two particles with the type of neighbourhood (direct: 0, diagonal: 1, diagonal over two dimensions: 2)
	 * 
	 * @param pairFunction Function to apply on the particle pair
	 */
	void applyToNeighbours( void (*pairFunction)(int, MembraneParticle&, MembraneParticle&) );
	
	/**
	 * @brief Checks if a particle is neighboured to this particle
	 * 
	 * @param particle Particles which neighbourhood is checked
	 * 
	 * @return true if the particles are neighboured
	 */
	bool isNeighbour( MembraneParticle* particle );
	
	/**
	 * @brief Checks if two particles are neighboured
	 * 
	 * @param p1 first particle
	 * @param p2 second particle
	 * 
	 * @return true if the particles are neighboured
	 */
	static bool neighboured( MembraneParticle* p1, MembraneParticle* p2 );
	
	/**
	 * @brief Checks if two particles are on the same membrane
	 * 
	 * @param p1 first particle
	 * @param p2 second particle
	 * 
	 * @return true if the particles are on the same membrane
	 */
	static bool sameMembrane( MembraneParticle* p1, MembraneParticle* p2 );
	
	/**
	 * @brief Cast a Particle to a MembraneParticle
	 * 
	 * @param p the particle to cast
	 * 
	 * @return casted pointer to the MembraneParticle (NULL-pointer if impossible)
	 */
	static MembraneParticle* castMembraneParticle( Particle& p );
};

inline double MembraneParticle::getStiffnessConstant()
{
	return membrane->stiffnessConstant;
}
	
inline double MembraneParticle::getAverageBondLength()
{
	return membrane->averageBondLength;
}

inline double MembraneParticle::getAverageBondLengthTyped( int type )
{
	return membrane->averageBondLengthTyped[type];
}

inline bool MembraneParticle::neighboured( MembraneParticle* p1, MembraneParticle* p2 )
{
	for ( int j = 0; j < 3; j++ )
	{
		for ( NeighbourList::iterator i = p1->neighbours[j].begin(); i != p1->neighbours[j].end(); i++ )
		{
			if ( p2 == *i )
				return true;
		}
	}
	
	return false;
}

inline bool MembraneParticle::sameMembrane( MembraneParticle* p1, MembraneParticle* p2 )
{
	return p1->membrane = p2->membrane;
}


inline MembraneParticle* MembraneParticle::castMembraneParticle( Particle& p )
{
	return dynamic_cast<MembraneParticle*>(&p);
}

#endif /* MEMBRANEPARTICLE_H_ */
