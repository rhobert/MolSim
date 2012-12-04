 

#ifndef CELL_H_
#define CELL_H_

#include <list>
#include "Particle.h"

using namespace std;


/**
* @brief This class offers the possibility to to iterate over single particles or pairs of particles without boundary conditions.
**/
class Cell {
	
public:
	/**
	* @brief Type of a list of single particles
	**/
	typedef list<Particle> SingleList;
	
private:
	
	/**
	 * @brief Logger for Cell class
	 */
	static log4cxx::LoggerPtr logger;

private:
	
	/**
	* @brief All particles in this cell
	**/
	SingleList particles;

public:
	Cell();
	
	/**
	* @brief Add particle to cell
	* 
	* @param p Particle to add to cell
	**/
	void addParticle( Particle& p );
	
	/**
	* @brief Count of all single particles
	* 
	* @return Particle count
	**/
	int size();
	
	friend class LinkedCellParticleContainer;
};

#endif /* CELL_H_ */

