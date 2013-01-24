

#ifndef CELL_H_
#define CELL_H_

#include <list>
#include <omp.h>
#include "../Particle.h"

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
	
	/**
	* @brief Type of a list of cells
	**/
	typedef list<Cell*> CellList;
	
private:
	
	/**
	 * @brief Logger for Cell class
	 */
	static log4cxx::LoggerPtr logger;
	
	/**
	* @brief All particles in this cell
	**/
	SingleList particles;
	
	/**
	* @brief Neighboured cells
	**/
	CellList neighbours;
	
	/**
	 * @brief 
	**/
	
	/**
	 * @brief Lock for access to particles
	**/
	omp_lock_t lock;
	
public:
	/**
	 * @brief Create an instance of Cell which contains no particles yet
	 */
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


