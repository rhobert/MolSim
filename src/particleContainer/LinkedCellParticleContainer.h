

#ifndef LINKEDCELLPARTICLECONTAINER_H_
#define LINKEDCELLPARTICLECONTAINER_H_

#include "Particle.h"
#include "ParticleContainer.h"
#include "Cell.h"
#include "utils/Vector.h"
#include <list>
#include <vector>

using namespace std;


/**
* @brief This class offers the possibility to to iterate over single particles or pairs of particles without boundary conditions.
**/
class LinkedCellParticleContainer: public ParticleContainer {
	
public:
	
	/**
	* @brief Type of a vector of cells
	**/
	typedef vector<Cell> CellList;
	
	/**
	* @brief Type of a vector of cell pairs
	**/
	typedef vector<pair<Cell*, Cell*> > CellPairList;
	
	
private:
	
	/**
	 * @brief Logger for ParticleContainer class
	 */
	static log4cxx::LoggerPtr logger;
	
	/**
	 * @brief Vector with all cells
	 */
	CellList cells;
	
	/**
	 * @brief Vector with indicies to of all boundary cells
	 */
	vector<int>  boundaryCells;
	
	/**
	 * @brief Vector with indicies to of all halo cells
	 */
	vector<int> haloCells;
	
	/**
	* @brief Pairs of cells which are neighboured
	**/
	CellPairList cellPairs;
	
	/**
	 * @brief Count of all cells
	 */
	int cellCount;
	
	/**
	 * @brief Side length of the cells
	 */
	double sideLength;
	
	/**
	 * @brief Vector of count of cells
	 */
	utils::Vector<int,3> cellDimensions;
	
	/**
	 * @brief Vector of domain size
	 */
	utils::Vector<double, 3> domainSize;
	
	/**
	 * @brief Count of all single particles
	 */
	int count;
	
	/**
	 * @brief Get cell of position
	 * 
	 * @param x Position
	 * 
	 * @return Cell id which position is maped to and -1 when out of domain and halo
	 */
	int getCell( utils::Vector<double,3> x );
	
	/**
	 * @brief Get cell of position
	 * 
	 * @param x Position
	 * 
	 * @return Cell id which position is maped to and -1 when out of domain and halo
	 */
	int getCell( utils::Vector<int,3> x );
	
public:
	/**
	* @brief Create an instance of LinkedCellParticleContainer which contains no particles yet
	* 
	* @param domainSize Domain size for simulation
	* @param cutoff Cutoff radius for particle pair iteration
	**/
	LinkedCellParticleContainer(utils::Vector<double, 3> domainSize, double cutoff);
	
	void addParticles( list<Particle> pList );
	
	void applyToSingleParticles( void (*singleFunction)(Particle&) );
	
	/**
	* @brief Iterate over all particle pairs with distance less equal the cutoff radius
	* 
	* @param pairFunction function to apply to particle pairs
	**/
	void applyToParticlePairs( void (*pairFunction)(Particle&, Particle&) );
	
	SingleList getParticles();
	
	int size();
	
	/**
	* @brief Iterate over all single boundary particles
	* 
	* @param singleFunction function to apply on all single boundary particles
	**/
	void applyToBoundaryParticles( void (*singleFunction)(Particle&) );
	
	/**
	* @brief Move particles to their destined cells
	**/
	void updateContainingCells();
	
	/**
	 * @brief Delete all particles in the halo
	 */
	void deleteHaloParticles();
	
	/**
	 * @brief Returns the domain size
	 * 
	 * @return Vector with the domain size
	 */
	utils::Vector<double, 3> getDomainSize();
	
};

#endif /* LINKEDCELLPARTICLECONTAINER_H_ */


