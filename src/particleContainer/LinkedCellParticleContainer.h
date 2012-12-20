

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
	vector<int> * boundaryCells[6];
	
	/**
	 * @brief Vector with indicies to of all halo cells
	 */
	vector<int> * haloCells[6];
	
	/**
	* @brief Pairs of cells which are neighboured
	**/
	CellPairList cellPairs;
	
	/**
	* @brief Pairs of neighboured boundary periodic over the boundary
	**/
	CellPairList * periodicCellPairs[6];
	
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
	 * @brief List of single particles, since last getParticles() call
	 */
	SingleList singleList;
	
	/**
	 * @brief Get cell of position
	 * 
	 * @param x Position
	 * 
	 * @return Cell id which position is maped to and -1 if out of domain and halo
	 */
	int getCell( utils::Vector<double,3> x );
	
	/**
	 * @brief Get cell of position
	 * 
	 * @param x Position
	 * 
	 * @return Cell id which position is maped to and -1 if out of domain and halo
	 */
	int getCell( utils::Vector<int,3> x );
	
	/**
	 * @brief Check if cell is boundary cell
	 * 
	 * @param x Position
	 * 
	 * @return True if cell is boundary cell
	 */
	bool isBoundary( utils::Vector<int,3> x );
	
	/**
	 * @brief Check if cell is halo cell
	 * 
	 * @param x Position
	 * 
	 * @return True if cell is halo cell
	 */
	bool isHalo( utils::Vector<int,3> x );
	
public:
	/**
	* @brief Create an instance of LinkedCellParticleContainer which contains no particles yet
	* 
	* @param domainSize Domain size for simulation
	* @param cutoff Cutoff radius for particle pair iteration
	* @param periodic Indicates whether particle pairs are constructed above domain-boundary
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
	
	SingleList& getParticles();
	
	int size();
	
	/**
	* @brief Iterate over all single boundary particles
	* 
	* @param boundary Specification of the boundary (in 0-5)
	* @param singleFunction function to apply on all single boundary particles
	**/
	void applyToBoundaryParticles( int boundary, void (*singleFunction)(Particle&) );
	
	/**
	* @brief Iterate over all periodic boundary particle pairs
	* 
	* @param boundary Specification of the boundary (in 0-5)
	* @param singleFunction function to apply on periodic boundary particles
	**/
	void applyToPeriodicBoundaryParticlePairs( int boundary, void (*pairFunction)(Particle&, Particle) );
	
	/**
	* @brief Iterate over all single halo particles
	* 
	* @param boundary Specification of the boundary (in 0-5)
	* @param singleFunction function to apply on all single halo particles
	**/
	void applyToHaloParticles( int boundary, void (*singleFunction)(Particle&) );
	
	/**
	* @brief Move particles to their destined cells
	**/
	void updateContainingCells();
	
	/**
	 * @brief Delete all particles in the halo
	 * 
	 * @param boundary Specification of the boundary (in 0-5)
	 */
	void deleteHaloParticles( int boundary );
	
	/**
	 * @brief Returns the domain size
	 * 
	 * @return Vector with the domain size
	 */
	utils::Vector<double, 3> getDomainSize();
	
};

#endif /* LINKEDCELLPARTICLECONTAINER_H_ */


