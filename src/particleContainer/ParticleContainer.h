

#ifndef PARTICLECONTAINER_H_
#define PARTICLECONTAINER_H_

#include "../Particle.h"
#include <list>
#include <vector>

using namespace std;



/** 
* \mainpage
* \ref efficiency
*/
/**
 * @brief This class offers the possibility to to iterate over single particles or pairs of particles.
 * 
 * \section efficiency Efficiency of the ParticleContainer implementations
 * 
 * \subsection efficiencyPrerequisites Prerequisites
 * 
 * \li efficiency-tests done without output
 * \li simulated single cuboids (particle distance 1.1225) with all particles in
 * \li LinkedCellParticleContainer initialized with cutoff radius 3 (=> ~7 particles per cell) without boundary conditions were neccessary
 * \li the runtime  given in milliseconds per iteration step (Intel Core i5-2410M @ 2.30GHz)
 * 
 * \subsection efficiencySimulation Simulation
 * 
 * \image html "charts.png"
 * \image latex "charts.eps"
 * 
 * \subsection efficiencyConclusions Conclusions
 * 
 * \li for almost all count of particles: LinkedCellParticleContainer much more efficient
 * \li runtime of SimpleParticleContainer scales quadratic, of LinkedCellParticleContainer linear with the number of particles
 * 
 **/
class ParticleContainer {

public:
	/**
	* @brief Type of a list of single particles
	**/
	typedef list<Particle> SingleList;
	
	/**
	* @brief Type of a list of particles pairs
	**/
	typedef vector<pair<Particle*, Particle*> > PairList;
	
private:
	
	/**
	 * @brief Logger for ParticleContainer class
	 */
	static log4cxx::LoggerPtr logger;
	
public:
	
	/**
	* @brief Destructor of the class ParticleContainer
	*/
	virtual ~ParticleContainer();
	
	/**
	* @brief Add a list of particles to the ParticleContainer and build inner datastructures
	* 
	* @param pList function to apply on all single particles
	**/
	virtual void addParticles( list<Particle> pList ) = 0;
	
	/**
	* @brief Iterate over all sinlge particles
	* 
	* @param singleFunction function to apply to single particles
	**/
	virtual void applyToSingleParticles( void (*singleFunction)(Particle&) ) = 0;

	/**
	* @brief Iterate over particle pairs
	* 
	* @param pairFunction function to apply to the particle pairs
	**/
	virtual void applyToParticlePairs( void (*pairFunction)(Particle&, Particle&) ) = 0;
	
	/**
	 * @brief Return a list of all single particles
	 * 
	 * @return list of all single particles
	 */
	virtual SingleList& getParticles() = 0;
	
	/**
	* @brief Count of all single particles
	* 
	* @return Particle count
	**/
	virtual int size() = 0;
};


#endif /* PARTICLECONTAINER_H_ */


