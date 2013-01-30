
#ifndef STATISTICS_H_
#define STATISTICS_H_

#include <list>

#include "Particle.h"
#include "particleContainer/ParticleContainer.h"

using namespace std;

/**
 * @brief This class represents statistcs of particles (movement, distribution)
 */
class Statistics
{
	/**
	 * @brief Logger for Statistics class
	 */
	static log4cxx::LoggerPtr logger;
	
	/**
	* @brief ParticleContainer to which the statistcs refers
	*/
	ParticleContainer* pc;
	
	/**
	* @brief Count of particle pairs for the RDF for the diffrent distances
	*/
	static int* rdfStatistics;
	
	/**
	* @brief The step size of the intervals in which the radius of the distance for the RDF are divided
	*/
	static double rdfDeltaR;
	
	/**
	* @brief Saves the current position of the particle
	* 
	* @param p Particle which position should be saved
	*/
	static void saveOldX(Particle& p);
	
	/**
	* @brief Updates rdfStatistics by calculating the distance of two particles
	* 
	* @param p1 First particle
	* @param p2 Second particle
	*/
	static void updateRDFStatistics(Particle& p1, Particle& p2);
	
public:
	
	/**
	* @brief Creates an instance of Statistics
	* 
	* @param pc ParticleContainer to which the statistcs refers
	*/
	Statistics(ParticleContainer& pc);
	
	/**
	* @brief Start to calculate the diffusion by saving the current position of all particles
	*/
	void beginCalcDiffusion();
	
	/**
	* @brief Finish the calculation of the diffusion (requires beginCalcDiffusion called before)
	* 
	* @return The diffusion
	*/
	double endCalcDiffusion();
	
	/**
	* @brief Calculates the RDF by dividing the cutoff radius in a certain count of intervals 
	* 
	* @param intervalCount Count of the intervals the cutoff radius should be divided
	* @param cutoff The cutoff radius (maximal radius which should be analyzed)
	* 
	* @return The RDF for all intervals
	*/
	double* calcRDF(int intervalCount, double cutoff);
};

#endif /* STATISTICS_H_ */


