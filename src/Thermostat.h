
#ifndef THERMOSTAT_H_
#define THERMOSTAT_H_

#define kB 1.0
#define thermostat_mass 1.0

#include "particleContainer/ParticleContainer.h"
#include "utils/Vector.h"
#include <cmath>

/**
 * @brief This class offers the possibility to simulate temperature adjustment control
 */
class Thermostat {

private:

    /**
	 * @brief Logger for Thermostat class
	 */
	static log4cxx::LoggerPtr logger;
	
	/**
	 * @brief ParticleContainer to regulate
	 */
	ParticleContainer * pc;
	
	/**
	 * @brief The desired temperature
	 */
	double targetT;
	
	/**
	 * @brief Count of dimensions for simulation 
	 */
	int dimensionCount;
	
	/**
	 * @brief Frequency to regulate temperature
	 */
	int regulationFrequency;
	
	/**
	 * @brief Step size in which temperature should be changed
	 */
	double diffT;
	
	/**
	 * @brief The number of timesteps after which temperature has to be changed
	 */
	int nMax; 
	
	/**
	 * @brief Kinetic energy which is needed to reach the disered temperature
	 */
	//double targetEnergy;
	
	/**
	 * @brief Kinetic energy of all particles at a particluar timestep 
	 */
	double currentEnergy;
	/**
	 * @brief Temperature at a particluar timestep
	 */
	double currentT;
	
	/**
	* @brief Apply Maxwell Boltzmann Distribution
	*
	* @param p Particle to distribute
	**/
	static void applyMaxwellBoltzmannDistribution( Particle& p );
	
	/**
	* @brief Scales veclocity of particle by static parameter beta
	* 
	* @param p Particle to scale velocity
	*/
	static void scaleVelocity( Particle& p );
	
	/**
	 * @brief Parameter for Maxwell Boltzmann Distribution
	 */
	static double meanVelocity;
	
	/**
	 * @brief Parameter for Maxwell Boltzmann Distribution
	 */
	static int dimensions;
	
	/**
	 * @brief Scale factor for velocity
	 */
	static double beta;
	
public:
	
	/**
	* @brief Create an instance of class Thermostat
	*
	* @param pc ParticleContainer which should be regulated
	* @param initialT Initial temperature for Maxwell-Boltzmann-Distribution
	* @param dimensionCount count of dimensions for simulation
	**/
	Thermostat( ParticleContainer& pc, double initialT, int dimensionCount );
	
	/**
	 * @brief Set the frequency to regulate temperature
	 * 
	 * @param frequency Frequency to regulate temperature
	 */
	void setFrequency ( int frequency );
	
	/**
	* @brief Checks if temperature should be regulated and do it if neccessary
	*
	* @param iteration Current iteration
	*/
	void apply ( int iteration );
	
	/**
	* @brief Regulates temperature
	*
	* @param targetT the disered temperature
	*/
	void regulateTemperature ( double targetT );
	
	/**
	 * @brief Returns the current energy
	 * 
	 * @return current energy
	 **/
	double getEnergy();
	
	/**
	 * @brief Returns the current temperature
	 * 
	 * @return current temperature
	 **/
	double getTemperature();
	
	/**
	* @brief Calculates the initial scaling value for particle velocities to reach initial temperature
	*
	* @param initialT initial temperature
	*
	* @return mean velocity for Maxwell-Boltzmann distribution
	*/
	static double initializeTemperature( double initialT );
};

#endif /* THERMOSTAT_H_ */
