
#ifndef THERMOSTAT_H_
#define THERMOSTAT_H_

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
	 * @brief The desired temperature
	 */
	double targetT;
	/**
	 * @brief Step size in which temperature should be changed
	 */
	double diffT;
	/**
	 * @brief The number of timesteps after which temperature has to be changed
	 */
	int nMax; 
	/**
	 * @brief Count of dimensions for simulation 
	 */
	int dimensionCount;
	/**
	 * @brief Kinetic energy of all particles at simulation start
	 */
	double initialEnergy;
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
	 * @brief Boltzmann-Konstante
	 */
	double k;
	/**
	 * @brief Scaling factor for particle velocities
	 */
	double beta;
	/**
	 * @brief Mass of the particles
	 */
	double m;

public:

	/**
	* @brief Create an instance of class Thermostat
	*
	* @param targetT the disered temperature
	* @param diffT step size in which temperature should be changed
	* @param nMax the number of timesteps after which temperature has to be changed
	* @param dimensionCount count of dimensions for simulation
	**/
	Thermostat( double targetT, double diffT, int nMax , int dimensionCount );

	/**
	* @brief Calculates the initial scaling value for particle velocities to reach initial temperature
	*
	* @param size count of all particles
	* @param initialT initila temperature
	*
	* @return mean velocity for Maxwell-Boltzmann distribution
	*/
	double initializeTemperature( int size , double initialT );

	/**
	* @brief Regulates temperature
	*
	* @param pc all particles
	* @param iteration current iteration
	* @param nThermostat number of iterations after which the thermostat is applied
	*/
	void regulateTemperature( ParticleContainer& pc, int iteration, int nThermostat );

};

#endif /* THERMOSTAT_H_ */
