/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "utils/Vector.h"
#include "log4cxx/logger.h"

/**
* @brief This class represents a particle
**/
class Particle {

private:

	/**
	 * @brief Logger for Particle class
	 */
	static log4cxx::LoggerPtr logger;
	
	/**
	* @brief The position of the particle.
	*/
	utils::Vector<double, 3> x;

	/**
	* @brief The velocity of the particle.
	*/
	utils::Vector<double, 3> v;

	/**
	* @brief The force effective on this particle.
     */
	utils::Vector<double, 3> f;

	/**
	* @brief The force wich was effective on this particle.
	*/
	utils::Vector<double, 3> old_f;

	/**
	* @brief The mass of this particle.
    */
	double m;

	/**
	* @brief Type of the particle. Use it for whatever you want (e.g. to seperate
     * molecules belonging to different bodies, matters, and so on).
    */
	int type;

public:
/**
* @brief Create an instance of Particle.
*/
	Particle(int type = 0);
/**
* @brief Copy-constructor of the class Particle.
*/
	Particle(const Particle& other);
/**
* @brief Create an instance of Particle using parameters containing position, velocity and mass.
*/
	Particle(
			// for visualization, we need always 3 coordinates
			// -> in case of 2d, we use only the first and the second
			utils::Vector<double, 3> x_arg,
			utils::Vector<double, 3> v_arg,
			double m_arg,
			int type = 0
	);

/**
* @brief Destructor of the class Particle.
*/
	virtual ~Particle();

/**
* @brief Returns the position of the particle.
*/
	utils::Vector<double, 3>& getX();

/**
* @brief Returns the force effective on the particle.
*/
	utils::Vector<double, 3>& getF();

/**
* @brief Returns the force which was effective on the particle.
*/
	utils::Vector<double, 3>& getOldF();

/**
* @brief Returns the velocity of the particle.
*/
	utils::Vector<double, 3>& getV();

/**
* @brief Changes the position of the particle.
*/
	void setX(utils::Vector<double, 3> other_x);

/**
* @brief Changes the velocity of the particle.
*/
	void setV(utils::Vector<double, 3> other_v);

/**
* @brief Changes the forces effective on the particle.
*/
	void setF(utils::Vector<double, 3> other_f);
	
/**
* @brief Changes the forces effective on the particle and saves the old forces effective on the particle.
*/
	void newF(utils::Vector<double, 3> other_f);

/**
* @brief Returns the mass of the particle.
*/
	double getM();

/**
* @brief Returns the type of the particle.
*/
	int getType();

/**
* @brief Compares two particles.
*/
	bool operator==(Particle& other);

/**
* @brief Creates a string-representation of the particle.
*/
	std::string toString();
};

std::ostream& operator<<(std::ostream& stream, Particle& p);

#endif /* PARTICLE_H_ */
