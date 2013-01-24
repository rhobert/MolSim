/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <vector>
#include "utils/Vector.h"
#include "log4cxx/logger.h"

#define PARTCLE_EPSILON 5
#define PARTCLE_SIGMA 1
#define PARTCLE_TYPE 0

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
	* @brief The force which static effective on this particle
     */
	utils::Vector<double, 3> static_f;
	
	/**
	* @brief The mass of this particle.
    */
	double m;

	/**
	* @brief
    */
	double sigma;

	/**
	* @brief
    */
	double epsilon;

	/**
	* @brief Type of the particle. Use it for whatever you want (e.g. to seperate
    * molecules belonging to different bodies, matters, and so on).
    */
	int type;

//	/**
//	* @brief Contains all neighboring particles of a particle
//    */
//	std::vector<Particle*> neighbors;

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
* @brief Create an instance of Particle using parameters containing position, velocity, mass, sigma, epsilon and type.
*/
	Particle(
			// for visualization, we need always 3 coordinates
			// -> in case of 2d, we use only the first and the second
			utils::Vector<double, 3> x_arg,
			utils::Vector<double, 3> v_arg,
			double m_arg,
			double sigma = PARTCLE_SIGMA,
			double epsilon = PARTCLE_EPSILON,
			int type = PARTCLE_TYPE
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
	void newF();

/**
* @brief Returns the mass of the particle.
*/
	double getM();

/**
* @brief
*/
	double getSigma();

/**
* @brief
*/
	double getEpsilon();

/**
* @brief Returns the type of the particle.
*/
	int getType();
	
/**
* @brief Returns the static force effective on the particle.
*/
	utils::Vector<double, 3>& getStaticF();
	
/**
* @brief Changes the static force effective on the particle
*/	
	 void setStaticF(utils::Vector<double, 3> other_f);
	 
/**
* @brief Returns the type of the particle.
*/


/**
	* @brief Contains all neighboring particles of a particle
    */
	std::vector<Particle*> neighbors;

/**
* @brief Returns neighboring particles
*/
//	std::vector<Particle&> getNeighbors();

/**
* @brief Add neighboring particle
*/
//	void addNeighbor(Particle& p);

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

inline  utils::Vector<double, 3>& Particle::getX() {
	return x;
}

inline utils::Vector<double, 3>& Particle::getV() {
	return v;
}

inline  utils::Vector<double, 3>& Particle::getF() {
	return f;
}

inline  utils::Vector<double, 3>& Particle::getOldF() {
	return old_f;
}

inline  utils::Vector<double, 3>& Particle::getStaticF() {
	return static_f;
}


inline  void Particle::setX(utils::Vector<double, 3> other_x) {
	x = other_x;
}

inline  void Particle::setV(utils::Vector<double, 3> other_v) {
	v = other_v;
}

inline  void Particle::setF(utils::Vector<double, 3> other_f) {
	f = other_f;
}

inline  void Particle::setStaticF(utils::Vector<double, 3> other_f) {
	static_f = other_f;
}

inline  void Particle::newF() {
	old_f = f;
	f = static_f;
}

inline  double Particle::getM() {
	return m;
}

inline  double Particle::getSigma() {
	return sigma;
}

inline  double Particle::getEpsilon() {
	return epsilon;
}

inline  int Particle::getType() {
	return type;
}

//inline std::vector<Particle&>& Particle::getNeighbors() {
//   return neighbors;
//}

//inline void Particle::addNeighbor(Particle& p) {
//    neighbors->push_back(p);
//}

inline bool Particle::operator ==(Particle& other) {
	if ( (x == other.x) && (v == other.v) && (f == other.f) &&
			(type == other.type) & (m == other.m) && (old_f == other.old_f)) {
		return true;
	}

	return false;
}

#endif /* PARTICLE_H_ */
