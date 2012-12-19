
#ifndef PHASESPACE_H_
#define PHASESPACE_H_

#include "Particle.h"
#include <list>
#include <iostream>

/**
* @class FileReader
*
* @brief This class implements the functionality to initialize the particle-list.
*/

class PhaseSpace {
	
	/**
	 * @brief Logger for PhaseSpace class
	 */
	static log4cxx::LoggerPtr logger;
	
public:

	/**
	* @brief Create an instance of PhaseSpace.
	*/
	PhaseSpace();

	/**
	* @brief Destructor of the class PhaseSpace.
	*/
	virtual ~PhaseSpace();

	/**
	* @brief Read in a particle-list using input stream
	* 
	* @param particles List to write read particles in
	* 
	* @param inStream Stream to read phase space from
	*/
	void readPhaseSpace(std::list<Particle>& particles, std::istream & inStream);
	
	/**
	* @brief Read in a particle-list using the input file.
	* 
	* @param particles List to write read particles in
	* 
	* @param filename File to read phase space from
	*/
	void readPhaseSpace(std::list<Particle>& particles, char* filename);
    
	/**
	* @brief Write out a particle-list using output stream
	* 
	* @param particles List to write read particles in
	* 
	* @param outStream Stream to write phase space in
	*/
	void writePhaseSpace(std::list<Particle>& particles, std::ostream & outStream);
	
	/**
	* @brief Write out a particle-list using the output file.
	* 
	* @param particles List to write read particles in
	* 
	* @param filename File to write phase space in
	*/
	void writePhaseSpace(std::list<Particle>& particles, char* filename);

};

#endif /* PHASESPACE_H_ */
