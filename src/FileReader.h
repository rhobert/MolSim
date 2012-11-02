/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#ifndef FILE_READER_H_
#define FILE_READER_H_

#include "Particle.h"
#include <list>

/**
* @brief This class implements the functionality to create a list of particles out of a file
**/
class FileReader {

public:
	/**
	* @brief Create an instance of FileReader
	**/
	FileReader();
	
	/**
	* @brief Destructor of the class ParticleContainer
	**/
	virtual ~FileReader();
	
	/**
	* @brief Read a file in and creates a list of particles out of it
	* 
	* @param particles List where particles are written in
	* 
	* @param filename Name of file where the particles information are stored
	**/
	void readFile(std::list<Particle>& particles, char* filename);

};

#endif /* FILE_READER_H_ */
