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
* @class FileReader
*
* @brief This class implements the functionality to initialize the particle-list.
*/

class FileReader {
	
	/**
	 * @brief Logger for FileReader class
	 */
	static log4cxx::LoggerPtr logger;
	
public:

	/**
	* @brief Create an instance of FileReader.
	*/
	FileReader();

	/**
	* @brief Destructor of the class FileReader.
	*/
	virtual ~FileReader();

	/**
	* @brief Initializes the particle-list using the input file.
	* 
	* @param particles List to write read particles in
	* 
	* @param filename File to read particles from
	*/
	void readFileList(std::list<Particle*>& particles, char* filename);
    
	/**
	* @brief Initializes the particle-list using the input file.
	* 
	* @param particles List to write read particles in
	* 
	* @param filename File to read cuboids from
	*/
	void readFileCuboid(std::list<Particle*>& particles, char* filename);

};

#endif /* FILE_READER_H_ */
