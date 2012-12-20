/*
 * XYZWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#ifndef XYZWRITER_H_
#define XYZWRITER_H_

#include "Particle.h"
#include <fstream>
#include <list>

namespace outputWriter {

/**
 * @class XYZWriter
 *
 * @brief This class implements the functionality to generate xyz output from particles.
 */

class XYZWriter {

public:
/**
* @brief Create an instance of XYZWriter.
*/
	XYZWriter();

/**
* @brief Destructor of the class VTKWriter.
*/
	virtual ~XYZWriter();

/**
* @brief plot type, mass, position, velocity and force of a particle.
*
* @param particles the list which contains the particles.
* @param filename the base name of the file to be written.
* @param iteration the number of the current iteration,
*        which is used to generate an unique filename
*/
	void plotParticles(std::list<Particle> particles, const std::string& filename, int iteration);

};

}

#endif /* XYZWRITER_H_ */