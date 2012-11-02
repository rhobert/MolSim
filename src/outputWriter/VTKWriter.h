/*
 * VTKWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#ifndef VTKWRITER_H_
#define VTKWRITER_H_

#include "Particle.h"
#include "outputWriter/vtk-unstructured.h"

#include <list>

namespace outputWriter {

/**
 * @class VTKWriter
 *
 * @brief This class implements the functionality to generate vtk output from particles.
 */
class VTKWriter {

public:
/**
* @brief Create an instance of VTKWriter.
*/
	VTKWriter();

/**
* @brief Destructor of the class VTKWriter.
*/
    virtual ~VTKWriter();

	/**
	 * @brief set up internal data structures and prepare to plot a particle.
	 */
	void initializeOutput(int numParticles);

	/**
	 * @brief plot type, mass, position, velocity and force of a particle.
	 *
	 * @note: initializeOutput() must have been called before.
	 */
	void plotParticle(Particle& p);

	/**
	 * @brief writes the final output file.
	 *
	 * @param filename the base name of the file to be written.
	 * @param iteration the number of the current iteration,
	 *        which is used to generate an unique filename
	 */
	void writeFile(const std::string& filename, int iteration);

private:
	VTKFile_t* vtkFile;
};

}

#endif /* VTKWRITER_H_ */
