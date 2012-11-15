/**
 * @file
 * @brief Simulation program
 */

#include "log4cxx/logger.h"
#include "log4cxx/propertyconfigurator.h"
#include "outputWriter/VTKWriter.h"
#include "FileReader.h"
#include "ParticleContainer.h"
#include "MaxwellBoltzmannDistribution.h"

#include <list>
#include <cstring>
#include <cstdlib>
#include <iostream>

using namespace std;

/**** forward declaration of the calculation functions ****/

#define ITERATION_PER_PLOT 10
#define EPSILON 5
#define SIGMA 1
#define BROWNIAN_MOTION 0.1

/**
 * @brief Calculate the force between two particles with the gravitational potential
 * 
 * @param p1 First particle
 * @param p2 Second particle
 * 
 * @return Force between particles
**/
utils::Vector<double, 3> gravitationalPotential(Particle& p1, Particle& p2);

/**
 * @brief Calculate the force between two particles with the Lenard-Jones potential
 * 
 * @param p1 First particle
 * @param p2 Second particle
 * 
 * @return Force between particles
**/
utils::Vector<double, 3> lenardJonesPotential(Particle& p1, Particle& p2);

/**
 * @brief Calculate the force for all particles
**/
void calculateF();

/**
 * @brief Calculate the position for all particles
**/
void calculateX();

/**
 * @brief Calculate the position for all particles
**/
void calculateV();

/**
 * @brief Plot the particles to a vtk-file
 * 
 * @param iteration Current iteration step 
**/
void plotParticles(int iteration);

/**
 * @brief Time step per iteration
**/
double delta_t;

/**
 * @brief Container with all particles and particle pairs
**/
ParticleContainer* particleContainer;

/**
 * @brief Function for force calculation
**/
utils::Vector<double, 3>  (*forceCalc)(Particle&, Particle&);

/**
 * @brief Program call syntax
**/
string molsim_usage = "\n"
	"Usage: ./MolSim END_T DELTA_T FILE_TYPE FILE POTENTIAL" "\n"
	"   END_T     - end time of simulation" "\n"
	"   DELTA_T   - timestep of simulation" "\n"
	"   FILE_TYPE - list or cuboid" "\n"
	"   FILE      - file with input data" "\n"
	"   POTENTIAL - gravitational or lenard_jones" "\n"
;

log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("MolSim"));


/**
 * @brief Simulation 
**/
int main(int argc, char* argsv[]) 
{
	log4cxx::PropertyConfigurator::configure("log4cxx.properties");
	
	LOG4CXX_INFO(logger, "Hello from MolSim for PSE!");
	
	// Check Parameters
	
	LOG4CXX_INFO(logger, "Check parameters");
	
	// Check count of parameters
	if ( (argc == 2 || argc == 3) && string(argsv[1]).compare("-test") == 0 )
	{
		LOG4CXX_INFO(logger, "Test option was passed");
		
		// Start all tests
		if ( argc == 2 )
		{
			LOG4CXX_INFO(logger, "Run all tests ...");
		}
		// Start sinlge test
		else
		{
			string test (argsv[2]);
			LOG4CXX_INFO(logger, "Run test " << test << " ...");
		}
		
		LOG4CXX_INFO(logger, "Finish tests");
		return EXIT_SUCCESS;
	}
	else if (argc != 6) 
	{
		LOG4CXX_FATAL(logger, "Errounous programme call: Wrong count of parameters!");
		cout << molsim_usage << endl;
		return EXIT_FAILURE;
	}
	
	LOG4CXX_INFO(logger, "Reading in parameters");
	
	// Init variables by parameters
	double end_time = atof(argsv[1]);
	delta_t = atof(argsv[2]);
	string file_type (argsv[3]);
	char* file_name = argsv[4];
	string potentialName (argsv[5]);
	
	LOG4CXX_INFO(logger, "Check END_T and DELTA_T");
	
	// Check end_time and delta_t
	if (end_time < 0 || delta_t < 0 || end_time < delta_t) 
	{
		LOG4CXX_FATAL(logger, "Errounous programme call: END_T and DELTA_T should be greater 0 and END_T greater DELTA_T!");
		cout << molsim_usage << endl;
		return EXIT_FAILURE;
	}
	
	LOG4CXX_INFO(logger, "END_T and DELTA_T are ok");
	
	// Read particles from file to Particle list and build ParticleContainer
	FileReader fileReader;
	list<Particle> particles;
	
	
	LOG4CXX_INFO(logger, "Check FILE_TYPE");
	
	// Check file_type
	if ( file_type.compare("list") == 0)
	{
		LOG4CXX_INFO(logger, "FILE_TYPE is set to list");
		fileReader.readFileList(particles, file_name);
		particleContainer = new ParticleContainer(particles);
	}
	else if ( file_type.compare("cuboid") == 0)
	{
		LOG4CXX_INFO(logger, "FILE_TYPE is set to cuboid");
		fileReader.readFileCuboid(particles, file_name);
		particleContainer = new ParticleContainer(particles);	
	}
	else
	{
		LOG4CXX_FATAL(logger, "Errounous programme call: For FILE_TYPE you have to specify list or cuboid!");
		cout << molsim_usage << endl;
		return EXIT_FAILURE;
	}
	
	LOG4CXX_INFO(logger, "Check POTENTIAL");
	
	// Check potentialName and set force calculation
	if ( potentialName.compare("gravitational") == 0)
	{
		LOG4CXX_INFO(logger, "force calculation for gravitational potential set");
		forceCalc = gravitationalPotential;
	}
	else if ( potentialName.compare("lenard_jones") == 0)
	{
		LOG4CXX_INFO(logger, "force calculation for Lenard-Jones potential set");
		forceCalc = lenardJonesPotential;
		
		LOG4CXX_INFO(logger, "Superpose velocity of particles with Brownian motion");
		// Superpose velocity of particles with Brownian motion
		for ( ParticleContainer::SingleList::iterator iterator = particleContainer->beginSingle();
			 iterator != particleContainer->endSingle();
			 iterator++ ) 
		{
			Particle& p = *iterator;
			MaxwellBoltzmannDistribution(p, BROWNIAN_MOTION, 3);
		}	
	}
	else
	{
		LOG4CXX_FATAL(logger, "Errounous programme call: For POTENTIAL you have to specify gravitational or lenard_jones");
		cout << molsim_usage << endl;
		return EXIT_FAILURE;
	}
	
	// Begin Simulation
	
	// the forces are needed to calculate x, but are not given in the input file.
	calculateF();
	
	LOG4CXX_INFO(logger, "Start simulation");
	
	// Init iteration variables
	double start_time = 0;
	double current_time;
	int iteration = 0;

	 // for this loop, we assume: current x, current f and current v are known
	for ( current_time = start_time; 
		 current_time < end_time; 
		 current_time += delta_t ) 
	{
		if (iteration % ITERATION_PER_PLOT == 0) {
			plotParticles(iteration);
		}
		
		// calculate new x
		calculateX();
		// calculate new f
		calculateF();
		// calculate new v
		calculateV();
		
		iteration++;
		LOG4CXX_INFO(logger, "Iteration " << iteration << " finished");
	}

	LOG4CXX_INFO(logger, "End simulation");
	
	return 0;
}

utils::Vector<double, 3> gravitationalPotential(Particle& p1, Particle& p2)
{
	utils::Vector<double, 3> x1_x2;
	utils::Vector<double, 3> F1_F2;
	
	//difference between coordinates of p1 and p2
	x1_x2 = p1.getX() - p2.getX();

	//force between p1 and p2
	F1_F2 = p1.getM() * p2.getM() / pow(x1_x2.L2Norm(),3) * (-1) * x1_x2;
	
	return F1_F2;
}

utils::Vector<double, 3> lenardJonesPotential(Particle& p1, Particle& p2)
{
	utils::Vector<double, 3> x1_x2;
	utils::Vector<double, 3> F1_F2;
	double l2Norm;
	double temp;
	
	//difference between coordinates of p1 and p2
	x1_x2 = p1.getX() - p2.getX();
	l2Norm = x1_x2.L2Norm();
	
	//force between p1 and p2
	temp = pow(SIGMA / l2Norm, 6);
	F1_F2 = 24*EPSILON / (l2Norm*l2Norm) * (temp - 2 * temp*temp) * (-1) * x1_x2;
	
	return F1_F2;
}

void calculateF() 
{
	utils::Vector<double, 3> x1_x2;
	utils::Vector<double, 3> F1_F2;
	utils::Vector<double, 3> F;
	
	F = 0.0;
	
	// Init forces effective on particles with 0.0 and save old force
	for ( ParticleContainer::SingleList::iterator iterator = particleContainer->beginSingle();
		 iterator != particleContainer->endSingle();
		 iterator++ ) 
	{
		Particle& p = *iterator;
		
		p.newF(F);
	}
	
	// Iterate over all particle pairs and calculate and sum effective forces
	for ( ParticleContainer::PairList::iterator iterator = particleContainer->beginPair();
		 iterator != particleContainer->endPair();
		 iterator++ ) 
	{		
		Particle& p1 = *(iterator->first);
		Particle& p2 = *(iterator->second);

		//force between p1 and p2
		F1_F2 = forceCalc(p1, p2);
		
		//sum forces
		F = p1.getF() + F1_F2;
		p1.setF(F);
		
		F = p2.getF() - F1_F2;
		p2.setF(F);
	}
}


void calculateX() 
{
	for ( ParticleContainer::SingleList::iterator iterator = particleContainer->beginSingle();
		 iterator != particleContainer->endSingle();
		 iterator++ ) 
	{
		Particle& p = *iterator;
		
		utils::Vector<double,3> x =
			p.getX() + 
			delta_t * p.getV() + 
			delta_t * delta_t / (2 * p.getM()) * p.getF();
			
		p.setX(x);
	}
}


void calculateV() 
{
	for ( ParticleContainer::SingleList::iterator iterator = particleContainer->beginSingle();
		 iterator != particleContainer->endSingle();
		 iterator++ ) 
	{
		Particle& p = *iterator;
		
		utils::Vector<double,3> v =
			p.getV() +
			delta_t / (2 * p.getM()) * (p.getF() + p.getOldF());
			
		p.setV(v);
	}
}

void plotParticles(int iteration) 
{
	LOG4CXX_INFO(logger, "Plot particles of Iteration " << iteration);
	string out_name("MD_vtk");
	outputWriter::VTKWriter vtk_writer;

	vtk_writer.initializeOutput(particleContainer->size());

	for ( ParticleContainer::SingleList::iterator iterator = particleContainer->beginSingle();
		 iterator != particleContainer->endSingle();
		 iterator++ ) 
	{
		vtk_writer.plotParticle(*iterator);	
	}

	vtk_writer.writeFile(out_name, iteration);
}