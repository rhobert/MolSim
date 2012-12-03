/**
 * @file
 * @brief Simulation program
 */

#include <list>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include <log4cxx/logger.h>
#include <log4cxx/propertyconfigurator.h>

#include "input/InputParameters.h"
#include "FileReader.h"
#include "outputWriter/VTKWriter.h"

#include "particleContainer/ParticleContainer.h"
#include "particleContainer/SimpleParticleContainer.h"

#include "ParticleGenerator.h"

#include "MaxwellBoltzmannDistribution.h"

#include "test/TestSettings.h"

using namespace std;

/**** forward declaration of the calculation functions ****/

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
 * @brief Apply Maxwell Boltzmann Distribution
 * 
 * @param p Particle to distribute
**/
void applyMaxwellBoltzmannDistribution( Particle& p );

/**
 * @brief Calculate and set the force between two particles
 * 
 * @param p1 First particle of the pair
 * @param p2 Second particle of the pair
**/
void calculateF( Particle& p1, Particle& p2 );

/**
 * @brief Set force effictive on a particle to 0
 * 
 * @param p Particle to modify
**/
void setNewForce( Particle& p );

/**
 * @brief Calculate the position of a particle
 * 
 * @param p Particle to modify
**/
void calculateX( Particle& p );

/**
 * @brief Calculate the position of a particle
 * 
 * @param p Particle to modify
**/
void calculateV( Particle& p );

/**
 * @brief Plot the particles to a vtk-file
 * 
 * @param iteration Current iteration step 
**/
void plotParticles(int iteration);

/**
 * @brief Plot particle to current VTKWriter vtk_writer
 * 
 * @param p Particle to modify
**/
void plotParticle( Particle& p );

/**
 * @brief Current output VTKWriter
**/
outputWriter::VTKWriter vtk_writer;

/**
 * @brief Time step per iteration
**/
double delta_t;

/**
 * @brief Destionation file to write output
**/
string outputFileName;

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
string molsim_usage = 
	"\n"
	"Usage: ./MolSim PARAMETER_FILE" "\n"
	"   PARAMETER_FILE - xml file with simulation parameters" "\n"
	"\n"
	"Usage: ./MolSim -test [TEST_NAME]" "\n"
	"   TEST_NAME      - run only TEST_NAME" "\n"
;

log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("MolSim"));


/**
 * @brief Simulation 
**/
int main(int argc, char* argsv[]) 
{	
	log4cxx::PropertyConfigurator::configure("log4cxx.properties");
	
	LOG4CXX_INFO(logger, "Hello from MolSim for PSE!");
	
	// Add available tests to usage info
	LOG4CXX_INFO(logger, "Detect available tests");
	
	TestSettings ts;
	vector<string> testNames = ts.getTestNames();
	
	molsim_usage += "Available tests:\n";
	
	for (int i; i < testNames.size(); i++)
	{	
		molsim_usage = molsim_usage + "   " + testNames[i] + "\n";
		LOG4CXX_TRACE(logger, "Detected " << testNames[i] << " as test");
	}
	
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
			
			ts.runTest();
		}
		// Start sinlge test
		else
		{
			string test_name (argsv[2]);
			LOG4CXX_INFO(logger, "Run test " << test_name << " ...");
			
			if ( ts.runTest(test_name) == 0 )
			{
				LOG4CXX_ERROR(logger, "Test " << test_name << " was not found");
				cout << molsim_usage << endl;
			}
		}
		
		
		LOG4CXX_INFO(logger, "Finish tests");
		return EXIT_SUCCESS;
	}
	else if (argc != 2) 
	{
		LOG4CXX_FATAL(logger, "Errounous programme call: Wrong count of parameters!");
		cout << molsim_usage << endl;
		return EXIT_FAILURE;
	}
	
	LOG4CXX_INFO(logger, "Reading parameter-file");
	
	std::ifstream file( argsv[1] );
	
	if ( !file.is_open() )
	{
		LOG4CXX_FATAL(logger, "Parameter-file can't be accessed!");
		cout << molsim_usage << endl;
		return EXIT_FAILURE;
	}
	
	LOG4CXX_INFO(logger, "Checking parameter-file");
	
	auto_ptr<PSE_Molekulardynamik_WS12::simulation_t> simulation;
	
	try
	{
		simulation = PSE_Molekulardynamik_WS12::simulation ( file );
	}
	catch (const xml_schema::exception& e)
	{
		LOG4CXX_FATAL(logger, "Parameter-file is not valid!");
		cout << e << endl;
		cout << molsim_usage << endl;
		return EXIT_FAILURE;
	}
	
	LOG4CXX_INFO(logger, "Reading in parameters");
	
	// Init variables by parameters
	double end_time = simulation->t_end();
	delta_t = simulation->delta_t();
	int writeFrequency = simulation->writeFrequency();
	outputFileName = simulation->outputFile();
	string file_name;
	
	LOG4CXX_INFO(logger, "Check END_T and DELTA_T");
	
	// Check end_time and delta_t
	if (end_time < delta_t) 
	{
		LOG4CXX_FATAL(logger, "Errounous programme call: END_T should be greater than DELTA_T!");
		cout << molsim_usage << endl;
		return EXIT_FAILURE;
	}
	
	// Read particles from file to Particle list and build ParticleContainer
	FileReader fileReader;
	list<Particle> particles;
	
	LOG4CXX_INFO(logger, "Reading in input files");
	
	for ( PSE_Molekulardynamik_WS12::inputFiles_t::inputFile_const_iterator i = simulation->inputFiles().inputFile().begin(); 
		 i != simulation->inputFiles().inputFile().end(); 
		 i++ )
	{
		file_name = *i;
		
		LOG4CXX_INFO(logger, "Reading in " << file_name << " with type " << i->type());
		
		// Check file_type
		if ( i->type() == PSE_Molekulardynamik_WS12::inputType_t::list )
		{
			fileReader.readFileList(particles, (char*) file_name.c_str());
			particleContainer = new SimpleParticleContainer(particles);
		}
		else if ( i->type() == PSE_Molekulardynamik_WS12::inputType_t::cuboid )
		{
			fileReader.readFileCuboid(particles, (char*) file_name.c_str());
			particleContainer = new SimpleParticleContainer(particles);	
		}
	}
	
	LOG4CXX_INFO(logger, "Set POTENTIAL");
	
	// Check potentialName and set force calculation
	if ( simulation->potential() == PSE_Molekulardynamik_WS12::potential_t::gravitational )
	{
		LOG4CXX_INFO(logger, "force calculation for gravitational potential set");
		forceCalc = gravitationalPotential;
	}
	else if ( simulation->potential() == PSE_Molekulardynamik_WS12::potential_t::lenard_jones )
	{
		LOG4CXX_INFO(logger, "force calculation for Lenard-Jones potential set");
		forceCalc = lenardJonesPotential;
		
		LOG4CXX_INFO(logger, "Superpose velocity of particles with Brownian motion");
		
		// Superpose velocity of particles with Brownian motion
		particleContainer->applyToSingleParticles ( applyMaxwellBoltzmannDistribution );
	}
	
	// Begin Simulation
	
	// the forces are needed to calculate x, but are not given in the input file.
	particleContainer->applyToSingleParticles ( setNewForce );
	particleContainer->applyToParticlePairs ( calculateF );
	
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
		if (iteration % writeFrequency == 0) {
			plotParticles(iteration);
		}
		
		// calculate new x
		particleContainer->applyToSingleParticles ( calculateX );
		// calculate new f
		particleContainer->applyToSingleParticles ( setNewForce );
		particleContainer->applyToParticlePairs ( calculateF );
		// calculate new v
		particleContainer->applyToSingleParticles ( calculateV );
		
		iteration++;
		LOG4CXX_TRACE(logger, "Iteration " << iteration << " finished");
	}

	LOG4CXX_INFO(logger, "End simulation");
	
	return 0;
}

void applyMaxwellBoltzmannDistribution( Particle& p )
{
	MaxwellBoltzmannDistribution(p, BROWNIAN_MOTION, 3);
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

void setNewForce( Particle& p )
{
	p.setF( utils::Vector<double,3>(0.0) );
}

void calculateF( Particle& p1, Particle& p2 ) 
{
	utils::Vector<double, 3> x1_x2;
	utils::Vector<double, 3> F1_F2;
	utils::Vector<double, 3> F;

	//force between p1 and p2
	F1_F2 = forceCalc(p1, p2);
	
	//sum forces
	F = p1.getF() + F1_F2;
	p1.setF(F);
	
	F = p2.getF() - F1_F2;
	p2.setF(F);
}


void calculateX( Particle& p ) 
{
	utils::Vector<double,3> x =
		p.getX() + 
		delta_t * p.getV() + 
		delta_t * delta_t / (2 * p.getM()) * p.getF();
		
	p.setX(x);
}


void calculateV( Particle& p ) 
{
	utils::Vector<double,3> v =
		p.getV() +
		delta_t / (2 * p.getM()) * (p.getF() + p.getOldF());
		
	p.setV(v);
}

void plotParticle( Particle& p )
{
	vtk_writer.plotParticle(p);
}

void plotParticles(int iteration) 
{
	LOG4CXX_INFO(logger, "Plot particles of Iteration " << iteration);
	
	vtk_writer.initializeOutput(particleContainer->size());
	
	particleContainer->applyToSingleParticles( plotParticle );

	vtk_writer.writeFile(outputFileName, iteration);
}