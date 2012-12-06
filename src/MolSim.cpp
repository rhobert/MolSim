/** 
 * @brief Molecular dynamics simulation program
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
#include "particleContainer/LinkedCellParticleContainer.h"

#include "ParticleGenerator.h"

#include "MaxwellBoltzmannDistribution.h"

#include "test/TestSettings.h"

using namespace std;

/**** forward declaration of the calculation functions ****/

#define EPSILON 5
#define SIGMA 1

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
 * @brief Calcuation for the reflection boundary condition for lenard jones potential
 * 
 * @param p Particle to reflect
**/
void calcReflection (Particle& p);

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
 * @brief LinkedCellParticleContainer with all particles and particle pairs if exist
**/
LinkedCellParticleContainer* linkedCellParticleContainer = NULL;

/**
 * @brief Function for force calculation
**/
utils::Vector<double, 3>  (*forceCalc)(Particle&, Particle&);

/**
 * @brief Count of dimensions where Brownian Motion is effictive
**/
int dimensionsBrownianMotion = 3;

/**
 * @brief Count of dimensions where Brownian Motion is effictive
**/
utils::Vector<double,3> domainSize (0.0);

/**
 * @brief Mean velocity for Brownian Motion
**/
double meanVelocityBrownianMotion = 0.1;

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
	
	// Add available tests for usage info
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
	
	LOG4CXX_TRACE(logger, "Check parameters");
	
	// Check count of parameters
	if ( (argc == 2 || argc == 3) && string(argsv[1]).compare("-test") == 0 )
	{
		LOG4CXX_DEBUG(logger, "Test option was passed");
		
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
				LOG4CXX_FATAL(logger, "Test " << test_name << " was not found");
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
	
	LOG4CXX_DEBUG(logger, "Checking parameter-file");
	
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
	
	LOG4CXX_DEBUG(logger, "Reading in parameters");
	
	// Init variables by parameters
	double end_time = simulation->t_end();
	delta_t = simulation->delta_t();
	int writeFrequency = simulation->writeFrequency();
	outputFileName = simulation->outputFile();
	
	string file_name;
	PSE_Molekulardynamik_WS12::boundary_t boundary (PSE_Molekulardynamik_WS12::boundary_t::outflow);
	
	PSE_Molekulardynamik_WS12::domain_t* domain = NULL;
	double cutoff = 0;
	
	LOG4CXX_DEBUG(logger, "Check END_T and DELTA_T");
	
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
	
	LOG4CXX_DEBUG(logger, "Reading in input files");
	
	for ( PSE_Molekulardynamik_WS12::inputs_t::inputFile_const_iterator i = simulation->inputs().inputFile().begin(); 
		 i != simulation->inputs().inputFile().end(); 
		 i++ )
	{
		file_name = *i;
		
		LOG4CXX_INFO(logger, "Reading in " << file_name << " with type " << i->type());
		
		// Check file_type
		if ( i->type() == PSE_Molekulardynamik_WS12::inputType_t::list )
		{
			fileReader.readFileList(particles, (char*) file_name.c_str());
		}
		else if ( i->type() == PSE_Molekulardynamik_WS12::inputType_t::cuboid )
		{
			fileReader.readFileCuboid(particles, (char*) file_name.c_str());
		}
	}
	
	LOG4CXX_DEBUG(logger, "Reading in cuboids");
	
	for ( PSE_Molekulardynamik_WS12::inputs_t::cuboid_const_iterator i = simulation->inputs().cuboid().begin(); 
		 i != simulation->inputs().cuboid().end(); 
		 i++ )
	{
		PSE_Molekulardynamik_WS12::cuboid_t cuboid = *i;
		
		utils::Vector<double,3> position;
		utils::Vector<double,3> velocity;
		utils::Vector<int,3> dimensions;
		double mass;
		double distance;
		
		position[0] = cuboid.position().x();
		position[1] = cuboid.position().y();
		position[2] = cuboid.position().z();
		
		velocity[0] = cuboid.velocity().x();
		velocity[1] = cuboid.velocity().y();
		velocity[2] = cuboid.velocity().z();
		
		dimensions[0] = cuboid.dimensions().x();
		dimensions[1] = cuboid.dimensions().y();
		dimensions[2] = cuboid.dimensions().z();
		
		mass = cuboid.mass();
		distance = cuboid.distance();
		
		LOG4CXX_INFO(logger, "Reading in cuboid at " << position.toString() << " with velocity " << velocity.toString() 
			<< ", dimensions " << dimensions.toString() << ", mass " << mass << " and distance " << distance);
		
		generateCuboid(particles, position, velocity, dimensions, distance, mass);
	}
	
	LOG4CXX_DEBUG(logger, "Reading in spheres");
	
	for ( PSE_Molekulardynamik_WS12::inputs_t::sphere_const_iterator i = simulation->inputs().sphere().begin(); 
		 i != simulation->inputs().sphere().end(); 
		 i++ )
	{
		PSE_Molekulardynamik_WS12::sphere_t sphere = *i;
		
		utils::Vector<double,3> position;
		utils::Vector<double,3> velocity;
		int radiusDimension;
		int dimensionCount;
		double mass;
		double distance;
		
		position[0] = sphere.position().x();
		position[1] = sphere.position().y();
		position[2] = sphere.position().z();
		
		velocity[0] = sphere.velocity().x();
		velocity[1] = sphere.velocity().y();
		velocity[2] = sphere.velocity().z();
		
		radiusDimension = sphere.radiusDimension();
		dimensionCount = sphere.dimensionCount();
		mass = sphere.mass();
		distance = sphere.distance();	
		
		LOG4CXX_INFO(logger, "Reading in sphere at " << position.toString() << " with velocity " << velocity.toString() 
			<< ", radius dimension " << radiusDimension << ", mass " << mass << " and distance " << distance << " in " << dimensionCount << " dimensions" );
		
		generateSphere(particles, position, velocity, radiusDimension, dimensionCount, distance, mass);
	}
	
	// Read in domain, boundary condition and cutoff radius for LinkedCellParticleContainer 
	 
	LOG4CXX_DEBUG(logger, "Detect if domain size is specified");
	
	if ( simulation->domain().present() )
	{
		domain = &simulation->domain().get();
		
		domainSize[0] = domain->dimensions().x();
		domainSize[1] = domain->dimensions().y();
		domainSize[2] = domain->dimensions().z();
		cutoff = domain->cutoff();
		boundary = domain->boundary();
		
		LOG4CXX_INFO(logger, "Domain " << domainSize.toString() << " with cutoff-radius " << cutoff << " and " << boundary << " boundary condition");
		
		LOG4CXX_DEBUG(logger, "Create LinkedCellParticleContainer");
		
		linkedCellParticleContainer =  new LinkedCellParticleContainer(domainSize, cutoff );
		particleContainer = dynamic_cast<ParticleContainer*>( linkedCellParticleContainer );
		
		domainSize = linkedCellParticleContainer->getDomainSize();
		
		LOG4CXX_DEBUG(logger, "Real domain size is " << domainSize.toString() );
	}
	else
	{
		LOG4CXX_DEBUG(logger, "No domain and boundary conditions specified");
		LOG4CXX_DEBUG(logger, "Create SimpleParticleContainer");
		
		particleContainer = new SimpleParticleContainer();	
	}
	
	LOG4CXX_INFO(logger, "Add " << particles.size() << " particles to ParticleContainer");
	
	particleContainer->addParticles( particles );
	
	LOG4CXX_INFO(logger, "Set potential for force calulation to " << simulation->potential() );
	
	// Check potentialName and set force calculation
	if ( simulation->potential() == PSE_Molekulardynamik_WS12::potential_t::gravitational )
	{
		forceCalc = gravitationalPotential;
	}
	else if ( simulation->potential() == PSE_Molekulardynamik_WS12::potential_t::lenard_jones )
	{
		forceCalc = lenardJonesPotential;
		
		LOG4CXX_DEBUG(logger, "Detect if parameters for Brownian Motion are specified");
		
		if ( simulation->brownianMotion().present() )
		{
			meanVelocityBrownianMotion = simulation->brownianMotion().get().meanVelocity();
			dimensionsBrownianMotion = simulation->brownianMotion().get().dimensionCount();
			
			LOG4CXX_DEBUG(logger, "Parameters for Brownian Motion are specified");
		}
		else
		{
			LOG4CXX_DEBUG(logger, "No parameters for Brownian Motion are specified");
		}
		
		LOG4CXX_INFO(logger, "Superpose velocity of particles with Brownian motion with mean velocity " << meanVelocityBrownianMotion << " in " << dimensionsBrownianMotion << " dimensions");
		
		// Superpose velocity of particles with Brownian motion
		particleContainer->applyToSingleParticles ( applyMaxwellBoltzmannDistribution );
	}
	
	// Begin Simulation
	
	// the forces are needed to calculate x, but are not given in the input file.
	particleContainer->applyToSingleParticles ( setNewForce );
	particleContainer->applyToParticlePairs ( calculateF );
	
	LOG4CXX_INFO(logger, "Write vtk output to files with basename " << outputFileName << " with write frequency " << writeFrequency );
	
	// Init iteration variables
	double start_time = 0;
	double current_time;
	int iteration = 0;
	
	LOG4CXX_INFO(logger, "Start simulation from time " << start_time << " to " << end_time << " with time steps " << delta_t);
	LOG4CXX_INFO(logger, "There will be " << floor((end_time - start_time) / delta_t) << " simulation steps and " << ceil(floor((end_time - start_time) / delta_t) / writeFrequency) << " output files" );

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
		
		// Set up new force (Could be modified by boundary condition)
		particleContainer->applyToSingleParticles ( setNewForce );
		
		if ( linkedCellParticleContainer != NULL )
		{
			linkedCellParticleContainer->updateContainingCells();
			
			if (boundary == PSE_Molekulardynamik_WS12::boundary_t::outflow)
			{
				linkedCellParticleContainer->deleteHaloParticles();
			}
			else if (boundary == PSE_Molekulardynamik_WS12::boundary_t::reflecting)
			{
				linkedCellParticleContainer->applyToBoundaryParticles(calcReflection);
			}
		}
		
		// calculate new f
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
	MaxwellBoltzmannDistribution(p, meanVelocityBrownianMotion, dimensionsBrownianMotion);
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

void calcReflection (Particle& p)
{
	utils::Vector<double,3> x;
	utils::Vector<double,3> x1_x2;
	double d = pow(2,1/6) * SIGMA;
	
	
	LOG4CXX_TRACE(logger, "Reflection is checked for " << p.getX().toString() )
	
	for ( int i = 0; i < 3; i++ )
	{
		if ( domainSize[i] != 0 )
		{
			LOG4CXX_TRACE(logger, "Reflection is for " << i << " dimension" );
			
			x = p.getX();
			
			x[i] = 0;
			x1_x2 = p.getX() - x;
			
			if ( x1_x2.L2Norm() <= d  )
			{
				LOG4CXX_TRACE(logger, "CounterParticle " << x.toString() );
				Particle counterParticle ( x, utils::Vector<double,3>(0.0), p.getM() );
				calculateF( p, counterParticle );
			}
			
			x[i] = domainSize[i];
			x1_x2 = p.getX() - x;
			
			if ( x1_x2.L2Norm() <= d  )
			{
				LOG4CXX_TRACE(logger, "CounterParticle " << x.toString()  );
				Particle counterParticle ( x, utils::Vector<double,3>(0.0), p.getM() );
				calculateF( p, counterParticle );
			}
		}
	}
}

void setNewForce( Particle& p )
{
	p.newF( utils::Vector<double,3>(0.0) );
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