/**
 * @brief Molecular dynamics simulation program
 */

#include <list>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sys/time.h>

#include <log4cxx/logger.h>
#include <log4cxx/propertyconfigurator.h>

#include "input/InputParameters.h"
#include "PhaseSpace.h"
#include "FileReader.h"
#include "outputWriter/VTKWriter.h"

#include "particleContainer/ParticleContainer.h"
#include "particleContainer/SimpleParticleContainer.h"
#include "particleContainer/LinkedCellParticleContainer.h"

#include "ParticleGenerator.h"

#include "MaxwellBoltzmannDistribution.h"

#include "Thermostat.h"

#include "test/TestSettings.h"

using namespace std;

/**** forward declaration of the calculation functions ****/


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
 * @brief Calcuation for the periodic boundary condition (move halo particles to their destinied boundary cells)
 *
 * @param p Particle to move
**/
void calcPeriodicHalo(Particle &p);

/**
 * @brief Calcuation for the periodic boundary condition
 *
 * @param p1 Boundary particle to modfiy
 * @param p2 Boundary particle of the opposite boundary
**/
void calcPeriodicBoundary(Particle &p1, Particle p2);


/**
 * @brief Apply gravitation on partcile
 *
 * @param p Particle to modify
**/
void applyGravitation(Particle& p);

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
utils::Vector<double,3> domainSize (0.0);

/**
 * @brief Mean velocity for Maxwell-Boltzmann distribution
**/
utils::Vector<double,3> gravitation (0.0);

/**
 * @brief Mean velocity for Maxwell-Boltzmann distribution
**/
double meanVelocity;

/**
 * @brief Count of dimensions where the simulation takes place
**/
int dimensionCount;

/**
 * @brief Cutoff-radius for force calulation
**/
double cutoff = 0;

/**
 * @brief Thermostat for temperature regulation
**/
Thermostat* thermostat;

/**
 * @brief The number of timesteps after which the thermostat is applied
**/
int nThermostat;

/**
 * @brief Controls use of thermostat
**/
bool thermostatOn = false;

/**
 * @brief Boundary conditions for all boundaries
**/
PSE_Molekulardynamik_WS12::boundary_t boundary []  = {
		PSE_Molekulardynamik_WS12::boundary_t::outflow, PSE_Molekulardynamik_WS12::boundary_t::outflow, PSE_Molekulardynamik_WS12::boundary_t::outflow, 
		PSE_Molekulardynamik_WS12::boundary_t::outflow, PSE_Molekulardynamik_WS12::boundary_t::outflow, PSE_Molekulardynamik_WS12::boundary_t::outflow};

/**
 * @brief Program call syntax
**/
string molsim_usage =
	"\n" "Usage:" "\n"
	"    ./MolSim help" "\n"
	"          Prints this usage information on standard output." "\n" "\n"
	"    ./MolSim run PARAMETER_FILE" "\n"
	"          Start simulation with parameters specified in PARAMETER_FILE." "\n"
	"          File is in xml-format which is specified in ./src/input/InputParameters.xsd (description in input.xml too)" "\n" "\n"
	"    ./MolSim test [TEST_NAME]" "\n"
	"          Runs all tests or only TEST_NAME:" "\n"
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

	for (int i = 0; i < testNames.size(); i++)
	{
		molsim_usage = molsim_usage + "          - " + testNames[i] + "\n";
		LOG4CXX_TRACE(logger, "Detected " << testNames[i] << " as test");
	}

	// Check Parameters

	LOG4CXX_TRACE(logger, "Check parameters");

	// Check count of parameters
	if ( (argc == 2 || argc == 3) && (string(argsv[1]).compare("test") == 0 || string(argsv[1]).compare("-test") == 0 || string(argsv[1]).compare("--test") == 0) )
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
	else if ( argc == 2 && (string(argsv[1]).compare("help") == 0 || string(argsv[1]).compare("-help") == 0 || string(argsv[1]).compare("--help") == 0) )
	{
		LOG4CXX_INFO(logger, "Usage information ist printed");
		cout << molsim_usage << endl;
		return EXIT_FAILURE;
	}
	else if ( argc != 3 || (string(argsv[1]).compare("run") != 0 && string(argsv[1]).compare("-run") != 0 && string(argsv[1]).compare("--run") != 0) )
	{
		LOG4CXX_FATAL(logger, "Errounous programme call: Wrong parameters!");
		cout << molsim_usage << endl;
		return EXIT_FAILURE;
	}

	LOG4CXX_INFO(logger, "Reading parameter-file " << argsv[2] );

	std::ifstream file( argsv[2] );

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
	dimensionCount = simulation->dimensionCount();
	
	PhaseSpace phaseSpace;
	int writeFrequency = 0;
	string file_name = "";

	PSE_Molekulardynamik_WS12::domain_t* domain = NULL;

	LOG4CXX_DEBUG(logger, "Check END_T and DELTA_T");

	// Check end_time and delta_t
	if (end_time < delta_t)
	{
		LOG4CXX_FATAL(logger, "Errounous programme call: END_T should be greater than DELTA_T!");
		cout << molsim_usage << endl;
		return EXIT_FAILURE;
	}
	
	LOG4CXX_DEBUG(logger, "Check if output is specified");
	
	if ( simulation->output().present() )
	{
		outputFileName = simulation->output().get().file();
		writeFrequency = simulation->output().get().writeFrequency();
		
		LOG4CXX_INFO(logger, "Write vtk output to files with basename " << outputFileName << " with write frequency " << writeFrequency );
	}
	else
	{
		LOG4CXX_WARN(logger, "No output specified" );
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
			fileReader.readFileList( particles, (char*) file_name.c_str() );
		}
		else if ( i->type() == PSE_Molekulardynamik_WS12::inputType_t::cuboid )
		{
			fileReader.readFileCuboid( particles, (char*) file_name.c_str() );
		}
		else if ( i->type() == PSE_Molekulardynamik_WS12::inputType_t::phaseSpace )
		{
			phaseSpace.readPhaseSpace( particles, (char*) file_name.c_str() );
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
		double sigma;
		double epsilon;
		int type;

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

		sigma = cuboid.sigma();
		epsilon = cuboid.epsilon();
		type = cuboid.type();

		LOG4CXX_INFO(logger, "Reading in cuboid at " << position.toString() << " with velocity " << velocity.toString()
			<< ", dimensions " << dimensions.toString() << ", mass " << mass << ", distance " << distance << ", sigma " << sigma
			<< ", epsilon " << epsilon << " and type " << type);

		generateCuboid(particles, position, velocity, dimensions, distance, mass, sigma, epsilon, type);
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
		double mass;
		double distance;
		double sigma;
		double epsilon;
		int type;

		position[0] = sphere.position().x();
		position[1] = sphere.position().y();
		position[2] = sphere.position().z();

		velocity[0] = sphere.velocity().x();
		velocity[1] = sphere.velocity().y();
		velocity[2] = sphere.velocity().z();

		radiusDimension = sphere.radiusDimension();
		mass = sphere.mass();
		distance = sphere.distance();

		sigma = sphere.sigma();
		epsilon = sphere.epsilon();
		type = sphere.type();

		LOG4CXX_INFO(logger, "Reading in sphere at " << position.toString() << " with velocity " << velocity.toString()
			<< ", radius dimension " << radiusDimension << ", mass " << mass << ", distance " << distance << ", sigma " << sigma
			<< ", epsilon " << epsilon << " and type " << type << " in " << dimensionCount << " dimensions" );

		generateSphere(particles, position, velocity, radiusDimension, dimensionCount, distance, mass, sigma, epsilon, type);
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
		
		LOG4CXX_INFO(logger, "Domain " << domainSize.toString() << " with cutoff-radius " << cutoff);
		
		boundary[0] = domain->x().lower();
		boundary[1] = domain->x().upper();
		
		LOG4CXX_INFO(logger, "In x-dimension the lower boundary condition is " << boundary[0] << " and the upper " << boundary[1]);
		
		if ( dimensionCount > 1 )
		{
			if ( domain->y().present() )
			{
				boundary[2] = domain->y().get().lower();
				boundary[3] = domain->y().get().upper();
				
				LOG4CXX_INFO(logger, "In y-dimension the lower boundary condition is " << boundary[2] << " and the upper " << boundary[3]);
				
				if ( dimensionCount > 2 )
				{
					if ( domain->z().present() )
					{
						boundary[4] = domain->z().get().lower();
						boundary[5] = domain->z().get().upper();
						
						LOG4CXX_INFO(logger, "In z-dimension the lower boundary condition is " << boundary[4] << " and the upper " << boundary[5]);
					}
					else
					{
						LOG4CXX_ERROR(logger, "Dimension count is set to " << dimensionCount << ", but no boundary condition is set for z-dimension");
						return EXIT_FAILURE;
					}
				}
			}
			else
			{
				LOG4CXX_ERROR(logger, "Dimension count is set to " << dimensionCount << ", but no boundary condition is set for y-dimension");
				return EXIT_FAILURE;
			}
		}

		LOG4CXX_DEBUG(logger, "Create LinkedCellParticleContainer");

		linkedCellParticleContainer =  new LinkedCellParticleContainer( domainSize, cutoff );
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
	
	LOG4CXX_DEBUG(logger, "Detect if gravitation is specified");
	
	if ( simulation->gravitation().present() )
	{
		gravitation[0] = simulation->gravitation().get().x();
		gravitation[1] = simulation->gravitation().get().y();
		gravitation[2] = simulation->gravitation().get().z();
		
		LOG4CXX_INFO(logger, "Gravitation is set to " << gravitation.toString());
	}
	
	LOG4CXX_INFO(logger, "Set potential for force calulation to " << simulation->potential() );

	// Check potentialName and set force calculation
	if ( simulation->potential() == PSE_Molekulardynamik_WS12::potential_t::gravitational )
	{
		forceCalc = gravitationalPotential;
	}
	else if ( simulation->potential() == PSE_Molekulardynamik_WS12::potential_t::lenard_jones )
	{
		forceCalc = lenardJonesPotential;

		LOG4CXX_DEBUG(logger, "Detect if parameters for Thermostat are specified");

		if ( simulation->thermostat().present() )
		{
			double initialT = simulation->thermostat().get().initialT();
			
//			double targetT = simulation->thermostat().get().targetT();
//			double diffT = simulation->thermostat().get().diffT();
//			int nMax = simulation->thermostat().get().nMax();
//			nThermostat = simulation->thermostat().get().nThermostat();

			thermostatOn = true;

			LOG4CXX_DEBUG(logger, "Parameters for Thermostat are specified");

			thermostat = new Thermostat(*particleContainer, initialT, dimensionCount);
			
			if ( simulation->thermostat().get().frequency().present() )
			{
				int regulationFrequency = simulation->thermostat().get().frequency().get();
				
				thermostat->setFrequency( regulationFrequency );
				LOG4CXX_DEBUG(logger, "Regulation frequency for Thermostat is " << regulationFrequency );
			}
			
			LOG4CXX_INFO(logger, "Superpose velocity of particles to get initial temperature with temperature " << initialT << " in " << dimensionCount << " dimensions");
		}
		else
		{
			LOG4CXX_DEBUG(logger, "No parameters for Brownian Motion or Thermostat are specified");
		}
	}

	// Begin Simulation

	// the forces are needed to calculate x, but are not given in the input file.
	particleContainer->applyToParticlePairs ( calculateF );

	// Init iteration variables
	double start_time = 0;
	double current_time;
	int iteration = 0;
	
	double totalTime = 0;
	double timeDiff;
	timeval t0, t1; 
	
	
	LOG4CXX_INFO(logger, "Start simulation from time " << start_time << " to " << end_time << " with time steps " << delta_t);
	LOG4CXX_INFO(logger, "Expecting " << floor((end_time - start_time) / delta_t) << " simulation steps" );
	if ( writeFrequency > 0 )
		LOG4CXX_INFO(logger, "Expecting " << ceil(floor((end_time - start_time) / delta_t) / writeFrequency) << " output files" );
	
	 // for this loop, we assume: current x, current f and current v are known
	for ( current_time = start_time;
		 current_time < end_time;
		 current_time += delta_t )
	{
		gettimeofday(&t0, NULL);
		
		// apply thermostat
		if ( thermostatOn == true )
		{
			thermostat->apply( iteration );
		}

		
		if (writeFrequency != 0 && iteration % writeFrequency == 0) {
			plotParticles(iteration);
		}
		
		// calculate new x
		particleContainer->applyToSingleParticles ( calculateX );

		// Set up new force (Could be modified by boundary condition)
		particleContainer->applyToSingleParticles ( setNewForce );
		
		if ( linkedCellParticleContainer != NULL )
		{	
			linkedCellParticleContainer->updateContainingCells();
			
			for ( int i = 0; i < dimensionCount * 2; i++ )
			{					
				if (boundary[i] == PSE_Molekulardynamik_WS12::boundary_t::outflow)
				{
					linkedCellParticleContainer->deleteHaloParticles(i);
				}
				else if (boundary[i] == PSE_Molekulardynamik_WS12::boundary_t::reflecting)
				{
					linkedCellParticleContainer->applyToBoundaryParticles(i, calcReflection);
				}
				else if (boundary[i] == PSE_Molekulardynamik_WS12::boundary_t::periodic)
				{
					linkedCellParticleContainer->applyToPeriodicBoundaryParticlePairs(i, calcPeriodicBoundary);
					linkedCellParticleContainer->applyToHaloParticles(i, calcPeriodicHalo);
					linkedCellParticleContainer->updateContainingCells();
				}
			}
		}
		
		// calculate new f
		particleContainer->applyToParticlePairs ( calculateF );
		
		// apply gravitation
		particleContainer->applyToSingleParticles ( applyGravitation );
		
		// calculate new v
		particleContainer->applyToSingleParticles ( calculateV );
		
		
		gettimeofday(&t1, NULL);
		timeDiff = (t1.tv_sec - t0.tv_sec) + (t1.tv_usec - t0.tv_usec) * 1e-9;
		totalTime += timeDiff;
		
		iteration++;
		LOG4CXX_TRACE(logger, "Iteration " << iteration << " finished in " << timeDiff << " seconds");
	}

	LOG4CXX_INFO(logger, "End simulation in " << totalTime << " seconds" );
	
	if ( iteration > 0 )
		LOG4CXX_INFO(logger, (totalTime / (double) iteration) << " per iteration");
	
	if ( simulation->outputPhaseSpace().present() )
	{
		string phaseSpaceFilename = simulation->outputPhaseSpace().get();
		
		LOG4CXX_INFO(logger, "Write phase space to " << phaseSpaceFilename );
		
		list<Particle> pList = particleContainer->getParticles();
		phaseSpace.writePhaseSpace( *&pList, const_cast<char*> ( phaseSpaceFilename.c_str() ) );
	}

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
	double sigma;
	double epsilon;

	if ( p1.getType() == p2.getType() )
	{
		sigma = p1.getSigma();
		epsilon = p1.getEpsilon();
	}
	else
	{
		sigma = ( p1.getSigma() + p2.getSigma() ) / 2;
		epsilon = sqrt( p1.getEpsilon() * p2.getEpsilon() );
	}

	//difference between coordinates of p1 and p2
	x1_x2 = p1.getX() - p2.getX();
	l2Norm = x1_x2.L2Norm();

	//force between p1 and p2
	temp = pow(sigma / l2Norm, 6);
	F1_F2 = 24*epsilon / (l2Norm*l2Norm) * (temp - 2 * temp*temp) * (-1) * x1_x2;

	return F1_F2;
}

void calcReflection (Particle& p)
{
	utils::Vector<double,3> x;
	utils::Vector<double,3> x1_x2;
	double d = pow(2,1/6) * p.getSigma();

	LOG4CXX_TRACE(logger, "Reflection is checked for " << p.getX().toString() );

	for ( int i = 0; i < 3; i++ )
	{
		if ( domainSize[i] != 0 )
		{
			LOG4CXX_TRACE(logger, "Reflection is set for " << i << " dimension" );

			x = p.getX();
			
			if ( boundary[i*2] == PSE_Molekulardynamik_WS12::boundary_t::reflecting )
			{
				x[i] = 0;
				x1_x2 = p.getX() - x;
				
				if ( x1_x2.L2Norm() <= d  )
				{
					LOG4CXX_TRACE(logger, "CounterParticle " << x.toString() );
					Particle counterParticle ( x, utils::Vector<double,3>(0.0), p.getM() );
					calculateF( p, counterParticle );
				}
			}
			
			if ( boundary[i*2+1] == PSE_Molekulardynamik_WS12::boundary_t::reflecting )
			{
			
				x[i] = domainSize[i];
				x1_x2 = p.getX() - x;

				if ( x1_x2.L2Norm() <= d  )
				{
					LOG4CXX_TRACE(logger, "CounterParticle " << x.toString() );
					Particle counterParticle ( x, utils::Vector<double,3>(0.0), p.getM() );
					calculateF( p, counterParticle );	
				}
			}
		}
	}
}

void calcPeriodicHalo(Particle &p)
{
	utils::Vector<double,3> h;
	
	LOG4CXX_TRACE(logger, "Move particle from " << p.getX().toString() );
	
	if(p.getX()[0] < 0){
		h[0] = domainSize[0];
		h[1] = 0;
		h[2] = 0;
		p.setX(p.getX()+h);
	}
	else if(p.getX()[0] > domainSize[0]){
		h[0] = domainSize[0];
		h[1] = 0;
		h[2] = 0;
		p.setX(p.getX()-h);
	}
	
	if(p.getX()[1] < 0){
		h[1] = domainSize[1];
		h[0] = 0;
		h[2] = 0;
		p.setX(p.getX()+h);
	}
	else if(p.getX()[1] > domainSize[1]){
		h[1] = domainSize[1];
		h[0] = 0;
		h[2] = 0;
		p.setX(p.getX()-h);
	}
	
	if(p.getX()[2] < 0){
		h[2] = domainSize[2];
		h[0] = 0;
		h[1] = 0;
		p.setX(p.getX()+h);
	}
	else if(p.getX()[2] > domainSize[2]){
		h[2] = domainSize[2];
		h[0] = 0;
		h[1] = 0;
		p.setX(p.getX()-h);
	}
	
	LOG4CXX_TRACE(logger, "Moved particle to " << p.getX().toString() );
}

void calcPeriodicBoundary(Particle &p1, Particle p2)
{	
	LOG4CXX_TRACE(logger, "Virtual particle-copy created at " << p2.getX().toString() << " in pair with " << p1.getX().toString() );
	calculateF(p1,p2);
}

void applyGravitation( Particle& p )
{
	p.setF( p.getF() + p.getM() * gravitation );
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
	LOG4CXX_DEBUG(logger, "Plot particles of Iteration " << iteration);

	vtk_writer.initializeOutput(particleContainer->size());

	particleContainer->applyToSingleParticles( plotParticle );

	vtk_writer.writeFile(outputFileName, iteration);
}
