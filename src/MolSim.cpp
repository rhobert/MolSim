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

#include <omp.h>
#include <papi.h>

#include "input/InputParameters.h"
#include "PhaseSpace.h"
#include "FileReader.h"
#include "outputWriter/VTKWriter.h"

#include "particleContainer/ParticleContainer.h"
#include "particleContainer/SimpleParticleContainer.h"
#include "particleContainer/LinkedCellParticleContainer.h"

#include "ParticleGenerator.h"
#include "Membrane.h"
#include "MembraneParticle.h"

#include "Thermostat.h"
#include "Statistics.h"

#include "test/TestSettings.h"

#define root6of2 1.122462048
#define root3of2 1.25992105

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
 * @brief Calculate the force between two particles with the Lennard-Jones potential
 *
 * @param p1 First particle
 * @param p2 Second particle
 *
 * @return Force between particles
**/
utils::Vector<double, 3> lennardJonesPotential(Particle& p1, Particle& p2);

/**
 * @brief Calculate the force between two particles with the smoothed Lennard-Jones potential
 *
 * @param p1 First particle
 * @param p2 Second particle
 *
 * @return Force between particles
**/
utils::Vector<double, 3> smoothedLennardJonesPotential(Particle& p1, Particle& p2);

/**
 * @brief Calculate and set the membrane force between this particle and its neighbours (only if the particle is a MembranParticle)
 *
 * @param p Particle to modfiy
**/
void calcMembraneForces( Particle& p );

/**
 * @brief Calculate and set the membrane force between the particle and on of his neighbours
 *
 * @param type type of neighbourship
 * @param p1 Particle to modfiy
 * @param p2 Neighbour particle
**/
void calcMembraneForce( int type, MembraneParticle& p1, MembraneParticle& p2 );

/**
 * @brief Calcuation for the reflection boundary condition for lennard jones potential
 *
 * @param boundary Boundary on which boundary conditions are applied
 * @param p Particle to reflect
**/
void calcReflection (int boundary, Particle& p);

/**
 * @brief Calcuation for the periodic boundary condition (move halo particles to their destinied boundary cells)
 *
 * @param boundary Boundary on which boundary conditions are applied
 * @param p Particle to move
**/
void calcPeriodicHalo(int boundary, Particle &p);

/**
 * @brief Calcuation for the periodic boundary condition
 *
 * @param p1 Boundary particle to modfiy
 * @param p2 Boundary particle of the opposite boundary
**/
void calcPeriodicBoundary(Particle &p1, Particle p2);


/**
 * @brief Set gravitation as static force on this particle
 *
 * @param p Particle to modify
**/
void setGravitation(Particle& p);

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

template <class Type1, class Type2, class Type3>
struct Triple
{
	Type1 first;
	Type2 second;
	Type3 third;
};

struct SmoothedLennardJonesPotentialParameter
{
	double cutoff;
	double cutoff3_rl;
	double cutoff_rl_exp3;
} sLJparamter;

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
 * @brief Parameter for smoothed Lennard-Jones potential
**/ 
double r_l = 0;

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
* @brief Controls use of membrane force calculation
**/
bool membraneOn = false;

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
	
	list< Triple<double, list<MembraneParticle*>*, utils::Vector<double,3> > > staticForcesList;

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
	list<Particle*> particles;

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
		
		if ( cuboid.membrane().present() )
		{
			membraneOn = true;
			PSE_Molekulardynamik_WS12::membrane_t membrane = cuboid.membrane().get();
			
			double stiffnessConstant = membrane.stiffnessConstant();
			double averageBondLength = membrane.averageBondLength();
			list<utils::Vector<int,3> > positions;
			list<MembraneParticle*>* positionParticles = new list<MembraneParticle*>;
			
			utils::Vector<double,3> F (0.0);
			
			if ( membrane.staticForce().present() )
			{
				PSE_Molekulardynamik_WS12::staticForce_t force = membrane.staticForce().get();
					
				F[0] = force.F().x();
				F[1] = force.F().y();
				F[2] = force.F().z();
				
				for ( PSE_Molekulardynamik_WS12::nonNegativeIntegerVectorList_t::position_const_iterator i = force.positions().position().begin();
					i != force.positions().position().end();
					i++ )
				{
					utils::Vector<int,3> position;
					
					position[0] = i->x();
					position[1] = i->y();
					position[2] = i->z();
					
					positions.push_back( position );
				}
				
				if ( force.timeEffective().present() )
				{
					Triple<double, list<MembraneParticle*>*, utils::Vector<double,3> > staticForceReset;
					
					staticForceReset.first = force.timeEffective().get();
					staticForceReset.second = positionParticles;
					staticForceReset.third = -1.0*F;
					
					staticForcesList.push_back( staticForceReset );
				}
			}
			
			generateCuboidMembrane(particles, position, velocity, dimensions, distance, mass, sigma, epsilon, type, stiffnessConstant, averageBondLength, positions, *positionParticles);
			
			for ( list<MembraneParticle*>::iterator i = positionParticles->begin();
				 i != positionParticles->end();
				 i++ )
			{
				(*i)->setStaticF(F);	
			}
			
			LOG4CXX_INFO(logger, "Cuboid is membrane with stiffness constant " << stiffnessConstant << " and average bond length " << averageBondLength)
		}
		else
		{
			generateCuboid(particles, position, velocity, dimensions, distance, mass, sigma, epsilon, type);
		}
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
	
	Statistics statistics (*particleContainer);
	
	LOG4CXX_DEBUG(logger, "Detect if gravitation is specified");
	
	if ( simulation->gravitation().present() )
	{
		gravitation[0] = simulation->gravitation().get().x();
		gravitation[1] = simulation->gravitation().get().y();
		gravitation[2] = simulation->gravitation().get().z();
		
		LOG4CXX_INFO(logger, "Gravitation is set to " << gravitation.toString());
	}
	
	LOG4CXX_DEBUG(logger, "Detect if parameters for Thermostat are specified");
	
	if ( simulation->thermostat().present() )
	{
		thermostat = new Thermostat( *particleContainer, dimensionCount );
		thermostatOn = true;
		
		if ( simulation->thermostat().get().initialT().present() )
		{
			double initialT = simulation->thermostat().get().initialT().get();
			thermostat->initializeTemperature( initialT );
			
			LOG4CXX_INFO(logger, "Superpose velocity of particles to get initial temperature with temperature " << initialT << " in " << dimensionCount << " dimensions");
		}

		if ( simulation->thermostat().get().frequency().present() )
		{
			double frequency = simulation->thermostat().get().frequency().get();
			double targetT = 0;
			double deltaT = 0;
			
			if ( simulation->thermostat().get().targetT().present() )
			{
				targetT = simulation->thermostat().get().targetT().get();
			}
			else if ( simulation->thermostat().get().initialT().present() )
			{
				targetT = simulation->thermostat().get().initialT().get();
			}
			
			if ( simulation->thermostat().get().deltaT().present() )
			{
				deltaT = simulation->thermostat().get().deltaT().get();
			}
			
			thermostat->regulate(frequency, targetT, deltaT);
			LOG4CXX_DEBUG(logger, "Regulation frequency for Thermostat is " << frequency );
		}
		
		
	}
	else
	{
		LOG4CXX_DEBUG(logger, "No parameters for Brownian Motion or Thermostat are specified");
	}

	// Check potentialName and set force calculation
	if ( simulation->potential().gravitational().present() )
	{
		forceCalc = gravitationalPotential;
		LOG4CXX_INFO(logger, "Set potential for force calulation to gravitational potential" );
	}
	else if ( simulation->potential().lennardJones().present() )
	{
		forceCalc = lennardJonesPotential;
		LOG4CXX_INFO(logger, "Set potential for force calulation to Lennard-Jones potential" );
	}
	else if ( simulation->potential().smoothedLennardJones().present() )
	{
		forceCalc = smoothedLennardJonesPotential;
		r_l = simulation->potential().smoothedLennardJones().get().rl();
		LOG4CXX_INFO(logger, "Set potential for force calulation to smoothed Lennard-Jones potential" );
		
		sLJparamter.cutoff = cutoff;
		sLJparamter.cutoff3_rl = 3.0*cutoff - r_l;
		sLJparamter.cutoff_rl_exp3 = pow(cutoff - r_l, 3.0);
	}

	// Begin Simulation
	
	
	// set up static gravitation
	particleContainer->applyToSingleParticles ( setGravitation );
	
	// the forces are needed to calculate x, but are not given in the input file.
	particleContainer->applyToParticlePairs ( calculateF );
	
	// Init iteration variables
	double start_time = 0;
	double current_time;
	int iteration = 0;
	int end_iteration = (int) ( (end_time - start_time) / delta_t );
	
	double totalTime = 0;
	double timeDiff;
	timeval t0, t1; 
	
	double doneTime = totalTime;
	double doneIteration = iteration;
	
	
	#define NUM_EVENTS 4
	int events[NUM_EVENTS] = { PAPI_TOT_CYC, PAPI_L1_TCM, PAPI_L2_TCM, PAPI_L3_TCM  };
	long_long values[NUM_EVENTS];
	int retval;
	char EventCodeStr[PAPI_MAX_STR_LEN];
	
	retval = PAPI_library_init(PAPI_VER_CURRENT);
	
	if (retval != PAPI_VER_CURRENT) 
	{
		fprintf(stderr, "PAPI library init error!\n");
		exit(1);
	}
	
	
	LOG4CXX_INFO(logger, "Start simulation from time " << start_time << " to " << end_time << " with time steps " << delta_t);
	LOG4CXX_INFO(logger, "Expecting " << end_iteration << " simulation steps" );
	if ( writeFrequency > 0 )
		LOG4CXX_INFO(logger, "Expecting " << ceil( ( (double) end_iteration ) / writeFrequency ) << " output files" );
	
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
		
		// check if we should reset static forces
		for ( list< Triple<double, list<MembraneParticle*>*, utils::Vector<double,3> > >::iterator i = staticForcesList.begin();
			 i != staticForcesList.end();
			 i++ )
		{
			if ( i->first <= current_time )
			{
				for ( list<MembraneParticle*>::iterator j = i->second->begin();
					 j != i->second->end();
					 j++ )
				{
					MembraneParticle& mp = **j;
					
					mp.setStaticF( mp.getStaticF() + i->third );
				}
				
				LOG4CXX_DEBUG(logger, "Reset static forces");
				
				i = staticForcesList.erase(i);
			}
		}
		
		if (writeFrequency != 0 && iteration % writeFrequency == 0) {
			plotParticles(iteration);
		}
		
		if ( iteration % 1000 == 0 )
		{
			statistics.beginCalcDiffusion();
		}
		
		if ( iteration % 1000 == 1 )
		{
			double diffusion = statistics.endCalcDiffusion();
			double* rdf = statistics.calcRDF(50, cutoff);
			
			LOG4CXX_INFO(logger, "Diffusion " << current_time <<  ": " << diffusion );
			
			string rdfOut = "";
			stringstream ss;
			
			for ( int i = 0; i < 10; i++ )
			{
				ss << "(" << i  << ", " << rdf[i] << ") ";
			}
			
			LOG4CXX_INFO(logger, "RDF " << current_time << ": " << ss.str() );
		}
		
		// calculate new x
		particleContainer->applyToSingleParticles ( calculateX );

		// Set up new force (Could be modified by boundary conditions)
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
		
		if ( membraneOn )
			particleContainer->applyToSingleParticles ( calcMembraneForces );
		
//		PAPI_start_counters(events, NUM_EVENTS);
		
		// calculate new v
		particleContainer->applyToSingleParticles ( calculateV );
		
//		PAPI_stop_counters(values, NUM_EVENTS);
/*		
		for ( int i = 0; i < NUM_EVENTS; i++ )
		{
			for ( int i = 0; i < NUM_EVENTS; i++ )
			{
				if (PAPI_event_code_to_name(events[i], EventCodeStr) == PAPI_OK)
				{
//					LOG4CXX_DEBUG(logger, EventCodeStr << ": " << values[i]);
					printf ( "%s: %lld \n", EventCodeStr, values[i] );
				}
			}
		}
*/	
		gettimeofday(&t1, NULL);
		timeDiff = (t1.tv_sec - t0.tv_sec) + (t1.tv_usec - t0.tv_usec) * 1e-6;
		totalTime += timeDiff;
		
		iteration++;
		
		LOG4CXX_TRACE(logger, "Iteration " << iteration << " finished in " << timeDiff << " seconds");
		
		if ( totalTime - doneTime >= 15 || iteration == 1 )
		{
			double remainingIterations = (double) (end_iteration - iteration);
			double deltaIteration = (double) (iteration - doneIteration);
			
			LOG4CXX_DEBUG(logger, "Expected remaining runtime " << ( remainingIterations / deltaIteration * (totalTime - doneTime) ) << "s"  
				<< " (iteration " << iteration << "/" << end_iteration << ")");
			
			doneTime = totalTime;
			doneIteration = iteration;
		}
	}

	LOG4CXX_INFO(logger, "End simulation in " << totalTime << " seconds" );
	
	if ( iteration > 0 )
		LOG4CXX_INFO(logger, (totalTime / (double) iteration) << " per iteration");
	
	if ( simulation->outputPhaseSpace().present() )
	{
		string phaseSpaceFilename = simulation->outputPhaseSpace().get();
		
		LOG4CXX_INFO(logger, "Write phase space to " << phaseSpaceFilename );
		
		list<Particle*> pList = particleContainer->getParticles();
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
/*
utils::Vector<double, 3> lennardJonesPotential(Particle& p1, Particle& p2)
{
	double sigma;
	double epsilon;

	if ( p1.getType() == p2.getType() )
	{
		sigma = p1.getSigma();
		epsilon = p1.getEpsilon();
	}
	else
	{
		sigma = ( p1.getSigma() + p2.getSigma() ) / 2.0;
		epsilon = sqrt( p1.getEpsilon() * p2.getEpsilon() );
	}
	
	double sigmaExp3 = sigma*sigma*sigma;
	double sigmaExp6 = sigmaExp3 * sigmaExp3;
	
	//difference between coordinates of p1 and p2
	utils::Vector<double, 3> x1_x2 = p2.getX() - p1.getX();
	
	double squareSumExp_1 = 1 / ( x1_x2.innerProduct() );
	double squareSumExp_3 = squareSumExp_1 * squareSumExp_1 * squareSumExp_1;
	
	double temp = sigmaExp6 * squareSumExp_3;
	
	//force between p1 and p2
	utils::Vector<double, 3> F1_F2 = 
		24.0 * epsilon * temp * squareSumExp_1 * 
		(1.0 - 2.0 * temp) * 
		x1_x2
	;

	return F1_F2;
}
*/

utils::Vector<double, 3> lennardJonesPotential(Particle& p1, Particle& p2)
{
	//difference between coordinates of p1 and p2
	utils::Vector<double, 3> x1_x2 = p2.getX() - p1.getX();
	double squareSum = x1_x2.innerProduct();
	
	if ( membraneOn )
	{
		MembraneParticle* mp1 = MembraneParticle::castMembraneParticle(p1);
	
		if ( mp1 != NULL )
		{
			MembraneParticle* mp2 = MembraneParticle::castMembraneParticle(p2);
			
			if ( mp2 != NULL )
			{
				if ( MembraneParticle::sameMembrane(mp1,mp2) )
				{
					double d = root6of2 * p1.getSigma();
					d = d*d;
					
					if ( squareSum > d || MembraneParticle::neighboured(mp1,mp2) )
						return utils::Vector<double, 3>(0.0);
				}
			}
		}
	}
	
	double sigma;
	double epsilon;
	
	if ( p1.getSigma() == p2.getSigma() )
		sigma = p1.getSigma();
	else
		sigma = ( p1.getSigma() + p2.getSigma() ) / 2.0;
	
	if ( p1.getEpsilon() == p2.getEpsilon() )
		epsilon = p1.getEpsilon();
	else
		epsilon = sqrt( p1.getEpsilon() * p2.getEpsilon() );
		
	
	double sigmaExp3 = sigma*sigma*sigma;
	double sigmaExp6 = sigmaExp3 * sigmaExp3;
	
	double squareSumExp_1 = 1 / squareSum;
	double squareSumExp_3 = squareSumExp_1 * squareSumExp_1 * squareSumExp_1;
	
	double temp = sigmaExp6 * squareSumExp_3;
	
	//force between p1 and p2
	utils::Vector<double, 3> F1_F2 = 
		24.0 * epsilon * temp * squareSumExp_1 * 
		(1.0 - 2.0 * temp) * 
		x1_x2
	;

	return F1_F2;
}

utils::Vector<double, 3> smoothedLennardJonesPotential(Particle& p1, Particle& p2)
{
	//difference between coordinates of p1 and p2
	utils::Vector<double, 3> x1_x2 = p2.getX() - p1.getX();
	double squareSum = x1_x2.innerProduct();
	
	if ( membraneOn )
	{
		MembraneParticle* mp1 = MembraneParticle::castMembraneParticle(p1);
	
		if ( mp1 != NULL )
		{
			MembraneParticle* mp2 = MembraneParticle::castMembraneParticle(p2);
			
			if ( mp2 != NULL )
			{
				if ( MembraneParticle::sameMembrane(mp1,mp2) )
				{
					double d = root6of2 * p1.getSigma();
					d = d*d;
					
					if ( squareSum > d || MembraneParticle::neighboured(mp1,mp2) )
						return utils::Vector<double, 3>(0.0);
				}
			}
		}
	}
	
	double sigma;
	double epsilon;
	
	if ( p1.getSigma() == p2.getSigma() )
		sigma = p1.getSigma();
	else
		sigma = ( p1.getSigma() + p2.getSigma() ) / 2.0;
	
	if ( p1.getEpsilon() == p2.getEpsilon() )
		epsilon = p1.getEpsilon();
	else
		epsilon = sqrt( p1.getEpsilon() * p2.getEpsilon() );
	
	double l2Norm = sqrt(squareSum);
	
	double sigmaExp3 = sigma*sigma*sigma;
	double sigmaExp6 = sigmaExp3 * sigmaExp3;
	
	double squareSumExp_1 = 1 / squareSum;
	double squareSumExp_3 = squareSumExp_1 * squareSumExp_1 * squareSumExp_1;
	
	double temp = sigmaExp6 * squareSumExp_3;
	
	/* begin smoothing */
	
	if ( l2Norm >= r_l )
	{
		double S1_S2 = 
			1.0 - (l2Norm - r_l)*(l2Norm - r_l) *
			sLJparamter.cutoff3_rl - 2.0 * l2Norm /	sLJparamter.cutoff_rl_exp3;
		;
		
		temp = temp * S1_S2;
	}
	
	/* end smoothing */
	
	//force between p1 and p2
	utils::Vector<double, 3> F1_F2 = 
		24.0 * epsilon * temp * squareSumExp_1 * 
		(1.0 - 2.0 * temp) * 
		x1_x2
	;

	return F1_F2;
}

void calcMembraneForces( Particle& p )
{
	MembraneParticle* mp = MembraneParticle::castMembraneParticle(p);
	
	if (mp != NULL)
	{
		mp->applyToNeighbours( calcMembraneForce );
	}
}

void calcMembraneForce( int type, MembraneParticle& p1, MembraneParticle& p2 )
{
	utils::Vector<double, 3> x1_x2 = p2.getX() - p1.getX();
	utils::Vector<double, 3> F1_F2 = p1.getStiffnessConstant() * ( x1_x2.L2Norm() - p1.getAverageBondLengthTyped(type) ) / x1_x2.L2Norm() * x1_x2;

	//sum forces
	p1.setF( p1.getF() + F1_F2 );
}

void calcReflection (int cBoundary, Particle& p)
{
	utils::Vector<double,3> x;
	double d = root6of2 * p.getSigma();
	int boundaryDimension = cBoundary / 2;
	
//	LOG4CXX_TRACE(logger, "Reflection is checked for " << p.getX().toString() );
		
	if ( domainSize[boundaryDimension] != 0 )
	{
		if ( boundary[cBoundary] == PSE_Molekulardynamik_WS12::boundary_t::reflecting )
		{
			if ( p.getX()[boundaryDimension] <= d  )
			{
				x = p.getX();
				x[boundaryDimension] = 0;
				
				LOG4CXX_TRACE(logger, "CounterParticle for " << p.getX().toString() << " is " << x.toString() );
				Particle counterParticle ( x, utils::Vector<double,3>(0.0), p.getM(), p.getSigma(), p.getEpsilon(), p.getType() );
				calculateF( p, counterParticle );
			}
		}
		
		if ( boundary[cBoundary] == PSE_Molekulardynamik_WS12::boundary_t::reflecting )
		{
			if ( domainSize[boundaryDimension] - p.getX()[boundaryDimension] <= d  )
			{
				x = p.getX();
				x[boundaryDimension] = domainSize[boundaryDimension];
				
				LOG4CXX_TRACE(logger, "CounterParticle " << p.getX().toString() << " is " << x.toString() );
				Particle counterParticle ( x, utils::Vector<double,3>(0.0), p.getM(), p.getSigma(), p.getEpsilon(), p.getType() );
				calculateF( p, counterParticle );	
			}
		}
	}
}

void calcPeriodicHalo(int cBoundary, Particle &p)
{
	int boundaryDimension = cBoundary / 2;
	utils::Vector<double,3> x = p.getX();
//	utils::Vector<double,3> h;
	
	LOG4CXX_TRACE(logger, "Move particle from " << p.getX().toString() );
	
	if ( cBoundary % 2 == 0 )
	{
		x[boundaryDimension] += domainSize[boundaryDimension];
		p.setX(x);
	}
	else
	{
		x[boundaryDimension] -= domainSize[boundaryDimension];
		p.setX(x);
	}
	
	LOG4CXX_TRACE(logger, "Moved particle to " << p.getX().toString() );
}

void calcPeriodicBoundary(Particle &p1, Particle p2)
{
	LOG4CXX_TRACE(logger, "Virtual particle-copy created at " << p2.getX().toString() << " in pair with " << p1.getX().toString() );
	
	calculateF(p1,p2);
}

// We don't need this anymore -> void setNewForce(Particle& p)
void setGravitation( Particle& p )
{
	p.setStaticF( p.getStaticF() + p.getM() * gravitation );
}

void setNewForce( Particle& p )
{
	p.newF();
}

void calculateF( Particle& p1, Particle& p2 )
{
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
		delta_t * delta_t / (2.0 * p.getM()) * p.getF();
		
	p.setX(x);
}


void calculateV( Particle& p )
{
	utils::Vector<double,3> v =
		p.getV() +
		delta_t / (2.0 * p.getM()) * (p.getF() + p.getOldF());

	p.setV(v);
}

void plotParticle( Particle& p )
{
	#pragma omp critical
	vtk_writer.plotParticle(p);
}

void plotParticles(int iteration)
{
	LOG4CXX_DEBUG(logger, "Plot particles of Iteration " << iteration);

	vtk_writer.initializeOutput(particleContainer->size());
	
/*	LinkedCellParticleContainer::SingleList& particles = particleContainer->getParticles();
	
	for ( LinkedCellParticleContainer::SingleList::iterator i = particles.begin();
		 i != particles.end();
		 i++ )
	{
		plotParticle(*i);
	}
*/	
	particleContainer->applyToSingleParticles( plotParticle );

	vtk_writer.writeFile(outputFileName, iteration);
}
