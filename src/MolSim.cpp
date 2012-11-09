/**
 * @file
 * @brief Simulation program
 */

#include "outputWriter/VTKWriter.h"
#include "FileReader.h"
#include "ParticleContainer.h"
#include "ParticleContainer.cpp"

#include <list>
#include <cstring>
#include <cstdlib>
#include <iostream>

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
 * @brief Simulation 
**/
int main(int argc, char* argsv[]) 
{
	cout << "Hello from MolSim for PSE!" << endl;
	
	if (argc != 5) 
	{
		cout << "Errounous programme call! " << endl;
		cout << "./MolSim filename t_end delta_t force" << endl;
		return EXIT_FAILURE;
	}
	
	// Init time variables
	double start_time = 0;
	double end_time = atof(argsv[2]);
	delta_t = atof(argsv[3]);
	
	if (end_time < 0 || delta_t < 0) 
	{
		cout << "Errounous programme call! " << endl;
		cout << "./MolSim filename t_end delta_t force" << endl;
		return EXIT_FAILURE;
	}
	
	string potentialName (argsv[4]);
	
	if ( potentialName.compare("gravitational") == 0)
	{
		forceCalc = gravitationalPotential;
	}
	else if ( potentialName.compare("lenard_jones") == 0)
	{
		forceCalc = lenardJonesPotential;
	}
	else
	{
		cout << "Errounous programme call! " << endl;
		cout << "./MolSim filename t_end delta_t force" << endl;
		return EXIT_FAILURE;
	}
		

	// Read particles from file to Particle list and build ParticleContainer
	FileReader fileReader;
	list<Particle> particles;
	
	fileReader.readFile(particles, argsv[1]);
	particleContainer = new ParticleContainer(particles);
	
	// the forces are needed to calculate x, but are not given in the input file.
	calculateF();

	double current_time;
	int iteration = 0;

	 // for this loop, we assume: current x, current f and current v are known
	for ( current_time = start_time; 
		 current_time < end_time; 
		 current_time += delta_t ) 
	{
		// calculate new x
		calculateX();
		// calculate new f
		calculateF();
		// calculate new v
		calculateV();
		
		iteration++;
		cout << "Iteration " << iteration << " finished." << endl;
		
		if (iteration % 10 == 0) {
			plotParticles(iteration);
		}
	}

	cout << "output written. Terminating..." << endl;
	
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