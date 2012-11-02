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

/**
 * calculate the force for all particles
 */
void calculateF();

/**
 * calculate the position for all particles
 */
void calculateX();

/**
 * calculate the position for all particles
 */
void calculateV();

/**
 * plot the particles to a xyz-file
 */
void plotParticles(int iteration);

double delta_t;

ParticleContainer* particleContainer;

int main(int argc, char* argsv[]) 
{
	cout << "Hello from MolSim for PSE!" << endl;
	
	if (argc != 4) 
	{
		cout << "Errounous programme call! " << endl;
		cout << "./molsym filename t_end delta_t" << endl;
	}

	double start_time = 0;
	double end_time = atof(argsv[2]);
	delta_t = atof(argsv[3]);

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


void calculateF() 
{
	utils::Vector<double, 3> x1_x2;
	utils::Vector<double, 3> F1_F2;
	utils::Vector<double, 3> F;
	
	F = 0.0;
	
	for ( ParticleContainer::SingleList::iterator iterator = particleContainer->beginSingle();
		 iterator != particleContainer->endSingle();
		 iterator++ ) 
	{
		Particle& p = *iterator;
		
		p.newF(F);
	}
	
	
	for ( ParticleContainer::PairList::iterator iterator = particleContainer->beginPair();
		 iterator != particleContainer->endPair();
		 iterator++ ) 
	{		
		Particle& p1 = *(iterator->first);
		Particle& p2 = *(iterator->second);

		//difference between coordinates of p1 and p2
		x1_x2 = p1.getX() - p2.getX();

		//force between p1 and p2
		F1_F2 = p1.getM() * p2.getM() / pow(x1_x2.L2Norm(),3) * (-1) * x1_x2;
		
		//sum forces
		F = p1.getF() + F1_F2;
		p1.setF(F);
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