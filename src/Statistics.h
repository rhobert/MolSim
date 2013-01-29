
#ifndef STATISTICS_H_
#define STATISTICS_H_

#include <list>

#include "Particle.h"
#include "particleContainer/ParticleContainer.h"

using namespace std;

double calcDiffusion( ParticleContainer* pc );

//list<double> calcDiffusion( ParticleContainer* container );

class Statistics
{
	/**
	 * @brief Logger for Statistics class
	 */
	static log4cxx::LoggerPtr logger;
	
	ParticleContainer* pc;
	
	static int* rdfStatistics;
	
	static double rdfDeltaR;
	
	static void saveOldX(Particle& p);
	
	static void updateRDFStatistics(Particle& p1, Particle& p2);
	
public:
	
	Statistics(ParticleContainer& pc);
	
	void beginCalcDiffusion();
	
	double endCalcDiffusion();
	
	double* calcRDF(int intervalCount, double cutoff);
};

#endif /* STATISTICS_H_ */


