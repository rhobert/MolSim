
#include "Statistics.h"

#define PI4_3 4.18879

using namespace std;

log4cxx::LoggerPtr Statistics::logger(log4cxx::Logger::getLogger("Statistics"));

double Statistics::rdfDeltaR = 0.0;
int* Statistics::rdfStatistics = NULL;

Statistics::Statistics( ParticleContainer& pc )
{
	this->pc = &pc;
}

void Statistics::beginCalcDiffusion()
{
	pc->applyToSingleParticles( Statistics::saveOldX );
}

double Statistics::endCalcDiffusion()
{
	ParticleContainer::SingleList& pList = pc->getParticles();
	double size = (double) pList.size();
	
	double var = 0;
	
	for (ParticleContainer::SingleList::iterator i = pList.begin(); i != pList.end(); i++)
	{
		Particle& p = **i;

		var += (p.getX() - p.getOldX()).L2Norm();
	}
	
	var = var / size;
	
	return var;
}

double* Statistics::calcRDF(int intervalCount, double cutoff)
{
	double* statistics = new double[intervalCount];
	Statistics::rdfStatistics = new int[intervalCount];
	Statistics::rdfDeltaR = cutoff / (double) intervalCount;
	
	for ( int i = 0; i < intervalCount; i++ )
		Statistics::rdfStatistics[i] = 0;
	
	pc->applyToParticlePairs(Statistics::updateRDFStatistics);
	
	double ri = 0;
	
	for ( int i = 0; i < intervalCount; i++ )
	{ 
		statistics[i] = Statistics::rdfStatistics[i] / ( PI4_3 * ( pow(ri + Statistics::rdfDeltaR, 3.0) - pow(ri, 3.0) ) );
		ri += Statistics::rdfDeltaR;
	}
	
	return statistics;
}

void Statistics::saveOldX(Particle& p)
{
	p.saveX();
}

void Statistics::updateRDFStatistics(Particle& p1, Particle& p2)
{
	int container = (int) ((p1.getX() - p2.getX()).L2Norm() / rdfDeltaR);
	
	#ifdef _OPENMP
		#pragma omp critical
	#endif
	rdfStatistics[container] += 2;
}