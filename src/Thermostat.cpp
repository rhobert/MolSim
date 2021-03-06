
#include "Thermostat.h"
#include "MaxwellBoltzmannDistribution.h"

log4cxx::LoggerPtr Thermostat::logger(log4cxx::Logger::getLogger("Thermostat"));

Thermostat::Thermostat( ParticleContainer& pc, int dimensionCount )
{
	this->pc = &pc;
	this->dimensionCount = dimensionCount;
	this->deltaT = 0;
	this->targetT = 0;
	this->currentTargetT = 0;
	this->regulationFrequency = 0;

	Thermostat::dimensions = dimensionCount;
}

void Thermostat::regulate( int frequency, double targetT, double deltaT )
{
	this->regulationFrequency = frequency;
	this->targetT = targetT;
	this->currentTargetT = getTemperature();
	this->deltaT = (currentTargetT < targetT) ? deltaT : -deltaT;
}

void Thermostat::apply ( int currentIteration )
{
	LOG4CXX_TRACE(logger, "Check if temperature should be regulated" );
	
	if ( regulationFrequency > 0 && currentIteration % regulationFrequency == 0 )
	{
		double regulateT = targetT;
		
		if ( deltaT != 0.0 )
		{
			if ( abs(currentT - targetT) > abs(deltaT) && abs(currentTargetT - targetT) > abs(deltaT) )
			{
				currentTargetT += deltaT;
				regulateT = currentTargetT; 
			}
		}

		regulateTemperature(regulateT);
	}
}

void Thermostat::regulateTemperature( double targetT )
{	
	double currentT = getTemperature();

	LOG4CXX_DEBUG(logger, "Regulate temperature from " << currentT << " to " << targetT);
	
	Thermostat::beta = sqrt( targetT / currentT );
	pc->applyToSingleParticles( scaleVelocity );
}

double Thermostat::getEnergy()
{
	ParticleContainer::SingleList& pList = pc->getParticles();
	double size = (double) pList.size();
	
	double currentEnergy = 0;
	
	for (ParticleContainer::SingleList::iterator i = pList.begin(); i != pList.end(); i++)
	{
		Particle& p = **i;

		currentEnergy += p.getM() / 2.0 * p.getV().innerProduct();
	}
	
	return currentEnergy;
}
	
double Thermostat::getTemperature()
{
	double currentT = ( getEnergy() * 2.0 ) / ( ((double) ( dimensionCount * pc->size() )) * kB );
	
	return currentT;
}

double Thermostat::meanVelocity = 0.1;
int Thermostat::dimensions = 3;
double Thermostat::beta = 1;

void Thermostat::applyMaxwellBoltzmannDistribution( Particle& p )
{
	MaxwellBoltzmannDistribution(p, Thermostat::meanVelocity, Thermostat::dimensions);
}

void Thermostat::scaleVelocity( Particle& p )
{
	p.setV( Thermostat::beta * p.getV() );
}

void Thermostat::initializeTemperature( double initialT )
{
	Thermostat::meanVelocity = sqrt( kB * initialT / thermostat_mass );
	
	LOG4CXX_INFO(logger, "Mean velocity for Maxwell Boltzmann Distribution is " << Thermostat::meanVelocity );
	
	this->pc->applyToSingleParticles( Thermostat::applyMaxwellBoltzmannDistribution );
}
