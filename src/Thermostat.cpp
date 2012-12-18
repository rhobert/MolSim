
#include "Thermostat.h"

log4cxx::LoggerPtr Thermostat::logger(log4cxx::Logger::getLogger("Thermostat"));

Thermostat::Thermostat( double targetT, double diffT, int nMax, int dimensionCount )
{
	this->targetT = targetT;
	this->diffT = diffT;
	this->nMax = nMax;
	this->dimensionCount = dimensionCount;

    k = 1.3806488e-23;
	m = 1.0;

}

double Thermostat::initializeTemperature( int size , double initialT )
{
	initialEnergy = dimensionCount / 2.0 * size * k * initialT;

	return sqrt( 2 * initialEnergy / dimensionCount * size * m );
}

void Thermostat::regulateTemperature( ParticleContainer& pc, int iteration, int nThermostat )
{
		if( iteration > 1 && iteration %  nThermostat == 1 && iteration <= nMax )
		{	
			//targetEnergy = dimensionCount / 2 * (&pc)->size() * k * targetT;

			for (ParticleContainer::SingleList::iterator i = pc.getParticles().begin(); i != pc.getParticles().end(); i++)
			{
				Particle& p = *i;

				currentEnergy += p.getM() / 2 * p.getV().innerProduct();
			}
				
			currentT = ( currentEnergy * 2 ) / ( dimensionCount * (&pc)->size() * k );

			if ( currentT != targetT ) 
			{
				beta = sqrt( ( currentT + diffT ) / currentT );

				for (ParticleContainer::SingleList::iterator j = pc.getParticles().begin(); j != pc.getParticles().end(); j++)
				{
					Particle& p = *j;

					p.setV( p.getV() * beta );
				}
			}
		}
		else
		{
			return;
		}
}

double Thermostat::getEnergy()
{
	return currentEnergy;
}
	
double Thermostat::getTemperature()
{
	return currentT;
}