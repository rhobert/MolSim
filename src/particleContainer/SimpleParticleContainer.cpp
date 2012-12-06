
#include "SimpleParticleContainer.h"
#include <list>
#include <vector>

using namespace std;

log4cxx::LoggerPtr SimpleParticleContainer::logger(log4cxx::Logger::getLogger("SimpleParticleContainer"));

SimpleParticleContainer::SimpleParticleContainer() 
{
	
}

void SimpleParticleContainer::addParticles( list<Particle> pList )
{
	for ( list<Particle>::iterator it = pList.begin(); it != pList.end(); it++ )
	{	
		singleList.push_back(*it);
	}
	
	LOG4CXX_INFO(logger, "SingleList generated");
	pairList = createPairs(singleList);
	LOG4CXX_INFO(logger, "PairList generated");
}

SimpleParticleContainer::PairList SimpleParticleContainer::createPairs( SimpleParticleContainer::SingleList& sList ) 
{
	SimpleParticleContainer::PairList pList;
	
	SimpleParticleContainer::SingleList::iterator it1;
	SimpleParticleContainer::SingleList::iterator it2;
	
	for (it1 = sList.begin(); it1 != sList.end(); it1++){
		for (it2 = it1, it2++; it2 != sList.end(); it2++){
			pair<Particle*, Particle*> help (&(*it1), &(*it2));
			pList.push_back(help);
		}
	}
	
	return pList;
}


void SimpleParticleContainer::applyToSingleParticles( void (*singleFunction)(Particle&) )
{
	for ( SimpleParticleContainer::SingleList::iterator i = singleList.begin();
		 i != singleList.end();
		 i++ )
	{
		singleFunction( *i );
	}
}

void SimpleParticleContainer::applyToParticlePairs( void (*pairFunction)(Particle&, Particle&) ) 
{
	for ( SimpleParticleContainer::PairList::iterator i = pairList.begin();
		 i != pairList.end();
		 i++ )
	{
		pairFunction( *(i->first), *(i->second) );
	}
	
}

SimpleParticleContainer::SingleList SimpleParticleContainer::getParticles()
{
	return singleList;
}

int SimpleParticleContainer::size() 
{
	return singleList.size();
}