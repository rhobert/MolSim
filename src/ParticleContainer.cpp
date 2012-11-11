
#include "ParticleContainer.h"
#include <list>
#include <vector>

using namespace std;

log4cxx::LoggerPtr ParticleContainer::logger(log4cxx::Logger::getLogger("ParticleContainer"));

ParticleContainer::ParticleContainer(list<Particle> pList) 
{
	for ( list<Particle>::iterator it = pList.begin(); it != pList.end(); it++ )
	{	
		singleList.push_back(*it);
	}
	
	LOG4CXX_INFO(logger, "SingleList generated");
	pairList = createPairs(singleList);
	LOG4CXX_INFO(logger, "PairList generated");
}

ParticleContainer::~ParticleContainer() 
{
	LOG4CXX_INFO(logger, "ParticleContainer destructed");
}

ParticleContainer::PairList ParticleContainer::createPairs( ParticleContainer::SingleList& sList ) 
{
	ParticleContainer::PairList pList;
	
	ParticleContainer::SingleList::iterator it1;
	ParticleContainer::SingleList::iterator it2;
	
	for (it1 = sList.begin(); it1 != sList.end(); it1++){
		for (it2 = it1, it2++; it2 != sList.end(); it2++){
			if (it1 != it2){
				pair<Particle*, Particle*> help = pair<Particle*, Particle*>(&(*it1), &(*it2));
				pList.push_back(help);
			}
		}
	}
	
	return pList;
}

ParticleContainer::SingleList::iterator ParticleContainer::beginSingle() 
{
	return this->singleList.begin();
}

ParticleContainer::SingleList::iterator ParticleContainer::endSingle()
{
	return this->singleList.end();
}

ParticleContainer::PairList::iterator ParticleContainer::beginPair() 
{
	return this->pairList.begin();
}

ParticleContainer::PairList::iterator ParticleContainer::endPair() 
{
	return this->pairList.end();
}

int ParticleContainer::size() 
{
	return singleList.size();
}