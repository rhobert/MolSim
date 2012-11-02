
#include "ParticleContainer.h"
#include <list>

using namespace std;

ParticleContainer::ParticleContainer(list<Particle> list) 
{
	singleList = list;
	cout << "Single-list generated!" << endl;
	pairList = createPairs(singleList);
	cout << "Pair-list generated!" << endl;
}

ParticleContainer::~ParticleContainer() 
{
	cout << "Particle-list destructed!" << endl;
}

ParticleContainer::PairList ParticleContainer::createPairs( ParticleContainer::SingleList& sList ) 
{
	ParticleContainer::PairList pList;
	
	list<Particle>::iterator it1;
	list<Particle>::iterator it2;

	for (it1 = sList.begin(); it1 != sList.end(); it1++){
		for (it2 = sList.begin(); it2 != sList.end(); it2++){
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