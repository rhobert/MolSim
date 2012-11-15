/*
 * ParticleContainerTest.cpp
 *
 *  Created on: 09.11.2012
 *      Author: chris
 */

#include "ParticleContainerTest.h"

using namespace std;

void ParticleContainerTest::setUp()
{
	double x_[] = {0,0,0};
	double v_[] = {0,0,0};
	
	utils::Vector<double,3> x (x_);
	utils::Vector<double,3> v (v_);
	double m = 1;
	
	Particle p ( x, v, m );
	particle = p;
	
	particles.push_back(p);
	particles.push_back(p);
	particles.push_back(p);
	particles.push_back(p);
	particles.push_back(p);
	
	container = new ParticleContainer(particles);
}

void ParticleContainerTest::tearDown()
{
	delete container;
}

void ParticleContainerTest::testEmpty()
{
	list<Particle> empty_particles;
	ParticleContainer container (empty_particles);
	
	CPPUNIT_ASSERT(container.size() == 0);
	CPPUNIT_ASSERT(container.beginSingle() == container.endSingle());
	CPPUNIT_ASSERT(container.beginPair() == container.endPair());
}

void ParticleContainerTest::testOne()
{
	list<Particle> one_particles;
	one_particles.push_back(particle);
	ParticleContainer one_container (one_particles);
	
	CPPUNIT_ASSERT(one_container.size() == 1);
	CPPUNIT_ASSERT(++one_container.beginSingle() == one_container.endSingle());
	CPPUNIT_ASSERT(one_container.beginPair() == one_container.endPair());
}
	
void ParticleContainerTest::testTwo()
{
	list<Particle> two_particles;
	two_particles.push_back(particle);
	two_particles.push_back(particle);
	ParticleContainer container (two_particles);
	
	CPPUNIT_ASSERT(container.size() == 2);
	CPPUNIT_ASSERT(++++container.beginSingle() == container.endSingle());
	CPPUNIT_ASSERT(++container.beginPair() == container.endPair());
}
	
void ParticleContainerTest::testForDoublesPair()
{		
	pair<Particle*, Particle*> pair1;
	pair<Particle*, Particle*> pair2;

	for( ParticleContainer::PairList::iterator i = container->beginPair(); 
		i != container->endPair(); 
		i++)
	{
		pair1 = *i;
		
		for( ParticleContainer::PairList::iterator j = container->beginPair(); 
			j != container->endPair(); 
			j++)
		{
			pair2 = *j;
			
			CPPUNIT_ASSERT(pair1 == pair2 || pair1.first != pair2.first || pair1.second != pair2.second);
			CPPUNIT_ASSERT(pair1.first != pair2.second || pair1.second != pair2.first);
		}
	}
}

void ParticleContainerTest::testCompletenessSingle()
{
	CPPUNIT_ASSERT(particles.size() == container->size());
	
	int counter = 0;
	
	for( ParticleContainer::SingleList::iterator i = container->beginSingle(); 
		i != container->endSingle(); 
		i++)
	{
		counter++;
	}
	
	CPPUNIT_ASSERT(particles.size() == counter);
}

void ParticleContainerTest::testCompletenessPair()
{
	int req_size = (particles.size() - 1) * particles.size() /2;
	int counter = 0;
	
	for( ParticleContainer::PairList::iterator i = container->beginPair(); 
		i != container->endPair(); 
		i++)
	{
		counter++;
	}
	
	CPPUNIT_ASSERT(req_size == counter);
}
