/*
 * ParticleContainerTest.cpp
 *
 *  Created on: 09.11.2012
 *      Author: chris
 */

#include <stdlib.h>
#include <time.h>

#include "ParticleContainerTest.h"



using namespace std;

void ParticleContainerTest::setUp()
{
	emptyContainer = setUpParticleContainer();
	oneContainer = setUpParticleContainer();
	twoContainer = setUpParticleContainer();
	randomContainer = setUpParticleContainer();
	
	srand( time(NULL) );
	randomSize = rand() % (RANDOM_SIZE_MAX - RANDOM_SIZE_MIN + 1) + RANDOM_SIZE_MIN;
	
	utils::Vector<double,3> x (PARTICLE_POSITION_STD);
	utils::Vector<double,3> v (PARTICLE_VELOCITY_STD);
	double m = 1;
	
	particle = Particle( x, v, m );
	
	Particle p0 = particle;
	Particle p1 = particle;
	Particle p2 = particle;
	
	Particle * randomParticles = new Particle[randomSize];
	
	for (int i = 0; i < randomSize; i++)
		randomParticles[i] = particle;
	
	list<Particle> particles;
	
	emptyContainer->addParticles(particles);
	
	particles.push_back(p0);
	oneContainer->addParticles(particles);
	
	particles.clear();
	
	particles.push_back(p1);
	particles.push_back(p2);
	twoContainer->addParticles(particles);
	
	particles.clear();
	
	for (int i = 0; i < randomSize; i++)
		particles.push_back(randomParticles[i]);
	
	randomContainer->addParticles(particles);
}

void ParticleContainerTest::tearDown()
{
	delete emptyContainer;
	delete oneContainer;
	delete twoContainer;
	delete randomContainer;
}

void ParticleContainerTest::testSize()
{
	CPPUNIT_ASSERT( emptyContainer->size() == 0 );
	CPPUNIT_ASSERT( oneContainer->size() == 1 );
	CPPUNIT_ASSERT( twoContainer->size() == 2 );
	CPPUNIT_ASSERT( randomContainer->size() == randomSize );
}

void ParticleContainerTest::testParticleCount()
{
	CPPUNIT_ASSERT( countParticles( emptyContainer ) == 0 );
	CPPUNIT_ASSERT( countParticles( oneContainer ) == 1 );
	CPPUNIT_ASSERT( countParticles( twoContainer ) == 2 );
	CPPUNIT_ASSERT( countParticles( randomContainer ) == randomSize );
}

void ParticleContainerTest::testParticleValues()
{
	CPPUNIT_ASSERT( checkParticles( emptyContainer ) );
	CPPUNIT_ASSERT( checkParticles( oneContainer ) );
	CPPUNIT_ASSERT( checkParticles( twoContainer ) );
	CPPUNIT_ASSERT( checkParticles( randomContainer ) );
}

void ParticleContainerTest::testApplyToSingleParticles()
{
	emptyContainer->applyToSingleParticles( modifyParticle );
	oneContainer->applyToSingleParticles( modifyParticle );
	twoContainer->applyToSingleParticles( modifyParticle );
	randomContainer->applyToSingleParticles( modifyParticle );
	
	CPPUNIT_ASSERT( checkModfication( emptyContainer ) );
	CPPUNIT_ASSERT( checkModfication( oneContainer ) );
	CPPUNIT_ASSERT( checkModfication( twoContainer ) );
	CPPUNIT_ASSERT( checkModfication( randomContainer ) );
}

void ParticleContainerTest::testDistinctParticlePairs()
{
	emptyContainer->applyToParticlePairs( modifyNotDistinctParticlePair );
	oneContainer->applyToParticlePairs( modifyNotDistinctParticlePair );
	twoContainer->applyToParticlePairs( modifyNotDistinctParticlePair );
	randomContainer->applyToParticlePairs( modifyNotDistinctParticlePair );
	
	CPPUNIT_ASSERT( checkParticles( emptyContainer ) );
	CPPUNIT_ASSERT( checkParticles( oneContainer ) );
	CPPUNIT_ASSERT( checkParticles( twoContainer ) );
	CPPUNIT_ASSERT( checkParticles( randomContainer ) );
}

void ParticleContainerTest::testParticlePairCount()
{
	CPPUNIT_ASSERT( countParticlePairs( emptyContainer ) == 0 );
	CPPUNIT_ASSERT( countParticlePairs( oneContainer ) == 0 );
	CPPUNIT_ASSERT( countParticlePairs( twoContainer ) == 1 );
	CPPUNIT_ASSERT( countParticlePairs( randomContainer ) == randomSize * (randomSize - 1) / 2 );
}


int ParticleContainerTest::countParticles(ParticleContainer * container)
{
	ParticleContainer::SingleList particles = container->getParticles();
	
	int count = 0;
	
	for ( ParticleContainer::SingleList::iterator i = particles.begin(); i != particles.end(); i++ )
		count++;
	
	return count;
}

int ParticleContainerTest::countParticlePairs(ParticleContainer * container)
{
	container->applyToParticlePairs( ParticleContainerTest::incrementParticlePair );
	
	ParticleContainer::SingleList particles = container->getParticles();
	
	double v = 0;
	
	for ( ParticleContainer::SingleList::iterator i = particles.begin(); i != particles.end(); i++ )
	{
		Particle & p = *i;
		
		v += p.getV()[0] - PARTICLE_VELOCITY_STD;
		p.setV(PARTICLE_VELOCITY_STD);
	}
	
	int count = ( (int) ( v / PARTICLE_VELOCITY_STD ) ) / 2;
	
	return count;
}

bool ParticleContainerTest::checkParticles(ParticleContainer * container)
{
	ParticleContainer::SingleList particles = container->getParticles();
		
	for ( ParticleContainer::SingleList::iterator i = particles.begin(); i != particles.end(); i++ )
	{
		Particle & p = *i;
		
		if ( !( p == particle  ) )
			return false;
	}
	
	return true;
}

bool ParticleContainerTest::checkModfication(ParticleContainer * container)
{
	ParticleContainer::SingleList particles = container->getParticles();
	
	Particle modifiedParticle = particle;
	modifyParticle(modifiedParticle);
	
	for ( ParticleContainer::SingleList::iterator i = particles.begin(); i != particles.end(); i++ )
	{
		Particle & p = *i;
		
		if ( !( p == modifiedParticle ) )
			return false;
	}
	
	return true;
}

void ParticleContainerTest::modifyParticle(Particle & p)
{
	p.setV( 3 * p.getV() );
}

void ParticleContainerTest::modifyNotDistinctParticlePair(Particle & p1, Particle & p2)
{
	if ( &p1 == &p2 )
	{
		modifyParticle(p1);
		modifyParticle(p2);
	}
}

void ParticleContainerTest::incrementParticlePair(Particle & p1, Particle & p2)
{
	p1.setV( p1.getV() + utils::Vector<double,3>(PARTICLE_VELOCITY_STD) );
	p2.setV( p2.getV() + utils::Vector<double,3>(PARTICLE_VELOCITY_STD) );
}