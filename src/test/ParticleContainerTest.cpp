/*
 * ParticleContainerTest.cpp
 *
 *  Created on: 09.11.2012
 *      Author: chris
 */

#include "ParticleContainerTest.h"
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>

void ParticleContainerTest::setUp(){

}

void ParticleContainerTest::init(ParticleContainer& pc){
	cont = pc;
	std::pair<Particle*, Particle*> help_pair = *cont.beginPair();
	p1 = help_pair.first;
	p2 = help_pair.second;

}

void ParticleContainerTest::tearDown(){

}


/**
 * tests if the size of the ParticleContainer equals the expected size
 */
void ParticleContainerTest::testCompleteness(){
	int req_size = ((cont.size()-1) * cont.size())/2;
	int counter = 0;
	for( ParticleContainer::PairList::iterator iterator = cont.beginPair();
			 iterator != cont.endPair();
			 iterator++){
		counter++;
	}
	CPPUNIT_ASSERT(req_size == counter);
}


/**
 * tests if the ParticleContainer contains pairs {(p1,p2),(p2,p1)}
 */
void ParticleContainerTest::testForDoubles(){
	std::pair<Particle*, Particle*> pair1;
	std::pair<Particle*, Particle*> pair2;

	for(ParticleContainer::PairList::iterator iterator = cont.beginPair(); iterator != cont.endPair(); iterator++){
		pair1 = *iterator;
		for(ParticleContainer::PairList::iterator iterator2 = cont.beginPair(); iterator2 != cont.endPair(); iterator2++){
			pair2 = *iterator2;

			CPPUNIT_ASSERT((pair1.first == pair2.second) && (pair1.second = pair2.first));

		}
	}

}
