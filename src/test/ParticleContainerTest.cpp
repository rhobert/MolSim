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

void ParticleContainerTest::testEquality(){
	CPPUNIT_ASSERT(*p1 == *p1);
	CPPUNIT_ASSERT(!(*p1 == *p2));

}
