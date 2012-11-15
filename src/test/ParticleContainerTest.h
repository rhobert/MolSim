/*
 * ParticleContainerTest.h
 *
 *  Created on: 09.11.2012
 *      Author: chris
 */


#include "ParticleContainer.h"
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestFixture.h>

#ifndef PARTICLECONTAINERTEST_H_
#define PARTICLECONTAINERTEST_H_

class ParticleContainerTest : public CppUnit::TestFixture 
{

	CPPUNIT_TEST_SUITE( ParticleContainerTest );
	
	CPPUNIT_TEST( testEmpty );
	CPPUNIT_TEST( testOne );
	CPPUNIT_TEST( testTwo );
	CPPUNIT_TEST( testForDoublesPair );
	CPPUNIT_TEST( testCompletenessSingle );
	CPPUNIT_TEST( testCompletenessPair );
	
	CPPUNIT_TEST_SUITE_END();
	
	ParticleContainer * container;
	list<Particle> particles;
	Particle particle;

public:
	
	void setUp();
	
	void tearDown();
	
	void testEmpty();
	
	void testOne();
	
	void testTwo();
	
	void testForDoublesPair();
	
	void testCompletenessSingle();
	
	void testCompletenessPair();
	
};


#endif /* PARTICLECONTAINERTEST_H_ */
