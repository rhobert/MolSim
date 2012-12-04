/*
 * ParticleContainerTest.h
 *
 *  Created on: 09.11.2012
 *      Author: chris
 */


#include "particleContainer/ParticleContainer.h"
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestFixture.h>

#ifndef PARTICLECONTAINERTEST_H_
#define PARTICLECONTAINERTEST_H_

/**
* @class ParticleContainerTest
*
* @brief Tests for class ParticleContainer
*/
class ParticleContainerTest : public CppUnit::TestFixture 
{

	CPPUNIT_TEST_SUITE( ParticleContainerTest );
/*	
	CPPUNIT_TEST( testEmpty );
	CPPUNIT_TEST( testOne );
	CPPUNIT_TEST( testTwo );
	CPPUNIT_TEST( testForDoublesPair );
	CPPUNIT_TEST( testCompletenessSingle );
	CPPUNIT_TEST( testCompletenessPair );
*/	
	CPPUNIT_TEST_SUITE_END();
	
	/**
	* @brief Instance of ParticleContainer used for tests
	*/
	ParticleContainer * container;
	
	/**
	* @brief Single particle used for tests
	*/
	Particle particle;
	
	/**
	* @brief List of particle used for tests
	*/
	list<Particle> particles;

public:
	
	/**
	* @brief Set up a test
	*/	
//	void setUp();
	
	/**
	* @brief Delete test-variables
	*/
//	void tearDown();
	
	/**
	* @brief Tests if an empty ParticleContainer has no single particles and no particle pairs
	*/
//	void testEmpty();
	
	/**
	* @brief Tests if an one element ParticleContainer has one single particles and no particle pairs
	*/
//	void testOne();
	
	/**
	* @brief Tests if a two element ParticleContainer has two single particles and one particle pair
	*/
//	void testTwo();
	
	/**
	* @brief Tests if a ParticleContainer has no duplicate particle pairs
	*/
//	void testForDoublesPair();
	
	/**
	* @brief Tests if a ParticleContainer has so much single particle like it should
	*/
//	void testCompletenessSingle();
	
	/**
	* @brief Tests if a ParticleContainer has so much particle pairs like it should
	*/
//	void testCompletenessPair();

};


#endif /* PARTICLECONTAINERTEST_H_ */
