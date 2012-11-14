/*
 * ParticleContainerTest.h
 *
 *  Created on: 09.11.2012
 *      Author: chris
 */


#include "ParticleContainer.h"
//#include "ParticleContainer.cpp"
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>



#ifndef PARTICLECONTAINERTEST_H_
#define PARTICLECONTAINERTEST_H_

class ParticleContainerTest : public CppUnit::TestFixture {

	CPPUNIT_TEST_SUITE(ParticleContainerTest);

	CPPUNIT_TEST(testForDoubles);
	CPPUNIT_TEST(testCompleteness);

	CPPUNIT_TEST_SUITE_END();

private:
	ParticleContainer &cont;
	Particle *p1;
	Particle *p2;

public:

	ParticleContainerTest();

	void setUp();

	void init(ParticleContainer &pc);

	void tearDown();

	void testForDoubles();

	void testCompleteness();


};


#endif /* PARTICLECONTAINERTEST_H_ */
