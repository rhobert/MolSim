/*
 * GeneralTests.h
 *
 *  Created on: 13.11.2012
 *      Author: chris
 */

#ifndef GENERALTESTS_H_
#define GENERALTESTS_H_

#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>
#include "utils/Vector.h"
#include "Particle.h"

class GeneralTests : public CppUnit::TestFixture{

	CPPUNIT_TEST_SUITE(GeneralTests);

	CPPUNIT_TEST(testVandX);

	CPPUNIT_TEST_SUITE_END();


private:

	utils::Vector<double,3> F;
	utils::Vector<double,3> v;
	utils::Vector<double,3> x;
	Particle p;


public:
	GeneralTests();
	void init(Particle particle);
	void setUp();
	void tearDown();
	void testVandX();

};


#endif /* GENERALTESTS_H_ */
