/*
 * TestSettings.cpp
 *
 *  Created on: 12.11.2012
 *      Author: chris
 */

#include "TestSettings.h"
#include "ParticleContainerTest.h"
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/TestResult.h>

// Register Test-Suites
CPPUNIT_TEST_SUITE_REGISTRATION( ParticleContainerTest );

void TestSettings::runTest()
{
	CppUnit::TextUi::TestRunner runner;
	CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
	runner.addTest( registry.makeTest() );
	runner.run();
}

void TestSettings::runTest(string name)
{
	CppUnit::TextUi::TestRunner runner;
	CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
	runner.addTest( registry.makeTest() );
	runner.run(name);
}
