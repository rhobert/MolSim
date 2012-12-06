/*
 * TestSettings.cpp
 *
 *  Created on: 12.11.2012
 *      Author: chris
 */


#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/TestResult.h>
#include <cppunit/Test.h>

#include "TestSettings.h"
#include "SimpleParticleContainerTest.h"
#include "LinkedCellParticleContainerTest.h"

// Register Test-Suites
CPPUNIT_TEST_SUITE_REGISTRATION( SimpleParticleContainerTest );
CPPUNIT_TEST_SUITE_REGISTRATION( LinkedCellParticleContainerTest );

void TestSettings::runTest()
{
	CppUnit::TextUi::TestRunner runner;
	CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
	runner.addTest( registry.makeTest() );
	runner.run();
}

int TestSettings::runTest(string name)
{
	CppUnit::TextUi::TestRunner runner;
	CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
	CppUnit::TestSuite* suite = (CppUnit::TestSuite*) registry.makeTest();
	std::vector<CppUnit::Test*> tests = suite->getTests();
	
	for ( int i = 0; i < tests.size(); i++ )
	{	
		CppUnit::Test & test = *(tests[i]);
		
		if (test.getName().compare(name) == 0)
		{
			runner.addTest( registry.makeTest() );
			runner.run();
			
			return 1;
		}
	}
	
	return 0;
}

std::vector<string> TestSettings::getTestNames()
{
	CppUnit::TextUi::TestRunner runner;
	CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
	CppUnit::TestSuite* suite = (CppUnit::TestSuite*) registry.makeTest();
	std::vector<CppUnit::Test*> tests = suite->getTests();
	std::vector<string> testNames (tests.size(), "");
	
	for ( int i = 0; i < tests.size(); i++ )
	{	
		testNames[i] = tests[i]->getName();
	}
	
	return testNames;
}
