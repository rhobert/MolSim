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

TestSettings::TestSettings(){

}

void TestSettings::runTest()
{
	CppUnit::TextUi::TestRunner runner;
	CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
	runner.addTest( registry.makeTest() );
	runner.run();
}

void TestSettings::runTest(string name){

	// Informiert Test-Listener ueber Testresultate
	CppUnit :: TestResult testresult;
	
	// Listener zum Sammeln der Testergebnisse registrieren
	CppUnit :: TestResultCollector collectedresults;
	testresult.addListener (&collectedresults);

	// Listener zur Ausgabe der Ergebnisse einzelner Tests
	CppUnit :: BriefTestProgressListener progress;
	testresult.addListener (&progress);

	// Test-Suite ueber die Registry im Test-Runner einfuegen
	CppUnit :: TestRunner testrunner;
	testrunner.addTest (CppUnit :: TestFactoryRegistry :: getRegistry (name).makeTest ());
	testrunner.run (testresult);

	// Resultate im Compiler-Format ausgeben
	CppUnit :: CompilerOutputter compileroutputter (&collectedresults, std::cerr);
	compileroutputter.write ();
}
