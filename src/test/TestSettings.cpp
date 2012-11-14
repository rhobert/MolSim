/*
 * TestSettings.cpp
 *
 *  Created on: 12.11.2012
 *      Author: chris
 */

#include "TestSettings.h"
#include "ParticleContainerTest.h"
//#include "ParticleContainerTest.cpp"
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/TestResult.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestRunner.h>
#include <cppunit/BriefTestProgressListener.h>


TestSettings::TestSettings(){

}

int TestSettings::runTest(string name){

	// Informiert Test-Listener ueber Testresultate
	    CPPUNIT_NS :: TestResult testresult;

	    // Listener zum Sammeln der Testergebnisse registrieren
	    CPPUNIT_NS :: TestResultCollector collectedresults;
	    testresult.addListener (&collectedresults);

	    // Listener zur Ausgabe der Ergebnisse einzelner Tests
	    CPPUNIT_NS :: BriefTestProgressListener progress;
	    testresult.addListener (&progress);

	    // Test-Suite ueber die Registry im Test-Runner einfuegen
	    CPPUNIT_NS :: TestRunner testrunner;
	    testrunner.addTest (CPPUNIT_NS :: TestFactoryRegistry :: getRegistry (name).makeTest ());
	    testrunner.run (testresult);

	    // Resultate im Compiler-Format ausgeben
	    CPPUNIT_NS :: CompilerOutputter compileroutputter (&collectedresults, std::cerr);
	    compileroutputter.write ();

	    // Rueckmeldung, ob Tests erfolgreich waren
	    return collectedresults.wasSuccessful () ? 0 : 1;
}
