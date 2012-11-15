/*
 * TestSettings.h
 *
 *  Created on: 12.11.2012
 *      Author: chris
 */

#ifndef TESTSETTINGS_H_
#define TESTSETTINGS_H_


#include <string>

using namespace std;

class TestSettings{

public:
	
	/**
	* @brief Run this test by name
	* 
	* @param name Name of the test to run
	*/
	void runTest(string name);
	
	/**
	* @brief Run all tests
	*/
	void runTest();
};


#endif /* TESTSETTINGS_H_ */
