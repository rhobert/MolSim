/*
 * TestSettings.h
 *
 *  Created on: 12.11.2012
 *      Author: chris
 */

#ifndef TESTSETTINGS_H_
#define TESTSETTINGS_H_


#include <string>
#include <vector>

using namespace std;

class TestSettings{

public:
	
	/**
	* @brief Run this test by name
	* 
	* @param name Name of the test to run
	* 
	* @return 1 if the test was found, 0 if not
	*/
	int runTest(string name);
	
	/**
	* @brief Run all tests
	*/
	void runTest();
	
	/**
	* @brief Return names of all available tests
	* 
	* @return Vector with names of all available tests
	*/
	std::vector<string> getTestNames();
};


#endif /* TESTSETTINGS_H_ */
