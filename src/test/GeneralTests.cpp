/*
 * GeneralTests.cpp
 *
 *  Created on: 13.11.2012
 *      Author: chris
 */

#include "GeneralTests.h"

GeneralTests::GeneralTests(){

}

void GeneralTests::init(Particle particle){
	p = particle;
}

void GeneralTests::setUp(){

	F = p.getF();
	v = p.getV();
	x = p.getX();

}

void GeneralTests::tearDown(){

}

void GeneralTests::testVandX(){
	bool b;
	utils::Vector<double,3> null_vec;
	null_vec[0] = 0;
	null_vec[1] = 0;
	null_vec[2] = 0;
	if(v == null_vec){
		if(x == null_vec){
			b = true;
		}
		else{
			b = false;
		}
	}
	else{
		if(x == null_vec){
			b = false;
		}
		else{
			b = true;
		}
	}

	CPPUNIT_ASSERT(b);

}
