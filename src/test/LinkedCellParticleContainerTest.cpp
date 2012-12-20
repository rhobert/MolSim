
#include "LinkedCellParticleContainerTest.h"

#define INNER_MASS 1.0
#define BOUNDARY_MASS 2.0
#define HALO_MASS 3.0


LinkedCellParticleContainer * LinkedCellParticleContainerTest::setUpParticleContainer()
{
	return new LinkedCellParticleContainer(DOMAIN_SIZE_STD, CUTOFF_RADIUS_STD);
}

void LinkedCellParticleContainerTest::setUp()
{
	ParticleContainerTest::setUp();
	
	container = setUpParticleContainer();
	
	list<Particle> particles;
	
	utils::Vector<double,3> x (0.0);
	utils::Vector<double,3> v (PARTICLE_VELOCITY_STD);
	
	innerParticle = Particle( x, v, INNER_MASS );
	boundaryParticle =  Particle( x, v, BOUNDARY_MASS );
	haloParticle =  Particle( x, v, HALO_MASS );
	
	Particle innerParticle1 = innerParticle;
	Particle innerParticle2 = innerParticle;
	Particle innerParticle3 = innerParticle;
	Particle boundaryParticle1 = boundaryParticle;
	Particle boundaryParticle2 = boundaryParticle;
	Particle haloParticle1 = haloParticle;
	Particle haloParticle2 = haloParticle;
	
	innerParticle1.setX( utils::Vector<double,3>( 0.5  * container->getDomainSize() )  );
	innerParticle2.setX( utils::Vector<double,3>( 1.1  * CUTOFF_RADIUS_STD ) );
	innerParticle3.setX( utils::Vector<double,3>( 1.2 * CUTOFF_RADIUS_STD ) );
	
	boundaryParticle1.setX(  utils::Vector<double,3>( utils::Vector<double,3>( 0.8 * CUTOFF_RADIUS_STD ) ) );
	boundaryParticle2.setX(  utils::Vector<double,3>( container->getDomainSize() ) - utils::Vector<double,3>( 0.2 * CUTOFF_RADIUS_STD ) );
	
	haloParticle1.setX(  utils::Vector<double,3>( utils::Vector<double,3>( -0.8 * CUTOFF_RADIUS_STD ) ) );
	haloParticle2.setX(  utils::Vector<double,3>( container->getDomainSize() ) + utils::Vector<double,3>( 0.1 * CUTOFF_RADIUS_STD ) );
	
	particles.push_back(innerParticle1);
	particles.push_back(innerParticle2);
	particles.push_back(innerParticle3);
	particles.push_back(boundaryParticle1);
	particles.push_back(boundaryParticle2);
	particles.push_back(haloParticle1);
	particles.push_back(haloParticle2);
	/*
	std::cout << " domain: " << container->getDomainSize().toString() << std::endl;
	std::cout << innerParticle1 << std::endl;
	std::cout << innerParticle2 << std::endl;
	std::cout << innerParticle3 << std::endl;
	std::cout << boundaryParticle1 << std::endl;
	std::cout << boundaryParticle2 << std::endl;
	std::cout << haloParticle1 << std::endl;
	std::cout << haloParticle2 << std::endl;
	*/
	container->addParticles(particles);
}


void LinkedCellParticleContainerTest::testSize(){

	utils::Vector<double,3> dom_size (5.0);
	CPPUNIT_ASSERT(container->getDomainSize() == dom_size);

}

void LinkedCellParticleContainerTest::testParticleCount(){

	CPPUNIT_ASSERT(container->size() == 7);

}

void LinkedCellParticleContainerTest::testParticleValues(){

	LinkedCellParticleContainer::SingleList sList = container->getParticles();

	utils::Vector<double,3> v_ (2.0);
	utils::Vector<double,3> x_1 (-0.8);
	utils::Vector<double,3> x_2 (0.8);
	utils::Vector<double,3> x_3 (1.1);
	utils::Vector<double,3> x_4 (1.2);
	utils::Vector<double,3> x_5 (2.5);
	utils::Vector<double,3> x_6 (4.8);
	utils::Vector<double,3> x_7 (5.1);

	Particle p1 = sList.front();
	sList.pop_front();
	Particle p2 = sList.front();
	sList.pop_front();
	Particle p3 = sList.front();
	sList.pop_front();
	Particle p4 = sList.front();
	sList.pop_front();
	Particle p5 = sList.front();
	sList.pop_front();
	Particle p6 = sList.front();
	sList.pop_front();
	Particle p7 = sList.front();

	CPPUNIT_ASSERT(p1.getX() == x_1 && p1.getV() == v_);
	CPPUNIT_ASSERT(p2.getX() == x_2 && p2.getV() == v_);
	CPPUNIT_ASSERT(p3.getX() == x_3 && p3.getV() == v_);
	CPPUNIT_ASSERT(p4.getX() == x_4 && p4.getV() == v_);
	CPPUNIT_ASSERT(p5.getX() == x_5 && p5.getV() == v_);
	CPPUNIT_ASSERT(p6.getX() == x_6 && p6.getV() == v_);
	CPPUNIT_ASSERT(p7.getX() == x_7 && p7.getV() == v_);

}

void LinkedCellParticleContainerTest::testApplyToSingleParticles(){

	Particle p = container->getParticles().front();

	container->applyToSingleParticles(modifyParticle);

	Particle p_ = container->getParticles().front();

	CPPUNIT_ASSERT(p.getX() * 2 == p_.getX());

}


/*
void LinkedCellParticleContainerTest::testDeleteHaloParticles()
{
	container->deleteHaloParticles();
	ParticleContainer::SingleList particles = container->getParticles();
	
	int count = 0;
	
	for ( ParticleContainer::SingleList::iterator i = particles.begin(); i != particles.end(); i++ )
	{
		Particle & p = *i;
		
		CPPUNIT_ASSERT(p.getM() != HALO_MASS);
		count++;
	}
	
	CPPUNIT_ASSERT ( count == 5 );
}

void LinkedCellParticleContainerTest::testApplyToBoundaryParticles()
{
	container->applyToBoundaryParticles( LinkedCellParticleContainerTest::modifyParticle );
	Particle modifiedBoundaryParticle = particle;
	modifyParticle( modifiedBoundaryParticle );
	
	ParticleContainer::SingleList particles = container->getParticles();
	
	for ( ParticleContainer::SingleList::iterator i = particles.begin(); i != particles.end(); i++ )
	{
		Particle & p = *i;
		
		CPPUNIT_ASSERT ( p.getM() != BOUNDARY_MASS || p.getV() == modifiedBoundaryParticle.getV() );
	}
}

void LinkedCellParticleContainerTest::testCutOffRadius()
{
	CPPUNIT_ASSERT( countParticlePairs( container ) == 4 );
	
	container->applyToParticlePairs( LinkedCellParticleContainerTest::moveToHalo );
	
	container->updateContainingCells();
	container->deleteHaloParticles();
	
	ParticleContainer::SingleList particles = container->getParticles();
	
	int count = 0;
	
	for ( ParticleContainer::SingleList::iterator i = particles.begin(); i != particles.end(); i++ )
	{
		Particle & p = *i;
		
		count++;
	}
	
	CPPUNIT_ASSERT ( count == 2 );
}

void LinkedCellParticleContainerTest::moveToHalo( Particle& p )
{
	p.setX( utils::Vector<double,3>( -0.5  * CUTOFF_RADIUS_STD ) );
}

void LinkedCellParticleContainerTest::moveToHalo( Particle& p1, Particle& p2 )
{
	p1.setX( utils::Vector<double,3>( -0.25 * CUTOFF_RADIUS_STD ) );
	p2.setX( utils::Vector<double,3>( -0.75 * CUTOFF_RADIUS_STD ) );
}
*/
