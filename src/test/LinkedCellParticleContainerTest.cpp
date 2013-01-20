
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



void LinkedCellParticleContainerTest::testDeleteHaloParticles()
{
	for ( int i = 0; i < 6; i++ )
		container->deleteHaloParticles(i);
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
	for ( int i = 0; i < 6; i++ )
		container->applyToBoundaryParticles( i, LinkedCellParticleContainerTest::modifyBoundaryParticle );
	Particle modifiedBoundaryParticle = particle;
	modifyParticle( modifiedBoundaryParticle );
	
	ParticleContainer::SingleList particles = container->getParticles();
	
	for ( ParticleContainer::SingleList::iterator i = particles.begin(); i != particles.end(); i++ )
	{
		Particle & p = *i;
		
		CPPUNIT_ASSERT ( p.getM() != BOUNDARY_MASS || p.getV() == modifiedBoundaryParticle.getV() || p.getM() == 2 );
	}
}

void LinkedCellParticleContainerTest::testCutOffRadius()
{
	CPPUNIT_ASSERT( countParticlePairs( container ) == 4 );
	
	container->applyToParticlePairs( LinkedCellParticleContainerTest::moveToHalo );
	container->updateContainingCells();
	
	for ( int i = 0; i < 6; i++ )
	{
		container->deleteHaloParticles(i);
	}
	
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

void LinkedCellParticleContainerTest::modifyBoundaryParticle( int boundary, Particle& p )
{
	modifyParticle(p);
}
