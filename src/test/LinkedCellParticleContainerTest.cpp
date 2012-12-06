
#include "LinkedCellParticleContainerTest.h"


LinkedCellParticleContainer * LinkedCellParticleContainerTest::setUpParticleContainer()
{
	return new LinkedCellParticleContainer(DOMAIN_SIZE_STD, CUTOFF_RADIUS_STD);
}