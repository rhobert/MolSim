 
#include "Cell.h"


Cell::Cell ()
{
	particles = list<Particle>();
}

void Cell::addParticle( Particle& p )
{
	particles.push_back(p);
}

int Cell::size() 
{
	return particles.size();
}