 
#include "Cell.h"

Cell::Cell ()
{
	particles = Cell::SingleList();
	neighbours = Cell::CellList();
	omp_init_lock(&lock);
}

void Cell::addParticle( Particle& p )
{
	particles.push_back(p);
}

int Cell::size() 
{
	return particles.size();
}