 
#include "Cell.h"

using namespace std;

Cell::Cell ()
{
	particles = Cell::SingleList();
	neighbours = Cell::CellList();
	periodicNeighbours = list<pair<Cell*,bool*> >();
	#ifdef _OPENMP
		omp_init_lock(&lock);
	#endif
}

void Cell::addParticle( Particle& p )
{
	particles.push_back(&p);
}

int Cell::size() 
{
	return particles.size();
}