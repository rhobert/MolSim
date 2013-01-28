
#include "MembraneParticle.h"
#include <cassert>

MembraneParticle::MembraneParticle(
		utils::Vector<double, 3> x_arg,
		utils::Vector<double, 3> v_arg,
		double m_arg,
		double sigma,
		double epsilon,
		int type,
		Membrane* membrane
) : Particle( x_arg, v_arg, m_arg, sigma, epsilon, type )
{
	this->membrane = membrane;
}

MembraneParticle::~MembraneParticle()
{
	for ( int j = 0; j < 3; j++ )
	{
		for ( NeighbourList::iterator i = neighbours[j].begin(); i != neighbours[j].end(); i++ )
		{
			(*i)->removeNeighbour(j, this);
		}
	}
}

void MembraneParticle::addNeighbour( int type, MembraneParticle* particle )
{
	assert( type >= 0 && type < 3 );
	
	neighbours[type].push_back(particle);
}

void MembraneParticle::removeNeighbour( int type, MembraneParticle* particle )
{
	assert( type >= 0 && type < 3 );
	
	for ( NeighbourList::iterator i = neighbours[type].begin(); i != neighbours[type].end(); i++ )
	{
		if ( particle == *i )
		{
			neighbours[type].erase(i);
			break;
		}
	}
}

void MembraneParticle::applyToNeighbours( void (*pairFunction)(int, MembraneParticle&, MembraneParticle&) )
{
	for ( int j = 0; j < 3; j++ )
	{
		for ( NeighbourList::iterator i = neighbours[j].begin(); i != neighbours[j].end(); i++ )
		{
			pairFunction ( j, *this, **i );
		}
	}
}


