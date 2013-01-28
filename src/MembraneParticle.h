
#ifndef MEMBRANEPARTICLE_H_
#define MEMBRANEPARTICLE_H_

#include <vector>

#include "utils/Vector.h"

#include "Particle.h"
#include "Membrane.h"

using namespace std;

class MembraneParticle: public Particle 
{
	typedef vector<MembraneParticle*> NeighbourList;
	
	Membrane* membrane;
	NeighbourList neighbours[3];
	
public:
	
	MembraneParticle(
		utils::Vector<double, 3> x_arg,
		utils::Vector<double, 3> v_arg,
		double m_arg,
		double sigma,
		double epsilon,
		int type,
		Membrane* membrane
	);
	
	virtual ~MembraneParticle();
	
	double getStiffnessConstant();
	
	double getAverageBondLength();
	
	double getAverageBondLengthTyped( int type );
	
	void addNeighbour( int type, MembraneParticle* particle );
	
	void removeNeighbour( int type, MembraneParticle* particle );
	
	void applyToNeighbours( void (*pairFunction)(int, MembraneParticle&, MembraneParticle&) );
	
	bool isNeighbour( MembraneParticle* particle );
	
	static bool neighboured( MembraneParticle* p1, MembraneParticle* p2 );
	
	static bool sameMembrane( MembraneParticle* p1, MembraneParticle* p2 );
	
	static MembraneParticle* castMembraneParticle( Particle& p );
};

inline double MembraneParticle::getStiffnessConstant()
{
	return membrane->stiffnessConstant;
}
	
inline double MembraneParticle::getAverageBondLength()
{
	return membrane->averageBondLength;
}

inline double MembraneParticle::getAverageBondLengthTyped( int type )
{
	return membrane->averageBondLengthTyped[type];
}

inline bool MembraneParticle::neighboured( MembraneParticle* p1, MembraneParticle* p2 )
{
	for ( int j = 0; j < 3; j++ )
	{
		for ( NeighbourList::iterator i = p1->neighbours[j].begin(); i != p1->neighbours[j].end(); i++ )
		{
			if ( p2 == *i )
				return true;
		}
	}
	
	return false;
}

inline bool MembraneParticle::sameMembrane( MembraneParticle* p1, MembraneParticle* p2 )
{
	return p1->membrane = p2->membrane;
}


inline MembraneParticle* MembraneParticle::castMembraneParticle( Particle& p )
{
	return dynamic_cast<MembraneParticle*>(&p);
}

#endif /* MEMBRANEPARTICLE_H_ */
