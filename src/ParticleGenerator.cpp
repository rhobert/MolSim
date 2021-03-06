
#include "ParticleGenerator.h"
#include "MembraneParticle.h"
#include "Membrane.h"
#include "MaxwellBoltzmannDistribution.h"
#include <cassert>
#include <cmath>
#include <list>

void generateCuboid(
	std::list<Particle*> &particles,
	utils::Vector<double, 3> x,
	utils::Vector<double, 3> v,
	utils::Vector<int, 3> N,
	double h,
	double m,
	double sigma,
	double epsilon,
	int type
)
{
	assert(N[0] > 0 && N[1] > 0 && N[2] > 0);
	assert(h >= 0);
	assert(m > 0);

	double particlePosition[3];
	particlePosition[0] = x[0];

	for ( int j0 = 0; j0 < N[0]; j0++ )
	{
		particlePosition[1] = x[1];

		for ( int j1 = 0; j1 < N[1]; j1++ )
		{
			particlePosition[2] = x[2];

			for ( int j2 = 0; j2 < N[2]; j2++ )
			{
				Particle* p = new Particle(particlePosition, v, m, sigma, epsilon, type);
				particles.push_back(p);

				particlePosition[2] += h;
			}

			particlePosition[1] += h;
		}

		particlePosition[0] += h;
	}
}

void generateCuboidMembrane(
	std::list<Particle*> &particles, 
	utils::Vector<double, 3> x,
	utils::Vector<double, 3> v, 
	utils::Vector<int, 3> N,
	double h,
	double m,
	double sigma,
	double epsilon,
	int type,
	double stiffnessConstant,
	double averageBondLength,
	std::list<utils::Vector<int,3> >& positions,
	std::list<MembraneParticle*>& positionParticles
)
{
	assert(N[0] > 0 && N[1] > 0 && N[2] > 0);
	assert(h >= 0);
	assert(m > 0);

	Membrane* membrane = new Membrane(stiffnessConstant, averageBondLength);
	
	double particlePosition[3];
	particlePosition[0] = x[0];
	
	MembraneParticle* grid [N[0]] [N[1]] [N[2]];
	
	// Generate particles and grid
	
	for ( int i0 = 0; i0 < N[0]; i0++ )
	{
		particlePosition[1] = x[1];

		for ( int i1 = 0; i1 < N[1]; i1++ )
		{
			particlePosition[2] = x[2];

			for ( int i2 = 0; i2 < N[2]; i2++ )
			{
				// Initialize particle
				MembraneParticle* p = new MembraneParticle(particlePosition, v, m, sigma, epsilon, type, membrane);
				
				// Assign particle to grid position
				grid[i0][i1][i2] = p;
				
				// Add particle to list
				particles.push_back(p);
				
				particlePosition[2] += h;
			}

			particlePosition[1] += h;
		}

		particlePosition[0] += h;
	}
	
	// Assign neighbours to particles
	
	for ( int i0 = 0; i0 < N[0]; i0++ )
	{
		for ( int i1 = 0; i1 < N[1]; i1++ )
		{
			for ( int i2 = 0; i2 < N[2]; i2++ )
			{
				MembraneParticle& p1 = *grid[i0][i1][i2];
				
				// Check all possible neighbours
				
				for ( int j0 = -1; j0 <= 1; j0++ )
				{
					if ( i0+j0 < 0 || i0+j0 >= N[0] )
						continue;
					
					for ( int j1 = -1; j1 <= +1; j1++ )
					{
						if ( i1+j1 < 0 || i1+j1 >= N[1] )
							continue;
						
						for ( int j2 = -1; j2 <= 1; j2++ )
						{
							if ( i2+j2 < 0 || i2+j2 >= N[2] )
								continue;
							
							if ( j0 == 0 && j1 == 0 && j2 == 0 )
								continue;
							
							MembraneParticle& p2 = *grid[i0+j0][i1+j1][i2+j2];
							/*
							if ( j0 == 0 || j1 == 0 || j2 == 0 )
							{
								if ( j0 == 0 && (j1 != 0 && j2 != 0) || 
									j1 == 0 && (j0 != 0 && j2 != 0) || 
									j2 == 0 && (j0 != 0 && j1 != 0) )
								{
									p1.addNeighbour( 1, &p2 );
								}
								else
								{
									p1.addNeighbour( 0, &p2 );
								}
							}
							else
							{
								p1.addNeighbour( 2, &p2 );
							}
							*/
							
							// Check type of neighbourship
							
							if ( j0 != 0 && j1 != 0 && j2 != 0 )
							{
								p1.addNeighbour( 2, &p2 );
							}
							else if ( (j0 == 0 && j1 != 0 && j2 != 0) || 
									(j0 != 0 && j1 == 0 && j2 != 0) || 
									(j0 != 0 && j1 != 0 && j2 == 0) )
							{
								p1.addNeighbour( 1, &p2 );
							}
							else
							{
								p1.addNeighbour( 0, &p2 );
							}
						}
					}
				}
			}
		}
	}
	
	// Add static forces to particles
	
	for ( std::list<utils::Vector<int,3> >::iterator i = positions.begin();
		 i != positions.end();
		 i++ )
	{
		utils::Vector<int,3>& pos = *i;
		
		assert( pos[0] >= 0 && pos[0] < N[0] && pos[1] >= 0 && pos[1] < N[1] && pos[2] >= 0 && pos[2] < N[2] );
		
		positionParticles.push_back( grid[pos[0]][pos[1]][pos[2]] );
	}
}

void generateSphere(
	std::list<Particle*> &particles,
	utils::Vector<double, 3> x,
	utils::Vector<double, 3> v,
	int N,
	int d,
	double h,
	double m,
	double sigma,
	double epsilon,
	int type
)
{
	assert(N > 0);
	assert(d > 0);
	assert(h >= 0);
	assert(m > 0);

	std::list<Particle*> particlesSphere;
	utils::Vector<double, 3> c ( x + utils::Vector<double,3>(- h * (N-1)) );
	utils::Vector<int, 3> n (2 * N - 1);

	if ( d < 3 )
	{
		n[2] = 1;
		c[2] = 0;
	}
	if ( d < 2 )
	{
		n[1] = 1;
		c[1] = 0;
	}

	generateCuboid(particlesSphere, c, v, n, h, m, sigma, epsilon, type);

	utils::Vector<double, 3> x1_x2;

	for ( std::list<Particle*>::iterator i = particlesSphere.begin(); i != particlesSphere.end(); i++ )
	{
		x1_x2 = (*i)->getX() - x;

		if ( x1_x2.L2Norm() <= (N-1)*h )
		{
			particles.push_back(*i);
		}
	}
}
