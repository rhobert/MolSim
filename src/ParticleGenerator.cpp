
#include "ParticleGenerator.h"
#include "MaxwellBoltzmannDistribution.h"
#include <cassert>
#include <cmath>
#include <list>

void generateCuboid(
	std::list<Particle> &particles, 
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
				Particle p(particlePosition, v, m, sigma, epsilon, type);
				particles.push_back(p);
				
				particlePosition[2] += h;
			}
			
			particlePosition[1] += h;
		}
		
		particlePosition[0] += h;
	}
}

void generateSphere(
	std::list<Particle> &particles, 
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
	
	std::list<Particle> particlesSphere;
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
	
	for ( std::list<Particle>::iterator i = particlesSphere.begin(); i != particlesSphere.end(); i++ )
	{
		x1_x2 = (*i).getX() - x;
		
		if ( x1_x2.L2Norm() <= (N-1)*h )
		{
			particles.push_back(*i);
		}
	}
}