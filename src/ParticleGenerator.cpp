
#include "ParticleGenerator.h"
#include "MaxwellBoltzmannDistribution.h"
#include <cassert>

void generateCuboid(std::list<Particle> &particles, utils::Vector<double, 3> x,
				utils::Vector<double, 3> v, utils::Vector<int, 3> N,
				double h, double m) 
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
				Particle p(particlePosition, v, m);
				particles.push_back(p);
				
				particlePosition[2] += h;
			}
			
			particlePosition[1] += h;
		}
		
		particlePosition[0] += h;
	}
}