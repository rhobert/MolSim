

#include "ParticleGenerator.h"
#include "MaxwellBoltzmannDistribution.h"

void generateCuboid(std::list<Particle> &particles, utils::Vector<double, 3> x,
							utils::Vector<double, 3> v, utils::Vector<int, 3> n,
								double h, double m, double meanVelocity) {

	int it1;
	int it2;
	int it3;
	
	it3 = 0;
	while (it3 <n [2]) {
		
		it2 = 0;
		while (it2 < n[1]) {

			it1 = 0;
			while (it1 < n[0]) {

				utils::Vector<double, 3> particlePosition;
				
				particlePosition[0] = x[0] + it1 * h;
				particlePosition[1] = x[1] + it2 * h;
				particlePosition[2] = x[2] + it3 * h;

				Particle newParticle(particlePosition, v, m);

				MaxwellBoltzmannDistribution(newParticle, meanVelocity, 2);

				particles.push_back(newParticle);

				++it1;
			}
			++it2;
		}
		++it3;
	}
}