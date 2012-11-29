
#include "ParticleContainer.h"
#include "LinkedCellParticleContainer.h"
#include "utils/Vector"
#include <cassert>

log4cxx::LoggerPtr logger = log4cxx::Logger::getLogger("LinkedCellParticleContainer");

LinkedCellParticleContainer::LinkedCellParticleContainer(double rCutoff, std::vector<double, 3> domainSize)
{
	assert(rCutoff > 0.0);
	sideLength = rCutoff;
	createCells();
}

LinkedCellParticleContainer::~LinkedCellParticleContainer() {}

void LinkedCellParticleContainer::createCells()
{
	xCells = std::floor(domainSize[0] / sideLength);
	yCells = std::floor(domainSize[1] / sideLength);
	zCells = std::floor(domainSize[2] / sideLength);

	for (int i=0; i<xCells; ++i)
	{
		for (int j=0; j<yCells; ++j)
		{
			for (int k=0; k<zCells; ++k)
			{
				ParticleContainer pc;
				cells.push_back(pc);
			}
		}
	}
}


bool LinkedCellParticleContainer::isHalo(int x, int y, int z)
{

}

bool LinkedCellParticleContainer::isBoundary(int x, int y, int z)
{

}


ParticleContainer * LinkedCellParticleContainer::getCell(int x, int y, int z)
{	
	assert(x >= 0);
	assert(x < xCells);
	assert(y >= 0);
	assert(y < yCells);
	assert(z >= 0);
	assert(z < zCells);

	int index = x + xCells * (y + yCells * z);
	return &cells[index];
}

void LinkedCellParticleContainer::addParticle(Particle& p) 
{
	//determine dedicated cell
	utils::Vector<double, 3> x = p.getX();
	int x_ = floor(x[0] / sideLength);
	int y_ = floor(x[1] / sideLength);
	int z_ = floor(x[2] / sideLength);

	//add particle
	ParticleContainer *pc = getCell(x_, y_, z_);
	pc.singleList.push_back(p);
}

void LinkedCellParticleContainer::updateCellParticleList()
{	
	int index;
	utils::Vector<double, 3> coordinates;
	int x_;
	int y_;
	int z_;

	for(int x=0; i < xCells; x++)
	{
		for(int y=0; j < yCells; y++)
		{
			for(int z=0; z < zCells; z++) 
			{	
				index = x + xCells * (y + yCells * z);
				ParticleContainer &container = cells[index];
				for (int i=0; i < container.size(); i++)
				{
					Particle &particle = container[i];
					coordinates = particle.getX();
					x_ = floor(coordinates[0] / sideLength);
					y_ = floor(coordinates[1] / sideLength);
					z_ = floor(coordinates[2] / sideLength);
					ParticleContainer *pc = getCell(x_, y_, z_);

					if()
					{

					}
				}
			}
		}
	}
}