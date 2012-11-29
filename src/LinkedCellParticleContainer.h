#ifndef LINKEDCELLPARTICLECONTAINER_H_
#define LINKEDCELLPARTICLECONTAINER_H_

#include "ParticleContainer.h"
#include "Cell.h"
//#include "utils/Vector"
#include <vector>
#include <log4cxx/logger.h>

class LinkedCellParticleContainer {

private:
	static log4cxx::LoggerPtr logger;

	std::vector<ParticleContainer> cells;

	double sideLength;

	int xCells, yCells, zCells;

	bool isBoundary();

	bool isHalo();


public:
	LinkedCellParticleContainer();

	virtual ~LinkedCellParticleContainer();

	ParticleContainer * getCell(int x, int y, int z);

	void addParticle(Particle& p);

	void updateCellParticleList();

};

#endif /* LINKEDCELLPARTICLECONTAINER_H_ */