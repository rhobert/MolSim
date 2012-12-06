
#include "LinkedCellParticleContainer.h"
#include <list>
#include <vector>
#include <cassert>
#include <iostream>

using namespace std;

log4cxx::LoggerPtr LinkedCellParticleContainer::logger(log4cxx::Logger::getLogger("LinkedCellParticleContainer"));

LinkedCellParticleContainer::LinkedCellParticleContainer(utils::Vector<double, 3> domainSize, double cutoff) 
{
	LOG4CXX_INFO(logger, "Create LinkedCellParticleContainer");
	
	assert(domainSize[0] != 0 || domainSize[1] != 0|| domainSize[2] != 0);
	
	cells = CellList();
	boundaryCells = vector<int>();
	haloCells = vector<int>();
	cellPairs = CellPairList();
	sideLength = cutoff;
	this->domainSize = utils::Vector<double,3> (0.0);
	
	cellDimensions = utils::Vector<int,3>(0);
	cellDimensions[0] = (domainSize[0] != 0) ? (int) ceil ( domainSize[0] / sideLength ) + 2 : 1;
	cellDimensions[1] = (domainSize[1] != 0) ? (int) ceil ( domainSize[1] / sideLength ) + 2 : 1;
	cellDimensions[2] = (domainSize[2] != 0) ? (int) ceil ( domainSize[2] / sideLength ) + 2 : 1;
	
	cellCount = cellDimensions[0] * cellDimensions[1] * cellDimensions[2];
	
	this->domainSize[0] = ceil ( domainSize[0] / sideLength ) * sideLength;
	this->domainSize[1] = ceil ( domainSize[1] / sideLength ) * sideLength;
	this->domainSize[2] = ceil ( domainSize[2] / sideLength ) * sideLength;
	
	int cellId;
	int i[3];
	int j[3];
	
	LOG4CXX_INFO(logger, "Create " << cellCount << " cells");
	
	// Create cells
	
	for ( i[0] = 0; i[0] < cellDimensions[0]; i[0]++ )
	{	
		for ( i[1] = 0; i[1] < cellDimensions[1]; i[1]++ )
		{
			for ( i[2] = 0; i[2] < cellDimensions[2]; i[2]++ )
			{	
				Cell& cell = *(new Cell());
				cells.push_back( cell );
				
				// Halo cell
				if		( (i[0] == 0 &&  domainSize[0] != 0) || (i[1] == 0 &&  domainSize[1] != 0) || (i[2] == 0 &&  domainSize[2] != 0) || 
						  (i[0] == cellDimensions[0]-1 &&  domainSize[0] != 0) || (i[1] == cellDimensions[1]-1 &&  domainSize[1] != 0) || (i[2] == cellDimensions[2]-1 &&  domainSize[2] != 0) )
				{
					haloCells.push_back( getCell(i) );
//					cout << "H ";
				}
				// Boundary cell
				else if	( (i[0] == 1 &&  domainSize[0] != 0) || (i[1] == 1 &&  domainSize[1] != 0) || (i[2] == 1 &&  domainSize[2] != 0) || 
						  (i[0] == cellDimensions[0]-2 &&  domainSize[0] != 0) || (i[1] == cellDimensions[1]-2 &&  domainSize[1] != 0) || (i[2] == cellDimensions[2]-2 &&  domainSize[2] != 0) )
				{
					boundaryCells.push_back( getCell(i) );
//					cout << "B ";
				}
				else
				{
//					cout << "N ";
				}
			}
		}
//		cout << endl;
	}
	
	LOG4CXX_INFO(logger, "Created " << boundaryCells.size() << " boundary cells");
	LOG4CXX_INFO(logger, "Created " << haloCells.size() << " halo cells");
	
	LOG4CXX_INFO(logger, "Create cell pairs");
	
	// Create cell pairs
	
	vector<bool> visted (cellCount, false);
	int cellVisted;
	
	for ( i[0] = 0; i[0] < cellDimensions[0]; i[0]++ )
	{	
		for ( i[1] = 0; i[1] < cellDimensions[1]; i[1]++ )
		{
			for ( i[2] = 0; i[2] < cellDimensions[2]; i[2]++ )
			{
				cellId = getCell( utils::Vector<int,3>(i) );
				
				for ( j[0] = i[0]-1; j[0] <= i[0]+1; j[0]++ )
				{
					for ( j[1] = i[1]-1; j[1] <= i[1]+1; j[1]++ )
					{
						for ( j[2] = i[2]-1; j[2] <= i[2]+1; j[2]++ )
						{
							
							cellVisted = getCell( utils::Vector<int,3>(j) );
							
							if ( cellVisted >= 0  && !visted[cellVisted]  )
							{
								cellPairs.push_back ( pair<Cell*, Cell*>( &(cells[cellId]), &(cells[cellVisted] )) );
							}
						}
					}
				}
				
				visted[cellId] = true;
			}
		}
	}
	
	

}

void LinkedCellParticleContainer::addParticles( list<Particle> pList )
{
	LOG4CXX_INFO(logger, "Add particles to cells");
	
	int cellId;
	count = 0;
	
	for ( list<Particle>::iterator i = pList.begin(); i != pList.end(); i++ )
	{
		Particle& p = *i;
		
		cellId = getCell( p.getX() );
		
		if ( cellId >= 0 )
		{
			count++;
			cells[cellId].addParticle(p);
		}
	}
}

void LinkedCellParticleContainer::applyToSingleParticles( void (*singleFunction)(Particle&) )
{
	
	for ( LinkedCellParticleContainer::CellList::iterator i = cells.begin(); i != cells.end(); i++ )
	{
		Cell& cell = *i;
		
		for ( Cell::SingleList::iterator j = cell.particles.begin(); j != cell.particles.end(); j++  )
		{
			Particle& p = *j;
			
			singleFunction(p);
		}
	}
}

void LinkedCellParticleContainer::applyToBoundaryParticles( void (*singleFunction)(Particle&) )
{
	for ( vector<int>::iterator i = boundaryCells.begin(); i != boundaryCells.end(); i++ )
	{
		Cell& cell = cells[*i];
		
		for ( Cell::SingleList::iterator j = cell.particles.begin(); j != cell.particles.end(); j++ )
		{
			Particle& p = *j;
			
			singleFunction(p);
		}
	}
}

void LinkedCellParticleContainer::applyToParticlePairs( void (*pairFunction)(Particle&, Particle&) ) 
{	
	for ( LinkedCellParticleContainer::CellPairList::iterator i = cellPairs.begin(); i != cellPairs.end(); i++ )
	{
		Cell& cell1 = *(i->first);
		Cell& cell2 = *(i->second);
		
		Cell::SingleList::iterator j1;
		Cell::SingleList::iterator j2;
		
		if ( &cell1 == &cell2 )
		{	
			
			for (  j1 = cell1.particles.begin(); j1 != cell1.particles.end(); j1++  )
			{
				Particle& p1 = *j1;
				
				for ( j2 = j1, j2++; j2 != cell1.particles.end(); j2++ )
				{
					Particle& p2 = *j2;
					
					utils::Vector<double,3> x1_x2 = p1.getX() - p2.getX();
					
					if ( x1_x2.L2Norm() <= sideLength )
					{
						pairFunction(p1,p2);
					}
				}
			}
		}
		else
		{
			Cell& cell1 = *(i->first);
			Cell& cell2 = *(i->second);
			
			for ( j1 = cell1.particles.begin(); j1 != cell1.particles.end(); j1++  )
			{
				Particle& p1 = *j1;
				
				for ( j2 = cell2.particles.begin(); j2 != cell2.particles.end(); j2++  )
				{
					Particle& p2 = *j2;
					
					utils::Vector<double,3> x1_x2 = p1.getX() - p2.getX();
					
					if ( x1_x2.L2Norm() <= sideLength )
					{
						pairFunction(p1,p2);
					}
				}
			}
		}
	}
}

void LinkedCellParticleContainer::updateContainingCells()
{
	for ( LinkedCellParticleContainer::CellList::iterator i = cells.begin(); i != cells.end(); i++ )
	{
		Cell& cell = *i;
		
		for ( Cell::SingleList::iterator j = cell.particles.begin(); j != cell.particles.end(); j++  )
		{
			Particle& p = *j;
			int cellId = getCell( p.getX() );
			
			if ( cellId >= 0 )
			{
				if ( &cell != &(cells[cellId]) )
				{
					cells[cellId].addParticle(p);
					j = cell.particles.erase(j);
				}
			}
			else
			{
				j = cell.particles.erase(j);
				count--;
				LOG4CXX_DEBUG(logger, "Delete particle " << p.getX().toString() << " from outside of domain and halo");
			}
		}
	}
}

void LinkedCellParticleContainer::deleteHaloParticles()
{
	for ( vector<int>::iterator i = haloCells.begin(); i != haloCells.end(); i++ )
	{
		Cell& cell = cells[*i];
		
		if ( cell.size() > 0 )
		{
			LOG4CXX_DEBUG(logger, "Delete " << cell.size() <<  " particles from halo");
			count -= cell.size();
			cell.particles.clear();
		}
	}
}

int LinkedCellParticleContainer::size() 
{
	return count;
}

utils::Vector<double, 3> LinkedCellParticleContainer::getDomainSize()
{
	return domainSize;
}

int LinkedCellParticleContainer::getCell( utils::Vector<double,3> x )
{
	x[0] = x[0] / sideLength;
	x[1] = x[1] / sideLength;
	x[2] = x[2] / sideLength;
	
	int cellPos[] = {(int) ceil(x[0]), (int) ceil(x[1]), (int) ceil(x[2])};
	
	return getCell( utils::Vector<int,3>(cellPos) );
}

int LinkedCellParticleContainer::getCell( utils::Vector<int,3> x )
{	
	if ( x[0] < 0 || x[1] < 0 || x[2] < 0 || 
		x[0] >= cellDimensions[0] || x[1] >= cellDimensions[1] || x[2] >= cellDimensions[2] )
	{
		return -1;
	}
	
	int cell = x[0] * cellDimensions[1] * cellDimensions[2]  + x[1] * cellDimensions[2] + x[2];
	
	return cell;
}

LinkedCellParticleContainer::SingleList LinkedCellParticleContainer::getParticles()
{
	LinkedCellParticleContainer::SingleList sList;
	
	for ( LinkedCellParticleContainer::CellList::iterator i = cells.begin(); i != cells.end(); i++ )
	{
		Cell& cell = *i;
		
		for ( Cell::SingleList::iterator j = cell.particles.begin(); j != cell.particles.end(); j++  )
		{
			Particle& p = *j;
			
			sList.push_back(p);
		}
	}
	
	return sList;
}