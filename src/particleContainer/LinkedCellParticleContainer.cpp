
#include "LinkedCellParticleContainer.h"
#include <list>
#include <vector>
#include <cassert>
#include <iostream>

#define LINKED_CELL_THREAD_COUNT 4

using namespace std;

log4cxx::LoggerPtr LinkedCellParticleContainer::logger(log4cxx::Logger::getLogger("LinkedCellParticleContainer"));

LinkedCellParticleContainer::LinkedCellParticleContainer(utils::Vector<double, 3> domainSize, double cutoff) 
{
	LOG4CXX_INFO(logger, "Create LinkedCellParticleContainer");
	
	assert(domainSize[0] != 0 || domainSize[1] != 0|| domainSize[2] != 0);
	
	count = 0;
	cells = CellList();
	singleList = LinkedCellParticleContainer::SingleList();
	
	boundaryCells[0] = new vector<int>();
	boundaryCells[1] = new vector<int>();
	boundaryCells[2] = new vector<int>();
	boundaryCells[3] = new vector<int>();
	boundaryCells[4] = new vector<int>();
	boundaryCells[5] = new vector<int>();
	
	haloCells[0] = new vector<int>();
	haloCells[1] = new vector<int>();
	haloCells[2] = new vector<int>();
	haloCells[3] = new vector<int>();
	haloCells[4] = new vector<int>();
	haloCells[5] = new vector<int>();
	
	periodicCellPairs[0] = new CellPairList();
	periodicCellPairs[1] = new CellPairList();
	periodicCellPairs[2] = new CellPairList();
	periodicCellPairs[3] = new CellPairList();
	periodicCellPairs[4] = new CellPairList();
	periodicCellPairs[5] = new CellPairList();
	
	cellPairs = CellPairList();
	sideLength = cutoff;
	this->domainSize = utils::Vector<double,3> (0.0);
	
	cellDimensions = utils::Vector<int,3>(0);
	cellDimensions[0] = (domainSize[0] != 0) ? (int) ceil ( domainSize[0] / sideLength ) + 2 : 1;
	cellDimensions[1] = (domainSize[1] != 0) ? (int) ceil ( domainSize[1] / sideLength ) + 2 : 1;
	cellDimensions[2] = (domainSize[2] != 0) ? (int) ceil ( domainSize[2] / sideLength ) + 2 : 1;
	
	cellCount = cellDimensions[0] * cellDimensions[1] * cellDimensions[2];
	
	this->domainSize[0] = (domainSize[0] != 0) ? ceil ( domainSize[0] / sideLength ) * sideLength : 0.0;
	this->domainSize[1] = (domainSize[1] != 0) ? ceil ( domainSize[1] / sideLength ) * sideLength : 0.0;
	this->domainSize[2] = (domainSize[2] != 0) ? ceil ( domainSize[2] / sideLength ) * sideLength : 0.0;
	
	int cellId;
	int i[3];
	int j[3];
	int k[3];
	
	boundaryCellCount = 0;
	haloCellCount = 0;
	
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
				
				if ( isHalo( utils::Vector<int,3>(i) ) )
				{
					if	( i[0] == 0 &&  domainSize[0] != 0 )
					{
						haloCells[0]->push_back( getCell(i) );
					}
					if	( i[0] == cellDimensions[0]-1 &&  domainSize[0] != 0 )
					{
						haloCells[1]->push_back( getCell(i) );
					}
					
					if	( i[1] == 0 &&  domainSize[1] != 0 )
					{
						haloCells[2]->push_back( getCell(i) );
					}
					if	( i[1] == cellDimensions[1]-1 &&  domainSize[1] != 0 )
					{
						haloCells[3]->push_back( getCell(i) );
					}
					
					if	( i[2] == 0 &&  domainSize[2] != 0 )
					{
						haloCells[4]->push_back( getCell(i) );
					}
					if	( i[2] == cellDimensions[2]-1 &&  domainSize[2] != 0 )
					{
						haloCells[5]->push_back( getCell(i) );
					}
					
					haloCellCount++;
				}
				
				// Boundary cell
				
				if ( isBoundary( utils::Vector<int,3>(i) ) )
					{
					if	( i[0] == 1 &&  domainSize[0] != 0 )
					{
						boundaryCells[0]->push_back( getCell(i) );
					}
					if	( i[0] == cellDimensions[0]-2 &&  domainSize[0] != 0 )
					{
						boundaryCells[1]->push_back( getCell(i) );
					}
					
					if	( i[1] == 1 &&  domainSize[1] != 0 )
					{
						boundaryCells[2]->push_back( getCell(i) );
					}
					if	( i[1] == cellDimensions[1]-2 &&  domainSize[1] != 0 )
					{
						boundaryCells[3]->push_back( getCell(i) );
					}
					
					if	( i[2] == 1 &&  domainSize[2] != 0 )
					{
						boundaryCells[4]->push_back( getCell(i) );
					}
					if	( i[2] == cellDimensions[2]-2 &&  domainSize[2] != 0 )
					{
						boundaryCells[5]->push_back( getCell(i) );
					}
					
					boundaryCellCount++;
				}
			}
		}
	}
	
	LOG4CXX_DEBUG(logger, "Created " << boundaryCellCount << " boundary cells");
	LOG4CXX_DEBUG(logger, "Created " << haloCellCount << " halo cells");
	
	LOG4CXX_DEBUG(logger, "Create cell pairs");
	
	// Create cell pairs
	
	vector<bool> visted (cellCount, false);
	int cellVisted;
	int c = 0;
	int p = 0;
	
	for ( i[0] = 0; i[0] < cellDimensions[0]; i[0]++ )
	{	
		for ( i[1] = 0; i[1] < cellDimensions[1]; i[1]++ )
		{
			for ( i[2] = 0; i[2] < cellDimensions[2]; i[2]++ )
			{
				cellId = getCell( utils::Vector<int,3>(i) );
				vector<bool> visitedPeriodic (cellCount, false);
				visitedPeriodic[cellId] = true;
				
//				cout <<  utils::Vector<int,3>(i) << endl;
				
				for ( j[0] = i[0]-1; j[0] <= i[0]+1; j[0]++ )
				{
					for ( j[1] = i[1]-1; j[1] <= i[1]+1; j[1]++ )
					{
						for ( j[2] = i[2]-1; j[2] <= i[2]+1; j[2]++ )
						{			
							cellVisted = getCell( utils::Vector<int,3>(j) );
							
							if ( cellVisted >= 0  && !visted[cellVisted]  )
							{
								c++;
								cells[cellId].neighbours.push_back( & cells[cellVisted] );
								cellPairs.push_back ( pair<Cell*, Cell*>( &(cells[cellId]), &(cells[cellVisted] )) );
							}
							
							// Boundary cell?
							if ( cellVisted >= 0 && isBoundary( utils::Vector<int,3>(i) ) && isHalo( utils::Vector<int,3>(j) ) )
							{
								
								utils::Vector<int,3> cellPeriodic (j);
								
//								cout << "\t" << cellPeriodic << " => ";
								
								for ( int n = 0; n <= 2; n++ )
								{
									if ( cellDimensions[n] != 1 )
									{
										if ( j[n] == 0 )
										{
											cellPeriodic[n] += cellDimensions[n] - 2;
										}
										if ( j[n] == cellDimensions[n]-1 )
										{
											cellPeriodic[n] -= cellDimensions[n] - 2;
										}
									}
								}
								
								cellVisted = getCell( cellPeriodic );
								
//								cout << cellPeriodic << " - " << cellVisted << endl;
								
								if ( !( cellPeriodic == utils::Vector<int,3>(j) || visitedPeriodic[ getCell(cellPeriodic) ] || isHalo(cellPeriodic) ) )
								{														
									if	( cellPeriodic[0] == 1 && cellDimensions[0] != 1  )
									{
										periodicCellPairs[1]->push_back ( pair<Cell*, Cell*>( &(cells[cellId]), &(cells[cellVisted] )) );
									}
									if	( cellPeriodic[0] == cellDimensions[0]-2 && cellDimensions[0] != 1  )
									{
										periodicCellPairs[0]->push_back ( pair<Cell*, Cell*>( &(cells[cellId]), &(cells[cellVisted] )) );
									}
									
									if	( cellPeriodic[1] == 1 && cellDimensions[1] != 1  )
									{
										periodicCellPairs[3]->push_back ( pair<Cell*, Cell*>( &(cells[cellId]), &(cells[cellVisted] )) );
									}
									if	( cellPeriodic[1] == cellDimensions[1]-2 && cellDimensions[1] != 1  )
									{
										periodicCellPairs[2]->push_back ( pair<Cell*, Cell*>( &(cells[cellId]), &(cells[cellVisted] )) );
									}
									
									if	( cellPeriodic[2] == 1 && cellDimensions[2] != 1 )
									{
										periodicCellPairs[5]->push_back ( pair<Cell*, Cell*>( &(cells[cellId]), &(cells[cellVisted] )) );
									}
									if	( cellPeriodic[2] == cellDimensions[2]-2 && cellDimensions[2] != 1 )
									{
										periodicCellPairs[4]->push_back ( pair<Cell*, Cell*>( &(cells[cellId]), &(cells[cellVisted] )) );
									}
									
									visitedPeriodic[ getCell(cellPeriodic) ] = true;
									
									p++;
								}
							}
						}
					}
				}
				
				visted[cellId] = true;
//				cout << endl;
			}
		}
	}
	
	LOG4CXX_DEBUG(logger, "Created " << c << " cell pairs");
	LOG4CXX_DEBUG(logger, "Created " << p << " periodic cell pairs");	
}

void LinkedCellParticleContainer::addParticles( list<Particle> pList )
{
	LOG4CXX_TRACE(logger, "Add " << pList.size() << " particles to cells");
	
	int cellId;
	
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
	omp_set_num_threads(LINKED_CELL_THREAD_COUNT);
	
	#ifdef _OPENMP
	#pragma omp parallel for
	#endif
	for ( int i = 0; i < cellCount; i++ )
	{
		Cell& cell = cells[i];
		
		for ( Cell::SingleList::iterator j = cell.particles.begin(); j != cell.particles.end(); j++  )
		{
			Particle& p = *j;
			
			singleFunction(p);
		}
	}
}

void LinkedCellParticleContainer::applyToBoundaryParticles( int boundary, void (*singleFunction)(int, Particle&) )
{
	assert ( boundary >= 0 && boundary < 6 );
	
	omp_set_num_threads(LINKED_CELL_THREAD_COUNT);
	
	int len = boundaryCells[boundary]->size();
	
	#ifdef _OPENMP
	#pragma omp parallel for
	#endif
	for ( int i = 0; i < len; i++ )
	{
		Cell& cell = cells[ boundaryCells[boundary]->at(i) ];
		
		for ( Cell::SingleList::iterator j = cell.particles.begin(); j != cell.particles.end(); j++ )
		{
			Particle& p = *j;
			
			singleFunction(boundary,p);
		}
	}
}

void LinkedCellParticleContainer::applyToHaloParticles( int boundary, void (*singleFunction)(int, Particle&) )
{
	assert ( boundary >= 0 && boundary < 6 );
	
	omp_set_num_threads(LINKED_CELL_THREAD_COUNT);
	
	int len = haloCells[boundary]->size();
	
	#ifdef _OPENMP
	#pragma omp parallel for
	#endif
	for ( int i = 0; i < len; i++ )
	{
		Cell& cell = cells[ haloCells[boundary]->at(i) ];
		
		for ( Cell::SingleList::iterator j = cell.particles.begin(); j != cell.particles.end(); j++ )
		{
			Particle& p = *j;
			
			singleFunction(boundary,p);
		}
	}
}

void LinkedCellParticleContainer::applyToParticlePairs( void (*pairFunction)(Particle&, Particle&) ) 
{	
	omp_set_num_threads(LINKED_CELL_THREAD_COUNT);
	
	#ifdef _OPENMP
	#pragma omp parallel for
	#endif
	for ( int i = 0; i < cellCount; i++ )
	{
		Cell& cell1 = cells[i];
		omp_set_lock(&cell1.lock);
		
		for ( Cell::CellList::iterator j = cell1.neighbours.begin(); j != cell1.neighbours.end(); j++ )
		{
			Cell& cell2 = **j;
			
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
				omp_set_lock(&cell2.lock);
				
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
				
				omp_unset_lock(&cell2.lock);
			}
		}
		
		omp_unset_lock(&cell1.lock);
	}
}

void LinkedCellParticleContainer::applyToPeriodicBoundaryParticlePairs( int boundary, void (*pairFunction)(Particle&, Particle) )
{
	assert ( boundary >= 0 && boundary < 6 );
	
	omp_set_num_threads(LINKED_CELL_THREAD_COUNT);
	
	int len = periodicCellPairs[boundary]->size();
	
	#ifdef _OPENMP
	#pragma omp parallel for
	#endif
	for ( int i = 0; i < len; i++ )
	{
		Cell& cell1 = *(periodicCellPairs[boundary]->at(i).first);
		Cell& cell2 = *(periodicCellPairs[boundary]->at(i).second);
		
		for ( Cell::SingleList::iterator j1 = cell1.particles.begin(); j1 != cell1.particles.end(); j1++  )
		{
			Particle& p1 = *j1;
			
			for ( Cell::SingleList::iterator j2 = cell2.particles.begin(); j2 != cell2.particles.end(); j2++ )
			{
				Particle p2 (*j2);
				utils::Vector<double,3> x = p2.getX();
				int boundaryDimension = boundary / 2;
		
//				cout << "Boundary " << boundary << ": ";
//				cout << "Search partner for " << p1.getX().toString() << ", ";
//				cout << "test " << p2.getX().toString() << " => ";
				
				if ( boundary % 2 == 0 )
				{
					x[boundaryDimension] -= domainSize[boundaryDimension];
				}
				else
				{
					x[boundaryDimension] += domainSize[boundaryDimension];
				}
				
				p2.setX(x);
				
//				cout << p2.getX().toString() << endl;
				
				utils::Vector<double,3> x1_x2 = p1.getX() - p2.getX();
					
				if ( x1_x2.L2Norm() <= sideLength )
				{
					pairFunction(p1,p2);
				}
			}
		}
	}
}

void LinkedCellParticleContainer::updateContainingCells()
{
	omp_set_num_threads(LINKED_CELL_THREAD_COUNT);
	
	#ifdef _OPENMP
	#pragma omp parallel for
	#endif
	for ( int i = 0; i < cellCount; i++ )
	{
		Cell& cell = cells[i];
		
		omp_set_lock ( &cell.lock );
		
		for ( Cell::SingleList::iterator j = cell.particles.begin(); j != cell.particles.end(); j++  )
		{
			Particle& p = *j;
			int cellId = getCell( p.getX() );
			
			if ( cellId >= 0 )
			{
				if ( &cell != &(cells[cellId]) )
				{
					LOG4CXX_TRACE(logger, "Move particle " << p.getX().toString() << " to other cell");
					
					omp_set_lock ( &cells[cellId].lock );
					cells[cellId].addParticle(p);
					omp_unset_lock ( &cells[cellId].lock );
					
					j = cell.particles.erase(j);
					j--;
				}
			}
			else
			{
				j = cell.particles.erase(j);
				count--;
				LOG4CXX_DEBUG(logger, "Delete particle " << p.getX().toString() << " from outside of domain and halo");
			}
		}
		
		omp_unset_lock ( &cell.lock );
	}
}

void LinkedCellParticleContainer::deleteHaloParticles( int boundary )
{
	assert ( boundary >= 0 && boundary < 6 );
	
	for ( vector<int>::iterator i = haloCells[boundary]->begin(); i != haloCells[boundary]->end(); i++ )
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

bool LinkedCellParticleContainer::isHalo( utils::Vector<int,3> x )
{	
	if ( cellDimensions[0] != 1 && (x[0] == 0 || x[0] == cellDimensions[0]-1) ||
		cellDimensions[1] != 1 && (x[1] == 0 || x[1] == cellDimensions[1]-1) ||
		cellDimensions[2] != 1 && (x[2] == 0 || x[2] == cellDimensions[2]-1)
	)
		return true;
		
	return false;
}

bool LinkedCellParticleContainer::isBoundary( utils::Vector<int,3> x )
{	
	if ( !isHalo(x) && (
		cellDimensions[0] != 1 && (x[0] == 1 || x[0] == cellDimensions[0]-2) ||
		cellDimensions[1] != 1 && (x[1] == 1 || x[1] == cellDimensions[1]-2) ||
		cellDimensions[2] != 1 && (x[2] == 1 || x[2] == cellDimensions[2]-2)
	) )
		return true;
	
	return false;
}

LinkedCellParticleContainer::SingleList& LinkedCellParticleContainer::getParticles()
{
	singleList.clear();
	
	for ( LinkedCellParticleContainer::CellList::iterator i = cells.begin(); i != cells.end(); i++ )
	{
		Cell& cell = *i;
		
		for ( Cell::SingleList::iterator j = cell.particles.begin(); j != cell.particles.end(); j++  )
		{
			Particle& p = *j;
			
			singleList.push_back(p);
		}
	}
	
	return singleList;
}