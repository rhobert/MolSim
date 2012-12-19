
#include "PhaseSpace.h"
#include "utils/Vector.h"

#include <fstream>
#include <sstream>

using namespace std;

log4cxx::LoggerPtr PhaseSpace::logger(log4cxx::Logger::getLogger("PhaseSpace"));

PhaseSpace::PhaseSpace() 
{
	
}

PhaseSpace::~PhaseSpace() 
{
	
}

void PhaseSpace::readPhaseSpace(std::list<Particle>& particles, istream & inStream)
{
	string line;

	while ( !inStream.eof() )
	{
		getline(inStream, line);
		
		if ( line.size() == 0 || line[0] == '#' )
			continue;
		
		istringstream in (line);
		
		utils::Vector<double,3> x;
		utils::Vector<double,3> v;
		utils::Vector<double,3> F;
		utils::Vector<double,3> oldF;
		double m;
		double epsilon;
		double sigma;
		int type;
			
		for ( int j = 0; j < 3; j++ )
			in >> x[j];
		
		for ( int j = 0; j < 3; j++ )
			in >> v[j];
		
		for ( int j = 0; j < 3; j++ )
			in >> F[j];
		
		for ( int j = 0; j < 3; j++ )
			in >> oldF[j];;
		
		in >> m;
		in >> sigma;
		in >> epsilon;
		in >> type;

		Particle p ( x, v, m, sigma, epsilon, type );
		p.setF(oldF);
		p.newF(F);
		
		particles.push_back(p);
	}
}

void PhaseSpace::writePhaseSpace(std::list<Particle>& particles, ostream & outStream)
{
	outStream.setf(ios_base::showpoint);
		
	for ( list<Particle>::iterator i = particles.begin(); i != particles.end(); i++ )
	{
		Particle& p = *i;
		
		for ( int j = 0; j < 3; j++ )
			outStream << p.getX()[j] << " ";
		
		outStream << "\t";
		
		for ( int j = 0; j < 3; j++ )
			outStream << p.getV()[j] << " ";
		
		outStream << "\t";
		
		for ( int j = 0; j < 3; j++ )
			outStream << p.getF()[j] << " ";
		
		outStream << "\t";
		
		for ( int j = 0; j < 3; j++ )
			outStream << p.getOldF()[j] << " ";
		
		outStream << "\t";
		
		outStream << p.getM() << "\t";
		outStream << p.getSigma() << "\t";
		outStream << p.getEpsilon() << "\t";
		outStream << p.getType();
		
		outStream << endl;
	}
}

void PhaseSpace::readPhaseSpace(std::list<Particle>& particles, char* filename)
{
	LOG4CXX_INFO(logger, "Read phase space to " << filename );
	
	ifstream file;
	file.open(filename, ios::in);
	
	if ( file.is_open() )
	{
		readPhaseSpace(particles, file);
		
		LOG4CXX_DEBUG(logger, "Read phase space from " << filename );
	}
	else
	{
		LOG4CXX_ERROR(logger, "Couldn't open " << filename );
	}
	
	file.close();
}
    
void PhaseSpace::writePhaseSpace(std::list<Particle>& particles, char* filename)
{
	LOG4CXX_INFO(logger, "Write phase space to " << filename );
	
	ofstream file;
	file.open(filename, ios::out);
	
	if ( file.is_open() )
	{
		writePhaseSpace(particles, file);
		
		LOG4CXX_DEBUG(logger, "Wrote phase space to " << filename );
	}
	else
	{
		LOG4CXX_ERROR(logger, "Couldn't open " << filename );
	}
	
	file.close();
}