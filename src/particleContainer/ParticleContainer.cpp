
#include "ParticleContainer.h"
#include <list>
#include <vector>

using namespace std;

log4cxx::LoggerPtr ParticleContainer::logger(log4cxx::Logger::getLogger("ParticleContainer"));


ParticleContainer::~ParticleContainer() 
{
	LOG4CXX_INFO(logger, "ParticleContainer destructed");
}

