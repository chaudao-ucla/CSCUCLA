/*
 * ALCTEmulationTreeCreator.h
 * 
 * Created on 26 August 2019
 *      Author: Chau Dao
*/

#ifndef CSCPATTERNS_INCLUDE_ALCTEmulationTreeCreator_H_
#define CSCPATTERNS_INCLUDE_ALCTEmulationTreeCreator_H_

#include "../include/Processor.h"

class ALCTEmulationAnalyzer : public Processor {
    int run(std::string inputfile, std::string outputfile, int start = 0, int end = -1);
};

#endif