/* 
 * File:   MonteCarlo.h
 * Author: silvest
 *
 * Created on April 28, 2011, 5:49 PM
 */

#ifndef MONTECARLO_H
#define	MONTECARLO_H

#include <InputParser.h>
#include "MonteCarloEngine.h"

class MonteCarlo {
public:
    MonteCarlo(const std::string& ModelConf_i, const std::string& 
            MonteCarloConf_i, const std::string& OutFile_i, const std::string& JobTag_i);
//    virtual ~MonteCarlo();
    void Run();
private:
    MonteCarloEngine MCEngine;
    InputParser myInputParser;
    std::vector<ModelParameter> ModPars;
    std::vector<Observable> Obs;
    std::vector<Observable2D> Obs2D;   
    std::string ModelConf, MCMCConf, OutFile, JobTag;
    bool PrintAllMarginalized;
};

#endif	/* MONTECARLO_H */
