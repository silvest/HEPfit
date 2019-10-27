/*
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

/**
 * @example libmode_config.cpp
 * This is an example of how to compute observables from the input parameters
 * defined in a model configuration file.
 *
 */

#include <iostream>
#include <ComputeObservables.h>
#include <InputParameters.h>
#include "StandardModel.h"

class AlsTest : public ThObservable {
public:   
    /**
     * constructor
     * @param Flavour
     */
    AlsTest(const StandardModel& SM_i) : ThObservable(SM_i) {}
    
/**
     * 
     * @brief arXiv:1206.2634v2
     * @return theoretical value of |\f$ BR(B \rightarrow \tau \nu) \f$|
     */
    double computeThValue() {

        std::cout << "5 GeV - 5 flavour" << std::endl;
        std::cout << "LO: " << SM.Als(5., LO) - SM.AlsOLD(5., LO) << std::endl;
        std::cout << "NLO: " << SM.Als(5., FULLNLO) - SM.AlsOLD(5., FULLNLO) << std::endl;
        std::cout << "NNLO: " << SM.Als(5., FULLNNLO) - SM.AlsOLD(5., FULLNNLO) << std::endl;
        std::cout << "NNNLO: " << SM.Als(5., FULLNNNLO) << std::endl;

        std::cout << "3 GeV - 4 flavour" << std::endl;
        std::cout << "LO: " << SM.Als(3., LO) - SM.AlsOLD(3., LO) << std::endl;
        std::cout << "NLO: " << SM.Als(3., FULLNLO) - SM.AlsOLD(3., FULLNLO) << std::endl;
        std::cout << "NNLO: " << SM.Als(3., FULLNNLO) - SM.AlsOLD(3., FULLNNLO) << std::endl;
        std::cout << "NNNLO (no comparison): " << SM.Als(3., FULLNNNLO) << std::endl;

        std::cout << "1 GeV - 3 flavour" << std::endl;
        std::cout << "LO: " << SM.Als(1., LO) - SM.AlsOLD(1., LO) << std::endl;
        std::cout << "NLO: " << SM.Als(1., FULLNLO) - SM.AlsOLD(1., FULLNLO) << std::endl;
        std::cout << "NNLO: " << SM.Als(1., FULLNNLO) - SM.AlsOLD(1., FULLNNLO) << std::endl;
        std::cout << "NNNLO (no comparison): " << SM.Als(1., FULLNNNLO) << std::endl;

        std::cout << "200 GeV - 6 flavour" << std::endl;
        std::cout << "LO: " << SM.Als(200., LO) - SM.AlsOLD(200., LO) << std::endl;
        std::cout << "NLO: " << SM.Als(200., FULLNLO) - SM.AlsOLD(200., FULLNLO) << std::endl;
        std::cout << "NNLO: " << SM.Als(200., FULLNNLO) - SM.AlsOLD(200., FULLNNLO) << std::endl;
        std::cout << "NNNLO (no comparison): " << SM.Als(200., FULLNNNLO) << std::endl;

        return (0.);
    }
    
protected:
    
private:
    
};


int main(void) {
    /* Define the model configuration file.                        */
    /* Here it is passed as the first argument to the executable.  */
    /* The model configuration file provides the default values of */
    /* the mandatory model parameters.                             */
    std::string ModelConf = "StandardModel.conf";
    std::string ModelName = "StandardModel";

    /* Define a map for the parameters to be varied. */
    std::map<std::string, double> DPars;
        
    /* Create objects of the classes ModelFactory and ThObsFactory */
    ModelFactory ModelF;
    ThObsFactory ThObsF;

    ThObsF.addObsToFactory("AlsTest", boost::factory<AlsTest*>());
    
    /* register user-defined model named ModelName defined in class ModelClass using the following syntax: */
    /* ModelF.addModelToFactory(ModelName, boost::factory<ModelClass*>() ) */

    /* register user-defined ThObservable named ThObsName defined in class ThObsClass using the following syntax: */
    /* ThObsF.addObsToFactory(ThObsName, boost::factory<ThObsClass*>() )*/

    /* Create an object of the class ComputeObservables. */
    ComputeObservables CO(ModelF, ThObsF, ModelConf);
    
    std::map<std::string, std::string> DFlags;
    // DFlags["FLAG"] = "TRUE";
    CO.setFlags(DFlags);

    /* Get the map of observables if necessary. */
    std::map<std::string, double> DObs = CO.getObservables();

    /* Get the map of observables with the parameter values defined above. */
    DObs = CO.compute(DPars);


//    std::cout << "\nObservables[" << i + 1 << "]:" << std::endl;
    for (std::map<std::string, double>::iterator it = DObs.begin(); it != DObs.end(); it++) {
        std::cout << it->first << " = " << it->second << std::endl;
    }

//
//    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}