// Include the necessary header file.
#include <HEPfit.h>

// Define the model configuration file.
std::string ModelConf = "SomeModel.conf";

// Define a map for the observables.
std::map<std::string, double> DObs;

// Define a map for the parameters to be varied.
std::map<std::string, double> DPars;

// Initialize the observables to be returned.
DObs["Mw"] = 0;
DObs["GammaZ"] = 0.;
DObs["AFBbottom"] = 0.;

// Create and object of the class ComputeObservables.
ComputeObservables CO(ModelConf, DObs);

// Vary the parameters that need to be varied in the analysis.
DPars["Mz"] = 91.1875;
DPars["AlsMz"] = 0.1184;

// Get the map of observables with the parameter values defined above.
DObs = CO.compute(DPars);