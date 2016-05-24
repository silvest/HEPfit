// Include the necessary header file.
#include <HEPfit.h>

// Define a map for the observables.
std::map<std::string, double> DObs;

// Define a map for the mandatory model parameters used during initializing a model.
std::map<std::string, double> DPars_IN;

// Define a map for the parameters to be varied.
std::map<std::string, double> DPars;

// Define a map for the model flags.
std::map<std::string, std::string> DFlags;

// Define the name of the model to be used.
std::string ModelName = "NPZbbbar";

// Create and object of the class InputParameters.
InputParameters IP;

// Read a map for the mandatory model parameters. (Default values in InputParameters.h)
DPars_IN = IP.getInputParameters(ModelName);

// Change the default values of the mandatory model parameters if necessary.
// This can also be done with Dpars after creating an object of ComputeObservables
DPars_IN["mcharm"] = 1.3;
DPars_IN["mub"] = 4.2;

// Initialize the observables to be returned.
DObs["Mw"] = 0;
DObs["GammaZ"] = 0.;
DObs["AFBbottom"] = 0.;

// Initialize the model flags to be set.
DFlags["NPZbbbarLR"] = "TRUE";

// Create and object of the class ComputeObservables.
ComputeObservables CO(ModelName, DPars_IN, DObs);

// Set the flags for the model being used, if necessary.
CO.setFlags(DFlags);

// Vary the parameters that need to be varied in the analysis.
DPars["mtop"] = 170.0;
DPars["mHl"] = 126.0;

// Get the map of observables with the parameter values defined above.
DObs = CO.compute(DPars);