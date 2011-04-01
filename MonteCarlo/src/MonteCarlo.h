/* 
 * File:   MonteCarlo.h
 * Author: silvest
 *
 * Created on March 8, 2011, 3:18 PM
 */

#ifndef __MONTECARLO__H
#define __MONTECARLO__H

#include <BCModel.h>
#include <BCH1D.h>
#include <TH1D.h>
#include <TFile.h>
#include <Observable.h>
#include <ModelParameter.h>
#include <Model.h>
#include <map>

#define NBINS1D 100

// This is a MonteCarlo header file.
// Model source code is located in file MonteCarlo/MonteCarlo.cxx

// ---------------------------------------------------------
class MonteCarlo : public BCModel {
public:

    // Constructors and destructor
    MonteCarlo(const char * name, Model * Mod, std::vector<ModelParameter> ModPars, std::vector<Observable> Obs);
    ~MonteCarlo();

    // Methods to overload, see file MonteCarlo.cxx
    void DefineParameters();
    double LogLikelihood(std::vector <double> parameters);
    void PrintHistogram();
    //void MCMCIterationInterface();
    void MCMCCurrentPointInterface(std::vector <double> & point, int ichain, bool accepted);
private:
    std::vector<ModelParameter> ModPars;
    std::vector<Observable> Obs_MCMC, Obs_noMCMC;
    Model * Mod;
    std::map<std::string, double> DPars;
    std::map<std::string, BCH1D * > Histo1D;
    std::map<std::string, TH1D * > InHisto1D;
};
// ---------------------------------------------------------

#endif
