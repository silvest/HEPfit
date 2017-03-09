#include <mpi.h>
#include <math.h>
#include <iostream>
#include <string>
#include "interfaces.hpp"
#include <HEPfit.h>
#include <InputParameters.h>

#define NINPUT  2
#define NDERIV  1

double loglikelihood (double theta[], int nDims, double phi[], int nDerived);
void prior (double cube[], double theta[], int nDims);
double hepfit_call(double *theta, int nDims, double *phi, int nDerived, int index);

int main(){
  
  Settings settings;

  settings.nDims         = 2; 
  settings.nDerived      = 1;
  settings.nlive         = 500;
  settings.num_repeats   = settings.nDims*5;
  settings.do_clustering = false;

  settings.precision_criterion = 1e-3;

  settings.base_dir      = "chains";
  settings.file_root     = "h"; 

  settings.write_resume  = true;
  settings.read_resume   = true;
  settings.write_live    = true;
  settings.write_dead    = true;
  settings.write_stats   = true;

  settings.equals        = true;
  settings.posteriors    = true;
  settings.cluster_posteriors = false;

  settings.feedback      = 1;
  settings.update_files  = settings.nlive;

  settings.boost_posterior= 5.0;

  run_polychord(loglikelihood, prior, settings) ;
}

double loglikelihood (double theta[], int nDims, double phi[], int nDerived){
  // theta[]  input parameters
  // nDims    number of input parameters
  // phi[]    Output derivedParameters
  // nDerived number of output derived parameters

  //  double logZero = -DBL_MAX;
  
  int index = MPI::COMM_WORLD.Get_rank();
  //  std::cout << "index = " << index << std::endl;
  
  return hepfit_call(theta, nDims, phi, nDerived, index); 
}


double hepfit_call(double *theta, int nDims, double *phi, int nDerived, int index){

  std::string selected_pars[nDims] = {"mtop","AlsMz"};   
  std::string selected_obsvs[nDerived] = {"Mw"}; 

  std::string ModelName = "StandardModel";
        
  /* Create an object of the class InputParameters. */
  InputParameters IP;
        
  /* Read a map for the mandatory model parameters. (Default values in InputParameters.h) */ 
  std::map<std::string, double> DPars_IN = IP.getInputParameters(ModelName);
        
  /* Create objects of the classes ModelFactory and ThObsFactory */
  ModelFactory ModelF;
  ThObsFactory ThObsF;

  /* Create an object of the class ComputeObservables. */
  ComputeObservables CO(ModelF, ThObsF, ModelName, DPars_IN);
        
  /* Add the observables to be returned. */
  CO.AddObservable( selected_obsvs[0] );
        
  /* Set the flags for the model being used, if necessary. */
  // std::map<std::string, std::string> DFlags;
  // DFlags["epsilon2SM"] = "TRUE";
  // CO.setFlags(DFlags);
        
  /* Get the map of observables */
  std::map<std::string, double> DObs = CO.getObservables();
        
  /* Define a map for the parameters to be varied. */
  std::map<std::string, double> DPars;
        
  /* Vary the parameters that need to be varied in the analysis. */
  // to be obtained from sampler: PolyChord or MultiNest
  for(int i = 0; i < nDims; i++){
    DPars[ selected_pars[i] ] = theta[i];
    std::cout << "DPars[ "<< selected_pars[i] << " ] = " << DPars[ selected_pars[i] ];
  }

  /* Get the map of observables with the parameter values defined above. */
  DObs = CO.compute(DPars);

  for (int i = 0; i < nDerived; i++){
    phi[i] = DObs[ selected_obsvs[i] ];
    std::cout << "DObs[ "<< selected_obsvs[i] << " ] = " << DObs[ selected_obsvs[i] ] << std::endl;
  }

  //@@@@ HEPfit should have option for return log-likelihood ??

  // compute and retrun log likelihood
  // Observable: Mw = 80.3613 +- 0.03, say
  double Mw = 80.3613;
  double Mw_err = 0.03;
  
  double MwPredicted = phi[0];

  double chi = (MwPredicted - Mw) / ( sqrt(2.) * Mw_err );

  std::cout << "log-likelihood = " << - chi * chi - log( Mw_err * sqrt(2.0 * M_PI) );
  
  return - chi * chi - log( Mw_err * sqrt(2.0 * M_PI) );
}
