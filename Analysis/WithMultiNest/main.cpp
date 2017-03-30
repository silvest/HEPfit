#include <mpi.h>
#include <string.h>
#include <iostream>

#include"multinest.h"

#include <HEPfit.h>
#include <InputParameters.h>

#define NINPUT  2
#define NDERIV  1

double spriorran[NINPUT][2];

extern"C" {
  double priors_mp_uniformprior_(double *r, double *x1, double *x2);
  double priors_mp_logprior_(double *r, double *x1, double *x2); 
}

void setpriors(void);
void slikelihood(double *theta, int &ndim, int &npars, double &slhood, void *id);

void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &logZerr, void *context);
double hepfit_call(double *theta, int nDims, double *phi, int nDerived, int index);

int main(){

  int i;

                        // set the MultiNest sampling parameters
  int IS = 0;					// do Nested Importance Sampling?
  int ndims = NINPUT;   // dimensionality (no. of free parameters)
  int nPar = NINPUT + NDERIV;       // total no. of parameters including free & derived parameters
  int nlive = 3500; //5623; //5196;     // number of live points
  int mmodal = 0;       // do mode separation? !whether to do multimodal sampling
  int ceff = 1;         // run in constant efficiency mode?   !whether to run in constant efficiency mode
  double efr = 0.8;     // set the required efficiency   !enlargement factor reduction parameter
  double tol = 0.5;     // tol, defines the stopping criteria     !evidence tolerance factor
  int nClsPar = 0;      // no. of parameters to do mode separation on !no.ofparameters which clustering shouldbe performed

  char root[100] = "chains/F-";// root for output files 
  int updInt = 8;       // after how many iterations feedback is required & the output files should be updated

  double Ztol = -1E90;  // all the modes with logZ < Ztol are ignored
                        // !null evidence (set it to very high negative no. if null evidence is unknown) 
  int maxModes = 10;    // expected max no. of modes (used only for memory allocation)
  int pWrap[ndims];     // which parameters to have periodic boundary conditions?
  for(i = 0; i < ndims; i++) pWrap[i] = 0;  // !parameters to wrap around (0 is F & non-zero T)
  int seed = -1;        // random no. generator seed, if < 0 then take the seed from system clock
  int fb = 1;           // need feedback on standard output?     !feedback on the sampling progress?
  int resume = 1;       // resume from a previous job? !whether resume a previous run   CHECK-POINTING, YES please == 1
  int outfile = 1;      // write output files?
  int initMPI = 1;      // initialize MPI routines?, relevant only if compiling with MPI
                        // set it to F if you want your main program to handle MPI initialization
  double logZero = -DBL_MAX;// points with loglike < logZero will be ignored by MultiNest
  int maxiter = 0;      // max no. of iterations,a value < 0 means infinity. MultiNest will terminate if either it 
                        // has done max no. iterations or convergence criterion (defined through tol) has been satisfied
  int context = 0;      // not required by MultiNest, any additional information user wants to pass

  setpriors();     

  // calling MultiNest
  nested::run(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, slikelihood, dumper, (void *)context); 

  return 0;
}

void slikelihood(double *theta, int &ndim, int &npars, double &slhood, void *id){
  //  double logZero = -DBL_MAX;
  
  double *phi = theta + ndim;
  int nDerived = NDERIV;
  int *index = (int *)id;
  *index = MPI::COMM_WORLD.Get_rank();
  //  std::cout << "index = " << *index << std::endl;

  for(int i=0; i < ndim; i++){ 
    theta[i] = priors_mp_uniformprior_(&theta[i], &spriorran[i][0], &spriorran[i][1]);
  } 
  slhood = hepfit_call(theta, ndim, phi, nDerived, *index); 
}


double hepfit_call(double *theta, int nDims, double *phi, int nDerived, int index){

  std::string selected_pars[NINPUT] = {"mtop","AlsMz"};   
  std::string selected_obsvs[NDERIV] = {"Mw"}; 

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

void setpriors(void){
  // Parameter "mtop": [169.54, 177.14]
  spriorran[0][0] = 169.54; 
  spriorran[0][1] = 177.14;

  // Parameter "AlsMz": [0.116, 0.121]
  spriorran[1][0] = 0.116;
  spriorran[1][1] = 0.121;
}

/************************************************* dumper routine ******************************************************/

// The dumper routine will be called every updInt*10 iterations
// MultiNest doesn not need to the user to do anything. User can use the arguments in whichever way he/she wants
//
//
// Arguments:
//
// nSamples 						= total number of samples in posterior distribution
// nlive 						= total number of live points
// nPar 						= total number of parameters (free + derived)
// physLive[1][nlive * (nPar + 1)] 			= 2D array containing the last set of live points (physical parameters plus derived parameters) along with their loglikelihood values
// posterior[1][nSamples * (nPar + 2)] 			= posterior distribution containing nSamples points. Each sample has nPar parameters (physical + derived) along with the their loglike value & posterior probability
// paramConstr[1][4*nPar]:
// paramConstr[0][0] to paramConstr[0][nPar - 1] 	= mean values of the parameters
// paramConstr[0][nPar] to paramConstr[0][2*nPar - 1] 	= standard deviation of the parameters
// paramConstr[0][nPar*2] to paramConstr[0][3*nPar - 1] = best-fit (maxlike) parameters
// paramConstr[0][nPar*4] to paramConstr[0][4*nPar - 1] = MAP (maximum-a-posteriori) parameters
// maxLogLike						= maximum loglikelihood value
// logZ							= log evidence value
// logZerr						= error on log evidence value
// context						void pointer, any additional information

void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &logZerr, void *context)
{
	// convert the 2D Fortran arrays to C arrays
	
	
	// the posterior distribution
	// postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns
	
	int i, j;
	
	double postdist[nSamples][nPar + 2];
	for( i = 0; i < nPar + 2; i++ )
		for( j = 0; j < nSamples; j++ )
			postdist[j][i] = posterior[0][i * nSamples + j];
	
	
	
	// last set of live points
	// pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column
	
	double pLivePts[nlive][nPar + 1];
	for( i = 0; i < nPar + 1; i++ )
		for( j = 0; j < nlive; j++ )
			pLivePts[j][i] = physLive[0][i * nlive + j];
}

/***********************************************************************************************************************/
