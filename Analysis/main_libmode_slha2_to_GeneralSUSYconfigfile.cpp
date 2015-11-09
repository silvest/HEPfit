/*
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

/**
 * @example libmode_config.cpp
 * This is an example of how to compute observables from the input parameters
 * defined in a model configuration file.
 *
 */

#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <iterator>

#include "SLHA.hpp"     // From LHPC, see https://lhpc.hepforge.org

#include <iostream>
#include <ComputeObservables.h>


std::string dtostring(double );
void slhaparsing( LHPC::SlhaTwoWithSpheno& );

int main(int argc, char** argv)
{
    try {
        
        if(argc != 2){
            std::cout << "\nThe name of 1 SLHA files to read is required\n" << std::endl;
            return EXIT_SUCCESS;
        }
        
	// make the configuration file from the SLHA2 file 
	std::string SLHAFileName( argv[1] );
	LHPC::SlhaParser ParserForSlhaBlockSet( true, false );
	LHPC::SlhaTwoWithSpheno SlhaSet( ParserForSlhaBlockSet, false );

	ParserForSlhaBlockSet.readFile( SLHAFileName );
	slhaparsing( SlhaSet );

	/* Define the model configuration file.                        */
        /* The model configuration file provides the default values of */
        /* the mandatory model parameters.                             */
        std::string ModelConf = "GeneralSUSY_SLHAparsed.conf"; 
        
        /* Define a map for the parameters to be varied. */
        std::map<std::string, double> DPars;
        
        /* Create objects of the classes ModelFactory and ThObsFactory */
        ModelFactory ModelF;
        ThObsFactory ThObsF;

        /* Create an object of the class ComputeObservables. */
        ComputeObservables CO(ModelF, ThObsF, ModelConf);
        
        /* Set the flags for the model being used, if necessary.                         */
        /* The flags have to correspond to the model specified in the model config file. */
        std::map<std::string, std::string> DFlags;
        // DFlags["FLAG"] = "TRUE";
        CO.setFlags(DFlags);
        
        /* Get the map of observables if necessary. */
        std::map<std::string, double> DObs = CO.getObservables();
        
	DPars = CO.getParameters();     // get the list of the parameters, READ from the config file
	std::cout << "\nParameters:"<< std::endl;
	for (std::map<std::string, double>::iterator it = DPars.begin(); it != DPars.end(); it++) {
	  std::cout << it->first << " = " << it->second << std::endl;
	}


	DObs = CO.compute(DPars);
	std::cout << "\nObservables:" << std::endl;
	for (std::map<std::string, double>::iterator it = DObs.begin(); it != DObs.end(); it++) {
	  std::cout << it->first << " = " << it->second << std::endl;
	}
        
        return EXIT_SUCCESS;
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}


void slhaparsing( LHPC::SlhaTwoWithSpheno& SlhaSet){

  std::vector<std::string> v;
  std::ifstream in("GeneralSUSY.conf");
  std::ofstream out("GeneralSUSY_SLHAparsed.conf");
  std::string s, line;
  std::string word;

  while( getline(in, line) ){           

    std::istringstream iss(line);
    std::vector<std::string> words;
    copy(std::istream_iterator<std::string>(iss),
	 std::istream_iterator<std::string>(),
	 back_inserter(words)
	 );

    if( words.size() >= 3 && words[0] == "ModelParameter" ){
      if(words[1] == "m1r") words[2] = dtostring( SlhaSet.EXTPAR( 1 ) );
      if(words[1] == "m1i") words[2] = dtostring( SlhaSet.IMEXTPAR( 1 ) );
      if(words[1] == "m2r") words[2] = dtostring( SlhaSet.EXTPAR( 2 ) );
      if(words[1] == "m2i") words[2] = dtostring( SlhaSet.IMEXTPAR( 2 ) );      
      if(words[1] == "m3") words[2] = dtostring( SlhaSet.EXTPAR( 3 ) );
      if(words[1] == "muHr") words[2] = dtostring( SlhaSet.EXTPAR( 23 ) );
      if(words[1] == "muHi") words[2] = dtostring( SlhaSet.IMEXTPAR( 23 ) );      
      //      if(words[1] == "mHptree") words[2] = dtostring( SlhaSet.EXTPAR( 27 ) );
      // using MA instead of mHptree despite the naming miss-match
      //      if(words[1] == "mHptree") words[2] = dtostring( SlhaSet.EXTPAR( 24 ) );
      if(words[1] == "mHptree") words[2] = dtostring( SlhaSet.EXTPAR( 26 ) );
      if(words[1] == "tanb") words[2] = dtostring( SlhaSet.EXTPAR( 25 ) );
      if(words[1] == "Q") words[2] = dtostring( SlhaSet.EXTPAR( 0 ) );


      if(words[1] == "msQhat2_11r") words[2] =  dtostring( SlhaSet.MSQ2IN( 1, 1 ) );
      if(words[1] == "msQhat2_22r") words[2] =  dtostring( SlhaSet.MSQ2IN( 2, 2 ) );
      if(words[1] == "msQhat2_33r") words[2] =  dtostring( SlhaSet.MSQ2IN( 3, 3 ) );
      if(words[1] == "msQhat2_12r") words[2] =  dtostring( SlhaSet.MSQ2IN( 1, 2 ) );
      if(words[1] == "msQhat2_13r") words[2] =  dtostring( SlhaSet.MSQ2IN( 1, 3 ) );
      if(words[1] == "msQhat2_23r") words[2] =  dtostring( SlhaSet.MSQ2IN( 2, 3 ) );
      if(words[1] == "msQhat2_12i") words[2] =  dtostring( SlhaSet.IMMSQ2IN( 1, 2 ) );
      if(words[1] == "msQhat2_13i") words[2] =  dtostring( SlhaSet.IMMSQ2IN( 1, 3 ) );
      if(words[1] == "msQhat2_23i") words[2] =  dtostring( SlhaSet.IMMSQ2IN( 2, 3 ) );


      if(words[1] == "msUhat2_11r") words[2] =  dtostring( SlhaSet.MSU2IN( 1, 1 ) );
      if(words[1] == "msUhat2_22r") words[2] =  dtostring( SlhaSet.MSU2IN( 2, 2 ) );
      if(words[1] == "msUhat2_33r") words[2] =  dtostring( SlhaSet.MSU2IN( 3, 3 ) );
      if(words[1] == "msUhat2_12r") words[2] =  dtostring( SlhaSet.MSU2IN( 1, 2 ) );
      if(words[1] == "msUhat2_13r") words[2] =  dtostring( SlhaSet.MSU2IN( 1, 3 ) );
      if(words[1] == "msUhat2_23r") words[2] =  dtostring( SlhaSet.MSU2IN( 2, 3 ) );
      if(words[1] == "msUhat2_12i") words[2] =  dtostring( SlhaSet.IMMSU2IN( 1, 2 ) );
      if(words[1] == "msUhat2_13i") words[2] =  dtostring( SlhaSet.IMMSU2IN( 1, 3 ) );
      if(words[1] == "msUhat2_23i") words[2] =  dtostring( SlhaSet.IMMSU2IN( 2, 3 ) );


      if(words[1] == "msDhat2_11r") words[2] =  dtostring( SlhaSet.MSD2IN( 1, 1 ) );
      if(words[1] == "msDhat2_22r") words[2] =  dtostring( SlhaSet.MSD2IN( 2, 2 ) );
      if(words[1] == "msDhat2_33r") words[2] =  dtostring( SlhaSet.MSD2IN( 3, 3 ) );
      if(words[1] == "msDhat2_12r") words[2] =  dtostring( SlhaSet.MSD2IN( 1, 2 ) );
      if(words[1] == "msDhat2_13r") words[2] =  dtostring( SlhaSet.MSD2IN( 1, 3 ) );
      if(words[1] == "msDhat2_23r") words[2] =  dtostring( SlhaSet.MSD2IN( 2, 3 ) );
      if(words[1] == "msDhat2_12i") words[2] =  dtostring( SlhaSet.IMMSD2IN( 1, 2 ) );
      if(words[1] == "msDhat2_13i") words[2] =  dtostring( SlhaSet.IMMSD2IN( 1, 3 ) );
      if(words[1] == "msDhat2_23i") words[2] =  dtostring( SlhaSet.IMMSD2IN( 2, 3 ) );


      if(words[1] == "msLhat2_11r") words[2] =  dtostring( SlhaSet.MSL2IN( 1, 1 ) );
      if(words[1] == "msLhat2_22r") words[2] =  dtostring( SlhaSet.MSL2IN( 2, 2 ) );
      if(words[1] == "msLhat2_33r") words[2] =  dtostring( SlhaSet.MSL2IN( 3, 3 ) );
      if(words[1] == "msLhat2_12r") words[2] =  dtostring( SlhaSet.MSL2IN( 1, 2 ) );
      if(words[1] == "msLhat2_13r") words[2] =  dtostring( SlhaSet.MSL2IN( 1, 3 ) );
      if(words[1] == "msLhat2_23r") words[2] =  dtostring( SlhaSet.MSL2IN( 2, 3 ) );
      if(words[1] == "msLhat2_12i") words[2] =  dtostring( SlhaSet.IMMSL2IN( 1, 2 ) );
      if(words[1] == "msLhat2_13i") words[2] =  dtostring( SlhaSet.IMMSL2IN( 1, 3 ) );
      if(words[1] == "msLhat2_23i") words[2] =  dtostring( SlhaSet.IMMSL2IN( 2, 3 ) );


      if(words[1] == "msEhat2_11r") words[2] =  dtostring( SlhaSet.MSE2IN( 1, 1 ) );
      if(words[1] == "msEhat2_22r") words[2] =  dtostring( SlhaSet.MSE2IN( 2, 2 ) );
      if(words[1] == "msEhat2_33r") words[2] =  dtostring( SlhaSet.MSE2IN( 3, 3 ) );
      if(words[1] == "msEhat2_12r") words[2] =  dtostring( SlhaSet.MSE2IN( 1, 2 ) );
      if(words[1] == "msEhat2_13r") words[2] =  dtostring( SlhaSet.MSE2IN( 1, 3 ) );
      if(words[1] == "msEhat2_23r") words[2] =  dtostring( SlhaSet.MSE2IN( 2, 3 ) );
      if(words[1] == "msEhat2_12i") words[2] =  dtostring( SlhaSet.IMMSE2IN( 1, 2 ) );
      if(words[1] == "msEhat2_13i") words[2] =  dtostring( SlhaSet.IMMSE2IN( 1, 3 ) );
      if(words[1] == "msEhat2_23i") words[2] =  dtostring( SlhaSet.IMMSE2IN( 2, 3 ) );


      if(words[1] == "TUhat_11r") words[2] =  dtostring( SlhaSet.TUIN( 1, 1 ) );
      if(words[1] == "TUhat_22r") words[2] =  dtostring( SlhaSet.TUIN( 2, 2 ) );
      if(words[1] == "TUhat_12r") words[2] =  dtostring( SlhaSet.TUIN( 1, 2 ) );
      if(words[1] == "TUhat_33r") words[2] =  dtostring( SlhaSet.TUIN( 3, 3 ) );
      if(words[1] == "TUhat_13r") words[2] =  dtostring( SlhaSet.TUIN( 1, 3 ) );
      if(words[1] == "TUhat_21r") words[2] =  dtostring( SlhaSet.TUIN( 2, 1 ) );
      if(words[1] == "TUhat_23r") words[2] =  dtostring( SlhaSet.TUIN( 2, 3 ) );
      if(words[1] == "TUhat_31r") words[2] =  dtostring( SlhaSet.TUIN( 3, 1 ) );
      if(words[1] == "TUhat_32r") words[2] =  dtostring( SlhaSet.TUIN( 3, 2 ) );
      if(words[1] == "TUhat_11i") words[2] =  dtostring( SlhaSet.IMTUIN( 1, 1 ) );
      if(words[1] == "TUhat_22i") words[2] =  dtostring( SlhaSet.IMTUIN( 2, 2 ) );
      if(words[1] == "TUhat_12i") words[2] =  dtostring( SlhaSet.IMTUIN( 1, 2 ) );
      if(words[1] == "TUhat_33i") words[2] =  dtostring( SlhaSet.IMTUIN( 3, 3 ) );
      if(words[1] == "TUhat_13i") words[2] =  dtostring( SlhaSet.IMTUIN( 1, 3 ) );
      if(words[1] == "TUhat_21i") words[2] =  dtostring( SlhaSet.IMTUIN( 2, 1 ) );
      if(words[1] == "TUhat_23i") words[2] =  dtostring( SlhaSet.IMTUIN( 2, 3 ) );
      if(words[1] == "TUhat_31i") words[2] =  dtostring( SlhaSet.IMTUIN( 3, 1 ) );
      if(words[1] == "TUhat_32i") words[2] =  dtostring( SlhaSet.IMTUIN( 3, 2 ) );


      if(words[1] == "TDhat_11r") words[2] =  dtostring( SlhaSet.TDIN( 1, 1 ) );
      if(words[1] == "TDhat_22r") words[2] =  dtostring( SlhaSet.TDIN( 2, 2 ) );
      if(words[1] == "TDhat_12r") words[2] =  dtostring( SlhaSet.TDIN( 1, 2 ) );
      if(words[1] == "TDhat_33r") words[2] =  dtostring( SlhaSet.TDIN( 3, 3 ) );
      if(words[1] == "TDhat_13r") words[2] =  dtostring( SlhaSet.TDIN( 1, 3 ) );
      if(words[1] == "TDhat_21r") words[2] =  dtostring( SlhaSet.TDIN( 2, 1 ) );
      if(words[1] == "TDhat_23r") words[2] =  dtostring( SlhaSet.TDIN( 2, 3 ) );
      if(words[1] == "TDhat_31r") words[2] =  dtostring( SlhaSet.TDIN( 3, 1 ) );
      if(words[1] == "TDhat_32r") words[2] =  dtostring( SlhaSet.TDIN( 3, 2 ) );
      if(words[1] == "TDhat_11i") words[2] =  dtostring( SlhaSet.IMTDIN( 1, 1 ) );
      if(words[1] == "TDhat_22i") words[2] =  dtostring( SlhaSet.IMTDIN( 2, 2 ) );
      if(words[1] == "TDhat_12i") words[2] =  dtostring( SlhaSet.IMTDIN( 1, 2 ) );
      if(words[1] == "TDhat_33i") words[2] =  dtostring( SlhaSet.IMTDIN( 3, 3 ) );
      if(words[1] == "TDhat_13i") words[2] =  dtostring( SlhaSet.IMTDIN( 1, 3 ) );
      if(words[1] == "TDhat_21i") words[2] =  dtostring( SlhaSet.IMTDIN( 2, 1 ) );
      if(words[1] == "TDhat_23i") words[2] =  dtostring( SlhaSet.IMTDIN( 2, 3 ) );
      if(words[1] == "TDhat_31i") words[2] =  dtostring( SlhaSet.IMTDIN( 3, 1 ) );
      if(words[1] == "TDhat_32i") words[2] =  dtostring( SlhaSet.IMTDIN( 3, 2 ) );


      if(words[1] == "TEhat_11r") words[2] =  dtostring( SlhaSet.TEIN( 1, 1 ) );
      if(words[1] == "TEhat_22r") words[2] =  dtostring( SlhaSet.TEIN( 2, 2 ) );
      if(words[1] == "TEhat_12r") words[2] =  dtostring( SlhaSet.TEIN( 1, 2 ) );
      if(words[1] == "TEhat_33r") words[2] =  dtostring( SlhaSet.TEIN( 3, 3 ) );
      if(words[1] == "TEhat_13r") words[2] =  dtostring( SlhaSet.TEIN( 1, 3 ) );
      if(words[1] == "TEhat_21r") words[2] =  dtostring( SlhaSet.TEIN( 2, 1 ) );
      if(words[1] == "TEhat_23r") words[2] =  dtostring( SlhaSet.TEIN( 2, 3 ) );
      if(words[1] == "TEhat_31r") words[2] =  dtostring( SlhaSet.TEIN( 3, 1 ) );
      if(words[1] == "TEhat_32r") words[2] =  dtostring( SlhaSet.TEIN( 3, 2 ) );
      if(words[1] == "TEhat_11i") words[2] =  dtostring( SlhaSet.IMTEIN( 1, 1 ) );
      if(words[1] == "TEhat_22i") words[2] =  dtostring( SlhaSet.IMTEIN( 2, 2 ) );
      if(words[1] == "TEhat_12i") words[2] =  dtostring( SlhaSet.IMTEIN( 1, 2 ) );
      if(words[1] == "TEhat_33i") words[2] =  dtostring( SlhaSet.IMTEIN( 3, 3 ) );
      if(words[1] == "TEhat_13i") words[2] =  dtostring( SlhaSet.IMTEIN( 1, 3 ) );
      if(words[1] == "TEhat_21i") words[2] =  dtostring( SlhaSet.IMTEIN( 2, 1 ) );
      if(words[1] == "TEhat_23i") words[2] =  dtostring( SlhaSet.IMTEIN( 2, 3 ) );
      if(words[1] == "TEhat_31i") words[2] =  dtostring( SlhaSet.IMTEIN( 3, 1 ) );
      if(words[1] == "TEhat_32i") words[2] =  dtostring( SlhaSet.IMTEIN( 3, 2 ) );


    }

    std::ostringstream oss;
    copy(words.begin(), words.end(),
	 std::ostream_iterator<std::string>(oss, " ")   
	 );

    line=oss.str();    // the modified line
    s += line + "\n";      
  }    
  out << s << std::endl; // write the file
}


std::string dtostring(double xd){
  std::ostringstream convert;
  if ( isnan(xd) != 0 ) xd = 0.0;
  convert << std::scientific << xd;	
  return convert.str();
}
