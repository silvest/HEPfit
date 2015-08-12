/*
 * FlhaOne.cpp
 *
 *  Created on: Feb 6, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      REFBAGME.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "SLHA.hpp"

namespace LHPC
{
  FlhaOne::FlhaOne( SlhaParser& fileParser,
                    bool const isVerbose ) :
      FCINFO( "FCINFO",
              isVerbose ),
      MODSEL( "MODSEL",
              BOL::UsefulStuff::notANumber,
              isVerbose ),
      FMODSEL( "FMODSEL",
               BOL::UsefulStuff::notANumber,
               isVerbose ),
      SMINPUTS( "SMINPUTS",
                BOL::UsefulStuff::notANumber,
                isVerbose ),
      VCKMIN( "VCKMIN",
              0.0,
              isVerbose ),
      UPMNSIN( "UPMNSIN",
               0.0,
               isVerbose ),
      VCKM( "VCKM",
            0.0,
            isVerbose ),
      IMVCKM( "IMVCKM",
            0.0,
            isVerbose ),
      UPMNS( "UPMNS",
            0.0,
            isVerbose ),
      IMUPMNS( "IMUPMNS",
               0.0,
               isVerbose ),
      FMASS( "FMASS",
             RunningConstant(),
             isVerbose ),
      FMASSERR( "FMASSERR",
                RunningConstantError(),
                isVerbose ),
      FLIFE( "FLIFE",
             BOL::UsefulStuff::notANumber,
             isVerbose ),
      FLIFEERR( "FLIFEERR",
                std::pair< double, double >( BOL::UsefulStuff::notANumber,
                                             BOL::UsefulStuff::notANumber ),
                isVerbose ),
      FCONST( "FCONST",
              BOL::UsefulStuff::notANumber,
              isVerbose ),
      FCONSTERR( "FCONSTERR",
                 std::pair< double, double >( BOL::UsefulStuff::notANumber,
                                              BOL::UsefulStuff::notANumber ),
                 isVerbose ),
      FCONSTRATIO( "FCONSTRATIO",
                   BOL::UsefulStuff::notANumber,
             isVerbose ),
      FCONSTRATIOERR( "FCONSTRATIOERR",
                     std::pair< double, double >( BOL::UsefulStuff::notANumber,
                                                BOL::UsefulStuff::notANumber ),
                      isVerbose ),
      FBAG( "FBAG",
            RunningConstant(),
            isVerbose ),
      FBAGERR( "FBAGERR",
               RunningConstantError(),
               isVerbose ),
      FWCOEF( "FWCOEF",
               BOL::UsefulStuff::notANumber,
          isVerbose ),
      FWCOEFERR( "FWCOEFERR",
                  std::pair< double, double >( BOL::UsefulStuff::notANumber,
                                               BOL::UsefulStuff::notANumber ),
                  isVerbose ),
      IMFWCOEF( "IMFWCOEF",
                 BOL::UsefulStuff::notANumber,
                 isVerbose ),
      IMFWCOEFERR( "IMFWCOEFERR",
                    std::pair< double, double >( BOL::UsefulStuff::notANumber,
                                                BOL::UsefulStuff::notANumber ),
                    isVerbose ),
      FOBS( "FOBS",
            FlavorObservable(),
            isVerbose ),
      FOBSERR( "FOBSERR",
               FlavorObservableError(),
               isVerbose ),
      FOBSSM( "FOBSSM",
              FlavorObservable(),
              isVerbose ),
      FOBSSMERR( "FOBSSMERR",
                 FlavorObservableError(),
                 isVerbose ),
      FPARAM( "FPARAM",
              FlavorObservable(),
              isVerbose ),
      FPARAMERR( "FPARAMERR",
                 FlavorObservableError(),
                 isVerbose ),
      fileParser( fileParser ),
      isVerbose( isVerbose )
  {
    fileParser.registerBlock( FCINFO );
    fileParser.registerBlock( MODSEL );
    fileParser.registerBlock( FMODSEL );
    fileParser.registerBlock( SMINPUTS );
    fileParser.registerBlock( VCKMIN );
    fileParser.registerBlock( UPMNSIN );
    fileParser.registerBlock( VCKM );
    fileParser.registerBlock( IMVCKM );
    fileParser.registerBlock( UPMNS );
    fileParser.registerBlock( IMUPMNS );
    fileParser.registerBlock( FMASS );
    fileParser.registerBlock( FMASSERR );
    fileParser.registerBlock( FLIFE );
    fileParser.registerBlock( FLIFEERR );
    fileParser.registerBlock( FCONST );
    fileParser.registerBlock( FCONSTERR );
    fileParser.registerBlock( FCONSTRATIO );
    fileParser.registerBlock( FCONSTRATIOERR );
    fileParser.registerBlock( FBAG );
    fileParser.registerBlock( FBAGERR );
    fileParser.registerBlock( FWCOEF );
    fileParser.registerBlock( FWCOEFERR );
    fileParser.registerBlock( IMFWCOEF );
    fileParser.registerBlock( IMFWCOEFERR );
    fileParser.registerBlock( FOBS );
    fileParser.registerBlock( FOBSERR );
    fileParser.registerBlock( FOBSSM );
    fileParser.registerBlock( FOBSSMERR );
    fileParser.registerBlock( FPARAM );
    fileParser.registerBlock( FPARAMERR );
  }

  FlhaOne::~FlhaOne()
  {
    // does nothing.
  }

}
