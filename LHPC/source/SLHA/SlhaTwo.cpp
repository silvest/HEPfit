/*
 * SlhaTwo.cpp
 *
 *  Created on: Feb 22, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "SLHA.hpp"

namespace LHPC
{
  SlhaTwo::SlhaTwo( SlhaParser& fileParser,
                    bool const isVerbose ) :
      SlhaOne( fileParser,
               isVerbose ),
      QEXTPAR( "QEXTPAR",
               BOL::UsefulStuff::notANumber,
               isVerbose ),
      IMMINPAR( "IMMINPAR",
                BOL::UsefulStuff::notANumber,
                isVerbose ),
      IMEXTPAR( "IMEXTPAR",
                BOL::UsefulStuff::notANumber,
                isVerbose ),
      IMNMIX( "IMNMIX",
              0.0,
              isVerbose ),
      IMUMIX( "IMUMIX",
              0.0,
              isVerbose ),
      IMVMIX( "IMVMIX",
              0.0,
              isVerbose ),
      IMSTOPMIX( "IMSTOPMIX",
                 0.0,
                 isVerbose ),
      IMSBOTMIX( "IMSBOTMIX",
                 0.0,
                 isVerbose ),
      IMSTAUMIX( "IMSTAUMIX",
                 0.0,
                 isVerbose ),
      IMALPHA( "IMALPHA",
               BOL::UsefulStuff::notANumber,
               isVerbose ),
      IMHMIX( "IMHMIX",
              BOL::UsefulStuff::notANumber,
              isVerbose ),
      IMGAUGE( "IMGAUGE",
                0.0,
                isVerbose ),
      IMMSOFT( "IMMSOFT",
               BOL::UsefulStuff::notANumber,
               isVerbose ),
      IMAU( "IMAU",
            0.0,
            isVerbose ),
      IMAD( "IMAD",
            0.0,
            isVerbose ),
      IMAE( "IMAE",
            0.0,
            isVerbose ),
      IMYU( "IMYU",
            0.0,
            isVerbose ),
      IMYD( "IMYD",
            0.0,
            isVerbose ),
      IMYE( "IMYE",
            0.0,
            isVerbose ),
      CVHMIX( "CVHMIX",
              0.0,
              isVerbose ),
      IMCVHMIX( "IMCVHMIX",
                0.0,
                isVerbose ),
      VCKMIN( "VCKMIN",
              0.0,
              isVerbose ),
      IMVCKMIN( "IMVCKMIN",
                0.0,
                isVerbose ),
      VCKM( "VCKM",
            0.0,
            isVerbose ),
      IMVCKM( "IMVCKM",
              0.0,
              isVerbose ),
      UPMNSIN( "UPMNSIN",
               0.0,
               isVerbose ),
      IMUPMNSIN( "IMUPMNSIN",
                 0.0,
                 isVerbose ),
      UPMNS( "UPMNS",
             0.0,
             isVerbose ),
      IMUPMNS( "IMUPMNS",
               0.0,
               isVerbose ),
      MSQ2IN( "MSQ2IN",
              0.0,
              isVerbose ),
      IMMSQ2IN( "IMMSQ2IN",
                0.0,
                isVerbose ),
      MSQ2( "MSQ2",
            0.0,
            isVerbose ),
      IMMSQ2( "IMMSQ2",
              0.0,
              isVerbose ),
      MSU2IN( "MSU2IN",
              0.0,
              isVerbose ),
      IMMSU2IN( "IMMSU2IN",
                0.0,
                isVerbose ),
      MSU2( "MSU2",
            0.0,
            isVerbose ),
      IMMSU2( "IMMSU2",
              0.0,
              isVerbose ),
      MSD2IN( "MSD2IN",
              0.0,
              isVerbose ),
      IMMSD2IN( "IMMSD2IN",
                0.0,
                isVerbose ),
      MSD2( "MSD2",
            0.0,
            isVerbose ),
      IMMSD2( "IMMSD2",
              0.0,
              isVerbose ),
      MSL2IN( "MSL2IN",
              0.0,
              isVerbose ),
      IMMSL2IN( "IMMSL2IN",
                0.0,
                isVerbose ),
      MSL2( "MSL2",
            0.0,
            isVerbose ),
      IMMSL2( "IMMSL2",
              0.0,
              isVerbose ),
      MSE2IN( "MSE2IN",
              0.0,
              isVerbose ),
      IMMSE2IN( "IMMSE2IN",
                0.0,
                isVerbose ),
      MSE2( "MSE2",
            0.0,
            isVerbose ),
      IMMSE2( "IMMSE2",
              0.0,
              isVerbose ),
      MSN2IN( "MSN2IN",
              0.0,
              isVerbose ),
      IMMSN2IN( "IMMSN2IN",
                0.0,
                isVerbose ),
      MSN2( "MSN2",
            0.0,
            isVerbose ),
      IMMSN2( "IMMSN2",
              0.0,
              isVerbose ),
      TUIN( "TUIN",
            0.0,
            isVerbose ),
      IMTUIN( "IMTUIN",
              0.0,
              isVerbose ),
      TU( "TU",
          0.0,
          isVerbose ),
      IMTU( "IMTU",
            0.0,
            isVerbose ),
      TDIN( "TDIN",
            0.0,
            isVerbose ),
      IMTDIN( "IMTDIN",
              0.0,
              isVerbose ),
      TD( "TD",
          0.0,
          isVerbose ),
      IMTD( "IMTD",
            0.0,
            isVerbose ),
      TEIN( "TEIN",
            0.0,
            isVerbose ),
      IMTEIN( "IMTEIN",
              0.0,
              isVerbose ),
      TE( "TE",
          0.0,
          isVerbose ),
      IMTE( "IMTE",
            0.0,
            isVerbose ),
      TNIN( "TNIN",
            0.0,
            isVerbose ),
      IMTNIN( "IMTNIN",
              0.0,
              isVerbose ),
      TN( "TN",
          0.0,
          isVerbose ),
      IMTN( "IMTN",
            0.0,
            isVerbose ),
      YN( "YN",
          0.0,
          isVerbose ),
      USQMIX( "USQMIX",
              0.0,
              isVerbose ),
      IMUSQMIX( "IMUSQMIX",
                0.0,
                isVerbose ),
      DSQMIX( "DSQMIX",
              0.0,
              isVerbose ),
      IMDSQMIX( "IMDSQMIX",
                0.0,
                isVerbose ),
      SELMIX( "SELMIX",
              0.0,
              isVerbose ),
      IMSELMIX( "IMSELMIX",
                0.0,
                isVerbose ),
      SNUMIX( "SNUMIX",
              0.0,
              isVerbose ),
      IMSNUMIX( "IMSNUMIX",
                0.0,
                isVerbose ),
      SNSMIX( "SNSMIX",
              0.0,
              isVerbose ),
      IMSNSMIX( "IMSNSMIX",
                0.0,
                isVerbose ),
      SNAMIX( "SNAMIX",
              0.0,
              isVerbose ),
      IMSNAMIX( "IMSNAMIX",
                0.0,
                isVerbose ),
      RVLAMLLEIN( "RVLAMLLEIN",
                  0.0,
                  isVerbose ),
      IMRVLAMLLEIN( "IMRVLAMLLEIN",
                    0.0,
                    isVerbose ),
      RVLAMLLE( "RVLAMLLE",
                0.0,
                isVerbose ),
      IMRVLAMLLE( "IMRVLAMLLE",
                  0.0,
                  isVerbose ),
      RVLAMLQDIN( "RVLAMLQDIN",
                  0.0,
                  isVerbose ),
      IMRVLAMLQDIN( "IMRVLAMLQDIN",
                    0.0,
                    isVerbose ),
      RVLAMLQD( "RVLAMLQD",
                0.0,
                isVerbose ),
      IMRVLAMLQD( "IMRVLAMLQD",
                  0.0,
                  isVerbose ),
      RVLAMUDDIN( "RVLAMUDDIN",
                  0.0,
                  isVerbose ),
      IMRVLAMUDDIN( "IMRVLAMUDDIN",
                    0.0,
                    isVerbose ),
      RVLAMUDD( "RVLAMUDD",
                0.0,
                isVerbose ),
      IMRVLAMUDD( "IMRVLAMUDD",
                  0.0,
                  isVerbose ),
      RVTLLEIN( "RVTLLEIN",
                0.0,
                isVerbose ),
      IMRVTLLEIN( "IMRVTLLEIN",
                  0.0,
                  isVerbose ),
      RVTLLE( "RVTLLE",
              0.0,
              isVerbose ),
      IMRVTLLE( "IMRVTLLE",
                0.0,
                isVerbose ),
      RVTLQDIN( "RVTLQDIN",
                0.0,
                isVerbose ),
      IMRVTLQDIN( "IMRVTLQDIN",
                  0.0,
                  isVerbose ),
      RVTLQD( "RVTLQD",
              0.0,
              isVerbose ),
      IMRVTLQD( "IMRVTLQD",
                0.0,
                isVerbose ),
      RVTUDDIN( "RVTUDDIN",
                0.0,
                isVerbose ),
      IMRVTUDDIN( "IMRVTUDDIN",
                  0.0,
                  isVerbose ),
      RVTUDD( "RVTUDD",
              0.0,
              isVerbose ),
      IMRVTUDD( "IMRVTUDD",
                0.0,
                isVerbose ),
      RVKAPPAIN( "RVKAPPAIN",
                 0.0,
                 isVerbose ),
      IMRVKAPPAIN( "IMRVKAPPAIN",
                   0.0,
                   isVerbose ),
      RVKAPPA( "RVKAPPA",
               0.0,
               isVerbose ),
      IMRVKAPPA( "IMRVKAPPA",
                 0.0,
                 isVerbose ),
      RVDIN( "RVDIN",
             0.0,
             isVerbose ),
      IMRVDIN( "IMRVDIN",
               0.0,
               isVerbose ),
      RVD( "RVD",
           0.0,
           isVerbose ),
      IMRVD( "IMRVD",
             0.0,
             isVerbose ),
      RVM2LH1( "RVM2LH1",
               0.0,
               isVerbose ),
      IMRVM2LH1( "IMRVM2LH1",
                 0.0,
                 isVerbose ),
      RVSNVEVIN( "RVSNVEVIN",
                 0.0,
                 isVerbose ),
      IMRVSNVEVIN( "IMRVSNVEVIN",
                   0.0,
                   isVerbose ),
      RVSNVEV( "RVSNVEV",
               0.0,
               isVerbose ),
      IMRVSNVEV( "IMRVSNVEV",
                 0.0,
                 isVerbose ),
      RVNMIX( "RVNMIX",
              0.0,
              isVerbose ),
      IMRVNMIX( "IMRVNMIX",
                0.0,
                isVerbose ),
      RVUMIX( "RVUMIX",
              0.0,
              isVerbose ),
      IMRVUMIX( "IMRVUMIX",
                0.0,
                isVerbose ),
      RVVMIX( "RVVMIX",
              0.0,
              isVerbose ),
      IMRVVMIX( "IMRVVMIX",
                0.0,
                isVerbose ),
      RVHMIX( "RVHMIX",
              0.0,
              isVerbose ),
      IMRVHMIX( "IMRVHMIX",
                0.0,
                isVerbose ),
      RVAMIX( "RVAMIX",
              0.0,
              isVerbose ),
      IMRVAMIX( "IMRVAMIX",
                0.0,
                isVerbose ),
      RVLMIX( "RVLMIX",
              0.0,
              isVerbose ),
      IMRVLMIX( "IMRVLMIX",
                0.0,
                isVerbose ),
      NMSSMRUN( "NMSSMRUN",
                BOL::UsefulStuff::notANumber,
                isVerbose ),
      NMHMIX( "NMHMIX",
              0.0,
              isVerbose ),
      IMNMHMIX( "IMNMHMIX",
                0.0,
                isVerbose ),
      NMAMIX( "NMAMIX",
              0.0,
              isVerbose ),
      IMNMAMIX( "IMNMAMIX",
                0.0,
                isVerbose ),
      NMNMIX( "NMNMIX",
              0.0,
              isVerbose ),
      IMNMNMIX( "IMNMNMIX",
            0.0,
            isVerbose )
  {
    fileParser.registerBlock( QEXTPAR );
    fileParser.registerBlock( IMMINPAR );
    fileParser.registerBlock( IMEXTPAR );
    fileParser.registerBlock( IMNMIX );
    fileParser.registerBlock( IMUMIX );
    fileParser.registerBlock( IMVMIX );
    fileParser.registerBlock( IMSTOPMIX );
    fileParser.registerBlock( IMSBOTMIX );
    fileParser.registerBlock( IMSTAUMIX );
    fileParser.registerBlock( IMALPHA );
    fileParser.registerBlock( IMHMIX );
    fileParser.registerBlock( IMGAUGE );
    fileParser.registerBlock( IMMSOFT );
    fileParser.registerBlock( IMAU );
    fileParser.registerBlock( IMAD );
    fileParser.registerBlock( IMAE );
    fileParser.registerBlock( IMYU );
    fileParser.registerBlock( IMYD );
    fileParser.registerBlock( IMYE );
    fileParser.registerBlock( CVHMIX );
    fileParser.registerBlock( IMCVHMIX );
    fileParser.registerBlock( VCKMIN );
    fileParser.registerBlock( IMVCKMIN );
    fileParser.registerBlock( VCKM );
    fileParser.registerBlock( IMVCKM );
    fileParser.registerBlock( UPMNSIN );
    fileParser.registerBlock( IMUPMNSIN );
    fileParser.registerBlock( UPMNS );
    fileParser.registerBlock( IMUPMNS );
    fileParser.registerBlock( MSQ2IN );
    fileParser.registerBlock( IMMSQ2IN );
    fileParser.registerBlock( MSQ2 );
    fileParser.registerBlock( IMMSQ2 );
    fileParser.registerBlock( MSU2IN );
    fileParser.registerBlock( IMMSU2IN );
    fileParser.registerBlock( MSU2 );
    fileParser.registerBlock( IMMSU2 );
    fileParser.registerBlock( MSD2IN );
    fileParser.registerBlock( IMMSD2IN );
    fileParser.registerBlock( MSD2 );
    fileParser.registerBlock( IMMSD2 );
    fileParser.registerBlock( MSL2IN );
    fileParser.registerBlock( IMMSL2IN );
    fileParser.registerBlock( MSL2 );
    fileParser.registerBlock( IMMSL2 );
    fileParser.registerBlock( MSE2IN );
    fileParser.registerBlock( IMMSE2IN );
    fileParser.registerBlock( MSE2 );
    fileParser.registerBlock( IMMSE2 );
    fileParser.registerBlock( MSN2IN );
    fileParser.registerBlock( IMMSN2IN );
    fileParser.registerBlock( MSN2 );
    fileParser.registerBlock( IMMSN2 );
    fileParser.registerBlock( TUIN );
    fileParser.registerBlock( IMTUIN );
    fileParser.registerBlock( TU );
    fileParser.registerBlock( IMTU );
    fileParser.registerBlock( TDIN );
    fileParser.registerBlock( IMTDIN );
    fileParser.registerBlock( TD );
    fileParser.registerBlock( IMTD );
    fileParser.registerBlock( TEIN );
    fileParser.registerBlock( IMTEIN );
    fileParser.registerBlock( TE );
    fileParser.registerBlock( IMTE );
    fileParser.registerBlock( TNIN );
    fileParser.registerBlock( IMTNIN );
    fileParser.registerBlock( TN );
    fileParser.registerBlock( IMTN );
    fileParser.registerBlock( YN );
    fileParser.registerBlock( USQMIX );
    fileParser.registerBlock( IMUSQMIX );
    fileParser.registerBlock( DSQMIX );
    fileParser.registerBlock( IMDSQMIX );
    fileParser.registerBlock( SELMIX );
    fileParser.registerBlock( IMSELMIX );
    fileParser.registerBlock( SNUMIX );
    fileParser.registerBlock( IMSNUMIX );
    fileParser.registerBlock( SNSMIX );
    fileParser.registerBlock( IMSNSMIX );
    fileParser.registerBlock( SNAMIX );
    fileParser.registerBlock( IMSNAMIX );
    fileParser.registerBlock( RVLAMLLEIN );
    fileParser.registerBlock( IMRVLAMLLEIN );
    fileParser.registerBlock( RVLAMLLE );
    fileParser.registerBlock( IMRVLAMLLE );
    fileParser.registerBlock( RVLAMLQDIN );
    fileParser.registerBlock( IMRVLAMLQDIN );
    fileParser.registerBlock( RVLAMLQD );
    fileParser.registerBlock( IMRVLAMLQD );
    fileParser.registerBlock( RVLAMUDDIN );
    fileParser.registerBlock( IMRVLAMUDDIN );
    fileParser.registerBlock( RVLAMUDD );
    fileParser.registerBlock( IMRVLAMUDD );
    fileParser.registerBlock( RVTLLEIN );
    fileParser.registerBlock( IMRVTLLEIN );
    fileParser.registerBlock( RVTLLE );
    fileParser.registerBlock( IMRVTLLE );
    fileParser.registerBlock( RVTLQDIN );
    fileParser.registerBlock( IMRVTLQDIN );
    fileParser.registerBlock( RVTLQD );
    fileParser.registerBlock( IMRVTLQD );
    fileParser.registerBlock( RVTUDDIN );
    fileParser.registerBlock( IMRVTUDDIN );
    fileParser.registerBlock( RVTUDD );
    fileParser.registerBlock( IMRVTUDD );
    fileParser.registerBlock( RVKAPPAIN );
    fileParser.registerBlock( IMRVKAPPAIN );
    fileParser.registerBlock( RVKAPPA );
    fileParser.registerBlock( IMRVKAPPA );
    fileParser.registerBlock( RVDIN );
    fileParser.registerBlock( IMRVDIN );
    fileParser.registerBlock( RVD );
    fileParser.registerBlock( IMRVD );
    fileParser.registerBlock( RVM2LH1 );
    fileParser.registerBlock( IMRVM2LH1 );
    fileParser.registerBlock( RVSNVEVIN );
    fileParser.registerBlock( IMRVSNVEVIN );
    fileParser.registerBlock( RVSNVEV );
    fileParser.registerBlock( IMRVSNVEV );
    fileParser.registerBlock( RVNMIX );
    fileParser.registerBlock( IMRVNMIX );
    fileParser.registerBlock( RVUMIX );
    fileParser.registerBlock( IMRVUMIX );
    fileParser.registerBlock( RVVMIX );
    fileParser.registerBlock( IMRVVMIX );
    fileParser.registerBlock( RVHMIX );
    fileParser.registerBlock( IMRVHMIX );
    fileParser.registerBlock( RVAMIX );
    fileParser.registerBlock( IMRVAMIX );
    fileParser.registerBlock( RVLMIX );
    fileParser.registerBlock( IMRVLMIX );
    fileParser.registerBlock( NMSSMRUN );
    fileParser.registerBlock( NMHMIX );
    fileParser.registerBlock( IMNMHMIX );
    fileParser.registerBlock( NMAMIX );
    fileParser.registerBlock( IMNMAMIX );
    fileParser.registerBlock( NMNMIX );
    fileParser.registerBlock( IMNMNMIX );
  }

  SlhaTwo::~SlhaTwo()
  {
    // does nothing.
  }

}
