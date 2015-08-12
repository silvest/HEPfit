/*
 * SlhaTwo.hpp
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

#ifndef SLHATWO_HPP_
#define SLHATWO_HPP_

#include "SlhaOne.hpp"

namespace LHPC
{
  /* this is just a collection of instances of classes derived from SlhaBlock
   * which covers all the blocks specified in SLHA2 (Comput. Phys. Commun. 180
   * (2009) 8 [arXiv:0801.0045 [hep-ph]]).
   */
  class SlhaTwo : public SlhaOne
  {
  public:
    SlhaTwo( SlhaParser& fileParser,
             bool const isVerbose = false );
    virtual
    ~SlhaTwo();

    SLHA::SparseSinglyIndexedBlock< double > QEXTPAR;
    SLHA::SparseSinglyIndexedBlock< double > IMMINPAR;
    SLHA::SparseSinglyIndexedBlock< double > IMEXTPAR;
    SLHA::DenseDoublyIndexedBlock< double > IMNMIX;
    SLHA::DenseDoublyIndexedBlock< double > IMUMIX;
    SLHA::DenseDoublyIndexedBlock< double > IMVMIX;
    SLHA::DenseDoublyIndexedBlock< double > IMSTOPMIX;
    SLHA::DenseDoublyIndexedBlock< double > IMSBOTMIX;
    SLHA::DenseDoublyIndexedBlock< double > IMSTAUMIX;
    SLHA::JustSingleValueBlock< double > IMALPHA;
    SLHA::SparseSinglyIndexedBlock< double > IMHMIX;
    SLHA::DenseSinglyIndexedBlock< double > IMGAUGE;
    SLHA::SparseSinglyIndexedBlock< double > IMMSOFT;
    SLHA::DenseDoublyIndexedBlock< double > IMAU;
    SLHA::DenseDoublyIndexedBlock< double > IMAD;
    SLHA::DenseDoublyIndexedBlock< double > IMAE;
    SLHA::DenseDoublyIndexedBlock< double > IMYU;
    SLHA::DenseDoublyIndexedBlock< double > IMYD;
    SLHA::DenseDoublyIndexedBlock< double > IMYE;
    SLHA::DenseDoublyIndexedBlock< double > CVHMIX;
    SLHA::DenseDoublyIndexedBlock< double > IMCVHMIX;
    SLHA::DenseSinglyIndexedBlock< double > VCKMIN;
    SLHA::DenseSinglyIndexedBlock< double > IMVCKMIN;
    SLHA::DenseSinglyIndexedBlock< double > VCKM;
    SLHA::DenseSinglyIndexedBlock< double > IMVCKM;
    SLHA::DenseSinglyIndexedBlock< double > UPMNSIN;
    SLHA::DenseSinglyIndexedBlock< double > IMUPMNSIN;
    SLHA::DenseSinglyIndexedBlock< double > UPMNS;
    SLHA::DenseSinglyIndexedBlock< double > IMUPMNS;
    SLHA::DenseDoublyIndexedBlock< double > MSQ2IN;
    SLHA::DenseDoublyIndexedBlock< double > IMMSQ2IN;
    SLHA::DenseDoublyIndexedBlock< double > MSQ2;
    SLHA::DenseDoublyIndexedBlock< double > IMMSQ2;
    SLHA::DenseDoublyIndexedBlock< double > MSU2IN;
    SLHA::DenseDoublyIndexedBlock< double > IMMSU2IN;
    SLHA::DenseDoublyIndexedBlock< double > MSU2;
    SLHA::DenseDoublyIndexedBlock< double > IMMSU2;
    SLHA::DenseDoublyIndexedBlock< double > MSD2IN;
    SLHA::DenseDoublyIndexedBlock< double > IMMSD2IN;
    SLHA::DenseDoublyIndexedBlock< double > MSD2;
    SLHA::DenseDoublyIndexedBlock< double > IMMSD2;
    SLHA::DenseDoublyIndexedBlock< double > MSL2IN;
    SLHA::DenseDoublyIndexedBlock< double > IMMSL2IN;
    SLHA::DenseDoublyIndexedBlock< double > MSL2;
    SLHA::DenseDoublyIndexedBlock< double > IMMSL2;
    SLHA::DenseDoublyIndexedBlock< double > MSE2IN;
    SLHA::DenseDoublyIndexedBlock< double > IMMSE2IN;
    SLHA::DenseDoublyIndexedBlock< double > MSE2;
    SLHA::DenseDoublyIndexedBlock< double > IMMSE2;
    SLHA::DenseDoublyIndexedBlock< double > MSN2IN;
    SLHA::DenseDoublyIndexedBlock< double > IMMSN2IN;
    SLHA::DenseDoublyIndexedBlock< double > MSN2;
    SLHA::DenseDoublyIndexedBlock< double > IMMSN2;
    SLHA::DenseDoublyIndexedBlock< double > TUIN;
    SLHA::DenseDoublyIndexedBlock< double > IMTUIN;
    SLHA::DenseDoublyIndexedBlock< double > TU;
    SLHA::DenseDoublyIndexedBlock< double > IMTU;
    SLHA::DenseDoublyIndexedBlock< double > TDIN;
    SLHA::DenseDoublyIndexedBlock< double > IMTDIN;
    SLHA::DenseDoublyIndexedBlock< double > TD;
    SLHA::DenseDoublyIndexedBlock< double > IMTD;
    SLHA::DenseDoublyIndexedBlock< double > TEIN;
    SLHA::DenseDoublyIndexedBlock< double > IMTEIN;
    SLHA::DenseDoublyIndexedBlock< double > TE;
    SLHA::DenseDoublyIndexedBlock< double > IMTE;
    SLHA::DenseDoublyIndexedBlock< double > TNIN;
    SLHA::DenseDoublyIndexedBlock< double > IMTNIN;
    SLHA::DenseDoublyIndexedBlock< double > TN;
    SLHA::DenseDoublyIndexedBlock< double > IMTN;
    SLHA::DenseDoublyIndexedBlock< double > YN;
    SLHA::DenseDoublyIndexedBlock< double > USQMIX;
    SLHA::DenseDoublyIndexedBlock< double > IMUSQMIX;
    SLHA::DenseDoublyIndexedBlock< double > DSQMIX;
    SLHA::DenseDoublyIndexedBlock< double > IMDSQMIX;
    SLHA::DenseDoublyIndexedBlock< double > SELMIX;
    SLHA::DenseDoublyIndexedBlock< double > IMSELMIX;
    SLHA::DenseDoublyIndexedBlock< double > SNUMIX;
    SLHA::DenseDoublyIndexedBlock< double > IMSNUMIX;
    SLHA::DenseDoublyIndexedBlock< double > SNSMIX;
    SLHA::DenseDoublyIndexedBlock< double > IMSNSMIX;
    SLHA::DenseDoublyIndexedBlock< double > SNAMIX;
    SLHA::DenseDoublyIndexedBlock< double > IMSNAMIX;
    SLHA::DenseTriplyIndexedBlock< double > RVLAMLLEIN;
    SLHA::DenseTriplyIndexedBlock< double > IMRVLAMLLEIN;
    SLHA::DenseTriplyIndexedBlock< double > RVLAMLLE;
    SLHA::DenseTriplyIndexedBlock< double > IMRVLAMLLE;
    SLHA::DenseTriplyIndexedBlock< double > RVLAMLQDIN;
    SLHA::DenseTriplyIndexedBlock< double > IMRVLAMLQDIN;
    SLHA::DenseTriplyIndexedBlock< double > RVLAMLQD;
    SLHA::DenseTriplyIndexedBlock< double > IMRVLAMLQD;
    SLHA::DenseTriplyIndexedBlock< double > RVLAMUDDIN;
    SLHA::DenseTriplyIndexedBlock< double > IMRVLAMUDDIN;
    SLHA::DenseTriplyIndexedBlock< double > RVLAMUDD;
    SLHA::DenseTriplyIndexedBlock< double > IMRVLAMUDD;
    SLHA::DenseTriplyIndexedBlock< double > RVTLLEIN;
    SLHA::DenseTriplyIndexedBlock< double > IMRVTLLEIN;
    SLHA::DenseTriplyIndexedBlock< double > RVTLLE;
    SLHA::DenseTriplyIndexedBlock< double > IMRVTLLE;
    SLHA::DenseTriplyIndexedBlock< double > RVTLQDIN;
    SLHA::DenseTriplyIndexedBlock< double > IMRVTLQDIN;
    SLHA::DenseTriplyIndexedBlock< double > RVTLQD;
    SLHA::DenseTriplyIndexedBlock< double > IMRVTLQD;
    SLHA::DenseTriplyIndexedBlock< double > RVTUDDIN;
    SLHA::DenseTriplyIndexedBlock< double > IMRVTUDDIN;
    SLHA::DenseTriplyIndexedBlock< double > RVTUDD;
    SLHA::DenseTriplyIndexedBlock< double > IMRVTUDD;
    SLHA::DenseSinglyIndexedBlock< double > RVKAPPAIN;
    SLHA::DenseSinglyIndexedBlock< double > IMRVKAPPAIN;
    SLHA::DenseSinglyIndexedBlock< double > RVKAPPA;
    SLHA::DenseSinglyIndexedBlock< double > IMRVKAPPA;
    SLHA::DenseSinglyIndexedBlock< double > RVDIN;
    SLHA::DenseSinglyIndexedBlock< double > IMRVDIN;
    SLHA::DenseSinglyIndexedBlock< double > RVD;
    SLHA::DenseSinglyIndexedBlock< double > IMRVD;
    SLHA::DenseSinglyIndexedBlock< double > RVM2LH1;
    SLHA::DenseSinglyIndexedBlock< double > IMRVM2LH1;
    SLHA::DenseSinglyIndexedBlock< double > RVSNVEVIN;
    SLHA::DenseSinglyIndexedBlock< double > IMRVSNVEVIN;
    SLHA::DenseSinglyIndexedBlock< double > RVSNVEV;
    SLHA::DenseSinglyIndexedBlock< double > IMRVSNVEV;
    SLHA::DenseDoublyIndexedBlock< double > RVNMIX;
    SLHA::DenseDoublyIndexedBlock< double > IMRVNMIX;
    SLHA::DenseDoublyIndexedBlock< double > RVUMIX;
    SLHA::DenseDoublyIndexedBlock< double > IMRVUMIX;
    SLHA::DenseDoublyIndexedBlock< double > RVVMIX;
    SLHA::DenseDoublyIndexedBlock< double > IMRVVMIX;
    SLHA::DenseDoublyIndexedBlock< double > RVHMIX;
    SLHA::DenseDoublyIndexedBlock< double > IMRVHMIX;
    SLHA::DenseDoublyIndexedBlock< double > RVAMIX;
    SLHA::DenseDoublyIndexedBlock< double > IMRVAMIX;
    SLHA::DenseDoublyIndexedBlock< double > RVLMIX;
    SLHA::DenseDoublyIndexedBlock< double > IMRVLMIX;
    SLHA::SparseSinglyIndexedBlock< double > NMSSMRUN;
    SLHA::DenseDoublyIndexedBlock< double > NMHMIX;
    SLHA::DenseDoublyIndexedBlock< double > IMNMHMIX;
    SLHA::DenseDoublyIndexedBlock< double > NMAMIX;
    SLHA::DenseDoublyIndexedBlock< double > IMNMAMIX;
    SLHA::DenseDoublyIndexedBlock< double > NMNMIX;
    SLHA::DenseDoublyIndexedBlock< double > IMNMNMIX;
  };

}

#endif /* SLHATWO_HPP_ */
