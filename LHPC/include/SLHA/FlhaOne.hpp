/*
 * FlhaOne.hpp
 *
 *  Created on: Apr 1, 2012 (really!)
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef FLHAONE_HPP_
#define FLHAONE_HPP_

#include "BlockTypes.hpp"
#include "SlhaParser.hpp"
#include "FlavorObservable.hpp"
#include "FlavorObservableError.hpp"

namespace LHPC
{
  /* this is just a collection of instances of classes derived from SlhaBlock
   * which covers all the blocks specified in FLHA1
   * (Comput.Phys.Commun. 183 (2012) 285-298 [arXiv:1008.0762]).
   */
  class FlhaOne
  {
  public:
    FlhaOne( SlhaParser& fileParser,
             bool const isVerbose = false );
    virtual
    ~FlhaOne();

    SLHA::LinesAsStringsBlock FCINFO;
    SLHA::SparseSinglyIndexedBlock< double > MODSEL;
    SLHA::SparseSinglyIndexedBlock< double > FMODSEL;
    /* the string value of entry 99 can be found by using
     * getLineWithoutComment( int const ), looking for a string that begins
     * with some whitespace then "99".
     */
    SLHA::SparseSinglyIndexedBlock< double > SMINPUTS;
    SLHA::DenseSinglyIndexedBlock< double > VCKMIN;
    SLHA::DenseSinglyIndexedBlock< double > UPMNSIN;
    SLHA::DenseSinglyIndexedBlock< double > VCKM;
    SLHA::DenseSinglyIndexedBlock< double > IMVCKM;
    SLHA::DenseSinglyIndexedBlock< double > UPMNS;
    SLHA::DenseSinglyIndexedBlock< double > IMUPMNS;
    SLHA::SinglyIndexedMultipleEntriesBlock< RunningConstant > FMASS;
    SLHA::SinglyIndexedMultipleEntriesBlock< RunningConstantError > FMASSERR;
    SLHA::SparseSinglyIndexedBlock< double > FLIFE;
    SLHA::SparseSinglyIndexedBlock< std::pair< double, double > > FLIFEERR;
    SLHA::SparseDoublyIndexedBlock< double > FCONST;
    SLHA::SparseSinglyIndexedBlock< std::pair< double, double > > FCONSTERR;
    SLHA::SparseQuadruplyIndexedBlock< double > FCONSTRATIO;
    SLHA::SparseQuadruplyIndexedBlock< std::pair< double, double > >
    FCONSTRATIOERR;
    SLHA::DoublyIndexedMultipleEntriesBlock< RunningConstant > FBAG;
    SLHA::DoublyIndexedMultipleEntriesBlock< RunningConstantError > FBAGERR;
    SLHA::SparseQuadruplyIndexedBlock< double > FWCOEF;
    SLHA::SparseQuadruplyIndexedBlock< std::pair< double, double > >
    FWCOEFERR;
    SLHA::SparseQuadruplyIndexedBlock< double > IMFWCOEF;
    SLHA::SparseQuadruplyIndexedBlock< std::pair< double, double > >
    IMFWCOEFERR;
    SLHA::DoublyIndexedMultipleEntriesBlock< FlavorObservable > FOBS;
    SLHA::DoublyIndexedMultipleEntriesBlock< FlavorObservableError > FOBSERR;
    SLHA::DoublyIndexedMultipleEntriesBlock< FlavorObservable > FOBSSM;
    SLHA::DoublyIndexedMultipleEntriesBlock< FlavorObservableError > FOBSSMERR;
    SLHA::DoublyIndexedMultipleEntriesBlock< FlavorObservable > FPARAM;
    SLHA::DoublyIndexedMultipleEntriesBlock< FlavorObservableError > FPARAMERR;


  protected:
    SlhaParser& fileParser;
    bool const isVerbose;
  };

}

#endif /* FLHAONE_HPP_ */
