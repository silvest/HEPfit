/*
 * NmssmExtraEwsbSpinZeroBosonSet.hpp
 *
 *  Created on: Jan 27, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef NMSSMEXTRAEWSBSPINZEROBOSONSET_HPP_
#define NMSSMEXTRAEWSBSPINZEROBOSONSET_HPP_

#include "../CodesAndDataForMassEigenstates.hpp"
#include "MssmExtraEwsbSpinZeroBosonSet.hpp"

namespace LHPC
{
  namespace MassSpectrumClass
  {
    class NmssmExtraEwsbSpinZeroBosonSet : public virtual MassSpectrum,
                                           public MssmExtraEwsbSpinZeroBosonSet
    {
    public:
      NmssmExtraEwsbSpinZeroBosonSet( bool const isVerbose = false,
                              std::vector< bool >* const defaultFlags = NULL );
      virtual
      ~NmssmExtraEwsbSpinZeroBosonSet();

      virtual MassEigenstate&
      getMediumHiggs(){ return getNeutralColorlessScalarTwo(); }
      virtual MassEigenstate const&
      getMediumHiggs() const{ return getNeutralColorlessScalarTwo(); }
      MassEigenstate&
      getNeutralColorlessScalarThree();
      MassEigenstate const&
      getNeutralColorlessScalarThree() const;
      MassEigenstate&
      getHiggsScalarThree(){ return getNeutralColorlessScalarThree(); }
      MassEigenstate const&
      getHiggsScalarThree() const{ return getNeutralColorlessScalarThree(); }
      virtual MassEigenstate&
      getHeavyHiggs(){ return getNeutralColorlessScalarThree(); }
      virtual MassEigenstate const&
      getHeavyHiggs() const{ return getNeutralColorlessScalarThree(); }
      MassEigenstate&
      getLightHiggsPseudoscalar(){ return
                                   getNeutralColorlessPseudoscalarOne(); }
      MassEigenstate const&
      getLightHiggsPseudoscalar() const{ return
                                        getNeutralColorlessPseudoscalarOne(); }
      MassEigenstate&
      getNeutralColorlessPseudoscalarTwo();
      MassEigenstate const&
      getNeutralColorlessPseudoscalarTwo() const;
      MassEigenstate&
      getHiggsPseudoscalarTwo(){ return getNeutralColorlessPseudoscalarTwo(); }
      MassEigenstate const&
      getHiggsPseudoscalarTwo() const{ return
                                       getNeutralColorlessPseudoscalarTwo(); }
      MassEigenstate&
      getHeavyHiggsPseudoscalar(){ return
                                   getNeutralColorlessPseudoscalarTwo(); }
      MassEigenstate const&
      getHeavyHiggsPseudoscalar() const{ return
                                        getNeutralColorlessPseudoscalarTwo(); }


    protected:
      MassEigenstate neutralColorlessScalarThree;
      MassEigenstate neutralColorlessPseudoscalarTwo;
    };



    inline MassEigenstate&
    NmssmExtraEwsbSpinZeroBosonSet::getNeutralColorlessScalarThree()
    {
      return neutralColorlessScalarThree;
    }

    inline MassEigenstate const&
    NmssmExtraEwsbSpinZeroBosonSet::getNeutralColorlessScalarThree() const
    {
      return neutralColorlessScalarThree;
    }

    inline MassEigenstate&
    NmssmExtraEwsbSpinZeroBosonSet::getNeutralColorlessPseudoscalarTwo()
    {
      return neutralColorlessPseudoscalarTwo;
    }

    inline MassEigenstate const&
    NmssmExtraEwsbSpinZeroBosonSet::getNeutralColorlessPseudoscalarTwo() const
    {
      return neutralColorlessPseudoscalarTwo;
    }

  }

}

#endif /* NMSSMEXTRAEWSBSPINZEROBOSONSET_HPP_ */
