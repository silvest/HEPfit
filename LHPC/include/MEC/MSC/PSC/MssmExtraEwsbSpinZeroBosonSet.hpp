/*
 * MssmExtraEwsbSpinZeroBosonSet.hpp
 *
 *  Created on: Jan 18, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef MSSMEXTRAEWSBSPINZEROBOSONSET_HPP_
#define MSSMEXTRAEWSBSPINZEROBOSONSET_HPP_

#include "../CodesAndDataForMassEigenstates.hpp"

namespace LHPC
{
  namespace MassSpectrumClass
  {
    class MssmExtraEwsbSpinZeroBosonSet : public virtual MassSpectrum
    {
    public:
      MssmExtraEwsbSpinZeroBosonSet( bool const isVerbose = false,
                              std::vector< bool >* const defaultFlags = NULL );
      virtual
      ~MssmExtraEwsbSpinZeroBosonSet();

      MassEigenstate&
      getNeutralColorlessScalarTwo();
      MassEigenstate const&
      getNeutralColorlessScalarTwo() const;
      MassEigenstate&
      getHiggsScalarTwo(){ return getNeutralColorlessScalarTwo(); }
      MassEigenstate const&
      getHiggsScalarTwo() const{ return getNeutralColorlessScalarTwo(); }
      virtual MassEigenstate&
      getHeavyHiggs(){ return getNeutralColorlessScalarTwo(); }
      virtual MassEigenstate const&
      getHeavyHiggs() const{ return getNeutralColorlessScalarTwo(); }
      MassEigenstate&
      getNeutralColorlessPseudoscalarOne();
      MassEigenstate const&
      getNeutralColorlessPseudoscalarOne() const;
      MassEigenstate&
      getHiggsPseudoscalarOne(){ return getNeutralColorlessPseudoscalarOne(); }
      MassEigenstate const&
      getHiggsPseudoscalarOne() const{ return
                                       getNeutralColorlessPseudoscalarOne(); }
      MassEigenstate&
      getHiggsPseudoscalar(){ return getNeutralColorlessPseudoscalarOne(); }
      MassEigenstate const&
      getHiggsPseudoscalar() const{ return
                                    getNeutralColorlessPseudoscalarOne(); }
      MassEigenstate&
      getPositiveColorlessSpinZeroBosonOne();
      MassEigenstate const&
      getPositiveColorlessSpinZeroBosonOne() const;
      MassEigenstate&
      getPositiveHiggsOne(){ return getPositiveColorlessSpinZeroBosonOne(); }
      MassEigenstate const&
      getPositiveHiggsOne() const{ return
                                   getPositiveColorlessSpinZeroBosonOne(); }
      MassEigenstate&
      getPositiveHiggs(){ return getPositiveColorlessSpinZeroBosonOne(); }
      MassEigenstate const&
      getPositiveHiggs() const{ return
                                getPositiveColorlessSpinZeroBosonOne(); }
      MassEigenstate&
      getNegativeColorlessSpinZeroBosonOne();
      MassEigenstate const&
      getNegativeColorlessSpinZeroBosonOne() const;
      MassEigenstate&
      getNegativeHiggsOne(){ return getNegativeColorlessSpinZeroBosonOne(); }
      MassEigenstate const&
      getNegativeHiggsOne() const{ return
                                   getNegativeColorlessSpinZeroBosonOne(); }
      MassEigenstate&
      getNegativeHiggs(){ return getNegativeColorlessSpinZeroBosonOne(); }
      MassEigenstate const&
      getNegativeHiggs() const{ return
                                getNegativeColorlessSpinZeroBosonOne(); }
      std::vector< MassEigenstate* >&
      getNeutralScalarsAndPseudoscalars();
      std::vector< MassEigenstate* > const&
      getNeutralScalarsAndPseudoscalars() const;
      std::vector< MassEigenstate* >&
      getChargedColorlessSpinZeroBosons();
      std::vector< MassEigenstate* > const&
      getChargedColorlessSpinZeroBosons() const;
      std::vector< MassEigenstate* >&
      getEwsbSpinZeroAndOneBosons();
      std::vector< MassEigenstate* > const&
      getEwsbSpinZeroAndOneBosons() const;
      std::vector< MassEigenstate* >&
      getNeutralEwsbSpinZeroAndOneBosons();
      std::vector< MassEigenstate* > const&
      getNeutralEwsbSpinZeroAndOneBosons() const;
      std::vector< MassEigenstate* >&
      getChargedEwsbSpinZeroAndOneBosons();
      std::vector< MassEigenstate* > const&
      getChargedEwsbSpinZeroAndOneBosons() const;


    protected:
      MassEigenstate neutralColorlessScalarTwo;
      MassEigenstate neutralColorlessPseudoscalarOne;
      MassEigenstate positiveColorlessSpinZeroBosonOne;
      MassEigenstate negativeColorlessSpinZeroBosonOne;
      std::vector< MassEigenstate* > neutralScalarsAndPseudoscalarPointers;
      std::vector< MassEigenstate* > chargedColorlessSpinZeroBosonPointers;
      std::vector< MassEigenstate* > ewsbSpinZeroAndOneBosonPointers;
      std::vector< MassEigenstate* > neutralEwsbSpinZeroAndOneBosonPointers;
      std::vector< MassEigenstate* > chargedEwsbSpinZeroAndOneBosonPointers;
    };



    inline MassEigenstate&
    MssmExtraEwsbSpinZeroBosonSet::getNeutralColorlessScalarTwo()
    {
      return neutralColorlessScalarTwo;
    }

    inline MassEigenstate const&
    MssmExtraEwsbSpinZeroBosonSet::getNeutralColorlessScalarTwo() const
    {
      return neutralColorlessScalarTwo;
    }

    inline MassEigenstate&
    MssmExtraEwsbSpinZeroBosonSet::getNeutralColorlessPseudoscalarOne()
    {
      return neutralColorlessPseudoscalarOne;
    }

    inline MassEigenstate const&
    MssmExtraEwsbSpinZeroBosonSet::getNeutralColorlessPseudoscalarOne() const
    {
      return neutralColorlessPseudoscalarOne;
    }

    inline MassEigenstate&
    MssmExtraEwsbSpinZeroBosonSet::getPositiveColorlessSpinZeroBosonOne()
    {
      return positiveColorlessSpinZeroBosonOne;
    }

    inline MassEigenstate const&
    MssmExtraEwsbSpinZeroBosonSet::getPositiveColorlessSpinZeroBosonOne() const
    {
      return positiveColorlessSpinZeroBosonOne;
    }

    inline MassEigenstate&
    MssmExtraEwsbSpinZeroBosonSet::getNegativeColorlessSpinZeroBosonOne()
    {
      return negativeColorlessSpinZeroBosonOne;
    }

    inline MassEigenstate const&
    MssmExtraEwsbSpinZeroBosonSet::getNegativeColorlessSpinZeroBosonOne() const
    {
      return negativeColorlessSpinZeroBosonOne;
    }

    inline std::vector< MassEigenstate* >&
    MssmExtraEwsbSpinZeroBosonSet::getNeutralScalarsAndPseudoscalars()
    {
      return neutralScalarsAndPseudoscalarPointers;
    }

    inline std::vector< MassEigenstate* > const&
    MssmExtraEwsbSpinZeroBosonSet::getNeutralScalarsAndPseudoscalars() const
    {
      return neutralScalarsAndPseudoscalarPointers;
    }

    inline std::vector< MassEigenstate* >&
    MssmExtraEwsbSpinZeroBosonSet::getChargedColorlessSpinZeroBosons()
    {
      return chargedColorlessSpinZeroBosonPointers;
    }

    inline std::vector< MassEigenstate* > const&
    MssmExtraEwsbSpinZeroBosonSet::getChargedColorlessSpinZeroBosons() const
    {
      return chargedColorlessSpinZeroBosonPointers;
    }

    inline std::vector< MassEigenstate* >&
    MssmExtraEwsbSpinZeroBosonSet::getEwsbSpinZeroAndOneBosons()
    {
      return ewsbSpinZeroAndOneBosonPointers;
    }

    inline std::vector< MassEigenstate* > const&
    MssmExtraEwsbSpinZeroBosonSet::getEwsbSpinZeroAndOneBosons() const
    {
      return ewsbSpinZeroAndOneBosonPointers;
    }

    inline std::vector< MassEigenstate* >&
    MssmExtraEwsbSpinZeroBosonSet::getNeutralEwsbSpinZeroAndOneBosons()
    {
      return neutralEwsbSpinZeroAndOneBosonPointers;
    }

    inline std::vector< MassEigenstate* > const&
    MssmExtraEwsbSpinZeroBosonSet::getNeutralEwsbSpinZeroAndOneBosons() const
    {
      return neutralEwsbSpinZeroAndOneBosonPointers;
    }

    inline std::vector< MassEigenstate* >&
    MssmExtraEwsbSpinZeroBosonSet::getChargedEwsbSpinZeroAndOneBosons()
    {
      return chargedEwsbSpinZeroAndOneBosonPointers;
    }

    inline std::vector< MassEigenstate* > const&
    MssmExtraEwsbSpinZeroBosonSet::getChargedEwsbSpinZeroAndOneBosons() const
    {
      return chargedEwsbSpinZeroAndOneBosonPointers;
    }

  }

}

#endif /* MSSMEXTRAEWSBSPINZEROBOSONSET_HPP_ */
