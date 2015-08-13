/*
 * StandardModel.hpp
 *
 *  Created on: Jan 8, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef STANDARDMODEL_HPP_
#define STANDARDMODEL_HPP_

#include "CodesAndDataForMassEigenstates.hpp"

namespace LHPC
{
  namespace MassSpectrumClass
  {
    // this is the spectrum of the Standard Model.
    class StandardModel : public virtual MassSpectrum
    {
    public:
      StandardModel( bool const isVerbose = false,
                     bool const neutrinosAreMajorana = false,
                     std::vector< bool >* const defaultFlags = NULL );
      virtual
      ~StandardModel();

      StandardModel&
      setNeutrinosToDirac();
      StandardModel&
      setNeutrinosToMajorana();
      MassEigenstate&
      getNeutralColorlessScalarOne();
      MassEigenstate const&
      getNeutralColorlessScalarOne() const;
      MassEigenstate&
      getHiggs(){ return getNeutralColorlessScalarOne(); }
      MassEigenstate const&
      getHiggs() const{ return getNeutralColorlessScalarOne(); }
      MassEigenstate&
      getHiggsScalarOne(){ return getNeutralColorlessScalarOne(); }
      MassEigenstate const&
      getHiggsScalarOne() const{ return getNeutralColorlessScalarOne(); }
      MassEigenstate&
      getLightHiggs(){ return getNeutralColorlessScalarOne(); }
      MassEigenstate const&
      getLightHiggs() const{ return getNeutralColorlessScalarOne(); }
      MassEigenstate&
      getPositronOne();
      MassEigenstate const&
      getPositronOne() const;
      MassEigenstate&
      getAntielectronOne(){ return getPositronOne(); }
      MassEigenstate const&
      getAntielectronOne() const{ return getPositronOne(); }
      MassEigenstate&
      getPositron(){ return getPositronOne(); }
      MassEigenstate const&
      getPositron() const{ return getPositronOne(); }
      MassEigenstate&
      getAntielectron(){ return getPositronOne(); }
      MassEigenstate const&
      getAntielectron() const{ return getPositronOne(); }
      MassEigenstate&
      getAntipositronOne();
      MassEigenstate const&
      getAntipositronOne() const;
      MassEigenstate&
      getElectronOne(){ return getAntipositronOne(); }
      MassEigenstate const&
      getElectronOne() const{ return getAntipositronOne(); }
      MassEigenstate&
      getElectron(){ return getAntipositronOne(); }
      MassEigenstate const&
      getElectron() const{ return getAntipositronOne(); }
      MassEigenstate&
      getPositronTwo();
      MassEigenstate const&
      getPositronTwo() const;
      MassEigenstate&
      getAntielectronTwo(){ return getPositronTwo(); }
      MassEigenstate const&
      getAntielectronTwo() const{ return getPositronTwo(); }
      MassEigenstate&
      getAntimuon(){ return getPositronTwo(); }
      MassEigenstate const&
      getAntimuon() const{ return getPositronTwo(); }
      MassEigenstate&
      getAntipositronTwo();
      MassEigenstate const&
      getAntipositronTwo() const;
      MassEigenstate&
      getElectronTwo(){ return getAntipositronTwo(); }
      MassEigenstate const&
      getElectronTwo() const{ return getAntipositronTwo(); }
      MassEigenstate&
      getMuon(){ return getAntipositronTwo(); }
      MassEigenstate const&
      getMuon() const{ return getAntipositronTwo(); }
      MassEigenstate&
      getPositronThree();
      MassEigenstate const&
      getPositronThree() const;
      MassEigenstate&
      getAntielectronThree(){ return getPositronThree(); }
      MassEigenstate const&
      getAntielectronThree() const{ return getPositronThree(); }
      MassEigenstate&
      getAntitau(){ return getPositronThree(); }
      MassEigenstate const&
      getAntitau() const{ return getPositronThree(); }
      MassEigenstate&
      getAntipositronThree();
      MassEigenstate const&
      getAntipositronThree() const;
      MassEigenstate&
      getElectronThree(){ return getAntipositronThree(); }
      MassEigenstate const&
      getElectronThree() const{ return getAntipositronThree(); }
      MassEigenstate&
      getTau(){ return getAntipositronThree(); }
      MassEigenstate const&
      getTau() const{ return getAntipositronThree(); }
      MassEigenstate&
      getAntineutrinoOne();
      MassEigenstate const&
      getAntineutrinoOne() const;
      MassEigenstate&
      getElectronAntineutrino(){ return getAntineutrinoOne(); }
      MassEigenstate const&
      getElectronAntineutrino() const{ return getAntineutrinoOne(); }
      MassEigenstate&
      getNeutrinoOne();
      MassEigenstate const&
      getNeutrinoOne() const;
      MassEigenstate&
      getElectronNeutrino(){ return getNeutrinoOne(); }
      MassEigenstate const&
      getElectronNeutrino() const{ return getNeutrinoOne(); }
      MassEigenstate&
      getAntineutrinoTwo();
      MassEigenstate const&
      getAntineutrinoTwo() const;
      MassEigenstate&
      getMuonAntineutrino(){ return getAntineutrinoTwo(); }
      MassEigenstate const&
      getMuonAntineutrino() const{ return getAntineutrinoTwo(); }
      MassEigenstate&
      getNeutrinoTwo();
      MassEigenstate const&
      getNeutrinoTwo() const;
      MassEigenstate&
      getMuonNeutrino(){ return getNeutrinoTwo(); }
      MassEigenstate const&
      getMuonNeutrino() const{ return getNeutrinoTwo(); }
      MassEigenstate&
      getAntineutrinoThree();
      MassEigenstate const&
      getAntineutrinoThree() const;
      MassEigenstate&
      getTauAntineutrino(){ return getAntineutrinoThree(); }
      MassEigenstate const&
      getTauAntineutrino() const{ return getAntineutrinoThree(); }
      MassEigenstate&
      getNeutrinoThree();
      MassEigenstate const&
      getNeutrinoThree() const;
      MassEigenstate&
      getTauNeutrino(){ return getTauNeutrino(); }
      MassEigenstate const&
      getTauNeutrino() const{ return getTauNeutrino(); }
      MassEigenstate&
      getAntidownOne();
      MassEigenstate const&
      getAntidownOne() const;
      MassEigenstate&
      getAntidown(){ return getAntidownOne(); }
      MassEigenstate const&
      getAntidown() const{ return getAntidownOne(); }
      MassEigenstate&
      getDownOne();
      MassEigenstate const&
      getDownOne() const;
      MassEigenstate&
      getDown(){ return getDownOne(); }
      MassEigenstate const&
      getDown() const{ return getDownOne(); }
      MassEigenstate&
      getAntidownTwo();
      MassEigenstate const&
      getAntidownTwo() const;
      MassEigenstate&
      getAntistrange(){ return getAntidownTwo(); }
      MassEigenstate const&
      getAntistrange() const{ return getAntidownTwo(); }
      MassEigenstate&
      getDownTwo();
      MassEigenstate const&
      getDownTwo() const;
      MassEigenstate&
      getStrange(){ return getDownTwo(); }
      MassEigenstate const&
      getStrange() const{ return getDownTwo(); }
      MassEigenstate&
      getAntidownThree();
      MassEigenstate const&
      getAntidownThree() const;
      MassEigenstate&
      getAntibottom(){ return getAntidownThree(); }
      MassEigenstate const&
      getAntibottom() const{ return getAntidownThree(); }
      MassEigenstate&
      getDownThree();
      MassEigenstate const&
      getDownThree() const;
      MassEigenstate&
      getBottom(){ return getDownThree(); }
      MassEigenstate const&
      getBottom() const{ return getDownThree(); }
      MassEigenstate&
      getUpOne();
      MassEigenstate const&
      getUpOne() const;
      MassEigenstate&
      getUp(){ return getUpOne(); }
      MassEigenstate const&
      getUp() const{ return getUpOne(); }
      MassEigenstate&
      getAntiupOne();
      MassEigenstate const&
      getAntiupOne() const;
      MassEigenstate&
      getAntiup(){ return getAntiupOne(); }
      MassEigenstate const&
      getAntiup() const{ return getAntiupOne(); }
      MassEigenstate&
      getUpTwo();
      MassEigenstate const&
      getUpTwo() const;
      MassEigenstate&
      getCharm(){ return getUpTwo(); }
      MassEigenstate const&
      getCharm() const{ return getUpTwo(); }
      MassEigenstate&
      getAntiupTwo();
      MassEigenstate const&
      getAntiupTwo() const;
      MassEigenstate&
      getAnticharm(){ return getAntiupTwo(); }
      MassEigenstate const&
      getAnticharm() const{ return getAntiupTwo(); }
      MassEigenstate&
      getUpThree();
      MassEigenstate const&
      getUpThree() const;
      MassEigenstate&
      getTop(){ return getUpThree(); }
      MassEigenstate const&
      getTop() const{ return getUpThree(); }
      MassEigenstate&
      getAntiupThree();
      MassEigenstate const&
      getAntiupThree() const;
      MassEigenstate&
      getAntitop(){ return getAntiupThree(); }
      MassEigenstate const&
      getAntitop() const{ return getAntiupThree(); }
      MassEigenstate&
      getPhoton();
      MassEigenstate const&
      getPhoton() const;
      MassEigenstate&
      getWPlusBosonOne();
      MassEigenstate const&
      getWPlusBosonOne() const;
      MassEigenstate&
      getWPlus(){ return getWPlusBosonOne(); }
      MassEigenstate const&
      getWPlus() const{ return getWPlusBosonOne(); }
      MassEigenstate&
      getWMinusBosonOne();
      MassEigenstate const&
      getWMinusBosonOne() const;
      MassEigenstate&
      getWMinus(){ return getWMinusBosonOne(); }
      MassEigenstate const&
      getWMinus() const{ return getWMinusBosonOne(); }
      MassEigenstate&
      getZBosonOne();
      MassEigenstate const&
      getZBosonOne() const;
      MassEigenstate&
      getZBoson(){ return getZBosonOne(); }
      MassEigenstate const&
      getZBoson() const{ return getZBosonOne(); }
      MassEigenstate&
      getZ(){ return getZBosonOne(); }
      MassEigenstate const&
      getZ() const{ return getZBosonOne(); }
      MassEigenstate&
      getGluon();
      MassEigenstate const&
      getGluon() const;
      std::vector< MassEigenstate* >&
      getPositiveLeptons();
      std::vector< MassEigenstate* > const&
      getPositiveLeptons() const;
      std::vector< MassEigenstate* >&
      getNegativeLeptons();
      std::vector< MassEigenstate* > const&
      getNegativeLeptons() const;
      std::vector< MassEigenstate* >&
      getAntineutrinos();
      std::vector< MassEigenstate* > const&
      getAntineutrinos() const;
      std::vector< MassEigenstate* >&
      getNeutrinos();
      std::vector< MassEigenstate* > const&
      getNeutrinos() const;
      std::vector< MassEigenstate* >&
      getDownAntiquarks();
      std::vector< MassEigenstate* > const&
      getDownAntiquarks() const;
      std::vector< MassEigenstate* >&
      getDownQuarks();
      std::vector< MassEigenstate* > const&
      getDownQuarks() const;
      std::vector< MassEigenstate* >&
      getUpQuarks();
      std::vector< MassEigenstate* > const&
      getUpQuarks() const;
      std::vector< MassEigenstate* >&
      getUpAntiquarks();
      std::vector< MassEigenstate* > const&
      getUpAntiquarks() const;


    protected:
      MassEigenstate neutralColorlessScalarOne;
      MassEigenstate positronOne;
      MassEigenstate antipositronOne;
      MassEigenstate positronTwo;
      MassEigenstate antipositronTwo;
      MassEigenstate positronThree;
      MassEigenstate antipositronThree;
      MassEigenstate antineutrinoOne;
      MassEigenstate neutrinoOne;
      MassEigenstate antineutrinoTwo;
      MassEigenstate neutrinoTwo;
      MassEigenstate antineutrinoThree;
      MassEigenstate neutrinoThree;
      MassEigenstate antidownOne;
      MassEigenstate downOne;
      MassEigenstate antidownTwo;
      MassEigenstate downTwo;
      MassEigenstate antidownThree;
      MassEigenstate downThree;
      MassEigenstate upOne;
      MassEigenstate antiupOne;
      MassEigenstate upTwo;
      MassEigenstate antiupTwo;
      MassEigenstate upThree;
      MassEigenstate antiupThree;
      MassEigenstate photonBoson;
      MassEigenstate zBosonOne;
      MassEigenstate wPlusBosonOne;
      MassEigenstate wMinusBosonOne;
      MassEigenstate gluonBoson;
      std::vector< MassEigenstate* > positiveLeptonPointers;
      std::vector< MassEigenstate* > negativeLeptonPointers;
      std::vector< MassEigenstate* > antineutrinoPointers;
      std::vector< MassEigenstate* > neutrinoPointers;
      std::vector< MassEigenstate* > downAntiquarkPointers;
      std::vector< MassEigenstate* > downQuarkPointers;
      std::vector< MassEigenstate* > upQuarkPointers;
      std::vector< MassEigenstate* > upAntiquarkPointers;
      MassEigenstate::MassEigenstatesPairedWithBr defaultDecayFiller;
    };
    typedef StandardModel SM;



    inline StandardModel&
    StandardModel::setNeutrinosToDirac()
    {
      // only the neutrinos which are charge-conjugates of the antineutrinos
      // are set back to having charge conjugates.
      MassEigenstate* neutrinoFinder;
      for( int neutrinoIndex( antineutrinoPointers.size() - 1 );
           0 <= neutrinoIndex;
           --neutrinoIndex )
      {
        neutrinoFinder = getMassEigenstate(
                            antineutrinoPointers[ neutrinoIndex ]->getCode() );
        if( NULL != neutrinoFinder )
        {
          neutrinoFinder->setToBeChargeConjugate(
                                       antineutrinoPointers[ neutrinoIndex ] );
        }
      }
      return *this;
    }

    inline StandardModel&
    StandardModel::setNeutrinosToMajorana()
    {
      for( int neutrinoIndex( neutrinoPointers.size() - 1 );
           0 <= neutrinoIndex;
           --neutrinoIndex )
      {
        neutrinoPointers[ neutrinoIndex ]->setToBeSelfConjugate();
      }
      return *this;
    }

    inline MassEigenstate&
    StandardModel::getNeutralColorlessScalarOne()
    {
      return neutralColorlessScalarOne;
    }

    inline MassEigenstate const&
    StandardModel::getNeutralColorlessScalarOne() const
    {
      return neutralColorlessScalarOne;
    }

    inline MassEigenstate&
    StandardModel::getPositronOne()
    {
      return positronOne;
    }

    inline MassEigenstate const&
    StandardModel::getPositronOne() const
    {
      return positronOne;
    }

    inline MassEigenstate&
    StandardModel::getAntipositronOne()
    {
      return antipositronOne;
    }

    inline MassEigenstate const&
    StandardModel::getAntipositronOne() const
    {
      return antipositronOne;
    }

    inline MassEigenstate&
    StandardModel::getPositronTwo()
    {
      return positronTwo;
    }

    inline MassEigenstate const&
    StandardModel::getPositronTwo() const
    {
      return positronTwo;
    }

    inline MassEigenstate&
    StandardModel::getAntipositronTwo()
    {
      return antipositronTwo;
    }

    inline MassEigenstate const&
    StandardModel::getAntipositronTwo() const
    {
      return antipositronTwo;
    }

    inline MassEigenstate&
    StandardModel::getPositronThree()
    {
      return positronThree;
    }

    inline MassEigenstate const&
    StandardModel::getPositronThree() const
    {
      return positronThree;
    }

    inline MassEigenstate&
    StandardModel::getAntipositronThree()
    {
      return antipositronThree;
    }

    inline MassEigenstate const&
    StandardModel::getAntipositronThree() const
    {
      return antipositronThree;
    }

    inline MassEigenstate&
    StandardModel::getAntineutrinoOne()
    {
      return antineutrinoOne;
    }

    inline MassEigenstate const&
    StandardModel::getAntineutrinoOne() const
    {
      return antineutrinoOne;
    }

    inline MassEigenstate&
    StandardModel::getNeutrinoOne()
    {
      return neutrinoOne;
    }

    inline MassEigenstate const&
    StandardModel::getNeutrinoOne() const
    {
      return neutrinoOne;
    }

    inline MassEigenstate&
    StandardModel::getAntineutrinoTwo()
    {
      return antineutrinoTwo;
    }

    inline MassEigenstate const&
    StandardModel::getAntineutrinoTwo() const
    {
      return antineutrinoTwo;
    }

    inline MassEigenstate&
    StandardModel::getNeutrinoTwo()
    {
      return neutrinoTwo;
    }

    inline MassEigenstate const&
    StandardModel::getNeutrinoTwo() const
    {
      return neutrinoTwo;
    }

    inline MassEigenstate&
    StandardModel::getAntineutrinoThree()
    {
      return antineutrinoThree;
    }

    inline MassEigenstate const&
    StandardModel::getAntineutrinoThree() const
    {
      return antineutrinoThree;
    }

    inline MassEigenstate&
    StandardModel::getNeutrinoThree()
    {
      return neutrinoThree;
    }

    inline MassEigenstate const&
    StandardModel::getNeutrinoThree() const
    {
      return neutrinoThree;
    }

    inline MassEigenstate&
    StandardModel::getAntidownOne()
    {
      return antidownOne;
    }

    inline MassEigenstate const&
    StandardModel::getAntidownOne() const
    {
      return antidownOne;
    }

    inline MassEigenstate&
    StandardModel::getDownOne()
    {
      return downOne;
    }

    inline MassEigenstate const&
    StandardModel::getDownOne() const
    {
      return downOne;
    }

    inline MassEigenstate&
    StandardModel::getAntidownTwo()
    {
      return antidownTwo;
    }

    inline MassEigenstate const&
    StandardModel::getAntidownTwo() const
    {
      return antidownTwo;
    }

    inline MassEigenstate&
    StandardModel::getDownTwo()
    {
      return downTwo;
    }

    inline MassEigenstate const&
    StandardModel::getDownTwo() const
    {
      return downTwo;
    }

    inline MassEigenstate&
    StandardModel::getAntidownThree()
    {
      return antidownThree;
    }

    inline MassEigenstate const&
    StandardModel::getAntidownThree() const
    {
      return antidownThree;
    }

    inline MassEigenstate&
    StandardModel::getDownThree()
    {
      return downThree;
    }

    inline MassEigenstate const&
    StandardModel::getDownThree() const
    {
      return downThree;
    }

    inline MassEigenstate&
    StandardModel::getUpOne()
    {
      return upOne;
    }

    inline MassEigenstate const&
    StandardModel::getUpOne() const
    {
      return upOne;
    }

    inline MassEigenstate&
    StandardModel::getAntiupOne()
    {
      return antiupOne;
    }

    inline MassEigenstate const&
    StandardModel::getAntiupOne() const
    {
      return antiupOne;
    }

    inline MassEigenstate&
    StandardModel::getUpTwo()
    {
      return upTwo;
    }

    inline MassEigenstate const&
    StandardModel::getUpTwo() const
    {
      return upTwo;
    }

    inline MassEigenstate&
    StandardModel::getAntiupTwo()
    {
      return antiupTwo;
    }

    inline MassEigenstate const&
    StandardModel::getAntiupTwo() const
    {
      return antiupTwo;
    }

    inline MassEigenstate&
    StandardModel::getUpThree()
    {
      return upThree;
    }

    inline MassEigenstate const&
    StandardModel::getUpThree() const
    {
      return upThree;
    }

    inline MassEigenstate&
    StandardModel::getAntiupThree()
    {
      return antiupThree;
    }

    inline MassEigenstate const&
    StandardModel::getAntiupThree() const
    {
      return antiupThree;
    }

    inline MassEigenstate&
    StandardModel::getPhoton()
    {
      return photonBoson;
    }

    inline MassEigenstate const&
    StandardModel::getPhoton() const
    {
      return photonBoson;
    }

    inline MassEigenstate&
    StandardModel::getWPlusBosonOne()
    {
      return wPlusBosonOne;
    }

    inline MassEigenstate const&
    StandardModel::getWPlusBosonOne() const
    {
      return wPlusBosonOne;
    }

    inline MassEigenstate&
    StandardModel::getWMinusBosonOne()
    {
      return wMinusBosonOne;
    }

    inline MassEigenstate const&
    StandardModel::getWMinusBosonOne() const
    {
      return wMinusBosonOne;
    }

    inline MassEigenstate&
    StandardModel::getZBosonOne()
    {
      return zBosonOne;
    }

    inline MassEigenstate const&
    StandardModel::getZBosonOne() const
    {
      return zBosonOne;
    }

    inline MassEigenstate&
    StandardModel::getGluon()
    {
      return gluonBoson;
    }

    inline MassEigenstate const&
    StandardModel::getGluon() const
    {
      return gluonBoson;
    }

    inline std::vector< MassEigenstate* >&
    StandardModel::getPositiveLeptons()
    {
      return positiveLeptonPointers;
    }

    inline std::vector< MassEigenstate* > const&
    StandardModel::getPositiveLeptons() const
    {
      return positiveLeptonPointers;
    }

    inline std::vector< MassEigenstate* >&
    StandardModel::getNegativeLeptons()
    {
      return negativeLeptonPointers;
    }

    inline std::vector< MassEigenstate* > const&
    StandardModel::getNegativeLeptons() const
    {
      return negativeLeptonPointers;
    }

    inline std::vector< MassEigenstate* >&
    StandardModel::getAntineutrinos()
    {
      return antineutrinoPointers;
    }

    inline std::vector< MassEigenstate* > const&
    StandardModel::getAntineutrinos() const
    {
      return antineutrinoPointers;
    }

    inline std::vector< MassEigenstate* >&
    StandardModel::getNeutrinos()
    {
      return neutrinoPointers;
    }

    inline std::vector< MassEigenstate* > const&
    StandardModel::getNeutrinos() const
    {
      return neutrinoPointers;
    }

    inline std::vector< MassEigenstate* >&
    StandardModel::getDownAntiquarks()
    {
      return downAntiquarkPointers;
    }

    inline std::vector< MassEigenstate* > const&
    StandardModel::getDownAntiquarks() const
    {
      return downAntiquarkPointers;
    }

    inline std::vector< MassEigenstate* >&
    StandardModel::getDownQuarks()
    {
      return downQuarkPointers;
    }

    inline std::vector< MassEigenstate* > const&
    StandardModel::getDownQuarks() const
    {
      return downQuarkPointers;
    }

    inline std::vector< MassEigenstate* >&
    StandardModel::getUpQuarks()
    {
      return upQuarkPointers;
    }

    inline std::vector< MassEigenstate* > const&
    StandardModel::getUpQuarks() const
    {
      return upQuarkPointers;
    }

    inline std::vector< MassEigenstate* >&
    StandardModel::getUpAntiquarks()
    {
      return upAntiquarkPointers;
    }

    inline std::vector< MassEigenstate* > const&
    StandardModel::getUpAntiquarks() const
    {
      return upAntiquarkPointers;
    }

  }
  typedef MassSpectrumClass::StandardModel StandardModelSpectrum;
  typedef MassSpectrumClass::StandardModel SmSpectrum;

}

#endif /* STANDARDMODEL_HPP_ */
