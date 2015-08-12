/*
 * SquarksOneToSix.hpp
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

#ifndef SQUARKSONETOSIX_HPP_
#define SQUARKSONETOSIX_HPP_

#include "../CodesAndDataForMassEigenstates.hpp"

namespace LHPC
{
  namespace MassSpectrumClass
  {
    class SquarksOneToSix : public virtual MassSpectrum
    {
    public:
      SquarksOneToSix( bool const isVerbose = false,
                       bool const flavorConserving = false,
                       std::vector< bool >* const defaultFlags = NULL );
      virtual
      ~SquarksOneToSix();

      MassEigenstate&
      getAntisdownOne();
      MassEigenstate const&
      getAntisdownOne() const;
      MassEigenstate&
      getAntisdownL(){ return getAntisdownOne(); }
      MassEigenstate const&
      getAntisdownL() const{ return getAntisdownOne(); }
      MassEigenstate&
      getSdownOne();
      MassEigenstate const&
      getSdownOne() const;
      MassEigenstate&
      getSdownL(){ return getSdownOne(); }
      MassEigenstate const&
      getSdownL() const{ return getSdownOne(); }
      MassEigenstate&
      getAntisdownTwo();
      MassEigenstate const&
      getAntisdownTwo() const;
      MassEigenstate&
      getAntisstrangeL(){ return getAntisdownTwo(); }
      MassEigenstate const&
      getAntisstrangeL() const{ return getAntisdownTwo(); }
      MassEigenstate&
      getSdownTwo();
      MassEigenstate const&
      getSdownTwo() const;
      MassEigenstate&
      getStrangeL(){ return getSdownTwo(); }
      MassEigenstate const&
      getStrangeL() const{ return getSdownTwo(); }
      MassEigenstate&
      getAntisdownThree();
      MassEigenstate const&
      getAntisdownThree() const;
      MassEigenstate&
      getAntisbottomOne(){ return getAntisdownThree(); }
      MassEigenstate const&
      getAntisbottomOne() const{ return getAntisdownThree(); }
      MassEigenstate&
      getSdownThree();
      MassEigenstate const&
      getSdownThree() const;
      MassEigenstate&
      getSbottomOne(){ return getSdownThree(); }
      MassEigenstate const&
      getSbottomOne() const{ return getSdownThree(); }
      MassEigenstate&
      getAntisdownFour();
      MassEigenstate const&
      getAntisdownFour() const;
      MassEigenstate&
      getAntisdownR(){ return getAntisdownFour(); }
      MassEigenstate const&
      getAntisdownR() const{ return getAntisdownFour(); }
      MassEigenstate&
      getSdownFour();
      MassEigenstate const&
      getSdownFour() const;
      MassEigenstate&
      getSdownR(){ return getSdownFour(); }
      MassEigenstate const&
      getSdownR() const{ return getSdownFour(); }
      MassEigenstate&
      getAntisdownFive();
      MassEigenstate const&
      getAntisdownFive() const;
      MassEigenstate&
      getAntisstrangeR(){ return getAntisdownFive(); }
      MassEigenstate const&
      getAntisstrangeR() const{ return getAntisdownFive(); }
      MassEigenstate&
      getSdownFive();
      MassEigenstate const&
      getSdownFive() const;
      MassEigenstate&
      getStrangeR(){ return getSdownFive(); }
      MassEigenstate const&
      getStrangeR() const{ return getSdownFive(); }
      MassEigenstate&
      getAntisdownSix();
      MassEigenstate const&
      getAntisdownSix() const;
      MassEigenstate&
      getAntisbottomTwo(){ return getAntisdownSix(); }
      MassEigenstate const&
      getAntisbottomTwo() const{ return getAntisdownSix(); }
      MassEigenstate&
      getSdownSix();
      MassEigenstate const&
      getSdownSix() const;
      MassEigenstate&
      getSbottomTwo(){ return getSdownSix(); }
      MassEigenstate const&
      getSbottomTwo() const{ return getSdownSix(); }
      MassEigenstate&
      getSupOne();
      MassEigenstate const&
      getSupOne() const;
      MassEigenstate&
      getSupL(){ return getSupOne(); }
      MassEigenstate const&
      getSupL() const{ return getSupOne(); }
      MassEigenstate&
      getAntisupOne();
      MassEigenstate const&
      getAntisupOne() const;
      MassEigenstate&
      getAntisupL(){ return getAntisupOne(); }
      MassEigenstate const&
      getAntisupL() const{ return getAntisupOne(); }
      MassEigenstate&
      getSupTwo();
      MassEigenstate const&
      getSupTwo() const;
      MassEigenstate&
      getScharmL(){ return getSupTwo(); }
      MassEigenstate const&
      getScharmL() const{ return getSupTwo(); }
      MassEigenstate&
      getAntisupTwo();
      MassEigenstate const&
      getAntisupTwo() const;
      MassEigenstate&
      getAntischarmL(){ return getAntisupTwo(); }
      MassEigenstate const&
      getAntischarmL() const{ return getAntisupTwo(); }
      MassEigenstate&
      getSupThree();
      MassEigenstate const&
      getSupThree() const;
      MassEigenstate&
      getStopOne(){ return getSupThree(); }
      MassEigenstate const&
      getStopOne() const{ return getSupThree(); }
      MassEigenstate&
      getAntisupThree();
      MassEigenstate const&
      getAntisupThree() const;
      MassEigenstate&
      getAntistopOne(){ return getAntisupThree(); }
      MassEigenstate const&
      getAntistopOne() const{ return getAntisupThree(); }
      MassEigenstate&
      getSupFour();
      MassEigenstate const&
      getSupFour() const;
      MassEigenstate&
      getSupR(){ return getSupFour(); }
      MassEigenstate const&
      getSupR() const{ return getSupFour(); }
      MassEigenstate&
      getAntisupFour();
      MassEigenstate const&
      getAntisupFour() const;
      MassEigenstate&
      getAntisupR(){ return getAntisupFour(); }
      MassEigenstate const&
      getAntisupR() const{ return getAntisupFour(); }
      MassEigenstate&
      getSupFive();
      MassEigenstate const&
      getSupFive() const;
      MassEigenstate&
      getScharmR(){ return getSupFive(); }
      MassEigenstate const&
      getScharmR() const{ return getSupFive(); }
      MassEigenstate&
      getAntisupFive();
      MassEigenstate const&
      getAntisupFive() const;
      MassEigenstate&
      getAntischarmR(){ return getAntisupFive(); }
      MassEigenstate const&
      getAntischarmR() const{ return getAntisupFive(); }
      MassEigenstate&
      getSupSix();
      MassEigenstate const&
      getSupSix() const;
      MassEigenstate&
      getStopTwo(){ return getSupSix(); }
      MassEigenstate const&
      getStopTwo() const{ return getSupSix(); }
      MassEigenstate&
      getAntisupSix();
      MassEigenstate const&
      getAntisupSix() const;
      MassEigenstate&
      getAntistopTwo(){ return getAntisupSix(); }
      MassEigenstate const&
      getAntistopTwo() const{ return getAntisupSix(); }
      std::vector< MassEigenstate* >&
      getDownAntisquarks();
      std::vector< MassEigenstate* > const&
      getDownAntisquarks() const;
      std::vector< MassEigenstate* >&
      getDownSquarks();
      std::vector< MassEigenstate* > const&
      getDownSquarks() const;
      std::vector< MassEigenstate* >&
      getUpSquarks();
      std::vector< MassEigenstate* > const&
      getUpSquarks() const;
      std::vector< MassEigenstate* >&
      getUpAntisquarks();
      std::vector< MassEigenstate* > const&
      getUpAntisquarks() const;
      std::vector< MassEigenstate* >&
      getSquarkPointers();
      std::vector< MassEigenstate* > const&
      getSquarkPointers() const;
      std::vector< MassEigenstate* >&
      getAntisquarkPointers();
      std::vector< MassEigenstate* > const&
      getAntisquarkPointers() const;


    protected:
      MassEigenstate antisdownOne;
      MassEigenstate sdownOne;
      MassEigenstate antisdownTwo;
      MassEigenstate sdownTwo;
      MassEigenstate antisdownThree;
      MassEigenstate sdownThree;
      MassEigenstate antisdownFour;
      MassEigenstate sdownFour;
      MassEigenstate antisdownFive;
      MassEigenstate sdownFive;
      MassEigenstate antisdownSix;
      MassEigenstate sdownSix;
      MassEigenstate supOne;
      MassEigenstate antisupOne;
      MassEigenstate supTwo;
      MassEigenstate antisupTwo;
      MassEigenstate supThree;
      MassEigenstate antisupThree;
      MassEigenstate supFour;
      MassEigenstate antisupFour;
      MassEigenstate supFive;
      MassEigenstate antisupFive;
      MassEigenstate supSix;
      MassEigenstate antisupSix;
      std::vector< MassEigenstate* > downAntisquarkPointers;
      std::vector< MassEigenstate* > downSquarkPointers;
      std::vector< MassEigenstate* > upSquarkPointers;
      std::vector< MassEigenstate* > upAntisquarkPointers;
      std::vector< MassEigenstate* > squarkPointers;
      std::vector< MassEigenstate* > antisquarkPointers;
    };



    inline MassEigenstate&
    SquarksOneToSix::getAntisdownOne()
    {
      return antisdownOne;
    }

    inline MassEigenstate const&
    SquarksOneToSix::getAntisdownOne() const
    {
      return antisdownOne;
    }

    inline MassEigenstate&
    SquarksOneToSix::getSdownOne()
    {
      return sdownOne;
    }

    inline MassEigenstate const&
    SquarksOneToSix::getSdownOne() const
    {
      return sdownOne;
    }

    inline MassEigenstate&
    SquarksOneToSix::getAntisdownTwo()
    {
      return antisdownTwo;
    }

    inline MassEigenstate const&
    SquarksOneToSix::getAntisdownTwo() const
    {
      return antisdownTwo;
    }

    inline MassEigenstate&
    SquarksOneToSix::getSdownTwo()
    {
      return sdownTwo;
    }

    inline MassEigenstate const&
    SquarksOneToSix::getSdownTwo() const
    {
      return sdownTwo;
    }

    inline MassEigenstate&
    SquarksOneToSix::getAntisdownThree()
    {
      return antisdownThree;
    }

    inline MassEigenstate const&
    SquarksOneToSix::getAntisdownThree() const
    {
      return antisdownThree;
    }

    inline MassEigenstate&
    SquarksOneToSix::getSdownThree()
    {
      return sdownThree;
    }

    inline MassEigenstate const&
    SquarksOneToSix::getSdownThree() const
    {
      return sdownThree;
    }

    inline MassEigenstate&
    SquarksOneToSix::getAntisdownFour()
    {
      return antisdownFour;
    }

    inline MassEigenstate const&
    SquarksOneToSix::getAntisdownFour() const
    {
      return antisdownFour;
    }

    inline MassEigenstate&
    SquarksOneToSix::getSdownFour()
    {
      return sdownFour;
    }

    inline MassEigenstate const&
    SquarksOneToSix::getSdownFour() const
    {
      return sdownFour;
    }

    inline MassEigenstate&
    SquarksOneToSix::getAntisdownFive()
    {
      return antisdownFive;
    }

    inline MassEigenstate const&
    SquarksOneToSix::getAntisdownFive() const
    {
      return antisdownFive;
    }

    inline MassEigenstate&
    SquarksOneToSix::getSdownFive()
    {
      return sdownFive;
    }

    inline MassEigenstate const&
    SquarksOneToSix::getSdownFive() const
    {
      return sdownFive;
    }

    inline MassEigenstate&
    SquarksOneToSix::getAntisdownSix()
    {
      return antisdownSix;
    }

    inline MassEigenstate const&
    SquarksOneToSix::getAntisdownSix() const
    {
      return antisdownSix;
    }

    inline MassEigenstate&
    SquarksOneToSix::getSdownSix()
    {
      return sdownSix;
    }

    inline MassEigenstate const&
    SquarksOneToSix::getSdownSix() const
    {
      return sdownSix;
    }

    inline MassEigenstate&
    SquarksOneToSix::getSupOne()
    {
      return supOne;
    }

    inline MassEigenstate const&
    SquarksOneToSix::getSupOne() const
    {
      return supOne;
    }

    inline MassEigenstate&
    SquarksOneToSix::getAntisupOne()
    {
      return antisupOne;
    }

    inline MassEigenstate const&
    SquarksOneToSix::getAntisupOne() const
    {
      return antisupOne;
    }

    inline MassEigenstate&
    SquarksOneToSix::getSupTwo()
    {
      return supTwo;
    }

    inline MassEigenstate const&
    SquarksOneToSix::getSupTwo() const
    {
      return supTwo;
    }

    inline MassEigenstate&
    SquarksOneToSix::getAntisupTwo()
    {
      return antisupTwo;
    }

    inline MassEigenstate const&
    SquarksOneToSix::getAntisupTwo() const
    {
      return antisupTwo;
    }

    inline MassEigenstate&
    SquarksOneToSix::getSupThree()
    {
      return supThree;
    }

    inline MassEigenstate const&
    SquarksOneToSix::getSupThree() const
    {
      return supThree;
    }

    inline MassEigenstate&
    SquarksOneToSix::getAntisupThree()
    {
      return antisupThree;
    }

    inline MassEigenstate const&
    SquarksOneToSix::getAntisupThree() const
    {
      return antisupThree;
    }

    inline MassEigenstate&
    SquarksOneToSix::getSupFour()
    {
      return supFour;
    }

    inline MassEigenstate const&
    SquarksOneToSix::getSupFour() const
    {
      return supFour;
    }

    inline MassEigenstate&
    SquarksOneToSix::getAntisupFour()
    {
      return antisupFour;
    }

    inline MassEigenstate const&
    SquarksOneToSix::getAntisupFour() const
    {
      return antisupFour;
    }

    inline MassEigenstate&
    SquarksOneToSix::getSupFive()
    {
      return supFive;
    }

    inline MassEigenstate const&
    SquarksOneToSix::getSupFive() const
    {
      return supFive;
    }

    inline MassEigenstate&
    SquarksOneToSix::getAntisupFive()
    {
      return antisupFive;
    }

    inline MassEigenstate const&
    SquarksOneToSix::getAntisupFive() const
    {
      return antisupFive;
    }

    inline MassEigenstate&
    SquarksOneToSix::getSupSix()
    {
      return supSix;
    }

    inline MassEigenstate const&
    SquarksOneToSix::getSupSix() const
    {
      return supSix;
    }

    inline MassEigenstate&
    SquarksOneToSix::getAntisupSix()
    {
      return antisupSix;
    }

    inline MassEigenstate const&
    SquarksOneToSix::getAntisupSix() const
    {
      return antisupSix;
    }

    inline std::vector< MassEigenstate* >&
    SquarksOneToSix::getDownAntisquarks()
    {
      return downAntisquarkPointers;
    }

    inline std::vector< MassEigenstate* > const&
    SquarksOneToSix::getDownAntisquarks() const
    {
      return downAntisquarkPointers;
    }

    inline std::vector< MassEigenstate* >&
    SquarksOneToSix::getDownSquarks()
    {
      return downSquarkPointers;
    }

    inline std::vector< MassEigenstate* > const&
    SquarksOneToSix::getDownSquarks() const
    {
      return downSquarkPointers;
    }

    inline std::vector< MassEigenstate* >&
    SquarksOneToSix::getUpSquarks()
    {
      return upSquarkPointers;
    }

    inline std::vector< MassEigenstate* > const&
    SquarksOneToSix::getUpSquarks() const
    {
      return upSquarkPointers;
    }

    inline std::vector< MassEigenstate* >&
    SquarksOneToSix::getUpAntisquarks()
    {
      return upAntisquarkPointers;
    }

    inline std::vector< MassEigenstate* > const&
    SquarksOneToSix::getUpAntisquarks() const
    {
      return upAntisquarkPointers;
    }

    inline std::vector< MassEigenstate* >&
    SquarksOneToSix::getSquarkPointers()
    {
      return squarkPointers;
    }

    inline std::vector< MassEigenstate* > const&
    SquarksOneToSix::getSquarkPointers() const
    {
      return squarkPointers;
    }

    inline std::vector< MassEigenstate* >&
    SquarksOneToSix::getAntisquarkPointers()
    {
      return antisquarkPointers;
    }

    inline std::vector< MassEigenstate* > const&
    SquarksOneToSix::getAntisquarkPointers() const
    {
      return antisquarkPointers;
    }

  }

}

#endif /* SQUARKSONETOSIX_HPP_ */
