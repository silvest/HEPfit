/*
 * ChargedSleptonsOneToSix.hpp
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

#ifndef CHARGEDSLEPTONSONETOSIX_HPP_
#define CHARGEDSLEPTONSONETOSIX_HPP_

#include "../CodesAndDataForMassEigenstates.hpp"

namespace LHPC
{
  namespace MassSpectrumClass
  {
    class ChargedSleptonsOneToSix : public virtual MassSpectrum
    {
    public:
      ChargedSleptonsOneToSix( bool const isVerbose = false,
                               bool const flavorConserving = false,
                              std::vector< bool >* const defaultFlags = NULL );
      virtual
      ~ChargedSleptonsOneToSix();

      MassEigenstate&
      getSpositronOne();
      MassEigenstate const&
      getSpositronOne() const;
      MassEigenstate&
      getAntiselectronOne(){ return getSpositronOne(); }
      MassEigenstate const&
      getAntiselectronOne() const{ return getSpositronOne(); }
      MassEigenstate&
      getSpositronL(){ return getSpositronOne(); }
      MassEigenstate const&
      getSpositronL() const{ return getSpositronOne(); }
      MassEigenstate&
      getAntiselectronL(){ return getSpositronOne(); }
      MassEigenstate const&
      getAntiselectronL() const{ return getSpositronOne(); }
      MassEigenstate&
      getAntispositronOne();
      MassEigenstate const&
      getAntispositronOne() const;
      MassEigenstate&
      getSelectronOne(){ return getAntispositronOne(); }
      MassEigenstate const&
      getSelectronOne() const{ return getAntispositronOne(); }
      MassEigenstate&
      getSelectronL(){ return getAntispositronOne(); }
      MassEigenstate const&
      getSelectronL() const{ return getAntispositronOne(); }
      MassEigenstate&
      getSpositronTwo();
      MassEigenstate const&
      getSpositronTwo() const;
      MassEigenstate&
      getAntiselectronTwo(){ return getSpositronTwo(); }
      MassEigenstate const&
      getAntiselectronTwo() const{ return getSpositronTwo(); }
      MassEigenstate&
      getAntismuonL(){ return getSpositronTwo(); }
      MassEigenstate const&
      getAntismuonL() const{ return getSpositronTwo(); }
      MassEigenstate&
      getAntispositronTwo();
      MassEigenstate const&
      getAntispositronTwo() const;
      MassEigenstate&
      getSelectronTwo(){ return getAntispositronTwo(); }
      MassEigenstate const&
      getSelectronTwo() const{ return getAntispositronTwo(); }
      MassEigenstate&
      getSmuonL(){ return getAntispositronTwo(); }
      MassEigenstate const&
      getSmuonL() const{ return getAntispositronTwo(); }
      MassEigenstate&
      getSpositronThree();
      MassEigenstate const&
      getSpositronThree() const;
      MassEigenstate&
      getAntiselectronThree(){ return getSpositronThree(); }
      MassEigenstate const&
      getAntiselectronThree() const{ return getSpositronThree(); }
      MassEigenstate&
      getAntistauOne(){ return getSpositronThree(); }
      MassEigenstate const&
      getAntistauOne() const{ return getSpositronThree(); }
      MassEigenstate&
      getAntispositronThree();
      MassEigenstate const&
      getAntispositronThree() const;
      MassEigenstate&
      getSelectronThree(){ return getAntispositronThree(); }
      MassEigenstate const&
      getSelectronThree() const{ return getAntispositronThree(); }
      MassEigenstate&
      getStauOne(){ return getAntispositronThree(); }
      MassEigenstate const&
      getStauOne() const{ return getAntispositronThree(); }
      MassEigenstate&
      getSpositronFour();
      MassEigenstate const&
      getSpositronFour() const;
      MassEigenstate&
      getAntiselectronFour(){ return getSpositronFour(); }
      MassEigenstate const&
      getAntiselectronFour() const{ return getSpositronFour(); }
      MassEigenstate&
      getSpositronR(){ return getSpositronFour(); }
      MassEigenstate const&
      getSpositronR() const{ return getSpositronFour(); }
      MassEigenstate&
      getAntiselectronR(){ return getSpositronFour(); }
      MassEigenstate const&
      getAntiselectronR() const{ return getSpositronFour(); }
      MassEigenstate&
      getAntispositronFour();
      MassEigenstate const&
      getAntispositronFour() const;
      MassEigenstate&
      getSelectronFour(){ return getAntispositronFour(); }
      MassEigenstate const&
      getSelectronFour() const{ return getAntispositronFour(); }
      MassEigenstate&
      getSelectronR(){ return getAntispositronFour(); }
      MassEigenstate const&
      getSelectronR() const{ return getAntispositronFour(); }
      MassEigenstate&
      getSpositronFive();
      MassEigenstate const&
      getSpositronFive() const;
      MassEigenstate&
      getAntiselectronFive(){ return getSpositronFive(); }
      MassEigenstate const&
      getAntiselectronFive() const{ return getSpositronFive(); }
      MassEigenstate&
      getAntismuonR(){ return getSpositronFive(); }
      MassEigenstate const&
      getAntismuonR() const{ return getSpositronFive(); }
      MassEigenstate&
      getAntispositronFive();
      MassEigenstate const&
      getAntispositronFive() const;
      MassEigenstate&
      getSelectronFive(){ return getAntispositronFive(); }
      MassEigenstate const&
      getSelectronFive() const{ return getAntispositronFive(); }
      MassEigenstate&
      getSmuonR(){ return getAntispositronFive(); }
      MassEigenstate const&
      getSmuonR() const{ return getAntispositronFive(); }
      MassEigenstate&
      getSpositronSix();
      MassEigenstate const&
      getSpositronSix() const;
      MassEigenstate&
      getAntiselectronSix(){ return getSpositronSix(); }
      MassEigenstate const&
      getAntiselectronSix() const{ return getSpositronSix(); }
      MassEigenstate&
      getAntistauTwo(){ return getSpositronSix(); }
      MassEigenstate const&
      getAntistauTwo() const{ return getSpositronSix(); }
      MassEigenstate&
      getAntispositronSix();
      MassEigenstate const&
      getAntispositronSix() const;
      MassEigenstate&
      getSelectronSix(){ return getAntispositronSix(); }
      MassEigenstate const&
      getSelectronSix() const{ return getAntispositronSix(); }
      MassEigenstate&
      getStauTwo(){ return getAntispositronSix(); }
      MassEigenstate const&
      getStauTwo() const{ return getAntispositronSix(); }
      std::vector< MassEigenstate* >&
      getPositiveSleptons();
      std::vector< MassEigenstate* > const&
      getPositiveSleptons() const;
      std::vector< MassEigenstate* >&
      getNegativeSleptons();
      std::vector< MassEigenstate* > const&
      getNegativeSleptons() const;


    protected:
      MassEigenstate spositronOne;
      MassEigenstate antispositronOne;
      MassEigenstate spositronTwo;
      MassEigenstate antispositronTwo;
      MassEigenstate spositronThree;
      MassEigenstate antispositronThree;
      MassEigenstate spositronFour;
      MassEigenstate antispositronFour;
      MassEigenstate spositronFive;
      MassEigenstate antispositronFive;
      MassEigenstate spositronSix;
      MassEigenstate antispositronSix;
      std::vector< MassEigenstate* > positiveSleptonPointers;
      std::vector< MassEigenstate* > negativeSleptonPointers;
    };



    inline MassEigenstate&
    ChargedSleptonsOneToSix::getSpositronOne()
    {
      return spositronOne;
    }

    inline MassEigenstate const&
    ChargedSleptonsOneToSix::getSpositronOne() const
    {
      return spositronOne;
    }

    inline MassEigenstate&
    ChargedSleptonsOneToSix::getAntispositronOne()
    {
      return antispositronOne;
    }

    inline MassEigenstate const&
    ChargedSleptonsOneToSix::getAntispositronOne() const
    {
      return antispositronOne;
    }

    inline MassEigenstate&
    ChargedSleptonsOneToSix::getSpositronTwo()
    {
      return spositronTwo;
    }

    inline MassEigenstate const&
    ChargedSleptonsOneToSix::getSpositronTwo() const
    {
      return spositronTwo;
    }

    inline MassEigenstate&
    ChargedSleptonsOneToSix::getAntispositronTwo()
    {
      return antispositronTwo;
    }

    inline MassEigenstate const&
    ChargedSleptonsOneToSix::getAntispositronTwo() const
    {
      return antispositronTwo;
    }

    inline MassEigenstate&
    ChargedSleptonsOneToSix::getSpositronThree()
    {
      return spositronThree;
    }

    inline MassEigenstate const&
    ChargedSleptonsOneToSix::getSpositronThree() const
    {
      return spositronThree;
    }

    inline MassEigenstate&
    ChargedSleptonsOneToSix::getAntispositronThree()
    {
      return antispositronThree;
    }

    inline MassEigenstate const&
    ChargedSleptonsOneToSix::getAntispositronThree() const
    {
      return antispositronThree;
    }

    inline MassEigenstate&
    ChargedSleptonsOneToSix::getSpositronFour()
    {
      return spositronFour;
    }

    inline MassEigenstate const&
    ChargedSleptonsOneToSix::getSpositronFour() const
    {
      return spositronFour;
    }

    inline MassEigenstate&
    ChargedSleptonsOneToSix::getAntispositronFour()
    {
      return antispositronFour;
    }

    inline MassEigenstate const&
    ChargedSleptonsOneToSix::getAntispositronFour() const
    {
      return antispositronFour;
    }

    inline MassEigenstate&
    ChargedSleptonsOneToSix::getSpositronFive()
    {
      return spositronFive;
    }

    inline MassEigenstate const&
    ChargedSleptonsOneToSix::getSpositronFive() const
    {
      return spositronFive;
    }

    inline MassEigenstate&
    ChargedSleptonsOneToSix::getAntispositronFive()
    {
      return antispositronFive;
    }

    inline MassEigenstate const&
    ChargedSleptonsOneToSix::getAntispositronFive() const
    {
      return antispositronFive;
    }

    inline MassEigenstate&
    ChargedSleptonsOneToSix::getSpositronSix()
    {
      return spositronSix;
    }

    inline MassEigenstate const&
    ChargedSleptonsOneToSix::getSpositronSix() const
    {
      return spositronSix;
    }

    inline MassEigenstate&
    ChargedSleptonsOneToSix::getAntispositronSix()
    {
      return antispositronSix;
    }

    inline MassEigenstate const&
    ChargedSleptonsOneToSix::getAntispositronSix() const
    {
      return antispositronSix;
    }

    inline std::vector< MassEigenstate* >&
    ChargedSleptonsOneToSix::getPositiveSleptons()
    {
      return positiveSleptonPointers;
    }

    inline std::vector< MassEigenstate* > const&
    ChargedSleptonsOneToSix::getPositiveSleptons() const
    {
      return positiveSleptonPointers;
    }

    inline std::vector< MassEigenstate* >&
    ChargedSleptonsOneToSix::getNegativeSleptons()
    {
      return negativeSleptonPointers;
    }

    inline std::vector< MassEigenstate* > const&
    ChargedSleptonsOneToSix::getNegativeSleptons() const
    {
      return negativeSleptonPointers;
    }

  }

}

#endif /* CHARGEDSLEPTONSONETOSIX_HPP_ */
