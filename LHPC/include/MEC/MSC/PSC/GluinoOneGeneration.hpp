/*
 * GluinoOneGeneration.hpp
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

#ifndef GLUINOONEGENERATION_HPP_
#define GLUINOONEGENERATION_HPP_

#include "../CodesAndDataForMassEigenstates.hpp"

namespace LHPC
{
  namespace MassSpectrumClass
  {
    class GluinoOneGeneration : public virtual MassSpectrum
    {
    public:
      GluinoOneGeneration( bool const isVerbose = false,
                           std::vector< bool >* const defaultFlags = NULL );
      virtual
      ~GluinoOneGeneration();

      MassEigenstate&
      getGluino();
      MassEigenstate const&
      getGluino() const;


    protected:
      MassEigenstate gluinoOne;
    };



    inline MassEigenstate&
    GluinoOneGeneration::getGluino()
    {
      return gluinoOne;
    }

    inline MassEigenstate const&
    GluinoOneGeneration::getGluino() const
    {
      return gluinoOne;
    }

  }

}

#endif /* GLUINOONEGENERATION_HPP_ */
