/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef AMPDS1_H
#define	AMPDS1_H

class StandardModel;
#include "gslpp_complex.h"
#include "OrderScheme.h"
#include "gslpp.h"
#include "math.h"

class AmpDS1 {
public:
    /**
     *
     * @brief compute the amplitude for \f$ K_L \f$ decay in 2 pion
     * @param Flavour
     */
    AmpDS1(const StandardModel& SM_i);

protected:
    /**
     *
     * @param order
     * @return the amplitude for \f$ K_L \f$ decay in 2 pion with 0 isospin change
     */
    gslpp::complex AmpDS1pp0(orders order);
    
     /**
     *
     * @param order
     * @return the amplitude (pure lattice) for \f$ K_L \f$ decay in 2 pion with 0 isospin change
     */
    gslpp::complex AmpDS1pp0pureLAT(orders order);

    /**
     *
     * @param order
     * @return the amplitude for \f$ K_L \f$ decay in 2 pion with double isospin change
     */
    gslpp::complex AmpDS1pp2(orders order);

    /**
     *
     * @return the real part of the amplitude for \f$ K_L \f$ decay in 2 pion with 0 isospin change
     * from experimental inputs
     */
     double getReA0();

     /**
      *
      * @return the real part of the amplitude for \f$ K_L \f$ decay in 2 pion with double
      * isospin change from experimental inputs
      */
      double getReA2();

      /**
       *
       * @return transformation matrix for the matrix elements in the chiral basis to the
       * ten operator basis for the isospin zero channel (cfr. eqn. 59 ArXiv:1104.4948)
       */
      gslpp::matrix<double> getChiralMatrixpp0() const;

      /**
       *
       * @return transformation matrix for the matrix elements in the chiral basis to the
       * ten operator basis for the isospin two channel (cfr. eqn. 69-70 ArXiv:1502.00263)
       */
      gslpp::matrix<double> getChiralMatrixpp2() const;

      /**
       *
       * @return renormalization matrix for the lattice matrix elements in the chiral basis
       * to the RI-SMOM scheme in the chiral basis
       */
      gslpp::matrix<double> getRIMatrixpp0() const;

      /**
       *
       * @return transformation matrix for the matrix elements in the chiral basis in the RI-SMOM
       * scheme to the MSbar scheme in the chiral basis (cfr. ArXiv:1104.4948)
       *
       */
      gslpp::matrix<double> getRISMOMTransMatrix(double mu, orders order) const;

private:

    const StandardModel& mySM;

};


#endif	/* AMPDS1_H */
