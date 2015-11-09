/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SUSYMASSINSERTIONMATCHING_H
#define	SUSYMASSINSERTIONMATCHING_H

#include <gslpp.h>
#include <gsl/gsl_sf_dilog.h>
#include <stdexcept>
#include <StandardModelMatching.h>

class SUSYMassInsertion;

    /**
     * @class SUSYMassInsertionMatching class
     * @ingroup SUSYMassInsertion
     * @brief A class for the matching of SUSY MIA.
     * @author HEPfit Collaboration
     * @copyright GNU General Public License
     * @details this class contains the Wilson Coefficients at the SUSY matching scale
     * for \f$ \Delta F = 1 \f$ and \f$ \Delta F = 2 \f$ processes in the mass insertion approximation;
     * references: 
     * \Delta F = 2 \f$ Wilson coefficients: M. Ciuchini et al, hep-ph/0606197
     * \Delta F = 1 \f$ Wilson coefficients: F.Gabbiani et al, hep-ph/9604387
     * @param SUSYMassInsertion_i object of type SUSYMassInsertion
     */
class SUSYMassInsertionMatching: public StandardModelMatching{
public :
    /**
     * @brief constructor
     * @param SUSYMassInsertion_i object of type SUSYMassInsertion
     */
    SUSYMassInsertionMatching(const SUSYMassInsertion& SUSYMassInsertion_i);
    
    
    //delta F=2 loop functions
    
    double f6(double x);
    double f6t(double x);
    
    /**
     * 
     * @brief \f$ \DeltaF = 2 \f$ loop functions RI-MOM scheme, LO term
     * @param x the square ratio between squark mass and gluino mass
     */
    double C0LO(double x);
    
    /**
     * 
     * @brief \f$ \DeltaF = 2 \f$ loop functions RI-MOM scheme, LO term
     * @param x the square ratio between squark mass and gluino mass
     */
    double C1LO(double x);
    
    /**
     * 
     * @brief \f$ \DeltaF = 2 \f$ loop functions RI-MOM scheme, LO term
     * @param x the square ratio between squark mass and gluino mass
     */
    double C2LO(double x);
    
    /**
     * 
     * @brief \f$ \DeltaF = 2 \f$ loop functions RI-MOM scheme, LO term
     * @param x the square ratio between squark mass and gluino mass
     */
    double C3LOA(double x);
    
    /**
     * 
     * @brief \f$ \DeltaF = 2 \f$ loop functions RI-MOM scheme, LO term
     * @param x the square ratio between squark mass and gluino mass
     */
    double C3LOB(double x);
    
    /**
     * 
     * @brief \f$ \DeltaF = 2 \f$ loop functions RI-MOM scheme, LO term
     * @param x the square ratio between squark mass and gluino mass
     */
    double C4LOA(double x);
    
    /**
     * 
     * @brief \f$ \DeltaF = 2 \f$ loop functions RI-MOM scheme, LO term
     * @param x the square ratio between squark mass and gluino mass
     */
    double C4LOB(double x);
    
    /**
     * 
     * @brief \f$ \DeltaF = 2 \f$ loop functions RI-MOM scheme, LO term
     * @param x the square ratio between squark mass and gluino mass
     * @param mumatch2 double for the square of the SUSY matching scale
     * @param Ms2 double for the square of the average squark mass 
     */
    double C0NLO(double x, double mumatch2, double Ms2);
    
    /**
     * 
     * @brief \f$ \DeltaF = 2 \f$ loop functions RI-MOM scheme, LO term
     * @param x the square ratio between squark mass and gluino mass
     * @param mumatch2 double for the square of the SUSY matching scale
     * @param Ms2 double for the square of the average squark mass 
     */
    double C1NLO(double x, double mumatch2, double Ms2);
    
    /**
     * 
     * @brief \f$ \DeltaF = 2 \f$ loop functions RI-MOM scheme, LO term
     * @param x the square ratio between squark mass and gluino mass
     * @param mumatch2 double for the square of the SUSY matching scale
     * @param Ms2 double for the square of the average squark mass 
     */
    double C2NLO(double x, double mumatch2, double Ms2);
    
    /**
     * 
     * @brief \f$ \DeltaF = 2 \f$ loop functions RI-MOM scheme, LO term
     * @param x the square ratio between squark mass and gluino mass
     * @param mumatch2 double for the square of the SUSY matching scale
     * @param Ms2 double for the square of the average squark mass 
     */
    double C3NLOA(double x, double mumatch2, double Ms2);
    
    /**
     * 
     * @brief \f$ \DeltaF = 2 \f$ loop functions RI-MOM scheme, LO term
     * @param x the square ratio between squark mass and gluino mass
     * @param mumatch2 double for the square of the SUSY matching scale
     * @param Ms2 double for the square of the average squark mass 
     */
    double C3NLOB(double x, double mumatch2, double Ms2);
    
    /**
     * 
     * @brief \f$ \DeltaF = 2 \f$ loop functions RI-MOM scheme, LO term
     * @param x the square ratio between squark mass and gluino mass
     * @param mumatch2 double for the square of the SUSY matching scale
     * @param Ms2 double for the square of the average squark mass 
     */
    double C4NLOA(double x, double mumatch2, double Ms2);
    
    /**
     * 
     * @brief \f$ \DeltaF = 2 \f$ loop functions RI-MOM scheme, LO term
     * @param x the square ratio between squark mass and gluino mass
     * @param mumatch2 double for the square of the SUSY matching scale
     * @param Ms2 double for the square of the average squark mass 
     */
    double C4NLOB(double x, double mumatch2, double Ms2);
        
    //delta F=1 loop functions
    
    /**
     * 
     * @brief \f$ \DeltaF = 1 \f$ loop functions, LO term
     * @param x the square ratio between squark mass and gluino mass
     */
    double B1(double x) const;
    
    /**
     * 
     * @brief \f$ \DeltaF = 1 \f$ loop functions, LO term
     * @param x the square ratio between squark mass and gluino mass
     */
    double B2(double x) const;
    
    /**
     * 
     * @brief \f$ \DeltaF = 1 \f$ loop functions, LO term
     * @param x the square ratio between squark mass and gluino mass
     */
    double P1(double x) const;
    
    /**
     * 
     * @brief \f$ \DeltaF = 1 \f$ loop functions, LO term
     * @param x the square ratio between squark mass and gluino mass
     */
    double P2(double x) const;
    
    /**
     * 
     * @brief \f$ \DeltaF = 1 \f$ loop functions, LO term
     * @param x the square ratio between squark mass and gluino mass
     */
    double M1(double x) const;
    
    /**
     * 
     * @brief \f$ \DeltaF = 1 \f$ loop functions, LO term
     * @param x the square ratio between squark mass and gluino mass
     */
    double M2(double x) const;
    
    /**
     * 
     * @brief \f$ \DeltaF = 1 \f$ loop functions, LO term
     * @param x the square ratio between squark mass and gluino mass
     */
    double M3(double x) const;
    
    /**
     * 
     * @brief \f$ \DeltaF = 1 \f$ loop functions, LO term
     * @param x the square ratio between squark mass and gluino mass
     */
    double M4(double x) const;
    
    
    // Wilson coefficients DeltaF=2
    
    /**
     * 
     * @brief \f$ \Delta C = 2 \f$
     * @return the vector of 8 wilson coefficients: SM + SusyMI
     */
    virtual  std::vector<WilsonCoefficient>& CMdd2();
    
    /**
     * 
     * @brief \f$ \Delta B = 2 \f$, \f$ B_{d} \f$
     * @return the vector of 8 wilson coefficients: SM + SusyMI
     */
    virtual  std::vector<WilsonCoefficient>& CMdbd2();
    
    /**
     * 
     * @brief \f$ \Delta B = 2 \f$, \f$ B_{s} \f$
     * @return the vector of 8 wilson coefficients: SM + SusyMI
     */
    virtual  std::vector<WilsonCoefficient>& CMdbs2();
    
    /**
     * 
     * @brief \f$ \Delta S = 2 \f$
     * @return the vector of 8 wilson coefficients: SM + SusyMI
     */
    virtual  std::vector<WilsonCoefficient>& CMdk2();
    
    /**
     * 
     * @brief \f$ \Delta F = 1 \f$
     * @return the vector of 8 wilson coefficients: SM + SusyMI
     */
    virtual  std::vector<WilsonCoefficient>& CMd1();
    
    /**
     * @brief Misiak et al, hep-ph/0005183
     * @return matrix for the change of basis from the standard to Misiak one
     */
    gslpp::matrix<double> RtoMisiak() const;
    
    /**
     * 
     * @brief the basis is defined, for example, in Chetyrkin et al hep-ph/9612313
     * @return matrix for the change to effective basis
     */
    gslpp::matrix<double> EffectiveBase() const;
    
    /**
     * 
     * @brief it changes renormalization scheme from LRI to NDR
     * @param i int flag: 1 for D, 2 for B_d, 3 for B_s, 4 for K
     */
    void LRItoNDR (int i);
    
private:
    const SUSYMassInsertion& SusyMI;
    double MuM2, Ms2;
    unsigned int Nf;
    gslpp::complex DLL, DLR, DRL, DRR;
    gslpp::matrix<double> drNDRLRI;
    
    WilsonCoefficient mcd2, mcd1, mcbd, mcbs, mck2;
    
};

#endif	/* SUSYMASSINSERTIONMATCHING_H */

