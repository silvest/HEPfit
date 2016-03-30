/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef UNITARITY_H
#define	UNITARITY_H

#include <stdexcept>
#include "ThObservable.h"
#include "THDM.h"
#include "THDMcache.h"

/**
 * @class unitarity
 * @ingroup THDM 
 * @brief An observable class for the requirement of tree level perturbative 
 * unitarity.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to require the unitarity for all the tree level 
 * scalar-scalar scattering amplitudes.
 * The 12 eigenvalues of the S-matrix can be found in @cite Akeroyd:2000wc,Ginzburg:2005dt.
 * They should be smaller than @f$8\pi@f$ in magnitude to preserve unitarity of the S-matrix,
 * however corresponding two-loop calculations in the SM show 
 * that a more reasonable upper bound on their absolute values would be @f$2\pi@f$
 * (see discussion in @cite Baglio:2014nea).
 */
class unitarity : public ThObservable {
public:
    /**
     * @brief Constructor.
     * @param[in] ?
     */
   unitarity(const StandardModel& SM_i);
     
    /**
     * @brief An empty constructor.
     * @return
     */
    double computeThValue();

    const THDM * myTHDM;
    
};

/**
 * @class unitarity1
 * @ingroup THDM 
 * @brief Eigenvalue of the @f$\Phi \Phi \to \Phi \Phi@f$ scattering matrix.
 * @details @f$\Lambda^{even}_{21+}@f$.
 */
class unitarity1: public unitarity {
public:

    /**
     * @brief Constructor for @f$\Lambda^{even}_{21+}@f$ from equations (9) in @cite Ginzburg:2005dt.
     */
    unitarity1(const StandardModel& SM_i);

    /**
     * @return @f$\Lambda^{even}_{21+}=\frac{1}{2}(\lambda_1+\lambda_2+\sqrt{(\lambda_1-\lambda_2)^2+4\lambda_5^2})@f$
     * @details Corresponds to @f$c_+@f$ from equation (2.8) in @cite Akeroyd:2000wc.
     */
    double computeThValue();
};

/**
 * @class unitarity2
 * @ingroup THDM 
 * @brief Eigenvalue of the @f$\Phi \Phi \to \Phi \Phi@f$ scattering matrix.
 * @details @f$\Lambda^{even}_{21-}@f$.
 */
class unitarity2: public unitarity {
public:

    /**
     * @brief Constructor for @f$\Lambda^{even}_{21-}@f$ from equations (9) in @cite Ginzburg:2005dt.
     */
    unitarity2(const StandardModel& SM_i);

    /**
     * @return @f$\Lambda^{even}_{21-}=\frac{1}{2}(\lambda_1+\lambda_2-\sqrt{(\lambda_1-\lambda_2)^2+4\lambda_5^2})@f$
     * @details Corresponds to @f$c_-@f$ from equation (2.8) in @cite Akeroyd:2000wc.
     */
    double computeThValue();
};

/**
 * @class unitarity3
 * @ingroup THDM 
 * @brief Eigenvalue of the @f$\Phi \Phi \to \Phi \Phi@f$ scattering matrix.
 * @details @f$\Lambda^{even}_{01+}@f$.
 */
class unitarity3: public unitarity {
public:

    /**
     * @brief Constructor for @f$\Lambda^{even}_{01+}@f$ from equations (9) in @cite Ginzburg:2005dt.
     */
    unitarity3(const StandardModel& SM_i);

    /**
     * @return @f$\Lambda^{even}_{01+}=\frac{1}{2}(\lambda_1+\lambda_2+\sqrt{(\lambda_1-\lambda_2)^2+4\lambda_4^2})@f$
     * @details Corresponds to @f$b_+@f$ from equation (2.7) in @cite Akeroyd:2000wc.
     */
    double computeThValue();
};

/**
 * @class unitarity4
 * @ingroup THDM 
 * @brief Eigenvalue of the @f$\Phi \Phi \to \Phi \Phi@f$ scattering matrix.
 * @details @f$\Lambda^{even}_{01-}@f$.
 */
class unitarity4: public unitarity {
public:

    /**
     * @brief Constructor for @f$\Lambda^{even}_{01-}@f$ from equations (9) in @cite Ginzburg:2005dt.
     */
    unitarity4(const StandardModel& SM_i);

    /**
     * @return @f$\Lambda^{even}_{01-}=\frac{1}{2}(\lambda_1+\lambda_2-\sqrt{(\lambda_1-\lambda_2)^2+4\lambda_4^2})@f$
     * @details Corresponds to @f$b_-@f$ from equation (2.7) in @cite Akeroyd:2000wc.
     */
    double computeThValue();
};

/**
 * @class unitarity5
 * @ingroup THDM 
 * @brief Eigenvalue of the @f$\Phi \Phi \to \Phi \Phi@f$ scattering matrix.
 * @details @f$\Lambda^{even}_{00+}@f$.
 */
class unitarity5: public unitarity {
public:

    /**
     * @brief Constructor for @f$\Lambda^{even}_{00+}@f$ from equations (9) in @cite Ginzburg:2005dt.
     */
    unitarity5(const StandardModel& SM_i);

    /**
     * @return @f$\Lambda^{even}_{00+}=\frac{3}{2}(\lambda_1+\lambda_2)+\frac{1}{2}\sqrt{9(\lambda_1-\lambda_2)^2+4(2\lambda_3+\lambda_4)^2}@f$
     * @details Corresponds to @f$a_+@f$ from equation (2.6) in @cite Akeroyd:2000wc.
     */
    double computeThValue();
};

/**
 * @class unitarity6
 * @ingroup THDM 
 * @brief Eigenvalue of the @f$\Phi \Phi \to \Phi \Phi@f$ scattering matrix.
 * @details @f$\Lambda^{even}_{00-}@f$.
 */
class unitarity6: public unitarity {
public:

    /**
     * @brief Constructor for @f$\Lambda^{even}_{00-}@f$ from equations (9) in @cite Ginzburg:2005dt.
     */
    unitarity6(const StandardModel& SM_i);

    /**
     * @return @f$\Lambda^{even}_{00-}=\frac{3}{2}(\lambda_1+\lambda_2)-\frac{1}{2}\sqrt{9(\lambda_1-\lambda_2)^2+4(2\lambda_3+\lambda_4)^2}@f$
     * @details Corresponds to @f$a_-@f$ from equation (2.6) in @cite Akeroyd:2000wc.
     */
    double computeThValue();
};

/**
 * @class unitarity7
 * @ingroup THDM 
 * @brief Eigenvalue of the @f$\Phi \Phi \to \Phi \Phi@f$ scattering matrix.
 * @details @f$\Lambda^{odd}_{21}@f$.
 */
class unitarity7: public unitarity {
public:

    /**
     * @brief Constructor for @f$\Lambda^{odd}_{21}@f$ from equations (9) in @cite Ginzburg:2005dt.
     */
    unitarity7(const StandardModel& SM_i);

    /**
     * @return @f$\Lambda^{odd}_{21}=\lambda_3+\lambda_4@f$
     * @details Corresponds to @f$f_1=f_2@f$ from equation (2.5) in @cite Akeroyd:2000wc
     */
    double computeThValue();
};

/**
 * @class unitarity8
 * @ingroup THDM 
 * @brief Eigenvalue of the @f$\Phi \Phi \to \Phi \Phi@f$ scattering matrix.
 * @details @f$\Lambda^{odd}_{20}@f$.
 */
class unitarity8: public unitarity {
public:

    /**
     * @brief Constructor for @f$\Lambda^{odd}_{20}@f$ from equations (9) in @cite Ginzburg:2005dt.
     */
    unitarity8(const StandardModel& SM_i);

    /**
     * @return @f$\Lambda^{odd}_{20}=\lambda_3-\lambda_4@f$
     * @details Corresponds to @f$p_1@f$ from equation (2.9) in @cite Akeroyd:2000wc.
     */
    double computeThValue();
};

/**
 * @class unitarity9
 * @ingroup THDM 
 * @brief Eigenvalue of the @f$\Phi \Phi \to \Phi \Phi@f$ scattering matrix.
 * @details @f$\Lambda^{odd}_{01+}@f$.
 */
class unitarity9: public unitarity {
public:

    /**
     * @brief Constructor for @f$\Lambda^{odd}_{01+}@f$ from equations (9) in @cite Ginzburg:2005dt.
     */
    unitarity9(const StandardModel& SM_i);

    /**
     * @return @f$\Lambda^{odd}_{01+}=\lambda_3+\lambda_5@f$
     * @details Corresponds to @f$f_-@f$ from equations (2.5) in @cite Akeroyd:2000wc.
     */
    double computeThValue();
};

/**
 * @class unitarity10
 * @ingroup THDM 
 * @brief Eigenvalue of the @f$\Phi \Phi \to \Phi \Phi@f$ scattering matrix.
 * @details @f$\Lambda^{odd}_{01-}@f$.
 */
class unitarity10: public unitarity {
public:

    /**
     * @brief Constructor for @f$\Lambda^{odd}_{01-}@f$ from equations (9) in @cite Ginzburg:2005dt.
     */
    unitarity10(const StandardModel& SM_i);

    /**
     * @return @f$\Lambda^{odd}_{01-}=\lambda_3-\lambda_5@f$
     * @details Corresponds to @f$e_2@f$ from equations (2.5) in @cite Akeroyd:2000wc.
     */
    double computeThValue();
};

/**
 * @class unitarity11
 * @ingroup THDM 
 * @brief Eigenvalue of the @f$\Phi \Phi \to \Phi \Phi@f$ scattering matrix.
 * @details @f$\Lambda^{odd}_{00+}@f$.
 */
class unitarity11: public unitarity {
public:

    /**
     * @brief Constructor for @f$\Lambda^{odd}_{00+}@f$ from equations (9) in @cite Ginzburg:2005dt.
     */
    unitarity11(const StandardModel& SM_i);

    /**
     * @return @f$\Lambda^{odd}_{00+}=\lambda_3+2\lambda_4+3\lambda_5@f$
     * @details Corresponds to @f$f_+@f$ from equations (2.5) in @cite Akeroyd:2000wc.
     */
    double computeThValue();
};

/**
 * @class unitarity12
 * @ingroup THDM 
 * @brief Eigenvalue of the @f$\Phi \Phi \to \Phi \Phi@f$ scattering matrix.
 * @details @f$\Lambda^{odd}_{00-}@f$.
 */
class unitarity12: public unitarity {
public:

    /**
     * @brief Constructor for @f$\Lambda^{odd}_{00-}@f$ from equations (9) in @cite Ginzburg:2005dt.
     */
    unitarity12(const StandardModel& SM_i);

    /**
     * @return @f$\Lambda^{odd}_{00-}=\lambda_3+2\lambda_4-3\lambda_5@f$
     * @details Corresponds to @f$e_1@f$ from equations (2.5) in @cite Akeroyd:2000wc.
     */
    double computeThValue();
};

/**
 * @class unitarityNLO1
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO1: public ThObservable {
public:

    /**
     * @brief unitarityNLO1 constructor.
     */
    unitarityNLO1(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO1@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO2
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO2: public ThObservable {
public:

    /**
     * @brief unitarityNLO2 constructor.
     */
    unitarityNLO2(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO2@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO3
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO3: public ThObservable {
public:

    /**
     * @brief unitarityNLO3 constructor.
     */
    unitarityNLO3(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO3@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO4
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO4: public ThObservable {
public:

    /**
     * @brief unitarityNLO4 constructor.
     */
    unitarityNLO4(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO4@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO5
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO5: public ThObservable {
public:

    /**
     * @brief unitarityNLO5 constructor.
     */
    unitarityNLO5(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO5@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO6
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO6: public ThObservable {
public:

    /**
     * @brief unitarityNLO6 constructor.
     */
    unitarityNLO6(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO6@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO7
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO7: public ThObservable {
public:

    /**
     * @brief unitarityNLO7 constructor.
     */
    unitarityNLO7(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO7@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO8
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO8: public ThObservable {
public:

    /**
     * @brief unitarityNLO8 constructor.
     */
    unitarityNLO8(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO8@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO9
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO9: public ThObservable {
public:

    /**
     * @brief unitarityNLO9 constructor.
     */
    unitarityNLO9(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO9@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO10
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO10: public ThObservable {
public:

    /**
     * @brief unitarityNLO10 constructor.
     */
    unitarityNLO10(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO10@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO11
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO11: public ThObservable {
public:

    /**
     * @brief unitarityNLO11 constructor.
     */
    unitarityNLO11(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO11@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO12
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO12: public ThObservable {
public:

    /**
     * @brief unitarityNLO12 constructor.
     */
    unitarityNLO12(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO12@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO13
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO13: public ThObservable {
public:

    /**
     * @brief unitarityNLO13 constructor.
     */
    unitarityNLO13(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO13@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO14
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO14: public ThObservable {
public:

    /**
     * @brief unitarityNLO14 constructor.
     */
    unitarityNLO14(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO14@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO15
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO15: public ThObservable {
public:

    /**
     * @brief unitarityNLO15 constructor.
     */
    unitarityNLO15(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO15@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO16
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO16: public ThObservable {
public:

    /**
     * @brief unitarityNLO16 constructor.
     */
    unitarityNLO16(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO16@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO17
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO17: public ThObservable {
public:

    /**
     * @brief unitarityNLO17 constructor.
     */
    unitarityNLO17(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO17@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO18
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO18: public ThObservable {
public:

    /**
     * @brief unitarityNLO18 constructor.
     */
    unitarityNLO18(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO18@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO19
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO19: public ThObservable {
public:

    /**
     * @brief unitarityNLO19 constructor.
     */
    unitarityNLO19(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO19@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO20
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO20: public ThObservable {
public:

    /**
     * @brief unitarityNLO20 constructor.
     */
    unitarityNLO20(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO20@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO21
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO21: public ThObservable {
public:

    /**
     * @brief unitarityNLO21 constructor.
     */
    unitarityNLO21(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO21@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO22
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO22: public ThObservable {
public:

    /**
     * @brief unitarityNLO22 constructor.
     */
    unitarityNLO22(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO22@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO23
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO23: public ThObservable {
public:

    /**
     * @brief unitarityNLO23 constructor.
     */
    unitarityNLO23(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO23@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO24
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO24: public ThObservable {
public:

    /**
     * @brief unitarityNLO24 constructor.
     */
    unitarityNLO24(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO24@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO25
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO25: public ThObservable {
public:

    /**
     * @brief unitarityNLO25 constructor.
     */
    unitarityNLO25(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO25@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO26
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO26: public ThObservable {
public:

    /**
     * @brief unitarityNLO26 constructor.
     */
    unitarityNLO26(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO26@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev1
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev1: public ThObservable {
public:

    /**
     * @brief unitarityNLOev1 constructor.
     */
    unitarityNLOev1(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev1@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev2
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev2: public ThObservable {
public:

    /**
     * @brief unitarityNLOev2 constructor.
     */
    unitarityNLOev2(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev2@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev3
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev3: public ThObservable {
public:

    /**
     * @brief unitarityNLOev3 constructor.
     */
    unitarityNLOev3(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev3@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev4
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev4: public ThObservable {
public:

    /**
     * @brief unitarityNLOev4 constructor.
     */
    unitarityNLOev4(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev4@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev5
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev5: public ThObservable {
public:

    /**
     * @brief unitarityNLOev5 constructor.
     */
    unitarityNLOev5(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev5@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev6
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev6: public ThObservable {
public:

    /**
     * @brief unitarityNLOev6 constructor.
     */
    unitarityNLOev6(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev6@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev7
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev7: public ThObservable {
public:

    /**
     * @brief unitarityNLOev7 constructor.
     */
    unitarityNLOev7(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev7@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev8
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev8: public ThObservable {
public:

    /**
     * @brief unitarityNLOev8 constructor.
     */
    unitarityNLOev8(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev8@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev9
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev9: public ThObservable {
public:

    /**
     * @brief unitarityNLOev9 constructor.
     */
    unitarityNLOev9(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev9@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev10
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev10: public ThObservable {
public:

    /**
     * @brief unitarityNLOev10 constructor.
     */
    unitarityNLOev10(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev10@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev11
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev11: public ThObservable {
public:

    /**
     * @brief unitarityNLOev11 constructor.
     */
    unitarityNLOev11(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev11@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev12
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev12: public ThObservable {
public:

    /**
     * @brief unitarityNLOev12 constructor.
     */
    unitarityNLOev12(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev12@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev13
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev13: public ThObservable {
public:

    /**
     * @brief unitarityNLOev13 constructor.
     */
    unitarityNLOev13(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev13@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev14
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev14: public ThObservable {
public:

    /**
     * @brief unitarityNLOev14 constructor.
     */
    unitarityNLOev14(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev14@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev15
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev15: public ThObservable {
public:

    /**
     * @brief unitarityNLOev15 constructor.
     */
    unitarityNLOev15(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev15@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev16
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev16: public ThObservable {
public:

    /**
     * @brief unitarityNLOev16 constructor.
     */
    unitarityNLOev16(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev16@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev17
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev17: public ThObservable {
public:

    /**
     * @brief unitarityNLOev17 constructor.
     */
    unitarityNLOev17(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev17@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev18
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev18: public ThObservable {
public:

    /**
     * @brief unitarityNLOev18 constructor.
     */
    unitarityNLOev18(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev18@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

#endif	/* UNITARITY_H */

