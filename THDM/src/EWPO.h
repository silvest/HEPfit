/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWPO_H
#define	EWPO_H

#include <ThObservable.h>

class THDM;
class THDMcache;


/**
 * @class EWPO
 * @ingroup THDM 
 * @brief An observable class to calculate the electroweak precision observables in the %THDM.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details An observable class to calculate EWPO in context of the 2HDM.
 */
class EWPO : public ThObservable {
public:

    /**
     * @brief EWPO constructor.
     * @param[in] SM. A reference to a SM object
     */
    EWPO(const StandardModel& SM_i);


    double computeThValue();
    double dDelta_r();
    void computeTHDMcouplings();

    const THDM * myTHDM;
    THDMcache * mycache;

    private:
};

class  AlTHDM: public EWPO {
public:

    /**
     * @brief AlTHDM constructor.
     */
    AlTHDM(const StandardModel& SM_i);

    /**
     * @brief The left-right asymmetry in @f$e^+e^-\to Z\to e \bar{e}@f$ at the
     * @f$Z@f$-pole in the THDM
     * @return @f$\mathcal{A}_e@f$
     */
    double computeThValue ();
private:
};

class  PpoltauTHDM: public EWPO {
public:

    /**
     * @brief PpoltauTHDM constructor.
     */
    PpoltauTHDM(const StandardModel& SM_i);

    /**
     * @brief The left-right asymmetry in @f$e^+e^-\to Z\to \tau \bar{\tau}@f$ at the
     * @f$Z@f$-pole in the THDM
     * @return @f$\mathcal{P}_\tau@f$

     */
    double computeThValue ();
private:
};

class  AcTHDM: public EWPO {
public:

    /**
     * @brief AcTHDM constructor.
     */
    AcTHDM(const StandardModel& SM_i);

     /**
     * @brief The left-right asymmetry in @f$e^+e^-\to Z\to c \bar{c}@f$ at the
     * @f$Z@f$-pole in the THDM
     * @return @f$\mathcal{A}_c@f$

     */
    double computeThValue ();
private:
};

class  AbTHDM: public EWPO {
public:

    /**
     * @brief AbTHDM constructor.
     */
    AbTHDM(const StandardModel& SM_i);

 /**
     * @brief The left-right asymmetry in @f$e^+e^-\to Z\to b \bar{b}@f$ at the
     * @f$Z@f$-pole in the THDM
     * @return @f$\mathcal{A}_b@f$
     */
    double computeThValue ();
private:
};

class  AFBl0THDM: public EWPO {
public:

    /**
     * @brief AFBl0THDM constructor.
     */
    AFBl0THDM(const StandardModel& SM_i);

    /**
     * @brief The forward-backward assymetry for electrons at the
     * @f$Z@f$-pole in the THDM
     * @return @f$\mathcal{A}_{FB e}@f$
     */
    double computeThValue ();
private:
};

class  AFBc0THDM: public EWPO {
public:

    /**
     * @brief AFBc0THDM constructor.
     */
    AFBc0THDM(const StandardModel& SM_i);

    /**
     * @brief The forward-backward assymetry for charm quarks at the
     * @f$Z@f$-pole in the THDM
     * @return @f$\mathcal{A}_{FB c}@f$
     */
    double computeThValue ();
private:
};

class  AFBb0THDM: public EWPO {
public:

    /**
     * @brief AFBb0THDM constructor.
     */
    AFBb0THDM(const StandardModel& SM_i);

    /**
     * @brief The forward-backward assymetry for bottom quarks at the
     * @f$Z@f$-pole in the THDM
     * @return @f$\mathcal{A}_{FB b}@f$
     */
    double computeThValue ();
private:
};

class  GammaZTHDM: public EWPO {
public:

    /**
     * @brief GammaZTHDM constructor.
     */
    GammaZTHDM(const StandardModel& SM_i);

    /**
     * @brief  The total decay width of the @f$Z@f$ boson, @f$\Gamma_Z@f$
     * in the THDM
     * @return @f$\Gamma_Z@f$ in GeV
     */
    double computeThValue ();
private:
};

class  Rl0THDM: public EWPO {
public:

    /**
     * @brief Rl0THDM constructor.
     */
    Rl0THDM(const StandardModel& SM_i);

    /**
     * @brief Ratio between the decay width of  @f$Z@f$ to hadrons and 
     * to electrons @f$R_\ell^0= \frac{\Gamma_h}{\Gamma_\ell}@f$ in the THDM.
     * @return @f$R_\ell^0 @f$
     */
    double computeThValue ();
private:
};

class  Rc0THDM: public EWPO {
public:

    /**
     * @brief Rc0THDM constructor.
     */
    Rc0THDM(const StandardModel& SM_i);

    /**
     * @brief Ratio between the decay width of  @f$Z@f$ to hadrons and 
     * to charm quarks @f$R_\ell^0= \frac{\Gamma_h}{\Gamma_c}@f$ in the THDM.
     * @return @f$R_c^0 @f$
     */
    double computeThValue ();
private:
};

class  Rb0THDM: public EWPO {
public:

    /**
     * @brief Rb0THDM constructor.
     */
    Rb0THDM(const StandardModel& SM_i);

    /**
     * @brief Ratio between the decay width of  @f$Z@f$ to hadrons and 
     * to bottom quarks @f$R_\ell^0= \frac{\Gamma_h}{\Gamma_b}@f$\ in the THDM.
     * @return @f$R_b^0 @f$
     */
    double computeThValue ();
private:
};

class  SigmahadTHDM: public EWPO {
public:

    /**
     * @brief SigmahadTHDM constructor.
     */
    SigmahadTHDM(const StandardModel& SM_i);

    /*
     * @brief The hadronic cross section for @f$e^+e^- \to Z \to \mathrm{hadrons}@f$
     * at the @f$Z@f$-pole, @f$\sigma_h^0@f$ in the THDM.
     * @return @f$\sigma_h^0@f$ in GeV@f$^{-2}@f$
     */
    double computeThValue ();
private:
};

class  GammaWTHDM: public EWPO {
public:

    /**
     * @brief GammaWTHDM constructor.
     */
    GammaWTHDM(const StandardModel& SM_i);

    /**
     * @brief The total width of the @f$W@f$ boson, @f$\Gamma_W@f$.
     * @f$\Gamma_W@f$ in GeV in the THDM.
     * @return @f$\Gamma_W@f$ in GeV
     */
    double computeThValue ();
private:
};

class  sinthetaeffl_2THDM: public EWPO {
public:

    /**
     * @brief sinthetaeffl_2THDM constructor.
     */
    sinthetaeffl_2THDM(const StandardModel& SM_i);

    /**
     * @brief  The effective weak mixing angle @f$\sin^2\theta_{\rm eff}^{\,\ell}@f$
     * for @f$Z\ell\bar{\ell}@f$ at the the @f$Z@f$-mass scale.
     * @return @f$\sin^2\theta_{\rm eff}^{\,\ell}@f$
     */
    double computeThValue ();
private:
};

class  MWTHDM: public EWPO {
public:

    /**
     * @brief MWTHDM constructor.
     */
    MWTHDM(const StandardModel& SM_i);

    /**
     * brief  @f$W@f$-boson mass in the on-shell scheme,
     * @f$M_{W,\mathrm{SM}}@f$ in the THDM.
     * @return  @f$M_{W,\mathrm{THDM}}@f$ in GeV
     **/
    double computeThValue ();
private:
};

#endif	/* EWPO_H */
