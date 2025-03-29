/*
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HIGGSTHOBSERVABLES_H
#define	HIGGSTHOBSERVABLES_H

#include "ThObservable.h"

class NPbase;

/**
 * @class muggH
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{ggH}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{ggH}@f$ between the gluon-gluon
 * fusion Higgs production cross-section in the current model and in the Standard Model.
 */
class muggH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muggH(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{ggH}@f$ in the current model.
     * @return @f$\mu_{ggH}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVBF
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{VBF}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{VBF}@f$ between the vector-boson
 * fusion Higgs production cross-section in the current model and in the Standard Model.
 */
class muVBF : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVBF(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{VBF}@f$ in the current model.
     * @return @f$\mu_{VBF}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVBFgamma
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{VBF+\gamma}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{VBF+\gamma}@f$ between the vector-boson
 * fusion Higgs production cross-section in association with a hard photon
 * in the current model and in the Standard Model.
 */
class muVBFgamma : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVBFgamma(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{VBF+\gamma}@f$ in the current model.
     * @return @f$\mu_{VBF+\gamma}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class mueeWBF
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{eeWBF}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{eeWBF}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeWBF : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeWBF(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{eeWBF}@f$ in the current model.
     * @return @f$\mu_{eeWBF}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class mueeWBFPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{eeWBF}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{eeWBF}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H @f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeWBFPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeWBFPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$\mu_{eeWBF}@f$ in the current model.
     * @return @f$\mu_{eeWBF}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeHvv
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$, 
 * excluding contributions from on-shell @f$Z\to \nu\bar{\nu} @f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$ between the 
 * @f$e^+e^- \to H\nu\bar{\nu}@f$ 
 * production cross-section in the current model and in the Standard Model.
 */
class mueeHvv : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeHvv(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class mueeHvvPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$, 
 * excluding contributions from on-shell @f$Z\to \nu\bar{\nu} @f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$ between the 
 * @f$e^+e^- \to H\nu\bar{\nu}@f$ 
 * production cross-section in the current model and in the Standard Model.
 */
class mueeHvvPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeHvvPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};

/**
 * @class mueeZBF
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^{+}e^{-}\to e^{+}e^{-} H}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^{+}e^{-}\to e^{+}e^{-} H}@f$ between the 
 * @f$ e^{+}e^{-}\to e^{+}e^{-} H @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeZBF : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeZBF(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^{+}e^{-}\to e^{+}e^{-} H}@f$ in the current model.
     * @return @f$\mu_{e^{+}e^{-}\to e^{+}e^{-} H}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeZBFPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^{+}e^{-}\to e^{+}e^{-} H}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^{+}e^{-}\to e^{+}e^{-} H}@f$ between the 
 * @f$e^{+}e^{-}\to e^{+}e^{-} H@f$ 
 * production cross-section in the current model and in the Standard Model.
 */
class mueeZBFPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeZBFPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^{+}e^{-}\to e^{+}e^{-} H}@f$ in the current model.
     * @return @f$\mu_{e^{+}e^{-}\to e^{+}e^{-} H}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};

/**
 * @class muepWBF
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{epWBF}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{epWBF}@f$ between the 
 * @f$ e^{-}p\to \nu j H @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class muepWBF : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muepWBF(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{epWBF}@f$ in the current model.
     * @return @f$\mu_{epWBF}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muepZBF
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{epZBF}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{epZBF}@f$ between the 
 * @f$ e^{-}p\to e^{-} j H @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class muepZBF : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muepZBF(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{epZBF}@f$ in the current model.
     * @return @f$\mu_{epZBF}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muWH
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{WH}@f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{WH}@f$ between the W Higgs 
 * associated production cross-section in the current model and in the Standard Model. 
 */
class muWH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muWH(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{WH}@f$ in the current model.
     * @return @f$\mu_{WH}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muWHpT250
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{WH,p_T>250}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{WH,p_T>250}@f$ between the WH 
 * associated production cross-section in the current model and in the Standard Model, with @f$p_{T,H}>250@f$ GeV.
 */
class muWHpT250 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muWHpT250(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{WH,p_T>250}@f$ in the current model.
     * @return @f$\mu_{WH,p_T>250}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muZH
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{ZH}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{ZH}@f$ between the Z Higgs 
 * associated production cross-section in the current model and in the Standard Model.
 */
class muZH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muZH(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{ZH}@f$ in the current model.
     * @return @f$\mu_{ZH}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muZHpT250
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{ZH,p_T>250}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{ZH,p_T>250}@f$ between the ZH 
 * associated production cross-section in the current model and in the Standard Model, with @f$p_{T,H}>250@f$ GeV.
 */
class muZHpT250 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muZHpT250(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{ZH,p_T>250}@f$ in the current model.
     * @return @f$\mu_{ZH,p_T>250}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeZHGen
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH}@f$ between the 
 * @f$e^+e^- \to ZH@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZHGen : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeZHGen(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};

/**
 * @class mueeZH
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH}@f$ between the 
 * @f$e^+e^- \to ZH@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeZH(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};

/**
 * @class mueeZllH
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, Z \to e^+ e^-, \mu^+ \mu^- }@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$ between the 
 * @f$e^+e^- \to ZH, Z \to e^+ e^-, \mu^+ \mu^-@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZllH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeZllH(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class mueeZllHPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$ between the 
 * @f$e^+e^- \to ZH, Z \to e^+ e^-, \mu^+ \mu^-@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZllHPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeZllHPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);
    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};

/**
 * @class mueeZqqH
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, Z \to q \bar{q}}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, Z \to q \bar{q}}@f$ between the 
 * @f$e^+e^- \to ZH, Z \to q \bar{q}@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZqqH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeZqqH(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, Z \to q \bar{q}}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, Z \to q \bar{q}}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class mueeZqqHPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, Z \to q \bar{q}}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, Z \to q \bar{q}}@f$ between the 
 * @f$e^+e^- \to ZH, Z \to q \bar{q}@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZqqHPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeZqqHPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, Z \to q \bar{q}}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, Z \to q \bar{q}}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class aPsk
 * @ingroup NewPhysics
 * @brief A class for computing the angular parameter @f$a@f$ from 
 * @f$\mu_{e^+e^- \to ZH}@f$ (arXiv:1708.09079 [hep-ph]).
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the angular parameter @f$a@f$ from @f$\mu_{e^+e^- \to ZH}@f$.
 */
class aPsk : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    aPsk(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of the angular parameter @f$a@f$ from @f$\mu_{e^+e^- \to ZH}@f$ in the current model.
     * @return @f$a_{e^+e^- \to ZH}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class bPsk
 * @ingroup NewPhysics
 * @brief A class for computing the angular parameter @f$b@f$ from 
 * @f$\mu_{e^+e^- \to ZH}@f$ (arXiv:1708.09079 [hep-ph]).
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the angular parameter @f$b@f$ from @f$\mu_{e^+e^- \to ZH}@f$.
 */
class bPsk : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    bPsk(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of the angular parameter @f$b@f$ from @f$\mu_{e^+e^- \to ZH}@f$ in the current model.
     * @return @f$a_{e^+e^- \to ZH}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class muVH
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{VH}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{VH}@f$ between the WH+ZH 
 * associated production cross-section in the current model and in the Standard Model.
 */
class muVH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVH(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{VH}@f$ in the current model.
     * @return @f$\mu_{VH}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVHpT250
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{VH,p_T>250}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{VH,p_T>250}@f$ between the WH+ZH 
 * associated production cross-section in the current model and in the Standard Model, with @f$p_{T,H}>250@f$ GeV.
 */
class muVHpT250 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVHpT250(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{VH,p_T>250}@f$ in the current model.
     * @return @f$\mu_{VH,p_T>250}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVBFpVH
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{VBF+VH}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{VBF+VH}@f$ between the sum of
 * VBF and WH+ZH 
 * associated production cross-section in the current model and in the Standard Model.
 */
class muVBFpVH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVBFpVH(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{VBF+VH}@f$ in the current model.
     * @return @f$\mu_{VBF+VH}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muttH
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{ttH}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{ttH}@f$ between the t-tbar-Higgs 
 * associated production cross-section in the current model and in the Standard Model.
 */
class muttH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muttH(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{ttH}@f$ in the current model. 
     * @return @f$\mu_{ttH}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mutHq
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{tHq}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{tHq}@f$ between the t-q-Higgs 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mutHq : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mutHq(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{tHq}@f$ in the current model. 
     * @return @f$\mu_{tHq}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muggHpttH
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{ggH+ttH}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{ggH+ttH}@f$ between the sum
 * of gluon-gluon fusion and t-tbar-Higgs associated production cross-section in 
 * the current model and in the Standard Model.
 */
class muggHpttH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muggHpttH(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{ggH+ttH}@f$ in the current model. 
     * @return @f$\mu_{ggH+ttH}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class mueettH
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{eettH}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{eettH}@f$ between the 
 * @f$ e^{+}e^{-}\to t\bar{t} H @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueettH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueettH(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{eettH}@f$ in the current model.
     * @return @f$\mu_{eettH}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueettHPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{eettH}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{eettH}@f$ between the 
 * @f$ e^{+}e^{-}\to t\bar{t} H @f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueettHPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueettHPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$\mu_{eettH}@f$ in the current model.
     * @return @f$\mu_{eettH}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mummH
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu\mu H}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu\mu H}@f$ between the @f$\sigma(\mu \mu \to H)}@f$
 * production cross-section in the current model and in the Standard Model.
 */
class mummH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummH(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu\mu H}@f$ in the current model.
     * @return @f$\mu_{\mu\mu H}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHNWA
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu\mu H}@f$ in the narrow width approximation.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu\mu H}@f$ between the @f$\sigma(\mu \mu \to H)}@f$
 * production cross-section in the current model and in the Standard Model, in the narrow width approximation.
 */
class mummHNWA : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHNWA(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu\mu H}@f$ in the current model.
     * @return @f$\mu_{\mu\mu H}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


// -----------------------------------------------------------------------------
// Decay widths
// -----------------------------------------------------------------------------

/**
 * @class GammaHtoggRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the @f$\Gamma(H\to gg)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the @f$\Gamma(H\to gg)@f$
 * in the current model and in the Standard Model.
 */
class GammaHtoggRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    GammaHtoggRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the @f$\Gamma(H\to gg)@f$
     * in the current model and in the Standard Model.
     * @return @f$\Gamma(H\to gg)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};

/**
 * @class GammaHtoWWRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the @f$\Gamma(H\to WW\to 4f)@f$ with @f$f@f$ any fermion.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the @f$\Gamma(H\to WW\to 4f)@f$.
 * in the current model and in the Standard Model
 */
class GammaHtoWWRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    GammaHtoWWRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the @f$\Gamma(H\to WW\to 4f)@f$
     * in the current model and in the Standard Model.
     * @return @f$\Gamma(H\to WW\to 4f)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class GammaHtoZZRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the @f$\Gamma(H\to ZZ \to 4f)@f$ with @f$f@f$ any fermion.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the @f$\Gamma(H\to ZZ\to 4f)@f$
 * in the current model and in the Standard Model.
 */
class GammaHtoZZRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    GammaHtoZZRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the @f$\Gamma(H\to ZZ\to 4f)@f$
     * in the current model and in the Standard Model.
     * @return @f$\Gamma(H\to ZZ\to 4f)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};

/**
 * @class GammaHtoZgaRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the @f$\Gamma(H\to Z\gamma)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the @f$\Gamma(H\to Z\gamma)@f$
 * in the current model and in the Standard Model.
 */
class GammaHtoZgaRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    GammaHtoZgaRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the @f$\Gamma(H\to Z\gamma)@f$
     * in the current model and in the Standard Model.
     * @return @f$\Gamma(H\to Z\gamma)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};

/**
 * @class GammaHtogagaRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the @f$\Gamma(H\to\gamma\gamma)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the @f$\Gamma(H\to\gamma\gamma)@f$
 * in the current model and in the Standard Model.
 */
class GammaHtogagaRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    GammaHtogagaRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the @f$\Gamma(H\to \gamma\gamma)@f$
     * in the current model and in the Standard Model.
     * @return @f$\Gamma(H\to \gamma\gamma)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};

/**
 * @class GammaHtomumuRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the @f$\Gamma(H\to \mu^+\mu^-)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the @f$\Gamma(H\to \mu^+\mu^-)@f$
 * in the current model and in the Standard Model.
 */
class GammaHtomumuRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    GammaHtomumuRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the @f$\Gamma(H\to \mu^+\mu^-)@f$
     * in the current model and in the Standard Model.
     * @return @f$\Gamma(H\to \mu^+\mu^-)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};

/**
 * @class GammaHtotautauRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the @f$\Gamma(H\to \tau^+\tau^-)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the @f$\Gamma(H\to \tau^+\tau^-)@f$
 * in the current model and in the Standard Model.
 */
class GammaHtotautauRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    GammaHtotautauRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the @f$\Gamma(H\to \tau^+\tau^-)@f$
     * in the current model and in the Standard Model.
     * @return @f$\Gamma(H\to \tau^+\tau^-)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};

/**
 * @class GammaHtossRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the @f$\Gamma(H\to s\bar{s})@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the @f$\Gamma(H\to s\bar{s})@f$
 * in the current model and in the Standard Model.
 */
class GammaHtossRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    GammaHtossRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the @f$\Gamma(H\to s\bar{s})@f$
     * in the current model and in the Standard Model.
     * @return @f$\Gamma(H\to s\bar{s})@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};

/**
 * @class GammaHtoccRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the @f$\Gamma(H\to c\bar{c})@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the @f$\Gamma(H\to c\bar{c})@f$
 * in the current model and in the Standard Model.
 */
class GammaHtoccRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    GammaHtoccRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the @f$\Gamma(H\to c\bar{c})@f$
     * in the current model and in the Standard Model.
     * @return @f$\Gamma(H\to c\bar{c})@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};

/**
 * @class GammaHtobbRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the @f$\Gamma(H\to b\bar{b})@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the @f$\Gamma(H\to b\bar{b})@f$
 * in the current model and in the Standard Model.
 */
class GammaHtobbRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    GammaHtobbRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the @f$\Gamma(H\to b\bar{b})@f$
     * in the current model and in the Standard Model.
     * @return @f$\Gamma(H\to b\bar{b})@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class GammaHRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\Gamma_{H}/\Gamma_{H}^{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Higgs width
 * in the current model and in the Standard Model.
 */
class GammaHRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    GammaHRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Higgs width
     * in the current model and in the Standard Model.
     * @return @f$\Gamma_{H}/\Gamma_{H}^{SM}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};

// -----------------------------------------------------------------------------
// Branching Ratios
// -----------------------------------------------------------------------------

/**
 * @class BrHtoinvRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to invisible)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to invisible)@f$
 * in the current model and in the Standard Model.
 */
class BrHtoinvRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtoinvRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to invisible)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to invisible)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHinvisible
 * @ingroup NewPhysics
 * @brief A class for computing the branching ratio of Higgs decays into 
 * invisible particles.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the branching ratio Br@f$(H\to invisible)@f$.
 */
class BrHinvisible : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHinvisible(const StandardModel& SM_i);

    /**
     * @brief A method to compute the branching ratio of Higgs decays into
     * invisible particles.
     * @return Br@f$(H\to invisible)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHinvisibleNP
 * @ingroup NewPhysics
 * @brief A class for computing the branching ratio of Higgs decays into 
 * invisible particles (only decays into new invisible particles).
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the branching ratio Br@f$(H\to invisible, NP)@f$.
 */
class BrHinvisibleNP : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHinvisibleNP(const StandardModel& SM_i);

    /**
     * @brief A method to compute the branching ratio of Higgs decays into
     * invisible particles (only decays into new invisible particles).
     * @return Br@f$(H\to invisible, NP)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHexotic
 * @ingroup NewPhysics
 * @brief A class for computing the branching ratio of Higgs decays into 
 * exotics (invisible or not).
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the branching ratio Br@f$(H\to exotics)@f$.
 */
class BrHexotic : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHexotic(const StandardModel& SM_i);

    /**
     * @brief A method to compute the branching ratio of Higgs decays into
     * exotics (invisible or not).
     * @return Br@f$(H\to exotic)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtovisRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to visible)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to visible)@f$
 * in the current model and in the Standard Model.
 */
class BrHtovisRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtovisRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to visible)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to visible)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtoggRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to gg)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to gg)@f$
 * in the current model and in the Standard Model.
 */
class BrHtoggRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtoggRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to gg)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to gg)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtoWWRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to WW\to 4f)@f$ with @f$f@f$ any fermion.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to WW\to 4f)@f$.
 * in the current model and in the Standard Model
 */
class BrHtoWWRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtoWWRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to WW\to 4f)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to WW\to 4f)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtoZZRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to ZZ \to 4f)@f$ with @f$f@f$ any fermion.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to ZZ\to 4f)@f$
 * in the current model and in the Standard Model.
 */
class BrHtoZZRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    BrHtoZZRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to ZZ\to 4f)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to ZZ\to 4f)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtoVVRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to VV \to 4f)@f$ with @f$f@f$ any fermion.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to VV\to 4f)@f$
 * in the current model and in the Standard Model.
 */
class BrHtoVVRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    BrHtoVVRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to VV\to 4f)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to VV\to 4f)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtoZgaRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to Z\gamma)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to Z\gamma)@f$
 * in the current model and in the Standard Model.
 */
class BrHtoZgaRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtoZgaRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to Z\gamma)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to Z\gamma)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtoZgallRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to Z\gamma\to ll\gamma)@f$ with @f$l=e,\mu@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to Z\gamma\to ll\gamma)@f$
 * in the current model and in the Standard Model.
 */
class BrHtoZgallRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtoZgallRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to Z\gamma\to ll\gamma)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to Z\gamma\to ll\gamma)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtoZgaeeRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to Z\gamma\to ee\gamma)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to Z\gamma\to ee\gamma)@f$
 * in the current model and in the Standard Model.
 */
class BrHtoZgaeeRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtoZgaeeRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to Z\gamma\to ee\gamma)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to Z\gamma\to ee\gamma)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtoZgamumuRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to Z\gamma\to \mu\mu\gamma)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to Z\gamma\to \mu\mu\gamma)@f$
 * in the current model and in the Standard Model.
 */
class BrHtoZgamumuRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtoZgamumuRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to Z\gamma\to \mu\mu\gamma)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to Z\gamma\to \mu\mu\gamma)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtogagaRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to\gamma\gamma)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to\gamma\gamma)@f$
 * in the current model and in the Standard Model.
 */
class BrHtogagaRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtogagaRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to \gamma\gamma)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \gamma\gamma)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtomumuRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to \mu^+\mu^-)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to \mu^+\mu^-)@f$
 * in the current model and in the Standard Model.
 */
class BrHtomumuRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtomumuRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to \mu^+\mu^-)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \mu^+\mu^-)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtotautauRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to \tau^+\tau^-)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to \tau^+\tau^-)@f$
 * in the current model and in the Standard Model.
 */
class BrHtotautauRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtotautauRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to \tau^+\tau^-)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \tau^+\tau^-)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtoccRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to c\bar{c})@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to c\bar{c})@f$
 * in the current model and in the Standard Model.
 */
class BrHtoccRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtoccRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to c\bar{c})@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to c\bar{c})@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtobbRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to b\bar{b})@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to b\bar{b})@f$
 * in the current model and in the Standard Model.
 */
class BrHtobbRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtobbRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to b\bar{b})@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to b\bar{b})@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


// -----------------------------------------------------------------------------
// More 4 fermion decays
// -----------------------------------------------------------------------------

/**
 * @class BrHto2l2vRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to 2l 2\nu)@f$ with @f$l=e,\mu@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to 2l 2\nu)@f$
 * in the current model and in the Standard Model.
 */
class BrHto2l2vRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    BrHto2l2vRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to 2l 2\nu)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to 2l 2\nu)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtoevmuvRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to e\nu \mu\nu)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to e\nu \mu\nu)@f$.
 * in the current model and in the Standard Model
 */
class BrHtoevmuvRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtoevmuvRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to e\nu \mu\nu)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to e\nu \mu\nu)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHto2e2vRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to 2e 2\nu)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to 2e 2\nu)@f$
 * in the current model and in the Standard Model.
 */
class BrHto2e2vRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    BrHto2e2vRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to 2e 2\nu)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to 2e 2\nu)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHto2mu2vRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to 2\mu 2\nu)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to 2\mu 2\nu)@f$
 * in the current model and in the Standard Model.
 */
class BrHto2mu2vRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    BrHto2mu2vRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to 2\mu 2\nu)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to 2\mu 2\nu)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHto4lRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to 4l)@f$ with @f$l=e,\mu@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to 4l)@f$
 * in the current model and in the Standard Model.
 */
class BrHto4lRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    BrHto4lRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to 4l)@f$ with @f$l=e,\mu@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to 4l)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHto4eRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to 4e)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to 4e)@f$
 * in the current model and in the Standard Model.
 */
class BrHto4eRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    BrHto4eRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to 4e)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to 4e)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHto4muRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to 4\mu)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to 4\mu)@f$
 * in the current model and in the Standard Model.
 */
class BrHto4muRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    BrHto4muRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to 4\mu)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to 4\mu)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHto2e2muRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to 2 e 2\mu)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to 2 e 2\mu)@f$
 * in the current model and in the Standard Model.
 */
class BrHto2e2muRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    BrHto2e2muRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to 2 e 2\mu)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to 2 e 2\mu)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtolljjRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to l l j j)@f$, @f$l=e,\mu,~~j\not=b)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to l l j j)@f$, @f$l=e,\mu,~~j\not=b)@f$
 * in the current model and in the Standard Model.
 */
class BrHtolljjRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    BrHtolljjRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to l l j j)@f$, @f$l=e,\mu,~~j\not=b)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to l l j j)@f$, @f$l=e,\mu,~~j\not=b)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtolvjjRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to l \nu j j)@f$, @f$l=e,\mu,~~j\not=b)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to l \nu j j)@f$, @f$l=e,\mu,~~j\not=b)@f$
 * in the current model and in the Standard Model.
 */
class BrHtolvjjRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    BrHtolvjjRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to l \nu j j)@f$, @f$l=e,\mu,~~j\not=b)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to l \nu j j)@f$, @f$l=e,\mu,~~j\not=b)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtolv_lvorjjRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to l \nu l \nu, l \nu j j)@f$, @f$l=e,\mu,~~j\not=b)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to l \nu l \nu, l \nu j j)@f$, @f$l=e,\mu,~~j\not=b)@f$
 * in the current model and in the Standard Model.
 */
class BrHtolv_lvorjjRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    BrHtolv_lvorjjRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to l \nu l \nu, l \nu j j)@f$, @f$l=e,\mu,~~j\not=b)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to l \nu l \nu, l \nu j j)@f$, @f$l=e,\mu,~~j\not=b)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtoll_vvorjjRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to l l \nu\nu, l l j j)@f$, @f$l=e,\mu,~~j\not=b)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to l l \nu\nu, l l j j)@f$, @f$l=e,\mu,~~j\not=b)@f$
 * in the current model and in the Standard Model.
 */
class BrHtoll_vvorjjRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    BrHtoll_vvorjjRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to l l \nu\nu, l l j j)@f$, @f$l=e,\mu,~~j\not=b)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to l l \nu\nu, l l j j)@f$, @f$l=e,\mu,~~j\not=b)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


// -----------------------------------------------------------------------------
// Ratios of BR (ratios with SM)
// -----------------------------------------------------------------------------

/**
 * @class BrHtogaga_over_mumu_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to \gamma\gamma)/@f$Br@f$(H\to \mu\mu)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to \gamma\gamma)/@f$Br@f$(H\to \mu\mu)@f$
 * in the current model and in the Standard Model.
 */
class BrHtogaga_over_mumu_Ratio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtogaga_over_mumu_Ratio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to \gamma\gamma)/@f$Br@f$(H\to \mu\mu)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \gamma\gamma)/@f$Br@f$(H\to \mu\mu)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtoZga_over_mumu_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to Z\gamma)/@f$Br@f$(H\to \mu\mu)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to Z\gamma)/@f$Br@f$(H\to \mu\mu)@f$
 * in the current model and in the Standard Model.
 */
class BrHtoZga_over_mumu_Ratio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtoZga_over_mumu_Ratio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to Z\gamma)/@f$Br@f$(H\to \mu\mu)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to Z\gamma)/@f$Br@f$(H\to \mu\mu)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtoZmumuga_over_mumu_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to Z\gamma\to\mu\mu \gamma)/@f$Br@f$(H\to \mu\mu)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to Z\gamma\to\mu\mu \gamma)/@f$Br@f$(H\to \mu\mu)@f$
 * in the current model and in the Standard Model.
 */
class BrHtoZmumuga_over_mumu_Ratio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtoZmumuga_over_mumu_Ratio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to Z\gamma\to\mu\mu \gamma)/@f$Br@f$(H\to \mu\mu)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to Z\gamma\to\mu\mu \gamma)/@f$Br@f$(H\to \mu\mu)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtoZga_over_4mu_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to Z\gamma)/@f$Br@f$(H\to 4\mu)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to Z\gamma)/@f$Br@f$(H\to 4\mu)@f$
 * in the current model and in the Standard Model.
 */
class BrHtoZga_over_4mu_Ratio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtoZga_over_4mu_Ratio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to Z\gamma)/@f$Br@f$(H\to 4\mu)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to Z\gamma)/@f$Br@f$(H\to 4\mu)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtoZmumuga_over_4mu_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to Z\gamma\to\mu\mu \gamma)/@f$Br@f$(H\to 4\mu)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to Z\gamma\to\mu\mu \gamma)/@f$Br@f$(H\to 4\mu)@f$
 * in the current model and in the Standard Model.
 */
class BrHtoZmumuga_over_4mu_Ratio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtoZmumuga_over_4mu_Ratio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to Z\gamma\to\mu\mu \gamma)/@f$Br@f$(H\to 4\mu)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to Z\gamma\to\mu\mu \gamma)/@f$Br@f$(H\to 4\mu)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtogaga_over_4l_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to \gamma\gamma)/@f$Br@f$(H\to 4\ell)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to \gamma\gamma)/@f$Br@f$(H\to 4\ell)@f$
 * in the current model and in the Standard Model.
 */
class BrHtogaga_over_4l_Ratio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtogaga_over_4l_Ratio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to \gamma\gamma)/@f$Br@f$(H\to 4\ell)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \gamma\gamma)/@f$Br@f$(H\to 4\ell)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};



/**
 * @class BrHtobb_over_4l_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to bb)/@f$Br@f$(H\to 4\ell)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to bb)/@f$Br@f$(H\to 4\ell)@f$
 * in the current model and in the Standard Model.
 */
class BrHtobb_over_4l_Ratio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtobb_over_4l_Ratio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to bb)/@f$Br@f$(H\to 4\ell)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to bb)/@f$Br@f$(H\to 4\ell)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};



/**
 * @class BrHto2l2v_over_4l_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to 2\ell 2\nu)/@f$Br@f$(H\to 4\ell)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to 2\ell 2\nu)/@f$Br@f$(H\to 4\ell)@f$
 * in the current model and in the Standard Model.
 */
class BrHto2l2v_over_4l_Ratio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHto2l2v_over_4l_Ratio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to 2\ell 2\nu)/@f$Br@f$(H\to 4\ell)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to 2\ell 2\nu)/@f$Br@f$(H\to 4\ell)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};



/**
 * @class BrHtotautau_over_4l_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to \tau\tau)/@f$Br@f$(H\to 4\ell)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to \tau\tau)/@f$Br@f$(H\to 4\ell)@f$
 * in the current model and in the Standard Model.
 */
class BrHtotautau_over_4l_Ratio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtotautau_over_4l_Ratio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to \tau\tau)/@f$Br@f$(H\to 4\ell)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \tau\tau)/@f$Br@f$(H\to 4\ell)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};



/**
 * @class BrHtogaga_over_2e2mu_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to \gamma\gamma)/@f$Br@f$(H\to 2e 2\mu)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to \gamma\gamma)/@f$Br@f$(H\to 2e 2\mu)@f$
 * in the current model and in the Standard Model.
 */
class BrHtogaga_over_2e2mu_Ratio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtogaga_over_2e2mu_Ratio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to \gamma\gamma)/@f$Br@f$(H\to 2e 2\mu)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \gamma\gamma)/@f$Br@f$(H\to 2e 2\mu)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtoZga_over_4l_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to Z\gamma)/@f$Br@f$(H\to 4\ell)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to Z\gamma)/@f$Br@f$(H\to 4\ell)@f$
 * in the current model and in the Standard Model.
 */
class BrHtoZga_over_4l_Ratio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtoZga_over_4l_Ratio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to Z\gamma)/@f$Br@f$(H\to 4\ell)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to Z\gamma)/@f$Br@f$(H\to 4\ell)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};



/**
 * @class BrHtomumu_over_4l_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to \mu\mu)/@f$Br@f$(H\to 4\ell)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to \mu\mu)/@f$Br@f$(H\to 4\ell)@f$
 * in the current model and in the Standard Model.
 */
class BrHtomumu_over_4l_Ratio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtomumu_over_4l_Ratio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to \mu\mu)/@f$Br@f$(H\to 4\ell)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \mu\mu)/@f$Br@f$(H\to 4\ell)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtomumu_over_4mu_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to \mu\mu)/@f$Br@f$(H\to 4\mu)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to \mu\mu)/@f$Br@f$(H\to 4\mu)@f$
 * in the current model and in the Standard Model.
 */
class BrHtomumu_over_4mu_Ratio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtomumu_over_4mu_Ratio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to \mu\mu)/@f$Br@f$(H\to 4\mu)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \mu\mu)/@f$Br@f$(H\to 4\mu)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};



/**
 * @class BrHto4l_over_gaga_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to 4\ell)/@f$Br@f$(H\to \gamma\gamma)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to 4\ell)/@f$Br@f$(H\to \gamma\gamma)@f$
 * in the current model and in the Standard Model.
 */
class BrHto4l_over_gaga_Ratio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHto4l_over_gaga_Ratio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to 4\ell)/@f$Br@f$(H\to \gamma\gamma)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to 4\ell)/@f$Br@f$(H\to \gamma\gamma)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtoZga_over_gaga_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to Z \gamma)/@f$Br@f$(H\to \gamma\gamma)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to Z \gamma)/@f$Br@f$(H\to \gamma\gamma)@f$
 * in the current model and in the Standard Model.
 */
class BrHtoZga_over_gaga_Ratio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtoZga_over_gaga_Ratio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to Z \gamma)/@f$Br@f$(H\to \gamma\gamma)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to Z \gamma)/@f$Br@f$(H\to \gamma\gamma)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtomumu_over_gaga_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to \mu\mu)/@f$Br@f$(H\to \gamma\gamma)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to \mu\mu)/@f$Br@f$(H\to \gamma\gamma)@f$
 * in the current model and in the Standard Model.
 */
class BrHtomumu_over_gaga_Ratio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtomumu_over_gaga_Ratio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to \mu\mu)/@f$Br@f$(H\to \gamma\gamma)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \mu\mu)/@f$Br@f$(H\to \gamma\gamma)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHto2l2v_over_gaga_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio Br@f$(H\to 2\ell 2\nu)/@f$Br@f$(H\to\gamma\gamma)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio Br@f$(H\to 2\ell 2\nu)/@f$Br@f$(H\to\gamma\gamma)@f$
 * in the current model and in the Standard Model.
 */
class BrHto2l2v_over_gaga_Ratio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHto2l2v_over_gaga_Ratio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio Br@f$(H\to 2\ell 2\nu)/@f$Br@f$(H\to\gamma\gamma)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to 2\ell 2\nu)/@f$Br@f$(H\to\gamma\gamma)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtobb_over_cc_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to bb)/@f$Br@f$(H\to cc)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to bb)/@f$Br@f$(H\to cc)@f$
 * in the current model and in the Standard Model.
 */
class BrHtobb_over_cc_Ratio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtobb_over_cc_Ratio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to bb)/@f$Br@f$(H\to cc)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to bb)/@f$Br@f$(H\to cc)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};



/**
 * @class BrHtogaga_over_gg_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to \gamma\gamma)/@f$Br@f$(H\to gg)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to \gamma\gamma)/@f$Br@f$(H\to gg)@f$
 * in the current model and in the Standard Model.
 */
class BrHtogaga_over_gg_Ratio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtogaga_over_gg_Ratio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to \gamma\gamma)/@f$Br@f$(H\to gg)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \gamma\gamma)/@f$Br@f$(H\to gg)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};



/**
 * @class BrHtogg_over_bb_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to gg)/@f$Br@f$(H\to bb)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to gg)/@f$Br@f$(H\to bb)@f$
 * in the current model and in the Standard Model.
 */
class BrHtogg_over_bb_Ratio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtogg_over_bb_Ratio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to gg)/@f$Br@f$(H\to bb)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to gg)/@f$Br@f$(H\to bb)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};



/**
 * @class BrHtogg_over_cc_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to gg)/@f$Br@f$(H\to cc)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to gg)/@f$Br@f$(H\to cc)@f$
 * in the current model and in the Standard Model.
 */
class BrHtogg_over_cc_Ratio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtogg_over_cc_Ratio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to gg)/@f$Br@f$(H\to cc)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to gg)/@f$Br@f$(H\to cc)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class muggHgaga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muggHgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muggHgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muggHgagaInt
 * @ingroup NewPhysics
 * @brief To be used ONLY in Higgs Observables and for the diphoton decay. Replaces the 
 * narrow width approximation for the result including finite width effects of interference
 * with the background. Do NOT use for models with linearized expressions (e.g. NPSMEFTd6)
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muggHgagaInt : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muggHgagaInt(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muVBFHgaga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muVBFHgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVBFHgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muZHgaga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muZHgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muZHgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muWHgaga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muWHgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muWHgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVHgaga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muVHgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVHgaga(const StandardModel& SM_i, const double sqrt_s_i);
    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muttHgaga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muttHgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muttHgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class mutHgaga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class mutHgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mutHgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};                           //AG:added

/**
 * @class muggHZga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muggHZga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muggHZga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muggHpbbH_Hgaga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muggHpbbH_Hgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muggHpbbH_Hgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};                    //AG:added


/**
 * @class muggHZgamumu
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muggHZgamumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muggHZgamumu(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muVBFHZga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muVBFHZga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVBFHZga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muZHZga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muZHZga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muZHZga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();
private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muWHZga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muWHZga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muWHZga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVHZga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muVHZga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVHZga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muttHZga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muttHZga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muttHZga(const StandardModel& SM_i, const double sqrt_s_i);
    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muggHZZ
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muggHZZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muggHZZ(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVBFHZZ
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muVBFHZZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVBFHZZ(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muZHZZ
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muZHZZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muZHZZ(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muWHZZ
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muWHZZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muWHZZ(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVHZZ
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muVHZZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVHZZ(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muttHZZ
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muttHZZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muttHZZ(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muttHptH_HZZ
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muttHptH_HZZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muttHptH_HZZ(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};                       //AG:added

/**
 * @class muttHptH_Hgaga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muttHptH_Hgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muttHptH_Hgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};                     //AG:added

/**
 * @class muttHptH_Hmumu
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muttHptH_Hmumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muttHptH_Hmumu(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};                     //AG:added

/**
 * @class muggHpbbH_HZZ
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muggHpbbH_HZZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muggHpbbH_HZZ(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};                      //AG:added



/**
 * @class muggHZZ4l
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muggHZZ4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muggHZZ4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};



/**
 * @class muggHZZ4mu
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muggHZZ4mu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muggHZZ4mu(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muVBFHZZ4l
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muVBFHZZ4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVBFHZZ4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muZHZZ4l
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muZHZZ4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muZHZZ4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muWHZZ4l
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muWHZZ4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muWHZZ4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVHZZ4l
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muVHZZ4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVHZZ4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muttHZZ4l
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muttHZZ4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muttHZZ4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muggHWW
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muggHWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muggHWW(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVBFHWW
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muVBFHWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVBFHWW(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muZHWW
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muZHWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muZHWW(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muWHWW
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muWHWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muWHWW(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVHWW
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muVHWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVHWW(const StandardModel& SM_i, const double sqrt_s_i);
    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muttHWW
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muttHWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muttHWW(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muttHptH_HWW
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muttHptH_HWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muttHptH_HWW(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};                       //AG:added

/**
 * @class muggHpbbH_HWW
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muggHpbbH_HWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muggHpbbH_HWW(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};                      //AG:added


/**
 * @class muggHWW2l2v
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muggHWW2l2v : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muggHWW2l2v(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVBFHWW2l2v
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muVBFHWW2l2v : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVBFHWW2l2v(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muZHWW2l2v
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muZHWW2l2v : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muZHWW2l2v(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muWHWW2l2v
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muWHWW2l2v : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muWHWW2l2v(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVHWW2l2v
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muVHWW2l2v : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVHWW2l2v(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muttHWW2l2v
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muttHWW2l2v : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muttHWW2l2v(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muttHVV
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muttHVV : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muttHVV(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muggHmumu
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muggHmumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muggHmumu(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVBFHmumu
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muVBFHmumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVBFHmumu(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muZHmumu
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muZHmumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muZHmumu(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muWHmumu
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muWHmumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muWHmumu(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVHmumu
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muVHmumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVHmumu(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muttHmumu
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muttHmumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muttHmumu(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muggHpttHptHpbbH_Hmumu
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muggHpttHptHpbbH_Hmumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muggHpttHptHpbbH_Hmumu(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};             //AG:added

/**
 * @class muVBFpVH_Hmumu
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muVBFpVH_Hmumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVBFpVH_Hmumu(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};                     //AG:added


/**
 * @class muggHtautau
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muggHtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muggHtautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVBFHtautau
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muVBFHtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVBFHtautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muVBFpVHtautau
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muVBFpVHtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVBFpVHtautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muZHtautau
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muZHtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muZHtautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muWHtautau
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muWHtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muWHtautau(const StandardModel& SM_i, const double sqrt_s_i);
    
    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVHtautau
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muVHtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVHtautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muttHtautau
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muttHtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muttHtautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muttHptH_Htautau
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muttHptH_Htautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muttHptH_Htautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};                   //AG:added

/**
 * @class muggHpbbH_Htautau
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muggHpbbH_Htautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muggHpbbH_Htautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};                  //AG:added


/**
 * @class muggHbb
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muggHbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muggHbb(const StandardModel& SM_i, const double sqrt_s_i);
    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVBFHbb
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muVBFHbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVBFHbb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muZHbb
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muZHbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muZHbb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muWHbb
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muWHbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muWHbb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVHbb
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muVHbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVHbb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muttHbb
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muttHbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muttHbb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muttHptH_Hbb
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muttHptH_Hbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muttHptH_Hbb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};                        //AG:added

/**
 * @class muggHpVBFpbbH_Hbb
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muggHpVBFpbbH_Hbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muggHpVBFpbbH_Hbb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};                  //AG:added



/**
 * @class muVHcc
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muVHcc : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVHcc(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};



/**
 * @class muVBFBRinv
 * @ingroup NewPhysics
 * @brief A class for computing the quantity @f$\mu_{VBF}\times BR(H \to invisible}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the quantity @f$\mu_{VBF}\times BR(H \to invisible}@f$, i.e.
 * the ratio between the @f$pp \to jjH@f$ 
 * production cross-section in the current model and 
 * the SM, times the total invisible branching ratio.
 */
class muVBFBRinv : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVBFBRinv(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{VBF, H \to invisible}@f$ in the current model.
     * @return @f$\mu_{VBF, H \to invisible}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muVBFHinv
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{VBF, H \to invisible}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{VBF, H \to invisible}@f$ between the 
 * @f$pp \to jjH, H \to invisible@f$ 
 * production cross-section in the current model and the SM.
 */
class muVBFHinv : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVBFHinv(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{VBF, H \to invisible}@f$ in the current model.
     * @return @f$\mu_{VBF, H \to invisible}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muVHBRinv
 * @ingroup NewPhysics
 * @brief A class for computing the quantity @f$\mu_{pp \to VH, H \to invisible}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the quantity @f$\mu_{pp \to VH, H \to invisible}@f$, i.e.
 * the ratio between the @f$pp \to VH@f$ 
 * associated production cross-section in the current model and 
 * the SM, times the total invisible branching ratio.
 */
class muVHBRinv : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVHBRinv(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{pp \to VH}\times BR(H \to invisible}@f$ in the current model.
     * @return @f$\mu_{pp \to VH}\times BR(H \to invisible}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muVHinv
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{pp \to VH, H \to invisible}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{pp \to VH, H \to invisible}@f$ between the 
 * @f$pp \to VH, H \to invisible@f$ 
 * associated production cross-section in the current model and the SM.
 */
class muVHinv : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVHinv(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{pp \to VH, H \to invisible}@f$ in the current model.
     * @return @f$\mu_{pp \to VH, H \to invisible}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muppHmumu
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muppHmumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muppHmumu(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muppHZga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muppHZga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muppHZga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muggHH2ga2b
 * @ingroup NewPhysics
 * @brief A class for computing the signal strength for di-Higgs production via
 * gluon fusion in the @f$\gamma\gamma b b@f$ channel.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the signal strength for di-Higgs production via
 * gluon fusion in the @f$\gamma\gamma b b@f$ channel in the current model.
 */
class muggHH2ga2b : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muggHH2ga2b(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{ggH}@f$ in the current model.
     * @return @f$\mu_{ggH}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muttHZbbboost
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\sigma(ttH)/\sigma(ttZ)@f$ 
 * in the @f$H,Z\to b\bar{b}@f$ channel in the boosted region. Normalized to the SM.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\sigma(ttH)/\sigma(ttZ)@f$ 
 * in the @f$H,Z\to b\bar{b}@f$ channel in the boosted region in the current model.
 */
class muttHZbbboost : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muttHZbbboost(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\sigma(ttH)/\sigma(ttZ)@f$ 
     * in the @f$H,Z\to b\bar{b}@f$ channel in the current model.
     * @return @f$\sigma(ttH)/\sigma(ttZ)@f$ normalized to the SM.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muttHgagaZeeboost
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\sigma(ttH)/\sigma(ttZ)@f$ 
 * in the @f$H\to b\bar{b}@f$, @f$Z\to e^+e^-@f$ channel in the boosted region. Normalized to the SM.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\sigma(ttH)/\sigma(ttZ)@f$ 
 * in the @f$H\to b\bar{b}@f$, @f$Z\to e^+e^-@f$ channel in the boosted region in the current model.
 */
class muttHgagaZeeboost : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muttHgagaZeeboost(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\sigma(ttH)/\sigma(ttZ)@f$ 
     * in the @f$H,Z\to b\bar{b}@f$ channel in the current model.
     * @return @f$\sigma(ttH)/\sigma(ttZ)@f$ normalized to the SM.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


//AG:begin
/**
 * @class ggHgaga
 * @ingroup NewPhysics
 * @brief A class for computing the value of @f$\sigma(ggH) Br(\gamma\gamma)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class ggHgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    ggHgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\sigma(ggH) Br(\gamma\gamma)@f$.
     * @return @f$\sigma(ggH) Br(\gamma\gamma)@f$.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class ggHZZ
 * @ingroup NewPhysics
 * @brief A class for computing the value of @f$\sigma(ggH) Br(ZZ)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class ggHZZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    ggHZZ(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\sigma(ggH) Br(ZZ)@f$.
     * @return @f$\sigma(ggH) Br(ZZ)@f$.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class ggHWW
 * @ingroup NewPhysics
 * @brief A class for computing the value of @f$\sigma(ggH) Br(WW)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class ggHWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    ggHWW(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\sigma(ggH) Br(WW)@f$.
     * @return @f$\sigma(ggH) Br(WW)@f$.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class ggHtautau
 * @ingroup NewPhysics
 * @brief A class for computing the value of @f$\sigma(ggH) Br(\tau\tau)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class ggHtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    ggHtautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\sigma(ggH) Br(\tau\tau)@f$.
     * @return @f$\sigma(ggH) Br(\tau\tau)@f$.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

//--------------------------------------------------

/**
 * @class VBFHgaga
 * @ingroup NewPhysics
 * @brief A class for computing the value of @f$\sigma(VBF) Br(\gamma\gamma)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class VBFHgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    VBFHgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\sigma(VBF) Br(\gamma\gamma)@f$.
     * @return @f$\sigma(VBF) Br(\gamma\gamma)@f$.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class VBFHZZ
 * @ingroup NewPhysics
 * @brief A class for computing the value of @f$\sigma(VBF) Br(ZZ)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class VBFHZZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    VBFHZZ(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\sigma(VBF) Br(ZZ)@f$.
     * @return @f$\sigma(VBF) Br(ZZ)@f$.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class VBFHWW
 * @ingroup NewPhysics
 * @brief A class for computing the value of @f$\sigma(VBF) Br(WW)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class VBFHWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    VBFHWW(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\sigma(VBF) Br(WW)@f$.
     * @return @f$\sigma(VBF) Br(WW)@f$.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class VBFHtautau
 * @ingroup NewPhysics
 * @brief A class for computing the value of @f$\sigma(VBF) Br(\tau\tau)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class VBFHtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    VBFHtautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\sigma(VBF) Br(\tau\tau)@f$.
     * @return @f$\sigma(VBF) Br(\tau\tau)@f$.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

//--------------------------------------------------

/**
 * @class WHgaga
 * @ingroup NewPhysics
 * @brief A class for computing the value of @f$\sigma(WH) Br(\gamma\gamma)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class WHgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    WHgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\sigma(WH) Br(\gamma\gamma)@f$.
     * @return @f$\sigma(WH) Br(\gamma\gamma)@f$.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class WHWW
 * @ingroup NewPhysics
 * @brief A class for computing the value of @f$\sigma(WH) Br(WW)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class WHWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    WHWW(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\sigma(WH) Br(WW)@f$.
     * @return @f$\sigma(WH) Br(WW)@f$.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class ggHtautau
 * @ingroup NewPhysics
 * @brief A class for computing the value of @f$\sigma(ggH) Br(\tau\tau)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class WHtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    WHtautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\sigma(WH) Br(\tau\tau)@f$.
     * @return @f$\sigma(WH) Br(\tau\tau)@f$.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class WHbb
 * @ingroup NewPhysics
 * @brief A class for computing the value of @f$\sigma(WH) Br(bb)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class WHbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    WHbb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\sigma(WH) Br(bb)@f$.
     * @return @f$\sigma(WH) Br(bb)@f$.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

//--------------------------------------------------

/**
 * @class ZHgaga
 * @ingroup NewPhysics
 * @brief A class for computing the value of @f$\sigma(ZH) Br(\gamma\gamma)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class ZHgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    ZHgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\sigma(ZH) Br(\gamma\gamma)@f$.
     * @return @f$\sigma(ZH) Br(\gamma\gamma)@f$.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class ZHWW
 * @ingroup NewPhysics
 * @brief A class for computing the value of @f$\sigma(ZH) Br(WW)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class ZHWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    ZHWW(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\sigma(ZH) Br(WW)@f$.
     * @return @f$\sigma(ZH) Br(WW)@f$.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class ZHtautau
 * @ingroup NewPhysics
 * @brief A class for computing the value of @f$\sigma(ZH) Br(\tau\tau)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class ZHtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    ZHtautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\sigma(ZH) Br(\tau\tau)@f$.
     * @return @f$\sigma(ZH) Br(\tau\tau)@f$.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class ZHbb
 * @ingroup NewPhysics
 * @brief A class for computing the value of @f$\sigma(ZH) Br(bb)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class ZHbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    ZHbb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\sigma(ZH) Br(bb)@f$.
     * @return @f$\sigma(ZH) Br(bb)@f$.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

//--------------------------------------------------

/**
 * @class ttHgaga
 * @ingroup NewPhysics
 * @brief A class for computing the value of @f$\sigma(ttH) Br(\gamma\gamma)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class ttHgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    ttHgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\sigma(ttH) Br(\gamma\gamma)@f$.
     * @return @f$\sigma(ttH) Br(\gamma\gamma)@f$.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class ttHWW
 * @ingroup NewPhysics
 * @brief A class for computing the value of @f$\sigma(ttH) Br(WW)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class ttHWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    ttHWW(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\sigma(ttH) Br(WW)@f$.
     * @return @f$\sigma(ttH) Br(WW)@f$.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class ttHtautau
 * @ingroup NewPhysics
 * @brief A class for computing the value of @f$\sigma(ttH) Br(\tau\tau)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class ttHtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    ttHtautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\sigma(ttH) Br(\tau\tau)@f$.
     * @return @f$\sigma(ttH) Br(\tau\tau)@f$.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class ttHbb
 * @ingroup NewPhysics
 * @brief A class for computing the value of @f$\sigma(ttH) Br(bb)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class ttHbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    ttHbb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\sigma(ttH) Br(bb)@f$.
     * @return @f$\sigma(ttH) Br(bb)@f$.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};
//AG:end

/**
 * @class UpperLimit_ppHZgammaA
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class UpperLimit_ppHZgammaA : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    UpperLimit_ppHZgammaA(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute 
     * @return 
     */
    double computeThValue();
    
private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class UpperLimit_ppHZgammaA13
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class UpperLimit_ppHZgammaA13 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    UpperLimit_ppHZgammaA13(const StandardModel& SM_i, const double sqrt_s_i);
    
    /**
     * @brief A method to compute 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class UpperLimit_ppHZgammaC13
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class UpperLimit_ppHZgammaC13 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    UpperLimit_ppHZgammaC13(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class UpperLimit_ppHZgammaC
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class UpperLimit_ppHZgammaC : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    UpperLimit_ppHZgammaC(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class cg_plus_ct
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class cg_plus_ct : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    cg_plus_ct(const StandardModel& SM_i);

    /**
     * @brief A method to compute 
     * @return 
     */
    double computeThValue();
private:
    const NPbase* myNPbase;
};

/**
 * @class cga_plus_ct
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class cga_plus_ct : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    cga_plus_ct(const StandardModel& SM_i);

    /**
     * @brief A method to compute 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};

/**
 * @class cg_minus_cga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class cg_minus_cga : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    cg_minus_cga(const StandardModel& SM_i);

    /**
     * @brief A method to compute 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};

/**
 * @class cV_plus_cb
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class cV_plus_cb : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    cV_plus_cb(const StandardModel& SM_i);

    /**
     * @brief A method to compute 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};

/**
 * @class cV_plus_ctau
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class cV_plus_ctau : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    cV_plus_ctau(const StandardModel& SM_i);

    /**
     * @brief A method to compute 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};

/**
 * @class cb_minus_cc
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class cb_minus_cc : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    cb_minus_cc(const StandardModel& SM_i);

    /**
     * @brief A method to compute 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};

/**
 * @class cb_minus_ctau
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class cb_minus_ctau : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    cb_minus_ctau(const StandardModel& SM_i);

    /**
     * @brief A method to compute 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};

/**
 * @class cc_minus_ctau
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class cc_minus_ctau : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    cc_minus_ctau(const StandardModel& SM_i);

    /**
     * @brief A method to compute 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


//  Full signal strengths at e+ e- colliders
//  ----------------------------------------

/**
 * @class mueeZHbb
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to bb}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to bb}@f$ between the 
 * @f$e^+e^- \to ZH, H \to bb@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZHbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeZHbb(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to bb}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to bb}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeZHcc
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to cc}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to cc}@f$ between the 
 * @f$e^+e^- \to ZH, H \to cc@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZHcc : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeZHcc(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to cc}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to cc}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};

/**
 * @class mueeZHss
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to ss}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to ss}@f$ between the 
 * @f$e^+e^- \to ZH, H \to ss@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZHss : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeZHss(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to ss}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to ss}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeZHgg
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to gg}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to gg}@f$ between the 
 * @f$e^+e^- \to ZH, H \to gg@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZHgg : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeZHgg(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to gg}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to gg}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeZHWW
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to WW}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to WW}@f$ between the 
 * @f$e^+e^- \to ZH, H \to WW@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZHWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeZHWW(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);
    
    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to WW}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to WW}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeZHtautau
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to \tau\tau}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to \tau\tau}@f$ between the 
 * @f$e^+e^- \to ZH, H \to \tau\tau@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZHtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeZHtautau(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to \tau\tau}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to \tau\tau}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeZHZZ
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to ZZ}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to ZZ}@f$ between the 
 * @f$e^+e^- \to ZH, H \to ZZ@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZHZZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeZHZZ(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to ZZ}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to ZZ}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeZHZga
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to Z\gamma}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to Z\gamma}@f$ between the 
 * @f$e^+e^- \to ZH, H \to Z\gamma@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZHZga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeZHZga(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to Z\gamma}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to Z\gamma}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeZHgaga
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to \gamma\gamma}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to \gamma\gamma}@f$ between the 
 * @f$e^+e^- \to ZH, H \to \gamma\gamma@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZHgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeZHgaga(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to \gamma\gamma}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to \gamma\gamma}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeZHmumu
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to \mu\mu}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to \mu\mu}@f$ between the 
 * @f$e^+e^- \to ZH, H \to \mu\mu@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZHmumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeZHmumu(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);
    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to \mu\mu}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to \mu\mu}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeZHBRinv
 * @ingroup NewPhysics
 * @brief A class for computing the quantity @f$\mu_{e^+e^- \to ZH, H \to invisible)}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the quantity @f$\mu_{e^+e^- \to ZH, H \to invisible}@f$, i.e.
 * the ratio between the @f$e^+e^- \to ZH@f$ 
 * associated production cross-section in the current model and 
 * the SM, times the total invisible branching ratio.
 */
class mueeZHBRinv : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeZHBRinv(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to invisible}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to invisible}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeZHinv
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to invisible}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to invisible}@f$ between the 
 * @f$e^+e^- \to ZH, H \to invisible@f$ 
 * associated production cross-section in the current model and the SM.
 */
class mueeZHinv : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeZHinv(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to invisible}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to invisible}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};

/**
 * @class mueeWBFbb
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to bb}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to bb}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to bb @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeWBFbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeWBFbb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to bb@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to bb@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeWBFbbPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to bb}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to bb}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to bb @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeWBFbbPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeWBFbbPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to bb@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to bb@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeWBFcc
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to cc}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to cc}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to cc @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeWBFcc : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeWBFcc(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to cc@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to cc@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeWBFgg
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to gg}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to gg}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to gg @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeWBFgg : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeWBFgg(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to gg@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to gg@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeWBFWW
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to WW}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to WW}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to WW @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeWBFWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeWBFWW(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to WW@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to WW@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeWBFtautau
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to \tau\tau}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to \tau\tau}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to \tau\tau @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeWBFtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeWBFtautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to \tau\tau@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to \tau\tau@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeWBFZZ
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to ZZ}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to ZZ}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to ZZ @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeWBFZZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeWBFZZ(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to ZZ@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to ZZ@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeWBFZga
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to Z\gamma}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to Z\gamma}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to Z\gamma @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeWBFZga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeWBFZga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to Z\gamma@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to Z\gamma@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeWBFgaga
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to \gamma\gamma}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to \gamma\gamma}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to \gamma\gamma @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeWBFgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeWBFgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to \gamma\gamma@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to \gamma\gamma@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeWBFmumu
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to \mu\mu}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to \mu\mu}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to \mu\mu @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeWBFmumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeWBFmumu(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to \mu\mu@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to \mu\mu@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeHvvbb
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to bb}@f$,
 * excluding contributions from on-shell @f$Z\to \nu\bar{\nu} @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to bb}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to bb @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeHvvbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeHvvbb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to bb@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to bb@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeHvvbbPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to bb}@f$,
 * excluding contributions from on-shell @f$Z\to \nu\bar{\nu} @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to bb}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to bb @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeHvvbbPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeHvvbbPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to bb@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to bb@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeHvvcc
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to cc}@f$,
 * excluding contributions from on-shell @f$Z\to \nu\bar{\nu} @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to cc}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to cc @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeHvvcc : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeHvvcc(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to cc@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to cc@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeHvvccPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to cc}@f$,
 * excluding contributions from on-shell @f$Z\to \nu\bar{\nu} @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to cc}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to cc @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeHvvccPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeHvvccPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to cc@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to cc@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeHvvss
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to ss}@f$,
 * excluding contributions from on-shell @f$Z\to \nu\bar{\nu} @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to ss}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to ss@f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeHvvss : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeHvvss(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to ss@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to ss@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeHvvssPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to ss@f$,
 * excluding contributions from on-shell @f$Z\to \nu\bar{\nu} @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to ss@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to ss @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeHvvssPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeHvvssPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to ss@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to ss@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeHvvgg
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to gg}@f$,
 * excluding contributions from on-shell @f$Z\to \nu\bar{\nu} @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to gg}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to gg @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeHvvgg : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeHvvgg(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to gg@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to gg@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeHvvggPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to gg}@f$,
 * excluding contributions from on-shell @f$Z\to \nu\bar{\nu} @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to gg}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to gg @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeHvvggPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeHvvggPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to gg@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to bb@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeHvvWW
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to WW}@f$,
 * excluding contributions from on-shell @f$Z\to \nu\bar{\nu} @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to WW}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to WW @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeHvvWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeHvvWW(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to WW@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to WW@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeHvvWWPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to WW}@f$,
 * excluding contributions from on-shell @f$Z\to \nu\bar{\nu} @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to WW}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to WW @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeHvvWWPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeHvvWWPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to WW@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to WW@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeHvvtautau
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to \tau\tau}@f$,
 * excluding contributions from on-shell @f$Z\to \nu\bar{\nu} @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to \tau\tau}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to \tau\tau @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeHvvtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeHvvtautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to \tau\tau@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to \tau\tau@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeHvvtautauPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to \tau\tau}@f$,
 * excluding contributions from on-shell @f$Z\to \nu\bar{\nu} @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to \tau\tau}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to \tau\tau @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeHvvtautauPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeHvvtautauPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to \tau\tau@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to \tau\tau@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeHvvZZ
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to ZZ}@f$,
 * excluding contributions from on-shell @f$Z\to \nu\bar{\nu} @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to ZZ}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to ZZ @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeHvvZZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeHvvZZ(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to ZZ@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to ZZ@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeHvvZZPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to ZZ}@f$,
 * excluding contributions from on-shell @f$Z\to \nu\bar{\nu} @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to ZZ}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to ZZ @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeHvvZZPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeHvvZZPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to ZZ@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to ZZ@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeHvvZga
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to Z\gamma}@f$,
 * excluding contributions from on-shell @f$Z\to \nu\bar{\nu} @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to Z\gamma}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to Z\gamma @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeHvvZga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeHvvZga(const StandardModel& SM_i, const double sqrt_s_i);
    
    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to Z\gamma@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to Z\gamma@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeHvvZgaPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to Z\gamma}@f$,
 * excluding contributions from on-shell @f$Z\to \nu\bar{\nu} @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to Z\gamma}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to Z\gamma @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeHvvZgaPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeHvvZgaPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to Z\gamma@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to Z\gamma@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeHvvgaga
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to \gamma\gamma}@f$,
 * excluding contributions from on-shell @f$Z\to \nu\bar{\nu} @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to \gamma\gamma}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to \gamma\gamma @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeHvvgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeHvvgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to \gamma\gamma@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to \gamma\gamma@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeHvvgagaPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to \gamma\gamma}@f$,
 * excluding contributions from on-shell @f$Z\to \nu\bar{\nu} @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to \gamma\gamma}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to \gamma\gamma @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeHvvgagaPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeHvvgagaPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to \gamma\gamma@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to \gamma\gamma@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeHvvmumu
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to \mu\mu}@f$,
 * excluding contributions from on-shell @f$Z\to \nu\bar{\nu} @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to \mu\mu}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to \mu\mu @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeHvvmumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeHvvmumu(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to \mu\mu@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to \mu\mu@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeHvvmumuPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to \mu\mu}@f$,
 * excluding contributions from on-shell @f$Z\to \nu\bar{\nu} @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to \mu\mu}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to \mu\mu @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeHvvmumuPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeHvvmumuPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to \mu\mu@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to \mu\mu@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};



/**
 * @class mueeZBFbb
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to bb}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to bb}@f$ between the 
 * @f$e^+e^- \to ZH, H \to bb@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZBFbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeZBFbb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to bb}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to bb}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeZBFbbPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to bb}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to bb}@f$ between the 
 * @f$e^+e^- \to ZH, H \to bb@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZBFbbPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueeZBFbbPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to bb}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to bb}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueettHbb
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to bb}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to bb}@f$ between the 
 * @f$e^+e^- \to ZH, H \to bb@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueettHbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueettHbb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to bb}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to bb}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueettHbbPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to bb}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to bb}@f$ between the 
 * @f$e^+e^- \to ZH, H \to bb@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueettHbbPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     * @param[in] Pol_em_i polarization of the electron
     * @param[in] Pol_ep_i polarization of the positron
     */
    mueettHbbPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to bb}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to bb}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


//  Production signal strengths at mu+ mu- colliders
//  ------------------------------------------------

/**
 * @class mummZH
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+\mu^- \to ZH}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+\mu^- \to ZH}@f$ between the 
 * @f$\mu^+\mu^- \to ZH@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummZH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummZH(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+\mu^- \to ZH}@f$ in the current model.
     * @return @f$\mu_{\mu^+\mu^- \to ZH}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class mummHvv
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+\mu^- \to H\nu\bar{\nu}}@f$, 
 * excluding contributions from on-shell @f$Z\to \nu\bar{\nu} @f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+\mu^- \to H\nu\bar{\nu}}@f$ between the 
 * @f$\mu^+\mu^- \to H\nu\bar{\nu}@f$ 
 * production cross-section in the current model and in the Standard Model.
 */
class mummHvv : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHvv(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+\mu^- \to H\nu\bar{\nu}}@f$ in the current model.
     * @return @f$\mu_{\mu^+\mu^- \to H\nu\bar{\nu}}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHmm
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+\mu^- \to H\mu^+\mu^-}@f$, 
 * excluding contributions from on-shell @f$Z\to \mu^+\mu^- @f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+\mu^- \to H\mu^+\mu^-}@f$ between the 
 * @f$\mu^+\mu^- \to H\mu^+\mu^-@f$ 
 * production cross-section in the current model and in the Standard Model.
 */
class mummHmm : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHmm(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+\mu^- \to H\mu^+\mu^-}@f$ in the current model.
     * @return @f$\mu_{\mu^+\mu^- \to H\mu^+\mu^-}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummttH
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{mmttH}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{mmttH}@f$ between the 
 * @f$ \mu^{+}\mu^{-}\to t\bar{t} H @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mummttH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummttH(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{mmttH}@f$ in the current model.
     * @return @f$\mu_{mmttH}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};



//  Full signal strengths at mu+ mu- colliders
//  -------------------------------------------

/**
 * @class mummHbb
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to bb}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to bb}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to bb@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHbb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to bb}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to bb}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHcc
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to cc}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to cc}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to cc@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHcc : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHcc(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to cc}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to cc}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHgg
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to gg}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to gg}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to gg@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHgg : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHgg(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to gg}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to gg}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHWW
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to WW}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to WW}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to WW@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHWW(const StandardModel& SM_i, const double sqrt_s_i);
    
    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to WW}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to WW}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHtautau
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \tau\tau}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \tau\tau}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to \tau\tau@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHtautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to \tau\tau}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to \tau\tau}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHZZ
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to ZZ}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to ZZ}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to ZZ@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHZZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHZZ(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to ZZ}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to ZZ}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHZga
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to Z\gamma}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to Z\gamma}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to Z\gamma@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHZga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHZga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to Z\gamma}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to Z\gamma}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHgaga
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \gamma\gamma}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \gamma\gamma}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to \gamma\gamma@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to \gamma\gamma}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to \gamma\gamma}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHmumu
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \mu\mu}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \mu\mu}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to \mu\mu@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHmumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHmumu(const StandardModel& SM_i, const double sqrt_s_i);
    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to \mu\mu}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to \mu\mu}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHbbNWA
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to bb}@f$, in the narrow width approximation.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to bb}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to bb@f$ 
 * associated production cross-section in the current model and in the Standard Model, in the narrow width approximation.
 */
class mummHbbNWA : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHbbNWA(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to bb}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to bb}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHccNWA
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to cc}@f$, in the narrow width approximation.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to cc}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to cc@f$ 
 * associated production cross-section in the current model and in the Standard Model, in the narrow width approximation.
 */
class mummHccNWA : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHccNWA(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to cc}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to cc}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHggNWA
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to gg}@f$, in the narrow width approximation.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to gg}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to gg@f$ 
 * associated production cross-section in the current model and in the Standard Model, in the narrow width approximation.
 */
class mummHggNWA : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHggNWA(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to gg}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to gg}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHWWNWA
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to WW}@f$, in the narrow width approximation.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to WW}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to WW@f$ 
 * associated production cross-section in the current model and in the Standard Model, in the narrow width approximation.
 */
class mummHWWNWA : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHWWNWA(const StandardModel& SM_i, const double sqrt_s_i);
    
    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to WW}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to WW}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHtautauNWA
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \tau\tau}@f$, in the narrow width approximation.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \tau\tau}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to \tau\tau@f$ 
 * associated production cross-section in the current model and in the Standard Model, in the narrow width approximation.
 */
class mummHtautauNWA : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHtautauNWA(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to \tau\tau}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to \tau\tau}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHZZNWA
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to ZZ}@f$, in the narrow width approximation.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to ZZ}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to ZZ@f$ 
 * associated production cross-section in the current model and in the Standard Model, in the narrow width approximation.
 */
class mummHZZNWA : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHZZNWA(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to ZZ}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to ZZ}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHZgaNWA
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to Z\gamma}@f$, in the narrow width approximation.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to Z\gamma}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to Z\gamma@f$ 
 * associated production cross-section in the current model and in the Standard Model, in the narrow width approximation.
 */
class mummHZgaNWA : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHZgaNWA(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to Z\gamma}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to Z\gamma}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHgagaNWA
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \gamma\gamma}@f$, in the narrow width approximation.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \gamma\gamma}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to \gamma\gamma@f$ 
 * associated production cross-section in the current model and in the Standard Model, in the narrow width approximation.
 */
class mummHgagaNWA : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHgagaNWA(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to \gamma\gamma}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to \gamma\gamma}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHmumuNWA
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \mu\mu}@f$, in the narrow width approximation.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \mu\mu}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to \mu\mu@f$ 
 * associated production cross-section in the current model and in the Standard Model, in the narrow width approximation.
 */
class mummHmumuNWA : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHmumuNWA(const StandardModel& SM_i, const double sqrt_s_i);
    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to \mu\mu}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to \mu\mu}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummZHbb
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to bb}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to bb}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to bb@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummZHbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummZHbb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to bb}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to bb}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummZHcc
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to cc}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to cc}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to cc@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummZHcc : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummZHcc(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to cc}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to cc}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummZHgg
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to gg}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to gg}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to gg@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummZHgg : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummZHgg(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to gg}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to gg}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummZHWW
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to WW}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to WW}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to WW@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummZHWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummZHWW(const StandardModel& SM_i, const double sqrt_s_i);
    
    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to WW}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to WW}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummZHtautau
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \tau\tau}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \tau\tau}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to \tau\tau@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummZHtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummZHtautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to \tau\tau}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to \tau\tau}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummZHZZ
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to ZZ}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to ZZ}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to ZZ@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummZHZZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummZHZZ(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to ZZ}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to ZZ}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummZHZga
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to Z\gamma}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to Z\gamma}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to Z\gamma@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummZHZga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummZHZga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to Z\gamma}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to Z\gamma}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummZHgaga
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \gamma\gamma}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \gamma\gamma}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to \gamma\gamma@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummZHgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummZHgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to \gamma\gamma}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to \gamma\gamma}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummZHmumu
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \mu\mu}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \mu\mu}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to \mu\mu@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummZHmumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummZHmumu(const StandardModel& SM_i, const double sqrt_s_i);
    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to \mu\mu}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to \mu\mu}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};



/**
 * @class mummHvvbb
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to bb}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to bb}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to bb@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHvvbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHvvbb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to bb}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to bb}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHvvcc
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to cc}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to cc}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to cc@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHvvcc : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHvvcc(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to cc}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to cc}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHvvgg
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to gg}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to gg}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to gg@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHvvgg : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHvvgg(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to gg}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to gg}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHvvWW
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to WW}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to WW}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to WW@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHvvWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHvvWW(const StandardModel& SM_i, const double sqrt_s_i);
    
    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to WW}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to WW}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHvvtautau
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \tau\tau}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \tau\tau}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to \tau\tau@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHvvtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHvvtautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to \tau\tau}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to \tau\tau}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHvvZZ
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to ZZ}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to ZZ}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to ZZ@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHvvZZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHvvZZ(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to ZZ}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to ZZ}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHvvZga
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to Z\gamma}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to Z\gamma}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to Z\gamma@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHvvZga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHvvZga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to Z\gamma}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to Z\gamma}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHvvgaga
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \gamma\gamma}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \gamma\gamma}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to \gamma\gamma@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHvvgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHvvgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to \gamma\gamma}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to \gamma\gamma}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHvvmumu
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \mu\mu}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \mu\mu}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to \mu\mu@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHvvmumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHvvmumu(const StandardModel& SM_i, const double sqrt_s_i);
    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to \mu\mu}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to \mu\mu}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};



/**
 * @class mummHmmbb
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to bb}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to bb}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to bb@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHmmbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHmmbb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to bb}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to bb}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHmmcc
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to cc}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to cc}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to cc@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHmmcc : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHmmcc(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to cc}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to cc}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHmmgg
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to gg}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to gg}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to gg@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHmmgg : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHmmgg(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to gg}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to gg}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHmmWW
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to WW}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to WW}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to WW@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHmmWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHmmWW(const StandardModel& SM_i, const double sqrt_s_i);
    
    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to WW}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to WW}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHmmtautau
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \tau\tau}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \tau\tau}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to \tau\tau@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHmmtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHmmtautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to \tau\tau}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to \tau\tau}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHmmZZ
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to ZZ}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to ZZ}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to ZZ@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHmmZZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHmmZZ(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to ZZ}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to ZZ}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHmmZga
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to Z\gamma}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to Z\gamma}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to Z\gamma@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHmmZga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHmmZga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to Z\gamma}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to Z\gamma}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHmmgaga
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \gamma\gamma}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \gamma\gamma}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to \gamma\gamma@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHmmgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHmmgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to \gamma\gamma}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to \gamma\gamma}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummHmmmumu
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \mu\mu}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \mu\mu}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to \mu\mu@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummHmmmumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummHmmmumu(const StandardModel& SM_i, const double sqrt_s_i);
    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to \mu\mu}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to \mu\mu}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummttHbb
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to bb}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to bb}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to bb@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummttHbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummttHbb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to bb}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to bb}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummttHcc
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to cc}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to cc}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to cc@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummttHcc : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummttHcc(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to cc}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to cc}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummttHgg
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to gg}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to gg}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to gg@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummttHgg : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummttHgg(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to gg}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to gg}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummttHWW
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to WW}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to WW}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to WW@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummttHWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummttHWW(const StandardModel& SM_i, const double sqrt_s_i);
    
    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to WW}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to WW}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummttHtautau
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \tau\tau}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \tau\tau}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to \tau\tau@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummttHtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummttHtautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to \tau\tau}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to \tau\tau}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummttHZZ
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to ZZ}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to ZZ}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to ZZ@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummttHZZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummttHZZ(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to ZZ}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to ZZ}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummttHZga
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to Z\gamma}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to Z\gamma}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to Z\gamma@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummttHZga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummttHZga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to Z\gamma}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to Z\gamma}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummttHgaga
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \gamma\gamma}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \gamma\gamma}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to \gamma\gamma@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummttHgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummttHgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to \gamma\gamma}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to \gamma\gamma}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummttHmumu
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \mu\mu}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{\mu^+ \mu^- \to H, H \to \mu\mu}@f$ between the 
 * @f$\mu^+ \mu^- \to H, H \to \mu\mu@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mummttHmumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mummttHmumu(const StandardModel& SM_i, const double sqrt_s_i);
    /**
     * @brief A method to compute the value of @f$\mu_{\mu^+ \mu^- \to H, H \to \mu\mu}@f$ in the current model.
     * @return @f$\mu_{\mu^+ \mu^- \to H, H \to \mu\mu}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};




//  Full signal strengths at ep colliders
//  -------------------------------------

/**
 * @class muepWBFbb
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^- p \to H (WBF), H \to bb}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^-p \to H (WBF), H \to bb}@f$ between the 
 * @f$ e^{-}p\to \nu j H, H \to bb @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class muepWBFbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muepWBFbb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^- p \to H (WBF), H \to bb}@f$ in the current model.
     * @return @f$\mu_{e^- p \to H (WBF), H \to bb}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muepWBFcc
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^- p \to H (WBF), H \to cc}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^-p \to H (WBF), H \to cc}@f$ between the 
 * @f$ e^{-}p\to \nu j H, H \to cc @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class muepWBFcc : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muepWBFcc(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^- p \to H (WBF), H \to cc}@f$ in the current model.
     * @return @f$\mu_{e^- p \to H (WBF), H \to cc}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muepWBFgg
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^- p \to H (WBF), H \to gg}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^-p \to H (WBF), H \to gg}@f$ between the 
 * @f$ e^{-}p\to \nu j H, H \to gg @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class muepWBFgg : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muepWBFgg(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^- p \to H (WBF), H \to gg}@f$ in the current model.
     * @return @f$\mu_{e^- p \to H (WBF), H \to gg}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muepWBFWW2l2v
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^- p \to H (WBF), H \to WW\to 2l2v}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^-p \to H (WBF), H \to WW\to 2l2v}@f$ between the 
 * @f$ e^{-}p\to \nu j H, H \to WW\to 2l2v @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class muepWBFWW2l2v : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muepWBFWW2l2v(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^- p \to H (WBF), H \to WW\to 2l2v}@f$ in the current model.
     * @return @f$\mu_{e^- p \to H (WBF), H \to WW\to 2l2v}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muepWBFZZ4l
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^- p \to H (WBF), H \to ZZ\to 4l}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^-p \to H (WBF), H \to ZZ\to 4l}@f$ between the 
 * @f$ e^{-}p\to \nu j H, H \to ZZ\to 4l @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class muepWBFZZ4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muepWBFZZ4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^- p \to H (WBF), H \to ZZ\to 4l}@f$ in the current model.
     * @return @f$\mu_{e^- p \to H (WBF), H \to ZZ\to 4l}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muepWBFgaga
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^- p \to H (WBF), H \to \gamma\gamma}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^-p \to H (WBF), H \to \gamma\gamma}@f$ between the 
 * @f$ e^{-}p\to \nu j H, H \to \gamma\gamma @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class muepWBFgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muepWBFgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^- p \to H (WBF), H \to \gamma\gamma}@f$ in the current model.
     * @return @f$\mu_{e^- p \to H (WBF), H \to \gamma\gamma}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muepWBFtautau
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^- p \to H (WBF), H \to \tau\tau}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^-p \to H (WBF), H \to \tau\tau}@f$ between the 
 * @f$ e^{-}p\to \nu j H, H \to \tau\tau @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class muepWBFtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muepWBFtautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^- p \to H (WBF), H \to \tau\tau}@f$ in the current model.
     * @return @f$\mu_{e^- p \to H (WBF), H \to \tau\tau}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muepZBFbb
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^- p \to H (ZBF), H \to bb}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^- p \to H (ZBF), H \to bb}@f$ between the 
 * @f$ e^{-}p\to e^{-} j H, H \to bb @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class muepZBFbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muepZBFbb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^- p \to H (ZBF), H \to bb}@f$ in the current model.
     * @return @f$\mu_{e^- p \to H (ZBF), H \to bb}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muepZBFcc
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^- p \to H (ZBF), H \to cc}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^- p \to H (ZBF), H \to cc}@f$ between the 
 * @f$ e^{-}p\to e^{-} j H, H \to cc @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class muepZBFcc : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muepZBFcc(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^- p \to H (ZBF), H \to cc}@f$ in the current model.
     * @return @f$\mu_{e^- p \to H (ZBF), H \to cc}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muepZBFgg
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^- p \to H (ZBF), H \to gg}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^- p \to H (ZBF), H \to gg}@f$ between the 
 * @f$ e^{-}p\to e^{-} j H, H \to gg @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class muepZBFgg : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muepZBFgg(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^- p \to H (ZBF), H \to gg}@f$ in the current model.
     * @return @f$\mu_{e^- p \to H (ZBF), H \to gg}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muepZBFWW2l2v
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^- p \to H (ZBF), H \to WW\to 2l2v}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^- p \to H (ZBF), H \to WW\to 2l2v}@f$ between the 
 * @f$ e^{-}p\to e^{-} j H, H \to WW\to 2l2v @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class muepZBFWW2l2v : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muepZBFWW2l2v(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^- p \to H (ZBF), H \to WW\to 2l2v}@f$ in the current model.
     * @return @f$\mu_{e^- p \to H (ZBF), H \to WW\to 2l2v}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muepZBFZZ4l
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^- p \to H (ZBF), H \to ZZ\to 4l}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^- p \to H (ZBF), H \to ZZ\to 4l}@f$ between the 
 * @f$ e^{-}p\to e^{-} j H, H \to ZZ\to 4l @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class muepZBFZZ4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muepZBFZZ4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^- p \to H (ZBF), H \to ZZ\to 4l}@f$ in the current model.
     * @return @f$\mu_{e^- p \to H (ZBF), H \to ZZ\to 4l}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muepZBFgaga
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^- p \to H (ZBF), H \to \gamma\gamma}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^- p \to H (ZBF), H \to \gamma\gamma}@f$ between the 
 * @f$ e^{-}p\to e^{-} j H, H \to \gamma\gamma @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class muepZBFgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muepZBFgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^- p \to H (ZBF), H \to \gamma\gamma}@f$ in the current model.
     * @return @f$\mu_{e^- p \to H (ZBF), H \to \gamma\gamma}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muepZBFtautau
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^- p \to H (ZBF), H \to \tau\tau}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^- p \to H (ZBF), H \to \tau\tau}@f$ between the 
 * @f$ e^{-}p\to e^{-} j H, H \to \tau\tau @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class muepZBFtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muepZBFtautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^- p \to H (ZBF), H \to \tau\tau}@f$ in the current model.
     * @return @f$\mu_{e^- p \to H (ZBF), H \to \tau\tau}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


// -----------------------------------------------------------------------------
// STXS bins
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// Stage 0
// -----------------------------------------------------------------------------

/**
 * @class STXS_0_qqH
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$pp \to H qq@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$pp \to H qq@f$.
 */
class STXS_0_qqH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS_0_qqH(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


// -----------------------------------------------------------------------------
// Stage 1
// -----------------------------------------------------------------------------

/**
 * @class STXSggH_VBFtopo_j3v_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$gg \to H@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$gg \to H@f$.
 */
class STXSggH_VBFtopo_j3v_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSggH_VBFtopo_j3v_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSggH_VBFtopo_j3_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$gg \to H@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$gg \to H@f$.
 */
class STXSggH_VBFtopo_j3_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSggH_VBFtopo_j3_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSggH0j4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$gg \to H@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$gg \to H@f$.
 */
class STXSggH0j4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSggH0j4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSggH1j_pTH_0_60_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$gg \to H@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$gg \to H@f$.
 */
class STXSggH1j_pTH_0_60_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSggH1j_pTH_0_60_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSggH1j_pTH_60_120_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$gg \to H@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$gg \to H@f$.
 */
class STXSggH1j_pTH_60_120_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSggH1j_pTH_60_120_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSggH1j_pTH_120_200_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$gg \to H@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$gg \to H@f$.
 */
class STXSggH1j_pTH_120_200_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSggH1j_pTH_120_200_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSggH1j_pTH_200_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$gg \to H@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$gg \to H@f$.
 */
class STXSggH1j_pTH_200_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSggH1j_pTH_200_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSggH2j_pTH_0_200_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$gg \to H@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$gg \to H@f$.
 */
class STXSggH2j_pTH_0_200_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSggH2j_pTH_0_200_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSggH2j_pTH_0_60_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$gg \to H@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$gg \to H@f$.
 */
class STXSggH2j_pTH_0_60_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSggH2j_pTH_0_60_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSggH2j_pTH_60_120_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$gg \to H@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$gg \to H@f$.
 */
class STXSggH2j_pTH_60_120_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSggH2j_pTH_60_120_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSggH2j_pTH_120_200_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$gg \to H@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$gg \to H@f$.
 */
class STXSggH2j_pTH_120_200_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSggH2j_pTH_120_200_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSggH2j_pTH_200_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$gg \to H@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$gg \to H@f$.
 */
class STXSggH2j_pTH_200_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSggH2j_pTH_200_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};



/**
 * @class STXSqqHqq_VBFtopo_Rest_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H qq@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H qq@f$.
 */
class STXSqqHqq_VBFtopo_Rest_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHqq_VBFtopo_Rest_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};



/**
 * @class STXSqqHqq_VBFtopo_j3v_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H qq@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H qq@f$.
 */
class STXSqqHqq_VBFtopo_j3v_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHqq_VBFtopo_j3v_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSqqHqq_VBFtopo_j3_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H qq@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H qq@f$.
 */
class STXSqqHqq_VBFtopo_j3_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHqq_VBFtopo_j3_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSqqHqq_nonVHtopo_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H qq@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H qq@f$.
 */
class STXSqqHqq_nonVHtopo_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHqq_nonVHtopo_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSqqHqq_VHtopo_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H qq@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H qq@f$.
 */
class STXSqqHqq_VHtopo_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHqq_VHtopo_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSqqHqq_Rest_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H qq@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H qq@f$.
 */
class STXSqqHqq_Rest_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHqq_Rest_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSqqHqq_pTj_200_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H qq@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H qq@f$.
 */
class STXSqqHqq_pTj_200_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHqq_pTj_200_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSqqHlv_pTV_0_250_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H \ell \nu@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H \ell \nu@f$.
 */
class STXSqqHlv_pTV_0_250_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHlv_pTV_0_250_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSqqHlv_pTV_0_150_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H \ell \nu@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H \ell \nu@f$.
 */
class STXSqqHlv_pTV_0_150_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHlv_pTV_0_150_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSqqHlv_pTV_150_250_0j_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H \ell \nu@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H \ell \nu@f$.
 */
class STXSqqHlv_pTV_150_250_0j_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHlv_pTV_150_250_0j_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSqqHlv_pTV_150_250_1j_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H \ell \nu@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H \ell \nu@f$.
 */
class STXSqqHlv_pTV_150_250_1j_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHlv_pTV_150_250_1j_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSqqHlv_pTV_250_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H \ell \nu@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H \ell \nu@f$.
 */
class STXSqqHlv_pTV_250_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHlv_pTV_250_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSqqHll_pTV_0_150_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H \ell \ell@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H \ell \ell@f$.
 */
class STXSqqHll_pTV_0_150_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHll_pTV_0_150_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSqqHll_pTV_150_250_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H \ell \ell@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H \ell \ell@f$.
 */
class STXSqqHll_pTV_150_250_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHll_pTV_150_250_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSqqHll_pTV_150_250_0j_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H \ell \ell@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H \ell \ell@f$.
 */
class STXSqqHll_pTV_150_250_0j_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHll_pTV_150_250_0j_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};



/**
 * @class STXSqqHll_pTV_150_250_1j_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H \ell \ell@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H \ell \ell@f$.
 */
class STXSqqHll_pTV_150_250_1j_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHll_pTV_150_250_1j_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};



/**
 * @class STXSqqHll_pTV_250_4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H \ell \ell@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H \ell \ell@f$.
 */
class STXSqqHll_pTV_250_4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHll_pTV_250_4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};



/**
 * @class STXSttHtH4l
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ttH + tH@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ttH + tH@f$.
 */
class STXSttHtH4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSttHtH4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};



/**
 * @class STXSqqHlv_pTV_0_250_bb
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H \ell \nu@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H \ell \nu@f$.
 */
class STXSqqHlv_pTV_0_250_bb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHlv_pTV_0_250_bb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSqqHlv_pTV_0_150_bb
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H \ell \nu@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H \ell \nu@f$.
 */
class STXSqqHlv_pTV_0_150_bb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHlv_pTV_0_150_bb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSqqHlv_pTV_150_250_0j_bb
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H \ell \nu@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H \ell \nu@f$.
 */
class STXSqqHlv_pTV_150_250_0j_bb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHlv_pTV_150_250_0j_bb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSqqHlv_pTV_150_250_1j_bb
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H \ell \nu@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H \ell \nu@f$.
 */
class STXSqqHlv_pTV_150_250_1j_bb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHlv_pTV_150_250_1j_bb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSqqHlv_pTV_250_bb
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H \ell \nu@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H \ell \nu@f$.
 */
class STXSqqHlv_pTV_250_bb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHlv_pTV_250_bb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSqqHll_pTV_0_150_bb
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H \ell \ell@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H \ell \ell@f$.
 */
class STXSqqHll_pTV_0_150_bb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHll_pTV_0_150_bb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSqqHll_pTV_150_250_bb
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H \ell \ell@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H \ell \ell@f$.
 */
class STXSqqHll_pTV_150_250_bb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHll_pTV_150_250_bb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSqqHll_pTV_150_250_0j_bb
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H \ell \ell@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H \ell \ell@f$.
 */
class STXSqqHll_pTV_150_250_0j_bb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHll_pTV_150_250_0j_bb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};



/**
 * @class STXSqqHll_pTV_150_250_1j_bb
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H \ell \ell@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H \ell \ell@f$.
 */
class STXSqqHll_pTV_150_250_1j_bb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHll_pTV_150_250_1j_bb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};



/**
 * @class STXSqqHll_pTV_250_bb
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to H \ell \ell@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to H \ell \ell@f$.
 */
class STXSqqHll_pTV_250_bb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSqqHll_pTV_250_bb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSWHqqHqq_VBFtopo_j3v_2b
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to WH  \to H qq@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to WH  \to H qq@f$.
 */
class STXSWHqqHqq_VBFtopo_j3v_2b : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSWHqqHqq_VBFtopo_j3v_2b(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};



/**
 * @class STXSWHqqHqq_VBFtopo_j3_2b
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to WH  \to H qq@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to WH  \to H qq@f$.
 */
class STXSWHqqHqq_VBFtopo_j3_2b : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSWHqqHqq_VBFtopo_j3_2b(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};



/**
 * @class STXSWHqqHqq_VH2j_2b
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to WH  \to H qq@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to WH  \to H qq@f$.
 */
class STXSWHqqHqq_VH2j_2b : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSWHqqHqq_VH2j_2b(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};



/**
 * @class STXSWHqqHqq_Rest_2b
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to WH  \to H qq@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to WH  \to H qq@f$.
 */
class STXSWHqqHqq_Rest_2b : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSWHqqHqq_Rest_2b(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};



/**
 * @class STXSWHqqHqq_pTj1_200_2b
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to WH  \to H qq@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to WH  \to H qq@f$.
 */
class STXSWHqqHqq_pTj1_200_2b : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSWHqqHqq_pTj1_200_2b(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSZHqqHqq_VBFtopo_j3v_2b
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to ZH  \to H qq@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to ZH  \to H qq@f$.
 */
class STXSZHqqHqq_VBFtopo_j3v_2b : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSZHqqHqq_VBFtopo_j3v_2b(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};



/**
 * @class STXSZHqqHqq_VBFtopo_j3_2b
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to ZH  \to H qq@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to ZH  \to H qq@f$.
 */
class STXSZHqqHqq_VBFtopo_j3_2b : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSZHqqHqq_VBFtopo_j3_2b(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSZHqqHqq_VH2j_2b
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to ZH  \to H qq@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to ZH  \to H qq@f$.
 */
class STXSZHqqHqq_VH2j_2b : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSZHqqHqq_VH2j_2b(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class STXSZHqqHqq_Rest_2b
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to ZH  \to H qq@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to ZH  \to H qq@f$.
 */
class STXSZHqqHqq_Rest_2b : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSZHqqHqq_Rest_2b(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};



/**
 * @class STXSZHqqHqq_pTj1_200_2b
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$qq \to ZH  \to H qq@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$qq \to ZH  \to H qq@f$.
 */
class STXSZHqqHqq_pTj1_200_2b : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXSZHqqHqq_pTj1_200_2b(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};



// -----------------------------------------------------------------------------
// Stage 1.2
// -----------------------------------------------------------------------------


/**
 * @class STXS12_ggH_pTH200_300_Nj01
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_pTH200_300_Nj01 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_pTH200_300_Nj01(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};


/**
 * @class STXS12_ggH_pTH300_450_Nj01
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_pTH300_450_Nj01 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_pTH300_450_Nj01(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};


/**
 * @class STXS12_ggH_pTH450_650_Nj01
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_pTH450_650_Nj01 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_pTH450_650_Nj01(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};


/**
 * @class STXS12_ggH_pTH650_Inf_Nj01
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_pTH650_Inf_Nj01 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_pTH650_Inf_Nj01(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};



/**
 * @class STXS12_ggH_pTH10_Inf_Nj0
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_pTH10_Inf_Nj0 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_pTH10_Inf_Nj0(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};





/**
 * @class STXS12_ggH_pTH0_10_Nj0
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_pTH0_10_Nj0 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_pTH0_10_Nj0(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};


/**
 * @class STXS12_ggH_pTH10_200_Nj0
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_pTH10_200_Nj0 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_pTH10_200_Nj0(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};


/**
 * @class STXS12_ggH_pTH0_200_Nj0
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_pTH0_200_Nj0 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_pTH0_200_Nj0(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};





/**
 * @class STXS12_ggH_mjj0_350_pTH0_60_Nj1
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_mjj0_350_pTH0_60_Nj1 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_mjj0_350_pTH0_60_Nj1(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};




/**
 * @class STXS12_ggH_pTH0_60_Nj1
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_pTH0_60_Nj1 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_pTH0_60_Nj1(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};


/**
 * @class STXS12_ggH_pTH60_120_Nj1
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_pTH60_120_Nj1 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_pTH60_120_Nj1(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};


/**
 * @class STXS12_ggH_pTH120_200_Nj1
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_pTH120_200_Nj1 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_pTH120_200_Nj1(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};


/**
 * @class STXS12_ggH_mjj0_350_pTH0_60_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_mjj0_350_pTH0_60_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_mjj0_350_pTH0_60_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};


/**
 * @class STXS12_ggH_mjj0_350_pTH60_120_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_mjj0_350_pTH60_120_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_mjj0_350_pTH60_120_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};




/**
 * @class STXS12_ggH_mjj0_350_pTH0_120_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_mjj0_350_pTH0_120_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_mjj0_350_pTH0_120_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};




/**
 * @class STXS12_ggH_mjj0_350_pTH120_200_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_mjj0_350_pTH120_200_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_mjj0_350_pTH120_200_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};


/**
 * @class STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};


/**
 * @class STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};


/**
 * @class STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};


/**
 * @class STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};





/**
 * @class STXS12_ggH_mjj350_Inf_pTH0_200_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_mjj350_Inf_pTH0_200_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_mjj350_Inf_pTH0_200_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
}; //AG:added




/**
 * @class STXS12_ggH_pTH0_200_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_pTH0_200_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_pTH0_200_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
}; //VM:added




/**
 * @class STXS12_ggH_pTH200_Inf
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_pTH200_Inf : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_pTH200_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
}; //VM:added



/**
 * @class STXS12_ggH_pTH300_Inf
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_pTH300_Inf : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_pTH300_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
}; //VM:added





/**
 * @class STXS12_ggH_pTH200_300
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_pTH200_300 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_pTH200_300(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
}; //VM:added





/**
 * @class STXS12_ggH_pTH300_450
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_pTH300_450 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_pTH300_450(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
}; //VM:added




/**
 * @class STXS12_ggH_pTH450_600
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_pTH450_650 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_pTH450_650(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
}; //VM:added




/**
 * @class STXS12_ggH_pTH650_Inf
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_pTH650_Inf : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_pTH650_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
}; //VM:added




/**
 * @class STXS12_ggH_pTH450_Inf
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggH_pTH450_Inf : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggH_pTH450_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
}; //VM:added




/**
 * @class STXS12_ggHpttH
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggHpttH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggHpttH(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
}; //AG:added


/**
 * @class STXS12_ggHll_pTV0_75
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggHll_pTV0_75 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggHll_pTV0_75(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};


/**
 * @class STXS12_ggHll_pTV75_150
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggHll_pTV75_150 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggHll_pTV75_150(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};


/**
 * @class STXS12_ggHll_pTV150_250_Nj0
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggHll_pTV150_250_Nj0 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggHll_pTV150_250_Nj0(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};


/**
 * @class STXS12_ggHll_pTV150_250_Nj1
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggHll_pTV150_250_Nj1 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggHll_pTV150_250_Nj1(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};


/**
 * @class STXS12_ggHll_pTV250_Inf
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ggHll_pTV250_Inf : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ggHll_pTV250_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};



/**
 * @class STXS12_qqHqq_Nj0
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHqq_Nj0 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHqq_Nj0(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};

/**
 * @class STXS12_qqHqq_Nj1
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHqq_Nj1 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHqq_Nj1(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};

/**
 * @class STXS12_qqHqq_mjj0_60_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHqq_mjj0_60_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHqq_mjj0_60_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};


/**
 * @class STXS12_qqHqq_mjj60_120_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHqq_mjj60_120_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHqq_mjj60_120_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};





/**
 * @class STXS12_qqHqq_VH_veto_Nj01
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHqq_VH_veto_Nj01 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHqq_VH_veto_Nj01(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};




/**
 * @class STXS12_qqHqq_VH_had_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHqq_VH_had_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHqq_VH_had_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};

/**
 * @class STXS12_qqHqq_mjj120_350_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHqq_mjj120_350_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHqq_mjj120_350_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};

/**
 * @class STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};

/**
 * @class STXS12_qqHqq_mjj350_1000_pTH200_Inf_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHqq_mjj350_1000_pTH200_Inf_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHqq_mjj350_1000_pTH200_Inf_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};      //AG:added

/**
 * @class STXS12_qqHqq_mjj1000_Inf_pTH200_Inf_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHqq_mjj1000_Inf_pTH200_Inf_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHqq_mjj1000_Inf_pTH200_Inf_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};      //AG:added

/**
 * @class STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};

/**
 * @class STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};

/**
 * @class STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};

/**
 * @class STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};

/**
 * @class STXS12_qqHqq_mjj350_700_pTH0_200_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHqq_mjj350_700_pTH0_200_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHqq_mjj350_700_pTH0_200_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};      //AG:added

/**
 * @class STXS12_qqHqq_mjj700_1000_pTH0_200_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHqq_mjj700_1000_pTH0_200_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHqq_mjj700_1000_pTH0_200_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};      //AG:added

/**
 * @class STXS12_qqHqq_mjj1000_1500_pTH0_200_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHqq_mjj1000_1500_pTH0_200_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHqq_mjj1000_1500_pTH0_200_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};      //AG:added

/**
 * @class STXS12_qqHqq_mjj1500_Inf_pTH0_200_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHqq_mjj1500_Inf_pTH0_200_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHqq_mjj1500_Inf_pTH0_200_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};      //AG:added

/**
 * @class STXS12_qqHqq_mjj1000_Inf_pTH0_200_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHqq_mjj1000_Inf_pTH0_200_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHqq_mjj1000_Inf_pTH0_200_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};      //AG:added

/**
 * @class STXS12_qqHqq_mjj350_Inf_Nj2
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHqq_mjj350_Inf_Nj2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHqq_mjj350_Inf_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};      //AG:added



/**
 * @class STXS12_qqHlv_pTV0_75
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHlv_pTV0_75 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHlv_pTV0_75(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};

/**
 * @class STXS12_qqHlv_pTV75_150
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHlv_pTV75_150 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHlv_pTV75_150(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};

/**
 * @class STXS12_qqHlv_pTV150_250_Nj0
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHlv_pTV150_250_Nj0 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHlv_pTV150_250_Nj0(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};

/**
 * @class STXS12_qqHlv_pTV150_250_Nj1
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHlv_pTV150_250_Nj1 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHlv_pTV150_250_Nj1(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};

/**
 * @class STXS12_qqHlv_pTV250_Inf
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHlv_pTV250_Inf : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHlv_pTV250_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};

/**
 * @class STXS12_qqHlv_pTV250_Inf
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHlv_pTV0_150 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHlv_pTV0_150(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};              //AG:added

/**
 * @class STXS12_qqHlv_pTV250_400
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHlv_pTV250_400 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHlv_pTV250_400(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};            //AG:added

/**
 * @class STXS12_qqHlv_pTV400_Inf
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHlv_pTV400_Inf : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHlv_pTV400_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};            //AG:added

/**
 * @class STXS12_qqHlv_pTV150_Inf
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHlv_pTV150_Inf : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHlv_pTV150_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};            //AG:added



/**
 * @class STXS12_qqHll_pTV0_75
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHll_pTV0_75 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHll_pTV0_75(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};

/**
 * @class STXS12_qqHll_pTV75_150
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHll_pTV75_150 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHll_pTV75_150(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};

/**
 * @class STXS12_qqHll_pTV0_150
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHll_pTV0_150 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHll_pTV0_150(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};              //AG:added

/**
 * @class STXS12_qqHll_pTV150_250_Nj0
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHll_pTV150_250_Nj0 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHll_pTV150_250_Nj0(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};

/**
 * @class STXS12_qqHll_pTV150_250_Nj1
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHll_pTV150_250_Nj1 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHll_pTV150_250_Nj1(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};

/**
 * @class STXS12_qqHll_pTV250_Inf
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHll_pTV250_Inf : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHll_pTV250_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};

/**
 * @class STXS12_qqHll_pTV250_400
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHll_pTV250_400 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHll_pTV250_400(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};            //AG:added

/**
 * @class STXS12_qqHll_pTV400_Inf
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHll_pTV400_Inf : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHll_pTV400_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};            //AG:added

/**
 * @class STXS12_qqHll_pTV150_Inf
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_qqHll_pTV150_Inf : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_qqHll_pTV150_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};            //AG:added


/**
 * @class STXS12_VHlep
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_VHlep : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_VHlep(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};            //AG:added



/**
 * @class STXS12_ttH_pTH0_60
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ttH_pTH0_60 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ttH_pTH0_60(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};

/**
 * @class STXS12_ttH_pTH60_120
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ttH_pTH60_120 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ttH_pTH60_120(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};

/**
 * @class STXS12_ttH_pTH0_120
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ttH_pTH0_120 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ttH_pTH0_120(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};

/**
 * @class STXS12_ttH_pTH120_200
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ttH_pTH120_200 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ttH_pTH120_200(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};

/**
 * @class STXS12_ttH_pTH200_300
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ttH_pTH200_300 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ttH_pTH200_300(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};

/**
 * @class STXS12_ttH_pTH300_Inf
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ttH_pTH300_Inf : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ttH_pTH300_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};

/**
 * @class STXS12_ttH_pTH300_450
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ttH_pTH300_450 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ttH_pTH300_450(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};              //AG:added

/**
 * @class STXS12_ttH_pTH450_Inf
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ttH_pTH450_Inf : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ttH_pTH450_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};              //AG:added

/**
 * @class STXS12_ttH_pTH300_Inf_add
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ttH_pTH300_Inf_add : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ttH_pTH300_Inf_add(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};              //AG:added

/**
 * @class STXS12_ttH
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_ttH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_ttH(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};              //AG:added



/**
 * @class STXS12_tH
 * @ingroup NewPhysics
 * @brief A class for computing the STXS bin @f$ @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the STXS bin @f$ @f$.
 */
class STXS12_tH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    STXS12_tH(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i);

    /**
     * @brief A method to compute the value of the STXS bin in the current model.
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    unsigned int fstate;
};



//-----------------------------------------------------------------------------------------
//-- Special Hadron collider signal strengths with separate full TH unc U(prod x decay) ---
//-----------------------------------------------------------------------------------------

/**
 * @class muTHUggHgaga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUggHgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUggHgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUVBFHgaga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUVBFHgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUVBFHgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUZHgaga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUZHgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUZHgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUWHgaga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUWHgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUWHgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUVHgaga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUVHgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUVHgaga(const StandardModel& SM_i, const double sqrt_s_i);
    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUttHgaga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUttHgaga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUttHgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUggHZga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUggHZga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUggHZga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};



/**
 * @class muTHUggHZgamumu
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUggHZgamumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUggHZgamumu(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muTHUVBFHZga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUVBFHZga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUVBFHZga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUZHZga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUZHZga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUZHZga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();
private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUWHZga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUWHZga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUWHZga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUVHZga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUVHZga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUVHZga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUttHZga
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUttHZga : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUttHZga(const StandardModel& SM_i, const double sqrt_s_i);
    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUggHZZ
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUggHZZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUggHZZ(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUVBFHZZ
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUVBFHZZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUVBFHZZ(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUZHZZ
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUZHZZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUZHZZ(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUWHZZ
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUWHZZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUWHZZ(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUVHZZ
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUVHZZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUVHZZ(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUttHZZ
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUttHZZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUttHZZ(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muTHUggHZZ4l
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUggHZZ4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUggHZZ4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};



/**
 * @class muTHUggHZZ4mu
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUggHZZ4mu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUggHZZ4mu(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muTHUVBFHZZ4l
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUVBFHZZ4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUVBFHZZ4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUZHZZ4l
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUZHZZ4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUZHZZ4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUWHZZ4l
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUWHZZ4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUWHZZ4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUVHZZ4l
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUVHZZ4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUVHZZ4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUttHZZ4l
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUttHZZ4l : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUttHZZ4l(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muTHUggHWW
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUggHWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUggHWW(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUVBFHWW
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUVBFHWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUVBFHWW(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUZHWW
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUZHWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUZHWW(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUWHWW
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUWHWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUWHWW(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUVHWW
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUVHWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUVHWW(const StandardModel& SM_i, const double sqrt_s_i);
    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUttHWW
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUttHWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUttHWW(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muTHUggHWW2l2v
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUggHWW2l2v : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUggHWW2l2v(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUVBFHWW2l2v
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUVBFHWW2l2v : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUVBFHWW2l2v(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUZHWW2l2v
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUZHWW2l2v : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUZHWW2l2v(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUWHWW2l2v
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUWHWW2l2v : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUWHWW2l2v(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUVHWW2l2v
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUVHWW2l2v : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUVHWW2l2v(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUttHWW2l2v
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUttHWW2l2v : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUttHWW2l2v(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muTHUggHmumu
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUggHmumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUggHmumu(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUVBFHmumu
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUVBFHmumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUVBFHmumu(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUZHmumu
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUZHmumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUZHmumu(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUWHmumu
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUWHmumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUWHmumu(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUVHmumu
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUVHmumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUVHmumu(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUttHmumu
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUttHmumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUttHmumu(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUggHtautau
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUggHtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUggHtautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUVBFHtautau
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUVBFHtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUVBFHtautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUZHtautau
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUZHtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUZHtautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUWHtautau
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUWHtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUWHtautau(const StandardModel& SM_i, const double sqrt_s_i);
    
    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUVHtautau
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUVHtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUVHtautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUttHtautau
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUttHtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUttHtautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUggHbb
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUggHbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUggHbb(const StandardModel& SM_i, const double sqrt_s_i);
    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUVBFHbb
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUVBFHbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUVBFHbb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUZHbb
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUZHbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUZHbb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUWHbb
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUWHbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUWHbb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUVHbb
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUVHbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUVHbb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muTHUttHbb
 * @ingroup NewPhysics
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class muTHUttHbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUttHbb(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief 
     * @return 
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muTHUVBFBRinv
 * @ingroup NewPhysics
 * @brief A class for computing the quantity @f$\mu_{VBF, H \to invisible}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the quantity @f$\mu_{VBF, H \to invisible}@f$, i.e.
 * the ratio between the @f$pp \to jjH@f$ 
 * production cross-section in the current model and 
 * the SM, times the total invisible branching ratio.
 */
class muTHUVBFBRinv : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUVBFBRinv(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{VBF}\times BR(H \to invisible}@f$ in the current model.
     * @return @f$\mu_{VBF}\times BR(H \to invisible}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muTHUVBFHinv
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{VBF, H \to invisible}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{VBF, H \to invisible}@f$ between the 
 * @f$pp \to jjH, H \to invisible@f$ 
 * production cross-section in the current model and the SM.
 */
class muTHUVBFHinv : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUVBFHinv(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{VBF, H \to invisible}@f$ in the current model.
     * @return @f$\mu_{VBF, H \to invisible}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muTHUVHBRinv
 * @ingroup NewPhysics
 * @brief A class for computing the quantity @f$\mu_{pp \to VH, H \to invisible}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the quantity @f$\mu_{pp \to VH, H \to invisible}@f$, i.e.
 * the ratio between the @f$pp \to VH@f$ 
 * associated production cross-section in the current model and 
 * the SM, times the total invisible branching ratio.
 */
class muTHUVHBRinv : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUVHBRinv(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{pp \to VH}\times BR(H \to invisible}@f$ in the current model.
     * @return @f$\mu_{pp \to VH}\times BR(H \to invisible}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muTHUVHinv
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{pp \to VH, H \to invisible}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{pp \to VH, H \to invisible}@f$ between the 
 * @f$pp \to VH, H \to invisible@f$ 
 * associated production cross-section in the current model and the SM.
 */
class muTHUVHinv : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muTHUVHinv(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{pp \to VH, H \to invisible}@f$ in the current model.
     * @return @f$\mu_{pp \to VH, H \to invisible}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};







#endif	/* HIGGSTHOBSERVABLES_H */

