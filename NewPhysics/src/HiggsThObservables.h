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
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
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
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
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
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
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
     */
    mueeZH(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class mueeZHPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH}@f$ between the 
 * @f$e^+e^- \to ZH@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZHPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
     */
    mueeZHPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

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
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
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
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
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
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
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
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
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
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
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
 * @class BrHtoWW2l2vRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to WW\to 2l2\nu)@f$ with @f$l=e,\mu@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to WW\to 2l2\nu)@f$.
 * in the current model and in the Standard Model
 */
class BrHtoWW2l2vRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtoWW2l2vRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to WW\to 2l2\nu)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to WW\to 2l2\nu)@f$
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
 * @class BrHtoZZ4lRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to ZZ\to 4l)@f$ with @f$l=e,\mu@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to ZZ\to 4l)@f$
 * in the current model and in the Standard Model.
 */
class BrHtoZZ4lRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    BrHtoZZ4lRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to ZZ\to 4l)@f$ with @f$l=e,\mu@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to ZZ\to 4l)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtoZZ4eRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to ZZ\to 4e)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to ZZ\to 4e)@f$
 * in the current model and in the Standard Model.
 */
class BrHtoZZ4eRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    BrHtoZZ4eRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to ZZ\to 4e)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to ZZ\to 4e)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtoZZ2e2muRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to ZZ\to 2e 2\mu)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to ZZ\to 2e 2\mu)@f$
 * in the current model and in the Standard Model.
 */
class BrHtoZZ2e2muRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    BrHtoZZ2e2muRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to ZZ\to 2e 2\mu)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to ZZ\to 2e 2\mu)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtoZZ4muRatio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to ZZ\to 4\mu)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to ZZ\to 4\mu)@f$
 * in the current model and in the Standard Model.
 */
class BrHtoZZ4muRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    BrHtoZZ4muRatio(const StandardModel& SM_i);

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to ZZ\to 4\mu)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to ZZ\to 4\mu)@f$
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

/**
 * @class BrHtogaga_over_mumu_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to \gamma\gamma)@f$@f/@f$Br@f$(H\to \mu\mu)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to \gamma\gamma)@f$@f/@f$Br@f$(H\to \mu\mu)@f$
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
     * @brief A method to compute the the ratio of the Br@f$(H\to \gamma\gamma)@f$@f/@f$Br@f$(H\to \mu\mu)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \gamma\gamma)@f$@f/@f$Br@f$(H\to \mu\mu)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtoZga_over_mumu_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to Z\gamma)@f$@f/@f$Br@f$(H\to \mu\mu)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to Z\gamma)@f$@f/@f$Br@f$(H\to \mu\mu)@f$
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
     * @brief A method to compute the the ratio of the Br@f$(H\to Z\gamma)@f$@f/@f$Br@f$(H\to \mu\mu)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to Z\gamma)@f$@f/@f$Br@f$(H\to \mu\mu)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtoZmumuga_over_mumu_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to Z\gamma\to\mu\mu \gamma)@f$@f/@f$Br@f$(H\to \mu\mu)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to Z\gamma\to\mu\mu \gamma)@f$@f/@f$Br@f$(H\to \mu\mu)@f$
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
     * @brief A method to compute the the ratio of the Br@f$(H\to Z\gamma\to\mu\mu \gamma)@f$@f/@f$Br@f$(H\to \mu\mu)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to Z\gamma\to\mu\mu \gamma)@f$@f/@f$Br@f$(H\to \mu\mu)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtogaga_over_4l_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to \gamma\gamma)@f$@f/@f$Br@f$(H\to 4\ell)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to \gamma\gamma)@f$@f/@f$Br@f$(H\to 4\ell)@f$
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
     * @brief A method to compute the the ratio of the Br@f$(H\to \gamma\gamma)@f$@f/@f$Br@f$(H\to 4\ell)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \gamma\gamma)@f$@f/@f$Br@f$(H\to 4\ell)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};



/**
 * @class BrHtobb_over_4l_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to bb)@f$@f/@f$Br@f$(H\to 4\ell)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to bb)@f$@f/@f$Br@f$(H\to 4\ell)@f$
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
     * @brief A method to compute the the ratio of the Br@f$(H\to bb)@f$@f/@f$Br@f$(H\to 4\ell)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to bb)@f$@f/@f$Br@f$(H\to 4\ell)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};



/**
 * @class BrHto2l2v_over_4l_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to 2\ell 2\nu)@f$@f/@f$Br@f$(H\to 4\ell)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to 2\ell 2\nu)@f$@f/@f$Br@f$(H\to 4\ell)@f$
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
     * @brief A method to compute the the ratio of the Br@f$(H\to 2\ell 2\nu)@f$@f/@f$Br@f$(H\to 4\ell)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to 2\ell 2\nu)@f$@f/@f$Br@f$(H\to 4\ell)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};



/**
 * @class BrHtotautau_over_4l_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to \tau\tau)@f$@f/@f$Br@f$(H\to 4\ell)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to \tau\tau)@f$@f/@f$Br@f$(H\to 4\ell)@f$
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
     * @brief A method to compute the the ratio of the Br@f$(H\to \tau\tau)@f$@f/@f$Br@f$(H\to 4\ell)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \tau\tau)@f$@f/@f$Br@f$(H\to 4\ell)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};



/**
 * @class BrHtogaga_over_2e2mu_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to \gamma\gamma)@f$@f/@f$Br@f$(H\to 2e 2\mu)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to \gamma\gamma)@f$@f/@f$Br@f$(H\to 2e 2\mu)@f$
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
     * @brief A method to compute the the ratio of the Br@f$(H\to \gamma\gamma)@f$@f/@f$Br@f$(H\to 2e 2\mu)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \gamma\gamma)@f$@f/@f$Br@f$(H\to 2e 2\mu)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtoZga_over_4l_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to Z\gamma)@f$@f/@f$Br@f$(H\to 4\ell)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to Z\gamma)@f$@f/@f$Br@f$(H\to 4\ell)@f$
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
     * @brief A method to compute the the ratio of the Br@f$(H\to Z\gamma)@f$@f/@f$Br@f$(H\to 4\ell)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to Z\gamma)@f$@f/@f$Br@f$(H\to 4\ell)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};



/**
 * @class BrHtomumu_over_4l_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to \mu\mu)@f$@f/@f$Br@f$(H\to 4\ell)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to \mu\mu)@f$@f/@f$Br@f$(H\to 4\ell)@f$
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
     * @brief A method to compute the the ratio of the Br@f$(H\to \mu\mu)@f$@f/@f$Br@f$(H\to 4\ell)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \mu\mu)@f$@f/@f$Br@f$(H\to 4\ell)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtomumu_over_4mu_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to \mu\mu)@f$@f/@f$Br@f$(H\to 4\mu)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to \mu\mu)@f$@f/@f$Br@f$(H\to 4\mu)@f$
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
     * @brief A method to compute the the ratio of the Br@f$(H\to \mu\mu)@f$@f/@f$Br@f$(H\to 4\mu)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \mu\mu)@f$@f/@f$Br@f$(H\to 4\mu)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};



/**
 * @class BrHto4l_over_gaga_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to 4\ell)@f$@f/@f$Br@f$(H\to \gamma\gamma)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to 4\ell)@f$@f/@f$Br@f$(H\to \gamma\gamma)@f$
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
     * @brief A method to compute the the ratio of the Br@f$(H\to 4\ell)@f$@f/@f$Br@f$(H\to \gamma\gamma)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to 4\ell)@f$@f/@f$Br@f$(H\to \gamma\gamma)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtoZga_over_gaga_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to Z \gamma)@f$@f/@f$Br@f$(H\to \gamma\gamma)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to Z \gamma)@f$@f/@f$Br@f$(H\to \gamma\gamma)@f$
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
     * @brief A method to compute the the ratio of the Br@f$(H\to Z \gamma)@f$@f/@f$Br@f$(H\to \gamma\gamma)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to Z \gamma)@f$@f/@f$Br@f$(H\to \gamma\gamma)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtomumu_over_gaga_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio of the Br@f$(H\to \mu\mu)@f$@f/@f$Br@f$(H\to \gamma\gamma)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to \mu\mu)@f$@f/@f$Br@f$(H\to \gamma\gamma)@f$
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
     * @brief A method to compute the the ratio of the Br@f$(H\to \mu\mu)@f$@f/@f$Br@f$(H\to \gamma\gamma)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \mu\mu)@f$@f/@f$Br@f$(H\to \gamma\gamma)@f$
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
     * @brief A method to compute the value of @f$\mu_{VBF}\times BR(H \to invisible}@f$ in the current model.
     * @return @f$\mu_{VBF}\times BR(H \to invisible}@f$
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
 * @brief A class for computing the quantity @f$\mu_{pp \to VH}\times BR(H \to invisible}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the quantity @f$\mu_{pp \to VH}\times BR(H \to invisible}@f$, i.e.
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
 * in the @f$H,Z\to b\bar{b}@f$ channelin the boosted region in the current model.
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


//  Full signal strengths at lepton colliders
//  -----------------------------------------

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
     */
    mueeZHbb(const StandardModel& SM_i, const double sqrt_s_i);

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
     */
    mueeZHcc(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to cc}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to cc}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
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
     */
    mueeZHgg(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to gg}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to gg}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
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
     */
    mueeZHWW(const StandardModel& SM_i, const double sqrt_s_i);
    
    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to WW}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to WW}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
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
     */
    mueeZHtautau(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to \tau\tau}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to \tau\tau}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
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
     */
    mueeZHZZ(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to ZZ}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to ZZ}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
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
     */
    mueeZHZga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to Z\gamma}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to Z\gamma}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
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
     */
    mueeZHgaga(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to \gamma\gamma}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to \gamma\gamma}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
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
     */
    mueeZHmumu(const StandardModel& SM_i, const double sqrt_s_i);
    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to \mu\mu}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to \mu\mu}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeZHBRinv
 * @ingroup NewPhysics
 * @brief A class for computing the quantity @f$\mu_{e^+e^- \to ZH}\times BR(H \to invisible}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the quantity @f$\mu_{e^+e^- \to ZH}\times BR(H \to invisible}@f$, i.e.
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
     */
    mueeZHBRinv(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH}\times BR(H \to invisible}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH}\times BR(H \to invisible}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
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
     */
    mueeZHinv(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to invisible}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to invisible}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeZHbbPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to bb}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to bb}@f$ between the 
 * @f$e^+e^- \to ZH, H \to bb@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZHbbPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
     */
    mueeZHbbPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

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
 * @class mueeZHccPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to cc}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to cc}@f$ between the 
 * @f$e^+e^- \to ZH, H \to cc@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZHccPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
     */
    mueeZHccPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

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
 * @class mueeZHggPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to gg}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to gg}@f$ between the 
 * @f$e^+e^- \to ZH, H \to gg@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZHggPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
     */
    mueeZHggPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

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
 * @class mueeZHWWPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to WW}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to WW}@f$ between the 
 * @f$e^+e^- \to ZH, H \to WW@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZHWWPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
     */
    mueeZHWWPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

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
 * @class mueeZHtautauPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to \tau\tau}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to \tau\tau}@f$ between the 
 * @f$e^+e^- \to ZH, H \to \tau\tau@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZHtautauPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
     */
    mueeZHtautauPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

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
 * @class mueeZHZZPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to ZZ}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to ZZ}@f$ between the 
 * @f$e^+e^- \to ZH, H \to ZZ@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZHZZPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
     */
    mueeZHZZPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

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
 * @class mueeZHZgaPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to Z\gamma}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to Z\gamma}@f$ between the 
 * @f$e^+e^- \to ZH, H \to Z\gamma@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZHZgaPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
     */
    mueeZHZgaPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

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
 * @class mueeZHgagaPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to \gamma\gamma}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to \gamma\gamma}@f$ between the 
 * @f$e^+e^- \to ZH, H \to \gamma\gamma@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZHgagaPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
     */
    mueeZHgagaPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

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
 * @class mueeZHmumuPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to \mu\mu}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to \mu\mu}@f$ between the 
 * @f$e^+e^- \to ZH, H \to \mu\mu@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZHmumuPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
     */
    mueeZHmumuPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

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
 * @class mueeZHBRinvPol
 * @ingroup NewPhysics
 * @brief A class for computing the quantity @f$\mu_{e^+e^- \to ZH}\times BR(H \to invisible}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the quantity @f$\mu_{e^+e^- \to ZH}\times BR(H \to invisible}@f$, i.e.
 * the ratio between the @f$e^+e^- \to ZH@f$ 
 * associated production cross-section in the current model and 
 * the SM, times the total invisible branching ratio.
 */
class mueeZHBRinvPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
     */
    mueeZHBRinvPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH}\times BR(H \to invisible}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH}\times BR(H \to invisible}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeZHinvPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to invisible}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to invisible}@f$ between the 
 * @f$e^+e^- \to ZH, H \to invisible@f$ 
 * associated production cross-section in the current model and the SM.
 */
class mueeZHinvPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
     */
    mueeZHinvPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to bb}@f$@f$.
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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to bb}@f$@f$.
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
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to cc}@f$@f$.
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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to gg}@f$@f$.
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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to WW}@f$@f$.
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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to \tau\tau}@f$@f$.
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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to ZZ}@f$@f$.
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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to Z\gamma}@f$@f$.
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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to \gamma\gamma}@f$@f$.
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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to \mu\mu}@f$@f$.
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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to bb}@f$@f$,
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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to bb}@f$@f$,
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
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to cc}@f$@f$,
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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to cc}@f$@f$,
 * excluding contributions from on-shell @f$Z\to \nu\bar{\nu} @f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to \nu\nu H, H \to cc}@f$ between the 
 * @f$ e^{+}e^{-}\to \nu\bar{\nu} H, H \to bb @f$ production
 * cross-section in the current model and in the Standard Model.
 */
class mueeHvvccPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
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
 * @class mueeHvvgg
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to gg}@f$@f$,
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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to gg}@f$@f$,
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
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to WW}@f$@f$,
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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to WW}@f$@f$,
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
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to \tau\tau}@f$@f$,
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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to \tau\tau}@f$@f$,
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
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to ZZ}@f$@f$,
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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to ZZ}@f$@f$,
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
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to Z\gamma}@f$@f$,
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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to Z\gamma}@f$@f$,
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
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to \gamma\gamma}@f$@f$,
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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to \gamma\gamma}@f$@f$,
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
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to \mu\mu}@f$@f$,
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
 * @brief A class for computing the ratio @f$@f$\mu_{e^+e^- \to \nu\nu H, H \to \mu\mu}@f$@f$,
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
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
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
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
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
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
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
 * @brief A class for computing the quantity @f$\mu_{VBF}\times BR(H \to invisible}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the quantity @f$\mu_{VBF}\times BR(H \to invisible}@f$, i.e.
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
 * @brief A class for computing the quantity @f$\mu_{pp \to VH}\times BR(H \to invisible}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the quantity @f$\mu_{pp \to VH}\times BR(H \to invisible}@f$, i.e.
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



/**
 * @class BrHto2l2v_over_gaga_Ratio
 * @ingroup NewPhysics
 * @brief A class for computing the ratio Br@f$(H\to 2\ell 2\nu)@f$@f/@f$Br@f$(H\to\gamma\gamma)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio Br@f$(H\to 2\ell 2\nu)@f$@f/@f$Br@f$(H\to\gamma\gamma)@f$
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
     * @brief A method to compute the the ratio Br@f$(H\to 2\ell 2\nu)@f$@f/@f$Br@f$(H\to\gamma\gamma)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to 2\ell 2\nu)@f$@f/@f$Br@f$(H\to\gamma\gamma)@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
};



#endif	/* HIGGSTHOBSERVABLES_H */

