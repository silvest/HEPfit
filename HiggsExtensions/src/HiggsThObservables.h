/*
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HIGGSTHOBSERVABLES_H
#define	HIGGSTHOBSERVABLES_H

#include "ThObservable.h"
#include "NPbase.h"

/**
 * @class muggH
 * @ingroup HiggsExtensions
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
    muggH(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muggH called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{ggH}@f$ in the current model.
     * @return @f$\mu_{ggH}@f$
     */
    double computeThValue()
    {
        return myNPbase->muggH(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVBF
 * @ingroup HiggsExtensions
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
    muVBF(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muVBF called with a class whose parent is not NPbase");

    }

    /**
     * @brief A method to compute the value of @f$\mu_{VBF}@f$ in the current model.
     * @return @f$\mu_{VBF}@f$
     */
    double computeThValue()
    {
        return myNPbase->muVBF(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVBFgamma
 * @ingroup HiggsExtensions
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
    muVBFgamma(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muVBFgamma called with a class whose parent is not NPbase");

    }

    /**
     * @brief A method to compute the value of @f$\mu_{VBF+\gamma}@f$ in the current model.
     * @return @f$\mu_{VBF+\gamma}@f$
     */
    double computeThValue()
    {
        return myNPbase->muVBFgamma(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class mueeWBF
 * @ingroup HiggsExtensions
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
    mueeWBF(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeWBF called with a class whose parent is not NPbase");

    }

    /**
     * @brief A method to compute the value of @f$\mu_{eeWBF}@f$ in the current model.
     * @return @f$\mu_{eeWBF}@f$
     */
    double computeThValue()
    {
        return myNPbase->mueeWBF(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeZBF
 * @ingroup HiggsExtensions
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
    mueeZBF(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeZBF called with a class whose parent is not NPbase");

    }

    /**
     * @brief A method to compute the value of @f$\mu_{e^{+}e^{-}\to e^{+}e^{-} H}@f$ in the current model.
     * @return @f$\mu_{e^{+}e^{-}\to e^{+}e^{-} H}@f$
     */
    double computeThValue()
    {
        return myNPbase->mueeZBF(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeZBFPol
 * @ingroup HiggsExtensions
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
    mueeZBFPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeZBFPol called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{e^{+}e^{-}\to e^{+}e^{-} H}@f$ in the current model.
     * @return @f$\mu_{e^{+}e^{-}\to e^{+}e^{-} H}@f$
     */
    double computeThValue()
    {
        return myNPbase->mueeZBFPol(sqrt_s,Pol_em, Pol_ep);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};

/**
 * @class muepWBF
 * @ingroup HiggsExtensions
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
    muepWBF(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muepWBF called with a class whose parent is not NPbase");

    }

    /**
     * @brief A method to compute the value of @f$\mu_{epWBF}@f$ in the current model.
     * @return @f$\mu_{eeWBF}@f$
     */
    double computeThValue()
    {
        return myNPbase->muepWBF(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muepZBF
 * @ingroup HiggsExtensions
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
    muepZBF(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muepZBF called with a class whose parent is not NPbase");

    }

    /**
     * @brief A method to compute the value of @f$\mu_{epZBF}@f$ in the current model.
     * @return @f$\mu_{eeZBF}@f$
     */
    double computeThValue()
    {
        return myNPbase->muepZBF(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muWH
 * @ingroup HiggsExtensions
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
    muWH(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muWH called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{WH}@f$ in the current model.
     * @return @f$\mu_{WH}@f$
     */
    double computeThValue()
    {
        return myNPbase->muWH(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muZH
 * @ingroup HiggsExtensions
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
    muZH(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muZH called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{ZH}@f$ in the current model.
     * @return @f$\mu_{ZH}@f$
     */
    double computeThValue()
    {
        return myNPbase->muZH(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class mueeZH
 * @ingroup HiggsExtensions
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
    mueeZH(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeZH called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH}@f$
     */
    double computeThValue()
    {
        return myNPbase->mueeZH(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class mueeZHPol
 * @ingroup HiggsExtensions
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
    mueeZHPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeZH called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH}@f$
     */
    double computeThValue()
    {
        return myNPbase->mueeZHPol(sqrt_s,Pol_em, Pol_ep);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};

/**
 * @class muVH
 * @ingroup HiggsExtensions
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
    muVH(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muVH called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{VH}@f$ in the current model.
     * @return @f$\mu_{VH}@f$
     */
    double computeThValue()
    {
        return myNPbase->muVH(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVBFpVH
 * @ingroup HiggsExtensions
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
    muVBFpVH(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muVBFpVH called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{VBF+VH}@f$ in the current model.
     * @return @f$\mu_{VBF+VH}@f$
     */
    double computeThValue()
    {
        return myNPbase->muVBFpVH(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muttH
 * @ingroup HiggsExtensions
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
    muttH(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muttH called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{ttH}@f$ in the current model. 
     * @return @f$\mu_{ttH}@f$
     */
    double computeThValue()
    {
        return myNPbase->muttH(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muggHpttH
 * @ingroup HiggsExtensions
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
    muggHpttH(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muggHpttH called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{ggH+ttH}@f$ in the current model. 
     * @return @f$\mu_{ggH+ttH}@f$
     */
    double computeThValue()
    {
        return myNPbase->muggHpttH(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class mueettH
 * @ingroup HiggsExtensions
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
    mueettH(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueettH called with a class whose parent is not NPbase");

    }

    /**
     * @brief A method to compute the value of @f$\mu_{eettH}@f$ in the current model.
     * @return @f$\mu_{eettH}@f$
     */
    double computeThValue()
    {
        return myNPbase->mueettH(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mummH
 * @ingroup HiggsExtensions
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
    mummH(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mummH called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{\mu\mu H}@f$ in the current model.
     * @return @f$\mu_{\mu\mu H}@f$
     */
    double computeThValue()
    {
        return myNPbase->mummH(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class GammaHRatio
 * @ingroup HiggsExtensions
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
    GammaHRatio(const StandardModel& SM_i)
    : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("GammaHRatio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Higgs width
     * in the current model and in the Standard Model.
     * @return @f$\Gamma_{H}/\Gamma_{H}^{SM}@f$
     */
    double computeThValue()
    {
        return myNPbase->computeGammaTotalRatio();
    }

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHinvisible
 * @ingroup HiggsExtensions
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
    BrHinvisible(const StandardModel& SM_i)
    : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHinvisible called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the branching ratio of Higgs decays into
     * invisible partciles.
     * @return Br@f$(H\to invisible)@f$
     */
    double computeThValue()
    {
        return myNPbase->Br_H_inv();
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHexotic
 * @ingroup HiggsExtensions
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
    BrHexotic(const StandardModel& SM_i)
    : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHexotic called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the branching ratio of Higgs decays into
     * exotics (invisible or not).
     * @return Br@f$(H\to exotic)@f$
     */
    double computeThValue()
    {
        return myNPbase->Br_H_exo();
    }

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtovisRatio
 * @ingroup HiggsExtensions
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
    BrHtovisRatio(const StandardModel& SM_i)
    : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtovisRatio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to visible)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to visible)@f$
     */
    double computeThValue()
    {
        return myNPbase->BrHvisRatio();
    }

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtoggRatio
 * @ingroup HiggsExtensions
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
    BrHtoggRatio(const StandardModel& SM_i)
    : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtoggRatio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to gg)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to gg)@f$
     */
    double computeThValue()
    {
        return myNPbase->BrHggRatio();
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtoWWRatio
 * @ingroup HiggsExtensions
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
    BrHtoWWRatio(const StandardModel& SM_i)
    : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtoWWRatio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to WW\to 4f)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to WW\to 4f)@f$
     */
    double computeThValue()
    {
        return myNPbase->BrHWWRatio();
    }

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtoWW2l2vRatio
 * @ingroup HiggsExtensions
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
    BrHtoWW2l2vRatio(const StandardModel& SM_i)
    : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtoWW2l2vRatio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to WW\to 2l2\nu)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to WW\to 2l2\nu)@f$
     */
    double computeThValue()
    {
        return myNPbase->BrHWW2l2vRatio();
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtoZZRatio
 * @ingroup HiggsExtensions
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
    BrHtoZZRatio(const StandardModel& SM_i)
    : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtoZZRatio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to ZZ\to 4f)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to ZZ\to 4f)@f$
     */
    double computeThValue()
    {
        return myNPbase->BrHZZRatio();
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtoZZ4lRatio
 * @ingroup HiggsExtensions
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
    BrHtoZZ4lRatio(const StandardModel& SM_i)
    : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtoZZ4lRatio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to ZZ\to 4l)@f$ with @f$l=e,\mu@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to ZZ\to 4l)@f$
     */
    double computeThValue()
    {
        return myNPbase->BrHZZ4lRatio();
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtoZgaRatio
 * @ingroup HiggsExtensions
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
    BrHtoZgaRatio(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtoZgaRatio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to Z\gamma)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to Z\gamma)@f$
     */
    double computeThValue()
    {
        return myNPbase->BrHZgaRatio();
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtoZgallRatio
 * @ingroup HiggsExtensions
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
    BrHtoZgallRatio(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtoZgallRatio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to Z\gamma\to ll\gamma)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to Z\gamma\to ll\gamma)@f$
     */
    double computeThValue()
    {
        return myNPbase->BrHZgallRatio();
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtogagaRatio
 * @ingroup HiggsExtensions
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
    BrHtogagaRatio(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtogagaRatio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to \gamma\gamma)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \gamma\gamma)@f$
     */
    double computeThValue()
    {
        return myNPbase->BrHgagaRatio();
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtomumuRatio
 * @ingroup HiggsExtensions
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
    BrHtomumuRatio(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtomumuRatio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to \mu^+\mu^-)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \mu^+\mu^-)@f$
     */
    double computeThValue()
    {
        return myNPbase->BrHmumuRatio();
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtotautauRatio
 * @ingroup HiggsExtensions
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
    BrHtotautauRatio(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtotautauRatio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to \tau^+\tau^-)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \tau^+\tau^-)@f$
     */
    double computeThValue()
    {
        return myNPbase->BrHtautauRatio();
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtoccRatio
 * @ingroup HiggsExtensions
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
    BrHtoccRatio(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtoccRatio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to c\bar{c})@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to c\bar{c})@f$
     */
    double computeThValue()
    {
        return myNPbase->BrHccRatio();
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtobbRatio
 * @ingroup HiggsExtensions
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
    BrHtobbRatio(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtobbRatio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to b\bar{b})@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to b\bar{b})@f$
     */
    double computeThValue()
    {
        return myNPbase->BrHbbRatio();
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtogaga_over_mumu_Ratio
 * @ingroup HiggsExtensions
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
    BrHtogaga_over_mumu_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtogaga_over_mumu_Ratio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to \gamma\gamma)@f$@f/@f$Br@f$(H\to \mu\mu)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \gamma\gamma)@f$@f/@f$Br@f$(H\to \mu\mu)@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return (1.0 + (myNPbase->BrHgagaRatio()) - (myNPbase->BrHmumuRatio()) );
        } else {
            return (myNPbase->BrHgagaRatio())/(myNPbase->BrHmumuRatio());
        }
    }

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtogaga_over_4l_Ratio
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio of the Br@f$(H\to \gamma\gamma)@f$@f/@f$Br@f$(H\to 4\ell)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to \gamma\gamma)@f$@f/@f$Br@f$(H\to 4\ell)@f$
 * in the current model and in the Standard Model (neglects new physics in Z decays).
 */
class BrHtogaga_over_4l_Ratio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtogaga_over_4l_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtogaga_over_4l_Ratio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to \gamma\gamma)@f$@f/@f$Br@f$(H\to 4\ell)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \gamma\gamma)@f$@f/@f$Br@f$(H\to 4\ell)@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return (1.0 + (myNPbase->BrHgagaRatio()) - (myNPbase->BrHZZ4lRatio()) );
        } else {
            return (myNPbase->BrHgagaRatio())/(myNPbase->BrHZZ4lRatio());
        }
    }

private:
    const NPbase* myNPbase;
};



/**
 * @class BrHtoZga_over_4l_Ratio
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio of the Br@f$(H\to Z\gamma)@f$@f/@f$Br@f$(H\to 4\ell)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to Z\gamma)@f$@f/@f$Br@f$(H\to 4\ell)@f$
 * in the current model and in the Standard Model (neglects new physics in Z decays).
 */
class BrHtoZga_over_4l_Ratio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtoZga_over_4l_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtoZga_over_4l_Ratio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to Z\gamma)@f$@f/@f$Br@f$(H\to 4\ell)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to Z\gamma)@f$@f/@f$Br@f$(H\to 4\ell)@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return (1.0 + (myNPbase->BrHZgaRatio()) - (myNPbase->BrHZZ4lRatio()) );
        } else {
            return (myNPbase->BrHZgaRatio())/(myNPbase->BrHZZ4lRatio());
        }
    }

private:
    const NPbase* myNPbase;
};



/**
 * @class BrHtomumu_over_4l_Ratio
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio of the Br@f$(H\to \mu\mu)@f$@f/@f$Br@f$(H\to 4\ell)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to \mu\mu)@f$@f/@f$Br@f$(H\to 4\ell)@f$
 * in the current model and in the Standard Model (neglects new physics in Z decays).
 */
class BrHtomumu_over_4l_Ratio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtomumu_over_4l_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtogaga_over_4l_Ratio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to \mu\mu)@f$@f/@f$Br@f$(H\to 4\ell)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \mu\mu)@f$@f/@f$Br@f$(H\to 4\ell)@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return (1.0 + (myNPbase->BrHmumuRatio()) - (myNPbase->BrHZZ4lRatio()) );
        } else {
            return (myNPbase->BrHmumuRatio())/(myNPbase->BrHZZ4lRatio());
        }
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHto4l_over_gaga_Ratio
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio of the Br@f$(H\to 4\ell)@f$@f/@f$Br@f$(H\to \gamma\gamma)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to 4\ell)@f$@f/@f$Br@f$(H\to \gamma\gamma)@f$
 * in the current model and in the Standard Model (neglects new physics in Z decays).
 */
class BrHto4l_over_gaga_Ratio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHto4l_over_gaga_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHto4l_over_gaga_Ratio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to 4\ell)@f$@f/@f$Br@f$(H\to \gamma\gamma)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to 4\ell)@f$@f/@f$Br@f$(H\to \gamma\gamma)@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return (1.0 + (myNPbase->BrHZZ4lRatio()) - (myNPbase->BrHgagaRatio()) );
        } else {
            return (myNPbase->BrHZZ4lRatio())/(myNPbase->BrHgagaRatio());
        }
    }

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtoZga_over_gaga_Ratio
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio of the Br@f$(H\to Z \gamma)@f$@f/@f$Br@f$(H\to \gamma\gamma)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to Z \gamma)@f$@f/@f$Br@f$(H\to \gamma\gamma)@f$
 * in the current model and in the Standard Model (neglects new physics in Z decays).
 */
class BrHtoZga_over_gaga_Ratio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtoZga_over_gaga_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtoZga_over_gaga_Ratio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to Z \gamma)@f$@f/@f$Br@f$(H\to \gamma\gamma)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to Z \gamma)@f$@f/@f$Br@f$(H\to \gamma\gamma)@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return (1.0 + (myNPbase->BrHZgaRatio()) - (myNPbase->BrHgagaRatio()) );
        } else {
            return (myNPbase->BrHZgaRatio())/(myNPbase->BrHgagaRatio());
        }
    }

private:
    const NPbase* myNPbase;
};


/**
 * @class BrHtomumu_over_gaga_Ratio
 * @ingroup HiggsExtensions
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
    BrHtomumu_over_gaga_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtomumu_over_gaga_Ratio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to \mu\mu)@f$@f/@f$Br@f$(H\to \gamma\gamma)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \mu\mu)@f$@f/@f$Br@f$(H\to \gamma\gamma)@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return (1.0 + (myNPbase->BrHmumuRatio()) - (myNPbase->BrHgagaRatio()) );
        } else {
            return (myNPbase->BrHmumuRatio())/(myNPbase->BrHgagaRatio());
        }
    }

private:
    const NPbase* myNPbase;
};


/**
 * @class muggHgaga
 * @ingroup HiggsExtensions
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
    muggHgaga(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muggHgaga called with a class whose parent is not NPbase");
    }

    /**
     * @brief 
     * @return 
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ( -1.0 + (myNPbase->muggH(sqrt_s)) + (myNPbase->BrHgagaRatio()) );
        } else {
            return myNPbase->muggHgaga(sqrt_s);
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muggHgagaInt
 * @ingroup HiggsExtensions
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
    muggHgagaInt(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muggHgagaInt called with a class whose parent is not NPbase");
    }

    /**
     * @brief 
     * @return 
     */
    double computeThValue()
    {
        return (myNPbase->muggHgagaInt(sqrt_s))/(myNPbase->BrHgagaRatio());
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muVBFHgaga
 * @ingroup HiggsExtensions
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
    muVBFHgaga(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muVBFHgaga called with a class whose parent is not NPbase");
    }

    /**
     * @brief 
     * @return 
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ( -1.0 + (myNPbase->muVBF(sqrt_s)) + (myNPbase->BrHgagaRatio()) );
        } else {
            return myNPbase->muVBFHgaga(sqrt_s);
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVHgaga
 * @ingroup HiggsExtensions
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
    muVHgaga(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muVHgaga called with a class whose parent is not NPbase");
    }

    /**
     * @brief 
     * @return 
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ( -1.0 + (myNPbase->muVH(sqrt_s)) + (myNPbase->BrHgagaRatio()) );
        } else {
            return myNPbase->muVHgaga(sqrt_s);
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muttHgaga
 * @ingroup HiggsExtensions
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
    muttHgaga(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muttHgaga called with a class whose parent is not NPbase");
    }

    /**
     * @brief 
     * @return 
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ( -1.0 + (myNPbase->muttH(sqrt_s)) + (myNPbase->BrHgagaRatio()) );
        } else {
            return myNPbase->muttHgaga(sqrt_s);
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muggHZZ
 * @ingroup HiggsExtensions
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
    muggHZZ(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muggHZZ called with a class whose parent is not NPbase");
    }

    /**
     * @brief 
     * @return 
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ( -1.0 + (myNPbase->muggH(sqrt_s)) + (myNPbase->BrHZZRatio()) );
        } else {
            return myNPbase->muggHZZ(sqrt_s);
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVBFHZZ
 * @ingroup HiggsExtensions
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
    muVBFHZZ(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muVBFHZZ called with a class whose parent is not NPbase");
    }

    /**
     * @brief 
     * @return 
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ( -1.0 + (myNPbase->muVBF(sqrt_s)) + (myNPbase->BrHZZRatio()) );
        } else {
            return myNPbase->muVBFHZZ(sqrt_s);
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVHZZ
 * @ingroup HiggsExtensions
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
    muVHZZ(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muVHZZ called with a class whose parent is not NPbase");
    }

    /**
     * @brief 
     * @return 
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ( -1.0 + (myNPbase->muVH(sqrt_s)) + (myNPbase->BrHZZRatio()) );
        } else {
            return myNPbase->muVHZZ(sqrt_s);
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muttHZZ
 * @ingroup HiggsExtensions
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
    muttHZZ(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muttHZZ called with a class whose parent is not NPbase");
    }

    /**
     * @brief 
     * @return 
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ( -1.0 + (myNPbase->muttH(sqrt_s)) + (myNPbase->BrHZZRatio()) );
        } else {
            return myNPbase->muttHZZ(sqrt_s);
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muggHWW
 * @ingroup HiggsExtensions
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
    muggHWW(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muggHWW called with a class whose parent is not NPbase");
    }

    /**
     * @brief 
     * @return 
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ( -1.0 + (myNPbase->muggH(sqrt_s)) + (myNPbase->BrHWWRatio()) );
        } else {
            return myNPbase->muggHWW(sqrt_s);
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVBFHWW
 * @ingroup HiggsExtensions
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
    muVBFHWW(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muVBFHWW called with a class whose parent is not NPbase");
    }

    /**
     * @brief 
     * @return 
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ( -1.0 + (myNPbase->muVBF(sqrt_s)) + (myNPbase->BrHWWRatio()) );
        } else {
            return myNPbase->muVBFHWW(sqrt_s);
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVHWW
 * @ingroup HiggsExtensions
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
    muVHWW(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muVHWW called with a class whose parent is not NPbase");
    }

    /**
     * @brief 
     * @return 
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ( -1.0 + (myNPbase->muVH(sqrt_s)) + (myNPbase->BrHWWRatio()) );
        } else {
            return myNPbase->muVHWW(sqrt_s);
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muttHWW
 * @ingroup HiggsExtensions
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
    muttHWW(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muttHWW called with a class whose parent is not NPbase");
    }

    /**
     * @brief 
     * @return 
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ( -1.0 + (myNPbase->muttH(sqrt_s)) + (myNPbase->BrHWWRatio()) );
        } else {
            return myNPbase->muttHWW(sqrt_s);
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muggHtautau
 * @ingroup HiggsExtensions
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
    muggHtautau(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muggHtautau called with a class whose parent is not NPbase");
    }

    /**
     * @brief 
     * @return 
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ( -1.0 + (myNPbase->muggH(sqrt_s)) + (myNPbase->BrHtautauRatio()) );
        } else {
            return myNPbase->muggHtautau(sqrt_s);
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVBFHtautau
 * @ingroup HiggsExtensions
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
    muVBFHtautau(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muVBFHtautau called with a class whose parent is not NPbase");
    }

    /**
     * @brief 
     * @return 
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ( -1.0 + (myNPbase->muVBF(sqrt_s)) + (myNPbase->BrHtautauRatio()) );
        } else {
            return myNPbase->muVBFHtautau(sqrt_s);
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVHtautau
 * @ingroup HiggsExtensions
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
    muVHtautau(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muVHtautau called with a class whose parent is not NPbase");
    }

    /**
     * @brief 
     * @return 
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ( -1.0 + (myNPbase->muVH(sqrt_s)) + (myNPbase->BrHtautauRatio()) );
        } else {
            return myNPbase->muVHtautau(sqrt_s);
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muttHtautau
 * @ingroup HiggsExtensions
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
    muttHtautau(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muttHtautau called with a class whose parent is not NPbase");
    }

    /**
     * @brief 
     * @return 
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ( -1.0 + (myNPbase->muttH(sqrt_s)) + (myNPbase->BrHtautauRatio()) );
        } else {
            return myNPbase->muttHtautau(sqrt_s);
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muggHbb
 * @ingroup HiggsExtensions
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
    muggHbb(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muggHbb called with a class whose parent is not NPbase");
    }

    /**
     * @brief 
     * @return 
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ( -1.0 + (myNPbase->muggH(sqrt_s)) + (myNPbase->BrHbbRatio()) );
        } else {
            return myNPbase->muggHbb(sqrt_s);
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVBFHbb
 * @ingroup HiggsExtensions
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
    muVBFHbb(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muVBFHbb called with a class whose parent is not NPbase");
    }

    /**
     * @brief 
     * @return 
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ( -1.0 + (myNPbase->muVBF(sqrt_s)) + (myNPbase->BrHbbRatio()) );
        } else {
            return myNPbase->muVBFHbb(sqrt_s);
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVHbb
 * @ingroup HiggsExtensions
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
    muVHbb(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muVHbb called with a class whose parent is not NPbase");
    }

    /**
     * @brief 
     * @return 
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ( -1.0 + (myNPbase->muVH(sqrt_s)) + (myNPbase->BrHbbRatio()) );
        } else {
            return myNPbase->muVHbb(sqrt_s);
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muttHbb
 * @ingroup HiggsExtensions
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
    muttHbb(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muttHbb called with a class whose parent is not NPbase");
    }

    /**
     * @brief 
     * @return 
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ( -1.0 + (myNPbase->muttH(sqrt_s)) + (myNPbase->BrHbbRatio()) );
        } else {
            return myNPbase->muttHbb(sqrt_s);
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muppHmumu
 * @ingroup HiggsExtensions
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
    muppHmumu(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muppHmumu called with a class whose parent is not NPbase");
    }

    /**
     * @brief 
     * @return 
     */
    double computeThValue()
    {
        return myNPbase->muppHmumu(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muppHZga
 * @ingroup HiggsExtensions
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
    muppHZga(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muppHZga called with a class whose parent is not NPbase");
    }

    /**
     * @brief 
     * @return 
     */
    double computeThValue()
    {
        return myNPbase->muppHZga(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class muggHH2ga2b
 * @ingroup HiggsExtensions
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
    muggHH2ga2b(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muggH called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{ggH}@f$ in the current model.
     * @return @f$\mu_{ggH}@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return (-2.0 + (myNPbase->muggHH(sqrt_s)) + (myNPbase->BrHgagaRatio()) + (myNPbase->BrHbbRatio()) );
        } else {
            return (myNPbase->muggHH(sqrt_s))*(myNPbase->BrHgagaRatio())*(myNPbase->BrHbbRatio());
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muttHZbbboost
 * @ingroup HiggsExtensions
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
    muttHZbbboost(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muttHZbbboost called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\sigma(ttH)/\sigma(ttZ)@f$ 
     * in the @f$H,Z\to b\bar{b}@f$ channel in the current model.
     * @return @f$\sigma(ttH)/\sigma(ttZ)@f$ normalized to the SM.
     */
    double computeThValue()
    {
        return (myNPbase->muttHZbbboost(sqrt_s));
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class UpperLimit_ppHZgammaA
 * @ingroup HiggsExtensions
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
    UpperLimit_ppHZgammaA(const StandardModel& SM_i, const double sqrt_s_i) : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("UpperLimit_ppHZgammaA called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute 
     * @return 
     */
    double computeThValue()
    {
        return myNPbase->UpperLimitZgammaA(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class UpperLimit_ppHZgammaA13
 * @ingroup HiggsExtensions
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
    UpperLimit_ppHZgammaA13(const StandardModel& SM_i, const double sqrt_s_i) : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("UpperLimit_ppHZgammaA13 called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute 
     * @return 
     */
    double computeThValue()
    {
        return myNPbase->UpperLimitZgammaA13(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class UpperLimit_ppHZgammaC13
 * @ingroup HiggsExtensions
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
    UpperLimit_ppHZgammaC13(const StandardModel& SM_i, const double sqrt_s_i) : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("UpperLimit_ppHZgammaC13 called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute 
     * @return 
     */
    double computeThValue()
    {
        return myNPbase->UpperLimitZgammaC13(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class UpperLimit_ppHZgammaC
 * @ingroup HiggsExtensions
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
    UpperLimit_ppHZgammaC(const StandardModel& SM_i, const double sqrt_s_i) : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("UpperLimit_ppHZgammaC called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute 
     * @return 
     */
    double computeThValue()
    {
        return myNPbase->UpperLimitZgammaC(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class cg_plus_ct
 * @ingroup HiggsExtensions
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
    cg_plus_ct(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("cg_plus_ct called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute 
     * @return 
     */
    double computeThValue()
    {
        return myNPbase->cgplusct();
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class cga_plus_ct
 * @ingroup HiggsExtensions
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
    cga_plus_ct(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("cga_plus_ct called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute 
     * @return 
     */
    double computeThValue()
    {
        return myNPbase->cgaplusct();
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class cg_minus_cga
 * @ingroup HiggsExtensions
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
    cg_minus_cga(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("cg_minus_cga called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute 
     * @return 
     */
    double computeThValue()
    {
        return myNPbase->cgminuscga();
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class cV_plus_cb
 * @ingroup HiggsExtensions
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
    cV_plus_cb(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("cV_plus_cb called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute 
     * @return 
     */
    double computeThValue()
    {
        return myNPbase->cVpluscb();
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class cV_plus_ctau
 * @ingroup HiggsExtensions
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
    cV_plus_ctau(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("cV_plus_ctau called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute 
     * @return 
     */
    double computeThValue()
    {
        return myNPbase->cVplusctau();
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class cb_minus_cc
 * @ingroup HiggsExtensions
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
    cb_minus_cc(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("cb_minus_cc called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute 
     * @return 
     */
    double computeThValue()
    {
        return myNPbase->cbminuscc();
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class cb_minus_ctau
 * @ingroup HiggsExtensions
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
    cb_minus_ctau(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("cb_minus_ctau called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute 
     * @return 
     */
    double computeThValue()
    {
        return myNPbase->cbminusctau();
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class cc_minus_ctau
 * @ingroup HiggsExtensions
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
    cc_minus_ctau(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("cc_minus_ctau called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute 
     * @return 
     */
    double computeThValue()
    {
        return myNPbase->ccminusctau();
    }

private:
    const NPbase* myNPbase;
};


//  Full signal strengths at lepton colliders
//  -----------------------------------------

/**
 * @class mueeZHbb
 * @ingroup HiggsExtensions
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
    mueeZHbb(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeZHbb called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to bb}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to bb}@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ((myNPbase->mueeZH(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
        } else {
            return (myNPbase->mueeZH(sqrt_s))*(myNPbase->BrHbbRatio());
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeZHcc
 * @ingroup HiggsExtensions
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
    mueeZHcc(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeZHcc called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to cc}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to cc}@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ((myNPbase->mueeZH(sqrt_s)) + (myNPbase->BrHccRatio()) - 1.0);
        } else {
            return (myNPbase->mueeZH(sqrt_s))*(myNPbase->BrHccRatio());
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeZHgg
 * @ingroup HiggsExtensions
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
    mueeZHgg(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeZHgg called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to gg}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to gg}@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ((myNPbase->mueeZH(sqrt_s)) + (myNPbase->BrHggRatio()) - 1.0);
        } else {
            return (myNPbase->mueeZH(sqrt_s))*(myNPbase->BrHggRatio());
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeZHWW
 * @ingroup HiggsExtensions
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
    mueeZHWW(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeZHcc called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to WW}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to WW}@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ((myNPbase->mueeZH(sqrt_s)) + (myNPbase->BrHWWRatio()) - 1.0);
        } else {
            return (myNPbase->mueeZH(sqrt_s))*(myNPbase->BrHWWRatio());
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeZHtautau
 * @ingroup HiggsExtensions
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
    mueeZHtautau(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeZHtautau called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to \tau\tau}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to \tau\tau}@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ((myNPbase->mueeZH(sqrt_s)) + (myNPbase->BrHtautauRatio()) - 1.0);
        } else {
            return (myNPbase->mueeZH(sqrt_s))*(myNPbase->BrHtautauRatio());
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeZHZZ
 * @ingroup HiggsExtensions
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
    mueeZHZZ(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeZHZZ called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to ZZ}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to ZZ}@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ((myNPbase->mueeZH(sqrt_s)) + (myNPbase->BrHZZRatio()) - 1.0);
        } else {
            return (myNPbase->mueeZH(sqrt_s))*(myNPbase->BrHZZRatio());
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeZHgaga
 * @ingroup HiggsExtensions
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
    mueeZHgaga(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeZHgaga called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to \gamma\gamma}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to \gamma\gamma}@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ((myNPbase->mueeZH(sqrt_s)) + (myNPbase->BrHgagaRatio()) - 1.0);
        } else {
            return (myNPbase->mueeZH(sqrt_s))*(myNPbase->BrHgagaRatio());
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeZHmumu
 * @ingroup HiggsExtensions
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
    mueeZHmumu(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeZHmumu called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to \mu\mu}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to \mu\mu}@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ((myNPbase->mueeZH(sqrt_s)) + (myNPbase->BrHmumuRatio()) - 1.0);
        } else {
            return (myNPbase->mueeZH(sqrt_s))*(myNPbase->BrHmumuRatio());
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeZHinv
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to invisible}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to invisible}@f$ between the 
 * @f$e^+e^- \to ZH, H \to invisible@f$ 
 * associated production cross-section in the current model and 
 * the SM @f$e^+e^- \to ZH@f total cross section.
 */
class mueeZHinv : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeZHinv(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeZHinv called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to invisible}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to invisible}@f$
     */
    double computeThValue()
    {
        
        return (myNPbase->mueeZH(sqrt_s))*(myNPbase->Br_H_inv());
            
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeZHbbPol
 * @ingroup HiggsExtensions
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
    mueeZHbbPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeZHbbPol called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to bb}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to bb}@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ((myNPbase->mueeZHPol(sqrt_s,Pol_em, Pol_ep)) + (myNPbase->BrHbbRatio()) - 1.0);
        } else {
            return (myNPbase->mueeZHPol(sqrt_s,Pol_em, Pol_ep))*(myNPbase->BrHbbRatio());
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeZHccPol
 * @ingroup HiggsExtensions
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
    mueeZHccPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeZHccPol called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to cc}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to cc}@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ((myNPbase->mueeZHPol(sqrt_s,Pol_em, Pol_ep)) + (myNPbase->BrHccRatio()) - 1.0);
        } else {
            return (myNPbase->mueeZHPol(sqrt_s,Pol_em, Pol_ep))*(myNPbase->BrHccRatio());
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeZHggPol
 * @ingroup HiggsExtensions
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
    mueeZHggPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeZHggPol called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to gg}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to gg}@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ((myNPbase->mueeZHPol(sqrt_s,Pol_em, Pol_ep)) + (myNPbase->BrHggRatio()) - 1.0);
        } else {
            return (myNPbase->mueeZHPol(sqrt_s,Pol_em, Pol_ep))*(myNPbase->BrHggRatio());
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeZHWWPol
 * @ingroup HiggsExtensions
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
    mueeZHWWPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeZHccPol called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to WW}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to WW}@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ((myNPbase->mueeZHPol(sqrt_s,Pol_em, Pol_ep)) + (myNPbase->BrHWWRatio()) - 1.0);
        } else {
            return (myNPbase->mueeZHPol(sqrt_s,Pol_em, Pol_ep))*(myNPbase->BrHWWRatio());
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeZHtautauPol
 * @ingroup HiggsExtensions
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
    mueeZHtautauPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeZHtautauPol called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to \tau\tau}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to \tau\tau}@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ((myNPbase->mueeZHPol(sqrt_s,Pol_em, Pol_ep)) + (myNPbase->BrHtautauRatio()) - 1.0);
        } else {
            return (myNPbase->mueeZHPol(sqrt_s,Pol_em, Pol_ep))*(myNPbase->BrHtautauRatio());
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeZHZZPol
 * @ingroup HiggsExtensions
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
    mueeZHZZPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeZHZZPol called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to ZZ}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to ZZ}@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ((myNPbase->mueeZHPol(sqrt_s,Pol_em, Pol_ep)) + (myNPbase->BrHZZRatio()) - 1.0);
        } else {
            return (myNPbase->mueeZHPol(sqrt_s,Pol_em, Pol_ep))*(myNPbase->BrHZZRatio());
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeZHgagaPol
 * @ingroup HiggsExtensions
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
    mueeZHgagaPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeZHgagaPol called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to \gamma\gamma}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to \gamma\gamma}@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ((myNPbase->mueeZHPol(sqrt_s,Pol_em, Pol_ep)) + (myNPbase->BrHgagaRatio()) - 1.0);
        } else {
            return (myNPbase->mueeZHPol(sqrt_s,Pol_em, Pol_ep))*(myNPbase->BrHgagaRatio());
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeZHmumuPol
 * @ingroup HiggsExtensions
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
    mueeZHmumuPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeZHmumuPol called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to \mu\mu}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to \mu\mu}@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ((myNPbase->mueeZHPol(sqrt_s,Pol_em, Pol_ep)) + (myNPbase->BrHmumuRatio()) - 1.0);
        } else {
            return (myNPbase->mueeZHPol(sqrt_s,Pol_em, Pol_ep))*(myNPbase->BrHmumuRatio());
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeZHinvPol
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to invisible}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH, H \to invisible}@f$ between the 
 * @f$e^+e^- \to ZH, H \to invisible@f$ 
 * associated production cross-section in the current model and 
 * the SM @f$e^+e^- \to ZH@f total cross section.
 */
class mueeZHinvPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
     */
    mueeZHinvPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeZHinvPol called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH, H \to invisible}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH, H \to invisible}@f$
     */
    double computeThValue()
    {  
        return (myNPbase->mueeZHPol(sqrt_s,Pol_em, Pol_ep))*(myNPbase->Br_H_inv());
            
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeWBFbb
 * @ingroup HiggsExtensions
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
    mueeWBFbb(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeWBFbb called with a class whose parent is not NPbase");

    }

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to bb@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to bb@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ((myNPbase->mueeWBF(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
        } else {
            return (myNPbase->mueeWBF(sqrt_s))*(myNPbase->BrHbbRatio());
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeWBFbbPol
 * @ingroup HiggsExtensions
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
    mueeWBFbbPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeWBFbbPol called with a class whose parent is not NPbase");

    }

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to bb@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to bb@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ((myNPbase->mueeWBFPol(sqrt_s,Pol_em, Pol_ep)) + (myNPbase->BrHbbRatio()) - 1.0);
        } else {
            return (myNPbase->mueeWBFPol(sqrt_s,Pol_em, Pol_ep))*(myNPbase->BrHbbRatio());
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};


/**
 * @class mueeWBFcc
 * @ingroup HiggsExtensions
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
    mueeWBFcc(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeWBFcc called with a class whose parent is not NPbase");

    }

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to cc@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to cc@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ((myNPbase->mueeWBF(sqrt_s)) + (myNPbase->BrHccRatio()) - 1.0);
        } else {
            return (myNPbase->mueeWBF(sqrt_s))*(myNPbase->BrHccRatio());
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeWBFgg
 * @ingroup HiggsExtensions
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
    mueeWBFgg(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeWBFgg called with a class whose parent is not NPbase");

    }

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to gg@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to gg@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ((myNPbase->mueeWBF(sqrt_s)) + (myNPbase->BrHggRatio()) - 1.0);
        } else {
            return (myNPbase->mueeWBF(sqrt_s))*(myNPbase->BrHggRatio());
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeWBFWW
 * @ingroup HiggsExtensions
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
    mueeWBFWW(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeWBFWW called with a class whose parent is not NPbase");

    }

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to WW@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to WW@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ((myNPbase->mueeWBF(sqrt_s)) + (myNPbase->BrHWWRatio()) - 1.0);
        } else {
            return (myNPbase->mueeWBF(sqrt_s))*(myNPbase->BrHWWRatio());
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeWBFtautau
 * @ingroup HiggsExtensions
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
    mueeWBFtautau(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeWBFtautau called with a class whose parent is not NPbase");

    }

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to \tau\tau@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to \tau\tau@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ((myNPbase->mueeWBF(sqrt_s)) + (myNPbase->BrHtautauRatio()) - 1.0);
        } else {
            return (myNPbase->mueeWBF(sqrt_s))*(myNPbase->BrHtautauRatio());
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeWBFZZ
 * @ingroup HiggsExtensions
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
    mueeWBFZZ(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeWBFZZ called with a class whose parent is not NPbase");

    }

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to ZZ@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to ZZ@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ((myNPbase->mueeWBF(sqrt_s)) + (myNPbase->BrHZZRatio()) - 1.0);
        } else {
            return (myNPbase->mueeWBF(sqrt_s))*(myNPbase->BrHZZRatio());
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeWBFgaga
 * @ingroup HiggsExtensions
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
    mueeWBFgaga(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeWBFgaga called with a class whose parent is not NPbase");

    }

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to \gamma\gamma@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to \gamma\gamma@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ((myNPbase->mueeWBF(sqrt_s)) + (myNPbase->BrHgagaRatio()) - 1.0);
        } else {
            return (myNPbase->mueeWBF(sqrt_s))*(myNPbase->BrHgagaRatio());
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class mueeWBFmumu
 * @ingroup HiggsExtensions
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
    mueeWBFmumu(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeWBFmumu called with a class whose parent is not NPbase");

    }

    /**
     * @brief A method to compute the value of @f$e^+e^- \to \nu\nu H, H \to \mu\mu@f$ in the current model.
     * @return @f$e^+e^- \to \nu\nu H, H \to \mu\mu@f$
     */
    double computeThValue()
    {
        if ( (this->getModel()).isModelLinearized() ) {
            return ((myNPbase->mueeWBF(sqrt_s)) + (myNPbase->BrHmumuRatio()) - 1.0);
        } else {
            return (myNPbase->mueeWBF(sqrt_s))*(myNPbase->BrHmumuRatio());
        }
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


#endif	/* HIGGSTHOBSERVABLES_H */

