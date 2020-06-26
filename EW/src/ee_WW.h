/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EEWW_H
#define EEWW_H

#include <stdexcept>
#include <ThObservable.h>
#include "NPbase.h"


/**
 * @class xseeWW
 * @brief A class for computing the cross section @f$e^+ e^- \to W^+ W^-@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the cross section @f$e^+ e^- \to W^+ W^-@f$.
 */
class xseeWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in GeV
     */
    xseeWW(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("xseeWW called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of the cross section @f$e^+ e^- \to W^+ W^-@f$ in the current model.
     * @return @f$\sigma(e^+ e^- \to W^+ W^-)@f$
     */
    double computeThValue()
    {
        return myNPbase->xseeWW(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class dxseeWWdcosBin
 * @brief A class for computing the integral of the differential cross section 
 * for @f$e^+ e^- \to W^+ W^-@f$ in a given @f$\cos{\theta}@f$ bin.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the integral of the differential cross section
 * for @f$e^+ e^- \to W^+ W^-@f$ in a given @f$\cos{\theta}@f$ bin.
 */
class dxseeWWdcosBin : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in GeV
     */
    dxseeWWdcosBin(const StandardModel& SM_i, const double sqrt_s_i, const double cos1_i, const double cos2_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i), cos1(cos1_i), cos2(cos2_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("dxseeWWdcosBin called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the integral of the differential cross section 
     * for @f$e^+ e^- \to W^+ W^-@f$ in a given @f$\cos{\theta}@f$ bin in the current model.
     * @return @f$\int_{\cos{\theta_1}}^{\cos{\theta_2}} d\sigma(e^+ e^- \to W^+ W^-)/d\cos{\theta}@f$
     */
    double computeThValue()
    {
        return myNPbase->dxseeWWdcosBin(sqrt_s, cos1, cos2);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    const double cos1, cos2;
};


/**
 * @class xseeWWlept
 * @brief A class for computing the cross section @f$e^+ e^- \to W^+ W^- \to \ell \nu \ell \nu@f$
 * in the dim-6 SMEFT, as in arXiv: 1606.06693 [hep-ph].
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the cross section @f$e^+ e^- \to W^+ W^- \to \ell \nu \ell \nu@f$.
 */
class xseeWWlept : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in GeV
     */
    xseeWWlept(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("xseeWWlept called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of the cross section @f$e^+ e^- \to W^+ W^- \to \ell \nu \ell \nu@f$ in the current model.
     * @return @f$\sigma(e^+ e^- \to W^+ W^-\to \ell \nu \ell \nu)@f$
     */
    double computeThValue()
    {
        return myNPbase->xseeWWleptLEP2(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class deltaxseeWWlept
 * @brief A class for computing the NP contribution to the cross section @f$e^+ e^- \to W^+ W^- \to \ell \nu \ell \nu@f$
 * in the dim-6 SMEFT, as in arXiv: 1606.06693 [hep-ph].
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the NP contribution to the cross section @f$e^+ e^- \to W^+ W^- \to \ell \nu \ell \nu@f$.
 */
class deltaxseeWWlept : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in GeV
     */
    deltaxseeWWlept(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("deltaxseeWWlept called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of the NP contribution to the cross section @f$e^+ e^- \to W^+ W^- \to \ell \nu \ell \nu@f$ in the current model.
     * @return @f$\delta \sigma(e^+ e^- \to W^+ W^-\to \ell \nu \ell \nu)@f$
     */
    double computeThValue()
    {
        return myNPbase->deltaxseeWWleptLEP2(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class xseeWWsemil
 * @brief A class for computing the cross section @f$e^+ e^- \to W^+ W^- \to \ell \nu j j @f$
 * in the dim-6 SMEFT, as in arXiv: 1606.06693 [hep-ph].
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the cross section @f$e^+ e^- \to W^+ W^-\to \ell \nu j j @f$.
 */
class xseeWWsemil : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in GeV
     */
    xseeWWsemil(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("xseeWWsemil called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of the cross section @f$e^+ e^- \to W^+ W^-\to \ell \nu j j @f$ in the current model.
     * @return @f$\sigma(e^+ e^- \to W^+ W^-\to \ell \nu j j )@f$
     */
    double computeThValue()
    {
        return myNPbase->xseeWWsemilLEP2(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class deltaxseeWWsemil
 * @brief A class for computing the NP contribution to the cross section @f$e^+ e^- \to W^+ W^- \to \ell \nu j j @f$
 * in the dim-6 SMEFT, as in arXiv: 1606.06693 [hep-ph].
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the NP contribution to the cross section @f$e^+ e^- \to W^+ W^-\to \ell \nu j j @f$.
 */
class deltaxseeWWsemil : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in GeV
     */
    deltaxseeWWsemil(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("deltaxseeWWsemil called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of the NP contribution to the cross section @f$e^+ e^- \to W^+ W^-\to \ell \nu j j @f$ in the current model.
     * @return @f$\delta \sigma(e^+ e^- \to W^+ W^-\to \ell \nu j j )@f$
     */
    double computeThValue()
    {
        return myNPbase->deltaxseeWWsemilLEP2(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class xseeWWhad
 * @brief A class for computing the cross section @f$e^+ e^- \to W^+ W^- \to j j j j@f$
 * in the dim-6 SMEFT, as in arXiv: 1606.06693 [hep-ph].
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the cross section @f$e^+ e^- \to W^+ W^- \to j j j j@f$.
 */
class xseeWWhad : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in GeV
     */
    xseeWWhad(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("xseeWWhad called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of the cross section @f$e^+ e^- \to W^+ W^- \to j j j j@f$ in the current model.
     * @return @f$\sigma(e^+ e^- \to W^+ W^- \to j j j j)@f$
     */
    double computeThValue()
    {
        return myNPbase->xseeWWhadLEP2(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class deltaxseeWWhad
 * @brief A class for computing the NP contribution to the cross section @f$e^+ e^- \to W^+ W^- \to j j j j@f$
 * in the dim-6 SMEFT, as in arXiv: 1606.06693 [hep-ph].
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the NP contribution to the cross section @f$e^+ e^- \to W^+ W^- \to j j j j@f$.
 */
class deltaxseeWWhad : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in GeV
     */
    deltaxseeWWhad(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("deltaxseeWWhad called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of the NP contribution to the cross section @f$e^+ e^- \to W^+ W^- \to j j j j@f$ in the current model.
     * @return @f$\delta \sigma(e^+ e^- \to W^+ W^- \to j j j j)@f$
     */
    double computeThValue()
    {
        return myNPbase->deltaxseeWWhadLEP2(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class xseeWWtot
 * @brief A class for computing the total cross section @f$e^+ e^- \to W^+ W^-@f$
 * in the dim-6 SMEFT, as in arXiv: 1606.06693 [hep-ph].
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the cross section @f$e^+ e^- \to W^+ W^-@f$.
 */
class xseeWWtot : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in GeV
     */
    xseeWWtot(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("xseeWWtot called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of the cross section @f$e^+ e^- \to W^+ W^-@f$ in the current model.
     * @return @f$\sigma(e^+ e^- \to W^+ W^-)@f$
     */
    double computeThValue()
    {
        return myNPbase->xseeWWtotLEP2(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class deltaxseeWWtot
 * @brief A class for computing the NP contribution to the total cross section @f$e^+ e^- \to W^+ W^-@f$
 * in the dim-6 SMEFT, as in arXiv: 1606.06693 [hep-ph].
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the NP contribution to the cross section @f$e^+ e^- \to W^+ W^-@f$.
 */
class deltaxseeWWtot : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in GeV
     */
    deltaxseeWWtot(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("deltaxseeWWtot called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of the NP contribution to the cross section @f$e^+ e^- \to W^+ W^-@f$ in the current model.
     * @return @f$\delta \sigma(e^+ e^- \to W^+ W^-)@f$
     */
    double computeThValue()
    {
        return myNPbase->deltaxseeWWtotLEP2(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class dxseeWWLEP2Bin
 * @brief A class for computing the differential cross section 
 * for @f$e^+ e^- \to W^+ W^- \to l v jj@f$, with @f$ l= e,\mu @f$,
 * for the 4 @f$ cos{\theta}@f$ bins defined in arXiv: 1606.06693 [hep-ph].
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the integral of the differential cross section
 * for @f$e^+ e^- \to W^+ W^- \to l v jj@f$ in a given @f$\cos{\theta}@f$ bin.
 */
class dxseeWWLEP2Bin : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in GeV
     * @param[in] bin_i the bin nnumber: 1-4
     */
    dxseeWWLEP2Bin(const StandardModel& SM_i, const double sqrt_s_i, const int bin_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i), bin(bin_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("dxseeWWLEP2Bin called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the integral of the differential cross section 
     * for @f$e^+ e^- \to W^+ W^-@f$ in a given @f$\cos{\theta}@f$ bin in the current model.
     * @return @f$\int_{bin_i} d\sigma(e^+ e^- \to W^+ W^-)/d\cos{\theta}@f$
     */
    double computeThValue()
    {
        return myNPbase->dxsdcoseeWWlvjjLEP2(sqrt_s, bin);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    const int bin;
};


/**
 * @class deltadxseeWWLEP2Bin
 * @brief A class for computing the NP contribution to the differential cross section 
 * for @f$e^+ e^- \to W^+ W^- \to l v jj@f$, with @f$ l= e,\mu @f$,
 * for the 4 @f$ cos{\theta}@f$ bins defined in arXiv: 1606.06693 [hep-ph].
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the integral of the NP contribution to the differential cross section
 * for @f$e^+ e^- \to W^+ W^- \to l v jj@f$ in a given @f$\cos{\theta}@f$ bin.
 */
class deltadxseeWWLEP2Bin : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in GeV
     * @param[in] bin_i the bin nnumber: 1-4
     */
    deltadxseeWWLEP2Bin(const StandardModel& SM_i, const double sqrt_s_i, const int bin_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i), bin(bin_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("deltadxseeWWLEP2Bin called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the integral of the NP contribution to the differential cross section 
     * for @f$e^+ e^- \to W^+ W^-@f$ in a given @f$\cos{\theta}@f$ bin in the current model.
     * @return @f$\int_{bin_i} \delta d\sigma(e^+ e^- \to W^+ W^-)/d\cos{\theta}@f$
     */
    double computeThValue()
    {
        return myNPbase->deltadxsdcoseeWWlvjjLEP2(sqrt_s, bin);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    const int bin;
};


#endif	/* EEWW_H */

