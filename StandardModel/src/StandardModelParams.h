/*
 * Copyright (C) 2013-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef STANDARDMODELPARAMS_H
#define	STANDARDMODELPARAMS_H

#include <string>
#include <ThObservable.h>
#include <ThObsType.h>

/**
 * @class StandardModelParams
 * @ingroup StandardModel
 * @brief A class for retrieving parameters in QCD and StandardModel.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is responsible for retrieving the values of model
 * parameters in QCD and StandardModel classes. The available parameters are
 * listed in the description of the constructor StandardModelParams().
 */
class StandardModelParams : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] ObsType a reference to an object of type ThObsType
     * @param[in] name_i the name of the parameter to be retrieved:
     * @li "AlsMz":&nbsp; the fine-structure constant @f$\alpha@f$,
     * @li "dAle5Mz":&nbsp; the five-flavour hadronic contribution to the electromagnetic coupling,
     * @f$\Delta\alpha_{\mathrm{had}}^{(5)}(M_Z^2)@f$, 
     * @li "Mz":&nbsp; the Fermi constant @f$G_\mu@f$ in @f$\mathrm{GeV}^{-2}@f$,
     * @li "mtop":&nbsp; the top-quark mass @f$m_t@f$ in GeV,
     * @li "mHl":&nbsp; the Higgs mass @f$m_h@f$ in GeV,
     * @li "delMw":&nbsp; the theoretical uncertainty in @f$M_W@f$ in GeV,
     * @li "delSin2th_l":&nbsp; the theoretical uncertainty in @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$,
     * @li "delGammaZ":&nbsp; the theoretical uncertainty in @f$\Gamma_Z@f$ in GeV.
     */
    StandardModelParams(const ThObsType& ObsType, const std::string name_i)
    : ThObservable(ObsType), name(name_i)
    {
    };

    /**
     * @brief A method to retrieve the value of the desired parameter specified
     * with the constructor. 
     * @return the value of the parameter
     */
    double computeThValue();

private:
    const std::string name;///< The name of the parmeter to be retrieved.

};

#endif	/* STANDARDMODELPARAMS_H */

