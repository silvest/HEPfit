/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MTMSBAR_H
#define MTMSBAR_H

#include "ThObservable.h"

/**
 * @class MtMSbar
 * @ingroup StandardModel
 * @brief A class for @f$m_t(m_t)@f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the value of
 * @f$m_t(m_t)@f$ in the @f$\bar{\mathrm{MS}}@f$ scheme from the top pole mass, using the three-loop
 * (FULLNNNLO) pole @f$\leftrightarrow\bar{\mathrm{MS}}@f$ relation. Note that the internal conversion
 * between the top input mass and the pole mass in QCD::PostUpdate() uses FULLNNLO, so the value differs
 * from the internally used @f$m_t(m_t)@f$ by the @f$O(\alpha_s^3)@f$ term (about -0.5 GeV).
 */
class MtMSbar : public ThObservable {
public:
    /**
     * @brief Constructor. 
     * @param[in] SM_i
     */
    MtMSbar(const StandardModel& SM_i) : ThObservable(SM_i) {};

    /**
     *
     * @return theoretical value of @f$\bar{\mathrm{MS}}@f$
     */
    virtual double computeThValue();

    private:

};

#endif /* MTMSBAR_H */

