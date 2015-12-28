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
 * @f$m_t(m_t)@f$ in the @f$\bar{\mathrm{MS}}@f$ scheme.
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

