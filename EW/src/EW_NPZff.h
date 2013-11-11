/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EW_NPZFF_H
#define	EW_NPZFF_H

#include <StandardModel.h>

/**
 * @class EW_NPZff
 * @ingroup EW
 * @brief A class for ...
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class EW_NPZff {
public:

    EW_NPZff(const StandardModel& SM_i);

    ////////////////////////////////////////////////////////////////////////

    double GammaZ(const double GammaZ_SM) const;
    double sigmaHadron(const double sigmaHadron_SM) const;
    double sin2thetaEff(const double sin2thetaEff_SM) const;
    double PtauPol(const double PtauPol_SM) const;
    double Alepton(const double Alepton_SM) const;
    double Acharm(const double Acharm_SM) const;
    double Abottom(const double Abottom_SM) const;
    double AFBlepton(const double AFBlepton_SM) const;
    double AFBcharm(const double AFBcharm_SM) const;
    double AFBbottom(const double AFBbottom_SM) const;
    double Rlepton(const double Rlepton_SM) const;
    double Rcharm(const double Rcharm_SM) const;
    double Rbottom(const double Rbottom_SM) const;
    
private:
    const StandardModel& SM;
    
};

#endif	/* EW_NPZFF_H */

