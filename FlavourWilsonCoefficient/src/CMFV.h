/* 
 * Copyright (C) 2012 SusyFit Collaboration
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CMFV_H
#define CMFV_H

#include "StandardModel.h"
#include "CMFVMatching.h"

class CMFV : public StandardModel {
public:
    
    static const int NCMFVvars = 1;

    static const std::string CMFVvars[NCMFVvars];
    
    /**
     * @brief CMFV constructor
     */
    CMFV();
    
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);
    
    virtual CMFVMatching& getMatching() const
    {
        return CMFVM.getObj();
    }

    virtual void setMatching(CMFVMatching& FWCMr)
    {
        CMFVM.setObj(FWCMr);
    }
    
    double getFtt() const
    {
        return Ftt;
    }


protected: 
    
    virtual void setParameter(const std::string, const double&);
    mutable Matching<CMFVMatching,CMFV> CMFVM;

private:
    double Ftt;

};

#endif /* CMFV_H */

