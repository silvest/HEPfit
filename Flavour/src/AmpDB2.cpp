/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AmpDB2.h"
#include "StandardModel.h"
#include "EvolDF2.h"
#include "HeffDF2.h"

AmpDB2::AmpDB2(const StandardModel& SM_i) 
: mySM(SM_i) 
{
    mySM.initializeBParameter("BBs");
    mySM.initializeBParameter("BBd");
}

gslpp::complex AmpDB2::AmpBd(orders order) 
{
    if (mySM.getFlavour().getHDF2().getCoeffBd().getOrder() < order % 3)
        throw std::runtime_error("DmBd::computeThValue(): requires cofficient of order not computed"); 

    gslpp::vector<gslpp::complex> ** allcoeff = mySM.getFlavour().ComputeCoeffBd( 
            mySM.getBBd().getMu(),
        mySM.getBBd().getScheme());
    
    gslpp::vector<double> me(mySM.getBBd().getBpars());
    double MBd = mySM.getMesons(QCD::B_D).getMass();
    double Mb = mySM.Mrun(mySM.getBBd().getMu(),
                mySM.getQuarks(QCD::BOTTOM).getMass_scale(),
                mySM.getQuarks(QCD::BOTTOM).getMass(), FULLNNLO);
    double Md = mySM.Mrun(mySM.getBBd().getMu(),
                mySM.getQuarks(QCD::DOWN).getMass_scale(),
                mySM.getQuarks(QCD::DOWN).getMass(), FULLNNLO);
    double KBd = MBd/(Mb+Md)*MBd/(Mb+Md);
    double Fb = mySM.getMesons(QCD::B_D).getDecayconst();
    me(0) *= 1./3.*MBd*Fb*Fb;
    me(1) *= -5./24.*KBd*MBd*Fb*Fb;
    me(2) *= 1./24.*KBd*MBd*Fb*Fb;
    me(3) *= 1./4.*KBd*MBd*Fb*Fb;
    me(4) *= 1./12.*KBd*MBd*Fb*Fb;
    
#if SUSYFIT_DEBUG & 1
    std::cout << "Bd: me(0) = " << me(0)  << std::endl;
#endif
#if SUSYFIT_DEBUG & 2
    std::cout << "coefficient Bd: " << (*(allcoeff[LO]) + *(allcoeff[NLO]))(0) << std::endl;
    std::cout << "M: " << me << std::endl;
    std::cout << "mu : " << mySM.getBBd().getMu() << ", mut: " << mySM.getMut() << ", scheme: " << mySM.getBBd().getScheme() << ", B par.: " <<  mySM.getBBd().getBpars()(0) << std::endl;
    std::cout << "U (mut): " << (mySM.getFlavour().getHDF2().getUDF2().Df2Evol(mySM.getBBd().getMu(),mySM.getMut(),LO)(0,0) +  
            mySM.getFlavour().getHDF2().getUDF2().Df2Evol(mySM.getBBd().getMu(),mySM.getMut(),NLO)(0,0))<< std::endl;
#endif

    switch(order) {
        case FULLNLO:
            return((*(allcoeff[LO]) + *(allcoeff[NLO])) * me / HCUT);
        case LO:
            return((*(allcoeff[LO])) * me / HCUT);
        default:
            throw std::runtime_error("AmpDB2::AmpBd(): order not implemented"); 
    }
}

gslpp::complex AmpDB2::AmpBs(orders order) 
{
    if (mySM.getFlavour().getHDF2().getCoeffBs().getOrder() < order % 3)
        throw std::runtime_error("DmBd::computeThValue(): requires cofficient of order not computed"); 

    gslpp::vector<gslpp::complex> ** allcoeff = mySM.getFlavour().ComputeCoeffBs(
            mySM.getBBs().getMu(),
            mySM.getBBs().getScheme());

    gslpp::vector<double> me(mySM.getBBs().getBpars());
    double MBs = mySM.getMesons(QCD::B_S).getMass();
    double Mb = mySM.getQuarks(QCD::BOTTOM).getMass();
    double Ms = mySM.Mrun(mySM.getBBs().getMu(),
                mySM.getQuarks(QCD::STRANGE).getMass_scale(),
                mySM.getQuarks(QCD::STRANGE).getMass(), FULLNNLO);
    double KBs = MBs/(Mb+Ms)*MBs/(Mb+Ms);
    double Fbs = mySM.getMesons(QCD::B_S).getDecayconst();
    me(0) *= 1./3.*MBs*Fbs*Fbs;
    me(1) *= -5./24.*KBs*MBs*Fbs*Fbs;
    me(2) *= 1./24.*KBs*MBs*Fbs*Fbs;
    me(3) *= 1./4.*KBs*MBs*Fbs*Fbs;
    me(4) *= 1./12.*KBs*MBs*Fbs*Fbs;
#if SUSYFIT_DEBUG & 1
    std::cout << "Bs: me(0) = " << me(0)  << std::endl;
#endif
    
    switch(order) {
        case FULLNLO:
            return((*(allcoeff[LO]) + *(allcoeff[NLO])) * me / HCUT);
        case LO:
            return((*(allcoeff[LO])) * me / HCUT);
        default:
            throw std::runtime_error("AmpDB2::AmpBs(): order not implemented"); 
    }
}

gslpp::complex AmpDB2::PBd()
{
    double mbpole = mySM.Mbar2Mp(mySM.getQuarks(QCD::BOTTOM).getMass());
    double Mw = mySM.Mw();
    double kappa = -2. * M_PI * mbpole * mbpole / 
    (3. * Mw * Mw * mySM.getFlavour().getHDF2().getUDF2().etabS0(mySM.getBBd().getMu()));
    
    double n[13] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
    
    n[0] = 0.1797;
    n[1] = 0.1391;
    n[5] = 1.0116;
    n[6] = 0.0455;
    n[8] = -0.0714;
    n[10] = -0.3331;
    
    double B1 = mySM.getBBd().getBpars()(0);
    double B2 = mySM.getBBd().getBpars()(1);
    
    gslpp::complex PBd = -2. * kappa / mySM.getCBd() * 
            (gslpp::complex(1,2.*mySM.getPhiBd(),true) * (n[0] + (n[5] * B2 + n[10])/B1)
            - gslpp::complex(1./mySM.getCKM().computeRt(),mySM.getCKM().computeBeta()+2.*mySM.getPhiBd(),true)
            * (n[1] + (n[6] * B2 + n[11])/B1)
            + gslpp::complex(1./mySM.getCKM().computeRt()/mySM.getCKM().computeRt(),2.*(mySM.getCKM().computeBeta()+mySM.getPhiBd()),true)
            * (n[2] + (n[7] * B2 + n[12])/B1));

    return PBd;
}

gslpp::complex AmpDB2::PBs()
{
    double mbpole = mySM.Mbar2Mp(mySM.getQuarks(QCD::BOTTOM).getMass());
    double Mw = mySM.Mw();
    double kappa = -2. * M_PI * mbpole * mbpole / 
    (3. * Mw * Mw * mySM.getFlavour().getHDF2().getUDF2().etabS0(mySM.getBBs().getMu()));
    
    double n[13] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
    
    n[0] = 0.1797;
    n[1] = 0.1391;
    n[5] = 1.0116;
    n[6] = 0.0455;
    n[8] = -0.0714;
    n[10] = -0.3331;
    
    double B1 = mySM.getBBs().getBpars()(0);
    double B2 = mySM.getBBs().getBpars()(1);
    
    gslpp::complex PBs = -2. * kappa / mySM.getCBs() * 
            (gslpp::complex(1,2.*mySM.getPhiBs(),true) * (n[0] + (n[5] * B2 + n[10])/B1)
            - gslpp::complex(1./mySM.getCKM().computeRts(),-mySM.getCKM().computeBetas()+2.*mySM.getPhiBs(),true)
            * (n[1] + (n[6] * B2 + n[11])/B1)
            + gslpp::complex(1./mySM.getCKM().computeRts()/mySM.getCKM().computeRts(),2.*(-mySM.getCKM().computeBetas()+mySM.getPhiBs()),true)
            * (n[2] + (n[7] * B2 + n[12])/B1));

    return PBs;
}

