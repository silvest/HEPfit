/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include <math.h>
#include <sstream>
#include <stdexcept>
#include "StandardModelMatching.h"
#include "SUSY.h"
#include "SUSYSpectrum.h"
#include "EWSUSY.h"


const std::string SUSY::SUSYvars[NSUSYvars] = {
    "m1r", "m1i", "m2r", "m2i", "m3" , "muHr", "muHi", "mHptree", "tanb", "Q_SUSY"
};

SUSY::SUSY()
: StandardModel(), 
        msQhat2(3,3,0.), msUhat2(3,3,0.), msDhat2(3,3,0.),msLhat2(3,3,0.), msNhat2(3,3,0.), msEhat2(3,3,0.),
        TUhat(3,3,0.), TDhat(3,3,0.), TNhat(3,3,0.), TEhat(3,3,0.),
        mch(2,0.), mneu(4,0.), m_su2(6,0.), m_sd2(6,0.), m_sdresum2(6,0.), m_sn2(6,0.), m_se2(6,0.),
        U(2,2,0.), V(2,2,0.), N(4,4,0.),
        Ru(6,6,0.), Rd(6,6,0.), Rdresum(6,6,0.), Rn(6,6,0.), Rl(6,6,0.),SUSYM(*this)
{
    SMM.setObj((StandardModelMatching&) SUSYM.getObj());
    ModelParamMap.insert(std::make_pair("m1r", std::cref(m1.real())));
    ModelParamMap.insert(std::make_pair("m1i", std::cref(m1.imag())));
    ModelParamMap.insert(std::make_pair("m2r", std::cref(m2.real())));
    ModelParamMap.insert(std::make_pair("m2i", std::cref(m2.imag())));
    ModelParamMap.insert(std::make_pair("m3", std::cref(m3)));
    ModelParamMap.insert(std::make_pair("muHr", std::cref(muH.real())));
    ModelParamMap.insert(std::make_pair("muHi", std::cref(muH.imag())));
    ModelParamMap.insert(std::make_pair("mHptree", std::cref(mHptree)));
    ModelParamMap.insert(std::make_pair("tanb", std::cref(tanb)));
    ModelParamMap.insert(std::make_pair("Q_SUSY", std::cref(Q_SUSY)));
    
    flag_fh = false;
    flag_h = false;
    flag_g = false;
    flag_ch = false;
    flag_ne = false;
}

SUSY::~SUSY() 
{
    if (IsModelInitialized()) {
        if (mySUSYSpectrum != NULL) delete(mySUSYSpectrum);
        if (myEWSUSY != NULL) delete(myEWSUSY);
    }
}
///////////////////////////////////////////////////////////////////////////
// Initialization

bool SUSY::InitializeModel()
{
    if (!flag_h) mySUSYSpectrum = new SUSYSpectrum(*this);
    else mySUSYSpectrum = NULL;
    myEWSUSY = new EWSUSY(*this);
    setFlagStr("Mw", "NORESUM");
    setModelInitialized(StandardModel::InitializeModel());
    setModelSUSY();
    return(true);
}

///////////////////////////////////////////////////////////////////////////
// Parameters 

bool SUSY::Init(const std::map<std::string, double>& DPars)
{
    return(StandardModel::Init(DPars));
}

bool SUSY::PreUpdate()
{
    if(!StandardModel::PreUpdate()) return (false);

    return (true);
}

bool SUSY::Update(const std::map<std::string, double>& DPars)
{
    if(!PreUpdate()) return (false);

    UpdateError = false;

    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);

    if (UpdateError) return (false);

    if(!PostUpdate()) return (false);

    return (true);
}

bool SUSY::PostUpdate()
{
    if(!StandardModel::PostUpdate()) return (false);
    /* Set the squark and slepton mass matrices and the trilinear-coupling matrices */
    SetSoftTerms();
    computeYukawas();

    /* use approximate GUT relation if M1 & M2 are zero */
    if(m1.abs() == 0. && m2.abs() == 0.) {
        m1.real() = m3/6.;
        m2.real() = m3/3.;
    }
    
    /* Compute Higgs and sparticle spectra with FeynHiggs */
    if (IsFlag_FH()) {
    }
    else {
    /* Compute Higgs and sparticle spectra without FeynHiggs */
        /* sfermions */
        for (int i = 0; i < 6; i++) {
            m_sn2(i) = 0.;// heavy decoupled masses for i=3-5
            m_se2(i) = 0.;
            m_su2(i) = 0.;
            m_sd2(i) = 0.;
            m_sdresum2(i) = 0.;
            for (int j = 0; j < 6; j++) {
                /* R: first (second) index for mass (gauge) eigenstates */
                /* UASf: second (third) index for gauge (mass) eigenstates */
                Rn.assign(i,j, 0.);
                Rl.assign(i,j, 0.);
                Ru.assign(i,j, 0.);
                Rd.assign(i,j, 0.);
                Rdresum.assign(i,j, 0.);
            }
        }

        /* charginos */
        for (int i = 0; i < 2; i++) {
            mch(i) = 0.;
            for (int j = 0; j < 2; j++) {
                /* U and V: first (second) index for mass (gauge) eigenstates */
                /* Ucha and VCha: first (second) index for gauge (mass) eigenstates */
                U.assign(i,j, 0.);
                V.assign(i,j, 0.);
            }
        }

        /* neutralinos */
        for (int i = 0; i < 4; i++) {
            mneu(i) = 0.;
            for (int j = 0; j < 4; j++)
                /* N: first (second) index for mass (gauge) eigenstates */
                /* Zneu: first (second) index for gauge (mass) eigenstates */
                N.assign(i,j, 0.);
        }
        
    if(!mySUSYSpectrum->CalcHiggs(mh,saeff)) return (false);
    if(!mySUSYSpectrum->CalcChargino(U,V,mch)) return (false);
    if(!mySUSYSpectrum->CalcNeutralino(N,mneu)) return (false);
    if(!mySUSYSpectrum->CalcSup(Ru,m_su2)) return (false);
    mySUSYSpectrum->SortSfermionMasses(m_su2, Ru);
    if(!mySUSYSpectrum->CalcSdown(Rd,m_sd2)) return (false);
    mySUSYSpectrum->SortSfermionMasses(m_sd2, Rd);
    if(!mySUSYSpectrum->CalcSneutrino(Rn,m_sn2)) return (false);
    mySUSYSpectrum->SortSfermionMasses(m_sn2, Rn);
    if(!mySUSYSpectrum->CalcSelectron(Rl,m_se2)) return (false);
    mySUSYSpectrum->SortSfermionMasses(m_se2, Rl);
    }

    /* Set the mass of the SM-like Higgs */
    /* allowed range for the use of EWSMApproximateFormulae class */
    if (mh[0] < 10. || mh[0] > 1000.) {
        std::stringstream out;
        out << mh[0];
        throw std::runtime_error("SUSY::PostUpdate(): mh=" + out.str() + " is out of range for EWSMApproximateFormulae");
        return (false);
    }
    mHl = mh[0];
    
    if( Q_SUSY == -1 || Q_SUSY == 0) Q_SUSY = sqrt( sqrt(m_su2(2) * m_su2(5)) );

    /* For EWSUSY class */
    myEWSUSY->SetRosiekParameters();

    /* Necessary for updating SUSY and SUSY-derived parameters in SUSYMatching */
    /* For SUSY Models only: The SM Matching needs to be updated here since StandardModel::PostUpdate() will not do it since the Higgs masses need to be computed. */
    SUSYM.getObj().updateSMParameters();
    SUSYM.getObj().updateSUSYParameters();

    SUSYM.getObj().Comp_mySUSYMQ();

    if (IsFlag_ne()) SUSYM.getObj().Comp_VdDNL(0);
    if (IsFlag_ne()) SUSYM.getObj().Comp_VdDNR(0);
    if (IsFlag_ch()) SUSYM.getObj().Comp_VdUCL();
    if (IsFlag_ch()) SUSYM.getObj().Comp_VdUCR(0);

    SUSYM.getObj().Comp_DeltaMd();
    SUSYM.getObj().Comp_DeltaDL();
    SUSYM.getObj().Comp_Eps_J();
    SUSYM.getObj().Comp_Lambda0EpsY();
    SUSYM.getObj().Comp_mySUSY_CKM();

    if (IsFlag_h()) {
        SUSYM.getObj().Comp_PHLR();
        SUSYM.getObj().Comp_VUDHH();
        SUSYM.getObj().Comp_PHRL();
    }
    if (IsFlag_ne()) {
        SUSYM.getObj().Comp_VdDNL(1);
        SUSYM.getObj().Comp_VdDNR(1);
        SUSYM.getObj().Comp_VuUN();
    }
    if (IsFlag_ch()) {
        SUSYM.getObj().Comp_VdUCR(1);
        SUSYM.getObj().Comp_VuDCL();
        SUSYM.getObj().Comp_VuDCR();
    }

    return (true);
}

void SUSY::setParameter(const std::string name, const double& value)
{
    if (name.compare("m1r") == 0)
        m1.real() = value;
    else if (name.compare("m1i") == 0)
        m1.imag() = value;
    else if (name.compare("m2r") == 0)
        m2.real() = value;
    else if (name.compare("m2i") == 0)
        m2.imag() = value;
    else if (name.compare("m3") == 0)
        m3 = value;
    else if (name.compare("muHr") == 0)
        muH.real() = value;
    else if (name.compare("muHi") == 0)
        muH.imag() = value;
    else if (name.compare("mHptree") == 0)
        mHptree = value;
    else if (name.compare("tanb") == 0)
        SetTanb(value);
    else if (name.compare("Q_SUSY") == 0)
        Q_SUSY = value;
    else
        StandardModel::setParameter(name, value);
}

bool SUSY::CheckParameters(const std::map<std::string, double>& DPars)
{
    for(int i=0; i<NSUSYvars; i++)
        if(DPars.find(SUSYvars[i])==DPars.end()) {
            std::cout << "ERROR: missing mandatory SUSY parameter " << SUSYvars[i] << std::endl;
            raiseMissingModelParameterCount();
            addMissingModelParameter(SUSYvars[i]);
        }
    return(StandardModel::CheckParameters(DPars));
}

void SUSY::SetTanb(const double tanb)
{
    this->tanb = tanb;
    if (tanb < 0.)
        throw std::runtime_error("SUSY::setTanb(): Negative tanb is not allowed");

    /* cosb and sinb are defined to be positive, corresponding to 0<=beta<=pi/2 */
    cosb = sqrt(1. / (1. + tanb * tanb));
    sinb = tanb * cosb;
}

void SUSY::computeYukawas()
{
    /* initializations */
    Yu = gslpp::matrix<gslpp::complex>::Id(3);
    Yd = gslpp::matrix<gslpp::complex>::Id(3);
    Ye = gslpp::matrix<gslpp::complex>::Id(3);
    Yn = gslpp::matrix<gslpp::complex>::Id(3);

    /* Convert the top-quark pole mass to the MSbar mass */
    double mtbar = Mp2Mbar(mtpole, FULLNLO);

    double Q_SUSY_temp = Q_SUSY;
    if( Q_SUSY == -1 || Q_SUSY == 0) Q_SUSY = sqrt( sqrt(msQhat2(2,2).abs() * msUhat2(2,2).abs() ) );

    for (int i = 0; i < 3; i++) {
        /* Run the quark masses to scale Q */
        if (i != 2)
            mu_Q[i] = Mrun(Q_SUSY, getQuarks((quark)(UP + 2 * i)).getMass_scale(),
                           getQuarks((quark)(UP + 2 * i)).getMass(), FULLNLO);
        else
            mu_Q[i] = Mrun(Q_SUSY, mtbar, FULLNLO);
        md_Q[i] = Mrun(Q_SUSY, getQuarks((quark)(DOWN + 2 * i)).getMass_scale(),
                       getQuarks((quark)(DOWN + 2 * i)).getMass(), FULLNLO);
        me_Q[i] = getLeptons((lepton)(ELECTRON + 2 * i)).getMass();
        mn_Q[i] = getLeptons((lepton)(NEUTRINO_1 + 2 * i)).getMass();

        /* Convert MSbar to DRbar */
        mu_Q[i] = MS2DRqmass(Q_SUSY, mu_Q[i]);
        md_Q[i] = MS2DRqmass(Q_SUSY, md_Q[i]);

        Yu.assign(i, i, mu_Q[i] / v2() * sqrt(2.));
        Yd.assign(i, i, md_Q[i] / v1() * sqrt(2.));
        Ye.assign(i, i, me_Q[i] / v1() * sqrt(2.));
        Yn.assign(i, i, mn_Q[i] / v2() * sqrt(2.));
    }

    Yu = myCKM.getCKM().transpose()*Yu;
    Yn = Yn * myPMNS.getPMNS().hconjugate();

    Q_SUSY = Q_SUSY_temp;
    
}

void SUSY::SetSoftTerms()
{
    // MsQ2, MsU2, MsD2, MsL2, MsN2, MsE2, TU, TD, TN and TE are set to 0 in the constructor.
    // See also GeneralSUSY::SetSoftTerms().
}


///////////////////////////////////////////////////////////////////////////
// Flags

bool SUSY::setFlag(const std::string name, const bool value)
{
    bool res = false;
    if(name.compare("Flag_H") == 0) {
        flag_h = value;
        res = true;
    }
    else if(name.compare("Flag_g") == 0) {
        flag_g = value;
        res = true;
    }
    else if(name.compare("Flag_Chi") == 0) {
        flag_ch = value;
        res = true;
    }
    else if(name.compare("Flag_Chi0") == 0) {
        flag_ne = value;
        res = true;
    }
    else if(name.compare("Flag_FH") == 0) {
        flag_fh = value;
        res = true;
    }
    else
        res = StandardModel::setFlag(name,value);

    return(res);
}


///////////////////////////////////////////////////////////////////////////

double SUSY::v1() const
{
    return (v() * cosb);
}

double SUSY::v2() const
{
    return (v() * sinb);
}

double SUSY::getMGl() const
{
    return m3;
}


///////////////////////////////////////////////////////////////////////////

double SUSY::Mw() const
{
    return StandardModel::Mw();
    //return myEWSUSY->Mw_MSSM();
}

double SUSY::Mw_dRho() const
{
    return 0.;
}
