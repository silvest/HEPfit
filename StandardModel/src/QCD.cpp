/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <math.h>
#include <map>
#include <stdexcept>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <TF1.h>
#include <Math/WrappedTF1.h>
#include <Math/BrentRootFinder.h>
#include <algorithm>
#include "QCD.h"
#include "gslpp_special_functions.h"

std::string QCD::QCDvars[NQCDvars] = {
    "AlsM", "MAls",
    "mup", "mdown", "mcharm", "mstrange", "mtop", "mbottom",
    "muc", "mub", "mut"};

QCD::QCD()
{
    setModelName("QCD");
    FlagCsi = true;
    std::cout << "Warning: use of pole top mass as input is now deprecated. By default, the MSbar mass is used."  << std::endl;
    std::cout << "To use the pole mass, set FlagMtPole to true. (to be removed in future releases)" << std::endl;
    FlagMtPole = false;
    FlagMpole2MbarNumeric = true;
    computeFBd = false;
    computeFBp = false;
    computeBd = false;
    computeBs = false;
    computemt = false;
    requireYu = false;
    requireYd = false;
    Nc = 3.;
    TF = 0.5;
    CF = Nc / 2. - 1. / (2. * Nc);
    CA = Nc;
    dFdA_NA = Nc * (Nc * Nc + 6.) / 48.;
    dAdA_NA = Nc * Nc * (Nc * Nc + 36.) / 24.;
    dFdF_NA = (Nc * Nc - 6. + 18. / Nc / Nc) / 96.;
    NA = Nc * Nc - 1.;

    //    Particle(std::string name, double mass, double mass_scale = 0., double width = 0., double charge = 0.,double isospin = 0.);
    quarks[UP] = Particle("UP", 0., 2., 0., 2. / 3., .5);
    quarks[CHARM] = Particle("CHARM", 0., 0., 0., 2. / 3., .5);
    quarks[TOP] = Particle("TOP", 0., 0., 0., 2. / 3., .5);
    quarks[DOWN] = Particle("DOWN", 0., 2., 0., -1. / 3., -.5);
    quarks[STRANGE] = Particle("STRANGE", 0., 2., 0., -1. / 3., -.5);
    quarks[BOTTOM] = Particle("BOTTOM", 0., 0., 0., -1. / 3., -.5);

    zeta2 = gslpp_special_functions::zeta(2);
    zeta3 = gslpp_special_functions::zeta(3);
    for (int i = 0; i < CacheSize; i++)
    {
        for (int j = 0; j < 9; j++)
            als_cache[j][i] = 0.;
        for (int j = 0; j < 4; j++)
            logLambda5_cache[j][i] = 0.;
        for (int j = 0; j < 10; j++)
            mrun_cache[j][i] = 0.;
        for (int j = 0; j < 6; j++)
            mp2mbar_cache[j][i] = 0.;
    }

    ModelParamMap.insert(std::make_pair("AlsM", std::cref(AlsM)));
    ModelParamMap.insert(std::make_pair("MAls", std::cref(MAls)));
    ModelParamMap.insert(std::make_pair("mup", std::cref(quarks[UP].getMass())));
    ModelParamMap.insert(std::make_pair("mdown", std::cref(quarks[DOWN].getMass())));
    ModelParamMap.insert(std::make_pair("mcharm", std::cref(quarks[CHARM].getMass())));
    ModelParamMap.insert(std::make_pair("mstrange", std::cref(quarks[STRANGE].getMass())));
    if (FlagMtPole)
        ModelParamMap.insert(std::make_pair("mtop", std::cref(mtpole)));
    else
        ModelParamMap.insert(std::make_pair("mtop", std::cref(quarks[TOP].getMass())));

    ModelParamMap.insert(std::make_pair("mbottom", std::cref(quarks[BOTTOM].getMass())));
    ModelParamMap.insert(std::make_pair("muc", std::cref(muc)));
    ModelParamMap.insert(std::make_pair("mub", std::cref(mub)));
    ModelParamMap.insert(std::make_pair("mut", std::cref(mut)));

    unknownParameterWarning = true;
    realorder = LO;
}

const std::string QCD::orderToString(const orders order) const
{
    switch (order)
    {
    case LO:
        return "LO";
    case NLO:
        return "NLO";
    case FULLNLO:
        return "FULLNLO";
    case NNLO:
        return "NNLO";
    case FULLNNLO:
        return "FULLNNLO";
    case NNNLO:
        return "NNNLO";
    case FULLNNNLO:
        return "FULLNNNLO";
    default:
        throw std::runtime_error("QCD::orderToString(): order not implemented.");
    }
}

////////////////////////////////////////////////////////////////////////

bool QCD::Init(const std::map<std::string, double> &DPars)
{
    bool check = CheckParameters(DPars);
    if (!check)
        return (check);
    check *= Update(DPars) * QCDsuccess;
    unknownParameterWarning = false;
    return (check);
}

bool QCD::PreUpdate()
{
    computemt = false;
    QCDsuccess = true;
    return (true);
}

bool QCD::Update(const std::map<std::string, double> &DPars)
{
    QCDsuccess = true;
    if (!PreUpdate())
        return (false);

    UpdateError = false;

    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
    {
        setParameter(it->first, it->second);
    }

    if (UpdateError || !QCDsuccess)
        return (false);

    if (!PostUpdate() || !QCDsuccess)
        return (false);

    return (true);
}

bool QCD::PostUpdate()
{
    if (computeFBd)
        mesonsMap.at(QCD::B_D).setDecayconst(mesonsMap.at(QCD::B_S).getDecayconst() / mesonsMap.at(QCD::B_D).getFBsoFBd());
    if (computeFBp)
        mesonsMap.at(QCD::B_P).setDecayconst(mesonsMap.at(QCD::B_S).getDecayconst() / mesonsMap.at(QCD::B_P).getFBsoFBd()); /**** FOR NOW FB+ = FBd ****/
    if (computeBs && FlagCsi)
    {
        BParameterMap.at("BBs").setBpars(0, BParameterMap.at("BBs").getFBsSqrtBBs1() * BParameterMap.at("BBs").getFBsSqrtBBs1() / mesonsMap.at(QCD::B_S).getDecayconst() / mesonsMap.at(QCD::B_S).getDecayconst());
        BParameterMap.at("BBs").setBpars(1, BParameterMap.at("BBs").getFBsSqrtBBs2() * BParameterMap.at("BBs").getFBsSqrtBBs2() / mesonsMap.at(QCD::B_S).getDecayconst() / mesonsMap.at(QCD::B_S).getDecayconst());
        BParameterMap.at("BBs").setBpars(2, BParameterMap.at("BBs").getFBsSqrtBBs3() * BParameterMap.at("BBs").getFBsSqrtBBs3() / mesonsMap.at(QCD::B_S).getDecayconst() / mesonsMap.at(QCD::B_S).getDecayconst());
        BParameterMap.at("BBs").setBpars(3, BParameterMap.at("BBs").getFBsSqrtBBs4() * BParameterMap.at("BBs").getFBsSqrtBBs4() / mesonsMap.at(QCD::B_S).getDecayconst() / mesonsMap.at(QCD::B_S).getDecayconst());
        BParameterMap.at("BBs").setBpars(4, BParameterMap.at("BBs").getFBsSqrtBBs5() * BParameterMap.at("BBs").getFBsSqrtBBs5() / mesonsMap.at(QCD::B_S).getDecayconst() / mesonsMap.at(QCD::B_S).getDecayconst());
    }
    if (computeBd)
    {
        if (FlagCsi)
        {
            BParameterMap.at("BBd").setBpars(0, mesonsMap.at(QCD::B_D).getFBsoFBd() * mesonsMap.at(QCD::B_D).getFBsoFBd() * BParameterMap.at("BBs").getBpars()(0) / BParameterMap.at("BBd").getcsi() / BParameterMap.at("BBd").getcsi());
            BParameterMap.at("BBd").setBBsoBBd(BParameterMap.at("BBs").getBpars()(0) / BParameterMap.at("BBd").getBpars()(0));
            BParameterMap.at("BBd").setBpars(1, BParameterMap.at("BBd").getFBdSqrtBBd2() * BParameterMap.at("BBd").getFBdSqrtBBd2() / mesonsMap.at(QCD::B_D).getDecayconst() / mesonsMap.at(QCD::B_D).getDecayconst());
            BParameterMap.at("BBd").setBpars(2, BParameterMap.at("BBd").getFBdSqrtBBd3() * BParameterMap.at("BBd").getFBdSqrtBBd3() / mesonsMap.at(QCD::B_D).getDecayconst() / mesonsMap.at(QCD::B_D).getDecayconst());
            BParameterMap.at("BBd").setBpars(3, BParameterMap.at("BBd").getFBdSqrtBBd4() * BParameterMap.at("BBd").getFBdSqrtBBd4() / mesonsMap.at(QCD::B_D).getDecayconst() / mesonsMap.at(QCD::B_D).getDecayconst());
            BParameterMap.at("BBd").setBpars(4, BParameterMap.at("BBd").getFBdSqrtBBd5() * BParameterMap.at("BBd").getFBdSqrtBBd5() / mesonsMap.at(QCD::B_D).getDecayconst() / mesonsMap.at(QCD::B_D).getDecayconst());
        }
        else
            BParameterMap.at("BBd").setBpars(0, BParameterMap.at("BBs").getBpars()(0) / BParameterMap.at("BBd").getBBsoBBd());
    }
    if (computemt)
    {
        if(FlagMtPole){
            quarks[TOP].setMass(163.); //Cannot use MSbar mass for Alphas thresholds since it is not yet updated
            quarks[TOP].setMass(Mp2Mbar(mtpole, TOP, FULLNNLO));
            quarks[TOP].setMass_scale(quarks[TOP].getMass());
        }
        else
            mtpole = Mbar2Mp(quarks[TOP].getMass(), TOP, FULLNNLO);
   }
    return (QCDsuccess);
}

void QCD::addParameters(std::vector<std::string> params_i)
{
    for (std::vector<std::string>::iterator it = params_i.begin(); it < params_i.end(); it++)
    {
        if (optionalParameters.find(*it) == optionalParameters.end())
        {
            optionalParameters[*it] = 0.;
            ModelParamMap.insert(std::make_pair(*it, std::cref(optionalParameters[*it])));
        }
    }
}

void QCD::initializeBParameter(std::string name_i) const
{
    if (BParameterMap.find(name_i) != BParameterMap.end())
        return;

    if (name_i.compare("BBs") == 0 || name_i.compare("BBd") == 0)
    {
        BParameterMap.insert(std::pair<std::string, BParameter>("BBs", BParameter(5, "BBs")));
        BParameterMap.at("BBs").setFlagCsi(FlagCsi);
        BParameterMap.at("BBs").ModelParameterMapInsert(ModelParamMap);
        computeBs = true;
        initializeMeson(QCD::B_S);

        BParameterMap.insert(std::pair<std::string, BParameter>("BBd", BParameter(5, "BBd")));
        BParameterMap.at("BBd").setFlagCsi(FlagCsi);
        BParameterMap.at("BBd").ModelParameterMapInsert(ModelParamMap);
        computeBd = true;
        initializeMeson(QCD::B_D);
        return;
    }
    if (name_i.compare("BBs_subleading") == 0)
    {
        BParameterMap.insert(std::pair<std::string, BParameter>(name_i, BParameter(2, name_i)));
        BParameterMap.at(name_i).ModelParameterMapInsert(ModelParamMap);
        return;
    }
    if (name_i.compare("BBd_subleading") == 0)
    {
        BParameterMap.insert(std::pair<std::string, BParameter>(name_i, BParameter(2, name_i)));
        BParameterMap.at(name_i).ModelParameterMapInsert(ModelParamMap);
        return;
    }
    if (name_i.compare("BD") == 0)
    {
        BParameterMap.insert(std::pair<std::string, BParameter>(name_i, BParameter(5, name_i)));
        BParameterMap.at(name_i).ModelParameterMapInsert(ModelParamMap);
        initializeMeson(QCD::D_0);
        return;
    }
    if (name_i.compare("BK") == 0)
    {
        BParameterMap.insert(std::pair<std::string, BParameter>(name_i, BParameter(5, name_i)));
        BParameterMap.at(name_i).ModelParameterMapInsert(ModelParamMap);
        initializeMeson(QCD::K_0);
        return;
    }
    if (name_i.compare("BKd1") == 0)
    {
        BParameterMap.insert(std::pair<std::string, BParameter>(name_i, BParameter(10, name_i)));
        BParameterMap.at(name_i).ModelParameterMapInsert(ModelParamMap);
        initializeMeson(QCD::K_0);
        initializeMeson(QCD::K_P);
        initializeMeson(QCD::K_S);
        initializeMeson(QCD::P_0);
        initializeMeson(QCD::P_P);
        return;
    }
    if (name_i.compare("BKd3") == 0)
    {
        BParameterMap.insert(std::pair<std::string, BParameter>(name_i, BParameter(10, name_i)));
        BParameterMap.at(name_i).ModelParameterMapInsert(ModelParamMap);
        initializeMeson(QCD::K_0);
        initializeMeson(QCD::K_P);
        initializeMeson(QCD::P_0);
        initializeMeson(QCD::P_P);
        return;
    }
}

void QCD::initializeMeson(const QCD::meson meson_i) const
{
    if (mesonsMap.find(meson_i) != mesonsMap.end())
        return;

    mesonsMap.insert(std::pair<const QCD::meson, Meson>(meson_i, Meson()));

    if (meson_i == QCD::P_0)
        mesonsMap.at(meson_i).setName("P_0");
    else if (meson_i == QCD::P_P)
        mesonsMap.at(meson_i).setName("P_P");
    else if (meson_i == QCD::K_0)
        mesonsMap.at(meson_i).setName("K_0");
    else if (meson_i == QCD::K_P)
        mesonsMap.at(meson_i).setName("K_P");
    else if (meson_i == QCD::K_S)
        mesonsMap.at(meson_i).setName("K_S");
    else if (meson_i == QCD::D_0)
        mesonsMap.at(meson_i).setName("D_0");
    else if (meson_i == QCD::D_P)
        mesonsMap.at(meson_i).setName("D_P");
    else if (meson_i == QCD::D_S)
        mesonsMap.at(meson_i).setName("D_S");
    else if (meson_i == QCD::B_D)
        mesonsMap.at(meson_i).setName("B_D");
    else if (meson_i == QCD::B_P)
        mesonsMap.at(meson_i).setName("B_P");
    else if (meson_i == QCD::B_S)
        mesonsMap.at(meson_i).setName("B_S");
    else if (meson_i == QCD::B_C)
        mesonsMap.at(meson_i).setName("B_C");
    else if (meson_i == QCD::PHI)
        mesonsMap.at(meson_i).setName("PHI");
    else if (meson_i == QCD::K_star)
        mesonsMap.at(meson_i).setName("K_star");
    else if (meson_i == QCD::K_star_P)
        mesonsMap.at(meson_i).setName("K_star_P");
    else if (meson_i == QCD::D_star_P)
        mesonsMap.at(meson_i).setName("D_star_P");
    else if (meson_i == QCD::RHO)
        mesonsMap.at(meson_i).setName("RHO");
    else if (meson_i == QCD::RHO_P)
        mesonsMap.at(meson_i).setName("RHO_P");
    else if (meson_i == QCD::OMEGA)
        mesonsMap.at(meson_i).setName("OMEGA");
    else
    {
        std::stringstream out;
        out << meson_i;
        throw std::runtime_error("QCD::initializeMeson() meson " + out.str() + " not implemented");
    }

    if (meson_i == QCD::B_D)
        computeFBd = true;
    if (meson_i == QCD::B_P)
        computeFBp = true;

    if ((computeFBd || computeFBp) && (mesonsMap.find(QCD::B_S) == mesonsMap.end()))
        initializeMeson(QCD::B_S);

    mesonsMap.at(meson_i).ModelParameterMapInsert(ModelParamMap);
}

void QCD::setParameter(const std::string name, const double &value)
{
    int notABparameter = 0;
    int notAMesonParameter = 0;

    if (name.compare("AlsM") == 0)
    {
        AlsM = value;
        computemt = true;
    }
    else if (name.compare("MAls") == 0)
    {
        MAls = value;
        computemt = true;
    }
    else if (name.compare("mup") == 0)
    {
        if (value < MEPS)
            UpdateError = true;
        quarks[UP].setMass(value);
    }
    else if (name.compare("mdown") == 0)
    {
        if (value < MEPS)
            UpdateError = true;
        quarks[DOWN].setMass(value);
    }
    else if (name.compare("mcharm") == 0)
    {
        quarks[CHARM].setMass(value);
        quarks[CHARM].setMass_scale(value);
    }
    else if (name.compare("mstrange") == 0)
    {
        if (value < MEPS)
            UpdateError = true;
        quarks[STRANGE].setMass(value);
    }
    else if (name.compare("mtop") == 0)
    {
        if (FlagMtPole)
            mtpole = value;
        else
        {
            quarks[TOP].setMass(value);
            quarks[TOP].setMass_scale(value);
        }
        computemt = true;
    }
    else if (name.compare("mbottom") == 0)
    {
        quarks[BOTTOM].setMass(value);
        quarks[BOTTOM].setMass_scale(value);
    }
    else if (name.compare("muc") == 0)
        muc = value;
    else if (name.compare("mub") == 0)
        mub = value;
    else if (name.compare("mut") == 0)
        mut = value;
    else if (optionalParameters.find(name) != optionalParameters.end())
        setOptionalParameter(name, value);
    else
    {
        if (!BParameterMap.empty())
            for (std::map<std::string, BParameter>::iterator it = BParameterMap.begin(); it != BParameterMap.end(); it++)
                if (it->second.setParameter(name, value))
                    notABparameter += 1;
        if (!mesonsMap.empty())
            for (std::map<const QCD::meson, Meson>::iterator it = mesonsMap.begin(); it != mesonsMap.end(); it++)
                if (it->second.setParameter(name, value))
                    notAMesonParameter += 1;
        if (unknownParameterWarning && !isSliced && notABparameter == 0 && notAMesonParameter == 0)
            if (std::find(unknownParameters.begin(), unknownParameters.end(), name) == unknownParameters.end())
                unknownParameters.push_back(name);
    }
}

bool QCD::CheckParameters(const std::map<std::string, double> &DPars)
{
    for (int i = 0; i < NQCDvars; i++)
        if (DPars.find(QCDvars[i]) == DPars.end())
        {
            std::cout << "ERROR: missing mandatory QCD parameter " << QCDvars[i] << std::endl;
            raiseMissingModelParameterCount();
            addMissingModelParameter(QCDvars[i]);
        }
    for (std::map<std::string, double>::iterator it = optionalParameters.begin(); it != optionalParameters.end(); it++)
    {
        if (DPars.find(it->first) == DPars.end())
        {
            std::cout << "ERROR: missing optional parameter " << it->first << std::endl;
            raiseMissingModelParameterCount();
            addMissingModelParameter(it->first);
        }
    }
    if (!BParameterMap.empty())
    {
        for (std::map<std::string, BParameter>::iterator it = BParameterMap.begin(); it != BParameterMap.end(); it++)
        {
            std::vector<std::string> parameters = it->second.parameterList(it->first);
            for (std::vector<std::string>::iterator it1 = parameters.begin(); it1 != parameters.end(); it1++)
            {
                if (DPars.find(*it1) == DPars.end())
                {
                    std::cout << "ERROR: missing parameter for " << it->first << ": " << *it1 << std::endl;
                    raiseMissingModelParameterCount();
                    addMissingModelParameter(*it1);
                }
            }
        }
    }
    if (!mesonsMap.empty())
    {
        for (std::map<const QCD::meson, Meson>::iterator it = mesonsMap.begin(); it != mesonsMap.end(); it++)
        {
            std::vector<std::string> parameters = it->second.parameterList(it->second.getName());
            for (std::vector<std::string>::iterator it1 = parameters.begin(); it1 != parameters.end(); it1++)
            {
                if (DPars.find(*it1) == DPars.end())
                {
                    std::cout << "ERROR: missing parameter for " << mesonsMap.at(it->first).getName() << ": " << *it1 << std::endl;
                    raiseMissingModelParameterCount();
                    addMissingModelParameter(*it1);
                }
            }
        }
    }
    if (getMissingModelParametersCount() > 0)
        return false;
    else
        return true;
}

////////////////////////////////////////////////////////////////////////

bool QCD::setFlag(const std::string name, const bool value)
{
    bool res = false;
    if (name.compare("FlagCsi") == 0)
    {
        FlagCsi = value;
        if (computeBs)
            BParameterMap.at("BBs").setFlagCsi(FlagCsi);
        if (computeBd)
            BParameterMap.at("BBd").setFlagCsi(FlagCsi);
        res = true;
    }
    else if (name.compare("FlagMtPole") == 0)
    {
        FlagMtPole = value;
        if (FlagMtPole) {
            ModelParamMap.at("mtop")=std::cref(mtpole);
            std::cout << "WARNING: use of pole top mass as input is now deprecated. This feature will be removed in future releases. " << std::endl;
            std::cout << "Using anyway the pole mass as requested." << std::endl;
        }
        else
            ModelParamMap.at("mtop")=std::cref(quarks[TOP].getMass());
        res = true;
    }
    else if (name.compare("FlagMpole2MbarNumeric") == 0)
    {
        FlagMpole2MbarNumeric = value;
        res = true;
    }
    return res;
}

bool QCD::setFlagStr(const std::string name, const std::string value)
{
    std::cout << "WARNING: wrong name or value for ModelFlag " << name << std::endl;
    return (false);
}

bool QCD::CheckFlags() const
{
    return (true);
}

////////////////////////////////////////////////////////////////////////

const double QCD::Thresholds(const int i) const
{
    if (!(mut > mub && mub > muc))
    {
        QCDsuccess = false;
        std::cout << "inverted thresholds in QCD::Thresholds()! QCDsuccess set to false" << std::endl;
    }

    switch (i)
    {
    case 0:
        return (1.0E10);
    case 1:
        return (mut);
    case 2:
        return (mub);
    case 3:
        return (muc);
    default:
        return (0.);
    }
}

const double QCD::AboveTh(const double mu) const
{
    int i;
    for (i = 4; i >= 0; i--)
        if (mu < Thresholds(i))
            return (Thresholds(i));

    QCDsuccess = false;
    std::cout << "Error in QCD::AboveTh()! QCDsuccess set to false" << std::endl;
    return (0.);
}

const double QCD::BelowTh(const double mu) const
{
    int i;
    for (i = 0; i < 5; i++)
        if (mu >= Thresholds(i))
            return (Thresholds(i));

    QCDsuccess = false;
    std::cout << "Error in QCD::BelowTh()! QCDsuccess set to false" << std::endl;
    return (0.);
}

const double QCD::Nf(const double mu) const
{
    int i;
    for (i = 1; i < 5; i++)
        if (mu >= Thresholds(i))
            return (7. - (double)i);

    QCDsuccess = false;
    std::cout << "Error in QCD::Nf(" << mu <<")! QCDsuccess set to false" << std::endl;
    return (0.);
}

void QCD::CacheShift(double cache[][CacheSize], int n) const
{
    int i, j;
    for (i = CacheSize - 1; i > 0; i--)
        for (j = 0; j < n; j++)
            cache[j][i] = cache[j][i - 1];
}

void QCD::CacheShift(int cache[][CacheSize], int n) const
{
    int i, j;
    for (i = CacheSize - 1; i > 0; i--)
        for (j = 0; j < n; j++)
            cache[j][i] = cache[j][i - 1];
}

////////////////////////////////////////////////////////////////////////

const double QCD::Beta0(const double nf) const
{
    return ((11. * CA - 4. * TF * nf) / 3.);
}

const double QCD::Beta1(const double nf) const
{
    return (34. / 3. * CA * CA - (20. / 3. * CA + 4. * CF) * TF * nf);
}

const double QCD::Beta2(const double nf) const
{
    return (2857. / 54. * CA * CA * CA - (1415. / 27. * CA * CA + 205. / 9. * CF * CA - 2. * CF * CF) * TF * nf +
            (158. / 27. * CA + 44. / 9. * CF) * TF * TF * nf * nf);
}

// Czakon, hep-ph/0411261, eq. (14)
const double QCD::Beta3(const double nf) const
{
    return (CA * CF * TF * TF * nf * nf * (17152. / 243. + 448. / 9. * zeta3) +
            CA * CF * CF * TF * nf * (-4204. / 27. + 352. / 9. * zeta3) +
            424. / 243. * CA * TF * TF * TF * nf * nf * nf + CA * CA * CF * TF * nf * (7073. / 243. - 656. / 9. * zeta3) + CA * CA * TF * TF * nf * nf * (7930. / 81. + 224. / 9. * zeta3) + 1232. / 243. * CF * TF * TF * TF * nf * nf * nf +
            CA * CA * CA * TF * nf * (-39143. / 81. + 136. / 3. * zeta3) + CA * CA * CA * CA * (150653. / 486. - 44. / 9. * zeta3) + CF * CF * TF * TF * nf * nf * (1352. / 27. - 704. / 9. * zeta3) + 46. * CF * CF * CF * TF * nf +
            nf * dFdA_NA * (512. / 9. - 1664. / 3. * zeta3) + nf * nf * dFdF_NA * (-704. / 9. + 512. / 3. * zeta3) + dAdA_NA * (-80. / 9. + 704. / 3. * zeta3));
}

const double QCD::AlsWithInit(const double mu, const double alsi, const double mu_i, const int nf,
                              const orders order) const
{

    double b1_b0 = Beta1((double)nf) / Beta0((double)nf);
    double v = 1. - Beta0((double)nf) * alsi / 2. / M_PI * log(mu_i / mu);
    double logv = log(v);

    switch (order)
    {
    case LO:
        return (alsi / v);
    case NLO:
        return (-alsi * alsi / 4. / M_PI / v / v * b1_b0 * logv);
    case NNLO:
        return (alsi * alsi * alsi / 4. / 4. / M_PI / M_PI / v / v / v * (Beta2((double)nf) / Beta0((double)nf) * (1. - v) + b1_b0 * b1_b0 * (logv * logv - logv + v - 1.)));
    case NNNLO:
        return (alsi * alsi * alsi * alsi / 4. / 4. / 4. / M_PI / M_PI / M_PI /
                v / v / v / v * (Beta3((double)nf) / Beta0((double)nf) * (1. - v * v) / 2. + b1_b0 * Beta2((double)nf) / Beta0((double)nf) * ((2. * v - 3.) * logv + v * v - v) + b1_b0 * b1_b0 * b1_b0 * (-logv * logv * logv + 2.5 * logv * logv + 2. * (1. - v) * logv - 0.5 * (v - 1.) * (v - 1.))));
    case FULLNLO:
        return (AlsWithInit(mu, alsi, mu_i, nf, LO) + AlsWithInit(mu, alsi, mu_i, nf, NLO));
    case FULLNNLO:
        return (AlsWithInit(mu, alsi, mu_i, nf, LO) + AlsWithInit(mu, alsi, mu_i, nf, NLO) + AlsWithInit(mu, alsi, mu_i, nf, NNLO));
    case FULLNNNLO:
        return (AlsWithInit(mu, alsi, mu_i, nf, LO) + AlsWithInit(mu, alsi, mu_i, nf, NLO) + AlsWithInit(mu, alsi, mu_i, nf, NNLO) + AlsWithInit(mu, alsi, mu_i, nf, NNNLO));
    default:
        throw std::runtime_error("QCD::AlsWithInit(): " + orderToString(order) + " is not implemented.");
    }
}

const double QCD::Als4(const double mu) const
{
    double v = 1. - Beta0(4.) * AlsM / 2. / M_PI * log(MAls / mu);
    return (AlsM / v * (1. - Beta1(4.) / Beta0(4.) * AlsM / 4. / M_PI * log(v) / v));
}

const double QCD::AlsWithLambda(const double mu, const double logLambda,
                                const orders order) const
{
    double nf = Nf(mu);
    double L = 2. * (log(mu) - logLambda);

    // LO contribution
    double b0 = Beta0(nf);
    double b0L = b0 * L;
    double alsLO = 4. * M_PI / b0L;
    if (order == LO)
        return alsLO;

    // NLO contribution
    double b1 = Beta1(nf);
    double log_L = log(L);
    double alsNLO = 4. * M_PI / b0L * (-b1 * log_L / b0 / b0L);
    if (order == NLO)
        return alsNLO;
    if (order == FULLNLO)
        return (alsLO + alsNLO);

    // NNLO contribution
    double b2 = Beta2(nf);
    double alsNNLO = 4. * M_PI / b0L * (1. / b0L / b0L * (b1 * b1 / b0 / b0 * (log_L * log_L - log_L - 1.) + b2 / b0));
    if (order == NNLO)
        return alsNNLO;
    if (order == FULLNNLO)
        return (alsLO + alsNLO + alsNNLO);

    // NNNLO contribution
    double b3 = Beta3(nf);
    double alsNNNLO = 4. * M_PI / b0L * (-1. / b0L / b0L / b0L * (b1 * b1 * b1 / b0 / b0 / b0 * (log_L * log_L * log_L - 5. / 2. * log_L * log_L - 2. * log_L - 0.5) + 3. * b1 * b2 * log_L / b0 / b0 - 0.5 * b3 / b0));
    if (order == NNNLO)
        return alsNNNLO;
    if (order == FULLNNNLO)
        return (alsLO + alsNLO + alsNNLO + alsNNNLO);

    throw std::runtime_error(orderToString(order) + " is not implemented in QCD::AlsWithLambda().");
}

const double QCD::AlsWithLambda(const double mu, const orders order) const
{
    return AlsWithLambda(mu, logLambda(Nf(mu), order), order);
}

const double QCD::NfThresholdCorrections(double mu, double M, double als, int nf, orders order) const
{
    double lmM = 2. * log(mu / M), res = 0.;
    switch (order)
    {
    case FULLNNNLO:
        res += als * als * als / M_PI / M_PI / M_PI * (-564731. / 124416. + 82043. / 27648. * gslpp_special_functions::zeta(3) + 2191. / 576. * lmM + 511. / 576. * lmM * lmM + lmM * lmM * lmM / 216. + ((double)nf - 1.) * (2633. / 31104. - 281. / 1728. * lmM));
    case FULLNNLO:
        res += als * als / M_PI / M_PI * (-11. / 72. + 19. / 24. * lmM + lmM * lmM / 36.);
    case FULLNLO:
        res += als / M_PI * lmM / 6.;
    case LO:
        break;
    default:
        throw std::runtime_error("QCD::ThresholdCorrections(): order not implemented.");
    }
    return (res);
}

const orders QCD::FullOrder(orders order) const
{
    switch (order)
    {
    case LO:
        return (LO);
    case NLO:
        return (FULLNLO);
    case NNLO:
        return (FULLNNLO);
    case NNNLO:
        return (FULLNNNLO);
    default:
        throw std::runtime_error("QCD::FullOrder(): " + orderToString(order) + " is not implemented.");
    }
}

const double QCD::MassOfNf(int nf) const
{
    switch (nf)
    {
    case 6:
        return (quarks[TOP].getMass());
    case 5:
        return (quarks[BOTTOM].getMass());
    case 4:
        return (quarks[CHARM].getMass());
    case 3:
        return (quarks[STRANGE].getMass());
    default:
        throw std::runtime_error("QCD::MassOfNf(): no running masses for light quarks");
    }
}

const double QCD::Als(const double mu, const orders order, const bool Nf_thr) const
{
    switch (order)
    {
    case LO:
        realorder = order;
        return AlsByOrder(mu, LO, Nf_thr);
    case FULLNLO:
        realorder = order;
        return (AlsByOrder(mu, LO, Nf_thr) + AlsByOrder(mu, NLO, Nf_thr));
    case FULLNNLO:
        realorder = order;
        return (AlsByOrder(mu, LO, Nf_thr) + AlsByOrder(mu, NLO, Nf_thr) + AlsByOrder(mu, NNLO, Nf_thr));
    case FULLNNNLO:
        realorder = order;
        return (AlsByOrder(mu, LO, Nf_thr) + AlsByOrder(mu, NLO, Nf_thr) + AlsByOrder(mu, NNLO, Nf_thr) + AlsByOrder(mu, NNNLO, Nf_thr));
    default:
        throw std::runtime_error("QCD::Als(): " + orderToString(order) + " is not implemented.");
    }
}

const double QCD::Als(const double mu, const int Nf, const orders order) const
{
    switch (order)
    {
    case LO:
        realorder = order;
        return AlsByOrder(mu, Nf, LO);
    case FULLNLO:
        realorder = order;
        return (AlsByOrder(mu, Nf, LO) + AlsByOrder(mu, Nf, NLO));
    case FULLNNLO:
        realorder = order;
        return (AlsByOrder(mu, Nf, LO) + AlsByOrder(mu, Nf, NLO) + AlsByOrder(mu, Nf, NNLO));
    case FULLNNNLO:
        realorder = order;
        return (AlsByOrder(mu, Nf, LO) + AlsByOrder(mu, Nf, NLO) + AlsByOrder(mu, Nf, NNLO) + AlsByOrder(mu, Nf, NNNLO));
    default:
        throw std::runtime_error("QCD::Als(): " + orderToString(order) + " is not implemented.");
    }
}

const double QCD::AlsByOrder(const double mu, const orders order, const bool Nf_thr) const
{
    int i, nfAls = (int)Nf(MAls), nfmu = Nf_thr ? (int)Nf(mu) : nfAls;
    double als, alstmp, mutmp;
    orders fullord;

    for (i = 0; i < CacheSize; ++i)
        if ((mu == als_cache[0][i]) && ((double)order == als_cache[1][i]) &&
            (AlsM == als_cache[2][i]) && (MAls == als_cache[3][i]) &&
            (mut == als_cache[4][i]) && (mub == als_cache[5][i]) &&
            (muc == als_cache[6][i]) && (Nf_thr == (bool)als_cache[7][i]))
            return als_cache[8][i];

    switch (order)
    {
    case LO:
    case NLO:
    case NNLO:
    case NNNLO:
        if (nfAls == nfmu)
            als = AlsWithInit(mu, AlsM, MAls, nfmu, order);
        fullord = FullOrder(order);
        if (nfAls > nfmu)
        {
            mutmp = BelowTh(MAls);
            alstmp = AlsWithInit(mutmp, AlsM, MAls, nfAls, realorder);
            alstmp *= (1. - NfThresholdCorrections(mutmp, MassOfNf(nfAls), alstmp, nfAls, fullord));
            for (i = nfAls - 1; i > nfmu; i--)
            {
                mutmp = BelowTh(mutmp - MEPS);
                alstmp = AlsWithInit(mutmp, alstmp, AboveTh(mutmp) - MEPS, i, realorder);
                alstmp *= (1. - NfThresholdCorrections(mutmp, MassOfNf(i), alstmp, i, fullord));
            }
            als = AlsWithInit(mu, alstmp, AboveTh(mu) - MEPS, nfmu, order);
        }

        if (nfAls < nfmu)
        {
            mutmp = AboveTh(MAls) - MEPS;
            alstmp = AlsWithInit(mutmp, AlsM, MAls, nfAls, realorder);
            alstmp *= (1. + NfThresholdCorrections(mutmp, MassOfNf(nfAls + 1), alstmp, nfAls + 1, fullord));
            for (i = nfAls + 1; i < nfmu; i++)
            {
                mutmp = AboveTh(mutmp) - MEPS;
                alstmp = AlsWithInit(mutmp, alstmp, BelowTh(mutmp) + MEPS, i, realorder);
                alstmp *= (1. + NfThresholdCorrections(mutmp, MassOfNf(i + 1), alstmp, i + 1, fullord));
            }
            als = AlsWithInit(mu, alstmp, BelowTh(mu) + MEPS, nfmu, order);
        }

        CacheShift(als_cache, 9);
        als_cache[0][0] = mu;
        als_cache[1][0] = (double)order;
        als_cache[2][0] = AlsM;
        als_cache[3][0] = MAls;
        als_cache[4][0] = mut;
        als_cache[5][0] = mub;
        als_cache[6][0] = muc;
        als_cache[7][0] = (int)Nf_thr;
        als_cache[8][0] = als;

        return als;
    default:
        throw std::runtime_error("QCD::Als(): " + orderToString(order) + " is not implemented.");
    }
}

const double QCD::AlsOLD(const double mu, const orders order) const
{
    int i;
    for (i = 0; i < CacheSize; ++i)
        if ((mu == als_cache[0][i]) && ((double)order == als_cache[1][i]) &&
            (AlsM == als_cache[2][i]) && (MAls == als_cache[3][i]) &&
            (mut == als_cache[4][i]) && (mub == als_cache[5][i]) &&
            (muc == als_cache[6][i]) && -1 == (int)als_cache[7][i])
            return als_cache[8][i];

    double nfmu = Nf(mu), nfz = Nf(MAls), mu_thre1, mu_thre2, Als_tmp, mf;
    double als;

    switch (order)
    {
    case LO:
    case FULLNLO:
    case NLO:
        if (nfmu == nfz)
            als = AlsWithInit(mu, AlsM, MAls, nfmu, order);
        else if (nfmu > nfz)
        {
            if (order == NLO)
                throw std::runtime_error("NLO is not implemented in QCD::Als(mu,order).");
            if (nfmu == nfz + 1.)
            {
                mu_thre1 = AboveTh(MAls); // mut
                Als_tmp = AlsWithInit(mu_thre1 - MEPS, AlsM, MAls, nfz, order);
                if (order == FULLNLO)
                {
                    mf = getQuarks(TOP).getMass(); // mf = mtpole;
                    Als_tmp = (1. + Als_tmp / M_PI * log(mu_thre1 / mf) / 3.) * Als_tmp;
                }
                als = AlsWithInit(mu, Als_tmp, mu_thre1 + MEPS, nfmu, order);
            }
            else
            {
                QCDsuccess = false;
                std::cout << "Error in QCD::Als(mu,order)! QCDsuccess set to false" << std::endl;
                als = 0.;
            }
        }
        else
        {
            if (order == NLO)
                throw std::runtime_error("NLO is not implemented in QCD::Als(mu,order).");
            if (nfmu == nfz - 1.)
            {
                mu_thre1 = BelowTh(MAls); // mub
                Als_tmp = AlsWithInit(mu_thre1 + MEPS, AlsM, MAls, nfz, order);
                if (order == FULLNLO)
                {
                    mf = getQuarks(BOTTOM).getMass();
                    Als_tmp = (1. - Als_tmp / M_PI * log(mu_thre1 / mf) / 3.) * Als_tmp;
                }
                als = AlsWithInit(mu, Als_tmp, mu_thre1 - MEPS, nfmu, order);
            }
            else if (nfmu == nfz - 2.)
            {
                mu_thre1 = BelowTh(MAls); // mub
                mu_thre2 = AboveTh(mu);   // muc
                Als_tmp = Als(mu_thre1 + MEPS, order);
                if (order == FULLNLO)
                {
                    mf = getQuarks(BOTTOM).getMass();
                    Als_tmp = (1. - Als_tmp / M_PI * log(mu_thre1 / mf) / 3.) * Als_tmp;
                }
                Als_tmp = AlsWithInit(mu_thre2, Als_tmp, mu_thre1 - MEPS, nfz-1, order);
                if (order == FULLNLO)
                {
                    mf = getQuarks(CHARM).getMass();
                    Als_tmp = (1. - Als_tmp / M_PI * log(mu_thre2 / mf) / 3.) * Als_tmp;
                }
                als = AlsWithInit(mu, Als_tmp, mu_thre2 - MEPS, nfmu, order);
            }
            else
            {
                QCDsuccess = false;
                std::cout << "Error in QCD::Als(mu,order)! QCDsuccess set to false" << std::endl;
                als = 0.;
            }
        }
        break;
    case NNLO:
        //        case NNNLO:
    case FULLNNLO:
        //        case FULLNNNLO:
        /* alpha_s(mu) computed with Lambda_QCD for Nf=nfmu */
        als = AlsWithLambda(mu, order);
        break;
    default:
        throw std::runtime_error(orderToString(order) + " is not implemented in QCD::Als(mu,order).");
    }

    CacheShift(als_cache, 9);
    als_cache[0][0] = mu;
    als_cache[1][0] = (double)order;
    als_cache[2][0] = AlsM;
    als_cache[3][0] = MAls;
    als_cache[4][0] = mut;
    als_cache[5][0] = mub;
    als_cache[6][0] = muc;
    als_cache[7][0] = -1;
    als_cache[8][0] = als;

    return als;
}

const double QCD::AlsByOrder(const double mu, const int Nf_in, const orders order) const
{
    int i, nfAls = (int)Nf(MAls), nfmu = Nf_in;
    double als, alstmp, mutmp;
    orders fullord;

    for (i = 0; i < CacheSize; ++i)
        if ((mu == als_cache[0][i]) && ((double)order == als_cache[1][i]) &&
            (AlsM == als_cache[2][i]) && (MAls == als_cache[3][i]) &&
            (mut == als_cache[4][i]) && (mub == als_cache[5][i]) &&
            (muc == als_cache[6][i]) && (Nf_in == (int)als_cache[7][i]))
            return als_cache[8][i];

    switch (order)
    {
    case LO:
    case NLO:
    case NNLO:
    case NNNLO:
        if (nfAls == nfmu)
            als = AlsWithInit(mu, AlsM, MAls, nfmu, order);
        fullord = FullOrder(order);
        if (nfAls > nfmu)
        {
            mutmp = BelowTh(MAls);
            alstmp = AlsWithInit(mutmp, AlsM, MAls, nfAls, realorder);
            alstmp *= (1. - NfThresholdCorrections(mutmp, MassOfNf(nfAls), alstmp, nfAls, fullord));
            for (i = nfAls - 1; i > nfmu; i--)
            {
                mutmp = BelowTh(mutmp - MEPS);
                alstmp = AlsWithInit(mutmp, alstmp, AboveTh(mutmp) - MEPS, i, realorder);
                alstmp *= (1. - NfThresholdCorrections(mutmp, MassOfNf(i), alstmp, i, fullord));
            }
            als = AlsWithInit(mu, alstmp, mutmp, nfmu, order);
        }

        if (nfAls < nfmu)
        {
            mutmp = AboveTh(MAls) - MEPS;
            alstmp = AlsWithInit(mutmp, AlsM, MAls, nfAls, realorder);
            alstmp *= (1. + NfThresholdCorrections(mutmp, MassOfNf(nfAls + 1), alstmp, nfAls + 1, fullord));
            for (i = nfAls + 1; i < nfmu; i++)
            {
                double oldmu = mutmp;
                mutmp = AboveTh(mutmp + 2. * MEPS) - MEPS;
                alstmp = AlsWithInit(mutmp, alstmp, oldmu, i, realorder);
                alstmp *= (1. + NfThresholdCorrections(mutmp, MassOfNf(i + 1), alstmp, i + 1, fullord));
            }
            als = AlsWithInit(mu, alstmp, mutmp + 2. * MEPS, nfmu, order);
        }

        CacheShift(als_cache, 9);
        als_cache[0][0] = mu;
        als_cache[1][0] = (double)order;
        als_cache[2][0] = AlsM;
        als_cache[3][0] = MAls;
        als_cache[4][0] = mut;
        als_cache[5][0] = mub;
        als_cache[6][0] = muc;
        als_cache[7][0] = Nf_in;
        als_cache[8][0] = als;

        return als;
    default:
        throw std::runtime_error("QCD::Als(): " + orderToString(order) + " is not implemented.");
    }
}

const double QCD::ZeroNf6NLO(double *logLambda6, double *logLambda5_in) const
{
    return (AlsWithLambda(mut + 1.e-10, *logLambda6, FULLNLO) - AlsWithLambda(mut - 1.e-10, *logLambda5_in, FULLNLO));
}

const double QCD::ZeroNf5(double *logLambda5, double *order) const
{
    return (AlsWithLambda(MAls, *logLambda5, (orders)*order) - AlsM);
}

const double QCD::ZeroNf4NLO(double *logLambda4, double *logLambda5_in) const
{
    return (AlsWithLambda(mub - 1.e-10, *logLambda4, FULLNLO) - AlsWithLambda(mub + 1.e-10, *logLambda5_in, FULLNLO));
}

const double QCD::ZeroNf3NLO(double *logLambda3, double *logLambda4_in) const
{
    return (AlsWithLambda(muc - 1.e-10, *logLambda3, FULLNLO) - AlsWithLambda(muc + 1.e-10, *logLambda4_in, FULLNLO));
}

const double QCD::logLambda5(orders order) const
{
    if (order == NLO)
        order = FULLNLO;
    if (order == NNLO)
        order = FULLNNLO;

    for (int i = 0; i < CacheSize; ++i)
        if ((AlsM == logLambda5_cache[0][i]) && (MAls == logLambda5_cache[1][i]) && ((double)order == logLambda5_cache[2][i]))
            return logLambda5_cache[3][i];

    CacheShift(logLambda5_cache, 4);
    logLambda5_cache[0][0] = AlsM;
    logLambda5_cache[1][0] = MAls;
    logLambda5_cache[2][0] = (double)order;

    if (order == LO)
        logLambda5_cache[3][0] = log(MAls) - 2. * M_PI / Beta0(5.) / AlsM;
    else
    {
        double xmin = -4., xmax = -0.2;
        TF1 f = TF1("f", this, &QCD::ZeroNf5, xmin, xmax, 1, "QCD", "zeroNf5");

        ROOT::Math::WrappedTF1 wf1(f);
        double ledouble = (double)order;
        wf1.SetParameters(&ledouble);

        ROOT::Math::BrentRootFinder brf;
        brf.SetFunction(wf1, xmin, xmax);

        if (brf.Solve())
            logLambda5_cache[3][0] = brf.Root();
        else
        {
            QCDsuccess = false;
            std::cout << "Error in QCD::logLambda5()! QCDsuccess set to false" << std::endl;
        }
    }
    return (logLambda5_cache[3][0]);
}

const double QCD::logLambdaNLO(const double nfNEW, const double nfORG,
                               const double logLambdaORG) const
{
    for (int i = 0; i < CacheSize; ++i)
        if ((AlsM == logLambdaNLO_cache[0][i]) && (MAls == logLambdaNLO_cache[1][i]) && (mut == logLambdaNLO_cache[2][i]) && (mub == logLambdaNLO_cache[3][i]) && (muc == logLambdaNLO_cache[4][i]) && (nfNEW == logLambdaNLO_cache[5][i]) && (nfORG == logLambdaNLO_cache[6][i]) && (logLambdaORG == logLambdaNLO_cache[7][i]))
            return logLambdaNLO_cache[8][i];

    CacheShift(logLambdaNLO_cache, 9);
    logLambdaNLO_cache[0][0] = AlsM;
    logLambdaNLO_cache[1][0] = MAls;
    logLambdaNLO_cache[2][0] = mut;
    logLambdaNLO_cache[3][0] = mub;
    logLambdaNLO_cache[4][0] = muc;
    logLambdaNLO_cache[5][0] = nfNEW;
    logLambdaNLO_cache[6][0] = nfORG;
    logLambdaNLO_cache[7][0] = logLambdaORG;

    double xmin = -4., xmax = -0.2;

    TF1 f;
    if (nfNEW == 6. && nfORG == 5.)
    {
        f = TF1("f", this, &QCD::ZeroNf6NLO, xmin, xmax, 1, "QCD", "zeroNf6NLO");
    }
    else if (nfNEW == 4. && nfORG == 5.)
    {
        f = TF1("f", this, &QCD::ZeroNf4NLO, xmin, xmax, 1, "QCD", "zeroNf4NLO");
    }
    else if (nfNEW == 3. && nfORG == 4.)
    {
        f = TF1("f", this, &QCD::ZeroNf3NLO, xmin, xmax, 1, "QCD", "zeroNf3NLO");
    }
    else
    {
        QCDsuccess = false;
        std::cout << "Error in QCD::logLambdaNLO()! QCDsuccess set to false" << std::endl;
    }

    ROOT::Math::WrappedTF1 wf1(f);
    wf1.SetParameters(&logLambdaORG);

    ROOT::Math::BrentRootFinder brf;
    brf.SetFunction(wf1, xmin, xmax);

    if (brf.Solve())
        logLambdaNLO_cache[8][0] = brf.Root();
    else
    {
        QCDsuccess = false;
        std::cout << "Error in QCD::logLambdaNLO()! QCDsuccess set to false" << std::endl;
    }

    return (logLambdaNLO_cache[8][0]);
}

const double QCD::logLambda(const double muMatching, const double mf,
                            const double nfNEW, const double nfORG,
                            const double logLambdaORG, orders order) const
{
    if (fabs(nfNEW - nfORG) != 1.)
    {
        QCDsuccess = false;
        std::cout << "Error in QCD::logLambda()! QCDsuccess set to false" << std::endl;
    }

    if (order == NLO)
        order = FULLNLO;
    if (order == NNLO)
        order = FULLNNLO;

    /* We do not use the codes below for FULLNLO, since threshold corrections
     * can be regarded as an NNLO effect as long as setting the matching scale
     * to be close to the mass scale of the decoupling quark. In order to use
     * the relation als^{nf+1} = als^{nf} exactly, we use logLambdaNLO method.
     */
    if (order == FULLNLO)
        return logLambdaNLO(nfNEW, nfORG, logLambdaORG);

    double logMuMatching = log(muMatching);
    double L = 2. * (logMuMatching - logLambdaORG);
    double rNEW = 0.0, rORG = 0.0, log_mu2_mf2 = 0.0, log_L = 0.0;
    double C1 = 0.0, C2 = 0.0; // threshold corrections
    double logLambdaNEW;

    // LO contribution
    logLambdaNEW = 1. / 2. / Beta0(nfNEW) * (Beta0(nfNEW) - Beta0(nfORG)) * L + logLambdaORG;

    // NLO contribution
    if (order == FULLNLO || order == FULLNNLO)
    {
        rNEW = Beta1(nfNEW) / Beta0(nfNEW);
        rORG = Beta1(nfORG) / Beta0(nfORG);
        log_mu2_mf2 = 2. * (logMuMatching - log(mf));
        log_L = log(L);
        if (nfNEW < nfORG)
            C1 = 2. / 3. * log_mu2_mf2;
        else
            C1 = -2. / 3. * log_mu2_mf2;
        logLambdaNEW += 1. / 2. / Beta0(nfNEW) * ((rNEW - rORG) * log_L - rNEW * log(Beta0(nfNEW) / Beta0(nfORG)) - C1);
    }

    // NNLO contribution
    if (order == FULLNNLO)
    {
        if (nfNEW == 5. && nfORG == 6.)
            C2 = -16. * (log_mu2_mf2 * log_mu2_mf2 / 36. - 19. / 24. * log_mu2_mf2 - 7. / 24.);
        else if (nfNEW == 6. && nfORG == 5.)
            C2 = -16. * (log_mu2_mf2 * log_mu2_mf2 / 36. + 19. / 24. * log_mu2_mf2 + 7. / 24.);
        else
        {
            if (nfNEW < nfORG)
                C2 = -16. * (log_mu2_mf2 * log_mu2_mf2 / 36. - 19. / 24. * log_mu2_mf2 + 11. / 72.);
            else
                C2 = -16. * (log_mu2_mf2 * log_mu2_mf2 / 36. + 19. / 24. * log_mu2_mf2 - 11. / 72.);
        }
        logLambdaNEW += 1. / 2. / Beta0(nfNEW) / Beta0(nfORG) / L * (rORG * (rNEW - rORG) * log_L + rNEW * rNEW - rORG * rORG - Beta2(nfNEW) / Beta0(nfNEW) + Beta2(nfORG) / Beta0(nfORG) + rNEW * C1 - C1 * C1 - C2);
    }

    return logLambdaNEW;
}

const double QCD::logLambda(const double nf, orders order) const
{
    if (order == NLO)
        order = FULLNLO;
    if (order == NNLO)
        order = FULLNNLO;

    double muMatching, mf, logLambdaORG, logLambdaNEW;
    if (nf == 5.)
        return logLambda5(order);
    else if (nf == 6.)
    {
        muMatching = Thresholds(1); // mut
        /* matching condition from Nf=5 to Nf=6 is given in terms of the top pole mass. */
        mf = mtpole; // top pole mass
        return logLambda(muMatching, mf, 6., 5., logLambda5(order), order);
    }
    else if (nf == 4. || nf == 3.)
    {
        muMatching = Thresholds(2);       // mub
        mf = getQuarks(BOTTOM).getMass(); // m_b(m_b)
        logLambdaORG = logLambda5(order);
        logLambdaNEW = logLambda(muMatching, mf, 4., 5., logLambdaORG, order);
        if (nf == 3.)
        {
            muMatching = Thresholds(3);      // muc
            mf = getQuarks(CHARM).getMass(); // m_c(m_c)
            logLambdaORG = logLambdaNEW;
            logLambdaNEW = logLambda(muMatching, mf, 3., 4., logLambdaORG, order);
        }
        return logLambdaNEW;
    }
    else
    {
        QCDsuccess = false;
        std::cout << "Error in QCD::logLambda()! QCDsuccess set to false" << std::endl;
    }
    return 0.;
}

////////////////////////////////////////////////////////////////////////

const double QCD::Gamma0(const double nf) const
{
    return (6. * CF);
}

const double QCD::Gamma1(const double nf) const
{
    return (CF * (3. * CF + 97. / 3. * Nc - 10. / 3. * nf));
}

const double QCD::Gamma2(const double nf) const
{
    return (129. * CF * CF * CF - 129. / 2. * CF * CF * Nc + 11413. / 54. * CF * Nc * Nc + CF * CF * nf * (-46. + 48. * zeta3) + CF * Nc * nf * (-556. / 27. - 48. * zeta3) - 70. / 27. * CF * nf * nf);
}

const double QCD::threCorrForMass(const double nf_f, const double nf_i) const
{
    if (fabs(nf_f - nf_i) != 1.)
    {
        QCDsuccess = false;
        std::cout << "Error in QCD::threCorrForMass()! QCDsuccess set to false" << std::endl;
        return 0.;
    }

    // hep-ph/9712201v2 eq. (48)
    double mu_threshold, mf, log_mu2_mf2, AlsoverPi;
    if (nf_f > nf_i)
    {
        if (nf_f == 6.)
        {
            mu_threshold = mut;
            mf = quarks[TOP].getMass(); // m_t(m_t)
        }
        else if (nf_f == 5.)
        {
            mu_threshold = mub;
            mf = quarks[BOTTOM].getMass(); // m_b(m_b)
        }
        else if (nf_f == 4.)
        {
            mu_threshold = muc;
            mf = quarks[CHARM].getMass(); // m_c(m_c)
        }
        else
        {
            QCDsuccess = false;
            std::cout << "Error in QCD::threCorrForMass()! QCDsuccess set to false" << std::endl;
            return 0.;
        }
        log_mu2_mf2 = 2. * log(mu_threshold / mf);
        AlsoverPi = Als(mu_threshold - MEPS, FULLNNLO) / M_PI;
        return (1. + AlsoverPi * AlsoverPi * (-log_mu2_mf2 * log_mu2_mf2 / 12. + 5. / 36. * log_mu2_mf2 - 89. / 432.));
    }
    else
    {
        if (nf_i == 6.)
        {
            mu_threshold = mut;
            mf = quarks[TOP].getMass(); // m_t(m_t)
        }
        else if (nf_i == 5.)
        {
            mu_threshold = mub;
            mf = quarks[BOTTOM].getMass(); // m_b(m_b)
        }
        else if (nf_i == 4.)
        {
            mu_threshold = muc;
            mf = quarks[CHARM].getMass(); // m_c(m_c)
        }
        else
        {
            QCDsuccess = false;
            std::cout << "Error in QCD::threCorrForMass()! QCDsuccess set to false" << std::endl;
            return 0.;
        }
        log_mu2_mf2 = 2. * log(mu_threshold / mf);
        AlsoverPi = Als(mu_threshold + MEPS, FULLNNLO) / M_PI;
        return (1. + AlsoverPi * AlsoverPi * (log_mu2_mf2 * log_mu2_mf2 / 12. - 5. / 36. * log_mu2_mf2 + 89. / 432.));
    }
}

const double QCD::Mrun(const double mu, const double m, const quark q, const orders order) const
{
    return Mrun(mu, m, m, q, order);
}

const double QCD::Mrun(const double mu_f, const double mu_i, const double m, const quark q,
                       const orders order) const
{
    // Note: When the scale evolves across a flavour threshold, the definitions
    //       of the outputs for "NLO" and "NNLO" become complicated.

    int i;
    for (i = 0; i < CacheSize; ++i)
    {
        if ((mu_f == mrun_cache[0][i]) && (mu_i == mrun_cache[1][i]) &&
            (m == mrun_cache[2][i]) && ((double)order == mrun_cache[3][i]) &&
            ((double)q == mrun_cache[4][i]) && 
            (AlsM == mrun_cache[5][i]) && (MAls == mrun_cache[6][i]) &&
            (mut == mrun_cache[7][i]) && (mub == mrun_cache[8][i]) &&
            (muc == mrun_cache[9][i]))
            return mrun_cache[10][i];
    }

    double nfi = Nf(mu_i), nff = Nf(mu_f);
    double nfq;
    if (q == TOP)
        nfq = 6.;
    else if (q == BOTTOM)
        nfq = 5.;
    else if (q == CHARM)
        nfq = 4.;
    else 
        nfq = 3.;

    nfi = (nfi > nfq) ? nfi : nfq;
    nff = (nff > nfq) ? nff : nfq;

    double mu_threshold, mu_threshold2, mu_threshold3, mrun;
    if (nff == nfi)
        mrun = MrunTMP(mu_f, mu_i, m, lround(nfi), order);
    else if (nff > nfi)
    {
        if (order == NLO || order == NNLO)
            throw std::runtime_error(orderToString(order) + " is not implemented in QCD::Mrun(mu_f,mu_i,m,q,order)");
        mu_threshold = AboveTh(mu_i);
        mrun = MrunTMP(mu_threshold - MEPS, mu_i, m, lround(nfi), order);
        if (order == FULLNNLO)
            mrun *= threCorrForMass(nfi + 1., nfi); // threshold corrections
        if (nff == nfi + 1.)
        {
            mrun = MrunTMP(mu_f, mu_threshold + MEPS, mrun, lround(nff), order);
        }
        else if (nff == nfi + 2.)
        {
            mu_threshold2 = BelowTh(mu_f);
            mrun = MrunTMP(mu_threshold2 - MEPS, mu_threshold + MEPS, mrun, lround(nfi) + 1, order);
            if (order == FULLNNLO)
                mrun *= threCorrForMass(nff, nfi + 1.); // threshold corrections
            mrun = MrunTMP(mu_f, mu_threshold2 + MEPS, mrun, lround(nff), order);
        }
        else if (nff == nfi + 3.)
        {
            mu_threshold2 = AboveTh(mu_threshold);
            mrun = MrunTMP(mu_threshold2 - MEPS, mu_threshold + MEPS, mrun, lround(nfi) + 1, order);
            if (order == FULLNNLO)
                mrun *= threCorrForMass(nfi + 2., nfi + 1.); // threshold corrections
            mu_threshold3 = BelowTh(mu_f);
            mrun = MrunTMP(mu_threshold3 - MEPS, mu_threshold2 + MEPS, mrun, lround(nff) - 1, order);
            if (order == FULLNNLO)
                mrun *= threCorrForMass(nff, nfi + 2.); // threshold corrections
            mrun = MrunTMP(mu_f, mu_threshold3 + MEPS, mrun, lround(nff), order);
        }
        else
        {
            QCDsuccess = false;
            std::cout << "Error in QCD::Mrun(mu_f,mu_i,m,order)! QCDsuccess set to false" << std::endl;
            mrun = 0.;
        }
    }
    else
    {
        if (order == NLO || order == NNLO)
            throw std::runtime_error(orderToString(order) + " is not implemented in QCD::Mrun(mu_f,mu_i,m,order)");
        mu_threshold = BelowTh(mu_i);
        mrun = MrunTMP(mu_threshold + MEPS, mu_i, m, lround(nfi), order);
        if (order == FULLNNLO)
            mrun *= threCorrForMass(nfi - 1., nfi); // threshold corrections
        if (nff == nfi - 1.)
            mrun = MrunTMP(mu_f, mu_threshold - MEPS, mrun, lround(nff), order);
        else if (nff == nfi - 2.)
        {
            mu_threshold2 = AboveTh(mu_f);
            mrun = MrunTMP(mu_threshold2 + MEPS, mu_threshold - MEPS, mrun, lround(nfi) - 1, order);
            if (order == FULLNNLO)
                mrun *= threCorrForMass(nff, nfi - 1.); // threshold corrections
            mrun = MrunTMP(mu_f, mu_threshold2 - MEPS, mrun, lround(nff), order);
        }
        else
        {
            QCDsuccess = false;
            std::cout << "Error in QCD::Mrun(mu_f,mu_i,m,order)! QCDsuccess set to false" << std::endl;
            mrun = 0.;
        }
    }

    if (mrun < 0.0)
    {
        std::stringstream out;
        out << "QCD::Mrun(): A quark mass becomes tachyonic in QCD::Mrun("
            << mu_f << ", " << mu_i << ", " << m << ", " << orderToString(order) << ")"
            << std::endl
            << "             Als(" << mu_i << ", " << orderToString(order) << ")/(4pi)="
            << Als(mu_i, order) / (4. * M_PI) << std::endl
            << "             Als(" << mu_f << ", " << orderToString(order) << ")/(4pi)="
            << Als(mu_f, order) / (4. * M_PI);
        QCDsuccess = false;
        std::cout << "Error in QCD::Mrun(mu_f,mu_i,m,order)! QCDsuccess set to false" << std::endl;
        mrun = 0.;
    }

    CacheShift(mrun_cache, 11);
    mrun_cache[0][0] = mu_f;
    mrun_cache[1][0] = mu_i;
    mrun_cache[2][0] = m;
    mrun_cache[3][0] = (double)order;
    mrun_cache[4][0] = (double)q;
    mrun_cache[5][0] = AlsM;
    mrun_cache[6][0] = MAls;
    mrun_cache[7][0] = mut;
    mrun_cache[8][0] = mub;
    mrun_cache[9][0] = muc;
    mrun_cache[10][0] = mrun;

    return mrun;
}

const double QCD::MrunTMP(const double mu_f, const double mu_i, const double m, const int nf,
                          const orders order) const
{

    // alpha_s/(4pi)
    orders orderForAls;
    if (order == LO)
        orderForAls = LO;
    if (order == NLO || order == FULLNLO)
        orderForAls = FULLNLO;
    if (order == NNLO || order == FULLNNLO)
        orderForAls = FULLNNLO;
    double ai = Als(mu_i, nf, orderForAls) / (4. * M_PI);
    double af = Als(mu_f, nf, orderForAls) / (4. * M_PI);

    // LO contribution
    double b0 = Beta0((double) nf), g0 = Gamma0((double) nf);
    double mLO = m * pow(af / ai, g0 / (2. * b0));
    if (order == LO)
        return mLO;

    // NLO contribution
    double b1 = Beta1((double) nf), g1 = Gamma1((double) nf);
    double A1 = g1 / (2. * b0) - b1 * g0 / (2. * b0 * b0);
    double mNLO = mLO * A1 * (af - ai);
    if (order == NLO)
        return mNLO;
    if (order == FULLNLO)
        return (mLO + mNLO);

    // NNLO contribution: Chetyrkin, hep-ph/9703278v1  eq. (15)
    double b2 = Beta2((double) nf), g2 = Gamma2((double) nf);
    double A2 = b1 * b1 * g0 / (2. * b0 * b0 * b0) - b2 * g0 / (2. * b0 * b0) - b1 * g1 / (2. * b0 * b0) + g2 / (2. * b0);
    double mNNLO = mLO * (A1 * A1 / 2. * (af - ai) * (af - ai) + A2 / 2. * (af * af - ai * ai));
    if (order == NNLO)
        return mNNLO;
    if (order == FULLNNLO)
        return (mLO + mNLO + mNNLO);

    throw std::runtime_error(orderToString(order) + " is not implemented in QCD::MrunTMP()");
}

const double QCD::Mrun4(const double mu_f, const double mu_i, const double m) const
{
    double nf = 4.;

    // alpha_s/(4pi)
    double ai = Als4(mu_i) / (4. * M_PI);
    double af = Als4(mu_f) / (4. * M_PI);

    // LO contribution
    double b0 = Beta0(nf), g0 = Gamma0(nf);
    double mLO = m * pow(af / ai, g0 / (2. * b0));

    // NLO contribution
    double b1 = Beta1(nf), g1 = Gamma1(nf);
    double A1 = g1 / (2. * b0) - b1 * g0 / (2. * b0 * b0);
    double mNLO = mLO * A1 * (af - ai);
    return (mLO + mNLO);
}

////////////////////////////////////////////////////////////////////////

const double QCD::Mbar2Mp(const double mbar, const quark q, const orders order) const
{
    int Nf_m;
    double nl;
    std::vector<double> x = {0.0};

    if (q == CHARM)
        Nf_m = 4;
    else if (q == BOTTOM)
        Nf_m = 5; 
    else if (q == TOP)
        Nf_m = 6;
    else
        throw std::runtime_error("QCD::Mbar2Mp() can convert only top, bottom and charm masses");

    nl = Nf_m-1;

    // LO contribution
    double MpLO = mbar;
    if (order == LO)
        return MpLO;

    // alpha_s(mbar)/pi
    orders orderForAls;
    if (order == NLO || order == FULLNLO)
        orderForAls = FULLNLO;
    if (order == NNLO || order == FULLNNLO)
        orderForAls = FULLNNLO;
    double a = Als(mbar, Nf_m, orderForAls) / M_PI;

    // NLO contribution
    double MpNLO = mbar * 4. / 3. * a;
    if (order == NLO)
        return MpNLO;
    if (order == FULLNLO)
        return (MpLO + MpNLO);

    // NNLO contribution

    x.push_back(quarks[STRANGE].getMass() / mbar); 
    if (nl >3)
        x.push_back(quarks[CHARM].getMass() / mbar);
    if (nl > 4)
        x.push_back(quarks[BOTTOM].getMass() / mbar);
    
    double Delta = 0.;
    for (int i = 0; i < x.size(); i++)
        Delta += DeltaMass(x[i]);

    double MpNNLO = mbar * (307. / 32. + 2. * zeta2 + 2. / 3. * zeta2 * log(2.0) - zeta3 / 6. - nl / 3. * (zeta2 + 71. / 48.) + 4. / 3. * Delta) * a * a;

    if (order == NNLO)
        return MpNNLO;
    if (order == FULLNNLO)
        return (MpLO + MpNLO + MpNNLO);

    throw std::runtime_error(orderToString(order) + " is not implemented in QCD::Mbar2Mp().");
}

double QCD::DeltaMass(double x) const
{
    return M_PI * M_PI / 8. * x - 0.597 * x * x + 0.230 * x * x * x;
}


const double QCD::Mp2MbarTMP(double *mu, double *params) const
{
    double mp = params[0];
    quark q = (quark)params[1];
    orders order = (orders)params[2];
    return (mp - Mbar2Mp(*mu, q, order));
}

const double QCD::Mp2Mbar_bar(const double mp, const quark q, const orders order) const
{
    if (order == NLO || order == NNLO)
        throw std::runtime_error(orderToString(order) + " is not implemented in QCD::Mp2Mbar().");

    int i;
    double ms = quarks[STRANGE].getMass(), mc = quarks[CHARM].getMass(),
        mb = quarks[BOTTOM].getMass();
    double alsmp = Als(mp, order);
    for (i = 0; i < CacheSize; ++i)
        if (alsmp == mp2mbar_cache[0][i] && ms == mp2mbar_cache[1][i] &&
            mc == mp2mbar_cache[2][i] && mb == mp2mbar_cache[3][i]
            && (double)order == mp2mbar_cache[4][i])
            return mp2mbar_cache[5][i];

    TF1 f("f", this, &QCD::Mp2MbarTMP, mp / 2., 2. * mp, 3, "QCD", "Mp2MbarTMP");

    ROOT::Math::WrappedTF1 wf1(f);
    double params[3];
    params[0] = mp;
    params[1] = (double)q; 
    params[2] = (double)order;
    wf1.SetParameters(params);

    ROOT::Math::BrentRootFinder brf;

    brf.SetNpx(5);
    brf.SetFunction(wf1, .5 * mp, 2. * mp);
    if (brf.Solve())
    {
        CacheShift(mp2mbar_cache, 6);
        mp2mbar_cache[0][0] = alsmp;
        mp2mbar_cache[1][0] = ms;
        mp2mbar_cache[2][0] = mc;
        mp2mbar_cache[3][0] = mb;
        mp2mbar_cache[4][0] = (double)order;
        mp2mbar_cache[5][0] = brf.Root();
    }
    else
    {
        QCDsuccess = false;
        std::cout << "Error in QCD::Mp2Mbar()! QCDsuccess set to false" << std::endl;
    }

    return (mp2mbar_cache[5][0]);
}

const double QCD::Mp2Mbar(const double mp, const quark q, const orders order) const
{
    if (q == TOP)
//        return Mp2Mbar_pole(mp, order);
        return Mp2Mbar_bar(mp, q, order);
    else if (q == BOTTOM || q == CHARM)
        return Mp2Mbar_bar(mp, q, order);
    else
        throw std::runtime_error("QCD::Mp2Mbar() can convert only top, bottom and charm masses");
}

const double QCD::Mp2Mbar_pole(const double mp, const orders order) const
{
    int Nfmp = 6;
        // LO contribution
        double mbarLO = mp;
        if (order == LO)
            return mbarLO;

        // alpha_s(mp)/pi
        orders orderForAls;
        if (order == NLO || order == FULLNLO)
            orderForAls = FULLNLO;
        if (order == NNLO || order == FULLNNLO)
            orderForAls = FULLNNLO;
        double a = Als(mp, Nfmp, orderForAls) / M_PI;

        // NLO contribution
        double mbarNLO = mp * (-4. / 3.) * a;
        if (order == NLO)
            return mbarNLO;
        if (order == FULLNLO)
            return (mbarLO + mbarNLO);

        double mbarNNLO = mp * ((-1111. * CA * CF) / 384. + (199. * CF * CF) / 128. + (143. * CF * TF) / 96. +
                                    (71. * CF * ((double)Nfmp - 1.) * TF) / 96. + (CA * CF * zeta2) / 2. - (15. * CF * CF * zeta2) / 8. -
                                    (3. * CA * CF * log(2.) * zeta2) / 2. + 3. * CF * CF * log(2.) * zeta2 - CF * TF * zeta2 +
                                    (CF * ((double)Nfmp - 1.) * TF * zeta2) / 2. + (3. * CA * CF * zeta3) / 8. - (3. * CF * CF * zeta3) / 4.) * a * a;
        
        //for the top, the correction from light quark masses is negligible
        if (order == NNLO)
            return mbarNNLO;
        if (order == FULLNNLO)
            return (mbarLO + mbarNLO + mbarNNLO);

    throw std::runtime_error(orderToString(order) + " is not implemented in QCD::Mp2Mbar_new().");
}

const double QCD::MS2DRqmass(const double MSbar) const
{
    return (MSbar / (1. + Als(MSbar, FULLNLO) / 4. / M_PI * CF));
}

const double QCD::MS2DRqmass(const double MSscale, const double MSbar) const
{
    return (MSbar / (1. + Als(MSscale, FULLNLO) / 4. / M_PI * CF));
}

const double QCD::Mofmu2MbarTMP(double *mu, double *params) const
{
    double mofmu = params[0];
    double muI = params[1];
    quark q = (quark)params[2];
    return (*mu - Mrun(*mu, muI, mofmu, q));
}

const double QCD::Mofmu2Mbar(const double m, const double mu, const quark q) const
{

    // First move to the right region

    double mutmp=0.;
    switch (q)
    {
    case TOP:
        mutmp = 165.;
        break;
    case BOTTOM:
        mutmp = 4.2;
        break;
    case CHARM:
        mutmp = 2.;
        break;
    default:
        QCDsuccess = false;
        std::cout << "Error in QCD::Mofmu2Mbar()! QCDsuccess set to false" << std::endl;
        break;
    }
    double mlow = Mrun(mutmp, mu, m, q);
    TF1 f("f", this, &QCD::Mofmu2MbarTMP, mlow / 2., 2. * mlow, 3, "QCD", "mofmu2mbara");

    ROOT::Math::WrappedTF1 wf1(f);

    double params[3];
    params[0] = mlow;
    params[1] = mutmp;
    params[2] = (double)q;
    wf1.SetParameters(params);

    ROOT::Math::BrentRootFinder brf;
    brf.SetNpx(5);
    brf.SetFunction(wf1, .3 * mlow, 3. * mlow);
    if (brf.Solve())
        mlow = brf.Root();
    else
        {
            QCDsuccess = false;
            std::cout << "Error in QCD::Mofmu2Mbar()! QCDsuccess set to false" << std::endl;
            mlow = 0.;
        }

    return (mlow);
}
