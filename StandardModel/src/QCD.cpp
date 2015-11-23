/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
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
#include <gsl/gsl_sf.h>
#include <TF1.h>
#include <Math/WrappedTF1.h>
#include <Math/BrentRootFinder.h>
#include "QCD.h"

const std::string QCD::QCDvars[NQCDvars] = {
    "AlsM", "MAls",
    "mup", "mdown", "mcharm", "mstrange", "mtop", "mbottom",
    "muc", "mub", "mut",
    "MK0", "MKp", "MD", "MBd", "MBp", "MBs", "MKstar", "Mphi",
    "tKl", "tKp", "tBd", "tBp", "tBs", "tKstar", "tphi",
    "DGs_Gs",
    "FK", "FD", "FBs", "FBsoFBd", "FKstar", "FKstarp", "Fphi",
    "BK1", "BK2", "BK3", "BK4", "BK5", "BKscale", "BKscheme",
    "BD1", "BD2", "BD3", "BD4", "BD5", "BDscale", "BDscheme",
    "BBsoBBd",
    "BBs1", "BBs2", "BBs3", "BBs4", "BBs5", "BBsscale", "BBsscheme",
    "BK(1/2)1", "BK(1/2)2", "BK(1/2)3", "BK(1/2)4", "BK(1/2)5",
    "BK(1/2)6", "BK(1/2)7", "BK(1/2)8", "BK(1/2)9", "BK(1/2)10",
    "BK(3/2)1", "BK(3/2)2", "BK(3/2)3", "BK(3/2)4", "BK(3/2)5",
    "BK(3/2)6", "BK(3/2)7", "BK(3/2)8", "BK(3/2)9", "BK(3/2)10",
    "BKd_scale", "BKd_scheme",
    "ReA2_Kd", "ReA0_Kd", "Omega_eta_etap",
    "Br_Kp_P0enu", "Br_Kp_munu", "Br_B_Xcenu", "DeltaP_cu", "IB_Kl", "IB_Kp",
    "a_0V", "a_1V", "a_2V", "MRV", "a_0A0", "a_1A0", "a_2A0", "MRA0", "a_0A1", "a_1A1", "a_2A1", "MRA1", "a_0A12", "a_1A12", "a_2A12", "MRA12",
    "a_0T1", "a_1T1", "a_2T1", "MRT1", "a_0T2", "a_1T2", "a_2T2", "MRT2", "a_0T23", "a_1T23", "a_2T23", "MRT23",
    "a_0Vphi", "a_1Vphi", "a_2Vphi", "MRVphi", "a_0A0phi", "a_1A0phi", "a_2A0phi", "MRA0phi", "a_0A1phi", "a_1A1phi", "a_2A1phi", "MRA1phi", "a_0A12phi", "a_1A12phi", "a_2A12phi",
    "MRA12phi", "a_0T1phi", "a_1T1phi", "a_2T1phi", "MRT1phi", "a_0T2phi", "a_1T2phi", "a_2T2phi", "MRT2phi", "a_0T23phi", "a_1T23phi", "a_2T23phi", "MRT23phi",
    "reh_0", "reh_p", "reh_m", "imh_0", "imh_p", "imh_m",
    "reh_0_1", "reh_p_1", "reh_m_1", "imh_0_1", "imh_p_1", "imh_m_1",
    "reh_0_2", "reh_p_2", "reh_m_2", "imh_0_2", "imh_p_2", "imh_m_2",
    "r_1_fplus", "r_2_fplus", "m_fit2_fplus", "r_1_fT", "r_2_fT", "m_fit2_fT", "r_2_f0", "m_fit2_f0",
    "reh_0_MP", "imh_0_MP", "reh_0_1_MP", "imh_0_1_MP",
    "bsgamma_E0", "BLNPcorr", "Gambino_mukin", "Gambino_BRsem", "Gambino_Mbkin", "Gambino_Mcatmuc", "Gambino_mupi2", "Gambino_rhoD3", "Gambino_muG2", "Gambino_rhoLS3",
    "lambdaB", "alpha1kst", "alpha2kst"
    //"r_2A0", "r_2T1", "r_2T2", "r_2A0phi", "r_2T1phi", "r_2T2phi"
};

QCD::QCD()
: BBs(5), BBd(5), BD(5), BK(5), BKd1(10), BKd3(10)
{
    Nc = 3.;
    CF = Nc / 2. - 1. / (2. * Nc);
    //    Particle(std::string name, double mass, double mass_scale = 0., double width = 0., double charge = 0.,double isospin = 0.);
    quarks[UP] = Particle("UP", 0., 2., 0., 2. / 3., .5);
    quarks[CHARM] = Particle("CHARM", 0., 0., 0., 2. / 3., .5);
    quarks[TOP] = Particle("TOP", 0., 0., 0., 2. / 3., .5);
    quarks[DOWN] = Particle("DOWN", 0., 2., 0., -1. / 3., -.5);
    quarks[STRANGE] = Particle("STRANGE", 0., 2., 0., -1. / 3., -.5);
    quarks[BOTTOM] = Particle("BOTTOM", 0., 0., 0., -1. / 3., -.5);
    zeta2 = gsl_sf_zeta_int(2);
    zeta3 = gsl_sf_zeta_int(3);
    for (int i = 0; i < CacheSize; i++) {
        for (int j = 0; j < 8; j++)
            als_cache[j][i] = 0.;
        for (int j = 0; j < 4; j++)
            logLambda5_cache[j][i] = 0.;
        for (int j = 0; j < 10; j++)
            mrun_cache[j][i] = 0.;
        for (int j = 0; j < 5; j++)
            mp2mbar_cache[j][i] = 0.;
    }

    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("AlsM", boost::cref(AlsM)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MAls", boost::cref(MAls)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mup", boost::cref(quarks[UP].getMass())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mdown", boost::cref(quarks[DOWN].getMass())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mcharm", boost::cref(quarks[CHARM].getMass())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mstrange", boost::cref(quarks[STRANGE].getMass())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mtop", boost::cref(mtpole)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mbottom", boost::cref(quarks[BOTTOM].getMass())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("muc", boost::cref(muc)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mub", boost::cref(mub)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mut", boost::cref(mut)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MK0", boost::cref(mesons[K_0].getMass())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MKp", boost::cref(mesons[K_P].getMass())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MD", boost::cref(mesons[D_0].getMass())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MBd", boost::cref(mesons[B_D].getMass())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MBp", boost::cref(mesons[B_P].getMass())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MBs", boost::cref(mesons[B_S].getMass())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MKstar", boost::cref(mesons[K_star].getMass())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Mphi", boost::cref(mesons[PHI].getMass())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("tKl", boost::cref(mesons[K_0].getWidth())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("tKp", boost::cref(mesons[K_P].getWidth())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("tBd", boost::cref(mesons[B_D].getWidth())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("tBp", boost::cref(mesons[B_P].getWidth())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("tBs", boost::cref(mesons[B_S].getWidth())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("DGs_Gs", boost::cref(mesons[B_S].getDgamma_gamma())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("tKstar", boost::cref(mesons[K_star].getWidth())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("tphi", boost::cref(mesons[PHI].getWidth())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("FK", boost::cref(mesons[K_0].getDecayconst())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("FD", boost::cref(mesons[D_0].getDecayconst())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("FBs", boost::cref(mesons[B_S].getDecayconst())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("FKstar", boost::cref(mesons[K_star].getDecayconst())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("FKstarp", boost::cref(FKstarp)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Fphi", boost::cref(mesons[PHI].getDecayconst())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("FBsoFBd", boost::cref(FBsoFBd)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BK1", boost::cref(BKB0)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BK2", boost::cref(BKB1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BK3", boost::cref(BKB2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BK4", boost::cref(BKB3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BK5", boost::cref(BKB4)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BKscale", boost::cref(BKscale)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BD1", boost::cref(BDB0)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BD2", boost::cref(BDB1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BD3", boost::cref(BDB2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BD4", boost::cref(BDB3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BD5", boost::cref(BDB4)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BDscale", boost::cref(BDscale)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BBsoBBd", boost::cref(BBsoBBd)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BBs1", boost::cref(BBsB0)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BBs2", boost::cref(BBsB1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BBs3", boost::cref(BBsB2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BBs4", boost::cref(BBsB3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BBs5", boost::cref(BBsB4)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BBsscale", boost::cref(BBsscale)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BK(1/2)1", boost::cref(BKd1B0)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BK(1/2)2", boost::cref(BKd1B1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BK(1/2)3", boost::cref(BKd1B2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BK(1/2)4", boost::cref(BKd1B3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BK(1/2)5", boost::cref(BKd1B4)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BK(1/2)6", boost::cref(BKd1B5)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BK(1/2)7", boost::cref(BKd1B6)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BK(1/2)8", boost::cref(BKd1B7)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BK(1/2)9", boost::cref(BKd1B8)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BK(1/2)10", boost::cref(BKd1B9)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BK(3/2)1", boost::cref(BKd3B0)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BK(3/2)2", boost::cref(BKd3B1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BK(3/2)3", boost::cref(BKd3B2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BK(3/2)4", boost::cref(BKd3B3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BK(3/2)5", boost::cref(BKd3B4)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BK(3/2)6", boost::cref(BKd3B5)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BK(3/2)7", boost::cref(BKd3B6)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BK(3/2)8", boost::cref(BKd3B7)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BK(3/2)9", boost::cref(BKd3B8)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BK(3/2)10", boost::cref(BKd3B9)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BKd_scale", boost::cref(BKd_scale)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ReA2_Kd", boost::cref(ReA2_Kd)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ReA0_Kd", boost::cref(ReA0_Kd)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Omega_eta_etap", boost::cref(Omega_eta_etap)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Br_Kp_P0enu", boost::cref(Br_Kp_P0enu)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Br_Kp_munu", boost::cref(Br_Kp_munu)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Br_B_Xcenu", boost::cref(Br_B_Xcenu)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("DeltaP_cu", boost::cref(DeltaP_cu)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("IB_Kl", boost::cref(IB_Kl)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("IB_Kp", boost::cref(IB_Kp)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("reh_0", boost::cref(reh_0)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("reh_p", boost::cref(reh_p)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("reh_m", boost::cref(reh_m)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("imh_0", boost::cref(imh_0)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("imh_p", boost::cref(imh_p)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("imh_m", boost::cref(imh_m)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("reh_0_1", boost::cref(reh_0_1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("reh_p_1", boost::cref(reh_p_1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("reh_m_1", boost::cref(reh_m_1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("imh_0_1", boost::cref(imh_0_1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("imh_p_1", boost::cref(imh_p_1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("imh_m_1", boost::cref(imh_m_1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("reh_0_2", boost::cref(reh_0_2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("reh_p_2", boost::cref(reh_p_2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("reh_m_2", boost::cref(reh_m_2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("imh_0_2", boost::cref(imh_0_2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("imh_p_2", boost::cref(imh_p_2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("imh_m_2", boost::cref(imh_m_2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("reh_0_MP", boost::cref(reh_0_MP)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("imh_0_MP", boost::cref(imh_0_MP)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("reh_0_1_MP", boost::cref(reh_0_1_MP)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("imh_0_1_MP", boost::cref(imh_0_1_MP)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_0V", boost::cref(a_0V)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_1V", boost::cref(a_1V)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_2V", boost::cref(a_2V)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MRV", boost::cref(MRV)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_0A0", boost::cref(a_0A0)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_1A0", boost::cref(a_1A0)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_2A0", boost::cref(a_2A0)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MRA0", boost::cref(MRA0)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_0A1", boost::cref(a_0A1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_1A1", boost::cref(a_1A1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_2A1", boost::cref(a_2A1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MRA1", boost::cref(MRA1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_0A12", boost::cref(a_0A12)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_1A12", boost::cref(a_1A12)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_2A12", boost::cref(a_2A12)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MRA12", boost::cref(MRA12)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_0T1", boost::cref(a_0T1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_1T1", boost::cref(a_1T1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_2T1", boost::cref(a_2T1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MRT1", boost::cref(MRT1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_0T2", boost::cref(a_0T2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_1T2", boost::cref(a_1T2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_2T2", boost::cref(a_2T2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MRT2", boost::cref(MRT2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_0T23", boost::cref(a_0T23)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_1T23", boost::cref(a_1T23)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_2T23", boost::cref(a_2T23)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MRT23", boost::cref(MRT23)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_0Vphi", boost::cref(a_0Vphi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_1Vphi", boost::cref(a_1Vphi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_2Vphi", boost::cref(a_2Vphi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MRVphi", boost::cref(MRVphi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_0A0phi", boost::cref(a_0A0phi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_1A0phi", boost::cref(a_1A0phi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_2A0phi", boost::cref(a_2A0phi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MRA0phi", boost::cref(MRA0phi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_0A1phi", boost::cref(a_0A1phi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_1A1phi", boost::cref(a_1A1phi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_2A1phi", boost::cref(a_2A1phi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MRA1phi", boost::cref(MRA1phi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_0A12phi", boost::cref(a_0A12phi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_1A12phi", boost::cref(a_1A12phi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_2A12phi", boost::cref(a_2A12phi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MRA12phi", boost::cref(MRA12phi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_0T1phi", boost::cref(a_0T1phi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_1T1phi", boost::cref(a_1T1phi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_2T1phi", boost::cref(a_2T1phi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MRT1phi", boost::cref(MRT1phi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_0T2phi", boost::cref(a_0T2phi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_1T2phi", boost::cref(a_1T2phi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_2T2phi", boost::cref(a_2T2phi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MRT2phi", boost::cref(MRT2phi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_0T23phi", boost::cref(a_0T23phi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_1T23phi", boost::cref(a_1T23phi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("a_2T23phi", boost::cref(a_2T23phi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MRT23phi", boost::cref(MRT23phi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("r_1_fplus", boost::cref(r_1_fplus)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("r_2_fplus", boost::cref(r_2_fplus)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("m_fit2_fplus", boost::cref(m_fit2_fplus)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("r_1_fT", boost::cref(r_1_fT)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("r_2_fT", boost::cref(r_2_fT)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("m_fit2_fT", boost::cref(m_fit2_fT)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("r_2_f0", boost::cref(r_2_f0)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("m_fit2_f0", boost::cref(m_fit2_f0)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("bsgamma_E0", boost::cref(bsgamma_E0)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BLNPcorr", boost::cref(BLNPcorr)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Gambino_mukin", boost::cref(Gambino_mukin)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Gambino_BRsem", boost::cref(Gambino_BRsem)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Gambino_Mbkin", boost::cref(Gambino_Mbkin)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Gambino_Mcatmuc", boost::cref(Gambino_Mcatmuc)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Gambino_mupi2", boost::cref(Gambino_mupi2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Gambino_rhoD3", boost::cref(Gambino_rhoD3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Gambino_muG2", boost::cref(Gambino_muG2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Gambino_rhoLS3", boost::cref(Gambino_rhoLS3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("lambdaB", boost::cref(mesons[B_D].getLambdaM())));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("alpha1kst", boost::cref(mesons[K_star].getGegenalpha(0))));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("alpha2kst", boost::cref(mesons[K_star].getGegenalpha(1))));

    unknownParameterWarning = true;
}

std::string QCD::orderToString(const orders order) const
{
    switch (order) {
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
        default:
            throw std::runtime_error(orderToString(order) + " is not implemented in QCD::orderToString().");
    }
}

////////////////////////////////////////////////////////////////////////

bool QCD::Init(const std::map<std::string, double>& DPars)
{
    bool check = CheckParameters(DPars);
    if (!check) return (check);
    check *= Update(DPars);
    unknownParameterWarning = false;
    return (check);
}

bool QCD::PreUpdate()
{
    requireYu = false;
    requireYd = false;
    computeBd = false;
    computeFBd = false;
    computemt = false;

    return (true);
}

bool QCD::Update(const std::map<std::string, double>& DPars)
{
    if (!PreUpdate()) return (false);

    UpdateError = false;

    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);

    if (UpdateError) return (false);

    if (!PostUpdate()) return (false);

    return (true);
}

bool QCD::PostUpdate()
{
    if (computeFBd)
        mesons[B_D].setDecayconst(mesons[B_S].getDecayconst() / FBsoFBd);
    if (computeBd)
        BBd.setBpars(0, BBs.getBpars()(0) / BBsoBBd);
    if (computemt) {
        quarks[TOP].setMass(Mp2Mbar(mtpole, FULLNNLO));
        quarks[TOP].setMass_scale(quarks[TOP].getMass());
    }

    myh_0 = gslpp::complex(reh_0,imh_0,true);
    myh_p = gslpp::complex(reh_p,imh_p,true);
    myh_m = gslpp::complex(reh_m,imh_m,true);
    myh_0_2 = (3.*myh_0 + gslpp::complex(reh_0_2,imh_0_2,true) - 4. * 
            gslpp::complex(reh_0_1,imh_0_1,true))/12.;
    myh_p_2 = (3.*myh_p + gslpp::complex(reh_p_2,imh_p_2,true) - 4. * 
            gslpp::complex(reh_p_1,imh_p_1,true))/12.;
    myh_m_2 = (3.*myh_m + gslpp::complex(reh_m_2,imh_m_2,true) - 4. * 
            gslpp::complex(reh_m_1,imh_m_1,true))/12.;
    myh_0_1 = -myh_0 - myh_0_2 + gslpp::complex(reh_0_1,imh_0_1,true);
    myh_p_1 = -myh_p - myh_p_2 + gslpp::complex(reh_p_1,imh_p_1,true);
    myh_m_1 = -myh_m - myh_m_2 + gslpp::complex(reh_m_1,imh_m_1,true);
    
    return (true);
}

void QCD::setParameter(const std::string name, const double& value)
{
    if (name.compare("AlsM") == 0) {
        AlsM = value;
        computemt = true;
        requireYu = true;
        requireYd = true;
    } else if (name.compare("MAls") == 0) {
        MAls = value;
        computemt = true;
        requireYu = true;
        requireYd = true;
    } else if (name.compare("mup") == 0) {
        if (value < MEPS) UpdateError = true;
        quarks[UP].setMass(value);
        requireYu = true;
    } else if (name.compare("mdown") == 0) {
        if (value < MEPS) UpdateError = true;
        quarks[DOWN].setMass(value);
        requireYd = true;
    } else if (name.compare("mcharm") == 0) {
        quarks[CHARM].setMass(value);
        quarks[CHARM].setMass_scale(value);
        requireYu = true;
    } else if (name.compare("mstrange") == 0) {
        if (value < MEPS) UpdateError = true;
        quarks[STRANGE].setMass(value);
        requireYd = true;
    } else if (name.compare("mtop") == 0) {
        mtpole = value;
        requireYu = true;
        computemt = true;
    } else if (name.compare("mbottom") == 0) {
        quarks[BOTTOM].setMass(value);
        quarks[BOTTOM].setMass_scale(value);
        requireYd = true;
    } else if (name.compare("muc") == 0)
        muc = value;
    else if (name.compare("mub") == 0)
        mub = value;
    else if (name.compare("mut") == 0)
        mut = value;
    else if (name.compare("MK0") == 0)
        mesons[K_0].setMass(value);
    else if (name.compare("MKp") == 0)
        mesons[K_P].setMass(value);
    else if (name.compare("MD") == 0)
        mesons[D_0].setMass(value);
    else if (name.compare("MBd") == 0)
        mesons[B_D].setMass(value);
    else if (name.compare("MBp") == 0)
        mesons[B_P].setMass(value);
    else if (name.compare("MBs") == 0)
        mesons[B_S].setMass(value);
    else if (name.compare("MKstar") == 0)
        mesons[K_star].setMass(value);
    else if (name.compare("Mphi") == 0)
        mesons[PHI].setMass(value);
    else if (name.compare("tKl") == 0)
        mesons[K_0].setLifetime(value);
    else if (name.compare("tKp") == 0)
        mesons[K_P].setLifetime(value);
    else if (name.compare("tBd") == 0)
        mesons[B_D].setLifetime(value);
    else if (name.compare("tBp") == 0)
        mesons[B_P].setLifetime(value);
    else if (name.compare("tBs") == 0)
        mesons[B_S].setLifetime(value);
    else if (name.compare("DGs_Gs") == 0)
        mesons[B_S].setDgamma_gamma(value);
    else if (name.compare("tKstar") == 0)
        mesons[K_star].setLifetime(value);
    else if (name.compare("tphi") == 0)
        mesons[PHI].setLifetime(value);
        //else if(name.compare("FP")==0) {
        //    mesons[P_0].setDecayconst(value);
        //    mesons[P_P].setDecayconst(value);
        //}
    else if (name.compare("FK") == 0)
        mesons[K_0].setDecayconst(value);
    else if (name.compare("FD") == 0)
        mesons[D_0].setDecayconst(value);
    else if (name.compare("FKstar") == 0)
        mesons[K_star].setDecayconst(value);
    else if (name.compare("FKstarp") == 0)
        FKstarp = value;
    else if (name.compare("Fphi") == 0)
        mesons[PHI].setDecayconst(value);
    else if (name.compare("FBs") == 0) {
        mesons[B_S].setDecayconst(value);
        computeFBd = true;
    } else if (name.compare("FBsoFBd") == 0) {
        FBsoFBd = value;
        computeFBd = true;
    } else if (name.compare("BK1") == 0) {
        BK.setBpars(0, value);
        BKB0 = value;
    } else if (name.compare("BK2") == 0) {
        BK.setBpars(1, value);
        BKB1 = value;
    } else if (name.compare("BK3") == 0) {
        BK.setBpars(2, value);
        BKB2 = value;
    } else if (name.compare("BK4") == 0) {
        BK.setBpars(3, value);
        BKB3 = value;
    } else if (name.compare("BK5") == 0) {
        BK.setBpars(4, value);
        BKB4 = value;
    } else if (name.compare("BKscale") == 0) {
        BK.setMu(value);
        BKscale = value;
    } else if (name.compare("BKscheme") == 0)
        BK.setScheme((schemes) value);
    else if (name.compare("BD1") == 0) {
        BD.setBpars(0, value);
        BDB0 = value;
    } else if (name.compare("BD2") == 0) {
        BD.setBpars(1, value);
        BDB1 = value;
    } else if (name.compare("BD3") == 0) {
        BD.setBpars(2, value);
        BDB2 = value;
    } else if (name.compare("BD4") == 0) {
        BD.setBpars(3, value);
        BDB3 = value;
    } else if (name.compare("BD5") == 0) {
        BD.setBpars(4, value);
        BDB4 = value;
    } else if (name.compare("BDscale") == 0) {
        BD.setMu(value);
        BDscale = value;
    } else if (name.compare("BDscheme") == 0)
        BD.setScheme((schemes) value);
    else if (name.compare("BBsoBBd") == 0) {
        BBsoBBd = value;
        computeBd = true;
    } else if (name.compare("BBs1") == 0) {
        BBs.setBpars(0, value);
        BBsB0 = value;
        computeBd = true;
    } else if (name.compare("BBs2") == 0) {
        BBd.setBpars(1, value);
        BBs.setBpars(1, value);
        BBsB1 = value;
        //BBdB1 = value;
    } else if (name.compare("BBs3") == 0) {
        BBd.setBpars(2, value);
        BBs.setBpars(2, value);
        BBsB2 = value;
        //BBdB2 = value;
    } else if (name.compare("BBs4") == 0) {
        BBd.setBpars(3, value);
        BBs.setBpars(3, value);
        BBsB3 = value;
        //BBdB3 = value;
    } else if (name.compare("BBs5") == 0) {
        BBd.setBpars(4, value);
        BBs.setBpars(4, value);
        BBsB4 = value;
        //BBdB4 = value;
    } else if (name.compare("BBsscale") == 0) {
        BBd.setMu(value);
        BBs.setMu(value);
        BBsscale = value;
        //BBdscale = value;
    } else if (name.compare("BBsscheme") == 0) {
        BBd.setScheme((schemes) value);
        BBs.setScheme((schemes) value);
    } else if (name.compare("BK(1/2)1") == 0) {
        BKd1.setBpars(0, value);
        BKd1B0 = value;
    } else if (name.compare("BK(1/2)2") == 0) {
        BKd1.setBpars(1, value);
        BKd1B1 = value;
    } else if (name.compare("BK(1/2)3") == 0) {
        BKd1.setBpars(2, value);
        BKd1B2 = value;
    } else if (name.compare("BK(1/2)4") == 0) {
        BKd1.setBpars(3, value);
        BKd1B3 = value;
    } else if (name.compare("BK(1/2)5") == 0) {
        BKd1.setBpars(4, value);
        BKd1B4 = value;
    } else if (name.compare("BK(1/2)6") == 0){
        BKd1.setBpars(5, value);
        BKd1B5 = value;
    } else if (name.compare("BK(1/2)7") == 0){
        BKd1.setBpars(6, value);
        BKd1B6 = value;
    } else if (name.compare("BK(1/2)8") == 0) {
        BKd1.setBpars(7, value);
        BKd1B7 = value;
    } else if (name.compare("BK(1/2)9") == 0) {
        BKd1.setBpars(8, value);
        BKd1B8 = value;
    } else if (name.compare("BK(1/2)10") == 0) {
        BKd1.setBpars(9, value);
        BKd1B9 = value;
    } else if (name.compare("BK(3/2)1") == 0) {
        BKd3.setBpars(0, value);
        BKd3B0 = value;
    } else if (name.compare("BK(3/2)2") == 0) {
        BKd3.setBpars(1, value);
        BKd3B1 = value;
    } else if (name.compare("BK(3/2)3") == 0) {
        BKd3.setBpars(2, value);
        BKd3B2 = value;
    } else if (name.compare("BK(3/2)4") == 0) {
        BKd3.setBpars(3, value);
        BKd3B3 = value;
    } else if (name.compare("BK(3/2)5") == 0) {
        BKd3.setBpars(4, value);
        BKd3B4 = value;
    } else if (name.compare("BK(3/2)6") == 0) {
        BKd3.setBpars(5, value);
        BKd3B5 = value;
    } else if (name.compare("BK(3/2)7") == 0) {
        BKd3.setBpars(6, value);
        BKd3B6 = value;
    } else if (name.compare("BK(3/2)8") == 0) {
        BKd3.setBpars(7, value);
        BKd3B7 = value;
    } else if (name.compare("BK(3/2)9") == 0) {
        BKd3.setBpars(8, value);
        BKd3B8 = value;
    } else if (name.compare("BK(3/2)10") == 0) {
        BKd3.setBpars(9, value);
        BKd3B9 = value;
    } else if (name.compare("BKd_scale") == 0) {
        BKd1.setMu(value);
        BKd3.setMu(value);
        BKd_scale = value;
    } else if (name.compare("BKd_scheme") == 0) {
        BKd1.setScheme((schemes) value);
        BKd3.setScheme((schemes) value);
    } else if (name.compare("ReA0_Kd") == 0)
        ReA0_Kd = value;
    else if (name.compare("ReA2_Kd") == 0)
        ReA2_Kd = value;
    else if (name.compare("Omega_eta_etap") == 0)
        Omega_eta_etap = value;
    else if (name.compare("Br_Kp_P0enu") == 0)
        Br_Kp_P0enu = value;
    else if (name.compare("Br_Kp_munu") == 0)
        Br_Kp_munu = value;
    else if (name.compare("Br_B_Xcenu") == 0)
        Br_B_Xcenu = value;
    else if (name.compare("DeltaP_cu") == 0)
        DeltaP_cu = value;
    else if (name.compare("IB_Kl") == 0)
        IB_Kl = value;
    else if (name.compare("IB_Kp") == 0)
        IB_Kp = value;
    else if (name.compare("reh_0") == 0)
        reh_0 = value;
    else if (name.compare("reh_p") == 0)
        reh_p = value;
    else if (name.compare("reh_m") == 0)
        reh_m = value;
    else if (name.compare("imh_0") == 0)
        imh_0 = value;
    else if (name.compare("imh_p") == 0)
        imh_p = value;
    else if (name.compare("imh_m") == 0)
        imh_m = value;
    else if (name.compare("reh_0_1") == 0)
        reh_0_1 = value;
    else if (name.compare("reh_p_1") == 0)
        reh_p_1 = value;
    else if (name.compare("reh_m_1") == 0)
        reh_m_1 = value;
    else if (name.compare("imh_0_1") == 0)
        imh_0_1 = value;
    else if (name.compare("imh_p_1") == 0)
        imh_p_1 = value;
    else if (name.compare("imh_m_1") == 0)
        imh_m_1 = value;
    else if (name.compare("reh_0_2") == 0)
        reh_0_2 = value;
    else if (name.compare("reh_p_2") == 0)
        reh_p_2 = value;
    else if (name.compare("reh_m_2") == 0)
        reh_m_2 = value;
    else if (name.compare("imh_0_2") == 0)
        imh_0_2 = value;
    else if (name.compare("imh_p_2") == 0)
        imh_p_2 = value;
    else if (name.compare("imh_m_2") == 0)
        imh_m_2 = value;
    else if (name.compare("reh_0_MP") == 0)
        reh_0_MP = value;
    else if (name.compare("imh_0_MP") == 0)
        imh_0_MP = value;
    else if (name.compare("reh_0_1_MP") == 0)
        reh_0_1_MP = value;
    else if (name.compare("imh_0_1_MP") == 0)
        imh_0_1_MP = value;
    else if (name.compare("a_0V") == 0)
        a_0V = value;
    else if (name.compare("a_1V") == 0)
        a_1V = value;
    else if (name.compare("a_2V") == 0)
        a_2V = value;
    else if (name.compare("MRV") == 0)
        MRV = value;
    else if (name.compare("a_0A0") == 0)
        a_0A0 = value;
    else if (name.compare("a_1A0") == 0)
        a_1A0 = value;
    else if (name.compare("a_2A0") == 0)
        a_2A0 = value;
    else if (name.compare("MRA0") == 0)
        MRA0 = value;
    else if (name.compare("a_0A1") == 0)
        a_0A1 = value;
    else if (name.compare("a_1A1") == 0)
        a_1A1 = value;
    else if (name.compare("a_2A1") == 0)
        a_2A1 = value;
    else if (name.compare("MRA1") == 0)
        MRA1 = value;
    else if (name.compare("a_0A12") == 0)
        a_0A12 = value;
    else if (name.compare("a_1A12") == 0)
        a_1A12 = value;
    else if (name.compare("a_2A12") == 0)
        a_2A12 = value;
    else if (name.compare("MRA12") == 0)
        MRA12 = value;
    else if (name.compare("a_0T1") == 0)
        a_0T1 = value;
    else if (name.compare("a_1T1") == 0)
        a_1T1 = value;
    else if (name.compare("a_2T1") == 0)
        a_2T1 = value;
    else if (name.compare("MRT1") == 0)
        MRT1 = value;
    else if (name.compare("a_0T2") == 0)
        a_0T2 = value;
    else if (name.compare("a_1T2") == 0)
        a_1T2 = value;
    else if (name.compare("a_2T2") == 0)
        a_2T2 = value;
    else if (name.compare("MRT2") == 0)
        MRT2 = value;
    else if (name.compare("a_0T23") == 0)
        a_0T23 = value;
    else if (name.compare("a_1T23") == 0)
        a_1T23 = value;
    else if (name.compare("a_2T23") == 0)
        a_2T23 = value;
    else if (name.compare("MRT23") == 0)
        MRT23 = value;
    else if (name.compare("a_0Vphi") == 0)
        a_0Vphi = value;
    else if (name.compare("a_1Vphi") == 0)
        a_1Vphi = value;
    else if (name.compare("a_2Vphi") == 0)
        a_2Vphi = value;
    else if (name.compare("MRVphi") == 0)
        MRVphi = value;
    else if (name.compare("a_0A0phi") == 0)
        a_0A0phi = value;
    else if (name.compare("a_1A0phi") == 0)
        a_1A0phi = value;
    else if (name.compare("a_2A0phi") == 0)
        a_2A0phi = value;
    else if (name.compare("MRA0phi") == 0)
        MRA0phi = value;
    else if (name.compare("a_0A1phi") == 0)
        a_0A1phi = value;
    else if (name.compare("a_1A1phi") == 0)
        a_1A1phi = value;
    else if (name.compare("a_2A1phi") == 0)
        a_2A1phi = value;
    else if (name.compare("MRA1phi") == 0)
        MRA1phi = value;
    else if (name.compare("a_0A12phi") == 0)
        a_0A12phi = value;
    else if (name.compare("a_1A12phi") == 0)
        a_1A12phi = value;
    else if (name.compare("a_2A12phi") == 0)
        a_2A12phi = value;
    else if (name.compare("MRA12phi") == 0)
        MRA12phi = value;
    else if (name.compare("a_0T1phi") == 0)
        a_0T1phi = value;
    else if (name.compare("a_1T1phi") == 0)
        a_1T1phi = value;
    else if (name.compare("a_2T1phi") == 0)
        a_2T1phi = value;
    else if (name.compare("MRT1phi") == 0)
        MRT1phi = value;
    else if (name.compare("a_0T2phi") == 0)
        a_0T2phi = value;
    else if (name.compare("a_1T2phi") == 0)
        a_1T2phi = value;
    else if (name.compare("a_2T2phi") == 0)
        a_2T2phi = value;
    else if (name.compare("MRT2phi") == 0)
        MRT2phi = value;
    else if (name.compare("a_0T23phi") == 0)
        a_0T23phi = value;
    else if (name.compare("a_1T23phi") == 0)
        a_1T23phi = value;
    else if (name.compare("a_2T23phi") == 0)
        a_2T23phi = value;
    else if (name.compare("MRT23phi") == 0)
        MRT23phi = value;
    else if (name.compare("r_1_fplus") == 0)
        r_1_fplus = value;
    else if (name.compare("r_2_fplus") == 0)
        r_2_fplus = value;
    else if (name.compare("m_fit2_fplus") == 0)
        m_fit2_fplus = value;
    else if (name.compare("r_1_fT") == 0)
        r_1_fT = value;
    else if (name.compare("r_2_fT") == 0)
        r_2_fT = value;
    else if (name.compare("m_fit2_fT") == 0)
        m_fit2_fT = value;
    else if (name.compare("r_2_f0") == 0)
        r_2_f0 = value;
    else if (name.compare("m_fit2_f0") == 0)
        m_fit2_f0 = value;
    else if (name.compare("bsgamma_E0") == 0)
        bsgamma_E0 = value;
    else if (name.compare("BLNPcorr") == 0)
        BLNPcorr = value;
    else if (name.compare("Gambino_mukin") == 0)
        Gambino_mukin = value;
    else if (name.compare("Gambino_BRsem") == 0)
        Gambino_BRsem = value;
    else if (name.compare("Gambino_Mbkin") == 0)
        Gambino_Mbkin = value;
    else if (name.compare("Gambino_Mcatmuc") == 0)
        Gambino_Mcatmuc = value;
    else if (name.compare("Gambino_mupi2") == 0)
        Gambino_mupi2 = value;
    else if (name.compare("Gambino_rhoD3") == 0)
        Gambino_rhoD3 = value;
    else if (name.compare("Gambino_muG2") == 0)
        Gambino_muG2 = value;
    else if (name.compare("Gambino_rhoLS3") == 0)
        Gambino_rhoLS3 = value;
    else if (name.compare("lambdaB") == 0)
        mesons[B_D].setLambdaM(value);
    else if (name.compare("alpha1kst") == 0)
        mesons[K_star].setGegenalpha(0,value);
    else if (name.compare("alpha2kst") == 0)
        mesons[K_star].setGegenalpha(1,value);
    else
        if (unknownParameterWarning)
        std::cout << "WARNING: unknown parameter " << name << " in model initialization" << std::endl;
}

bool QCD::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NQCDvars; i++)
        if (DPars.find(QCDvars[i]) == DPars.end()) {
            std::cout << "missing mandatory QCD parameter " << QCDvars[i] << std::endl;
            return false;
        }
    return true;
}

////////////////////////////////////////////////////////////////////////

bool QCD::setFlag(const std::string name, const bool value)
{
    std::cout << "WARNING: wrong name or value for ModelFlag " << name << std::endl;
    return (false);
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

double QCD::Thresholds(const int i) const
{
    if (!(mut > mub && mub > muc))
        throw std::runtime_error("inverted thresholds in QCD::Thresholds()!");

    switch (i) {
        case 0: return (1.0E10);
        case 1: return (mut);
        case 2: return (mub);
        case 3: return (muc);
        default: return (0.);
    }
}

double QCD::AboveTh(const double mu) const
{
    int i;
    for (i = 4; i >= 0; i--)
        if (mu < Thresholds(i)) return (Thresholds(i));

    throw std::runtime_error("Error in QCD::AboveTh()");
}

double QCD::BelowTh(const double mu) const
{
    int i;
    for (i = 0; i < 5; i++)
        if (mu >= Thresholds(i)) return (Thresholds(i));

    throw std::runtime_error("Error in QCD::BelowTh()");
}

double QCD::Nf(const double mu) const
{
    int i;
    for (i = 1; i < 5; i++)
        if (mu >= Thresholds(i))
            return (7. - (double) i);

    throw std::runtime_error("Error in QCD::Nf()");
}

void QCD::CacheShift(double cache[][CacheSize], int n) const
{
    int i, j;
    for (i = CacheSize - 1; i > 0; i--)
        for (j = 0; j < n; j++)
            cache[j][i] = cache[j][i - 1];
}

////////////////////////////////////////////////////////////////////////

double QCD::Beta0(const double nf) const
{
    return ( (11. * Nc - 2. * nf) / 3.);
}

double QCD::Beta1(const double nf) const
{
    return ( 34. / 3. * Nc * Nc - 10. / 3. * Nc * nf - 2. * CF * nf);
}

double QCD::Beta2(const double nf) const
{
    return ( 2857. / 54. * Nc * Nc * Nc + CF * CF * nf - 205. / 18. * CF * Nc * nf
            - 1415. / 54. * Nc * Nc * nf + 11. / 9. * CF * nf * nf + 79. / 54. * Nc * nf * nf);
}

double QCD::AlsWithInit(const double mu, const double alsi, const double mu_i,
        const orders order) const
{
    double nf = Nf(mu);
    if (nf != Nf(mu_i))
        throw std::runtime_error("Error in QCD::AlsWithInit().");

    double v = 1. - Beta0(nf) * alsi / 2. / M_PI * log(mu_i / mu);
    switch (order) {
        case LO:
            return (alsi / v);
        case FULLNLO:
            return (alsi / v * (1. - Beta1(nf) / Beta0(nf) * alsi / 4. / M_PI * log(v) / v));
        case NLO:
            return (alsi / v * (-Beta1(nf) / Beta0(nf) * alsi / 4. / M_PI * log(v) / v));
        default:
            throw std::runtime_error(orderToString(order) + " is not implemented in QCD::Als(mu,alsi,mi,order).");
    }
}

double QCD::Als4(const double mu) const
{
    double v = 1. - Beta0(4.) * AlsM / 2. / M_PI * log(MAls / mu);
    return (AlsM / v * (1. - Beta1(4.) / Beta0(4.) * AlsM / 4. / M_PI * log(v) / v));
}

double QCD::Alstilde5(const double mu) const
{
    double mu_0 = MAls;
    double alphatilde_e = alphaMz()/4./M_PI;
    double alphatilde_s = AlsM/4./M_PI;
    unsigned int nf = 5;

    double B00S = Beta0(nf), B10S = Beta1(nf), B20S = Beta2(nf), B30S = gsl_sf_zeta_int(3) * 352864./81. - 598391./1458,
            B01S = -22./9., B11S = -308./27., B02S = 4945./243.; 

    double B00E = 80./9., B01E = 176./9., B10E = 464./27.; 

    double B10soB00s = B10S / B00S;
    double B01soB00e = B01S/B00E;

    double vs= 1. + 2. * B00S * alphatilde_s * log(mu/ mu_0);
    double ve= 1. - 2. * B00E * alphatilde_e * log(mu/ mu_0);
    double ps= B00S * alphatilde_s /(B00S * alphatilde_s + B00E * alphatilde_e);

    double logve = log(ve);
    double logvs = log(vs);
    double logeos = log(ve/vs);
    double logsoe = log(vs/ve);
    double asovs = alphatilde_s/vs;
    double aeove = alphatilde_e/ve;

    double result = 0;

    result = asovs - pow(asovs, 2) * (logvs * B10soB00s - logve * B01soB00e) 
            +  pow(asovs, 3) * ((1. - vs) * B20S / B00S + B10soB00s * B10soB00s * (logvs * logvs - logvs
            + vs - 1.) + B01soB00e * B01soB00e * logve * logve + (-2. * logvs * logve 
            + ps * ve * logve) * B01S * B10S/(B00E * B00S)) 
            +  pow(asovs, 4) * (0.5 * B30S *(1. - vs * vs)/ B00S + ((2. * vs - 3.) * logvs + vs * vs 
            - vs) * B20S * B10soB00s /(B00S) + B10soB00s * B10soB00s * B10soB00s * (- pow(logvs,3) 
            + 5. * pow(logvs,2) / 2. + 2. * (1. - vs) * logvs - (vs - 1.) * (vs - 1.)* 0.5)) 
            + pow(asovs, 2) * (aeove) * ((ve - 1.) * B02S / B00E 
            + ps * ve * logeos * B11S /B00S +(logve - ve + 1.) * B01soB00e * B10E/(B00S) 
            + logvs * ve * ps * B01S * B10soB00s/(B00S) +(logsoe * ve * ps - logvs) * B01soB00e * B01E/( B00S));
    return (result);
}

double QCD::AlsWithLambda(const double mu, const double logLambda,
        const orders order) const
{
    double nf = Nf(mu);
    double L = 2. * (log(mu) - logLambda);

    // LO contribution
    double b0 = Beta0(nf);
    double b0L = b0*L;
    double alsLO = 4. * M_PI / b0L;
    if (order == LO) return alsLO;

    // NLO contribution
    double b1 = Beta1(nf);
    double log_L = log(L);
    double alsNLO = 4. * M_PI / b0L * (-b1 * log_L / b0 / b0L);
    if (order == NLO) return alsNLO;
    if (order == FULLNLO) return (alsLO + alsNLO);

    // NNLO contribution
    double b2 = Beta2(nf);
    double alsNNLO = 4. * M_PI / b0L * (1. / b0L / b0L
            * (b1 * b1 / b0 / b0 * (log_L * log_L - log_L - 1.) + b2 / b0));
    if (order == NNLO) return alsNNLO;
    if (order == FULLNNLO) return (alsLO + alsNLO + alsNNLO);

    throw std::runtime_error(orderToString(order) + " is not implemented in QCD::AlsWithLambda().");
}

double QCD::AlsWithLambda(const double mu, const orders order) const
{
    return AlsWithLambda(mu, logLambda(Nf(mu), order), order);
}

double QCD::Als(const double mu, const orders order) const
{
    int i;
    for (i = 0; i < CacheSize; ++i)
        if ((mu == als_cache[0][i]) && ((double) order == als_cache[1][i]) &&
                (AlsM == als_cache[2][i]) && (MAls == als_cache[3][i]) &&
                (mut == als_cache[4][i]) && (mub == als_cache[5][i]) &&
                (muc == als_cache[6][i]))
            return als_cache[7][i];

    double nfmu = Nf(mu), nfz = Nf(MAls), mu_thre1, mu_thre2, Als_tmp, mf;
    double als;

    switch (order) {
        case LO:
        case FULLNLO:
        case NLO:
            if (nfmu == nfz)
                als = AlsWithInit(mu, AlsM, MAls, order);
            else if (nfmu > nfz) {
                if (order == NLO)
                    throw std::runtime_error("NLO is not implemented in QCD::Als(mu,order).");
                if (nfmu == nfz + 1.) {
                    mu_thre1 = AboveTh(MAls); // mut
                    Als_tmp = AlsWithInit(mu_thre1 - MEPS, AlsM, MAls, order);
                    if (order == FULLNLO) {
                        mf = mtpole;
                        Als_tmp = (1. + Als_tmp / M_PI * log(mu_thre1 / mf) / 3.) * Als_tmp;
                    }
                    als = AlsWithInit(mu, Als_tmp, mu_thre1 + MEPS, order);
                } else
                    throw std::runtime_error("Error in QCD::Als(mu,order)");
            } else {
                if (order == NLO)
                    throw std::runtime_error("NLO is not implemented in QCD::Als(mu,order).");
                if (nfmu == nfz - 1.) {
                    mu_thre1 = BelowTh(MAls); // mub
                    Als_tmp = AlsWithInit(mu_thre1 + MEPS, AlsM, MAls, order);
                    if (order == FULLNLO) {
                        mf = getQuarks(BOTTOM).getMass();
                        Als_tmp = (1. - Als_tmp / M_PI * log(mu_thre1 / mf) / 3.) * Als_tmp;
                    }
                    als = AlsWithInit(mu, Als_tmp, mu_thre1 - MEPS, order);
                } else if (nfmu == nfz - 2.) {
                    mu_thre1 = BelowTh(MAls); // mub
                    mu_thre2 = AboveTh(mu); // muc
                    Als_tmp = Als(mu_thre1 + MEPS, order);
                    if (order == FULLNLO) {
                        mf = getQuarks(BOTTOM).getMass();
                        Als_tmp = (1. - Als_tmp / M_PI * log(mu_thre1 / mf) / 3.) * Als_tmp;
                    }
                    Als_tmp = AlsWithInit(mu_thre2, Als_tmp, mu_thre1 - MEPS, order);
                    if (order == FULLNLO) {
                        mf = getQuarks(CHARM).getMass();
                        Als_tmp = (1. - Als_tmp / M_PI * log(mu_thre2 / mf) / 3.) * Als_tmp;
                    }
                    als = AlsWithInit(mu, Als_tmp, mu_thre2 - MEPS, order);
                } else
                    throw std::runtime_error("Error in QCD::Als(mu,order)");
            }
            break;
        case FULLNNLO:
        case NNLO:
            /* alpha_s(mu) computed with Lambda_QCD for Nf=nfmu */
            als = AlsWithLambda(mu, order);
            break;
        default:
            throw std::runtime_error(orderToString(order) + " is not implemented in QCD::Als(mu,order).");
    }

    CacheShift(als_cache, 8);
    als_cache[0][0] = mu;
    als_cache[1][0] = (double) order;
    als_cache[2][0] = AlsM;
    als_cache[3][0] = MAls;
    als_cache[4][0] = mut;
    als_cache[5][0] = mub;
    als_cache[6][0] = muc;
    als_cache[7][0] = als;

    return als;
}

double QCD::ZeroNf6NLO(double *logLambda6, double *logLambda5_in) const
{
    return ( AlsWithLambda(mut + 1.e-10, *logLambda6, FULLNLO)
            - AlsWithLambda(mut - 1.e-10, *logLambda5_in, FULLNLO));
}

double QCD::ZeroNf5(double *logLambda5, double *order) const
{
    return ( AlsWithLambda(MAls, *logLambda5, (orders) * order) - AlsM);
}

double QCD::ZeroNf4NLO(double *logLambda4, double *logLambda5_in) const
{
    return ( AlsWithLambda(mub - 1.e-10, *logLambda4, FULLNLO)
            - AlsWithLambda(mub + 1.e-10, *logLambda5_in, FULLNLO));
}

double QCD::ZeroNf3NLO(double *logLambda3, double *logLambda4_in) const
{
    return ( AlsWithLambda(muc - 1.e-10, *logLambda3, FULLNLO)
            - AlsWithLambda(muc + 1.e-10, *logLambda4_in, FULLNLO));
}

double QCD::logLambda5(orders order) const
{
    if (order == NLO) order = FULLNLO;
    if (order == NNLO) order = FULLNNLO;

    for (int i = 0; i < CacheSize; ++i)
        if ((AlsM == logLambda5_cache[0][i])
                && (MAls == logLambda5_cache[1][i])
                && ((double) order == logLambda5_cache[2][i]))
            return logLambda5_cache[3][i];

    CacheShift(logLambda5_cache, 4);
    logLambda5_cache[0][0] = AlsM;
    logLambda5_cache[1][0] = MAls;
    logLambda5_cache[2][0] = (double) order;

    if (order == LO)
        logLambda5_cache[3][0] = log(MAls) - 2. * M_PI / Beta0(5.) / AlsM;
    else {
        double xmin = -4., xmax = -0.2;
        TF1 f = TF1("f", this, &QCD::ZeroNf5, xmin, xmax, 1, "QCD", "zeroNf5");

        ROOT::Math::WrappedTF1 wf1(f);
        double ledouble = (double) order;
        wf1.SetParameters(&ledouble);

        ROOT::Math::BrentRootFinder brf;
        brf.SetFunction(wf1, xmin, xmax);

        if (brf.Solve()) logLambda5_cache[3][0] = brf.Root();
        else
            throw std::runtime_error("Error in QCD::logLambda5()");
    }
    return ( logLambda5_cache[3][0]);
}

double QCD::logLambdaNLO(const double nfNEW, const double nfORG,
        const double logLambdaORG) const
{
    for (int i = 0; i < CacheSize; ++i)
        if ((AlsM == logLambdaNLO_cache[0][i])
                && (MAls == logLambdaNLO_cache[1][i])
                && (mut == logLambdaNLO_cache[2][i])
                && (mub == logLambdaNLO_cache[3][i])
                && (muc == logLambdaNLO_cache[4][i])
                && (nfNEW == logLambdaNLO_cache[5][i])
                && (nfORG == logLambdaNLO_cache[6][i])
                && (logLambdaORG == logLambdaNLO_cache[7][i]))
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
    if (nfNEW == 6. && nfORG == 5.) {
        f = TF1("f", this, &QCD::ZeroNf6NLO, xmin, xmax, 1, "QCD", "zeroNf6NLO");
    } else if (nfNEW == 4. && nfORG == 5.) {
        f = TF1("f", this, &QCD::ZeroNf4NLO, xmin, xmax, 1, "QCD", "zeroNf4NLO");
    } else if (nfNEW == 3. && nfORG == 4.) {
        f = TF1("f", this, &QCD::ZeroNf3NLO, xmin, xmax, 1, "QCD", "zeroNf3NLO");
    } else
        throw std::runtime_error("Error in QCD::logLambdaNLO()");

    ROOT::Math::WrappedTF1 wf1(f);
    wf1.SetParameters(&logLambdaORG);

    ROOT::Math::BrentRootFinder brf;
    brf.SetFunction(wf1, xmin, xmax);

    if (brf.Solve()) logLambdaNLO_cache[8][0] = brf.Root();
    else
        throw std::runtime_error("Error in QCD::logLambdaNLO()");

    return ( logLambdaNLO_cache[8][0]);
}

double QCD::logLambda(const double muMatching, const double mf,
        const double nfNEW, const double nfORG,
        const double logLambdaORG, orders order) const
{
    if (fabs(nfNEW - nfORG) != 1.)
        throw std::runtime_error("Error in QCD::logLambda()");
    if (order == NLO) order = FULLNLO;
    if (order == NNLO) order = FULLNNLO;

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
    logLambdaNEW = 1. / 2. / Beta0(nfNEW)
            *(Beta0(nfNEW) - Beta0(nfORG)) * L + logLambdaORG;

    // NLO contribution
    if (order == FULLNLO || order == FULLNNLO) {
        rNEW = Beta1(nfNEW) / Beta0(nfNEW);
        rORG = Beta1(nfORG) / Beta0(nfORG);
        log_mu2_mf2 = 2. * (logMuMatching - log(mf));
        log_L = log(L);
        if (nfNEW < nfORG)
            C1 = 2. / 3. * log_mu2_mf2;
        else
            C1 = -2. / 3. * log_mu2_mf2;
        logLambdaNEW += 1. / 2. / Beta0(nfNEW)
                *((rNEW - rORG) * log_L
                - rNEW * log(Beta0(nfNEW) / Beta0(nfORG)) - C1);
    }

    // NNLO contribution
    if (order == FULLNNLO) {
        if (nfNEW == 5. && nfORG == 6.)
            C2 = -16. * (log_mu2_mf2 * log_mu2_mf2 / 36. - 19. / 24. * log_mu2_mf2 - 7. / 24.);
        else if (nfNEW == 6. && nfORG == 5.)
            C2 = -16. * (log_mu2_mf2 * log_mu2_mf2 / 36. + 19. / 24. * log_mu2_mf2 + 7. / 24.);
        else {
            if (nfNEW < nfORG)
                C2 = -16. * (log_mu2_mf2 * log_mu2_mf2 / 36. - 19. / 24. * log_mu2_mf2 + 11. / 72.);
            else
                C2 = -16. * (log_mu2_mf2 * log_mu2_mf2 / 36. + 19. / 24. * log_mu2_mf2 - 11. / 72.);
        }
        logLambdaNEW += 1. / 2. / Beta0(nfNEW) / Beta0(nfORG) / L
                * (rORG * (rNEW - rORG) * log_L + rNEW * rNEW - rORG * rORG
                - Beta2(nfNEW) / Beta0(nfNEW) + Beta2(nfORG) / Beta0(nfORG)
                + rNEW * C1 - C1 * C1 - C2);
    }

    return logLambdaNEW;
}

double QCD::logLambda(const double nf, orders order) const
{
    if (order == NLO) order = FULLNLO;
    if (order == NNLO) order = FULLNNLO;

    double muMatching, mf, logLambdaORG, logLambdaNEW;
    if (nf == 5.)
        return logLambda5(order);
    else if (nf == 6.) {
        muMatching = Thresholds(1); // mut
        /* matching condition from Nf=5 to Nf=6 is given in terms of the top pole mass. */
        mf = mtpole; // top pole mass
        return logLambda(muMatching, mf, 6., 5., logLambda5(order), order);
    } else if (nf == 4. || nf == 3.) {
        muMatching = Thresholds(2); // mub
        mf = getQuarks(BOTTOM).getMass(); // m_b(m_b)
        logLambdaORG = logLambda5(order);
        logLambdaNEW = logLambda(muMatching, mf, 4., 5., logLambdaORG, order);
        if (nf == 3.) {
            muMatching = Thresholds(3); // muc
            mf = getQuarks(CHARM).getMass(); // m_c(m_c)
            logLambdaORG = logLambdaNEW;
            logLambdaNEW = logLambda(muMatching, mf, 3., 4., logLambdaORG, order);
        }
        return logLambdaNEW;
    } else
        throw std::runtime_error("Error in QCD::logLambda()");
}

////////////////////////////////////////////////////////////////////////

double QCD::Gamma0(const double nf) const
{
    return ( 6. * CF);
}

double QCD::Gamma1(const double nf) const
{
    return ( CF * (3. * CF + 97. / 3. * Nc - 10. / 3. * nf));
}

double QCD::Gamma2(const double nf) const
{
    return ( 129. * CF * CF * CF - 129. / 2. * CF * CF * Nc + 11413. / 54. * CF * Nc * Nc
            + CF * CF * nf * (-46. + 48. * zeta3) + CF * Nc * nf * (-556. / 27. - 48. * zeta3)
            - 70. / 27. * CF * nf * nf);
}

double QCD::threCorrForMass(const double nf_f, const double nf_i) const
{
    if (fabs(nf_f - nf_i) != 1.)
        throw std::runtime_error("Error in QCD::threCorrForMass()");

    double mu_threshold, mf, log_mu2_mf2;
    if (nf_f > nf_i) {
        if (nf_f == 6.) {
            mu_threshold = mut;
            mf = quarks[TOP].getMass(); // m_t(m_t)
        } else if (nf_f == 5.) {
            mu_threshold = mub;
            mf = quarks[BOTTOM].getMass(); // m_b(m_b)
        } else if (nf_f == 4.) {
            mu_threshold = muc;
            mf = quarks[CHARM].getMass(); // m_c(m_c)
        } else
            throw std::runtime_error("Error in QCD::threCorrForMass()");
        log_mu2_mf2 = 2. * log(mu_threshold / mf);
        return (1. + pow(Als(mu_threshold - MEPS, FULLNNLO) / M_PI, 2.)
                *(-log_mu2_mf2 * log_mu2_mf2 / 12. + 5. / 36. * log_mu2_mf2 - 89. / 432.));
    } else {
        if (nf_i == 6.) {
            mu_threshold = mut;
            mf = quarks[TOP].getMass(); // m_t(m_t)
        } else if (nf_i == 5.) {
            mu_threshold = mub;
            mf = quarks[BOTTOM].getMass(); // m_b(m_b)
        } else if (nf_i == 4.) {
            mu_threshold = muc;
            mf = quarks[CHARM].getMass(); // m_c(m_c)
        } else
            throw std::runtime_error("Error in QCD::threCorrForMass()");
        log_mu2_mf2 = 2. * log(mu_threshold / mf);
        return (1. + pow(Als(mu_threshold + MEPS, FULLNNLO) / M_PI, 2.)
                *(log_mu2_mf2 * log_mu2_mf2 / 12. - 5. / 36. * log_mu2_mf2 + 89. / 432.));
    }
}

double QCD::Mrun(const double mu, const double m, const orders order) const
{
    return Mrun(mu, m, m, order);
}

double QCD::Mrun(const double mu_f, const double mu_i, const double m,
        const orders order) const
{
    // Note: When the scale evolves across a flavour threshold, the definitions 
    //       of the outputs for "NLO" and "NNLO" become complicated. 

    int i;
    for (i = 0; i < CacheSize; ++i) {
        if ((mu_f == mrun_cache[0][i]) && (mu_i == mrun_cache[1][i]) &&
                (m == mrun_cache[2][i]) && ((double) order == mrun_cache[3][i]) &&
                (AlsM == mrun_cache[4][i]) && (MAls == mrun_cache[5][i]) &&
                (mut == mrun_cache[6][i]) && (mub == mrun_cache[7][i]) &&
                (muc == mrun_cache[8][i]))
            return mrun_cache[9][i];
    }

    double nfi = Nf(mu_i), nff = Nf(mu_f);
    double mu_threshold, mu_threshold2, mu_threshold3, mrun;
    if (nff == nfi)
        mrun = MrunTMP(mu_f, mu_i, m, order);
    else if (nff > nfi) {
        if (order == NLO || order == NNLO)
            throw std::runtime_error(orderToString(order) + " is not implemented in QCD::Mrun(mu_f,mu_i,m,order)");
        mu_threshold = AboveTh(mu_i);
        mrun = MrunTMP(mu_threshold - MEPS, mu_i, m, order);
        if (order == FULLNNLO)
            mrun *= threCorrForMass(nfi + 1., nfi); // threshold corrections
        if (nff == nfi + 1.) {
            mrun = MrunTMP(mu_f, mu_threshold + MEPS, mrun, order);
        } else if (nff == nfi + 2.) {
            mu_threshold2 = BelowTh(mu_f);
            mrun = MrunTMP(mu_threshold2 - MEPS, mu_threshold + MEPS, mrun, order);
            if (order == FULLNNLO)
                mrun *= threCorrForMass(nff, nfi + 1.); // threshold corrections
            mrun = MrunTMP(mu_f, mu_threshold2 + MEPS, mrun, order);
        } else if (nff == nfi + 3.) {
            mu_threshold2 = AboveTh(mu_threshold);
            mrun = MrunTMP(mu_threshold2 - MEPS, mu_threshold + MEPS, mrun, order);
            if (order == FULLNNLO)
                mrun *= threCorrForMass(nfi + 2., nfi + 1.); // threshold corrections
            mu_threshold3 = BelowTh(mu_f);
            mrun = MrunTMP(mu_threshold3 - MEPS, mu_threshold2 + MEPS, mrun, order);
            if (order == FULLNNLO)
                mrun *= threCorrForMass(nff, nfi + 2.); // threshold corrections
            mrun = MrunTMP(mu_f, mu_threshold3 + MEPS, mrun, order);
        } else
            throw std::runtime_error("Error in QCD::Mrun(mu_f,mu_i,m,order)");
    } else {
        if (order == NLO || order == NNLO)
            throw std::runtime_error(orderToString(order) + " is not implemented in QCD::Mrun(mu_f,mu_i,m,order)");
        mu_threshold = BelowTh(mu_i);
        mrun = MrunTMP(mu_threshold + MEPS, mu_i, m, order);
        if (order == FULLNNLO)
            mrun *= threCorrForMass(nfi - 1., nfi); // threshold corrections
        if (nff == nfi - 1.)
            mrun = MrunTMP(mu_f, mu_threshold - MEPS, mrun, order);
        else if (nff == nfi - 2.) {
            mu_threshold2 = AboveTh(mu_f);
            mrun = MrunTMP(mu_threshold2 + MEPS, mu_threshold - MEPS, mrun, order);
            if (order == FULLNNLO)
                mrun *= threCorrForMass(nff, nfi - 1.); // threshold corrections
            mrun = MrunTMP(mu_f, mu_threshold2 - MEPS, mrun, order);
        } else
            throw std::runtime_error("Error in QCD::Mrun(mu_f,mu_i,m,order)");
    }

    if (mrun < 0.0) {
        std::stringstream out;
        out << "QCD::Mrun(): A quark mass becomes tachyonic in QCD::Mrun("
                << mu_f << ", " << mu_i << ", " << m << ", " << orderToString(order) << ")"
                << std::endl
                << "             Als(" << mu_i << ", " << orderToString(order) << ")/(4pi)="
                << Als(mu_i, order) / (4. * M_PI) << std::endl
                << "             Als(" << mu_f << ", " << orderToString(order) << ")/(4pi)="
                << Als(mu_f, order) / (4. * M_PI);
        throw std::runtime_error(out.str());
    }

    CacheShift(mrun_cache, 10);
    mrun_cache[0][0] = mu_f;
    mrun_cache[1][0] = mu_i;
    mrun_cache[2][0] = m;
    mrun_cache[3][0] = (double) order;
    mrun_cache[4][0] = AlsM;
    mrun_cache[5][0] = MAls;
    mrun_cache[6][0] = mut;
    mrun_cache[7][0] = mub;
    mrun_cache[8][0] = muc;
    mrun_cache[9][0] = mrun;

    return mrun;
}

double QCD::MrunTMP(const double mu_f, const double mu_i, const double m,
        const orders order) const
{
    double nf = Nf(mu_f);
    if (nf != Nf(mu_i))
        throw std::runtime_error("Error in QCD::MrunTMP().");

    // alpha_s/(4pi)
    orders orderForAls;
    if (order == LO) orderForAls = LO;
    if (order == NLO || order == FULLNLO) orderForAls = FULLNLO;
    if (order == NNLO || order == FULLNNLO) orderForAls = FULLNNLO;
    double ai = Als(mu_i, orderForAls) / (4. * M_PI);
    double af = Als(mu_f, orderForAls) / (4. * M_PI);

    // LO contribution
    double b0 = Beta0(nf), g0 = Gamma0(nf);
    double mLO = m * pow(af / ai, g0 / (2. * b0));
    if (order == LO) return mLO;

    // NLO contribution
    double b1 = Beta1(nf), g1 = Gamma1(nf);
    double A1 = g1 / (2. * b0) - b1 * g0 / (2. * b0 * b0);
    double mNLO = mLO * A1 * (af - ai);
    if (order == NLO) return mNLO;
    if (order == FULLNLO) return (mLO + mNLO);

    // NNLO contribution    
    double b2 = Beta2(nf), g2 = Gamma2(nf);
    double A2 = b1 * b1 * g0 / (2. * b0 * b0 * b0) - b2 * g0 / (2. * b0 * b0) - b1 * g1 / (2. * b0 * b0) + g2 / (2. * b0);
    double mNNLO = mLO * (A1 * A1 / 2. * (af - ai)*(af - ai) + A2 / 2. * (af * af - ai * ai));
    if (order == NNLO) return mNNLO;
    if (order == FULLNNLO) return (mLO + mNLO + mNNLO);

    throw std::runtime_error(orderToString(order) + " is not implemented in QCD::MrunTMP()");
}

double QCD::Mrun4(const double mu_f, const double mu_i, const double m) const
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

double QCD::Mbar2Mp(const double mbar, const orders order) const
{
    // LO contribution
    double MpLO = mbar;
    if (order == LO) return MpLO;

    // alpha_s(mbar)/pi
    orders orderForAls;
    if (order == NLO || order == FULLNLO) orderForAls = FULLNLO;
    if (order == NNLO || order == FULLNNLO) orderForAls = FULLNNLO;
    double a = Als(mbar + MEPS, orderForAls) / M_PI;

    // NLO contribution 
    double MpNLO = mbar * 4. / 3. * a;
    if (order == NLO) return MpNLO;
    if (order == FULLNLO) return (MpLO + MpNLO);

    // NNLO contribution
    double nl, x;
    if (mbar < 3.)
        throw std::runtime_error("QCD::Mbar2Mp() can convert only top and bottom masses");
    else if (mbar < 50.) {
        // for the b quark
        nl = 4.;
        /* simply adding m_s(2 GeV) and m_c(m_c) */
        x = (quarks[STRANGE].getMass() + quarks[CHARM].getMass()) / mbar;
    } else {
        // for the top quark
        nl = 5.;
        /* simply adding m_s(2 GeV), m_c(m_c) and m_b(m_b) */
        x = (quarks[STRANGE].getMass() + quarks[CHARM].getMass()
                + quarks[BOTTOM].getMass()) / mbar;
    }
    double Delta = M_PI * M_PI / 8. * x - 0.597 * x * x + 0.230 * x * x*x;
    double MpNNLO = mbar * (307. / 32. + 2. * zeta2 + 2. / 3. * zeta2 * log(2.0) - zeta3 / 6.
            - nl / 3. * (zeta2 + 71. / 48.) + 4. / 3. * Delta) * a*a;
    if (order == NNLO) return MpNNLO;
    if (order == FULLNNLO) return (MpLO + MpNLO + MpNNLO);

    throw std::runtime_error(orderToString(order) + " is not implemented in QCD::Mbar2Mp().");
}

double QCD::Mp2MbarTMP(double *mu, double *params) const
{
    double mp = params[0];
    orders order = (orders) params[1];
    return (mp - Mbar2Mp(*mu, order));
}

double QCD::Mp2Mbar(const double mp, const orders order) const
{
    if (order == NLO || order == NNLO)
        throw std::runtime_error(orderToString(order) + " is not implemented in QCD::Mp2Mbar().");

    int i;
    double ms = quarks[STRANGE].getMass(), mc = quarks[CHARM].getMass();
    double alsmp = Als(mp, order);
    for (i = 0; i < CacheSize; ++i)
        if (alsmp == mp2mbar_cache[0][i] && ms == mp2mbar_cache[1][i] &&
                mc == mp2mbar_cache[2][i] && (double) order == mp2mbar_cache[3][i])
            return mp2mbar_cache[4][i];

    CacheShift(mp2mbar_cache, 5);
    mp2mbar_cache[0][0] = alsmp;
    mp2mbar_cache[1][0] = ms;
    mp2mbar_cache[2][0] = mc;
    mp2mbar_cache[3][0] = (double) order;

    TF1 f("f", this, &QCD::Mp2MbarTMP, mp / 2., 2. * mp, 2, "QCD", "mp2mbara");

    ROOT::Math::WrappedTF1 wf1(f);
    double params[2];
    params[0] = mp;
    params[1] = (double) order;
    wf1.SetParameters(params);

    ROOT::Math::BrentRootFinder brf;

    brf.SetFunction(wf1, .7 * mp, 1.3 * mp);
    if (brf.Solve())
        mp2mbar_cache[4][0] = brf.Root();
    else
        throw std::runtime_error("error in QCD::mp2mbar");

    return (mp2mbar_cache[4][0]);
}

double QCD::MS2DRqmass(const double MSbar) const
{
    return (MSbar / (1. + Als(MSbar, FULLNLO) / 4. / M_PI * CF));
}

double QCD::MS2DRqmass(const double MSscale, const double MSbar) const
{
    return (MSbar / (1. + Als(MSscale, FULLNLO) / 4. / M_PI * CF));
}
