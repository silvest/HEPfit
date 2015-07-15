/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "CPenguinBoxMu.h"

CPenguinBoxMu::CPenguinBoxMu(const StandardModel& model_i) 
: model(model_i), modelmatching(*(model_i.getMyMatching()))
{}

CPenguinBoxMu::~CPenguinBoxMu()
{}

double CPenguinBoxMu::C_NL()
{
    double muc = model.getMuc();
    double Mw =  model.getMuw();
    double mc = model.Mrun(muc, model.getQuarks(QCD::CHARM).getMass_scale(),
                           model.getQuarks(QCD::CHARM).getMass(), FULLNNLO);

    double K = model.Als(Mw)/model.Als(muc);
    double Kc = model.Als(muc)/model.Als(mc);
    double Kp = pow(K, 6./25.);
    double Km = pow(K, -12./25.);
    double K3 = pow(K, -1./25.);
    
    return( (mc*mc/Mw/Mw)/32.*pow(Kc,24./25.)*
            ((48./7.*Kp + 24./11.*Km - 696./77.*K3)*(4.*M_PI/model.Als(muc) + 15212./1875.*(1.-1./Kc))
            + (1.-log(muc*muc/mc/mc))*(16.*Kp-8.*Km) - 1176244./13125.*Kp - 2302./6875.*Km 
            + 3529184./48125.*K3 + K*(56248./4375.*Kp - 81448./6875.*Km + 
            4563698./144375*K3)));    
}

double CPenguinBoxMu::B_NL()
{
    double muc = model.getMuc();
    double Mw =  model.getMuw();
    double mc = model.Mrun(muc, model.getQuarks(QCD::CHARM).getMass_scale(),
                           model.getQuarks(QCD::CHARM).getMass(), FULLNNLO);

    double K = model.Als(Mw)/model.Als(muc);
    double Kc = model.Als(muc)/model.Als(mc);
    double K3 = pow(K, -1./25.);
    
    return( (mc*mc/Mw/Mw)*.25*pow(Kc,24./25.)*
            (3.*(1.-K3)*(4.*M_PI/model.Als(muc) + 15212./1875.*(1.-1./Kc))
            - -log(muc*muc/mc/mc) - 329./12. + 15212./625.*K3 + 30581./7500.*K*K3));
}

double CPenguinBoxMu::X_ch()
{
    
    double x = pow(model.Mrun(model.getMuw(), model.getQuarks(QCD::TOP).getMass_scale(), 
                              model.getQuarks(QCD::TOP).getMass(), FULLNNLO) / model.Mw_tree(), 2.);
    double a = model.computelamc().real()/model.getLambda()*(C_NL() - B_NL());
    double b = model.computelamt().real()/model.getLambda()*
              (modelmatching.Y0(x) + model.Als(model.getMuw())/4./M_PI*modelmatching.Y1(x, model.getMuw()));
    
    return( a*a + 2.*a*b );
}