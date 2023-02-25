/* 
 * Copyright (C) 2023 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Asl.h"
#include "StandardModelMatching.h"

#include "std_make_vector.h"

Asl::Asl(const StandardModel& SM_i, QCD::lepton lep_i)
: ThObservable(SM_i), mySM(SM_i) {
    lep = lep_i;
    //setParametersForObservable(make_vector<std::string>() << "Asl_test" );
}

Asl::~Asl() {
}

void Asl::updateParameters() {
    rhob = mySM.getCKM().getRhoBar();
    etab = mySM.getCKM().getEtaBar();
    Md = mySM.getQuarks(QCD::DOWN).getMass();
    Mb = mySM.getQuarks(QCD::BOTTOM).getMass();
    Mc = mySM.getQuarks(QCD::CHARM).getMass();
    Mt = mySM.getQuarks(QCD::TOP).getMass();
    MW = mySM.Mw();
    MB = mySM.getMesons(QCD::B_D).getMass();

    //wrong basis for now: Misiak instead of Buras
    allcoeff = mySM.getFlavour().ComputeCoeffBMll(Mb, lep);

    C_1 = (*(allcoeff[LO]))(0) + (*(allcoeff[NLO]))(0);
    C_2 = ((*(allcoeff[LO]))(1) + (*(allcoeff[NLO]))(1));
    C_3 = ((*(allcoeff[LO]))(2) + (*(allcoeff[NLO]))(2));
    C_4 = (*(allcoeff[LO]))(3) + (*(allcoeff[NLO]))(3);
    C_5 = ((*(allcoeff[LO]))(4) + (*(allcoeff[NLO]))(4));
    C_6 = ((*(allcoeff[LO]))(5) + (*(allcoeff[NLO]))(5));

//    std::cout.precision(3);
//    std::cout << "C_1" << C_1 << "\n";
//    std::cout << "C_2" << C_2 << "\n";
//    std::cout << "C_3" << C_3 << "\n";
//    std::cout << "C_4" << C_4 << "\n";
//    std::cout << "C_5" << C_5 << "\n";
//    std::cout << "C_6" << C_6 << "\n";

    K_1 = 3. * C_1 * C_1 + 2. * C_1 * C_2;
    K_2 = C_2 * C_2;
    K_1prime = 2. * (3. * C_1 * C_3 + C_1 * C_4 + C_2 * C_3);
    K_2prime = 2. * C_2 * C_4;
    K_3prime = 2 * (3 * C_1 * C_5 + C_1 * C_6 + C_2 * C_5 + C_2 * C_6);

//    std::cout << "K_1" << K_1 << "\n";
//    std::cout << "K_2" << K_2 << "\n";
//    std::cout << "K_1prime" << K_1prime << "\n";
//    std::cout << "K_2prime" << K_2prime << "\n";
//    std::cout << "K_3prime" << K_3prime << "\n";

    Mb2 = Mb * Mb;
    Mt2 = Mt * Mt;
    MB2 = MB * MB;
    MW2 = MW * MW;
    K12 = K_1 + K_2;

    z = Mc * Mc / Mb2;

    B = 0.87;
    B_s = 0.83;
    B_sprime = MB2 / pow(Mb + Md, 2) * B_s;
    
    eta_B = 0.55;
    eta_Bb = eta_B * pow(mySM.Als(MB), -6. / 23.) * (1. + mySM.Als(MB) / (4. * M_PI) * 5165. / 3174.);
    //Inami-lim: not correct
    S_0 = StandardModelMatching(mySM).S0(Mt2,MW2);
    //std::cout << S_0 << "\n";     //3023
    // S_0 = StandardModelMatching(mySM).S0(Mt,MW);
    //std::cout << S_0 << "\n";     //33
    
    kappa = 4. * M_PI * Mb2 / MW2 * K12 / (eta_Bb * S_0) * z;
    
    //StandardModel CP asymmetry in semileptonic B decay from hep-ph/0202010v2
    A_sl = -kappa * etab / ((1. - rhob) * (1. - rhob) + etab * etab)
            *
            (
            1. + z * (5. / 4. * B_sprime / B * (K_2 - K_1) / K12 - K_2 / K12)
            + z * z * (K_2 - K_1) / K12 * (1. / 3. - 5. / 6. * B_sprime / B)
            + 2. * z * z * (1. - etab) / ((1. - rhob) * (1. - rhob) + etab * etab) * (K_2 - K_1) / K12 * (5. / 2. * B_sprime / B - 1.)
            + z * (7. * K_1 + 3 * K_2) / (2. * K12 * B) * (MB2 - Mb2) / Mb2
            + (K_1prime + K_2prime - K_3prime) / K12
            );
    
    //etab = 0;
    std::cout << etab << "\n";
   
    return;
}
    
    double Asl::computeThValue(){
    //std::cout << SM.getOptionalParameter("Asl_test") << "\n";
    updateParameters();
    std::cout << "Asl: " << A_sl;
    return A_sl.real();
}