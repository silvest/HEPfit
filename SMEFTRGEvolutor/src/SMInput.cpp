//#include <gsl/gsl_linalg.h>

//#include "gsl/gsl_complex.h"
//#include "gsl/gsl_complex_math.h"

double RGESolver::GetCKMAngle(std::string name)
{
    return * (CKMAngles.at(name));
}

double RGESolver::GetCKMPhase()
{
    return CKM_delta;
}




//Extracts the 4 parameters from the CKM 
//affected by unphysical phases

void RGESolver::ExtractParametersFromCKM()
{

    s13 = CKM(0, 2).abs();
    c13 = sqrt(1. - s13 * s13);
    c12 = CKM(0, 0).abs() / c13;
    s12 = sqrt(1. - c12 * c12);
    s23 = CKM(1, 2).abs() / c13;
    c23 = sqrt(1. - s23 * s23);

    //See http://www.utfit.org/UTfit/Formalism
    //delta in [-Pi,Pi]

    if (
            s12 == 0. || c12 == 0. ||
            s13 == 0. || c13 == 0. ||
            s23 == 0. || c23 == 0.) {
        CKM_delta = 0.;
    } else {
        double gamma = (-(CKM(0, 0)*(CKM(0, 2)).conjugate()) /
                (CKM(1, 0)*(CKM(1, 2)).conjugate())).arg();

        double a = (c12 * s13 * s23) / (s12 * c23);

        double tan_g = tan(gamma);
        if (gamma < M_PI * 0.5 && gamma > -M_PI * 0.5) {
            CKM_delta = 2. * atan(
                    (1. - sqrt(1. - (a * a - 1.) * tan_g * tan_g)) /
                    ((a - 1.) * tan_g)
                    );
        } else {
            CKM_delta = 2. * atan(
                    (1. + sqrt(1. - (a * a - 1.) * tan_g * tan_g)) /
                    ((a - 1.) * tan_g)
                    );
        }
    }


}

void RGESolver::FromMassesToYukawas(std::string basis)
{
    double v = sqrt(0.5 * mh2 / lambda);
    UpdateCKM();

    gslpp::vector<double> yudiag(3), yddiag(3), yediag(3);
    double sqrt2ov = sqrt(2.) / v;

    yudiag(0) = sqrt2ov * mu;
    yudiag(1) = sqrt2ov * mc;
    yudiag(2) = sqrt2ov * mt;

    yddiag(0) = sqrt2ov * md;
    yddiag(1) = sqrt2ov * ms;
    yddiag(2) = sqrt2ov * mb;

    yediag(0) = sqrt2ov * mel;
    yediag(1) = sqrt2ov * mmu;
    yediag(2) = sqrt2ov * mtau;

    yuR = gslpp::matrix<double>(yudiag);
    yuI.reset();
    yeR = gslpp::matrix<double>(yediag);
    yeI.reset();
    ydR = gslpp::matrix<double>(yddiag);
    ydI.reset();

    if (basis == "UP") {

        gslpp::matrix<gslpp::complex> yd = ydR * CKM.hconjugate();
        ydR = yd.real();
        ydI = yd.imag();

    } else if (basis == "DOWN") {

        gslpp::matrix<gslpp::complex> yu = yuR * CKM;
        yuR = yu.real();
        yuI = yu.imag();

    } else
        std::cout << "WARNING: wrong basis choice: " << basis << ": Yukawa couplings not updated!!!" << std::endl;

}

void RGESolver::UpdateCKM()
{
    gslpp::complex PhaseFactor(1., CKM_delta,
            true);
    //First row
    CKM.assign(0, 0, c12 * c13);
    CKM.assign(0, 1, s12 * c13);
    CKM.assign(0, 2, s13 * PhaseFactor.conjugate());

    //Second row 
    CKM.assign(1, 0, -s12 * c23 - c12 * s23 * s13 * PhaseFactor);
    CKM.assign(1, 1, c12 * c23 - s12 * s23 * s13 * PhaseFactor);
    CKM.assign(1, 2, s23 * c13);

    //Third row
    CKM.assign(2, 0, s12 * s23 - c12 * c23 * s13 * PhaseFactor);
    CKM.assign(2, 1, -c12 * s23 - s12 * c23 * s13 * PhaseFactor);
    CKM.assign(2, 2, c23 * c13);
}

void RGESolver::GenerateSMInitialConditions(
        double mu, std::string basis,
        std::string method)
{
    /*if (method != "Numeric" && method != "Approximate") {
        std::cout << "WARNING : invalid method\n"
                "Available methods: Numeric, Approximate"
                << std::endl;
    }*/
    SetSMDefaultInput();
    //Before evolving, eventual changes 
    //in the input values of CKM angles must be 
    //translated 
    c12 = sqrt(1. - s12 * s12);
    c13 = sqrt(1. - s13 * s13);
    c23 = sqrt(1. - s23 * s23);

    //if (inputCKM == true) {
    UpdateCKM();
    FromMassesToYukawas(basis);
    // }

    EvolveSMOnly(method, InputScale_SM, mu);
    GoToBasisSMOnly(basis);
}

void RGESolver::GenerateSMInitialConditions(double muIn, double muFin, std::string basis, std::string method,
        double g1in, double g2in, double g3in, double lambdain, double mh2in,
        double Muin[3], double Mdin[3], double Mein[3],
        double s12in, double s13in, double s23in, double deltain)
{

        g1 = g1in;
        g2 = g2in;
        g3 = g3in;

        mh2 = mh2in;
        lambda = lambdain;

        mu = Muin[0];
        mc = Muin[1];
        mt = Muin[2];

        md = Mdin[0];
        ms = Mdin[1];
        mb = Mdin[2];

        mel = Mein[0];
        mmu = Mein[1];
        mtau = Mein[2];

        s12 = s12in;
        s13 = s13in;
        s23 = s23in;
        CKM_delta = deltain;


        c12 = sqrt(1. - s12 * s12);
        c13 = sqrt(1. - s13 * s13);
        c23 = sqrt(1. - s23 * s23);

        UpdateCKM();
        FromMassesToYukawas(basis);

        EvolveSMOnly(method, muIn, muFin);
        GoToBasisSMOnly(basis);
    
}

void RGESolver::EvolveSMOnly(std::string method, double muI, double muF)
{
    if (method != "Numeric" && method != "Approximate") {
        std::cout << "WARNING : invalid method\n"
                "Available methods: Numeric, Approximate"
                << std::endl;
    }


    if (muF != muI) {
        //Initial conditions are inserted 
        //in the array x
        InitSMOnly();

        //Numeric solution 
        if (method == "Numeric") {
            double tI = log(muI);
            double tF = log(muF);
            double ttmp = tI;
            gsl_odeiv2_control * cSMOnly
                    = gsl_odeiv2_control_y_new(epsabs_, epsrel_);

            //Initial step 
            double stepIn = (tF - tI)*0.01;
            if (tF > tI) {
                while (ttmp < tF) {
                    int status = gsl_odeiv2_evolve_apply(eSMOnly, cSMOnly, sSMOnly,
                            &sysSMOnly, &ttmp, tF, &stepIn, x);
                    //std::cout << "step : " << stepIn << std::endl;
                    if (status != GSL_SUCCESS) {
                        printf("error, return value=%d\n", status);
                        break;
                    }
                }
            }
            if (tF < tI) {
                while (ttmp > tF) {
                    int status = gsl_odeiv2_evolve_apply(eSMOnly, cSMOnly, sSMOnly,
                            &sysSMOnly, &ttmp, tF, &stepIn, x);
                    if (status != GSL_SUCCESS) {
                        printf("error, return value=%d\n", status);
                        break;
                    }
                }
            }

            gsl_odeiv2_evolve_reset(eSMOnly);
            gsl_odeiv2_step_reset(sSMOnly);
            gsl_odeiv2_control_free(cSMOnly);
        }


        //Approximate solution 
        //-------------------------------------------
        if (method == "Approximate") {
            double beta[59] = {0.};
            double Log_muF_over_muI = log(muF / muI);
            /*int status = */funcSMOnly(10., x, beta, NULL);
            for (int i = 0; i < 59; i++) {
                x[i] += beta[i] * Log_muF_over_muI;
            }
        }
        //-------------------------------------------

        UpdateSMOnly();
        //Evolved values from x are put 
        //back in the coefficients 
    }
}

void RGESolver::GoToBasisSMOnly(std::string basis)
{

    gslpp::matrix<gslpp::complex> Uu(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Vu(3, 3, 0.);
    gslpp::matrix<gslpp::complex> yuDiag(3, 3, 0.);
    gslpp::vector<double> Su(3, 0.);

    gslpp::matrix<gslpp::complex> Ud(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Vd(3, 3, 0.);
    gslpp::matrix<gslpp::complex> ydDiag(3, 3, 0.);
    gslpp::vector<double> Sd(3, 0.);

    gslpp::matrix<gslpp::complex> Re(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rl(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Redag(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rldag(3, 3, 0.);
    gslpp::matrix<gslpp::complex> yeDiag(3, 3, 0.);
    gslpp::vector<double> Se(3, 0.);

    gslpp::matrix<gslpp::complex> yu(3, 3, 0.);
    gslpp::matrix<gslpp::complex> yd(3, 3, 0.);
    gslpp::matrix<gslpp::complex> ye(3, 3, 0.);

    yu = yuR + gslpp::complex::i() * yuI;
    yd = ydR + gslpp::complex::i() * ydI;
    ye = yeR + gslpp::complex::i() * yeI;

    yu.singularvalue(Uu, Vu, Su);
    yd.singularvalue(Ud, Vd, Sd);
    ye.singularvalue(Re, Rl, Se);




    //Matrix to rotate fields
    gslpp::matrix<gslpp::complex> Ru(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rd(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rq(3, 3, 0.);

    //Computing the CKM 
    CKM = (Vu.hconjugate()) * Vd;
    //Extract the 4 parameters from the raw CKM
    ExtractParametersFromCKM();
    //Build the CKM with the 4 parameters
    UpdateCKM();

    yuR = gslpp::matrix<double>(Su);
    yuI.reset();
    yeR = gslpp::matrix<double>(Se);
    yeI.reset();
    ydR = gslpp::matrix<double>(Sd);
    ydI.reset();

    if (basis == "UP") {

        Ru = Uu;
        Rd = Ud;
        Rq = Vu;

        gslpp::matrix<gslpp::complex> yd = ydR * CKM.hconjugate();
        ydR = yd.real();
        ydI = yd.imag();

    } else if (basis == "DOWN") {

        Ru = Uu;
        Rd = Ud;
        Rq = Vd;

        gslpp::matrix<gslpp::complex> yu = yuR * CKM; // CKM*yuR;
        yuR = yu.real();
        yuI = yu.imag();

    } else
        std::cout << "WARNING: wrong basis choice: " << basis << ": Yukawa couplings not updated!!!" << std::endl;


}

void RGESolver::InitSMOnly()
{

    int n = 0;
    int i, j;
    x[0] = g2;
    x[1] = g1;
    x[2] = g3;
    n += 3;
    x[n] = lambda;
    n++;
    x[n] = mh2;
    n++;

    for (i = 0; i < NG; i++) {
        for (j = 0; j < NG; j++) {
            int count = 0;
            x[n + count * DF] = yuR(i, j);
            count++;
            x[n + count * DF] = yuI(i, j);
            count++;
            x[n + count * DF] = ydR(i, j);
            count++;
            x[n + count * DF] = ydI(i, j);
            count++;
            x[n + count * DF] = yeR(i, j);
            count++;
            x[n + count * DF] = yeI(i, j);
            count++;

            n++;
        }
    }
    n += (2. * Nyukawa - 1.) * DF;

}

void RGESolver::UpdateSMOnly()
{
    //Here we return to original structures

    int n = 0;
    int a, i, j;
    g2 = x[0];
    g1 = x[1];
    g3 = x[2];
    n += 3;
    lambda = x[n];
    n++;
    mh2 = x[n];
    n++;

    for (i = 0; i < NG; i++) {
        for (j = 0; j < NG; j++) {
            a = 0;
            yuR.assign(i, j, x[n + a * DF]);
            a++;
            yuI.assign(i, j, x[n + a * DF]);
            a++;
            ydR.assign(i, j, x[n + a * DF]);
            a++;
            ydI.assign(i, j, x[n + a * DF]);
            a++;
            yeR.assign(i, j, x[n + a * DF]);
            a++;
            yeI.assign(i, j, x[n + a * DF]);
            a++;
            n++;
        }
    }
    n += (2. * Nyukawa - 1.) * DF;
}

int RGESolver::funcSMOnly(double logmu, const double y[], double f[], void* params)
{
    int i, j, a, b;
    double loop_factor = 1.
            / (16. * M_PI * M_PI);
    int c = 0;
    //counter initialized at 0. This index reads inside the array y, where all 
    //independent parameters are stored. In this part of the function the 
    //variables are organized in the correct structures.

    double g1 = y[c + 1]; //gauge couplings 
    double g2 = y[c ];
    double g3 = y[c + 2];

    c += Ngauge;

    double lambda = y[c]; //Higgs sector
    double mh2 = y[c + 1];
    c += Nh;

    gslpp::matrix<double> yeR(NG, 0.), yeI(NG, 0.), //yukawas
            ydR(NG, 0.), ydI(NG, 0.),
            yuR(NG, 0.), yuI(NG, 0.);
    gslpp::matrix<double> yedagR(NG, 0.), yedagI(NG, 0.),
            yddagR(NG, 0.), yddagI(NG, 0.),
            yudagR(NG, 0.), yudagI(NG, 0.);

    for (i = 0; i < NG; i++) {
        for (j = 0; j < NG; j++) {
            a = 0;
            yuR.assign(i, j, y[c]);
            yudagR.assign(j, i, y[c]);
            a++;
            yuI.assign(i, j, y[c + DF]);
            yudagI.assign(j, i, -y[c + DF]);
            a++;
            ydR.assign(i, j, y[c + 2 * DF]);
            yddagR.assign(j, i, y[c + 2 * DF]);
            a++;
            ydI.assign(i, j, y[c + 3 * DF]);
            yddagI.assign(j, i, -y[c + 3 * DF]);
            a++;
            yeR.assign(i, j, y[c + 4 * DF]);
            yedagR.assign(j, i, y[c + 4 * DF]);
            a++;
            yeI.assign(i, j, y[c + 5 * DF]);
            yedagI.assign(j, i, -y[c + 5 * DF]);
            a++;
            c++;
        }
    }
    c += (2. * Nyukawa - 1.) * DF;

    //Auxiliary quantities: 

    //Powers and products of gauge couplings
    double g12 = g1*g1;
    double g22 = g2*g2;
    double g32 = g3*g3;
    //double g1g3 = g1*g3;
    //double g1g2 = g1*g2;
    //double g2g3 = g2*g3;

    double gammaH = 0.; //Higgs wavefunction normalization  
    double H = 0.; //scalar yukawa-trace dependent appearing in lambda RGE
    gslpp::matrix<double> yudyuR(NG, 0.), yudyuI(NG, 0.); //yu^dag yu 
    gslpp::matrix<double> yddydR(NG, 0.), yddydI(NG, 0.); //yd^dag yd 
    gslpp::matrix<double> yedyeR(NG, 0.), yedyeI(NG, 0.); //ye^dag ye 
    gslpp::matrix<double> ydyudR(NG, 0.), ydyudI(NG, 0.); //yd yu^dag

    yudyuR = yudagR * yuR - yudagI * yuI;
    yudyuI = yudagI * yuR + yudagR * yuI;
    yddydR = yddagR * ydR - yddagI * ydI;
    yddydI = yddagI * ydR + yddagR * ydI;
    yedyeR = yedagR * yeR - yedagI * yeI;
    yedyeI = yedagI * yeR + yedagR * yeI;
    ydyudR = ydR * yudagR - ydI * yudagI;
    ydyudI = ydI * yudagR + ydR * yudagI;





    //gammaH
    for (i = 0; i < NG; i++) {
        gammaH += yedyeR(i, i) + NC * (yudyuR(i, i) + yddydR(i, i));
        for (j = 0; j < NG; j++) {
            H += yedyeR(i, j) * yedyeR(j, i) - yedyeI(i, j) * yedyeI(j, i)
                    + NC * (yudyuR(i, j) * yudyuR(j, i) - yudyuI(i, j) * yudyuI(j, i)
                    + yddydR(i, j) * yddydR(j, i) - yddydI(i, j) * yddydI(j, i));
        }
    }



    //------------------------------------
    //---------------RGE------------------
    //------------------------------------

    c = 0; //counter restarts from 0 

    //---------RGE GAUGE/HIGGS --------

    //SMEFT contributes to SM beta functions are proportional to mh2.
    //They are in RGE 1.
    //SM beta functions for mh2 and lambda are in 
    //https://arxiv.org/abs/hep-ph/0207271 
    //mh2,lambda and Yukawas follow the conventions of RGE 1

    {
        //g2
        f[c] = (-b02 * g22 //SM
                ) * g2 * loop_factor;
        c++;
        //g1
        f[c] = (-b01 * g12) * g1 * loop_factor;
        c++;
        //g3
        f[c] = (-b03 * g32) * g3 * loop_factor;
        c++;

        //lambda
        f[c] = lambda * (24. * lambda - 3. * g12 - 9. * g22 + 4. * gammaH)
                + 0.375 * g12 * g12 + 0.75 * g12 * g22
                + 1.125 * g22 * g22 - 2. * H //SM
                ;
        f[c] *= loop_factor;
        c++;
        //mh2
        f[c] = mh2 * (12. * lambda + 2. * gammaH - 1.5 * g12 - 4.5 * g22 //SM  
                ) * loop_factor;
        c++;
    }

    //---------RGE YUKAWA --------
    //SM beta functions can be found in
    // https://arxiv.org/abs/hep-ph/0207271
    for (i = 0; i < NG; i++) {
        for (j = 0; j < NG; j++) {
            //Entries without matrix products:      
            //yuR  
            f[c ] = (gammaH - (17. / 12.) * g12 - 2.25 * g22 - 8. * g32) * yuR(i, j)
                    ;
            //yuI
            f[c + DF] = (gammaH - (17. / 12.) * g12 - 2.25 * g22 - 8. * g32) * yuI(i, j)
                    ;
            //ydR     
            f[c + 2 * DF] = (gammaH - (5. / 12.) * g12 - 2.25 * g22
                    - 8. * g32) * ydR(i, j)
                    ;
            //ydI  
            f[c + 3 * DF] = (gammaH - (5. / 12.) * g12 - 2.25 * g22
                    - 8. * g32) * ydI(i, j)
                    ;
            //yeR
            f[c + 4 * DF] = (gammaH - 2.25 * g22 - 3.75 * g12) * yeR(i, j)
                    ;
            //yeI
            f[c + 5 * DF] = (gammaH - 2.25 * g22 - 3.75 * g12) * yeI(i, j)
                    ;
            //Entries with 1 matrix product (1 summed index)
            for (b = 0; b < NG; b++) {
                //yuR
                f[c ] +=
                        +1.5 * (yuR(i, b)*(yudyuR(b, j) - yddydR(b, j))
                        - yuI(i, b)*(yudyuI(b, j) - yddydI(b, j)))
                        ;
                //yuI
                f[c + DF] += 1.5 * (yuI(i, b)*(yudyuR(b, j) - yddydR(b, j))
                        + yuR(i, b)*(yudyuI(b, j) - yddydI(b, j)))
                        ;
                //ydR       
                f[c + 2 * DF] += 1.5 * (ydR(i, b)*(yddydR(b, j) - yudyuR(b, j))
                        - ydI(i, b)*(yddydI(b, j) - yudyuI(b, j)))
                        ;
                //ydI 
                f[c + 3 * DF] += 1.5 * (
                        ydI(i, b)*(yddydR(b, j) - yudyuR(b, j))
                        + ydR(i, b)*(yddydI(b, j) - yudyuI(b, j)))
                        ;
                //yeR
                f[c + 4 * DF] += 1.5 * (yeR(i, b)*(yedyeR(b, j))
                        - yeI(i, b)*(yedyeI(b, j)))
                        ;
                //yeI 
                f[c + 5 * DF] += 1.5 * (yeI(i, b)*(yedyeR(b, j))
                        + yeR(i, b)*(yedyeI(b, j)))
                        ;

            }


            for (a = 0; a < (2 * Nyukawa); a++) {
                f[c + a * DF] *= loop_factor;
            }
            c++;
        }
    }
    c += (2 * Nyukawa - 1) * DF;
    return GSL_SUCCESS;
}


//Default input for SM 

void RGESolver::SetSMDefaultInput()
{
    //Default input scale for SM in GeV
    //Higgs and Gauge sector


    /*
    mh2 = 126. * 126.;
    lambda = 0.14065;
    g3 = 1.22;
    g2 = 0.6516;
    g1 = .3576;*/

    //Inputs in gauge sector
    double alpha_em = 0.007812;
    double alpha_s = 0.1074;
    double sinWeakSquared = 0.23147;

    double e = sqrt(4. * M_PI * alpha_em);

    g1 = e / sqrt(1. - sinWeakSquared);
    g2 = e / sqrt(sinWeakSquared);
    g3 = sqrt(4. * M_PI * alpha_s);


    //Inputs in the Higgs sector (to be updated)
    mh2 = 125.10 * 125.10;
    //double v = 247.;
    double mZ = 91.191; //Z mass in GeV
    double v = mZ * 2. * sqrt(1. - sinWeakSquared) / g2;
    lambda = 0.5 * mh2 / (v * v);



    //Fermion masses in GeV
    mu = 0.0012;
    mc = 0.64;
    mt = 162.;

    md = 0.0027;
    ms = 0.052;
    mb = 2.75;
    //Not evolved via RGEs, equal to their experimental values (PDG)
    mel = 0.0005110;
    mmu = 0.1057;
    mtau = 1.776;


    //CKM parameters in radians 
    s12 = 0.225;
    s13 = 0.042;
    s23 = 0.003675;
    CKM_delta = 1.167625;
    c12 = sqrt(1. - s12 * s12);
    c13 = sqrt(1. - s13 * s13);
    c23 = sqrt(1. - s23 * s23);



    //By default, Yukawas are aligned with CKM input
    //in up basis
    FromMassesToYukawas("UP");

    InputScale_SM = v / sqrt(2.);
   

}

void RGESolver::Reset()
{

    //Resets integration parameters to default 
    Resetepsabs();
    Resetepsrel();


    int a, i, j, k, l;
    /*g2 = 0;
    g1 = 0;
    g3 = 0;
    lambda = 0;
    mh2 = 0;

    yuR.reset();
    yuI.reset();
    ydR.reset();
    ydI.reset();
    yeR.reset();
    yeI.reset();*/
    SetSMDefaultInput();

    //    for (i = 0; i < NG; i ++) {
    //        for (j = 0; j < NG; j ++) {
    //            a = 0;
    //            yuR.assign(i,j,0);
    //            a ++;
    //            yuI.assign(i,j,0);
    //            a ++;
    //            ydR.assign(i,j,0);
    //            a ++;
    //            ydI.assign(i,j,0);
    //            a ++;
    //            yeR.assign(i,j,0);
    //            a ++;
    //            yeI.assign(i,j,0);
    //            a ++;
    //            //        }
    //    }

    cG = 0;
    cGT = 0;
    cW = 0;
    cWT = 0;
    cH = 0;
    cHBOX = 0;
    cHD = 0;

    cHG = 0;
    cHB = 0;
    cHW = 0;
    cHWB = 0;
    cHGT = 0;
    cHBT = 0;
    cHWT = 0;
    cHWBT = 0;

    //class 5
    for (i = 0; i < NG; i++) {
        for (j = 0; j < NG; j++) {
            a = 0;
            WC1_set(cuHR, i, j, 0);
            WC1_set(cuHI, i, j, 0);
            a++;
            WC1_set(cdHR, i, j, 0);
            WC1_set(cdHI, i, j, 0);
            a++;
            WC1_set(ceHR, i, j, 0);
            WC1_set(ceHI, i, j, 0);
            a++;
        }
    }

    //Class 6
    for (i = 0; i < NG; i++) {
        for (j = 0; j < NG; j++) {
            a = 0;
            WC1_set(ceWR, i, j, 0);
            a++;
            WC1_set(ceWI, i, j, 0);
            a++;
            WC1_set(ceBR, i, j, 0);
            a++;
            WC1_set(ceBI, i, j, 0);
            a++;
            WC1_set(cuGR, i, j, 0);
            a++;
            WC1_set(cuGI, i, j, 0);
            a++;
            WC1_set(cuWR, i, j, 0);
            a++;
            WC1_set(cuWI, i, j, 0);
            a++;
            WC1_set(cuBR, i, j, 0);
            a++;
            WC1_set(cuBI, i, j, 0);
            a++;
            WC1_set(cdGR, i, j, 0);
            a++;
            WC1_set(cdGI, i, j, 0);
            a++;
            WC1_set(cdWR, i, j, 0);
            a++;
            WC1_set(cdWI, i, j, 0);
            a++;
            WC1_set(cdBR, i, j, 0);
            a++;
            WC1_set(cdBI, i, j, 0);
            a++;
        }
    }

    //class 7
    {
        for (i = 0; i < DWC2R; i++) {
            WC2R_set(cHl1R, WC2R_indices[i][0], WC2R_indices[i][1], 0);
        }
        for (i = 0; i < DWC2I; i++) {
            WC2I_set(cHl1I, WC2I_indices[i][0], WC2I_indices[i][1], 0);
        }
        for (i = 0; i < DWC2R; i++) {
            WC2R_set(cHl3R, WC2R_indices[i][0], WC2R_indices[i][1], 0);
        }
        for (i = 0; i < DWC2I; i++) {
            WC2I_set(cHl3I, WC2I_indices[i][0], WC2I_indices[i][1], 0);
        }
        for (i = 0; i < DWC2R; i++) {
            WC2R_set(cHeR, WC2R_indices[i][0], WC2R_indices[i][1], 0);
        }
        for (i = 0; i < DWC2I; i++) {
            WC2I_set(cHeI, WC2I_indices[i][0], WC2I_indices[i][1], 0);
        }

        for (i = 0; i < DWC2R; i++) {
            WC2R_set(cHq1R, WC2R_indices[i][0], WC2R_indices[i][1], 0);
        }
        for (i = 0; i < DWC2I; i++) {
            WC2I_set(cHq1I, WC2I_indices[i][0], WC2I_indices[i][1], 0);
        }
        for (i = 0; i < DWC2R; i++) {
            WC2R_set(cHq3R, WC2R_indices[i][0], WC2R_indices[i][1], 0);
        }
        for (i = 0; i < DWC2I; i++) {
            WC2I_set(cHq3I, WC2I_indices[i][0], WC2I_indices[i][1], 0);
        }

        for (i = 0; i < DWC2R; i++) {
            WC2R_set(cHuR, WC2R_indices[i][0], WC2R_indices[i][1], 0);
        }
        for (i = 0; i < DWC2I; i++) {
            WC2I_set(cHuI, WC2I_indices[i][0], WC2I_indices[i][1], 0);
        }

        for (i = 0; i < DWC2R; i++) {
            WC2R_set(cHdR, WC2R_indices[i][0], WC2R_indices[i][1], 0);
        }
        for (i = 0; i < DWC2I; i++) {
            WC2I_set(cHdI, WC2I_indices[i][0], WC2I_indices[i][1], 0);
        }

        for (i = 0; i < NG; i++) {
            for (j = 0; j < NG; j++) {
                WC1_set(cHudR, i, j, 0);
            }
        }
        for (i = 0; i < NG; i++) {
            for (j = 0; j < NG; j++) {
                WC1_set(cHudI, i, j, 0);
            }
        }


    }
    //class 8_LLLL
    {
        for (a = 0; a < DWC6R; a++) {
            WC6R_set(cllR, WC6R_indices[a][0], WC6R_indices[a][1],
                    WC6R_indices[a][2], WC6R_indices[a][3], 0);
        }
        for (a = 0; a < DWC6I; a++) {
            WC6I_set(cllI, WC6I_indices[a][0], WC6I_indices[a][1],
                    WC6I_indices[a][2], WC6I_indices[a][3], 0);
        }
        for (a = 0; a < DWC6R; a++) {
            WC6R_set(cqq1R, WC6R_indices[a][0], WC6R_indices[a][1],
                    WC6R_indices[a][2], WC6R_indices[a][3], 0);
        }
        for (a = 0; a < DWC6I; a++) {
            WC6I_set(cqq1I, WC6I_indices[a][0], WC6I_indices[a][1],
                    WC6I_indices[a][2], WC6I_indices[a][3], 0);
        }
        for (a = 0; a < DWC6R; a++) {
            WC6R_set(cqq3R, WC6R_indices[a][0], WC6R_indices[a][1],
                    WC6R_indices[a][2], WC6R_indices[a][3], 0);
        }
        for (a = 0; a < DWC6I; a++) {
            WC6I_set(cqq3I, WC6I_indices[a][0], WC6I_indices[a][1],
                    WC6I_indices[a][2], WC6I_indices[a][3], 0);
        }

        for (a = 0; a < DWC7R; a++) {
            WC7R_set(clq1R, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
        }
        for (a = 0; a < DWC7I; a++) {
            WC7I_set(clq1I, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
        }
        for (a = 0; a < DWC7R; a++) {
            WC7R_set(clq3R, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
        }
        for (a = 0; a < DWC7I; a++) {
            WC7I_set(clq3I, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
        }

    }
    //Class 8_RRRR
    {
        for (a = 0; a < DWC8R; a++) {
            WC8R_set(ceeR, WC8R_indices[a][0], WC8R_indices[a][1],
                    WC8R_indices[a][2], WC8R_indices[a][3], 0);
        }
        for (a = 0; a < DWC8I; a++) {
            WC8I_set(ceeI, WC8I_indices[a][0], WC8I_indices[a][1],
                    WC8I_indices[a][2], WC8I_indices[a][3], 0);
        }
        for (a = 0; a < DWC6R; a++) {
            WC6R_set(cuuR, WC6R_indices[a][0], WC6R_indices[a][1],
                    WC6R_indices[a][2], WC6R_indices[a][3], 0);
        }
        for (a = 0; a < DWC6I; a++) {
            WC6I_set(cuuI, WC6I_indices[a][0], WC6I_indices[a][1],
                    WC6I_indices[a][2], WC6I_indices[a][3], 0);
        }
        for (a = 0; a < DWC6R; a++) {
            WC6R_set(cddR, WC6R_indices[a][0], WC6R_indices[a][1],
                    WC6R_indices[a][2], WC6R_indices[a][3], 0);
        }
        for (a = 0; a < DWC6I; a++) {
            WC6I_set(cddI, WC6I_indices[a][0], WC6I_indices[a][1],
                    WC6I_indices[a][2], WC6I_indices[a][3], 0);
        }
        for (a = 0; a < DWC7R; a++) {
            WC7R_set(ceuR, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
        }
        for (a = 0; a < DWC7I; a++) {
            WC7I_set(ceuI, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
        }


        for (a = 0; a < DWC7R; a++) {
            WC7R_set(cedR, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
        }
        for (a = 0; a < DWC7I; a++) {
            WC7I_set(cedI, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
        }

        for (a = 0; a < DWC7R; a++) {
            WC7R_set(cud1R, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
        }
        for (a = 0; a < DWC7I; a++) {
            WC7I_set(cud1I, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
        }
        for (a = 0; a < DWC7R; a++) {
            WC7R_set(cud8R, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
        }
        for (a = 0; a < DWC7I; a++) {
            WC7I_set(cud8I, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
        }


    }

    //Class 8_LLRR
    {
        for (a = 0; a < DWC7R; a++) {
            WC7R_set(cleR, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
        }
        for (a = 0; a < DWC7I; a++) {
            WC7I_set(cleI, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
        }
        for (a = 0; a < DWC7R; a++) {
            WC7R_set(cluR, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
        }
        for (a = 0; a < DWC7I; a++) {
            WC7I_set(cluI, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
        }

        for (a = 0; a < DWC7R; a++) {
            WC7R_set(cldR, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
        }
        for (a = 0; a < DWC7I; a++) {
            WC7I_set(cldI, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
        }
        for (a = 0; a < DWC7R; a++) {
            WC7R_set(cqeR, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
        }
        for (a = 0; a < DWC7I; a++) {
            WC7I_set(cqeI, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
        }

        for (a = 0; a < DWC7R; a++) {
            WC7R_set(cqu1R, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
        }
        for (a = 0; a < DWC7I; a++) {
            WC7I_set(cqu1I, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
        }
        for (a = 0; a < DWC7R; a++) {
            WC7R_set(cqu8R, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
        }
        for (a = 0; a < DWC7I; a++) {
            WC7I_set(cqu8I, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
        }



        for (a = 0; a < DWC7R; a++) {
            WC7R_set(cqd1R, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
        }
        for (a = 0; a < DWC7I; a++) {
            WC7I_set(cqd1I, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
        }
        for (a = 0; a < DWC7R; a++) {
            WC7R_set(cqd8R, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
        }
        for (a = 0; a < DWC7I; a++) {
            WC7I_set(cqd8I, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
        }

    }

    //Class 8_LRRL

    for (i = 0; i < NG; i++) {
        for (j = 0; j < NG; j++) {
            for (k = 0; k < NG; k++) {
                for (l = 0; l < NG; l++) {
                    WC5_set(cledqR, i, j, k, l, 0);
                    WC5_set(cledqI, i, j, k, l, 0);
                }
            }
        }
    }

    //Class 8_LRLR
    for (i = 0; i < NG; i++) {
        for (j = 0; j < NG; j++) {
            for (k = 0; k < NG; k++) {
                for (l = 0; l < NG; l++) {
                    a = 0;
                    WC5_set(cquqd1R, i, j, k, l, 0);
                    a++;
                    WC5_set(cquqd1I, i, j, k, l, 0);
                    a++;
                    WC5_set(cquqd8R, i, j, k, l, 0);
                    a++;
                    WC5_set(cquqd8I, i, j, k, l, 0);
                    a++;
                    WC5_set(clequ1R, i, j, k, l, 0);
                    a++;
                    WC5_set(clequ1I, i, j, k, l, 0);
                    a++;
                    WC5_set(clequ3R, i, j, k, l, 0);
                    a++;
                    WC5_set(clequ3I, i, j, k, l, 0);
                    a++;
                }
            }
        }
    }

}

void RGESolver::GoToBasis(std::string basis)
{

    gslpp::matrix<gslpp::complex> Uu(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Vu(3, 3, 0.);
    gslpp::matrix<gslpp::complex> yuDiag(3, 3, 0.);
    gslpp::vector<double> Su(3, 0.);

    gslpp::matrix<gslpp::complex> Ud(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Vd(3, 3, 0.);
    gslpp::matrix<gslpp::complex> ydDiag(3, 3, 0.);
    gslpp::vector<double> Sd(3, 0.);

    gslpp::matrix<gslpp::complex> Re(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rl(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Redag(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rldag(3, 3, 0.);
    gslpp::matrix<gslpp::complex> yeDiag(3, 3, 0.);
    gslpp::vector<double> Se(3, 0.);

    int i, j;
    gslpp::matrix<gslpp::complex> yu(3, 3, 0.);
    gslpp::matrix<gslpp::complex> yd(3, 3, 0.);
    gslpp::matrix<gslpp::complex> ye(3, 3, 0.);

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            yu.assign(i, j, gslpp::complex(yuR(i, j), yuI(i, j), false));
            yd.assign(i, j, gslpp::complex(ydR(i, j), ydI(i, j), false));
            ye.assign(i, j, gslpp::complex(yeR(i, j), yeI(i, j), false));
        }
    }

    using namespace std;

    yu.singularvalue(Uu, Vu, Su);
    yd.singularvalue(Ud, Vd, Sd);
    ye.singularvalue(Re, Rl, Se);


    //Matrix to rotate fields
    gslpp::matrix<gslpp::complex> Ru(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rudag(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rd(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rddag(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rq(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rqdag(3, 3, 0.);

    //Computing the CKM 
    CKM = (Vu.hconjugate()) * Vd;
    gslpp::matrix<gslpp::complex> CKMUnphys = CKM;
    //Extract the 4 parameters from the raw CKM
    ExtractParametersFromCKM();
    //Build the CKM with the 4 parameters
    UpdateCKM();

    

    double a11 = remainder(CKMUnphys(0, 0).arg() - CKM(0, 0).arg(), 2.*M_PI);
    double a12 = remainder(CKMUnphys(0, 1).arg() - CKM(0, 1).arg(), 2.*M_PI);
    double a13 = remainder(CKMUnphys(0, 2).arg() - CKM(0, 2).arg(), 2.*M_PI);


    //    double a23 = (gslpp::complex(CKMUnphys(1, 0) / CKM(1, 0))).arg() - a11 + a13;
    //    double a33 = (gslpp::complex(CKMUnphys(2, 0) / CKM(2, 0))).arg() - a11 + a13;
    double a23 = remainder(CKMUnphys(1, 0).arg() - CKM(1, 0).arg(), 2.*M_PI) - a11 + a13;
    double a33 = remainder(CKMUnphys(2, 0).arg() - CKM(2, 0).arg(), 2.*M_PI) - a11 + a13;

    gslpp::matrix<gslpp::complex> phi1(3, 3, 0.);
    phi1.assign(0, 0, 1.);
    phi1.assign(1, 1, gslpp::complex(1., a23 - a13, true));
    phi1.assign(2, 2, gslpp::complex(1., a33 - a13, true));

    gslpp::matrix<gslpp::complex> phi2dag(3, 3, 0.);
    phi2dag.assign(0, 0, gslpp::complex(1., -a11, true));
    phi2dag.assign(1, 1, gslpp::complex(1., -a12, true));
    phi2dag.assign(2, 2, gslpp::complex(1., -a13, true));

    gslpp::matrix<gslpp::complex> phie(3, 3, 0.);
    phie.assign(0, 0, gslpp::complex(1., -(Re(0, 0)).arg(), true));
    phie.assign(1, 1, gslpp::complex(1., -(Re(1, 1)).arg(), true));
    phie.assign(2, 2, gslpp::complex(1., -(Re(2, 2)).arg(), true));

    yuR = gslpp::matrix<double>(Su);
    yuI.reset();
    yeR = gslpp::matrix<double>(Se);
    yeI.reset();
    ydR = gslpp::matrix<double>(Sd);
    ydI.reset();

    Re = Re*phie;
    Rl = Rl*phie;

    if (basis == "UP") {

        Ru = Uu*phi1;
        Rd = Ud*phi2dag;
        Rq = Vu*phi1;


        gslpp::matrix<gslpp::complex> yd = ydR * CKM.hconjugate();
        ydR = yd.real();
        ydI = yd.imag();

    } else if (basis == "DOWN") {

        Ru = Uu*phi1;
        Rd = Ud*phi2dag;
        Rq = Vd*phi2dag;

        gslpp::matrix<gslpp::complex> yu = yuR*CKM;
        yuR = yu.real();
        yuI = yu.imag();

    } else
        std::cout << "WARNING: wrong basis choice: " << basis << ": Yukawa couplings not updated!!!" << std::endl;

    //I save the hermitian conjugates of the rotation 
    //matrix to be more efficient. 
    Rudag = Ru.hconjugate();
    Rddag = Rd.hconjugate();
    Rqdag = Rq.hconjugate();

    Redag = Re.hconjugate();
    Rldag = Rl.hconjugate();


    using namespace std;
    

    int a, b, c, d, p, r, s, t, n;

    //Auxiliary objects to perform the rotation
    gslpp::complex z, w, x, y, h, f;
    gslpp::complex le, lu, ld, qe, qu1, qu8, qd1, qd8;



    //Coefficients in the new basis

    //Class 5
    double cuHRp[3 * 3] = {0.};
    double cuHIp[3 * 3] = {0.};
    double cdHRp[3 * 3] = {0.};
    double cdHIp[3 * 3] = {0.};
    double ceHRp[3 * 3] = {0.};
    double ceHIp[3 * 3] = {0.};

    //Class 6
    double ceWRp[3 * 3] = {0.};
    double ceWIp[3 * 3] = {0.};
    double ceBRp[3 * 3] = {0.};
    double ceBIp[3 * 3] = {0.};

    double cuWRp[3 * 3] = {0.};
    double cuWIp[3 * 3] = {0.};
    double cuBRp[3 * 3] = {0.};
    double cuBIp[3 * 3] = {0.};
    double cuGRp[3 * 3] = {0.};
    double cuGIp[3 * 3] = {0.};

    double cdWRp[3 * 3] = {0.};
    double cdWIp[3 * 3] = {0.};
    double cdBRp[3 * 3] = {0.};
    double cdBIp[3 * 3] = {0.};
    double cdGRp[3 * 3] = {0.};
    double cdGIp[3 * 3] = {0.};

    //Class 7
    double cHl1Rp[3 * 3] = {0.};
    double cHl1Ip[3 * 3] = {0.};
    double cHl3Rp[3 * 3] = {0.};
    double cHl3Ip[3 * 3] = {0.};
    double cHeRp[3 * 3] = {0.};
    double cHeIp[3 * 3] = {0.};
    double cHq1Rp[3 * 3] = {0.};
    double cHq1Ip[3 * 3] = {0.};
    double cHq3Rp[3 * 3] = {0.};
    double cHq3Ip[3 * 3] = {0.};
    double cHdRp[3 * 3] = {0.};
    double cHdIp[3 * 3] = {0.};
    double cHuRp[3 * 3] = {0.};
    double cHuIp[3 * 3] = {0.};
    double cHudRp[3 * 3] = {0.};
    double cHudIp[3 * 3] = {0.};


    //Class 8
    double cuuRp[81] = {0.};
    double cuuIp[81] = {0.};
    double cddRp[81] = {0.};
    double cddIp[81] = {0.};
    double ceeRp[81] = {0.};
    double ceeIp[81] = {0.};

    double cllRp[81] = {0.};
    double cllIp[81] = {0.};
    double cqq1Rp[81] = {0.};
    double cqq1Ip[81] = {0.};
    double cqq3Rp[81] = {0.};
    double cqq3Ip[81] = {0.};

    double clq1Rp[81] = {0.};
    double clq1Ip[81] = {0.};
    double clq3Rp[81] = {0.};
    double clq3Ip[81] = {0.};

    double ceuRp[81] = {0.};
    double ceuIp[81] = {0.};
    double cedRp[81] = {0.};
    double cedIp[81] = {0.};
    double cud1Rp[81] = {0.};
    double cud1Ip[81] = {0.};
    double cud8Rp[81] = {0.};
    double cud8Ip[81] = {0.};

    double cleRp[81] = {0.};
    double cleIp[81] = {0.};
    double cluRp[81] = {0.};
    double cluIp[81] = {0.};
    double cldRp[81] = {0.};
    double cldIp[81] = {0.};
    double cqeRp[81] = {0.};
    double cqeIp[81] = {0.};

    double cqu1Rp[81] = {0.};
    double cqu1Ip[81] = {0.};
    double cqu8Rp[81] = {0.};
    double cqu8Ip[81] = {0.};
    double cqd1Rp[81] = {0.};
    double cqd1Ip[81] = {0.};
    double cqd8Rp[81] = {0.};
    double cqd8Ip[81] = {0.};


    double cledqRp[81] = {0.};
    double cledqIp[81] = {0.};

    double clequ1Rp[81] = {0.};
    double clequ1Ip[81] = {0.};
    double clequ3Rp[81] = {0.};
    double clequ3Ip[81] = {0.};

    double cquqd1Rp[81] = {0.};
    double cquqd1Ip[81] = {0.};
    double cquqd8Rp[81] = {0.};
    double cquqd8Ip[81] = {0.};

    //Class 5
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            z = gslpp::complex(0., 0.);
            for (a = 0; a < 3; a++) {
                for (b = 0; b < 3; b++) {
                    z += Rqdag(i, a) *
                            gslpp::complex(
                            WC1(cuHR, a, b), WC1(cuHI, a, b)) *
                            Ru(b, j);
                }
            }
            WC1_set(cuHRp, i, j, z.real());
            WC1_set(cuHIp, i, j, z.imag());
        }
    }
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            z = gslpp::complex(0., 0.);
            for (a = 0; a < 3; a++) {
                for (b = 0; b < 3; b++) {
                    z += Rqdag(i, a) *
                            gslpp::complex(
                            WC1(cdHR, a, b), WC1(cdHI, a, b)) *
                            Rd(b, j);
                }
            }
            WC1_set(cdHRp, i, j, z.real());
            WC1_set(cdHIp, i, j, z.imag());
        }
    }
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            z = gslpp::complex(0., 0.);
            for (a = 0; a < 3; a++) {
                for (b = 0; b < 3; b++) {
                    z += Rldag(i, a) *
                            gslpp::complex(
                            WC1(ceHR, a, b), WC1(ceHI, a, b)) *
                            Re(b, j);
                }
            }
            WC1_set(ceHRp, i, j, z.real());
            WC1_set(ceHIp, i, j, z.imag());
        }
    }



    //The coefficients are updated with the rotated ones
    std::copy(std::begin(cuHRp), std::end(cuHRp), std::begin(cuHR));
    std::copy(std::begin(cuHIp), std::end(cuHIp), std::begin(cuHI));
    std::copy(std::begin(cdHRp), std::end(cdHRp), std::begin(cdHR));
    std::copy(std::begin(cdHIp), std::end(cdHIp), std::begin(cdHI));
    std::copy(std::begin(ceHRp), std::end(ceHRp), std::begin(ceHR));
    std::copy(std::begin(ceHIp), std::end(ceHIp), std::begin(ceHI));


    //Class 6
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            z = gslpp::complex(0., 0.);
            for (a = 0; a < 3; a++) {
                for (b = 0; b < 3; b++) {
                    z += Rldag(i, a) *
                            gslpp::complex(
                            WC1(ceWR, a, b), WC1(ceWI, a, b)) *
                            Re(b, j);
                }
            }
            WC1_set(ceWRp, i, j, z.real());
            WC1_set(ceWIp, i, j, z.imag());
        }
    }
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            z = gslpp::complex(0., 0.);
            for (a = 0; a < 3; a++) {
                for (b = 0; b < 3; b++) {
                    z += Rldag(i, a) *
                            gslpp::complex(
                            WC1(ceBR, a, b), WC1(ceBI, a, b)) *
                            Re(b, j);
                }
            }
            WC1_set(ceBRp, i, j, z.real());
            WC1_set(ceBIp, i, j, z.imag());
        }
    }


    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            z = gslpp::complex(0., 0.);
            for (a = 0; a < 3; a++) {
                for (b = 0; b < 3; b++) {
                    z += Rqdag(i, a) *
                            gslpp::complex(
                            WC1(cuWR, a, b), WC1(cuWI, a, b)) *
                            Ru(b, j);
                }
            }
            WC1_set(cuWRp, i, j, z.real());
            WC1_set(cuWIp, i, j, z.imag());
        }
    }

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            z = gslpp::complex(0., 0.);
            for (a = 0; a < 3; a++) {
                for (b = 0; b < 3; b++) {
                    z += Rqdag(i, a) *
                            gslpp::complex(
                            WC1(cuBR, a, b), WC1(cuBI, a, b)) *
                            Ru(b, j);
                }
            }
            WC1_set(cuBRp, i, j, z.real());
            WC1_set(cuBIp, i, j, z.imag());
        }
    }

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            z = gslpp::complex(0., 0.);
            for (a = 0; a < 3; a++) {
                for (b = 0; b < 3; b++) {
                    z += Rqdag(i, a) *
                            gslpp::complex(
                            WC1(cuGR, a, b), WC1(cuGI, a, b)) *
                            Ru(b, j);
                }
            }
            WC1_set(cuGRp, i, j, z.real());
            WC1_set(cuGIp, i, j, z.imag());
        }
    }

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            z = gslpp::complex(0., 0.);
            for (a = 0; a < 3; a++) {
                for (b = 0; b < 3; b++) {
                    z += Rqdag(i, a) *
                            gslpp::complex(
                            WC1(cdWR, a, b), WC1(cdWI, a, b)) *
                            Rd(b, j);
                }
            }
            WC1_set(cdWRp, i, j, z.real());
            WC1_set(cdWIp, i, j, z.imag());
        }
    }

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            z = gslpp::complex(0., 0.);
            for (a = 0; a < 3; a++) {
                for (b = 0; b < 3; b++) {
                    z += Rqdag(i, a) *
                            gslpp::complex(
                            WC1(cdBR, a, b), WC1(cdBI, a, b)) *
                            Rd(b, j);
                }
            }
            WC1_set(cdBRp, i, j, z.real());
            WC1_set(cdBIp, i, j, z.imag());
        }
    }

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            z = gslpp::complex(0., 0.);
            for (a = 0; a < 3; a++) {
                for (b = 0; b < 3; b++) {
                    z += Rqdag(i, a) *
                            gslpp::complex(
                            WC1(cdGR, a, b), WC1(cdGI, a, b)) *
                            Rd(b, j);
                }
            }
            WC1_set(cdGRp, i, j, z.real());
            WC1_set(cdGIp, i, j, z.imag());
        }
    }


    std::copy(std::begin(ceWRp), std::end(ceWRp), std::begin(ceWR));
    std::copy(std::begin(ceWIp), std::end(ceWIp), std::begin(ceWI));
    std::copy(std::begin(ceBRp), std::end(ceBRp), std::begin(ceBR));
    std::copy(std::begin(ceBIp), std::end(ceBIp), std::begin(ceBI));

    std::copy(std::begin(cuWRp), std::end(cuWRp), std::begin(cuWR));
    std::copy(std::begin(cuWIp), std::end(cuWIp), std::begin(cuWI));
    std::copy(std::begin(cuBRp), std::end(cuBRp), std::begin(cuBR));
    std::copy(std::begin(cuBIp), std::end(cuBIp), std::begin(cuBI));
    std::copy(std::begin(cuGRp), std::end(cuGRp), std::begin(cuGR));
    std::copy(std::begin(cuGIp), std::end(cuGIp), std::begin(cuGI));

    std::copy(std::begin(cdWRp), std::end(cdWRp), std::begin(cdWR));
    std::copy(std::begin(cdWIp), std::end(cdWIp), std::begin(cdWI));
    std::copy(std::begin(cdBRp), std::end(cdBRp), std::begin(cdBR));
    std::copy(std::begin(cdBIp), std::end(cdBIp), std::begin(cdBI));
    std::copy(std::begin(cdGRp), std::end(cdGRp), std::begin(cdGR));
    std::copy(std::begin(cdGIp), std::end(cdGIp), std::begin(cdGI));




    //Class 7
    for (n = 0; n < DWC2R; n++) {
        i = WC2R_indices[n][0];
        j = WC2R_indices[n][1];

        z = gslpp::complex(0., 0.);
        for (a = 0; a < 3; a++) {
            for (b = 0; b < 3; b++) {
                z += Rldag(i, a) *
                        gslpp::complex(
                        WC2R(cHl1R, a, b), WC2I(cHl1I, a, b)) *
                        Rl(b, j);
            }
        }
        WC2R_set(cHl1Rp, i, j, z.real());
        WC2I_set(cHl1Ip, i, j, z.imag());
    }

    for (n = 0; n < DWC2R; n++) {
        i = WC2R_indices[n][0];
        j = WC2R_indices[n][1];

        z = gslpp::complex(0., 0.);
        for (a = 0; a < 3; a++) {
            for (b = 0; b < 3; b++) {
                z += Rldag(i, a) *
                        gslpp::complex(
                        WC2R(cHl3R, a, b), WC2I(cHl3I, a, b)) *
                        Rl(b, j);
            }
        }
        WC2R_set(cHl3Rp, i, j, z.real());
        WC2I_set(cHl3Ip, i, j, z.imag());
    }

    for (n = 0; n < DWC2R; n++) {
        i = WC2R_indices[n][0];
        j = WC2R_indices[n][1];

        z = gslpp::complex(0., 0.);
        for (a = 0; a < 3; a++) {
            for (b = 0; b < 3; b++) {
                z += Redag(i, a) *
                        gslpp::complex(
                        WC2R(cHeR, a, b), WC2I(cHeI, a, b)) *
                        Re(b, j);
            }
        }
        WC2R_set(cHeRp, i, j, z.real());
        WC2I_set(cHeIp, i, j, z.imag());
    }

    for (n = 0; n < DWC2R; n++) {
        i = WC2R_indices[n][0];
        j = WC2R_indices[n][1];

        z = gslpp::complex(0., 0.);
        for (a = 0; a < 3; a++) {
            for (b = 0; b < 3; b++) {
                z += Rqdag(i, a) *
                        gslpp::complex(
                        WC2R(cHq1R, a, b), WC2I(cHq1I, a, b)) *
                        Rq(b, j);
            }
        }
        WC2R_set(cHq1Rp, i, j, z.real());
        WC2I_set(cHq1Ip, i, j, z.imag());
    }

    for (n = 0; n < DWC2R; n++) {
        i = WC2R_indices[n][0];
        j = WC2R_indices[n][1];

        z = gslpp::complex(0., 0.);
        for (a = 0; a < 3; a++) {
            for (b = 0; b < 3; b++) {
                z += Rqdag(i, a) *
                        gslpp::complex(
                        WC2R(cHq3R, a, b), WC2I(cHq3I, a, b)) *
                        Rq(b, j);
            }
        }
        WC2R_set(cHq3Rp, i, j, z.real());
        WC2I_set(cHq3Ip, i, j, z.imag());
    }

    for (n = 0; n < DWC2R; n++) {
        i = WC2R_indices[n][0];
        j = WC2R_indices[n][1];

        z = gslpp::complex(0., 0.);
        for (a = 0; a < 3; a++) {
            for (b = 0; b < 3; b++) {
                z += Rddag(i, a) *
                        gslpp::complex(
                        WC2R(cHdR, a, b), WC2I(cHdI, a, b)) *
                        Rd(b, j);
            }
        }
        WC2R_set(cHdRp, i, j, z.real());
        WC2I_set(cHdIp, i, j, z.imag());
    }


    for (n = 0; n < DWC2R; n++) {
        i = WC2R_indices[n][0];
        j = WC2R_indices[n][1];

        z = gslpp::complex(0., 0.);
        for (a = 0; a < 3; a++) {
            for (b = 0; b < 3; b++) {
                z += Rudag(i, a) *
                        gslpp::complex(
                        WC2R(cHuR, a, b), WC2I(cHuI, a, b)) *
                        Ru(b, j);
            }
        }
        WC2R_set(cHuRp, i, j, z.real());
        WC2I_set(cHuIp, i, j, z.imag());
    }

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            z = gslpp::complex(0., 0.);
            for (a = 0; a < 3; a++) {
                for (b = 0; b < 3; b++) {
                    z += Rudag(i, a) *
                            gslpp::complex(
                            WC1(cHudR, a, b), WC1(cHudI, a, b)) *
                            Rd(b, j);
                }
            }
            WC1_set(cHudRp, i, j, z.real());
            WC1_set(cHudIp, i, j, z.imag());
        }
    }

    std::copy(std::begin(cHl1Rp), std::end(cHl1Rp), std::begin(cHl1R));
    std::copy(std::begin(cHl1Ip), std::end(cHl1Ip), std::begin(cHl1I));
    std::copy(std::begin(cHl3Rp), std::end(cHl3Rp), std::begin(cHl3R));
    std::copy(std::begin(cHl3Ip), std::end(cHl3Ip), std::begin(cHl3I));
    std::copy(std::begin(cHeRp), std::end(cHeRp), std::begin(cHeR));
    std::copy(std::begin(cHeIp), std::end(cHeIp), std::begin(cHeI));
    std::copy(std::begin(cHq1Rp), std::end(cHq1Rp), std::begin(cHq1R));
    std::copy(std::begin(cHq1Ip), std::end(cHq1Ip), std::begin(cHq1I));
    std::copy(std::begin(cHq3Rp), std::end(cHq3Rp), std::begin(cHq3R));
    std::copy(std::begin(cHq3Ip), std::end(cHq3Ip), std::begin(cHq3I));
    std::copy(std::begin(cHdRp), std::end(cHdRp), std::begin(cHdR));
    std::copy(std::begin(cHdIp), std::end(cHdIp), std::begin(cHdI));
    std::copy(std::begin(cHuRp), std::end(cHuRp), std::begin(cHuR));
    std::copy(std::begin(cHuIp), std::end(cHuIp), std::begin(cHuI));
    std::copy(std::begin(cHudRp), std::end(cHudRp), std::begin(cHudR));
    std::copy(std::begin(cHudIp), std::end(cHudIp), std::begin(cHudI));







    //Class 8 (WC6)
    for (n = 0; n < DWC6R; n++) {
        p = WC6R_indices[n][0];
        r = WC6R_indices[n][1];
        s = WC6R_indices[n][2];
        t = WC6R_indices[n][3];
        z = gslpp::complex(0., 0.);
        w = gslpp::complex(0., 0.);
        x = gslpp::complex(0., 0.);
        y = gslpp::complex(0., 0.);
        h = gslpp::complex(0., 0.);
        for (a = 0; a < 3; a++) {
            for (b = 0; b < 3; b++) {
                for (c = 0; c < 3; c++) {
                    for (d = 0; d < 3; d++) {
                        z += Rudag(p, a) * Rudag(s, c) *
                                gslpp::complex(
                                WC6R(cuuR, a, b, c, d), WC6I(cuuI, a, b, c, d))
                                * Ru(b, r) * Ru(d, t);
                        w += Rddag(p, a) * Rddag(s, c) *
                                gslpp::complex(
                                WC6R(cddR, a, b, c, d), WC6I(cddI, a, b, c, d))
                                * Rd(b, r) * Rd(d, t);
                        x += Rldag(p, a) * Rldag(s, c) *
                                gslpp::complex(
                                WC6R(cllR, a, b, c, d), WC6I(cllI, a, b, c, d))
                                * Rl(b, r) * Rl(d, t);
                        y += Rqdag(p, a) * Rqdag(s, c) *
                                gslpp::complex(
                                WC6R(cqq1R, a, b, c, d), WC6I(cqq1I, a, b, c, d))
                                * Rq(b, r) * Rq(d, t);
                        h += Rqdag(p, a) * Rqdag(s, c) *
                                gslpp::complex(
                                WC6R(cqq3R, a, b, c, d), WC6I(cqq3I, a, b, c, d))
                                * Rq(b, r) * Rq(d, t);

                    }
                }
            }
        }
        WC6R_set(cuuRp, p, r, s, t, z.real());
        WC6I_set(cuuIp, p, r, s, t, z.imag());
        WC6R_set(cddRp, p, r, s, t, w.real());
        WC6I_set(cddIp, p, r, s, t, w.imag());

        WC6R_set(cllRp, p, r, s, t, x.real());
        WC6I_set(cllIp, p, r, s, t, x.imag());
        WC6R_set(cqq1Rp, p, r, s, t, y.real());
        WC6I_set(cqq1Ip, p, r, s, t, y.imag());
        WC6R_set(cqq3Rp, p, r, s, t, h.real());
        WC6I_set(cqq3Ip, p, r, s, t, h.imag());
    }


    std::copy(std::begin(cuuRp), std::end(cuuRp), std::begin(cuuR));
    std::copy(std::begin(cuuIp), std::end(cuuIp), std::begin(cuuI));
    std::copy(std::begin(cddRp), std::end(cddRp), std::begin(cddR));
    std::copy(std::begin(cddIp), std::end(cddIp), std::begin(cddI));

    std::copy(std::begin(cllRp), std::end(cllRp), std::begin(cllR));
    std::copy(std::begin(cllIp), std::end(cllIp), std::begin(cllI));
    std::copy(std::begin(cqq1Rp), std::end(cqq1Rp), std::begin(cqq1R));
    std::copy(std::begin(cqq1Ip), std::end(cqq1Ip), std::begin(cqq1I));
    std::copy(std::begin(cqq3Rp), std::end(cqq3Rp), std::begin(cqq3R));
    std::copy(std::begin(cqq3Ip), std::end(cqq3Ip), std::begin(cqq3I));

    //Class 8 (WC7)
    for (n = 0; n < DWC7R; n++) {
        p = WC7R_indices[n][0];
        r = WC7R_indices[n][1];
        s = WC7R_indices[n][2];
        t = WC7R_indices[n][3];

        z = gslpp::complex(0., 0.);
        w = gslpp::complex(0., 0.);
        x = gslpp::complex(0., 0.);
        y = gslpp::complex(0., 0.);
        h = gslpp::complex(0., 0.);
        f = gslpp::complex(0., 0.);

        le = gslpp::complex(0., 0.);
        lu = gslpp::complex(0., 0.);
        ld = gslpp::complex(0., 0.);
        qe = gslpp::complex(0., 0.);

        qu1 = gslpp::complex(0., 0.);
        qu8 = gslpp::complex(0., 0.);

        qd1 = gslpp::complex(0., 0.);
        qd8 = gslpp::complex(0., 0.);

        for (a = 0; a < 3; a++) {
            for (b = 0; b < 3; b++) {
                for (c = 0; c < 3; c++) {
                    for (d = 0; d < 3; d++) {
                        z += Rldag(p, a) * Rqdag(s, c) *
                                gslpp::complex(
                                WC7R(clq1R, a, b, c, d), WC7I(clq1I, a, b, c, d))
                                * Rl(b, r) * Rq(d, t);
                        w += Rldag(p, a) * Rqdag(s, c) *
                                gslpp::complex(
                                WC7R(clq3R, a, b, c, d), WC7I(clq3I, a, b, c, d))
                                * Rl(b, r) * Rq(d, t);

                        x += Redag(p, a) * Rudag(s, c) *
                                gslpp::complex(
                                WC7R(ceuR, a, b, c, d), WC7I(ceuI, a, b, c, d))
                                * Re(b, r) * Ru(d, t);
                        y += Redag(p, a) * Rddag(s, c) *
                                gslpp::complex(
                                WC7R(cedR, a, b, c, d), WC7I(cedI, a, b, c, d))
                                * Re(b, r) * Rd(d, t);

                        h += Rudag(p, a) * Rddag(s, c) *
                                gslpp::complex(
                                WC7R(cud1R, a, b, c, d), WC7I(cud1I, a, b, c, d))
                                * Ru(b, r) * Rd(d, t);
                        f += Rudag(p, a) * Rddag(s, c) *
                                gslpp::complex(
                                WC7R(cud8R, a, b, c, d), WC7I(cud8I, a, b, c, d))
                                * Ru(b, r) * Rd(d, t);

                        le += Rldag(p, a) * Redag(s, c) *
                                gslpp::complex(
                                WC7R(cleR, a, b, c, d), WC7I(cleI, a, b, c, d))
                                * Rl(b, r) * Re(d, t);
                        lu += Rldag(p, a) * Rudag(s, c) *
                                gslpp::complex(
                                WC7R(cluR, a, b, c, d), WC7I(cluI, a, b, c, d))
                                * Rl(b, r) * Ru(d, t);
                        ld += Rldag(p, a) * Rddag(s, c) *
                                gslpp::complex(
                                WC7R(cldR, a, b, c, d), WC7I(cldI, a, b, c, d))
                                * Rl(b, r) * Rd(d, t);
                        qe += Rqdag(p, a) * Redag(s, c) *
                                gslpp::complex(
                                WC7R(cqeR, a, b, c, d), WC7I(cqeI, a, b, c, d))
                                * Rq(b, r) * Re(d, t);

                        qu1 += Rqdag(p, a) * Rudag(s, c) *
                                gslpp::complex(
                                WC7R(cqu1R, a, b, c, d), WC7I(cqu1I, a, b, c, d))
                                * Rq(b, r) * Ru(d, t);
                        qu8 += Rqdag(p, a) * Rudag(s, c) *
                                gslpp::complex(
                                WC7R(cqu8R, a, b, c, d), WC7I(cqu8I, a, b, c, d))
                                * Rq(b, r) * Ru(d, t);

                        qd1 += Rqdag(p, a) * Rddag(s, c) *
                                gslpp::complex(
                                WC7R(cqd1R, a, b, c, d), WC7I(cqd1I, a, b, c, d))
                                * Rq(b, r) * Rd(d, t);
                        qd8 += Rqdag(p, a) * Rddag(s, c) *
                                gslpp::complex(
                                WC7R(cqd8R, a, b, c, d), WC7I(cqd8I, a, b, c, d))
                                * Rq(b, r) * Rd(d, t);

                    }
                }
            }
        }
        WC7R_set(clq1Rp, p, r, s, t, z.real());
        WC7I_set(clq1Ip, p, r, s, t, z.imag());
        WC7R_set(clq3Rp, p, r, s, t, w.real());
        WC7I_set(clq3Ip, p, r, s, t, w.imag());

        WC7R_set(ceuRp, p, r, s, t, x.real());
        WC7I_set(ceuIp, p, r, s, t, x.imag());
        WC7R_set(cedRp, p, r, s, t, y.real());
        WC7I_set(cedIp, p, r, s, t, y.imag());

        WC7R_set(cud1Rp, p, r, s, t, h.real());
        WC7I_set(cud1Ip, p, r, s, t, h.imag());
        WC7R_set(cud8Rp, p, r, s, t, f.real());
        WC7I_set(cud8Ip, p, r, s, t, f.imag());

        WC7R_set(cleRp, p, r, s, t, le.real());
        WC7I_set(cleIp, p, r, s, t, le.imag());
        WC7R_set(cluRp, p, r, s, t, lu.real());
        WC7I_set(cluIp, p, r, s, t, lu.imag());
        WC7R_set(cldRp, p, r, s, t, ld.real());
        WC7I_set(cldIp, p, r, s, t, ld.imag());
        WC7R_set(cqeRp, p, r, s, t, qe.real());
        WC7I_set(cqeIp, p, r, s, t, qe.imag());

        WC7R_set(cqu1Rp, p, r, s, t, qu1.real());
        WC7I_set(cqu1Ip, p, r, s, t, qu1.imag());
        WC7R_set(cqu8Rp, p, r, s, t, qu8.real());
        WC7I_set(cqu8Ip, p, r, s, t, qu8.imag());
        WC7R_set(cqd1Rp, p, r, s, t, qd1.real());
        WC7I_set(cqd1Ip, p, r, s, t, qd1.imag());
        WC7R_set(cqd8Rp, p, r, s, t, qd8.real());
        WC7I_set(cqd8Ip, p, r, s, t, qd8.imag());

    }


    std::copy(std::begin(clq1Rp), std::end(clq1Rp), std::begin(clq1R));
    std::copy(std::begin(clq1Ip), std::end(clq1Ip), std::begin(clq1I));
    std::copy(std::begin(clq3Rp), std::end(clq3Rp), std::begin(clq3R));
    std::copy(std::begin(clq3Ip), std::end(clq3Ip), std::begin(clq3I));
    std::copy(std::begin(ceuRp), std::end(ceuRp), std::begin(ceuR));
    std::copy(std::begin(ceuIp), std::end(ceuIp), std::begin(ceuI));
    std::copy(std::begin(cedRp), std::end(cedRp), std::begin(cedR));
    std::copy(std::begin(cedIp), std::end(cedIp), std::begin(cedI));
    std::copy(std::begin(cud1Rp), std::end(cud1Rp), std::begin(cud1R));
    std::copy(std::begin(cud1Ip), std::end(cud1Ip), std::begin(cud1I));
    std::copy(std::begin(cud8Rp), std::end(cud8Rp), std::begin(cud8R));
    std::copy(std::begin(cud8Ip), std::end(cud8Ip), std::begin(cud8I));

    std::copy(std::begin(cleRp), std::end(cleRp), std::begin(cleR));
    std::copy(std::begin(cleIp), std::end(cleIp), std::begin(cleI));
    std::copy(std::begin(cldRp), std::end(cldRp), std::begin(cldR));
    std::copy(std::begin(cldIp), std::end(cldIp), std::begin(cldI));
    std::copy(std::begin(cluRp), std::end(cluRp), std::begin(cluR));
    std::copy(std::begin(cluIp), std::end(cluIp), std::begin(cluI));
    std::copy(std::begin(cqeRp), std::end(cqeRp), std::begin(cqeR));
    std::copy(std::begin(cqeIp), std::end(cqeIp), std::begin(cqeI));

    std::copy(std::begin(cqu1Rp), std::end(cqu1Rp), std::begin(cqu1R));
    std::copy(std::begin(cqu1Ip), std::end(cqu1Ip), std::begin(cqu1I));
    std::copy(std::begin(cqu8Rp), std::end(cqu8Rp), std::begin(cqu8R));
    std::copy(std::begin(cqu8Ip), std::end(cqu8Ip), std::begin(cqu8I));
    std::copy(std::begin(cqd1Rp), std::end(cqd1Rp), std::begin(cqd1R));
    std::copy(std::begin(cqd1Ip), std::end(cqd1Ip), std::begin(cqd1I));
    std::copy(std::begin(cqd8Rp), std::end(cqd8Rp), std::begin(cqd8R));
    std::copy(std::begin(cqd8Ip), std::end(cqd8Ip), std::begin(cqd8I));

    //Class 8 (WC8)
    for (n = 0; n < DWC8R; n++) {
        p = WC8R_indices[n][0];
        r = WC8R_indices[n][1];
        s = WC8R_indices[n][2];
        t = WC8R_indices[n][3];

        z = gslpp::complex(0., 0.);
        for (a = 0; a < 3; a++) {
            for (b = 0; b < 3; b++) {
                for (c = 0; c < 3; c++) {
                    for (d = 0; d < 3; d++) {
                        z += Redag(p, a) * Redag(s, c) *
                                gslpp::complex(
                                WC8R(ceeR, a, b, c, d), WC8I(ceeI, a, b, c, d))
                                * Re(b, r) * Re(d, t);
                    }
                }
            }
        }
        WC8R_set(ceeRp, p, r, s, t, z.real());
        WC8I_set(ceeIp, p, r, s, t, z.imag());
    }

    std::copy(std::begin(ceeRp), std::end(ceeRp), std::begin(ceeR));
    std::copy(std::begin(ceeIp), std::end(ceeIp), std::begin(ceeI));








    //Class 8 (WC5)
    for (p = 0; p < 3; p++) {
        for (r = 0; r < 3; r++) {
            for (s = 0; s < 3; s++) {
                for (t = 0; t < 3; t++) {
                    z = gslpp::complex(0., 0.);
                    x = gslpp::complex(0., 0.);
                    y = gslpp::complex(0., 0.);
                    h = gslpp::complex(0., 0.);
                    w = gslpp::complex(0., 0.);


                    for (a = 0; a < 3; a++) {
                        for (b = 0; b < 3; b++) {
                            for (c = 0; c < 3; c++) {
                                for (d = 0; d < 3; d++) {
                                    h += Rldag(p, a) * Rddag(s, c) *
                                            gslpp::complex(
                                            WC5(cledqR, a, b, c, d), WC5(cledqI, a, b, c, d))
                                            * Re(b, r) * Rq(d, t);

                                    z += Rldag(p, a) * Rqdag(s, c) *
                                            gslpp::complex(
                                            WC5(clequ1R, a, b, c, d), WC5(clequ1I, a, b, c, d))
                                            * Re(b, r) * Ru(d, t);
                                    w += Rldag(p, a) * Rqdag(s, c) *
                                            gslpp::complex(
                                            WC5(clequ3R, a, b, c, d), WC5(clequ3I, a, b, c, d))
                                            * Re(b, r) * Ru(d, t);

                                    x += Rqdag(p, a) * Rqdag(s, c) *
                                            gslpp::complex(
                                            WC5(cquqd1R, a, b, c, d), WC5(cquqd1I, a, b, c, d))
                                            * Ru(b, r) * Rd(d, t);
                                    y += Rqdag(p, a) * Rqdag(s, c) *
                                            gslpp::complex(
                                            WC5(cquqd8R, a, b, c, d), WC5(cquqd8I, a, b, c, d))
                                            * Ru(b, r) * Rd(d, t);



                                }
                            }
                        }
                    }

                    WC5_set(cledqRp, p, r, s, t, h.real());
                    WC5_set(cledqIp, p, r, s, t, h.imag());

                    WC5_set(clequ1Rp, p, r, s, t, z.real());
                    WC5_set(clequ1Ip, p, r, s, t, z.imag());
                    WC5_set(clequ3Rp, p, r, s, t, w.real());
                    WC5_set(clequ3Ip, p, r, s, t, w.imag());

                    WC5_set(cquqd1Rp, p, r, s, t, x.real());
                    WC5_set(cquqd1Ip, p, r, s, t, x.imag());
                    WC5_set(cquqd8Rp, p, r, s, t, y.real());
                    WC5_set(cquqd8Ip, p, r, s, t, y.imag());


                }
            }
        }
    }


    std::copy(std::begin(cledqRp), std::end(cledqRp), std::begin(cledqR));
    std::copy(std::begin(cledqIp), std::end(cledqIp), std::begin(cledqI));

    std::copy(std::begin(clequ1Rp), std::end(clequ1Rp), std::begin(clequ1R));
    std::copy(std::begin(clequ1Ip), std::end(clequ1Ip), std::begin(clequ1I));
    std::copy(std::begin(clequ3Rp), std::end(clequ3Rp), std::begin(clequ3R));
    std::copy(std::begin(clequ3Ip), std::end(clequ3Ip), std::begin(clequ3I));

    std::copy(std::begin(cquqd1Rp), std::end(cquqd1Rp), std::begin(cquqd1R));
    std::copy(std::begin(cquqd1Ip), std::end(cquqd1Ip), std::begin(cquqd1I));
    std::copy(std::begin(cquqd8Rp), std::end(cquqd8Rp), std::begin(cquqd8R));
    std::copy(std::begin(cquqd8Ip), std::end(cquqd8Ip), std::begin(cquqd8I));



}

void RGESolver::EvolveToBasis(
        std::string method, double muI, double muF,
        std::string basis)
{
    Evolve(method, muI, muF);
    GoToBasis(basis);

}



