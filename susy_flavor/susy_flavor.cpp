#include<iostream>
/*This is a special type of statement called a preprocessor directive. 
Preprocessor directives tell the compiler to perform a special task. 
In this case, we are telling the compiler that we would like to use 
the iostream library. The iostream library contains code that tells 
the compiler what cout and endl do. In other words, we need to include 
the iostream library in order to write to the screen.*/
#include <fstream>
#include <stdlib.h>
#include<iomanip>
#include<stdio.h>
#include<complex> /* Standard Library for Complex Numbers */
using namespace std;
/*As we learned before, cout and endl live inside the iostream library. 
However, within iostream, they live inside a special compartment named 
std (short for standard). This using statement tells the compiler to 
look inside a compartment named std if it can’t find cout or endl 
defined anywhere else. In other words, this statement is also necessary 
so that the compiler can find cout and endl.*/

// Forward Declaration ============================================              

extern "C"
{ 
void vpar_update_(double&, double&, double&);
void lam_fit_(double&);
void lam_fit_nlo_(double&);
void ckm_wolf_(double&, double&, double&, double&);
void init_fermion_sector_(double&, double&, double&, double&);
int  init_higgs_sector_(double&, double&, complex<double>&, int&);
int  init_ino_sector_(complex<double>&, complex<double>&, double&,
                      complex<double>&, double&, int&);
int  init_slepton_sector_(double sll[3],double slr[3],complex<double> slmi_l[3],
                          complex<double> slmi_r[3],complex<double> slmi_lr[3][3],
                          complex<double> slmi_lrp[3][3], int&);

int init_squark_sector_(double amsq[3],double amsu[3],double amsd[3],
    complex<double> sqmi_l[3],complex<double> sumi_r[3],complex<double> sdmi_r[3],
    complex<double> sumi_lr[3][3],complex<double> sdmi_lr[3][3],
    complex<double> sumi_lrp[3][3],complex<double> sdmi_lrp[3][3],int&);
void reset_phys_data_();
int fcorr_epa_(double &,double &,double &,double &,double &,
               double &,double &,double &,double &,int &);
void print_mssm_par_();
void print_mssm_masses_();
int  set_resummation_level_(int &,int &);
double chiral_corr_size_(double corr_l[3],double corr_d[3], double corr_u[3],
double corr_ckm[3][3]);
double edm_l_(const int &);
double edm_n_();
double g_minus_2_anomaly_(const int &);
double br_llg_(const int &,const int &);
double k_pivv_(double &,double &);
double b_ll_(const int &, const int &, const int &, const int &);
double b_taunu_(double &, double &);
double bxg_nl_(double &,double &);
double dd_kaon_(double &,double &);
double uu_dmeson_(double &);
double dd_bmeson_(const int &, double &, double &, double &);
}  

// Access Fortran Common Blocks ===================================

typedef struct { double eps; int indx[3][3]; int iconv; } sf_cont_cpp;
/* typedef and struct are c++ keywords */



// Initialization

double var(string);
double v1();
double v2();

// The main program ================================================
int main() 
{ 

int  ilev, ierr;
double zm0,wm0,alpha_z,
       alpha_s,
       alam,apar,rhobar,etabar,
       top,top_scale,bot,bot_scale,
       pm,tanbe,
       amglu;

double abs_amue, abs_sdmi_lr_33, abs_sumi_lr_33;

complex<double> amg,amgg,amue;

double sll[3],slr[3],amsq[3],amsu[3],amsd[3],
       corr_l[3],corr_d[3],corr_u[3],corr_ckm[3][3];

complex<double> slmi_l[3], slmi_r[3], slmi_lr[3][3], slmi_lrp[3][3],
                sqmi_l[3], sdmi_r[3], sumi_r[3], sdmi_lr[3][3], 
                sumi_lr[3][3], sdmi_lrp[3][3], sumi_lrp[3][3];

double br_k0,br_kp,br_taunu,dtaunu_ratio,delb,amiu_b,
eps_k,delta_mk,delta_md,delta_mbd,dmb_re,dmb_im,delta_mbs;

// Access Fortran Common Blocks ===================================
extern sf_cont_cpp sf_cont_;

// Input parameters convention choice ===============================
   sf_cont_.iconv = 1;    // SLHA2 input conventions
// sf_cont_.iconv = 2;    // hep-ph/9511250 input conventions

// fixes the treatment of enhanced chiral correction resummation ====
// ilev = 0;  // No resummation, SUSY corrections strictly 1-loop
// ilev = 1;  // Resummation using the decoupling limit
   ilev = 2;  // Exact iterative solution, may not always converge

// SM basic input initialization ====================================
   zm0 = var("M_z");          // M_Z
   wm0 = 80.398e0;           // M_W
   alpha_z = 1/127.934e0;    // alpha_em(M_Z)

   vpar_update_(zm0,wm0,alpha_z);

// QCD parameters ===================================================
   alpha_s = 0.1172e0;      // alpha_s(MZ)
    
   lam_fit_(alpha_s);     // fits Lambda_QCD at 3 loop level
   lam_fit_nlo_(alpha_s); // fits Lambda_QCD at NLO level

// CKM matrix initialization ========================================
   alam = var("lambda");       // lambda
   apar = var("A");        // A
   rhobar = var("rhob");      // rho bar
   etabar = var("etab");      // eta bar
    
   ckm_wolf_(alam,apar,rhobar,etabar);

// Fermion mass initialization, input: MSbar running quark masses
   top_scale = 163.2e0;
   top = 163.2e0;            // m_t(top_scale)
   bot_scale = 4.17e0;
   bot = 4.17e0;             // m_b(bot_scale)
    
   init_fermion_sector_(top,top_scale,bot,bot_scale);

// Higgs sector parameters
   pm    = 200;            // M_A
   tanbe = var("tanb");             // tan(beta)
   amue= complex<double>(var("muHr"),var("muHi")); // mu

   init_higgs_sector_(pm,tanbe,amue,ierr);

   if (ierr != 0) 
      {
        cout << "negative tree level Higgs mass^2?" << endl;
        return 0;
      } 

// Gaugino sector parameters. CAUTION: if M1 is set to 0 here then
// program sets M1 and M2 GUT-related, i.e. M1 = 5/3 s_W^2/c_W^2*M2
    amgg  = complex<double>(var("m1r"),var("m1i"));  // M1 (bino mass)
    amg   = complex<double>(var("m2r"),var("m2i"));  // M2 (wino mass)
    amglu = var("m3");                           // M3 (gluino mass)
    
    init_ino_sector_(amgg,amg,amglu,amue,tanbe,ierr);
    
    if (ierr != 0) 
       {
        cout << "-ino mass below M_Z/2?" << endl;
       } 

// Slepton diagonal soft breaking parameters
      sll[0] = 300.0e0;           // left selectron mass scale
      sll[1] = 300.0e0;           // left smuon mass scale
      sll[2] = 300.0e0;           // left stau mass scale
      slr[0] = 300.0e0;           // right selectron mass scale
      slr[1] = 300.0e0;           // right smuon mass scale
      slr[2] = 300.0e0;           // right stau mass scale

// Slepton LL and RR mass insertions (hermitian matrices)
// slmi_x(1),slmi_x(2),slmi_x(3) are 12,23,31 entries respectively

      for (int i=0; i<3; i++)
           {
             slmi_l[i] = complex<double>(0.0e0,0.0e0); // slepton LL mass insertion
             slmi_r[i] = complex<double>(0.0e0,0.0e0); // slepton RR mass insertion
           } 

       slmi_l[0] = complex<double>(2.0e-2,1.0e-2); // non-vanishing LL 12 entry

// Slepton LR mass insertions, non-hermitian in general
// All entries dimensionless (normalized to diagonal masses)
      for (int i=0; i<3; i++) 
          {     
          for (int j=0; j<3; j++)
              {
              // holomorphic LR mixing terms
              slmi_lr[i][j] =complex<double>(0.0e0,0.0e0);
              // non-holomorphic LR mixing terms
              slmi_lrp[i][j] =complex<double>(0.0e0,0.0e0);
              }
          } 

// Example: diagonal entries normalized to Y_l as in SUGRA  
      slmi_lr[0][0] = complex<double>(1.5e-4,0.0e0); // A_e
      slmi_lr[1][1] = complex<double>(3.0e-2,0.0e0); // A_mu
      slmi_lr[2][2] = complex<double>(5.0e-1,0.0e0); // A_tau
      slmi_lr[2][1] = complex<double>(1.0e-2,2.0e-2); // non-vanishing LR 23 entry


// Calculate physical masses and mixing angles

      init_slepton_sector_(sll,slr,slmi_l,slmi_r,slmi_lr,slmi_lrp,ierr);
  
   if (ierr != 0) 
      {
        cout << "Negative tree level slepton mass^2?" << endl;
        return 0;
      } 

// Squark diagonal soft breaking parameters

      amsq[0] = sqrt(var("msQ2_11r"));          // left squark mass, 1st generation       // means sqrt(msQ2_11r)
      amsq[1] = sqrt(var("msQ2_22r"));          // left squark mass, 2nd generation       // means sqrt(msQ2_22r)  
      amsq[2] = sqrt(var("msQ2_33r"));          // left squark mass, 3rd generation       // means sqrt(msQ2_33r)
      amsd[0] = sqrt(var("msD2_11r"));          // right down squark mass                 // means sqrt(msD2_11r)
      amsd[1] = sqrt(var("msD2_22r"));          // right strange squark mass
      amsd[2] = sqrt(var("msD2_33r"));          // right sbottom mass
      amsu[0] = sqrt(var("msU2_11r"));          // right up squark mass
      amsu[1] = sqrt(var("msU2_22r"));          // right charm squark mass
      amsu[2] = sqrt(var("msU2_33r"));          // right stop mass

// Squark LL and RR mass insertions (hermitian matrices)
// sqmi_l(1),sqmi_l(2), sqmi_l(3) are 12,23,31 entry respectively, etc.  // Maybe is wrong the entry for sqmi_l(3) is the 13 entry ???

    //  for (int i=0; i<3; i++) 
    //      {
    //      sqmi_l[i] = complex<double>(0.0e0,0.0e0); // squark LL mass insertion
    //      sumi_r[i] = complex<double>(0.0e0,0.0e0); // up-squark RR mass insertion
    //      sdmi_r[i] = complex<double>(0.0e0,0.0e0); // down-squark RR mass insertion
    //      } 
    // 
    //  sqmi_l[1] = complex<double>(2.0e-2,-1.0e-2); // non-vanishing LL 23 entry
    
    // test - matrix
    //
    //   sqmi_l[1] = complex<double>(2.0e-2,1.e-2);
    //
    // end - test 
    
    /// Initialization from .txt file
    
      sqmi_l[0] = complex<double>(var("msQ2_12r")/(amsq[0] * amsq[1]), var("msQ2_12i")/(amsq[0] *amsq[1]));     // sqmi_l[0] means delta_12 
      sqmi_l[1] = complex<double>(var("msQ2_23r")/(amsq[1] * amsq[2]), var("msQ2_23i")/(amsq[1] * amsq[2]));   // sqmi_l[1] means delta_23
      sqmi_l[2] = complex<double>(var("msQ2_13r")/(amsq[0] * amsq[2]), var("msQ2_13i")/(amsq[0] * amsq[2]));  // sqmi_l[2] means delta_31
     
      sumi_r[0] = complex<double>(var("msU2_12r")/(amsu[0] * amsu[1]), var("msU2_12i")/(amsu[0]*amsu[1]));
      sumi_r[1] = complex<double>(var("msU2_23r")/(amsu[1] * amsu[2]), var("msU2_23i")/(amsu[1]*amsu[2]));
      sumi_r[2] = complex<double>(var("msU2_13r")/(amsu[0] * amsu[2]), var("msU2_13i")/(amsu[0]*amsu[2]));
    
      sdmi_r[0] = complex<double>(var("msD2_12r")/(amsu[0] * amsu[1]), var("msD2_12i")/(amsu[0]*amsu[1]));
      sdmi_r[1] = complex<double>(var("msD2_23r")/(amsu[1] * amsu[2]), var("msD2_23i")/(amsu[1]*amsu[2]));
      sdmi_r[2] = complex<double>(var("msD2_13r")/(amsu[0] * amsu[2]), var("msD2_13i")/(amsu[0]*amsu[2]));
    
      

// Squark holomorphic LR mass insertions, non-hermitian in general
// All entries dimensionless (normalized to masses)

      for (int i=0; i<3; i++) 
          {     
          for (int j=0; i<3; i++)
              {
              // holomorphic LR mixing terms
              sumi_lr[i][j] = complex<double>(0.0e0,0.0e0); // up-squark 
              sdmi_lr[i][j] = complex<double>(0.0e0,0.0e0); // down-squark 
              // non-holomorphic LR mixing terms
              sumi_lrp[i][j] = complex<double>(0.0e0,0.0e0); // up-squark
              sdmi_lrp[i][j] = complex<double>(0.0e0,0.0e0); // down-squark
              }
          }

// Example: diagonal entries normalized to Y_d,Y_u as in SUGRA  

   //   sumi_lr[0][0] = complex<double>(1.0e-5,0.0e0);
   //   sumi_lr[1][1] = complex<double>(4.0e-3,0.0e0);
   //   sumi_lr[2][2] = complex<double>(1.0e0,0.0e0);
   //   sdmi_lr[0][0] = complex<double>(-1.0e-3,0.0e0);
   //   sdmi_lr[1][1] = complex<double>(-2.0e-2,0.0e0);
   //   sdmi_lr[2][2] = complex<double>(-8.0e-1,0.0e0);
   //   sumi_lr[2][1] = complex<double>(1.0e-2,2.0e-2);  // non-vanishing up LR 23 entry     // Maybe is the entry 32 ???
   //   sdmi_lr[2][1] = complex<double>(-3.0e-2,1.0e-2); // non-vanishing down LR 23 entry

    
    //
    // Initialization (the elements LR are a bit different because here the quark masses are non running)
    //
    // demand to Luca
    
    
    sumi_lr[0][0] = -v2()/sqrt(2) * complex<double>(var("TU_11r"),var("TU_11i"))/(sqrt(amsq[0] * amsu[0]));  
    sumi_lr[1][0] = -v2()/sqrt(2) * complex<double>(var("TU_12r"),var("TU_12i"))/(sqrt(amsq[0] * amsu[1]));
    sumi_lr[2][0] = -v2()/sqrt(2) * complex<double>(var("TU_13r"),var("TU_13i"))/(sqrt(amsq[0] * amsu[2]));
    sumi_lr[0][1] = -v2()/sqrt(2) * complex<double>(var("TU_21r"),var("TU_21i"))/(sqrt(amsq[1] * amsu[0]));
    sumi_lr[1][1] = -v2()/sqrt(2) * complex<double>(var("TU_22r"),var("TU_22i"))/(sqrt(amsq[1] * amsu[1]));
    sumi_lr[2][1] = -v2()/sqrt(2) * complex<double>(var("TU_23r"),var("TU_23i"))/(sqrt(amsq[1] * amsu[2]));
    sumi_lr[0][2] = -v2()/sqrt(2) * complex<double>(var("TU_31r"),var("TU_31i"))/(sqrt(amsq[2] * amsu[0]));
    sumi_lr[1][2] = -v2()/sqrt(2) * complex<double>(var("TU_32"),var("TU_32i"))/(sqrt(amsq[2] * amsu[1]));
    sumi_lr[2][2] = -(v2()/sqrt(2) * complex<double>(var("TU_33r"),var("TU_33i")) - amue * top / tanbe)/(sqrt(amsq[2] * amsu[2]));
    
    sdmi_lr[0][0] = v1()/sqrt(2) * complex<double>(var("TD_11r"),var("TD_11i"))/(sqrt(amsq[0] * amsd[0]));
    sdmi_lr[1][0] = v1()/sqrt(2) * complex<double>(var("TD_12r"),var("TD_12i"))/(sqrt(amsq[0] * amsd[1]));
    sdmi_lr[2][0] = v1()/sqrt(2) * complex<double>(var("TD_13r"),var("TD_13i"))/(sqrt(amsq[0] * amsd[2]));
    sdmi_lr[0][1] = v1()/sqrt(2) * complex<double>(var("TD_21r"),var("TD_21i"))/(sqrt(amsq[1] * amsd[0]));
    sdmi_lr[1][1] = v1()/sqrt(2) * complex<double>(var("TD_22r"),var("TD_22i"))/(sqrt(amsq[1] * amsd[1]));
    sdmi_lr[2][1] = v1()/sqrt(2) * complex<double>(var("TD_23r"),var("TD_23i"))/(sqrt(amsq[1] * amsd[2]));
    sdmi_lr[0][2] = v1()/sqrt(2) * complex<double>(var("TD_31r"),var("TD_31i"))/(sqrt(amsq[2] * amsd[0]));
    sdmi_lr[1][2] = v1()/sqrt(2) * complex<double>(var("TD_32r"),var("TD_32i"))/(sqrt(amsq[2] * amsd[1]));
    sdmi_lr[2][2] = (v1()/sqrt(2) * complex<double>(var("TD_33r"),var("TD_33i")) - amue * bot * tanbe)/(sqrt(amsq[2] * amsd[2]));
    
    
// Calculate physical masses and mixing angles

      init_squark_sector_(amsq,amsu,amsd,sqmi_l,sumi_r,sdmi_r,
                          sumi_lr,sdmi_lr,sumi_lrp,sdmi_lrp,ierr);
  
    if (ierr != 0) 
      {
        cout << "Negative tree level squark mass^2?" << endl;
        return 0;
      } 

//  Reset status of physical Higgs mass after parameter changes
    reset_phys_data_();

//  Neutral CP-even Higgs masses in the 1-loop Effective Potential
//  Approximation. Only real mu, A_t, A_b allowed - replaced x->abs(x)

    abs_amue=abs(amue);
    abs_sdmi_lr_33=abs(sdmi_lr[2][2]);
    abs_sumi_lr_33=abs(sumi_lr[2][2]);

    fcorr_epa_(tanbe,pm,top,abs_amue,amsq[2],amsd[2],amsu[2],
              abs_sdmi_lr_33,abs_sumi_lr_33,ierr);
  

    if (ierr != 0) 
      {
        cout << "Negative 1-loop EPA CP-even Higgs mass^2?" << endl;
        return 0;
      }

// End of input section =============================================

    
    
    /*************************************************************/
    // In what follows I drop all the prints 
    
    
// Control output
   //cout << setw(50) << left << "************************************************" << endl;
   //cout << setw(50) << left << "MSSM Lagrangian parameters and tree level masses" << endl;
   //cout << setw(50) << left << "written on file mssm_data.txt" << endl;
   //cout << setw(50) << left << "************************************************" << endl;

   print_mssm_par_();    // Lagrangian parameters
   print_mssm_masses_(); // tree level physical masses
   
// perform resummation of chirally enhanced corrections
   set_resummation_level_(ilev,ierr);

   chiral_corr_size_(corr_l,corr_d,corr_u,corr_ckm);


  //cout << " " << endl;
  //printf ("Corrections to lepton Yukawa couplings:      %8.1f%% %8.1f%% %8.1f%% \n", 
  //         100*corr_l[0],100*corr_l[1],100*corr_l[2]);
  //printf ("Corrections to down-quark Yukawa couplings:  %8.1f%% %8.1f%% %8.1f%% \n", 
  //         100*corr_d[0],100*corr_d[1],100*corr_d[2]);
  //printf ("Corrections to up-quark Yukawa couplings:    %8.1f%% %8.1f%% %8.1f%% \n", 
  //         100*corr_u[0],100*corr_u[1],100*corr_u[2]);
  //cout << "Corrections to CKM matrix elements: " << endl;
  //for (int i=0; i<3; i++)
  //    {
  //     printf (" %8.1f%% %8.1f%% %8.1f%% \n", 
  //         100*corr_ckm[0][i],100*corr_ckm[1][i],100*corr_ckm[2][i]);
  //    }
  //cout << " " << endl;

//  Results for implemented observables:
     //cout <<" ******************************************************"<< endl;
 
     //printf ("Physical results for resummation level %d (error code %d)\n", 
     //       ilev,ierr);

      //cout << "   " << endl;
      //cout << "Electric dipole moments:" << endl;
      //printf ("Electron EDM = %11.4e \n",edm_l_(1));
      //printf ("Muon EDM     = %11.4e \n",edm_l_(2));
      //printf ("Tau EDM      = %11.4e \n",edm_l_(3));
      //printf ("Neutron EDM  = %11.4e \n",edm_n_());
      //cout << " " << endl;

      //cout << "g-2 anomaly, SUSY contribution:" << endl;
      //printf ("Electron (g-2)/2 = %11.4e \n",g_minus_2_anomaly_(1));
      //printf ("Muon (g-2)/2     = %11.4e \n",g_minus_2_anomaly_(2));
      //printf ("Tau  (g-2)/2     = %11.4e \n",g_minus_2_anomaly_(3));
      //cout << " " << endl;

      //cout << "l^J->l^I gamma decays:" << endl;
      //printf ("Br(mu-> e gamma)   = %11.4e \n",br_llg_(2,1));
      //printf ("Br(tau-> e gamma)  = %11.4e \n",br_llg_(3,1));
      //printf ("Br(tau-> mu gamma) = %11.4e \n",br_llg_(3,2));
      //cout << " " << endl;

      //cout << "Neutrino K decays:" << endl;
      k_pivv_(br_k0,br_kp);
      //printf ("BR(K_L^0 -> pi^0 vv) = %11.4e \n",br_k0);
      //printf ("BR(K^+   -> pi^+ vv) = %11.4e \n",br_kp);
      //cout << " " << endl;


      //cout << "Leptonic B decays:" << endl;
      //printf ("BR(B_d -> mu^+ mu^-) = %11.4e \n",b_ll_(3,1,2,2));
      //printf ("BR(B_s -> mu^+ mu^-) = %11.4e \n",b_ll_(3,2,2,2));
      //printf ("BR(B_s -> mu^+ e^-)  = %11.4e \n",b_ll_(3,2,2,1));

      b_taunu_(br_taunu,dtaunu_ratio);

      //printf ("BR(B_u -> tau nu)    = %11.4e \n",br_taunu);
      //printf ("BR(B_u -> D tau nu)/BR(B_u -> D l nu) = %11.4e \n",dtaunu_ratio);
      //cout << " " << endl;

//    Physical quantities for BR(B->X_s g) calculation
      delb = 0.99; //Photon energy infrared cutoff
      amiu_b= 4.8; //Renormalization scale miu_b
      //printf ("BR(B -> X_S gamma) = %11.4e \n",bxg_nl_(delb,amiu_b));
      //cout << " " << endl;

      //cout << "KK mixing:" << endl;
      dd_kaon_(eps_k,delta_mk);

      printf ("eps_K     = %11.4e \n",eps_k);
      //printf ("Delta m_K = %11.4e \n",delta_mk);
      //cout << " " << endl;

      //cout << "DD mixing:" << endl;
      uu_dmeson_(delta_md);
      printf ("Delta m_D = %11.4e \n",delta_md);
      //cout << " " << endl;

      //cout << "BB mixing:" << endl;
      dd_bmeson_(1,delta_mbd,dmb_re,dmb_im);

      printf ("Re(H_eff_Bd) = %11.4e \n",dmb_re);
      printf ("Im(H_eff_Bd) = %11.4e \n",dmb_im);
      printf ("Delta m_B_d  = %11.4e \n",delta_mbd / 6.58211899e-13);
      
      dd_bmeson_(2,delta_mbs,dmb_re,dmb_im);
      printf ("Delta m_B_s  = %11.4e \n",delta_mbs);
      //cout << " " << endl;

    
    /*************************************************************/
    
    //std::ofstream myfile("out_sflav.txt");
    //if (myfile.is_open())
    //{ myfile << 1;
    //}else std::cout << "out_sflav.txt";
    
return 0;
}// End of main()


// Initialization

double var(const string variable){
    
    const string com1 = "grep -r ";
    const string com2 ="| awk '{print $(NF-2)}' > var.txt \n";
    const string file ="GeneralSUSY_NEW.conf";
    string com;
    com.append(com1);
    com.append(variable);
    com.append(" ");
    com.append(file);
    com.append(com2);
    const char * c = com.c_str();
    
    // Test - lines 
    //cout << c << endl;
    // end - test
    
    
    system(c);
    
    
    double value;
    FILE* pointer;
    pointer = fopen("var.txt" , "r" );
    fscanf(pointer , "%lf" , &value);
    
    
    return (value);
    
}


// sinb = tanb * sqrt(1. / (1. + tanb * tanb));     from the project SUSY.cpp
// cosb = sqrt(1. / (1. + tanb * tanb));            from the project SUSY.cpp

double v1(){
    
    double tanb = var("tanb");
    return (1/(sqrt(sqrt(2) * var("GF"))) *  sqrt(1. / (1. + tanb * tanb)));   
    
}

double v2(){
    
    double tanb = var("tanb");
    return (1/(sqrt(sqrt(2) * var("GF"))) * tanb * sqrt(1. / (1. + tanb * tanb)));
}



