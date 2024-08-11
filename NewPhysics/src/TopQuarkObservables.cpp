#include "TopQuarkObservables.h"
#include "std_make_vector.h"


TopQuarkObservables::TopQuarkObservables(const NPSMEFTd6General& NP_i) : myNPSMEFTd6General(NP_i) {}


//// ALL THE UNITS OF THE WILSON COEFFICIENTS ARE IN TeV^(-2). THIS FACTOR IS CONSIDERED IN THE GETTER 




////// ttbar Forward Backward asymmetry at Tevatron ///////////////////////////////////


FB_asymmetry_Tevatron_tt_diff_mtt_LO::FB_asymmetry_Tevatron_tt_diff_mtt_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{

    setParametersForObservable(make_vector<std::string>() 
            << "SM_FB_asymmetry_bin_mtt_350_450" << "SM_FB_asymmetry_deno_bin_mtt_350_450"
            << "SM_FB_asymmetry_bin_mtt_450_550" << "SM_FB_asymmetry_deno_bin_mtt_450_550"
            << "SM_FB_asymmetry_bin_mtt_550_650" << "SM_FB_asymmetry_deno_bin_mtt_550_650"
            << "SM_FB_asymmetry_bin_mtt_650_750" << "SM_FB_asymmetry_deno_bin_mtt_650_750");
            
    
}


double FB_asymmetry_Tevatron_tt_diff_mtt_LO::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    bool   flag_Quadratic= false; //Needs to be properly defined
        
    
    //double C_tG = ewgc("CuGR",2,2);

    
    
    
    if(b_min == 350 && b_max == 450){
        
        
        double SM_FB_asymmetry_bin_mtt_350_450 = SM.getOptionalParameter("SM_FB_asymmetry_bin_mtt_350_450");
        double SM_FB_asymmetry_deno_bin_mtt_350_450 = SM.getOptionalParameter("SM_FB_asymmetry_deno_bin_mtt_350_450");
        double SM_FB_asymmetry_num_bin_mtt_350_450 = SM_FB_asymmetry_bin_mtt_350_450*SM_FB_asymmetry_deno_bin_mtt_350_450;
        
        double SM_sigma_pos_bin_mtt_350_450 =0.5*(SM_FB_asymmetry_num_bin_mtt_350_450+SM_FB_asymmetry_deno_bin_mtt_350_450);
        double SM_sigma_neg_bin_mtt_350_450 =0.5*(SM_FB_asymmetry_deno_bin_mtt_350_450-SM_FB_asymmetry_num_bin_mtt_350_450);

                
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            
            
            double sigma_pos_bin_mtt_350_450_madgraph = 1.675650;
            double sigma_pos_bin_mtt_350_450_NP = -0.006110*ewgc("CG")+0.015540*ewgc("Cqd8R",2,2,0,0)
            +0.304934*ewgc("Cqq1R",0,2,2,0)+0.419488*ewgc("Cqq3R",0,2,2,0)+0.004257*ewgc("Cqq3R",1,2,2,1)
            +0.087584*ewgc("Cqu8R",0,0,2,2)+0.073317*ewgc("Cqu8R",2,2,0,0)+0.014866*ewgc("Cud8R",2,2,0,0)
            -0.385947*ewgc("CuGR",2,2)+0.298571*ewgc("CuuR",0,2,2,0)
            +-0.000866*ewgc("Cqd1R",2,2,0,0)+-0.000177*ewgc("Cqd1R",2,2,1,1)+0.030033*ewgc("Cqq1R",0,0,2,2)
            +-0.000533*ewgc("Cqq1R",1,1,2,2)+0.040348*ewgc("Cqq3R",0,0,2,2)+-0.000078*ewgc("Cqq3R",1,1,2,2)
            +0.011868*ewgc("Cqu1R",0,0,2,2)+0.000311*ewgc("Cqu1R",1,1,2,2)+0.007404*ewgc("Cqu1R",2,2,0,0)
            +0.000365*ewgc("Cqu1R",2,2,1,1)+-0.000412*ewgc("Cud1R",2,2,0,0)+0.000135*ewgc("Cud1R",2,2,1,1)
            +0.014596*ewgc("CuuR",0,0,2,2)+-0.000309*ewgc("CuuR",1,1,2,2);
            double sigma_neg_bin_mtt_350_450_madgraph = 1.674510;
            double sigma_neg_bin_mtt_350_450_NP = -0.007566*ewgc("CG")+0.014887*ewgc("Cqd8R",2,2,0,0)
            +0.305105*ewgc("Cqq1R",0,2,2,0)+0.418602*ewgc("Cqq3R",0,2,2,0)+0.002955*ewgc("Cqq3R",1,2,2,1)
            +0.085906*ewgc("Cqu8R",0,0,2,2)+0.073521*ewgc("Cqu8R",2,2,0,0)+0.015697*ewgc("Cud8R",2,2,0,0)
            -0.383667*ewgc("CuGR",2,2)+0.300446*ewgc("CuuR",0,2,2,0)
            -0.001188*ewgc("Cqd1R",2,2,0,0)+0.000242*ewgc("Cqd1R",2,2,1,1)+0.028715*ewgc("Cqq1R",0,0,2,2)
            +0.000588*ewgc("Cqq1R",1,1,2,2)+0.041531*ewgc("Cqq3R",0,0,2,2)+0.000146*ewgc("Cqq3R",1,1,2,2)
            +0.013026*ewgc("Cqu1R",0,0,2,2)+-0.000294*ewgc("Cqu1R",1,1,2,2)+0.007109*ewgc("Cqu1R",2,2,0,0)
            +-0.000347*ewgc("Cqu1R",2,2,1,1)+-0.001721*ewgc("Cud1R",2,2,0,0)+-0.000067*ewgc("Cud1R",2,2,1,1)
            +0.016846*ewgc("CuuR",0,0,2,2)+0.000338*ewgc("CuuR",1,1,2,2);
            
            double sigma_pos_bin_mtt_350_450_NP_Corrected = SM_sigma_pos_bin_mtt_350_450*sigma_pos_bin_mtt_350_450_NP/sigma_pos_bin_mtt_350_450_madgraph;
            double sigma_neg_bin_mtt_350_450_NP_Corrected = SM_sigma_neg_bin_mtt_350_450*sigma_neg_bin_mtt_350_450_NP/sigma_neg_bin_mtt_350_450_madgraph;
            
            double NP_FB_asymmetry_num_bin_mtt_350_450 = sigma_pos_bin_mtt_350_450_NP_Corrected - sigma_neg_bin_mtt_350_450_NP_Corrected;
            double NP_FB_asymmetry_deno_bin_mtt_350_450 = sigma_pos_bin_mtt_350_450_NP_Corrected + sigma_neg_bin_mtt_350_450_NP_Corrected;;
            
            return  SM_FB_asymmetry_bin_mtt_350_450*(1+(NP_FB_asymmetry_num_bin_mtt_350_450-NP_FB_asymmetry_deno_bin_mtt_350_450)/SM_FB_asymmetry_deno_bin_mtt_350_450);            
        
        }
        
    } else if(b_min == 450 && b_max == 550){
        
        
        double SM_FB_asymmetry_bin_mtt_450_550 = SM.getOptionalParameter("SM_FB_asymmetry_bin_mtt_450_550");
        double SM_FB_asymmetry_deno_bin_mtt_450_550 = SM.getOptionalParameter("SM_FB_asymmetry_deno_bin_mtt_450_550");
        double SM_FB_asymmetry_num_bin_mtt_450_550 = SM_FB_asymmetry_bin_mtt_450_550*SM_FB_asymmetry_deno_bin_mtt_450_550;
        
        double SM_sigma_pos_bin_mtt_450_550 =0.5*(SM_FB_asymmetry_num_bin_mtt_450_550+SM_FB_asymmetry_deno_bin_mtt_450_550);
        double SM_sigma_neg_bin_mtt_450_550 =0.5*(SM_FB_asymmetry_deno_bin_mtt_450_550-SM_FB_asymmetry_num_bin_mtt_450_550);

        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            
            
            double sigma_pos_bin_mtt_450_550_madgraph = 0.683803;
            double sigma_pos_bin_mtt_450_550_NP = -0.003753*ewgc("CG")+0.007844*ewgc("Cqd8R",2,2,0,0)
            +0.201208*ewgc("Cqq1R",0,2,2,0)+0.257805*ewgc("Cqq3R",0,2,2,0)+0.000709*ewgc("Cqq3R",1,2,2,1)
            +0.056092*ewgc("Cqu8R",0,0,2,2)+0.047574*ewgc("Cqu8R",2,2,0,0)+0.008085*ewgc("Cud8R",2,2,0,0)
            -0.168317*ewgc("CuGR",2,2)+0.198755*ewgc("CuuR",0,2,2,0)
            +-0.000527*ewgc("Cqd1R",2,2,0,0)+0.000144*ewgc("Cqd1R",2,2,1,1)+0.022808*ewgc("Cqq1R",0,0,2,2)
            +0.000181*ewgc("Cqq1R",1,1,2,2)+0.030815*ewgc("Cqq3R",0,0,2,2)+0.000065*ewgc("Cqq3R",1,1,2,2)
            +0.005817*ewgc("Cqu1R",0,0,2,2)+0.000027*ewgc("Cqu1R",1,1,2,2)+0.003722*ewgc("Cqu1R",2,2,0,0)
            +0.000096*ewgc("Cqu1R",2,2,1,1)+0.000104*ewgc("Cud1R",2,2,0,0)+0.000087*ewgc("Cud1R",2,2,1,1)
            +0.013193*ewgc("CuuR",0,0,2,2)+0.000111*ewgc("CuuR",1,1,2,2);
            double sigma_neg_bin_mtt_450_550_madgraph = 0.685429;
            double sigma_neg_bin_mtt_450_550_NP = -0.004203*ewgc("CG")+0.008080*ewgc("Cqd8R",2,2,0,0)
            +0.201688*ewgc("Cqq1R",0,2,2,0)+0.257507*ewgc("Cqq3R",0,2,2,0)+0.000816*ewgc("Cqq3R",1,2,2,1)
            +0.056997*ewgc("Cqu8R",0,0,2,2)+0.047589*ewgc("Cqu8R",2,2,0,0)+0.007820*ewgc("Cud8R",2,2,0,0)
            -0.169079*ewgc("CuGR",2,2)+0.198609*ewgc("CuuR",0,2,2,0)
            +-0.000152*ewgc("Cqd1R",2,2,0,0)+-0.000219*ewgc("Cqd1R",2,2,1,1)+0.022948*ewgc("Cqq1R",0,0,2,2)
            +-0.000255*ewgc("Cqq1R",1,1,2,2)+0.030379*ewgc("Cqq3R",0,0,2,2)+-0.000024*ewgc("Cqq3R",1,1,2,2)
            +0.006338*ewgc("Cqu1R",0,0,2,2)+-0.000058*ewgc("Cqu1R",1,1,2,2)+0.004314*ewgc("Cqu1R",2,2,0,0)
            +-0.000126*ewgc("Cqu1R",2,2,1,1)+0.000286*ewgc("Cud1R",2,2,0,0)+-0.000164*ewgc("Cud1R",2,2,1,1)
            +0.013542*ewgc("CuuR",0,0,2,2)+-0.000141*ewgc("CuuR",1,1,2,2);
            
            double sigma_pos_bin_mtt_450_550_NP_Corrected = SM_sigma_pos_bin_mtt_450_550*sigma_pos_bin_mtt_450_550_NP/sigma_pos_bin_mtt_450_550_madgraph;
            double sigma_neg_bin_mtt_450_550_NP_Corrected = SM_sigma_neg_bin_mtt_450_550*sigma_neg_bin_mtt_450_550_NP/sigma_neg_bin_mtt_450_550_madgraph;
            
            double NP_FB_asymmetry_num_bin_mtt_450_550 = sigma_pos_bin_mtt_450_550_NP_Corrected - sigma_neg_bin_mtt_450_550_NP_Corrected;
            double NP_FB_asymmetry_deno_bin_mtt_450_550 = sigma_pos_bin_mtt_450_550_NP_Corrected + sigma_neg_bin_mtt_450_550_NP_Corrected;;
            
            return  SM_FB_asymmetry_bin_mtt_450_550*(1+(NP_FB_asymmetry_num_bin_mtt_450_550-NP_FB_asymmetry_deno_bin_mtt_450_550)/SM_FB_asymmetry_deno_bin_mtt_450_550);            
        
        }
        
    } else if(b_min == 550 && b_max == 650){
        
        
        double SM_FB_asymmetry_bin_mtt_550_650 = SM.getOptionalParameter("SM_FB_asymmetry_bin_mtt_550_650");
        double SM_FB_asymmetry_deno_bin_mtt_550_650 = SM.getOptionalParameter("SM_FB_asymmetry_deno_bin_mtt_550_650");
        double SM_FB_asymmetry_num_bin_mtt_550_650 = SM_FB_asymmetry_bin_mtt_550_650*SM_FB_asymmetry_deno_bin_mtt_550_650;
        
        double SM_sigma_pos_bin_mtt_550_650 =0.5*(SM_FB_asymmetry_num_bin_mtt_550_650+SM_FB_asymmetry_deno_bin_mtt_550_650);
        double SM_sigma_neg_bin_mtt_550_650 =0.5*(SM_FB_asymmetry_deno_bin_mtt_550_650-SM_FB_asymmetry_num_bin_mtt_550_650);

        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            
            
            double sigma_pos_bin_mtt_550_650_madgraph = 0.216436;
            double sigma_pos_bin_mtt_550_650_NP = -0.001345*ewgc("CG")+0.002993*ewgc("Cqd8R",2,2,0,0)
            +0.098574*ewgc("Cqq1R",0,2,2,0)+0.119856*ewgc("Cqq3R",0,2,2,0)+0.000356*ewgc("Cqq3R",1,2,2,1)
            +0.026186*ewgc("Cqu8R",0,0,2,2)+0.023259*ewgc("Cqu8R",2,2,0,0)+0.002973*ewgc("Cud8R",2,2,0,0)
            -0.057346*ewgc("CuGR",2,2)+0.095602*ewgc("CuuR",0,2,2,0)
            +0.000089*ewgc("Cqd1R",2,2,0,0)+0.000160*ewgc("Cqd1R",2,2,1,1)+0.012427*ewgc("Cqq1R",0,0,2,2)
            +0.000182*ewgc("Cqq1R",1,1,2,2)+0.015602*ewgc("Cqq3R",0,0,2,2)+0.000078*ewgc("Cqq3R",1,1,2,2)
            +0.002898*ewgc("Cqu1R",0,0,2,2)+0.000184*ewgc("Cqu1R",1,1,2,2)+0.001348*ewgc("Cqu1R",2,2,0,0)
            +0.000237*ewgc("Cqu1R",2,2,1,1)+-0.000559*ewgc("Cud1R",2,2,0,0)+0.000145*ewgc("Cud1R",2,2,1,1)
            +0.007087*ewgc("CuuR",0,0,2,2)+0.000149*ewgc("CuuR",1,1,2,2);
            double sigma_neg_bin_mtt_550_650_madgraph = 0.216949;
            double sigma_neg_bin_mtt_550_650_NP = -0.001414*ewgc("CG")+0.002846*ewgc("Cqd8R",2,2,0,0)
            +0.098612*ewgc("Cqq1R",0,2,2,0)+0.120103*ewgc("Cqq3R",0,2,2,0)+0.000216*ewgc("Cqq3R",1,2,2,1)
            +0.025800*ewgc("Cqu8R",0,0,2,2)+0.023527*ewgc("Cqu8R",2,2,0,0)+0.002835*ewgc("Cud8R",2,2,0,0)
            -0.057275*ewgc("CuGR",2,2)+0.095612*ewgc("CuuR",0,2,2,0)
            +-0.000225*ewgc("Cqd1R",2,2,0,0)+-0.000155*ewgc("Cqd1R",2,2,1,1)+0.012485*ewgc("Cqq1R",0,0,2,2)
            +-0.000175*ewgc("Cqq1R",1,1,2,2)+0.015723*ewgc("Cqq3R",0,0,2,2)+-0.000069*ewgc("Cqq3R",1,1,2,2)
            +0.002745*ewgc("Cqu1R",0,0,2,2)+-0.000178*ewgc("Cqu1R",1,1,2,2)+0.001266*ewgc("Cqu1R",2,2,0,0)
            +-0.000230*ewgc("Cqu1R",2,2,1,1)+-0.000703*ewgc("Cud1R",2,2,0,0)+-0.000139*ewgc("Cud1R",2,2,1,1)
            +0.006972*ewgc("CuuR",0,0,2,2)+-0.000143*ewgc("CuuR",1,1,2,2);
            
            double sigma_pos_bin_mtt_550_650_NP_Corrected = SM_sigma_pos_bin_mtt_550_650*sigma_pos_bin_mtt_550_650_NP/sigma_pos_bin_mtt_550_650_madgraph;
            double sigma_neg_bin_mtt_550_650_NP_Corrected = SM_sigma_neg_bin_mtt_550_650*sigma_neg_bin_mtt_550_650_NP/sigma_neg_bin_mtt_550_650_madgraph;
            
            double NP_FB_asymmetry_num_bin_mtt_550_650 = sigma_pos_bin_mtt_550_650_NP_Corrected - sigma_neg_bin_mtt_550_650_NP_Corrected;
            double NP_FB_asymmetry_deno_bin_mtt_550_650 = sigma_pos_bin_mtt_550_650_NP_Corrected + sigma_neg_bin_mtt_550_650_NP_Corrected;;
            
            return  SM_FB_asymmetry_bin_mtt_550_650*(1+(NP_FB_asymmetry_num_bin_mtt_550_650-NP_FB_asymmetry_deno_bin_mtt_550_650)/SM_FB_asymmetry_deno_bin_mtt_550_650);            
        
        }
        
    } else if(b_min == 650 && b_max == 750){
        
        
        double SM_FB_asymmetry_bin_mtt_650_750 = SM.getOptionalParameter("SM_FB_asymmetry_bin_mtt_650_750");
        double SM_FB_asymmetry_deno_bin_mtt_650_750 = SM.getOptionalParameter("SM_FB_asymmetry_deno_bin_mtt_650_750");
        double SM_FB_asymmetry_num_bin_mtt_650_750 = SM_FB_asymmetry_bin_mtt_650_750*SM_FB_asymmetry_deno_bin_mtt_650_750;
        
        double SM_sigma_pos_bin_mtt_650_750 =0.5*(SM_FB_asymmetry_num_bin_mtt_650_750+SM_FB_asymmetry_deno_bin_mtt_650_750);
        double SM_sigma_neg_bin_mtt_650_750 =0.5*(SM_FB_asymmetry_deno_bin_mtt_650_750-SM_FB_asymmetry_num_bin_mtt_650_750);

        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            
            
            double sigma_pos_bin_mtt_650_750_madgraph = 0.097727;
            double sigma_pos_bin_mtt_650_750_NP = -0.000271*ewgc("CG")+0.001360*ewgc("Cqd8R",2,2,0,0)
            +0.073375*ewgc("Cqq1R",0,2,2,0)+0.084010*ewgc("Cqq3R",0,2,2,0)+0.000196*ewgc("Cqq3R",1,2,2,1)
            +0.018728*ewgc("Cqu8R",0,0,2,2)+0.017419*ewgc("Cqu8R",2,2,0,0)+0.001357*ewgc("Cud8R",2,2,0,0)
            -0.026810*ewgc("CuGR",2,2)+0.071319*ewgc("CuuR",0,2,2,0)
            +-0.000015*ewgc("Cqd1R",2,2,0,0)+0.000006*ewgc("Cqd1R",2,2,1,1)+0.010420*ewgc("Cqq1R",0,0,2,2)
            +0.000021*ewgc("Cqq1R",1,1,2,2)+0.012135*ewgc("Cqq3R",0,0,2,2)+0.000029*ewgc("Cqq3R",1,1,2,2)
            +0.001508*ewgc("Cqu1R",0,0,2,2)+0.000036*ewgc("Cqu1R",1,1,2,2)+0.001196*ewgc("Cqu1R",2,2,0,0)
            +0.000014*ewgc("Cqu1R",2,2,1,1)+-0.000041*ewgc("Cud1R",2,2,0,0)+0.000016*ewgc("Cud1R",2,2,1,1)
            +0.005619*ewgc("CuuR",0,0,2,2)+0.000034*ewgc("CuuR",1,1,2,2);
            double sigma_neg_bin_mtt_650_750_madgraph = 0.097766;
            double sigma_neg_bin_mtt_650_750_NP = -0.000309*ewgc("CG")+0.001327*ewgc("Cqd8R",2,2,0,0)
            +0.073361*ewgc("Cqq1R",0,2,2,0)+0.083778*ewgc("Cqq3R",0,2,2,0)+0.000188*ewgc("Cqq3R",1,2,2,1)
            +0.018756*ewgc("Cqu8R",0,0,2,2)+0.017360*ewgc("Cqu8R",2,2,0,0)+0.001327*ewgc("Cud8R",2,2,0,0)
            -0.026852*ewgc("CuGR",2,2)+0.071283*ewgc("CuuR",0,2,2,0)
            +-0.000003*ewgc("Cqd1R",2,2,0,0)+0.000000*ewgc("Cqd1R",2,2,1,1)+0.010473*ewgc("Cqq1R",0,0,2,2)
            +-0.000013*ewgc("Cqq1R",1,1,2,2)+0.012141*ewgc("Cqq3R",0,0,2,2)+-0.000012*ewgc("Cqq3R",1,1,2,2)
            +0.001523*ewgc("Cqu1R",0,0,2,2)+-0.000029*ewgc("Cqu1R",1,1,2,2)+0.001164*ewgc("Cqu1R",2,2,0,0)
            +-0.000007*ewgc("Cqu1R",2,2,1,1)+-0.000073*ewgc("Cud1R",2,2,0,0)+-0.000009*ewgc("Cud1R",2,2,1,1)
            +0.005643*ewgc("CuuR",0,0,2,2)+-0.000027*ewgc("CuuR",1,1,2,2);
            
            double sigma_pos_bin_mtt_650_750_NP_Corrected = SM_sigma_pos_bin_mtt_650_750*sigma_pos_bin_mtt_650_750_NP/sigma_pos_bin_mtt_650_750_madgraph;
            double sigma_neg_bin_mtt_650_750_NP_Corrected = SM_sigma_neg_bin_mtt_650_750*sigma_neg_bin_mtt_650_750_NP/sigma_neg_bin_mtt_650_750_madgraph;
            
            double NP_FB_asymmetry_num_bin_mtt_650_750 = sigma_pos_bin_mtt_650_750_NP_Corrected - sigma_neg_bin_mtt_650_750_NP_Corrected;
            double NP_FB_asymmetry_deno_bin_mtt_650_750 = sigma_pos_bin_mtt_650_750_NP_Corrected + sigma_neg_bin_mtt_650_750_NP_Corrected;;
            
            return  SM_FB_asymmetry_bin_mtt_650_750*(1+(NP_FB_asymmetry_num_bin_mtt_650_750-NP_FB_asymmetry_deno_bin_mtt_650_750)/SM_FB_asymmetry_deno_bin_mtt_650_750);            
        
        }
        
    } else {
        throw std::runtime_error("\nERROR: Please specify a correct binning range for FB_asymmetry_Tevatron_tt_diff_mtt_LO.\n");
    }
    
    
}






////// ttbar differential cross section ///////////////////////////////////


//binning of CMS 1811.06625
sigma_tt_diff_mtt_LO_CMS_181106625::sigma_tt_diff_mtt_LO_CMS_181106625(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{

    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tt_bin_mtt_300_380_CMS_181106625"
            << "SM_sigma_tt_bin_mtt_380_470_CMS_181106625" << "SM_sigma_tt_bin_mtt_470_620_CMS_181106625"
            << "SM_sigma_tt_bin_mtt_620_820_CMS_181106625" << "SM_sigma_tt_bin_mtt_820_1100_CMS_181106625"
            << "SM_sigma_tt_bin_mtt_1100_1500_CMS_181106625" << "SM_sigma_tt_bin_mtt_1500_2500_CMS_181106625");
            
    
}

double sigma_tt_diff_mtt_LO_CMS_181106625::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
   
    
    if(b_min == 300 && b_max == 380){
        
        double SM_sigma_tt_bin_300_380 = SM.getOptionalParameter("SM_sigma_tt_bin_mtt_300_380_CMS_181106625");
        double sigma_tt_bin_300_380_madgraph = 54.003200;//pb
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tt_bin_300_380 + (+-1.080720*ewgc("CG")+0.015822*ewgc("Cqd1R",2,2,0,0)
            +0.027498*ewgc("Cqd1R",2,2,1,1)+0.239796*ewgc("Cqd8R",2,2,0,0)+0.073740*ewgc("Cqd8R",2,2,1,1)
            +0.089304*ewgc("Cqq1R",0,0,2,2)+1.256970*ewgc("Cqq1R",0,2,2,0)+0.033486*ewgc("Cqq1R",1,1,2,2)
            +0.108318*ewgc("Cqq1R",1,2,2,1)+0.230184*ewgc("Cqq3R",0,0,2,2)+2.777490*ewgc("Cqq3R",0,2,2,0)
            +0.067932*ewgc("Cqq3R",1,1,2,2)+0.457110*ewgc("Cqq3R",1,2,2,1)+0.071850*ewgc("Cqu1R",0,0,2,2)
            +0.029322*ewgc("Cqu1R",1,1,2,2)+0.080706*ewgc("Cqu1R",2,2,0,0)+0.034674*ewgc("Cqu1R",2,2,1,1)
            +0.535692*ewgc("Cqu8R",0,0,2,2)+0.091128*ewgc("Cqu8R",1,1,2,2)+0.347364*ewgc("Cqu8R",2,2,0,0)
            +0.057720*ewgc("Cqu8R",2,2,1,1)+0.014940*ewgc("Cud1R",2,2,0,0)+0.030534*ewgc("Cud1R",2,2,1,1)
            +0.241590*ewgc("Cud8R",2,2,0,0)+0.076716*ewgc("Cud8R",2,2,1,1)+-18.341500*ewgc("CuGR",2,2)
            +0.081258*ewgc("CuuR",0,0,2,2)+1.238740*ewgc("CuuR",0,2,2,0)+0.036636*ewgc("CuuR",1,1,2,2)
            +0.108228*ewgc("CuuR",1,2,2,1)
                    )*(SM_sigma_tt_bin_300_380/sigma_tt_bin_300_380_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
            
        }
        
    } else if(b_min == 380 && b_max == 470){
        
        double SM_sigma_tt_bin_380_470 = SM.getOptionalParameter("SM_sigma_tt_bin_mtt_380_470_CMS_181106625");
        double sigma_tt_bin_380_470_madgraph = 173.423000;//pb
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tt_bin_380_470 + (+-8.359080*ewgc("CG")+-0.042078*ewgc("Cqd1R",2,2,0,0)+0.004590*ewgc("Cqd1R",2,2,1,1)
            +0.649842*ewgc("Cqd8R",2,2,0,0)+0.144468*ewgc("Cqd8R",2,2,1,1)+0.260772*ewgc("Cqq1R",0,0,2,2)+4.262270*ewgc("Cqq1R",0,2,2,0)
            +-0.029244*ewgc("Cqq1R",1,1,2,2)+0.230286*ewgc("Cqq1R",1,2,2,1)+0.780282*ewgc("Cqq3R",0,0,2,2)+9.255620*ewgc("Cqq3R",0,2,2,0)
            +0.084678*ewgc("Cqq3R",1,1,2,2)+1.257790*ewgc("Cqq3R",1,2,2,1)+0.071382*ewgc("Cqu1R",0,0,2,2)+0.007008*ewgc("Cqu1R",1,1,2,2)
            +0.097896*ewgc("Cqu1R",2,2,0,0)+0.018918*ewgc("Cqu1R",2,2,1,1)+1.653840*ewgc("Cqu8R",0,0,2,2)+0.194502*ewgc("Cqu8R",1,1,2,2)
            +1.000630*ewgc("Cqu8R",2,2,0,0)+0.057138*ewgc("Cqu8R",2,2,1,1)+-0.038910*ewgc("Cud1R",2,2,0,0)+0.004434*ewgc("Cud1R",2,2,1,1)
            +0.648138*ewgc("Cud8R",2,2,0,0)+0.144918*ewgc("Cud8R",2,2,1,1)+-51.556900*ewgc("CuGR",2,2)+0.308838*ewgc("CuuR",0,0,2,2)
            +4.176160*ewgc("CuuR",0,2,2,0)+0.030888*ewgc("CuuR",1,1,2,2)+0.226812*ewgc("CuuR",1,2,2,1)
                    )*(SM_sigma_tt_bin_380_470/sigma_tt_bin_380_470_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 470 && b_max == 620){
        
        double SM_sigma_tt_bin_470_620 = SM.getOptionalParameter("SM_sigma_tt_bin_mtt_470_620_CMS_181106625");
        double sigma_tt_bin_470_620_madgraph = 152.893000;//pb
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tt_bin_470_620 + (+-12.015500*ewgc("CG")+-0.014436*ewgc("Cqd1R",2,2,0,0)+-0.010944*ewgc("Cqd1R",2,2,1,1)
            +0.738569*ewgc("Cqd8R",2,2,0,0)+0.127001*ewgc("Cqd8R",2,2,1,1)+0.294306*ewgc("Cqq1R",0,0,2,2)+4.783330*ewgc("Cqq1R",0,2,2,0)
            +-0.041070*ewgc("Cqq1R",1,1,2,2)+0.199307*ewgc("Cqq1R",1,2,2,1)+1.030540*ewgc("Cqq3R",0,0,2,2)+10.384200*ewgc("Cqq3R",0,2,2,0)
            +0.084312*ewgc("Cqq3R",1,1,2,2)+1.174270*ewgc("Cqq3R",1,2,2,1)+0.134514*ewgc("Cqu1R",0,0,2,2)+-0.007254*ewgc("Cqu1R",1,1,2,2)
            +0.095286*ewgc("Cqu1R",2,2,0,0)+0.000101*ewgc("Cqu1R",2,2,1,1)+1.890170*ewgc("Cqu8R",0,0,2,2)+0.177748*ewgc("Cqu8R",1,1,2,2)
            +1.147850*ewgc("Cqu8R",2,2,0,0)+0.048700*ewgc("Cqu8R",2,2,1,1)+-0.042132*ewgc("Cud1R",2,2,0,0)+-0.013392*ewgc("Cud1R",2,2,1,1)
            +0.739715*ewgc("Cud8R",2,2,0,0)+0.126949*ewgc("Cud8R",2,2,1,1)+-42.458000*ewgc("CuGR",2,2)+0.326724*ewgc("CuuR",0,0,2,2)
            +4.680410*ewgc("CuuR",0,2,2,0)+0.007338*ewgc("CuuR",1,1,2,2)+0.194285*ewgc("CuuR",1,2,2,1)
                    )*(SM_sigma_tt_bin_470_620/sigma_tt_bin_470_620_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 620 && b_max == 820){
        
        double SM_sigma_tt_bin_620_820 = SM.getOptionalParameter("SM_sigma_tt_bin_mtt_620_820_CMS_181106625");
        double sigma_tt_bin_620_820_madgraph = 66.954300;//pb
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tt_bin_620_820 + (+-6.708550*ewgc("CG")+-0.017060*ewgc("Cqd1R",2,2,0,0)+0.003028*ewgc("Cqd1R",2,2,1,1)
            +0.511089*ewgc("Cqd8R",2,2,0,0)+0.076736*ewgc("Cqd8R",2,2,1,1)+0.224652*ewgc("Cqq1R",0,0,2,2)+3.446780*ewgc("Cqq1R",0,2,2,0)
            +-0.023712*ewgc("Cqq1R",1,1,2,2)+0.115517*ewgc("Cqq1R",1,2,2,1)+0.804056*ewgc("Cqq3R",0,0,2,2)+7.363740*ewgc("Cqq3R",0,2,2,0)
            +0.059872*ewgc("Cqq3R",1,1,2,2)+0.658982*ewgc("Cqq3R",1,2,2,1)+0.068678*ewgc("Cqu1R",0,0,2,2)+0.006047*ewgc("Cqu1R",1,1,2,2)
            +0.049867*ewgc("Cqu1R",2,2,0,0)+0.005234*ewgc("Cqu1R",2,2,1,1)+1.331380*ewgc("Cqu8R",0,0,2,2)+0.102676*ewgc("Cqu8R",1,1,2,2)
            +0.823487*ewgc("Cqu8R",2,2,0,0)+0.029953*ewgc("Cqu8R",2,2,1,1)+-0.042581*ewgc("Cud1R",2,2,0,0)+-0.001945*ewgc("Cud1R",2,2,1,1)
            +0.511913*ewgc("Cud8R",2,2,0,0)+0.077032*ewgc("Cud8R",2,2,1,1)+-18.574200*ewgc("CuGR",2,2)+0.262504*ewgc("CuuR",0,0,2,2)
            +3.359280*ewgc("CuuR",0,2,2,0)+0.010662*ewgc("CuuR",1,1,2,2)+0.112864*ewgc("CuuR",1,2,2,1)
                    )*(SM_sigma_tt_bin_620_820/sigma_tt_bin_620_820_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 820 && b_max == 1100){
        
        double SM_sigma_tt_bin_820_1100 = SM.getOptionalParameter("SM_sigma_tt_bin_mtt_820_1100_CMS_181106625");
        double sigma_tt_bin_820_1100_madgraph = 24.163600;//pb
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tt_bin_820_1100 + (-2.802570*ewgc("CG")+-0.006901*ewgc("Cqd1R",2,2,0,0)+0.000237*ewgc("Cqd1R",2,2,1,1)
            +0.333880*ewgc("Cqd8R",2,2,0,0)+0.037791*ewgc("Cqd8R",2,2,1,1)+0.161210*ewgc("Cqq1R",0,0,2,2)+2.286240*ewgc("Cqq1R",0,2,2,0)
            +-0.012073*ewgc("Cqq1R",1,1,2,2)+0.059007*ewgc("Cqq1R",1,2,2,1)+0.566288*ewgc("Cqq3R",0,0,2,2)+4.802780*ewgc("Cqq3R",0,2,2,0)
            +0.033808*ewgc("Cqq3R",1,1,2,2)+0.333648*ewgc("Cqq3R",1,2,2,1)+0.042755*ewgc("Cqu1R",0,0,2,2)+0.002326*ewgc("Cqu1R",1,1,2,2)
            +0.035514*ewgc("Cqu1R",2,2,0,0)+0.001443*ewgc("Cqu1R",2,2,1,1)+0.875141*ewgc("Cqu8R",0,0,2,2)+0.051098*ewgc("Cqu8R",1,1,2,2)
            +0.542868*ewgc("Cqu8R",2,2,0,0)+0.013234*ewgc("Cqu8R",2,2,1,1)+-0.023700*ewgc("Cud1R",2,2,0,0)+-0.001051*ewgc("Cud1R",2,2,1,1)
            +0.332700*ewgc("Cud8R",2,2,0,0)+0.037501*ewgc("Cud8R",2,2,1,1)-6.912790*ewgc("CuGR",2,2)+0.178817*ewgc("CuuR",0,0,2,2)
            +2.227770*ewgc("CuuR",0,2,2,0)+0.006347*ewgc("CuuR",1,1,2,2)+0.057377*ewgc("CuuR",1,2,2,1)
                    )*(SM_sigma_tt_bin_820_1100/sigma_tt_bin_820_1100_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 1100 && b_max == 1500){
        
        double SM_sigma_tt_bin_1100_1500 = SM.getOptionalParameter("SM_sigma_tt_bin_mtt_1100_1500_CMS_181106625");
        double sigma_tt_bin_1100_1500_madgraph = 6.983800;//pb
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tt_bin_1100_1500 + (+-0.863200*ewgc("CG")+-0.004510*ewgc("Cqd1R",2,2,0,0)+0.000022*ewgc("Cqd1R",2,2,1,1)
            +0.191513*ewgc("Cqd8R",2,2,0,0)+0.016414*ewgc("Cqd8R",2,2,1,1)+0.103366*ewgc("Cqq1R",0,0,2,2)+1.351660*ewgc("Cqq1R",0,2,2,0)
            +-0.005690*ewgc("Cqq1R",1,1,2,2)+0.026285*ewgc("Cqq1R",1,2,2,1)+0.355618*ewgc("Cqq3R",0,0,2,2)+2.795780*ewgc("Cqq3R",0,2,2,0)
            +0.014960*ewgc("Cqq3R",1,1,2,2)+0.147160*ewgc("Cqq3R",1,2,2,1)+0.024812*ewgc("Cqu1R",0,0,2,2)+0.001067*ewgc("Cqu1R",1,1,2,2)
            +0.017626*ewgc("Cqu1R",2,2,0,0)+0.000632*ewgc("Cqu1R",2,2,1,1)+0.509667*ewgc("Cqu8R",0,0,2,2)+0.022537*ewgc("Cqu8R",1,1,2,2)
            +0.318022*ewgc("Cqu8R",2,2,0,0)+0.006561*ewgc("Cqu8R",2,2,1,1)+-0.016528*ewgc("Cud1R",2,2,0,0)+-0.001045*ewgc("Cud1R",2,2,1,1)
            +0.191209*ewgc("Cud8R",2,2,0,0)+0.016353*ewgc("Cud8R",2,2,1,1)+-2.105990*ewgc("CuGR",2,2)+0.112485*ewgc("CuuR",0,0,2,2)
            +1.312690*ewgc("CuuR",0,2,2,0)+0.002397*ewgc("CuuR",1,1,2,2)+0.025537*ewgc("CuuR",1,2,2,1)
                    )*(SM_sigma_tt_bin_1100_1500/sigma_tt_bin_1100_1500_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 1500 && b_max == 2500){
        
        double SM_sigma_tt_bin_1500_2500 = SM.getOptionalParameter("SM_sigma_tt_bin_mtt_1500_2500_CMS_181106625");
        double sigma_tt_bin_1500_2500_madgraph = 1.821290;//pb
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tt_bin_1500_2500 + (+-0.231926*ewgc("CG")+-0.003089*ewgc("Cqd1R",2,2,0,0)+-0.000047*ewgc("Cqd1R",2,2,1,1)
            +0.130652*ewgc("Cqd8R",2,2,0,0)+0.007679*ewgc("Cqd8R",2,2,1,1)+0.078637*ewgc("Cqq1R",0,0,2,2)+0.983935*ewgc("Cqq1R",0,2,2,0)
            +-0.002654*ewgc("Cqq1R",1,1,2,2)+0.014376*ewgc("Cqq1R",1,2,2,1)+0.252179*ewgc("Cqq3R",0,0,2,2)+1.968040*ewgc("Cqq3R",0,2,2,0)
            +0.007838*ewgc("Cqq3R",1,1,2,2)+0.071898*ewgc("Cqq3R",1,2,2,1)+0.017883*ewgc("Cqu1R",0,0,2,2)+0.000657*ewgc("Cqu1R",1,1,2,2)
            +0.011678*ewgc("Cqu1R",2,2,0,0)+0.000359*ewgc("Cqu1R",2,2,1,1)+0.363022*ewgc("Cqu8R",0,0,2,2)+0.011088*ewgc("Cqu8R",1,1,2,2)
            +0.232613*ewgc("Cqu8R",2,2,0,0)+0.003450*ewgc("Cqu8R",2,2,1,1)+-0.011469*ewgc("Cud1R",2,2,0,0)+-0.000549*ewgc("Cud1R",2,2,1,1)
            +0.130711*ewgc("Cud8R",2,2,0,0)+0.007725*ewgc("Cud8R",2,2,1,1)+-0.607749*ewgc("CuGR",2,2)+0.081595*ewgc("CuuR",0,0,2,2)
            +0.955959*ewgc("CuuR",0,2,2,0)+0.001335*ewgc("CuuR",1,1,2,2)+0.014025*ewgc("CuuR",1,2,2,1)
                    )*(SM_sigma_tt_bin_1500_2500/sigma_tt_bin_1500_2500_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } 
    
}






//binning of CMS 2108.02803
sigma_tt_diff_mtt_CMS_LO::sigma_tt_diff_mtt_CMS_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{

    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tt_CMS_bin_mtt_250_400" << "SM_sigma_tt_CMS_bin_mtt_400_480" <<
            "SM_sigma_tt_CMS_bin_mtt_480_560" << "SM_sigma_tt_CMS_bin_mtt_560_640" << "SM_sigma_tt_CMS_bin_mtt_640_720" << 
            "SM_sigma_tt_CMS_bin_mtt_720_800" << "SM_sigma_tt_CMS_bin_mtt_800_900" << "SM_sigma_tt_CMS_bin_mtt_900_1000" << 
            "SM_sigma_tt_CMS_bin_mtt_1000_1150" << "SM_sigma_tt_CMS_bin_mtt_1150_1300" << "SM_sigma_tt_CMS_bin_mtt_1300_1500" <<
            "SM_sigma_tt_CMS_bin_mtt_1500_1700" << "SM_sigma_tt_CMS_bin_mtt_1700_2000" << "SM_sigma_tt_CMS_bin_mtt_2000_2300" <<
            "SM_sigma_tt_CMS_bin_mtt_2300_3500");
            //"SM_sigma_tt_bin_2300_2600" << "SM_sigma_tt_bin_2600_3000" << "SM_sigma_tt_bin_3000_3500" << "SM_sigma_tt_bin_3500_4000");
    
}

double sigma_tt_diff_mtt_CMS_LO::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    
    //double C_tG = ewgc("CuGR",2,2);

    
    
    
    if(b_min == 250 && b_max == 400){
        
        double SM_sigma_tt_bin_250_400 = SM.getOptionalParameter("SM_sigma_tt_CMS_bin_mtt_250_400");
        //double sigma_tt_bin_250_400_madgraph = 105600.0;//fb maybe over the bin width? Check!
        double sigma_tt_bin_250_400_madgraph = 95.612100;//pb
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            

            
            double total = SM_sigma_tt_bin_250_400 + (-2.422920*ewgc("CG")+0.401694*ewgc("Cqd8R",2,2,0,0)+0.123768*ewgc("Cqd8R",2,2,1,1)
                    +2.248510*ewgc("Cqq1R",0,2,2,0)+0.184272*ewgc("Cqq1R",1,2,2,1)+4.951130*ewgc("Cqq3R",0,2,2,0)
                    +0.780684*ewgc("Cqq3R",1,2,2,1)+0.926298*ewgc("Cqu8R",0,0,2,2)+0.152874*ewgc("Cqu8R",1,1,2,2)
                    +0.575544*ewgc("Cqu8R",2,2,0,0)+0.079620*ewgc("Cqu8R",2,2,1,1)+0.404088*ewgc("Cud8R",2,2,0,0)
                    +0.127740*ewgc("Cud8R",2,2,1,1)-31.370400*ewgc("CuGR",2,2)+2.216290*ewgc("CuuR",0,2,2,0)
                    +0.180834*ewgc("CuuR",1,2,2,1)+
                    0.007578*ewgc("Cqd1R",2,2,0,0)+0.042222*ewgc("Cqd1R",2,2,1,1)+0.159660*ewgc("Cqq1R",0,0,2,2)
                    +0.024942*ewgc("Cqq1R",1,1,2,2)+0.397038*ewgc("Cqq3R",0,0,2,2)+0.088098*ewgc("Cqq3R",1,1,2,2)
                    +0.100374*ewgc("Cqu1R",0,0,2,2)+0.039066*ewgc("Cqu1R",1,1,2,2)+0.110028*ewgc("Cqu1R",2,2,0,0)
                    +0.051414*ewgc("Cqu1R",2,2,1,1)+0.023658*ewgc("Cud1R",2,2,0,0)+0.044238*ewgc("Cud1R",2,2,1,1)
                    +0.146250*ewgc("CuuR",0,0,2,2)+0.059376*ewgc("CuuR",1,1,2,2)
                    )*(SM_sigma_tt_bin_250_400/sigma_tt_bin_250_400_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
            
        }
        
    } else if(b_min == 400 && b_max == 480){
        
        double SM_sigma_tt_bin_400_480 = SM.getOptionalParameter("SM_sigma_tt_CMS_bin_mtt_400_480");
        double sigma_tt_bin_400_480_madgraph = 147.218000;//pb
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tt_bin_400_480 + (-8.063800*ewgc("CG")+0.537774*ewgc("Cqd8R",2,2,0,0)+0.102036*ewgc("Cqd8R",2,2,1,1)
                    +3.659260*ewgc("Cqq1R",0,2,2,0)+0.168744*ewgc("Cqq1R",1,2,2,1)+7.940650*ewgc("Cqq3R",0,2,2,0)
                    +1.036730*ewgc("Cqq3R",1,2,2,1)+1.420910*ewgc("Cqu8R",0,0,2,2)+0.147192*ewgc("Cqu8R",1,1,2,2)
                    +0.859650*ewgc("Cqu8R",2,2,0,0)+0.036930*ewgc("Cqu8R",2,2,1,1)+0.537786*ewgc("Cud8R",2,2,0,0)
                    +0.101004*ewgc("Cud8R",2,2,1,1)-42.879900*ewgc("CuGR",2,2)+3.585550*ewgc("CuuR",0,2,2,0)
                    +0.166566*ewgc("CuuR",1,2,2,1)
                    -0.030360*ewgc("Cqd1R",2,2,0,0)+-0.014292*ewgc("Cqd1R",2,2,1,1)+0.202368*ewgc("Cqq1R",0,0,2,2)
                    +-0.028410*ewgc("Cqq1R",1,1,2,2)+0.682230*ewgc("Cqq3R",0,0,2,2)+0.060234*ewgc("Cqq3R",1,1,2,2)
                    +0.056178*ewgc("Cqu1R",0,0,2,2)+-0.007908*ewgc("Cqu1R",1,1,2,2)+0.073824*ewgc("Cqu1R",2,2,0,0)
                    +-0.001914*ewgc("Cqu1R",2,2,1,1)+-0.048564*ewgc("Cud1R",2,2,0,0)+-0.013968*ewgc("Cud1R",2,2,1,1)
                    +0.267738*ewgc("CuuR",0,0,2,2)+0.001314*ewgc("CuuR",1,1,2,2)
                    )*(SM_sigma_tt_bin_400_480/sigma_tt_bin_400_480_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 480 && b_max == 560){
        
        double SM_sigma_tt_bin_480_560 = SM.getOptionalParameter("SM_sigma_tt_CMS_bin_mtt_480_560");
        double sigma_tt_bin_480_560_madgraph = 93.593900;//pb
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tt_bin_480_560+(-7.025760*ewgc("CG")+0.449310*ewgc("Cqd8R",2,2,0,0)+0.079662*ewgc("Cqd8R",2,2,1,1)
                    +2.769450*ewgc("Cqq1R",0,2,2,0)+0.129120*ewgc("Cqq1R",1,2,2,1)+6.058480*ewgc("Cqq3R",0,2,2,0)
                    +0.711270*ewgc("Cqq3R",1,2,2,1)+1.107320*ewgc("Cqu8R",0,0,2,2)+0.113052*ewgc("Cqu8R",1,1,2,2)
                    +0.676488*ewgc("Cqu8R",2,2,0,0)+0.029088*ewgc("Cqu8R",2,2,1,1)+0.448860*ewgc("Cud8R",2,2,0,0)
                    +0.080076*ewgc("Cud8R",2,2,1,1)-26.003200*ewgc("CuGR",2,2)+2.714080*ewgc("CuuR",0,2,2,0)
                    +0.128328*ewgc("CuuR",1,2,2,1)
                    -0.010302*ewgc("Cqd1R",2,2,0,0)+-0.001752*ewgc("Cqd1R",2,2,1,1)+0.150552*ewgc("Cqq1R",0,0,2,2)
                    +-0.018192*ewgc("Cqq1R",1,1,2,2)+0.575952*ewgc("Cqq3R",0,0,2,2)+0.060726*ewgc("Cqq3R",1,1,2,2)
                    +0.086460*ewgc("Cqu1R",0,0,2,2)+0.002484*ewgc("Cqu1R",1,1,2,2)+0.062640*ewgc("Cqu1R",2,2,0,0)
                    +0.001020*ewgc("Cqu1R",2,2,1,1)+-0.025974*ewgc("Cud1R",2,2,0,0)+-0.005118*ewgc("Cud1R",2,2,1,1)
                    +0.175404*ewgc("CuuR",0,0,2,2)+0.010116*ewgc("CuuR",1,1,2,2)
                    )*(SM_sigma_tt_bin_480_560/sigma_tt_bin_480_560_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 560 && b_max == 640){
        
        double SM_sigma_tt_bin_560_640 = SM.getOptionalParameter("SM_sigma_tt_CMS_bin_mtt_560_640");
        double sigma_tt_bin_560_640_madgraph = 55.041700;//pb
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tt_bin_560_640+(-4.951430*ewgc("CG")+0.305342*ewgc("Cqd8R",2,2,0,0)+0.048212*ewgc("Cqd8R",2,2,1,1)
                    +2.089870*ewgc("Cqq1R",0,2,2,0)+0.073246*ewgc("Cqq1R",1,2,2,1)+4.469060*ewgc("Cqq3R",0,2,2,0)
                    +0.460417*ewgc("Cqq3R",1,2,2,1)+0.799277*ewgc("Cqu8R",0,0,2,2)+0.062208*ewgc("Cqu8R",1,1,2,2)
                    +0.489781*ewgc("Cqu8R",2,2,0,0)+0.018157*ewgc("Cqu8R",2,2,1,1)+0.304480*ewgc("Cud8R",2,2,0,0)
                    +0.048640*ewgc("Cud8R",2,2,1,1)-15.174600*ewgc("CuGR",2,2)+2.032110*ewgc("CuuR",0,2,2,0)
                    +0.070646*ewgc("CuuR",1,2,2,1)
                    -0.015037*ewgc("Cqd1R",2,2,0,0)+-0.004797*ewgc("Cqd1R",2,2,1,1)+0.156666*ewgc("Cqq1R",0,0,2,2)
                    +-0.019126*ewgc("Cqq1R",1,1,2,2)+0.483078*ewgc("Cqq3R",0,0,2,2)+0.035353*ewgc("Cqq3R",1,1,2,2)
                    +0.040200*ewgc("Cqu1R",0,0,2,2)+-0.004762*ewgc("Cqu1R",1,1,2,2)+0.026975*ewgc("Cqu1R",2,2,0,0)
                    +0.000032*ewgc("Cqu1R",2,2,1,1)+-0.024794*ewgc("Cud1R",2,2,0,0)+-0.003859*ewgc("Cud1R",2,2,1,1)
                    +0.150933*ewgc("CuuR",0,0,2,2)+0.002163*ewgc("CuuR",1,1,2,2)
                    )*(SM_sigma_tt_bin_560_640/sigma_tt_bin_560_640_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 640 && b_max == 720){
        
        double SM_sigma_tt_bin_640_720 = SM.getOptionalParameter("SM_sigma_tt_CMS_bin_mtt_640_720");
        double sigma_tt_bin_640_720_madgraph = 32.530700;//pb
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tt_bin_640_720+(-3.219570*ewgc("CG")+0.232042*ewgc("Cqd8R",2,2,0,0)+0.037100*ewgc("Cqd8R",2,2,1,1)
                    +1.550940*ewgc("Cqq1R",0,2,2,0)+0.053066*ewgc("Cqq1R",1,2,2,1)+3.301500*ewgc("Cqq3R",0,2,2,0)
                    +0.307033*ewgc("Cqq3R",1,2,2,1)+0.604172*ewgc("Cqu8R",0,0,2,2)+0.049607*ewgc("Cqu8R",1,1,2,2)
                    +0.375725*ewgc("Cqu8R",2,2,0,0)+0.014703*ewgc("Cqu8R",2,2,1,1)+0.232665*ewgc("Cud8R",2,2,0,0)
                    +0.036439*ewgc("Cud8R",2,2,1,1)-8.998020*ewgc("CuGR",2,2)+1.513220*ewgc("CuuR",0,2,2,0)
                    +0.052060*ewgc("CuuR",1,2,2,1)
                    -0.002628*ewgc("Cqd1R",2,2,0,0)+0.001228*ewgc("Cqd1R",2,2,1,1)+0.109033*ewgc("Cqq1R",0,0,2,2)
                    +-0.011971*ewgc("Cqq1R",1,1,2,2)+0.357790*ewgc("Cqq3R",0,0,2,2)+0.027758*ewgc("Cqq3R",1,1,2,2)
                    +0.031403*ewgc("Cqu1R",0,0,2,2)+0.002604*ewgc("Cqu1R",1,1,2,2)+0.030589*ewgc("Cqu1R",2,2,0,0)
                    +0.004275*ewgc("Cqu1R",2,2,1,1)+-0.015530*ewgc("Cud1R",2,2,0,0)+-0.001006*ewgc("Cud1R",2,2,1,1)
                    +0.127340*ewgc("CuuR",0,0,2,2)+0.005749*ewgc("CuuR",1,1,2,2)
                    )*(SM_sigma_tt_bin_640_720/sigma_tt_bin_640_720_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
            
        }
        
    } else if(b_min == 720 && b_max == 800){
        
        double SM_sigma_tt_bin_720_800 = SM.getOptionalParameter("SM_sigma_tt_CMS_bin_mtt_720_800");
        double sigma_tt_bin_720_800_madgraph = 19.672800;//pb
        
        //PROBABLY THE OPERATORS PURELY FROM THE THIRD FAMILY NEED TO BE NEGLECTED SINCE THEIR CONTRIBUTION IS REALLY SMALL
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tt_bin_720_800+(-2.083310*ewgc("CG")+0.175894*ewgc("Cqd8R",2,2,0,0)+0.025732*ewgc("Cqd8R",2,2,1,1)
                    +1.182780*ewgc("Cqq1R",0,2,2,0)+0.038678*ewgc("Cqq1R",1,2,2,1)+2.528020*ewgc("Cqq3R",0,2,2,0)
                    +0.209536*ewgc("Cqq3R",1,2,2,1)+0.457521*ewgc("Cqu8R",0,0,2,2)+0.034597*ewgc("Cqu8R",1,1,2,2)
                    +0.282292*ewgc("Cqu8R",2,2,0,0)+0.011773*ewgc("Cqu8R",2,2,1,1)+0.176461*ewgc("Cud8R",2,2,0,0)
                    +0.026183*ewgc("Cud8R",2,2,1,1)-5.490750*ewgc("CuGR",2,2)+1.152010*ewgc("CuuR",0,2,2,0)
                    +0.037207*ewgc("CuuR",1,2,2,1)
                    -0.005403*ewgc("Cqd1R",2,2,0,0)+0.001667*ewgc("Cqd1R",2,2,1,1)+0.074296*ewgc("Cqq1R",0,0,2,2)
                    +-0.006382*ewgc("Cqq1R",1,1,2,2)+0.282862*ewgc("Cqq3R",0,0,2,2)+0.020603*ewgc("Cqq3R",1,1,2,2)
                    +0.024782*ewgc("Cqu1R",0,0,2,2)+0.003390*ewgc("Cqu1R",1,1,2,2)+0.019139*ewgc("Cqu1R",2,2,0,0)
                    +0.003552*ewgc("Cqu1R",2,2,1,1)+-0.011411*ewgc("Cud1R",2,2,0,0)+0.000128*ewgc("Cud1R",2,2,1,1)
                    +0.087248*ewgc("CuuR",0,0,2,2)+0.006031*ewgc("CuuR",1,1,2,2)
                    )*(SM_sigma_tt_bin_720_800/sigma_tt_bin_720_800_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 800 && b_max == 900){
        
        double SM_sigma_tt_bin_800_900 = SM.getOptionalParameter("SM_sigma_tt_CMS_bin_mtt_800_900");
        double sigma_tt_bin_800_900_madgraph = 14.503700;//pb
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tt_bin_800_900+(-1.630550*ewgc("CG")+0.163468*ewgc("Cqd8R",2,2,0,0)+0.019984*ewgc("Cqd8R",2,2,1,1)
                    +1.108320*ewgc("Cqq1R",0,2,2,0)+0.030141*ewgc("Cqq1R",1,2,2,1)+2.335190*ewgc("Cqq3R",0,2,2,0)
                    +0.177932*ewgc("Cqq3R",1,2,2,1)+0.423182*ewgc("Cqu8R",0,0,2,2)+0.026737*ewgc("Cqu8R",1,1,2,2)
                    +0.263074*ewgc("Cqu8R",2,2,0,0)+0.007495*ewgc("Cqu8R",2,2,1,1)+0.162380*ewgc("Cud8R",2,2,0,0)
                    +0.019683*ewgc("Cud8R",2,2,1,1)-4.096130*ewgc("CuGR",2,2)+1.080140*ewgc("CuuR",0,2,2,0)
                    +0.029855*ewgc("CuuR",1,2,2,1)
                    -0.003016*ewgc("Cqd1R",2,2,0,0)+-0.000014*ewgc("Cqd1R",2,2,1,1)+0.082508*ewgc("Cqq1R",0,0,2,2)
                    +-0.006467*ewgc("Cqq1R",1,1,2,2)+0.284167*ewgc("Cqq3R",0,0,2,2)+0.016937*ewgc("Cqq3R",1,1,2,2)
                    +0.024673*ewgc("Cqu1R",0,0,2,2)+0.001100*ewgc("Cqu1R",1,1,2,2)+0.017210*ewgc("Cqu1R",2,2,0,0)
                    +0.000256*ewgc("Cqu1R",2,2,1,1)+-0.012320*ewgc("Cud1R",2,2,0,0)+-0.001016*ewgc("Cud1R",2,2,1,1)
                    +0.092839*ewgc("CuuR",0,0,2,2)+0.002999*ewgc("CuuR",1,1,2,2)
                    )*(SM_sigma_tt_bin_800_900/sigma_tt_bin_800_900_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 900 && b_max == 1000){
        
        double SM_sigma_tt_bin_900_1000 = SM.getOptionalParameter("SM_sigma_tt_CMS_bin_mtt_900_1000");
        double sigma_tt_bin_900_1000_madgraph = 8.319110;
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tt_bin_900_1000+(-0.968369*ewgc("CG")+0.119704*ewgc("Cqd8R",2,2,0,0)+0.014122*ewgc("Cqd8R",2,2,1,1)
                    +0.816250*ewgc("Cqq1R",0,2,2,0)+0.021400*ewgc("Cqq1R",1,2,2,1)+1.722240*ewgc("Cqq3R",0,2,2,0)
                    +0.118511*ewgc("Cqq3R",1,2,2,1)+0.313127*ewgc("Cqu8R",0,0,2,2)+0.019055*ewgc("Cqu8R",1,1,2,2)
                    +0.193743*ewgc("Cqu8R",2,2,0,0)+0.005861*ewgc("Cqu8R",2,2,1,1)+0.118769*ewgc("Cud8R",2,2,0,0)
                    +0.014256*ewgc("Cud8R",2,2,1,1)-2.388600*ewgc("CuGR",2,2)+0.792650*ewgc("CuuR",0,2,2,0)
                    +0.020774*ewgc("CuuR",1,2,2,1)
                    -0.003293*ewgc("Cqd1R",2,2,0,0)+0.000536*ewgc("Cqd1R",2,2,1,1)+0.050591*ewgc("Cqq1R",0,0,2,2)
                    +-0.004069*ewgc("Cqq1R",1,1,2,2)+0.198784*ewgc("Cqq3R",0,0,2,2)+0.012664*ewgc("Cqq3R",1,1,2,2)
                    +0.014075*ewgc("Cqu1R",0,0,2,2)+0.000829*ewgc("Cqu1R",1,1,2,2)+0.012454*ewgc("Cqu1R",2,2,0,0)
                    +0.001399*ewgc("Cqu1R",2,2,1,1)+-0.009767*ewgc("Cud1R",2,2,0,0)+-0.000207*ewgc("Cud1R",2,2,1,1)
                    +0.062268*ewgc("CuuR",0,0,2,2)+0.002893*ewgc("CuuR",1,1,2,2)
                    )*(SM_sigma_tt_bin_900_1000/sigma_tt_bin_900_1000_madgraph);
            
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 1000 && b_max == 1150){
        
        double SM_sigma_tt_bin_1000_1150 = SM.getOptionalParameter("SM_sigma_tt_CMS_bin_mtt_1000_1150");
        double sigma_tt_bin_1000_1150_madgraph = 6.639310;
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tt_bin_1000_1150+(-0.804304*ewgc("CG")+0.123780*ewgc("Cqd8R",2,2,0,0)+0.012627*ewgc("Cqd8R",2,2,1,1)
                    +0.861600*ewgc("Cqq1R",0,2,2,0)+0.019534*ewgc("Cqq1R",1,2,2,1)+1.803180*ewgc("Cqq3R",0,2,2,0)
                    +0.110947*ewgc("Cqq3R",1,2,2,1)+0.330188*ewgc("Cqu8R",0,0,2,2)+0.016923*ewgc("Cqu8R",1,1,2,2)
                    +0.204389*ewgc("Cqu8R",2,2,0,0)+0.004738*ewgc("Cqu8R",2,2,1,1)+0.124652*ewgc("Cud8R",2,2,0,0)
                    +0.012582*ewgc("Cud8R",2,2,1,1)-1.942700*ewgc("CuGR",2,2)+0.840043*ewgc("CuuR",0,2,2,0)
                    +0.018953*ewgc("CuuR",1,2,2,1)
                    -0.003310*ewgc("Cqd1R",2,2,0,0)+-0.000267*ewgc("Cqd1R",2,2,1,1)+0.067012*ewgc("Cqq1R",0,0,2,2)
                    +-0.004102*ewgc("Cqq1R",1,1,2,2)+0.217415*ewgc("Cqq3R",0,0,2,2)+0.011236*ewgc("Cqq3R",1,1,2,2)
                    +0.013693*ewgc("Cqu1R",0,0,2,2)+0.001001*ewgc("Cqu1R",1,1,2,2)+0.009635*ewgc("Cqu1R",2,2,0,0)
                    +0.000624*ewgc("Cqu1R",2,2,1,1)+-0.011255*ewgc("Cud1R",2,2,0,0)+-0.000678*ewgc("Cud1R",2,2,1,1)
                    +0.068914*ewgc("CuuR",0,0,2,2)+0.001922*ewgc("CuuR",1,1,2,2)
                    )*(SM_sigma_tt_bin_1000_1150/sigma_tt_bin_1000_1150_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
            
        }
        
    } else if(b_min == 1150 && b_max == 1300){
        
        double SM_sigma_tt_bin_1150_1300 = SM.getOptionalParameter("SM_sigma_tt_CMS_bin_mtt_1150_1300");
        double sigma_tt_bin_1150_1300_madgraph = 3.246500;
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tt_bin_1150_1300+(-0.401966*ewgc("CG")+0.083879*ewgc("Cqd8R",2,2,0,0)+0.007322*ewgc("Cqd8R",2,2,1,1)
                    +0.586193*ewgc("Cqq1R",0,2,2,0)+0.011479*ewgc("Cqq1R",1,2,2,1)+1.212270*ewgc("Cqq3R",0,2,2,0)
                    +0.066123*ewgc("Cqq3R",1,2,2,1)+0.220114*ewgc("Cqu8R",0,0,2,2)+0.009786*ewgc("Cqu8R",1,1,2,2)
                    +0.138518*ewgc("Cqu8R",2,2,0,0)+0.002758*ewgc("Cqu8R",2,2,1,1)+0.083424*ewgc("Cud8R",2,2,0,0)
                    +0.007362*ewgc("Cud8R",2,2,1,1)-0.974239*ewgc("CuGR",2,2)+0.569482*ewgc("CuuR",0,2,2,0)
                    +0.011161*ewgc("CuuR",1,2,2,1)
                    -0.001702*ewgc("Cqd1R",2,2,0,0)+-0.000120*ewgc("Cqd1R",2,2,1,1)+0.040512*ewgc("Cqq1R",0,0,2,2)
                    +-0.002375*ewgc("Cqq1R",1,1,2,2)+0.151994*ewgc("Cqq3R",0,0,2,2)+0.007076*ewgc("Cqq3R",1,1,2,2)
                    +0.013447*ewgc("Cqu1R",0,0,2,2)+0.000342*ewgc("Cqu1R",1,1,2,2)+0.008418*ewgc("Cqu1R",2,2,0,0)
                    +0.000261*ewgc("Cqu1R",2,2,1,1)+-0.006759*ewgc("Cud1R",2,2,0,0)+-0.000548*ewgc("Cud1R",2,2,1,1)
                    +0.048311*ewgc("CuuR",0,0,2,2)+0.001036*ewgc("CuuR",1,1,2,2)
                    )*(SM_sigma_tt_bin_1150_1300/sigma_tt_bin_1150_1300_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 1300 && b_max == 1500){
        
        double SM_sigma_tt_bin_1300_1500 = SM.getOptionalParameter("SM_sigma_tt_CMS_bin_mtt_1300_1500");
        double sigma_tt_bin_1300_1500_madgraph = 2.043090;
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tt_bin_1300_1500+(-0.257644*ewgc("CG")+0.071800*ewgc("Cqd8R",2,2,0,0)+0.005583*ewgc("Cqd8R",2,2,1,1)
                    +0.514089*ewgc("Cqq1R",0,2,2,0)+0.009040*ewgc("Cqq1R",1,2,2,1)+1.057920*ewgc("Cqq3R",0,2,2,0)
                    +0.049771*ewgc("Cqq3R",1,2,2,1)+0.193778*ewgc("Cqu8R",0,0,2,2)+0.007659*ewgc("Cqu8R",1,1,2,2)
                    +0.121016*ewgc("Cqu8R",2,2,0,0)+0.002184*ewgc("Cqu8R",2,2,1,1)+0.071700*ewgc("Cud8R",2,2,0,0)
                    +0.005503*ewgc("Cud8R",2,2,1,1)-0.630898*ewgc("CuGR",2,2)+0.499461*ewgc("CuuR",0,2,2,0)
                    +0.008718*ewgc("CuuR",1,2,2,1)
                    -0.001690*ewgc("Cqd1R",2,2,0,0)+0.000018*ewgc("Cqd1R",2,2,1,1)+0.040412*ewgc("Cqq1R",0,0,2,2)
                    +-0.002151*ewgc("Cqq1R",1,1,2,2)+0.135654*ewgc("Cqq3R",0,0,2,2)+0.004878*ewgc("Cqq3R",1,1,2,2)
                    +0.008707*ewgc("Cqu1R",0,0,2,2)+0.000370*ewgc("Cqu1R",1,1,2,2)+0.005981*ewgc("Cqu1R",2,2,0,0)
                    +0.000084*ewgc("Cqu1R",2,2,1,1)+-0.006194*ewgc("Cud1R",2,2,0,0)+-0.000432*ewgc("Cud1R",2,2,1,1)
                    +0.043381*ewgc("CuuR",0,0,2,2)+0.000671*ewgc("CuuR",1,1,2,2)
                    )*(SM_sigma_tt_bin_1300_1500/sigma_tt_bin_1300_1500_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 1500 && b_max == 1700){
        
        double SM_sigma_tt_bin_1500_1700 = SM.getOptionalParameter("SM_sigma_tt_CMS_bin_mtt_1500_1700");
        double sigma_tt_bin_1500_1700_madgraph = 0.914445;
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tt_bin_1500_1700+(-0.117785*ewgc("CG")+0.045019*ewgc("Cqd8R",2,2,0,0)+0.003048*ewgc("Cqd8R",2,2,1,1)
                    +0.324342*ewgc("Cqq1R",0,2,2,0)+0.005084*ewgc("Cqq1R",1,2,2,1)+0.662237*ewgc("Cqq3R",0,2,2,0)
                    +0.028217*ewgc("Cqq3R",1,2,2,1)+0.121878*ewgc("Cqu8R",0,0,2,2)+0.004270*ewgc("Cqu8R",1,1,2,2)
                    +0.077000*ewgc("Cqu8R",2,2,0,0)+0.001206*ewgc("Cqu8R",2,2,1,1)+0.045112*ewgc("Cud8R",2,2,0,0)
                    +0.003088*ewgc("Cud8R",2,2,1,1)-0.292794*ewgc("CuGR",2,2)+0.315186*ewgc("CuuR",0,2,2,0)
                    +0.004977*ewgc("CuuR",1,2,2,1)
                    -0.001201*ewgc("Cqd1R",2,2,0,0)+-0.000002*ewgc("Cqd1R",2,2,1,1)+0.021836*ewgc("Cqq1R",0,0,2,2)
                    +-0.001153*ewgc("Cqq1R",1,1,2,2)+0.081725*ewgc("Cqq3R",0,0,2,2)+0.003047*ewgc("Cqq3R",1,1,2,2)
                    +0.006181*ewgc("Cqu1R",0,0,2,2)+0.000329*ewgc("Cqu1R",1,1,2,2)+0.004044*ewgc("Cqu1R",2,2,0,0)
                    +0.000197*ewgc("Cqu1R",2,2,1,1)+-0.003978*ewgc("Cud1R",2,2,0,0)+-0.000220*ewgc("Cud1R",2,2,1,1)
                    +0.026489*ewgc("CuuR",0,0,2,2)+0.000505*ewgc("CuuR",1,1,2,2)
                    )*(SM_sigma_tt_bin_1500_1700/sigma_tt_bin_1500_1700_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 1700 && b_max == 2000){
        
        double SM_sigma_tt_bin_1700_2000 = SM.getOptionalParameter("SM_sigma_tt_CMS_bin_mtt_1700_2000");
        double sigma_tt_bin_1700_2000_madgraph = 0.562282;
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tt_bin_1700_2000+(-0.071302*ewgc("CG")+0.039170*ewgc("Cqd8R",2,2,0,0)+0.002352*ewgc("Cqd8R",2,2,1,1)
                    +0.289186*ewgc("Cqq1R",0,2,2,0)+0.004238*ewgc("Cqq1R",1,2,2,1)+0.584013*ewgc("Cqq3R",0,2,2,0)
                    +0.021609*ewgc("Cqq3R",1,2,2,1)+0.107387*ewgc("Cqu8R",0,0,2,2)+0.003355*ewgc("Cqu8R",1,1,2,2)
                    +0.068348*ewgc("Cqu8R",2,2,0,0)+0.001051*ewgc("Cqu8R",2,2,1,1)+0.039148*ewgc("Cud8R",2,2,0,0)
                    +0.002369*ewgc("Cud8R",2,2,1,1)-0.187937*ewgc("CuGR",2,2)+0.281115*ewgc("CuuR",0,2,2,0)
                    +0.004135*ewgc("CuuR",1,2,2,1)
                    -0.001093*ewgc("Cqd1R",2,2,0,0)+0.000001*ewgc("Cqd1R",2,2,1,1)+0.027035*ewgc("Cqq1R",0,0,2,2)
                    +-0.000869*ewgc("Cqq1R",1,1,2,2)+0.077631*ewgc("Cqq3R",0,0,2,2)+0.002319*ewgc("Cqq3R",1,1,2,2)
                    +0.005833*ewgc("Cqu1R",0,0,2,2)+0.000155*ewgc("Cqu1R",1,1,2,2)+0.003598*ewgc("Cqu1R",2,2,0,0)
                    +0.000098*ewgc("Cqu1R",2,2,1,1)+-0.003566*ewgc("Cud1R",2,2,0,0)+-0.000132*ewgc("Cud1R",2,2,1,1)
                    +0.023764*ewgc("CuuR",0,0,2,2)+0.000360*ewgc("CuuR",1,1,2,2)
                    )*(SM_sigma_tt_bin_1700_2000/sigma_tt_bin_1700_2000_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
            
        }
        
    } else if(b_min == 2000 && b_max == 2300){
        
        double SM_sigma_tt_bin_2000_2300 = SM.getOptionalParameter("SM_sigma_tt_CMS_bin_mtt_2000_2300");
        double sigma_tt_bin_2000_2300_madgraph = 0.204084;
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tt_bin_2000_2300+(-0.025926*ewgc("CG")+0.020907*ewgc("Cqd8R",2,2,0,0)+0.001111*ewgc("Cqd8R",2,2,1,1)
                    +0.158780*ewgc("Cqq1R",0,2,2,0)+0.002134*ewgc("Cqq1R",1,2,2,1)+0.315561*ewgc("Cqq3R",0,2,2,0)
                    +0.010496*ewgc("Cqq3R",1,2,2,1)+0.058108*ewgc("Cqu8R",0,0,2,2)+0.001605*ewgc("Cqu8R",1,1,2,2)
                    +0.037313*ewgc("Cqu8R",2,2,0,0)+0.000506*ewgc("Cqu8R",2,2,1,1)+0.020857*ewgc("Cud8R",2,2,0,0)
                    +0.001103*ewgc("Cud8R",2,2,1,1)-0.072063*ewgc("CuGR",2,2)+0.154056*ewgc("CuuR",0,2,2,0)
                    +0.002073*ewgc("CuuR",1,2,2,1)
                    -0.000189*ewgc("Cqd1R",2,2,0,0)+-0.000024*ewgc("Cqd1R",2,2,1,1)+0.011560*ewgc("Cqq1R",0,0,2,2)
                    +-0.000319*ewgc("Cqq1R",1,1,2,2)+0.039980*ewgc("Cqq3R",0,0,2,2)+0.001166*ewgc("Cqq3R",1,1,2,2)
                    +0.002430*ewgc("Cqu1R",0,0,2,2)+0.000093*ewgc("Cqu1R",1,1,2,2)+0.001730*ewgc("Cqu1R",2,2,0,0)
                    +0.000036*ewgc("Cqu1R",2,2,1,1)+-0.001642*ewgc("Cud1R",2,2,0,0)+-0.000091*ewgc("Cud1R",2,2,1,1)
                    +0.013720*ewgc("CuuR",0,0,2,2)+0.000218*ewgc("CuuR",1,1,2,2)
                    )*(SM_sigma_tt_bin_2000_2300/sigma_tt_bin_2000_2300_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
            
        }
        
    } else if(b_min == 2300 && b_max == 3500){
        
        double SM_sigma_tt_bin_2300_3500 = SM.getOptionalParameter("SM_sigma_tt_CMS_bin_mtt_2300_3500");
        double sigma_tt_bin_2300_3500_madgraph = 0.140478;
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tt_bin_2300_3500+(-0.016913*ewgc("CG")+0.025555*ewgc("Cqd8R",2,2,0,0)+0.001168*ewgc("Cqd8R",2,2,1,1)
                    +0.211627*ewgc("Cqq1R",0,2,2,0)+0.002920*ewgc("Cqq1R",1,2,2,1)+0.406226*ewgc("Cqq3R",0,2,2,0)
                    +0.011576*ewgc("Cqq3R",1,2,2,1)+0.075650*ewgc("Cqu8R",0,0,2,2)+0.001858*ewgc("Cqu8R",1,1,2,2)
                    +0.049952*ewgc("Cqu8R",2,2,0,0)+0.000687*ewgc("Cqu8R",2,2,1,1)+0.025594*ewgc("Cud8R",2,2,0,0)
                    +0.001165*ewgc("Cud8R",2,2,1,1)-0.054955*ewgc("CuGR",2,2)+0.205602*ewgc("CuuR",0,2,2,0)
                    +0.002840*ewgc("CuuR",1,2,2,1)
                    -0.000606*ewgc("Cqd1R",2,2,0,0)+-0.000021*ewgc("Cqd1R",2,2,1,1)+0.018206*ewgc("Cqq1R",0,0,2,2)
                    +-0.000313*ewgc("Cqq1R",1,1,2,2)+0.052844*ewgc("Cqq3R",0,0,2,2)+0.001307*ewgc("Cqq3R",1,1,2,2)
                    +0.003438*ewgc("Cqu1R",0,0,2,2)+0.000081*ewgc("Cqu1R",1,1,2,2)+0.002306*ewgc("Cqu1R",2,2,0,0)
                    +0.000029*ewgc("Cqu1R",2,2,1,1)+-0.002283*ewgc("Cud1R",2,2,0,0)+-0.000106*ewgc("Cud1R",2,2,1,1)
                    +0.017622*ewgc("CuuR",0,0,2,2)+0.000252*ewgc("CuuR",1,1,2,2)
                    )*(SM_sigma_tt_bin_2300_3500/sigma_tt_bin_2300_3500_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
            
            
        }
        
    } else {
        throw std::runtime_error("\nERROR: Please specify a correct binning range for sigma_tt_diff_LO.\n");
    }

}



//binning of ATLAS 1908.07305

sigma_norm_tt_diff_mtt_ATLAS_LO::sigma_norm_tt_diff_mtt_ATLAS_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{

    setParametersForObservable(make_vector<std::string>() << "SM_sigma_norm_tt_ATLAS_bin_mtt_325_400" << 
            "SM_sigma_norm_tt_ATLAS_bin_mtt_400_480" << "SM_sigma_norm_tt_ATLAS_bin_mtt_480_580" << 
            "SM_sigma_norm_tt_ATLAS_bin_mtt_580_700" << "SM_sigma_norm_tt_ATLAS_bin_mtt_700_860" << 
            "SM_sigma_norm_tt_ATLAS_bin_mtt_860_1020" << "SM_sigma_norm_tt_ATLAS_bin_mtt_1020_1250" <<
            "SM_sigma_norm_tt_ATLAS_bin_mtt_1250_1500" << "SM_sigma_norm_tt_ATLAS_bin_mtt_1500_2000" );
                
}

double sigma_norm_tt_diff_mtt_ATLAS_LO::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    
    
    double sigma_tt_inclusive_madgraph = 480.242000;
    double sigma_tt_inclusive_NP = (-32.061600*ewgc("CG")+-0.072252*ewgc("Cqd1R",2,2,0,0)+0.024384*ewgc("Cqd1R",2,2,1,1)
    +2.795340*ewgc("Cqd8R",2,2,0,0)+0.483828*ewgc("Cqd8R",2,2,1,1)+1.212250*ewgc("Cqq1R",0,0,2,2)+18.371200*ewgc("Cqq1R",0,2,2,0)+
    -0.080958*ewgc("Cqq1R",1,1,2,2)+0.753096*ewgc("Cqq1R",1,2,2,1)+4.019140*ewgc("Cqq3R",0,0,2,2)+39.347700*ewgc("Cqq3R",0,2,2,0)
    +0.353400*ewgc("Cqq3R",1,1,2,2)+4.100860*ewgc("Cqq3R",1,2,2,1)+0.431874*ewgc("Cqu1R",0,0,2,2)+0.039174*ewgc("Cqu1R",1,1,2,2)
    +0.388572*ewgc("Cqu1R",2,2,0,0)+0.061362*ewgc("Cqu1R",2,2,1,1)+7.158910*ewgc("Cqu8R",0,0,2,2)+0.650778*ewgc("Cqu8R",1,1,2,2)
    +4.412830*ewgc("Cqu8R",2,2,0,0)+0.216756*ewgc("Cqu8R",2,2,1,1)+-0.160380*ewgc("Cud1R",2,2,0,0)+0.016986*ewgc("Cud1R",2,2,1,1)
    +2.795980*ewgc("Cud8R",2,2,0,0)+0.487194*ewgc("Cud8R",2,2,1,1)+-140.557000*ewgc("CuGR",2,2)+1.352220*ewgc("CuuR",0,0,2,2)
    +17.951000*ewgc("CuuR",0,2,2,0)+0.095604*ewgc("CuuR",1,1,2,2)+0.739128*ewgc("CuuR",1,2,2,1));
    
    
    if(b_min == 325 && b_max == 400){
        
        double SM_sigma_norm_tt_bin_325_400 = SM.getOptionalParameter("SM_sigma_norm_tt_ATLAS_bin_mtt_325_400");
        double sigma_tt_bin_325_400_madgraph = 95.612100;
        
        
        double sigma_tt_bin_325_400_NP_lin = (-2.422920*ewgc("CG")+0.007578*ewgc("Cqd1R",2,2,0,0)+0.042222*ewgc("Cqd1R",2,2,1,1)+
        0.401694*ewgc("Cqd8R",2,2,0,0)+0.123768*ewgc("Cqd8R",2,2,1,1)+0.159660*ewgc("Cqq1R",0,0,2,2)+2.248510*ewgc("Cqq1R",0,2,2,0)
        +0.024942*ewgc("Cqq1R",1,1,2,2)+0.184272*ewgc("Cqq1R",1,2,2,1)+0.397038*ewgc("Cqq3R",0,0,2,2)+4.951130*ewgc("Cqq3R",0,2,2,0)
        +0.088098*ewgc("Cqq3R",1,1,2,2)+0.780684*ewgc("Cqq3R",1,2,2,1)+0.100374*ewgc("Cqu1R",0,0,2,2)+0.039066*ewgc("Cqu1R",1,1,2,2)
        +0.110028*ewgc("Cqu1R",2,2,0,0)+0.051414*ewgc("Cqu1R",2,2,1,1)+0.926298*ewgc("Cqu8R",0,0,2,2)+0.152874*ewgc("Cqu8R",1,1,2,2)
        +0.575544*ewgc("Cqu8R",2,2,0,0)+0.079620*ewgc("Cqu8R",2,2,1,1)+0.023658*ewgc("Cud1R",2,2,0,0)+0.044238*ewgc("Cud1R",2,2,1,1)
        +0.404088*ewgc("Cud8R",2,2,0,0)+0.127740*ewgc("Cud8R",2,2,1,1)+-31.370400*ewgc("CuGR",2,2)+0.146250*ewgc("CuuR",0,0,2,2)
        +2.216290*ewgc("CuuR",0,2,2,0)+0.059376*ewgc("CuuR",1,1,2,2)+0.180834*ewgc("CuuR",1,2,2,1));
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total_norm_tt_bin_325_400 = SM_sigma_norm_tt_bin_325_400 + SM_sigma_norm_tt_bin_325_400*(
            sigma_tt_bin_325_400_NP_lin/sigma_tt_bin_325_400_madgraph - sigma_tt_inclusive_NP/sigma_tt_inclusive_madgraph);
            
            if (total_norm_tt_bin_325_400 < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total_norm_tt_bin_325_400;
            
        }
        
    } else if(b_min == 400 && b_max == 480){
        
        double SM_sigma_norm_tt_bin_400_480 = SM.getOptionalParameter("SM_sigma_norm_tt_ATLAS_bin_mtt_400_480");
        double sigma_tt_bin_400_480_madgraph = 147.218000;
        
        
        double sigma_tt_bin_400_480_NP_lin = (+-8.063800*ewgc("CG")+-0.030360*ewgc("Cqd1R",2,2,0,0)+-0.014292*ewgc("Cqd1R",2,2,1,1)
        +0.537774*ewgc("Cqd8R",2,2,0,0)+0.102036*ewgc("Cqd8R",2,2,1,1)+0.202368*ewgc("Cqq1R",0,0,2,2)+3.659260*ewgc("Cqq1R",0,2,2,0)
        +-0.028410*ewgc("Cqq1R",1,1,2,2)+0.168744*ewgc("Cqq1R",1,2,2,1)+0.682230*ewgc("Cqq3R",0,0,2,2)+7.940650*ewgc("Cqq3R",0,2,2,0)
        +0.060234*ewgc("Cqq3R",1,1,2,2)+1.036730*ewgc("Cqq3R",1,2,2,1)+0.056178*ewgc("Cqu1R",0,0,2,2)+-0.007908*ewgc("Cqu1R",1,1,2,2)
        +0.073824*ewgc("Cqu1R",2,2,0,0)+-0.001914*ewgc("Cqu1R",2,2,1,1)+1.420910*ewgc("Cqu8R",0,0,2,2)+0.147192*ewgc("Cqu8R",1,1,2,2)
        +0.859650*ewgc("Cqu8R",2,2,0,0)+0.036930*ewgc("Cqu8R",2,2,1,1)+-0.048564*ewgc("Cud1R",2,2,0,0)+-0.013968*ewgc("Cud1R",2,2,1,1)
        +0.537786*ewgc("Cud8R",2,2,0,0)+0.101004*ewgc("Cud8R",2,2,1,1)+-42.879900*ewgc("CuGR",2,2)+0.267738*ewgc("CuuR",0,0,2,2));
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total_norm_tt_bin_400_480 = SM_sigma_norm_tt_bin_400_480 + SM_sigma_norm_tt_bin_400_480*(
            sigma_tt_bin_400_480_NP_lin/sigma_tt_bin_400_480_madgraph - sigma_tt_inclusive_NP/sigma_tt_inclusive_madgraph);

            if (total_norm_tt_bin_400_480 < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total_norm_tt_bin_400_480;
            
        }
        
        
    } else if(b_min == 480 && b_max == 580){
        
        double SM_sigma_norm_tt_bin_480_580 = SM.getOptionalParameter("SM_sigma_norm_tt_ATLAS_bin_mtt_480_580");
        double sigma_tt_bin_480_580_madgraph = 110.230000;
        
        
        double sigma_tt_bin_480_580_NP_lin = (+-8.469500*ewgc("CG")+-0.009222*ewgc("Cqd1R",2,2,0,0)+-0.002448*ewgc("Cqd1R",2,2,1,1)
        +0.534144*ewgc("Cqd8R",2,2,0,0)+0.096336*ewgc("Cqd8R",2,2,1,1)+0.204966*ewgc("Cqq1R",0,0,2,2)+3.353950*ewgc("Cqq1R",0,2,2,0)
        +-0.024336*ewgc("Cqq1R",1,1,2,2)+0.148722*ewgc("Cqq1R",1,2,2,1)+0.719910*ewgc("Cqq3R",0,0,2,2)+7.301870*ewgc("Cqq3R",0,2,2,0)
        +0.071550*ewgc("Cqq3R",1,1,2,2)+0.848766*ewgc("Cqq3R",1,2,2,1)+0.100416*ewgc("Cqu1R",0,0,2,2)+0.002706*ewgc("Cqu1R",1,1,2,2)
        +0.071370*ewgc("Cqu1R",2,2,0,0)+0.005880*ewgc("Cqu1R",2,2,1,1)+1.334360*ewgc("Cqu8R",0,0,2,2)+0.130434*ewgc("Cqu8R",1,1,2,2)
        +0.813354*ewgc("Cqu8R",2,2,0,0)+0.035946*ewgc("Cqu8R",2,2,1,1)+-0.027816*ewgc("Cud1R",2,2,0,0)+-0.003090*ewgc("Cud1R",2,2,1,1)
        +0.529638*ewgc("Cud8R",2,2,0,0)+0.097110*ewgc("Cud8R",2,2,1,1)+-30.595400*ewgc("CuGR",2,2)+0.224958*ewgc("CuuR",0,0,2,2)
        +3.281180*ewgc("CuuR",0,2,2,0)+0.011376*ewgc("CuuR",1,1,2,2)+0.146868*ewgc("CuuR",1,2,2,1));
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{

            double total_norm_tt_bin_480_580 = SM_sigma_norm_tt_bin_480_580 + SM_sigma_norm_tt_bin_480_580*(
            sigma_tt_bin_480_580_NP_lin/sigma_tt_bin_480_580_madgraph - sigma_tt_inclusive_NP/sigma_tt_inclusive_madgraph);

            
            if (total_norm_tt_bin_480_580 < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total_norm_tt_bin_480_580;
            
        }
        
        
    } else if(b_min == 580 && b_max == 700){
        
        double SM_sigma_norm_tt_bin_580_700 = SM.getOptionalParameter("SM_sigma_norm_tt_ATLAS_bin_mtt_580_700");
        double sigma_tt_bin_580_700_madgraph = 64.291200;
        
        
        double sigma_tt_bin_580_700_NP_lin = (-6.041880*ewgc("CG")+-0.020183*ewgc("Cqd1R",2,2,0,0)+-0.001903*ewgc("Cqd1R",2,2,1,1)
        +0.402131*ewgc("Cqd8R",2,2,0,0)+0.062093*ewgc("Cqd8R",2,2,1,1)+0.184619*ewgc("Cqq1R",0,0,2,2)+2.713140*ewgc("Cqq1R",0,2,2,0)
        +-0.021160*ewgc("Cqq1R",1,1,2,2)+0.097214*ewgc("Cqq1R",1,2,2,1)+0.621278*ewgc("Cqq3R",0,0,2,2)+5.805340*ewgc("Cqq3R",0,2,2,0)
        +0.045730*ewgc("Cqq3R",1,1,2,2)+0.564125*ewgc("Cqq3R",1,2,2,1)+0.050442*ewgc("Cqu1R",0,0,2,2)+-0.000803*ewgc("Cqu1R",1,1,2,2)
        +0.043475*ewgc("Cqu1R",2,2,0,0)+0.002446*ewgc("Cqu1R",2,2,1,1)+1.042640*ewgc("Cqu8R",0,0,2,2)+0.085265*ewgc("Cqu8R",1,1,2,2)
        +0.647214*ewgc("Cqu8R",2,2,0,0)+0.025471*ewgc("Cqu8R",2,2,1,1)+-0.037187*ewgc("Cud1R",2,2,0,0)+-0.005738*ewgc("Cud1R",2,2,1,1)
        +0.405243*ewgc("Cud8R",2,2,0,0)+0.062329*ewgc("Cud8R",2,2,1,1)+-17.737100*ewgc("CuGR",2,2)+0.206038*ewgc("CuuR",0,0,2,2)
        +2.643790*ewgc("CuuR",0,2,2,0)+0.007879*ewgc("CuuR",1,1,2,2)+0.093899*ewgc("CuuR",1,2,2,1));
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total_norm_tt_bin_580_700 = SM_sigma_norm_tt_bin_580_700 + SM_sigma_norm_tt_bin_580_700*(
            sigma_tt_bin_580_700_NP_lin/sigma_tt_bin_580_700_madgraph - sigma_tt_inclusive_NP/sigma_tt_inclusive_madgraph);
            
            if (total_norm_tt_bin_580_700 < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total_norm_tt_bin_580_700;
            
        }
        
        
    } else if(b_min == 700 && b_max == 860){
        
        double SM_sigma_norm_tt_bin_700_860 = SM.getOptionalParameter("SM_sigma_norm_tt_ATLAS_bin_mtt_700_860");
        double sigma_tt_bin_700_860_madgraph = 35.990600;
        
        
        double sigma_tt_bin_700_860_NP_lin = (+-3.854860*ewgc("CG")+-0.006342*ewgc("Cqd1R",2,2,0,0)+0.000815*ewgc("Cqd1R",2,2,1,1)
        +0.330269*ewgc("Cqd8R",2,2,0,0)+0.044940*ewgc("Cqd8R",2,2,1,1)+0.146891*ewgc("Cqq1R",0,0,2,2)+2.230090*ewgc("Cqq1R",0,2,2,0)
        +-0.014938*ewgc("Cqq1R",1,1,2,2)+0.068073*ewgc("Cqq1R",1,2,2,1)+0.533029*ewgc("Cqq3R",0,0,2,2)+4.739040*ewgc("Cqq3R",0,2,2,0)
        +0.037619*ewgc("Cqq3R",1,1,2,2)+0.391305*ewgc("Cqq3R",1,2,2,1)+0.048161*ewgc("Cqu1R",0,0,2,2)+0.002162*ewgc("Cqu1R",1,1,2,2)
        +0.035290*ewgc("Cqu1R",2,2,0,0)+0.001030*ewgc("Cqu1R",2,2,1,1)+0.861595*ewgc("Cqu8R",0,0,2,2)+0.061229*ewgc("Cqu8R",1,1,2,2)
        +0.532168*ewgc("Cqu8R",2,2,0,0)+0.016974*ewgc("Cqu8R",2,2,1,1)+-0.021213*ewgc("Cud1R",2,2,0,0)+-0.001135*ewgc("Cud1R",2,2,1,1)
        +0.331528*ewgc("Cud8R",2,2,0,0)+0.044438*ewgc("Cud8R",2,2,1,1)+-10.060700*ewgc("CuGR",2,2)+0.167516*ewgc("CuuR",0,0,2,2)
        +2.172550*ewgc("CuuR",0,2,2,0)+0.007310*ewgc("CuuR",1,1,2,2)+0.067370*ewgc("CuuR",1,2,2,1)
        );
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
                        
            double total_norm_tt_bin_700_860 = SM_sigma_norm_tt_bin_700_860 + SM_sigma_norm_tt_bin_700_860*(
            sigma_tt_bin_700_860_NP_lin/sigma_tt_bin_700_860_madgraph - sigma_tt_inclusive_NP/sigma_tt_inclusive_madgraph);

            
            
            if (total_norm_tt_bin_700_860 < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total_norm_tt_bin_700_860;
            
        }
        
        
    } else if(b_min == 860 && b_max == 1020){
        
        double SM_sigma_norm_tt_bin_860_1020 = SM.getOptionalParameter("SM_sigma_norm_tt_ATLAS_bin_mtt_860_1020");
        double sigma_tt_bin_860_1020_madgraph = 14.350200;
        
        
        double sigma_tt_bin_860_1020_NP_lin = (+-1.661030*ewgc("CG")+-0.004506*ewgc("Cqd1R",2,2,0,0)+-0.000767*ewgc("Cqd1R",2,2,1,1)
        +0.197992*ewgc("Cqd8R",2,2,0,0)+0.022552*ewgc("Cqd8R",2,2,1,1)+0.098024*ewgc("Cqq1R",0,0,2,2)+1.356740*ewgc("Cqq1R",0,2,2,0)
        +-0.007518*ewgc("Cqq1R",1,1,2,2)+0.033961*ewgc("Cqq1R",1,2,2,1)+0.339810*ewgc("Cqq3R",0,0,2,2)+2.851490*ewgc("Cqq3R",0,2,2,0)
        +0.020041*ewgc("Cqq3R",1,1,2,2)+0.198170*ewgc("Cqq3R",1,2,2,1)+0.022757*ewgc("Cqu1R",0,0,2,2)+0.001121*ewgc("Cqu1R",1,1,2,2)
        +0.017951*ewgc("Cqu1R",2,2,0,0)+0.000351*ewgc("Cqu1R",2,2,1,1)+0.517550*ewgc("Cqu8R",0,0,2,2)+0.029952*ewgc("Cqu8R",1,1,2,2)
        +0.320654*ewgc("Cqu8R",2,2,0,0)+0.008393*ewgc("Cqu8R",2,2,1,1)+-0.016088*ewgc("Cud1R",2,2,0,0)+-0.002105*ewgc("Cud1R",2,2,1,1)
        +0.196345*ewgc("Cud8R",2,2,0,0)+0.022422*ewgc("Cud8R",2,2,1,1)+-4.106080*ewgc("CuGR",2,2)+0.108949*ewgc("CuuR",0,0,2,2)
        +1.321780*ewgc("CuuR",0,2,2,0)+0.002587*ewgc("CuuR",1,1,2,2)+0.032883*ewgc("CuuR",1,2,2,1)
        );
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            
            
            double total_norm_tt_bin_860_1020 = SM_sigma_norm_tt_bin_860_1020 + SM_sigma_norm_tt_bin_860_1020*(
            sigma_tt_bin_860_1020_NP_lin/sigma_tt_bin_860_1020_madgraph - sigma_tt_inclusive_NP/sigma_tt_inclusive_madgraph);

            
            
            if (total_norm_tt_bin_860_1020 < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total_norm_tt_bin_860_1020;
            
        }
        
        
    } else if(b_min == 1020 && b_max == 1250){
        
        double SM_sigma_norm_tt_bin_1020_1250 = SM.getOptionalParameter("SM_sigma_norm_tt_ATLAS_bin_mtt_1020_1250");
        double sigma_tt_bin_1020_1250_madgraph = 7.837970;
        
        
        double sigma_tt_bin_1020_1250_NP_lin = (+-0.952352*ewgc("CG")+-0.004072*ewgc("Cqd1R",2,2,0,0)+0.000929*ewgc("Cqd1R",2,2,1,1)
        +0.164245*ewgc("Cqd8R",2,2,0,0)+0.016778*ewgc("Cqd8R",2,2,1,1)+0.082114*ewgc("Cqq1R",0,0,2,2)+1.139650*ewgc("Cqq1R",0,2,2,0)
        +-0.004040*ewgc("Cqq1R",1,1,2,2)+0.025379*ewgc("Cqq1R",1,2,2,1)+0.292593*ewgc("Cqq3R",0,0,2,2)+2.374990*ewgc("Cqq3R",0,2,2,0)
        +0.015213*ewgc("Cqq3R",1,1,2,2)+0.140499*ewgc("Cqq3R",1,2,2,1)+0.024363*ewgc("Cqu1R",0,0,2,2)+0.001848*ewgc("Cqu1R",1,1,2,2)
        +0.016591*ewgc("Cqu1R",2,2,0,0)+0.001685*ewgc("Cqu1R",2,2,1,1)+0.434440*ewgc("Cqu8R",0,0,2,2)+0.022335*ewgc("Cqu8R",1,1,2,2)
        +0.269693*ewgc("Cqu8R",2,2,0,0)+0.007106*ewgc("Cqu8R",2,2,1,1)+-0.014099*ewgc("Cud1R",2,2,0,0)+-0.000052*ewgc("Cud1R",2,2,1,1)
        +0.164420*ewgc("Cud8R",2,2,0,0)+0.016868*ewgc("Cud8R",2,2,1,1)+-2.312160*ewgc("CuGR",2,2)+0.091177*ewgc("CuuR",0,0,2,2)
        +1.107260*ewgc("CuuR",0,2,2,0)+0.003389*ewgc("CuuR",1,1,2,2)+0.024817*ewgc("CuuR",1,2,2,1)
        );
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            
            
            double total_norm_tt_bin_1020_1250 = SM_sigma_norm_tt_bin_1020_1250 + SM_sigma_norm_tt_bin_1020_1250*(
            sigma_tt_bin_1020_1250_NP_lin/sigma_tt_bin_1020_1250_madgraph - sigma_tt_inclusive_NP/sigma_tt_inclusive_madgraph);

            
            
            if (total_norm_tt_bin_1020_1250 < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total_norm_tt_bin_1020_1250;
            
        }
        
        
    } else if(b_min == 1250 && b_max == 1500){
        
        double SM_sigma_norm_tt_bin_1250_1500 = SM.getOptionalParameter("SM_sigma_norm_tt_ATLAS_bin_mtt_1250_1500");
        double sigma_tt_bin_1250_1500_madgraph = 2.891000;
        
        
        double sigma_tt_bin_1250_1500_NP_lin = (+-0.363269*ewgc("CG")+-0.002056*ewgc("Cqd1R",2,2,0,0)+-0.000125*ewgc("Cqd1R",2,2,1,1)
        +0.096439*ewgc("Cqd8R",2,2,0,0)+0.007646*ewgc("Cqd8R",2,2,1,1)+0.054967*ewgc("Cqq1R",0,0,2,2)+0.685920*ewgc("Cqq1R",0,2,2,0)
        +-0.002843*ewgc("Cqq1R",1,1,2,2)+0.012355*ewgc("Cqq1R",1,2,2,1)+0.181075*ewgc("Cqq3R",0,0,2,2)+1.415120*ewgc("Cqq3R",0,2,2,0)
        +0.007077*ewgc("Cqq3R",1,1,2,2)+0.068674*ewgc("Cqq3R",1,2,2,1)+0.011300*ewgc("Cqu1R",0,0,2,2)+0.000325*ewgc("Cqu1R",1,1,2,2)
        +0.008365*ewgc("Cqu1R",2,2,0,0)+0.000111*ewgc("Cqu1R",2,2,1,1)+0.258096*ewgc("Cqu8R",0,0,2,2)+0.010409*ewgc("Cqu8R",1,1,2,2)
        +0.161941*ewgc("Cqu8R",2,2,0,0)+0.002867*ewgc("Cqu8R",2,2,1,1)+-0.007602*ewgc("Cud1R",2,2,0,0)+-0.000615*ewgc("Cud1R",2,2,1,1)
        +0.096217*ewgc("Cud8R",2,2,0,0)+0.007558*ewgc("Cud8R",2,2,1,1)+-0.887639*ewgc("CuGR",2,2)+0.058000*ewgc("CuuR",0,0,2,2)
        +0.666640*ewgc("CuuR",0,2,2,0)+0.001038*ewgc("CuuR",1,1,2,2)+0.011866*ewgc("CuuR",1,2,2,1)
        );
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{

            
            
            double total_norm_tt_bin_1250_1500 = SM_sigma_norm_tt_bin_1250_1500 + SM_sigma_norm_tt_bin_1250_1500*(
            sigma_tt_bin_1250_1500_NP_lin/sigma_tt_bin_1250_1500_madgraph - sigma_tt_inclusive_NP/sigma_tt_inclusive_madgraph);

            
            
            if (total_norm_tt_bin_1250_1500 < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total_norm_tt_bin_1250_1500;
            
        }
        
        
    } else if(b_min == 1500 && b_max == 2000){
        
        double SM_sigma_norm_tt_bin_1500_2000 = SM.getOptionalParameter("SM_sigma_norm_tt_ATLAS_bin_mtt_1500_2000");
        double sigma_tt_bin_1500_2000_madgraph = 1.821290;
        
        
        double sigma_tt_bin_1500_2000_NP_lin = (+-0.231926*ewgc("CG")+-0.003089*ewgc("Cqd1R",2,2,0,0)+-0.000047*ewgc("Cqd1R",2,2,1,1)
        +0.130652*ewgc("Cqd8R",2,2,0,0)+0.007679*ewgc("Cqd8R",2,2,1,1)+0.078637*ewgc("Cqq1R",0,0,2,2)+0.983935*ewgc("Cqq1R",0,2,2,0)
        +-0.002654*ewgc("Cqq1R",1,1,2,2)+0.014376*ewgc("Cqq1R",1,2,2,1)+0.252179*ewgc("Cqq3R",0,0,2,2)+1.968040*ewgc("Cqq3R",0,2,2,0)
        +0.007838*ewgc("Cqq3R",1,1,2,2)+0.071898*ewgc("Cqq3R",1,2,2,1)+0.017883*ewgc("Cqu1R",0,0,2,2)+0.000657*ewgc("Cqu1R",1,1,2,2)
        +0.011678*ewgc("Cqu1R",2,2,0,0)+0.000359*ewgc("Cqu1R",2,2,1,1)+0.363022*ewgc("Cqu8R",0,0,2,2)+0.011088*ewgc("Cqu8R",1,1,2,2)
        +0.232613*ewgc("Cqu8R",2,2,0,0)+0.003450*ewgc("Cqu8R",2,2,1,1)+-0.011469*ewgc("Cud1R",2,2,0,0)+-0.000549*ewgc("Cud1R",2,2,1,1)
        +0.130711*ewgc("Cud8R",2,2,0,0)+0.007725*ewgc("Cud8R",2,2,1,1)+-0.607749*ewgc("CuGR",2,2)+0.081595*ewgc("CuuR",0,0,2,2)
        +0.955959*ewgc("CuuR",0,2,2,0)+0.001335*ewgc("CuuR",1,1,2,2)+0.014025*ewgc("CuuR",1,2,2,1)
        );
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            
            double total_norm_tt_bin_1500_2000 = SM_sigma_norm_tt_bin_1500_2000 + SM_sigma_norm_tt_bin_1500_2000*(
            sigma_tt_bin_1500_2000_NP_lin/sigma_tt_bin_1500_2000_madgraph - sigma_tt_inclusive_NP/sigma_tt_inclusive_madgraph);

            
            if (total_norm_tt_bin_1500_2000 < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total_norm_tt_bin_1500_2000;
            
        }
        
        
    }
}






//// ttbar inclusive 13 TeV ////

sigma_tt_13_LO::sigma_tt_13_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tt_13");
}

double sigma_tt_13_LO::computeThValue()
{
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    
    double SM_sigma_tt_13 = SM.getOptionalParameter("SM_sigma_tt_13");
    double sigma_tt_13_madgraph = 480.242000;
    double sigma_tt_13_NP = (-32.061600*ewgc("CG")+-0.072252*ewgc("Cqd1R",2,2,0,0)+0.024384*ewgc("Cqd1R",2,2,1,1)
    +2.795340*ewgc("Cqd8R",2,2,0,0)+0.483828*ewgc("Cqd8R",2,2,1,1)+1.212250*ewgc("Cqq1R",0,0,2,2)+18.371200*ewgc("Cqq1R",0,2,2,0)+
    -0.080958*ewgc("Cqq1R",1,1,2,2)+0.753096*ewgc("Cqq1R",1,2,2,1)+4.019140*ewgc("Cqq3R",0,0,2,2)+39.347700*ewgc("Cqq3R",0,2,2,0)
    +0.353400*ewgc("Cqq3R",1,1,2,2)+4.100860*ewgc("Cqq3R",1,2,2,1)+0.431874*ewgc("Cqu1R",0,0,2,2)+0.039174*ewgc("Cqu1R",1,1,2,2)
    +0.388572*ewgc("Cqu1R",2,2,0,0)+0.061362*ewgc("Cqu1R",2,2,1,1)+7.158910*ewgc("Cqu8R",0,0,2,2)+0.650778*ewgc("Cqu8R",1,1,2,2)
    +4.412830*ewgc("Cqu8R",2,2,0,0)+0.216756*ewgc("Cqu8R",2,2,1,1)+-0.160380*ewgc("Cud1R",2,2,0,0)+0.016986*ewgc("Cud1R",2,2,1,1)
    +2.795980*ewgc("Cud8R",2,2,0,0)+0.487194*ewgc("Cud8R",2,2,1,1)+-140.557000*ewgc("CuGR",2,2)+1.352220*ewgc("CuuR",0,0,2,2)
    +17.951000*ewgc("CuuR",0,2,2,0)+0.095604*ewgc("CuuR",1,1,2,2)+0.739128*ewgc("CuuR",1,2,2,1)
    )*(SM_sigma_tt_13/sigma_tt_13_madgraph);
    
    
    if(flag_Quadratic){
        
        return  0.;
        
        }
    
        else{
        
            double total_sigma_tt_13_LO = SM_sigma_tt_13 + sigma_tt_13_NP;

            if (total_sigma_tt_13_LO < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total_sigma_tt_13_LO;
            
        }

}




//// ttbar inclusive 8 TeV over 7 TeV prediction ////

R_tt_8_o_7_LO::R_tt_8_o_7_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_R_tt_8_o_7");
}

double R_tt_8_o_7_LO::computeThValue()
{
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    
    double SM_R_tt_8_o_7 = SM.getOptionalParameter("SM_R_tt_8_o_7");
    double sigma_tt_7_madgraph = 103.104000;
    double sigma_tt_8_madgraph = 146.968000;
    
    double sigma_tt_7_NP_lin = +-5.648030*ewgc("CG")+-0.032856*ewgc("Cqd1R",2,2,0,0)+-0.001782*ewgc("Cqd1R",2,2,1,1)
    +0.878394*ewgc("Cqd8R",2,2,0,0)+0.094212*ewgc("Cqd8R",2,2,1,1)+0.401058*ewgc("Cqq1R",0,0,2,2)+6.068220*ewgc("Cqq1R",0,2,2,0)
    +-0.021348*ewgc("Cqq1R",1,1,2,2)+0.145452*ewgc("Cqq1R",1,2,2,1)+1.251470*ewgc("Cqq3R",0,0,2,2)+12.803700*ewgc("Cqq3R",0,2,2,0)
    +0.064734*ewgc("Cqq3R",1,1,2,2)+0.848508*ewgc("Cqq3R",1,2,2,1)+0.138378*ewgc("Cqu1R",0,0,2,2)+0.001884*ewgc("Cqu1R",1,1,2,2)
    +0.118542*ewgc("Cqu1R",2,2,0,0)+0.006222*ewgc("Cqu1R",2,2,1,1)+2.324460*ewgc("Cqu8R",0,0,2,2)+0.130314*ewgc("Cqu8R",1,1,2,2)
    +1.450870*ewgc("Cqu8R",2,2,0,0)+0.037656*ewgc("Cqu8R",2,2,1,1)+-0.060552*ewgc("Cud1R",2,2,0,0)+-0.004116*ewgc("Cud1R",2,2,1,1)
    +0.880116*ewgc("Cud8R",2,2,0,0)+0.094002*ewgc("Cud8R",2,2,1,1)+-29.721300*ewgc("CuGR",2,2)+0.418404*ewgc("CuuR",0,0,2,2)
    +5.935330*ewgc("CuuR",0,2,2,0)+0.014052*ewgc("CuuR",1,1,2,2)+0.141732*ewgc("CuuR",1,2,2,1);
    
    
    double sigma_tt_8_NP_lin = -8.487130*ewgc("CG")+-0.046260*ewgc("Cqd1R",2,2,0,0)+-0.005154*ewgc("Cqd1R",2,2,1,1)
    +1.147210*ewgc("Cqd8R",2,2,0,0)+0.132288*ewgc("Cqd8R",2,2,1,1)+0.516528*ewgc("Cqq1R",0,0,2,2)+7.896180*ewgc("Cqq1R",0,2,2,0)
    +-0.036054*ewgc("Cqq1R",1,1,2,2)+0.207276*ewgc("Cqq1R",1,2,2,1)+1.648780*ewgc("Cqq3R",0,0,2,2)+16.705700*ewgc("Cqq3R",0,2,2,0)
    +0.089058*ewgc("Cqq3R",1,1,2,2)+1.222520*ewgc("Cqq3R",1,2,2,1)+0.170862*ewgc("Cqu1R",0,0,2,2)+-0.003144*ewgc("Cqu1R",1,1,2,2)
    +0.146586*ewgc("Cqu1R",2,2,0,0)+0.005466*ewgc("Cqu1R",2,2,1,1)+3.028300*ewgc("Cqu8R",0,0,2,2)+0.181704*ewgc("Cqu8R",1,1,2,2)
    +1.881010*ewgc("Cqu8R",2,2,0,0)+0.049662*ewgc("Cqu8R",2,2,1,1)+-0.086898*ewgc("Cud1R",2,2,0,0)+-0.009084*ewgc("Cud1R",2,2,1,1)
    +1.148470*ewgc("Cud8R",2,2,0,0)+0.132612*ewgc("Cud8R",2,2,1,1)+-42.562100*ewgc("CuGR",2,2)+0.560760*ewgc("CuuR",0,0,2,2)
    +7.717930*ewgc("CuuR",0,2,2,0)+0.015432*ewgc("CuuR",1,1,2,2)+0.202656*ewgc("CuuR",1,2,2,1);
    
    if(flag_Quadratic){
        
        return  0.;
        
        }
    
        else{

            double total_R_tt_8_o_7_LO = SM_R_tt_8_o_7 + SM_R_tt_8_o_7*(
            sigma_tt_8_NP_lin/sigma_tt_8_madgraph - sigma_tt_7_NP_lin/sigma_tt_7_madgraph);

            
            if (total_R_tt_8_o_7_LO < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total_R_tt_8_o_7_LO;
            
        }

}



//// ttbar inclusive 13 TeV over 8 TeV ////

R_tt_13_o_8_LO::R_tt_13_o_8_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_R_tt_13_o_8");
}

double R_tt_13_o_8_LO::computeThValue()
{
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    double SM_R_tt_13_o_8 = SM.getOptionalParameter("SM_R_tt_13_o_8");
    double sigma_tt_8_madgraph = 146.968000;
    double sigma_tt_13_madgraph = 480.242000;
    
    double sigma_tt_8_NP_lin = -8.487130*ewgc("CG")+-0.046260*ewgc("Cqd1R",2,2,0,0)+-0.005154*ewgc("Cqd1R",2,2,1,1)
    +1.147210*ewgc("Cqd8R",2,2,0,0)+0.132288*ewgc("Cqd8R",2,2,1,1)+0.516528*ewgc("Cqq1R",0,0,2,2)+7.896180*ewgc("Cqq1R",0,2,2,0)
    +-0.036054*ewgc("Cqq1R",1,1,2,2)+0.207276*ewgc("Cqq1R",1,2,2,1)+1.648780*ewgc("Cqq3R",0,0,2,2)+16.705700*ewgc("Cqq3R",0,2,2,0)
    +0.089058*ewgc("Cqq3R",1,1,2,2)+1.222520*ewgc("Cqq3R",1,2,2,1)+0.170862*ewgc("Cqu1R",0,0,2,2)+-0.003144*ewgc("Cqu1R",1,1,2,2)
    +0.146586*ewgc("Cqu1R",2,2,0,0)+0.005466*ewgc("Cqu1R",2,2,1,1)+3.028300*ewgc("Cqu8R",0,0,2,2)+0.181704*ewgc("Cqu8R",1,1,2,2)
    +1.881010*ewgc("Cqu8R",2,2,0,0)+0.049662*ewgc("Cqu8R",2,2,1,1)+-0.086898*ewgc("Cud1R",2,2,0,0)+-0.009084*ewgc("Cud1R",2,2,1,1)
    +1.148470*ewgc("Cud8R",2,2,0,0)+0.132612*ewgc("Cud8R",2,2,1,1)+-42.562100*ewgc("CuGR",2,2)+0.560760*ewgc("CuuR",0,0,2,2)
    +7.717930*ewgc("CuuR",0,2,2,0)+0.015432*ewgc("CuuR",1,1,2,2)+0.202656*ewgc("CuuR",1,2,2,1);
    
    double sigma_tt_13_NP_lin = -32.061600*ewgc("CG")+-0.072252*ewgc("Cqd1R",2,2,0,0)+0.024384*ewgc("Cqd1R",2,2,1,1)
    +2.795340*ewgc("Cqd8R",2,2,0,0)+0.483828*ewgc("Cqd8R",2,2,1,1)+1.212250*ewgc("Cqq1R",0,0,2,2)+18.371200*ewgc("Cqq1R",0,2,2,0)+
    -0.080958*ewgc("Cqq1R",1,1,2,2)+0.753096*ewgc("Cqq1R",1,2,2,1)+4.019140*ewgc("Cqq3R",0,0,2,2)+39.347700*ewgc("Cqq3R",0,2,2,0)
    +0.353400*ewgc("Cqq3R",1,1,2,2)+4.100860*ewgc("Cqq3R",1,2,2,1)+0.431874*ewgc("Cqu1R",0,0,2,2)+0.039174*ewgc("Cqu1R",1,1,2,2)
    +0.388572*ewgc("Cqu1R",2,2,0,0)+0.061362*ewgc("Cqu1R",2,2,1,1)+7.158910*ewgc("Cqu8R",0,0,2,2)+0.650778*ewgc("Cqu8R",1,1,2,2)
    +4.412830*ewgc("Cqu8R",2,2,0,0)+0.216756*ewgc("Cqu8R",2,2,1,1)+-0.160380*ewgc("Cud1R",2,2,0,0)+0.016986*ewgc("Cud1R",2,2,1,1)
    +2.795980*ewgc("Cud8R",2,2,0,0)+0.487194*ewgc("Cud8R",2,2,1,1)+-140.557000*ewgc("CuGR",2,2)+1.352220*ewgc("CuuR",0,0,2,2)
    +17.951000*ewgc("CuuR",0,2,2,0)+0.095604*ewgc("CuuR",1,1,2,2)+0.739128*ewgc("CuuR",1,2,2,1);
    
    
    if(flag_Quadratic){
        
        return  0.;
        
        }
    
        else{

            double total_R_tt_13_o_8_LO = SM_R_tt_13_o_8 + SM_R_tt_13_o_8*(
            sigma_tt_13_NP_lin/sigma_tt_13_madgraph - sigma_tt_8_NP_lin/sigma_tt_8_madgraph);

            
            if (total_R_tt_13_o_8_LO < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total_R_tt_13_o_8_LO;
            
        }

}






////// ttbar charge asymmetry differential ///////////////////////////////////

charge_asymmetry_tt_diff_mtt_LO::charge_asymmetry_tt_diff_mtt_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{

    setParametersForObservable(make_vector<std::string>() << "SM_charge_asymmetry_bin_mtt_0_500" << "SM_charge_asymmetry_deno_bin_mtt_0_500"
            << "SM_charge_asymmetry_bin_mtt_500_750" << "SM_charge_asymmetry_deno_bin_mtt_500_750" << "SM_charge_asymmetry_bin_mtt_750_1000" 
            << "SM_charge_asymmetry_deno_bin_mtt_750_1000" << "SM_charge_asymmetry_bin_mtt_1000_1500" << "SM_charge_asymmetry_deno_bin_mtt_1000_1500"
            << "SM_charge_asymmetry_bin_mtt_1500_3000" << "SM_charge_asymmetry_deno_bin_mtt_1500_3000");
            //"SM_sigma_tt_bin_2300_2600" << "SM_sigma_tt_bin_2600_3000" << "SM_sigma_tt_bin_3000_3500" << "SM_sigma_tt_bin_3500_4000");
    
}


double charge_asymmetry_tt_diff_mtt_LO::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    


    
    
    
    if(b_min == 0 && b_max == 500){
        
        double SM_charge_asymmetry_bin_mtt_0_500 = SM.getOptionalParameter("SM_charge_asymmetry_bin_mtt_0_500");
        double SM_charge_asymmetry_deno_bin_mtt_0_500 = SM.getOptionalParameter("SM_charge_asymmetry_deno_bin_mtt_0_500");
        double SM_charge_asymmetry_num_bin_mtt_0_500 = SM_charge_asymmetry_bin_mtt_0_500*SM_charge_asymmetry_deno_bin_mtt_0_500;
        
        double SM_sigma_pos_bin_mtt_0_500 =0.5*(SM_charge_asymmetry_num_bin_mtt_0_500+SM_charge_asymmetry_deno_bin_mtt_0_500);
        double SM_sigma_neg_bin_mtt_0_500 =0.5*(SM_charge_asymmetry_deno_bin_mtt_0_500-SM_charge_asymmetry_num_bin_mtt_0_500);

        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            
                    
            double sigma_pos_bin_mtt_0_500_madgraph = 135.533000;
            double sigma_pos_bin_mtt_0_500_NP = -6.478050*ewgc("CG")+0.364812*ewgc("Cqd8R",2,2,0,0)
            +0.065870*ewgc("Cqd8R",2,2,1,1)+4.159410*ewgc("Cqq1R",0,2,2,0)+0.184853*ewgc("Cqq1R",1,2,2,1)
            +8.461010*ewgc("Cqq3R",0,2,2,0)+0.893405*ewgc("Cqq3R",1,2,2,1)+1.126020*ewgc("Cqu8R",0,0,2,2)
            -0.010419*ewgc("Cqu8R",1,1,2,2)+0.588354*ewgc("Cqu8R",2,2,0,0)+-0.002524*ewgc("Cqu8R",2,2,1,1)
            +0.455566*ewgc("Cud8R",2,2,0,0)+0.131078*ewgc("Cud8R",2,2,1,1)-39.064300*ewgc("CuGR",2,2)
            +4.072180*ewgc("CuuR",0,2,2,0)+0.097783*ewgc("CuuR",1,2,2,1)
            +0.022580*ewgc("Cqd1R",2,2,0,0)+0.004491*ewgc("Cqd1R",2,2,1,1)+0.359419*ewgc("Cqq1R",0,0,2,2)
            +-0.028395*ewgc("Cqq1R",1,1,2,2)+0.861834*ewgc("Cqq3R",0,0,2,2)+0.132996*ewgc("Cqq3R",1,1,2,2)
            +-0.067518*ewgc("Cqu1R",0,0,2,2)+0.020733*ewgc("Cqu1R",1,1,2,2)+-0.041823*ewgc("Cqu1R",2,2,0,0)
            +-0.014534*ewgc("Cqu1R",2,2,1,1)+0.134580*ewgc("Cud1R",2,2,0,0)+-0.004280*ewgc("Cud1R",2,2,1,1)
            +0.096208*ewgc("CuuR",0,0,2,2)+0.038326*ewgc("CuuR",1,1,2,2);
            double sigma_neg_bin_mtt_0_500_madgraph = 135.300000;
            double sigma_neg_bin_mtt_0_500_NP = -6.358540*ewgc("CG")+0.411874*ewgc("Cqd8R",2,2,0,0)
            +-0.006404*ewgc("Cqd8R",2,2,1,1)+2.642460*ewgc("Cqq1R",0,2,2,0)+0.063511*ewgc("Cqq1R",1,2,2,1)
            +5.970760*ewgc("Cqq3R",0,2,2,0)+0.714859*ewgc("Cqq3R",1,2,2,1)+1.615130*ewgc("Cqu8R",0,0,2,2)
            -0.079004*ewgc("Cqu8R",1,1,2,2)+0.944519*ewgc("Cqu8R",2,2,0,0)+0.028645*ewgc("Cqu8R",2,2,1,1)
            +0.301135*ewgc("Cud8R",2,2,0,0)+-0.078071*ewgc("Cud8R",2,2,1,1)-39.340100*ewgc("CuGR",2,2)
            +2.559810*ewgc("CuuR",0,2,2,0)+-0.060753*ewgc("CuuR",1,2,2,1)
            +-0.015693*ewgc("Cqd1R",2,2,0,0)+-0.008030*ewgc("Cqd1R",2,2,1,1)+0.123273*ewgc("Cqq1R",0,0,2,2)
            +0.052135*ewgc("Cqq1R",1,1,2,2)+0.387403*ewgc("Cqq3R",0,0,2,2)+-0.093296*ewgc("Cqq3R",1,1,2,2)
            +-0.012692*ewgc("Cqu1R",0,0,2,2)+-0.019275*ewgc("Cqu1R",1,1,2,2)+-0.015188*ewgc("Cqu1R",2,2,0,0)
            +0.021649*ewgc("Cqu1R",2,2,1,1)+0.001446*ewgc("Cud1R",2,2,0,0)+0.011419*ewgc("Cud1R",2,2,1,1)
            +-0.088095*ewgc("CuuR",0,0,2,2)+-0.031337*ewgc("CuuR",1,1,2,2);
            
            double sigma_pos_bin_mtt_0_500_NP_Corrected = SM_sigma_pos_bin_mtt_0_500*sigma_pos_bin_mtt_0_500_NP/sigma_pos_bin_mtt_0_500_madgraph;
            double sigma_neg_bin_mtt_0_500_NP_Corrected = SM_sigma_neg_bin_mtt_0_500*sigma_neg_bin_mtt_0_500_NP/sigma_neg_bin_mtt_0_500_madgraph;
            
            double NP_charge_asymmetry_num_bin_mtt_0_500 = sigma_pos_bin_mtt_0_500_NP_Corrected - sigma_neg_bin_mtt_0_500_NP_Corrected;
            double NP_charge_asymmetry_deno_bin_mtt_0_500 = sigma_pos_bin_mtt_0_500_NP_Corrected + sigma_neg_bin_mtt_0_500_NP_Corrected;;
            
            return  SM_charge_asymmetry_bin_mtt_0_500*(1+(NP_charge_asymmetry_num_bin_mtt_0_500-NP_charge_asymmetry_deno_bin_mtt_0_500)/SM_charge_asymmetry_deno_bin_mtt_0_500);            
        
        }
        
    } else if(b_min == 500 && b_max == 750){
        
        double SM_charge_asymmetry_bin_mtt_500_750 = SM.getOptionalParameter("SM_charge_asymmetry_bin_mtt_500_750");
        double SM_charge_asymmetry_deno_bin_mtt_500_750 = SM.getOptionalParameter("SM_charge_asymmetry_deno_bin_mtt_500_750");
        double SM_charge_asymmetry_num_bin_mtt_500_750 = SM_charge_asymmetry_bin_mtt_500_750*SM_charge_asymmetry_deno_bin_mtt_500_750;
        
        
        double SM_sigma_pos_bin_mtt_500_750 =0.5*(SM_charge_asymmetry_num_bin_mtt_500_750+SM_charge_asymmetry_deno_bin_mtt_500_750);
        double SM_sigma_neg_bin_mtt_500_750 =0.5*(SM_charge_asymmetry_deno_bin_mtt_500_750-SM_charge_asymmetry_num_bin_mtt_500_750);

        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            
            
            
            double sigma_pos_bin_mtt_500_750_madgraph = 80.589400;
            double sigma_pos_bin_mtt_500_750_NP = -6.999250*ewgc("CG")+0.386926*ewgc("Cqd8R",2,2,0,0)
            +0.078147*ewgc("Cqd8R",2,2,1,1)+4.115120*ewgc("Cqq1R",0,2,2,0)+0.207149*ewgc("Cqq1R",1,2,2,1)
            +8.367150*ewgc("Cqq3R",0,2,2,0)+0.726869*ewgc("Cqq3R",1,2,2,1)+0.814576*ewgc("Cqu8R",0,0,2,2)
            +0.182552*ewgc("Cqu8R",1,1,2,2)+0.471441*ewgc("Cqu8R",2,2,0,0)+0.020673*ewgc("Cqu8R",2,2,1,1)
            +0.603346*ewgc("Cud8R",2,2,0,0)+0.063159*ewgc("Cud8R",2,2,1,1)-20.819700*ewgc("CuGR",2,2)
            +3.980290*ewgc("CuuR",0,2,2,0)+0.169850*ewgc("CuuR",1,2,2,1)
            +0.047499*ewgc("Cqd1R",2,2,0,0)+0.051676*ewgc("Cqd1R",2,2,1,1)+0.322849*ewgc("Cqq1R",0,0,2,2)
            +0.058428*ewgc("Cqq1R",1,1,2,2)+0.909203*ewgc("Cqq3R",0,0,2,2)+0.028468*ewgc("Cqq3R",1,1,2,2)
            +0.083136*ewgc("Cqu1R",0,0,2,2)+0.009166*ewgc("Cqu1R",1,1,2,2)+0.095883*ewgc("Cqu1R",2,2,0,0)
            +0.046037*ewgc("Cqu1R",2,2,1,1)+-0.029524*ewgc("Cud1R",2,2,0,0)+0.062062*ewgc("Cud1R",2,2,1,1)
            +0.373091*ewgc("CuuR",0,0,2,2)+0.067338*ewgc("CuuR",1,1,2,2);
            double sigma_neg_bin_mtt_500_750_madgraph = 80.891300;
            double sigma_neg_bin_mtt_500_750_NP = -7.057070*ewgc("CG")+0.573031*ewgc("Cqd8R",2,2,0,0)
            -0.004786*ewgc("Cqd8R",2,2,1,1)+2.052850*ewgc("Cqq1R",0,2,2,0)+0.125976*ewgc("Cqq1R",1,2,2,1)
            +4.836950*ewgc("Cqq3R",0,2,2,0)+0.711224*ewgc("Cqq3R",1,2,2,1)+1.470530*ewgc("Cqu8R",0,0,2,2)
            +0.069589*ewgc("Cqu8R",1,1,2,2)+0.938294*ewgc("Cqu8R",2,2,0,0)+-0.045546*ewgc("Cqu8R",2,2,1,1)
            +0.347501*ewgc("Cud8R",2,2,0,0)+0.011230*ewgc("Cud8R",2,2,1,1)+-20.768900*ewgc("CuGR",2,2)
            +1.942090*ewgc("CuuR",0,2,2,0)+0.161845*ewgc("CuuR",1,2,2,1)
            +-0.030151*ewgc("Cqd1R",2,2,0,0)+-0.051728*ewgc("Cqd1R",2,2,1,1)+0.012765*ewgc("Cqq1R",0,0,2,2)
            +-0.097104*ewgc("Cqq1R",1,1,2,2)+0.424513*ewgc("Cqq3R",0,0,2,2)+-0.015512*ewgc("Cqq3R",1,1,2,2)
            +0.014810*ewgc("Cqu1R",0,0,2,2)+-0.010594*ewgc("Cqu1R",1,1,2,2)+0.012532*ewgc("Cqu1R",2,2,0,0)
            +-0.046820*ewgc("Cqu1R",2,2,1,1)+-0.078305*ewgc("Cud1R",2,2,0,0)+-0.072728*ewgc("Cud1R",2,2,1,1)
            +0.132974*ewgc("CuuR",0,0,2,2)+-0.067169*ewgc("CuuR",1,1,2,2);
            
            double sigma_pos_bin_mtt_500_750_NP_Corrected = SM_sigma_pos_bin_mtt_500_750*sigma_pos_bin_mtt_500_750_NP/sigma_pos_bin_mtt_500_750_madgraph;
            double sigma_neg_bin_mtt_500_750_NP_Corrected = SM_sigma_neg_bin_mtt_500_750*sigma_neg_bin_mtt_500_750_NP/sigma_neg_bin_mtt_500_750_madgraph;
            
            double NP_charge_asymmetry_num_bin_mtt_500_750 = sigma_pos_bin_mtt_500_750_NP_Corrected - sigma_neg_bin_mtt_500_750_NP_Corrected;
            double NP_charge_asymmetry_deno_bin_mtt_500_750 = sigma_pos_bin_mtt_500_750_NP_Corrected + sigma_neg_bin_mtt_500_750_NP_Corrected;;
            
            return  SM_charge_asymmetry_bin_mtt_500_750*(1+(NP_charge_asymmetry_num_bin_mtt_500_750-NP_charge_asymmetry_deno_bin_mtt_500_750)/SM_charge_asymmetry_deno_bin_mtt_500_750);            
        
            
            
            
            
        }
        
    }  else if(b_min == 750 && b_max == 1000){
        
        
        double SM_charge_asymmetry_deno_bin_mtt_750_1000 = SM.getOptionalParameter("SM_charge_asymmetry_deno_bin_mtt_750_1000");
        double SM_charge_asymmetry_bin_mtt_750_1000 = SM.getOptionalParameter("SM_charge_asymmetry_bin_mtt_750_1000");
        double SM_charge_asymmetry_num_bin_mtt_750_1000 = SM_charge_asymmetry_bin_mtt_750_1000*SM_charge_asymmetry_deno_bin_mtt_750_1000;
        
        double SM_sigma_pos_bin_mtt_750_1000 =0.5*(SM_charge_asymmetry_num_bin_mtt_750_1000+SM_charge_asymmetry_deno_bin_mtt_750_1000);
        double SM_sigma_neg_bin_mtt_750_1000 =0.5*(SM_charge_asymmetry_deno_bin_mtt_750_1000-SM_charge_asymmetry_num_bin_mtt_750_1000);

        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            
            
            
            double sigma_pos_bin_mtt_750_1000_madgraph = 16.986900;
            double sigma_pos_bin_mtt_750_1000_NP = -1.924360*ewgc("CG")+0.144619*ewgc("Cqd8R",2,2,0,0)
            +0.015897*ewgc("Cqd8R",2,2,1,1)+1.863020*ewgc("Cqq1R",0,2,2,0)+0.031168*ewgc("Cqq1R",1,2,2,1)
            +3.676260*ewgc("Cqq3R",0,2,2,0)+0.245826*ewgc("Cqq3R",1,2,2,1)+0.316581*ewgc("Cqu8R",0,0,2,2)
            +0.022092*ewgc("Cqu8R",1,1,2,2)+0.190067*ewgc("Cqu8R",2,2,0,0)+0.021485*ewgc("Cqu8R",2,2,1,1)
            +0.258731*ewgc("Cud8R",2,2,0,0)+0.013146*ewgc("Cud8R",2,2,1,1)-4.307000*ewgc("CuGR",2,2)
            +1.792720*ewgc("CuuR",0,2,2,0)+0.019910*ewgc("CuuR",1,2,2,1)
            +-0.016893*ewgc("Cqd1R",2,2,0,0)+-0.008497*ewgc("Cqd1R",2,2,1,1)+0.141702*ewgc("Cqq1R",0,0,2,2)
            +-0.012729*ewgc("Cqq1R",1,1,2,2)+0.438894*ewgc("Cqq3R",0,0,2,2)+0.019810*ewgc("Cqq3R",1,1,2,2)
            +0.021260*ewgc("Cqu1R",0,0,2,2)+-0.004833*ewgc("Cqu1R",1,1,2,2)+0.011461*ewgc("Cqu1R",2,2,0,0)
            +-0.003083*ewgc("Cqu1R",2,2,1,1)+-0.036348*ewgc("Cud1R",2,2,0,0)+-0.008414*ewgc("Cud1R",2,2,1,1)
            +0.152046*ewgc("CuuR",0,0,2,2)+-0.006108*ewgc("CuuR",1,1,2,2);
            
            double sigma_neg_bin_mtt_750_1000_madgraph = 16.965400;
            double sigma_neg_bin_mtt_750_1000_NP = -1.909310*ewgc("CG")+0.253028*ewgc("Cqd8R",2,2,0,0)
            +0.026868*ewgc("Cqd8R",2,2,1,1)+0.775801*ewgc("Cqq1R",0,2,2,0)+0.021341*ewgc("Cqq1R",1,2,2,1)
            +1.850950*ewgc("Cqq3R",0,2,2,0)+0.208119*ewgc("Cqq3R",1,2,2,1)+0.666539*ewgc("Cqu8R",0,0,2,2)
            +0.038662*ewgc("Cqu8R",1,1,2,2)+0.437622*ewgc("Cqu8R",2,2,0,0)+0.027679*ewgc("Cqu8R",2,2,1,1)
            +0.143565*ewgc("Cud8R",2,2,0,0)+0.028559*ewgc("Cud8R",2,2,1,1)+-4.304670*ewgc("CuGR",2,2)
            +0.748939*ewgc("CuuR",0,2,2,0)+0.029156*ewgc("CuuR",1,2,2,1)
            +-0.016259*ewgc("Cqd1R",2,2,0,0)+0.007655*ewgc("Cqd1R",2,2,1,1)+0.028036*ewgc("Cqq1R",0,0,2,2)
            +-0.000705*ewgc("Cqq1R",1,1,2,2)+0.186361*ewgc("Cqq3R",0,0,2,2)+0.044037*ewgc("Cqq3R",1,1,2,2)
            +0.037482*ewgc("Cqu1R",0,0,2,2)+0.005345*ewgc("Cqu1R",1,1,2,2)+0.040132*ewgc("Cqu1R",2,2,0,0)
            +0.002865*ewgc("Cqu1R",2,2,1,1)+-0.012822*ewgc("Cud1R",2,2,0,0)+0.007088*ewgc("Cud1R",2,2,1,1)
            +0.055812*ewgc("CuuR",0,0,2,2)+0.007647*ewgc("CuuR",1,1,2,2);
            
            double sigma_pos_bin_mtt_750_1000_NP_Corrected = SM_sigma_pos_bin_mtt_750_1000*sigma_pos_bin_mtt_750_1000_NP/sigma_pos_bin_mtt_750_1000_madgraph;
            double sigma_neg_bin_mtt_750_1000_NP_Corrected = SM_sigma_neg_bin_mtt_750_1000*sigma_neg_bin_mtt_750_1000_NP/sigma_neg_bin_mtt_750_1000_madgraph;
            
            double NP_charge_asymmetry_num_bin_mtt_750_1000 = sigma_pos_bin_mtt_750_1000_NP_Corrected - sigma_neg_bin_mtt_750_1000_NP_Corrected;
            double NP_charge_asymmetry_deno_bin_mtt_750_1000 = sigma_pos_bin_mtt_750_1000_NP_Corrected + sigma_neg_bin_mtt_750_1000_NP_Corrected;;
            
            return  SM_charge_asymmetry_bin_mtt_750_1000*(1+(NP_charge_asymmetry_num_bin_mtt_750_1000-NP_charge_asymmetry_deno_bin_mtt_750_1000)/SM_charge_asymmetry_deno_bin_mtt_750_1000);            
        
            
            
            
            
        }
        
    }  else if(b_min == 1000 && b_max == 1500){
        
        
        double SM_charge_asymmetry_deno_bin_mtt_1000_1500 = SM.getOptionalParameter("SM_charge_asymmetry_deno_bin_mtt_1000_1500");
        double SM_charge_asymmetry_bin_mtt_1000_1500 = SM.getOptionalParameter("SM_charge_asymmetry_bin_mtt_1000_1500");
        double SM_charge_asymmetry_num_bin_mtt_1000_1500 = SM_charge_asymmetry_bin_mtt_1000_1500*SM_charge_asymmetry_deno_bin_mtt_1000_1500;
        
        
        double SM_sigma_pos_bin_mtt_1000_1500 =0.5*(SM_charge_asymmetry_num_bin_mtt_1000_1500+SM_charge_asymmetry_deno_bin_mtt_1000_1500);
        double SM_sigma_neg_bin_mtt_1000_1500 =0.5*(SM_charge_asymmetry_deno_bin_mtt_1000_1500-SM_charge_asymmetry_num_bin_mtt_1000_1500);

        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            
            double sigma_pos_bin_mtt_1000_1500_madgraph = 5.965330;
            double sigma_pos_bin_mtt_1000_1500_NP = -0.732765*ewgc("CG")+0.095774*ewgc("Cqd8R",2,2,0,0)
            +0.009861*ewgc("Cqd8R",2,2,1,1)+1.407820*ewgc("Cqq1R",0,2,2,0)+0.023055*ewgc("Cqq1R",1,2,2,1)
            +2.789740*ewgc("Cqq3R",0,2,2,0)+0.124751*ewgc("Cqq3R",1,2,2,1)+0.226047*ewgc("Cqu8R",0,0,2,2)
            +0.012025*ewgc("Cqu8R",1,1,2,2)+0.132009*ewgc("Cqu8R",2,2,0,0)+0.002076*ewgc("Cqu8R",2,2,1,1)
            +0.176967*ewgc("Cud8R",2,2,0,0)+0.014820*ewgc("Cud8R",2,2,1,1)-1.511020*ewgc("CuGR",2,2)
            +1.381670*ewgc("CuuR",0,2,2,0)+0.022990*ewgc("CuuR",1,2,2,1)
            +0.000480*ewgc("Cqd1R",2,2,0,0)+-0.000506*ewgc("Cqd1R",2,2,1,1)+0.119062*ewgc("Cqq1R",0,0,2,2)
            +-0.007047*ewgc("Cqq1R",1,1,2,2)+0.356119*ewgc("Cqq3R",0,0,2,2)+0.012398*ewgc("Cqq3R",1,1,2,2)
            +0.012939*ewgc("Cqu1R",0,0,2,2)+0.000035*ewgc("Cqu1R",1,1,2,2)+0.003252*ewgc("Cqu1R",2,2,0,0)
            +-0.002001*ewgc("Cqu1R",2,2,1,1)+-0.013101*ewgc("Cud1R",2,2,0,0)+-0.000576*ewgc("Cud1R",2,2,1,1)
            +0.117162*ewgc("CuuR",0,0,2,2)+0.001190*ewgc("CuuR",1,1,2,2);
            double sigma_neg_bin_mtt_1000_1500_madgraph = 5.963920;
            double sigma_neg_bin_mtt_1000_1500_NP = -0.734104*ewgc("CG")+0.177167*ewgc("Cqd8R",2,2,0,0)
            +0.016263*ewgc("Cqd8R",2,2,1,1)+0.546178*ewgc("Cqq1R",0,2,2,0)+0.024062*ewgc("Cqq1R",1,2,2,1)
            +1.297260*ewgc("Cqq3R",0,2,2,0)+0.100018*ewgc("Cqq3R",1,2,2,1)+0.517319*ewgc("Cqu8R",0,0,2,2)
            +0.016952*ewgc("Cqu8R",1,1,2,2)+0.336176*ewgc("Cqu8R",2,2,0,0)+0.003311*ewgc("Cqu8R",2,2,1,1)
            +0.095957*ewgc("Cud8R",2,2,0,0)+0.011525*ewgc("Cud8R",2,2,1,1)-1.511200*ewgc("CuGR",2,2)
            +0.534920*ewgc("CuuR",0,2,2,0)+0.024264*ewgc("CuuR",1,2,2,1)
            +0.000672*ewgc("Cqd1R",2,2,0,0)+0.000477*ewgc("Cqd1R",2,2,1,1)+0.023176*ewgc("Cqq1R",0,0,2,2)
            +-0.006814*ewgc("Cqq1R",1,1,2,2)+0.147804*ewgc("Cqq3R",0,0,2,2)+0.009854*ewgc("Cqq3R",1,1,2,2)
            +0.029905*ewgc("Cqu1R",0,0,2,2)+-0.000840*ewgc("Cqu1R",1,1,2,2)+0.013292*ewgc("Cqu1R",2,2,0,0)
            +0.001980*ewgc("Cqu1R",2,2,1,1)+-0.004801*ewgc("Cud1R",2,2,0,0)+0.001036*ewgc("Cud1R",2,2,1,1)
            +0.045432*ewgc("CuuR",0,0,2,2)+0.000336*ewgc("CuuR",1,1,2,2);
            
            
            double sigma_pos_bin_mtt_1000_1500_NP_Corrected = SM_sigma_pos_bin_mtt_1000_1500*sigma_pos_bin_mtt_1000_1500_NP/sigma_pos_bin_mtt_1000_1500_madgraph;
            double sigma_neg_bin_mtt_1000_1500_NP_Corrected = SM_sigma_neg_bin_mtt_1000_1500*sigma_neg_bin_mtt_1000_1500_NP/sigma_neg_bin_mtt_1000_1500_madgraph;
            
            double NP_charge_asymmetry_num_bin_mtt_1000_1500 = sigma_pos_bin_mtt_1000_1500_NP_Corrected - sigma_neg_bin_mtt_1000_1500_NP_Corrected;
            double NP_charge_asymmetry_deno_bin_mtt_1000_1500 = sigma_pos_bin_mtt_1000_1500_NP_Corrected + sigma_neg_bin_mtt_1000_1500_NP_Corrected;;
            
            return  SM_charge_asymmetry_bin_mtt_1000_1500*(1+(NP_charge_asymmetry_num_bin_mtt_1000_1500-NP_charge_asymmetry_deno_bin_mtt_1000_1500)/SM_charge_asymmetry_deno_bin_mtt_1000_1500);            
        
            
            
            
            
        }
        
    }     else if(b_min == 1500 && b_max == 3000){
        
        double SM_charge_asymmetry_bin_mtt_1500_3000 = SM.getOptionalParameter("SM_charge_asymmetry_bin_mtt_1500_3000");
        double SM_charge_asymmetry_deno_bin_mtt_1500_3000 = SM.getOptionalParameter("SM_charge_asymmetry_deno_bin_mtt_1500_3000");
        double SM_charge_asymmetry_num_bin_mtt_1500_3000 = SM_charge_asymmetry_bin_mtt_1500_3000*SM_charge_asymmetry_deno_bin_mtt_1500_3000;
        
        double SM_sigma_pos_bin_mtt_1500_3000 =0.5*(SM_charge_asymmetry_num_bin_mtt_1500_3000+SM_charge_asymmetry_deno_bin_mtt_1500_3000);
        double SM_sigma_neg_bin_mtt_1500_3000 =0.5*(SM_charge_asymmetry_deno_bin_mtt_1500_3000-SM_charge_asymmetry_num_bin_mtt_1500_3000);
        
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            
            
            
            double sigma_pos_bin_mtt_1500_3000_madgraph = 0.910259;
            double sigma_pos_bin_mtt_1500_3000_NP = -0.116441*ewgc("CG")+0.045409*ewgc("Cqd8R",2,2,0,0)
            +0.003096*ewgc("Cqd8R",2,2,1,1)+0.711655*ewgc("Cqq1R",0,2,2,0)+0.007175*ewgc("Cqq1R",1,2,2,1)
            +1.350330*ewgc("Cqq3R",0,2,2,0)+0.039100*ewgc("Cqq3R",1,2,2,1)+0.110131*ewgc("Cqu8R",0,0,2,2)
            +0.005191*ewgc("Cqu8R",1,1,2,2)+0.064112*ewgc("Cqu8R",2,2,0,0)+0.002451*ewgc("Cqu8R",2,2,1,1)
            +0.084721*ewgc("Cud8R",2,2,0,0)+0.004134*ewgc("Cud8R",2,2,1,1)-0.231007*ewgc("CuGR",2,2)
            +0.691483*ewgc("CuuR",0,2,2,0)+0.006844*ewgc("CuuR",1,2,2,1)
            +-0.001536*ewgc("Cqd1R",2,2,0,0)+0.000002*ewgc("Cqd1R",2,2,1,1)+0.063201*ewgc("Cqq1R",0,0,2,2)
            +-0.001288*ewgc("Cqq1R",1,1,2,2)+0.178152*ewgc("Cqq3R",0,0,2,2)+0.004233*ewgc("Cqq3R",1,1,2,2)
            +0.005775*ewgc("Cqu1R",0,0,2,2)+0.000328*ewgc("Cqu1R",1,1,2,2)+0.003439*ewgc("Cqu1R",2,2,0,0)
            +-0.000036*ewgc("Cqu1R",2,2,1,1)+-0.007397*ewgc("Cud1R",2,2,0,0)+-0.000466*ewgc("Cud1R",2,2,1,1)
            +0.059131*ewgc("CuuR",0,0,2,2)+0.000437*ewgc("CuuR",1,1,2,2);
            double sigma_neg_bin_mtt_1500_3000_madgraph = 0.909569;
            double sigma_neg_bin_mtt_1500_3000_NP = -0.116979*ewgc("CG")+0.084970*ewgc("Cqd8R",2,2,0,0)
            +0.004067*ewgc("Cqd8R",2,2,1,1)+0.274254*ewgc("Cqq1R",0,2,2,0)+0.007000*ewgc("Cqq1R",1,2,2,1)
            +0.618997*ewgc("Cqq3R",0,2,2,0)+0.032143*ewgc("Cqq3R",1,2,2,1)+0.253632*ewgc("Cqu8R",0,0,2,2)
            +0.006071*ewgc("Cqu8R",1,1,2,2)+0.167569*ewgc("Cqu8R",2,2,0,0)+0.002373*ewgc("Cqu8R",2,2,1,1)
            +0.045316*ewgc("Cud8R",2,2,0,0)+0.003020*ewgc("Cud8R",2,2,1,1)-0.231248*ewgc("CuGR",2,2)
            +0.265616*ewgc("CuuR",0,2,2,0)+0.006959*ewgc("CuuR",1,2,2,1)
            +-0.002375*ewgc("Cqd1R",2,2,0,0)+-0.000039*ewgc("Cqd1R",2,2,1,1)+0.015305*ewgc("Cqq1R",0,0,2,2)
            +-0.001083*ewgc("Cqq1R",1,1,2,2)+0.075644*ewgc("Cqq3R",0,0,2,2)+0.003511*ewgc("Cqq3R",1,1,2,2)
            +0.011743*ewgc("Cqu1R",0,0,2,2)+0.000548*ewgc("Cqu1R",1,1,2,2)+0.008199*ewgc("Cqu1R",2,2,0,0)
            +0.000069*ewgc("Cqu1R",2,2,1,1)+-0.003571*ewgc("Cud1R",2,2,0,0)+-0.000568*ewgc("Cud1R",2,2,1,1)
            +0.022219*ewgc("CuuR",0,0,2,2)+0.000486*ewgc("CuuR",1,1,2,2);
            
            
            
            double sigma_pos_bin_mtt_1500_3000_NP_Corrected = SM_sigma_pos_bin_mtt_1500_3000*sigma_pos_bin_mtt_1500_3000_NP/sigma_pos_bin_mtt_1500_3000_madgraph;
            double sigma_neg_bin_mtt_1500_3000_NP_Corrected = SM_sigma_neg_bin_mtt_1500_3000*sigma_neg_bin_mtt_1500_3000_NP/sigma_neg_bin_mtt_1500_3000_madgraph;
            
            double NP_charge_asymmetry_num_bin_mtt_1500_3000 = sigma_pos_bin_mtt_1500_3000_NP_Corrected - sigma_neg_bin_mtt_1500_3000_NP_Corrected;
            double NP_charge_asymmetry_deno_bin_mtt_1500_3000 = sigma_pos_bin_mtt_1500_3000_NP_Corrected + sigma_neg_bin_mtt_1500_3000_NP_Corrected;;
            
            return  SM_charge_asymmetry_bin_mtt_1500_3000*(1+(NP_charge_asymmetry_num_bin_mtt_1500_3000-NP_charge_asymmetry_deno_bin_mtt_1500_3000)/SM_charge_asymmetry_deno_bin_mtt_1500_3000);            
            //return 0.;

        }
        
    }
    
    
    else {
        throw std::runtime_error("\nERROR: Please specify a correct binning range for charge_asymmetry_tt_diff_mtt_LO.\n");
    }

}



sigma_tta_diff_LO_ATLAS_emu_200706946::sigma_tta_diff_LO_ATLAS_emu_200706946(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{

    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tta_bin_20_25_ATLAS_emu" << "SM_sigma_tta_bin_25_30_ATLAS_emu"
            << "SM_sigma_tta_bin_30_35_ATLAS_emu" << "SM_sigma_tta_bin_35_40_ATLAS_emu" << "SM_sigma_tta_bin_40_47_ATLAS_emu" 
            << "SM_sigma_tta_bin_47_55_ATLAS_emu" << "SM_sigma_tta_bin_55_70_ATLAS_emu" << "SM_sigma_tta_bin_70_85_ATLAS_emu"
            << "SM_sigma_tta_bin_85_132_ATLAS_emu" << "SM_sigma_tta_bin_132_180_ATLAS_emu" << "SM_sigma_tta_bin_180_300_ATLAS_emu");

    
}


double sigma_tta_diff_LO_ATLAS_emu_200706946::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    

    if(b_min == 20 && b_max == 25){
        
        double SM_sigma_tta_bin_20_25 = SM.getOptionalParameter("SM_sigma_tta_bin_20_25_ATLAS_emu");
        double sigma_tta_bin_20_25_madgraph = 0.247297; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            //There is no dependence on ewgc("Cqd8R",2,2,2,2) neither on ewgc("Cud8R",2,2,2,2) with the precision considered
        
            double total = SM_sigma_tta_bin_20_25 +(-0.014297*ewgc("CG")+0.002013*ewgc("Cqd8R",2,2,0,0)
                    +0.000242*ewgc("Cqd8R",2,2,1,1)+0.035750*ewgc("Cqq1R",0,2,2,0)
                    +0.001121*ewgc("Cqq1R",1,2,2,1)+0.062345*ewgc("Cqq3R",0,2,2,0)
                    +0.004594*ewgc("Cqq3R",1,2,2,1)+0.015299*ewgc("Cqu8R",0,0,2,2)
                    +0.000694*ewgc("Cqu8R",1,1,2,2)+0.013298*ewgc("Cqu8R",2,2,0,0)
                    +0.000390*ewgc("Cqu8R",2,2,1,1)+-0.002389*ewgc("CuBR",2,2)
                    +0.003507*ewgc("Cud8R",2,2,0,0)+0.000447*ewgc("Cud8R",2,2,1,1)
                    -0.065862*ewgc("CuGR",2,2)+0.034943*ewgc("CuuR",0,2,2,0)
                    +0.001094*ewgc("CuuR",1,2,2,1)+-0.001764*ewgc("CuWR",2,2)
                    +-0.026156*ewgc("CHD")+-0.015017*ewgc("CHl3R",0,0)+-0.015017*ewgc("CHl3R",1,1)
                    +-0.056084*ewgc("CHWB")+0.015096*ewgc("CllR",0,1,1,0)+-0.000133*ewgc("Cqd1R",2,2,0,0)
                    +-0.000015*ewgc("Cqd1R",2,2,1,1)+0.002878*ewgc("Cqq1R",0,0,2,2)+-0.000104*ewgc("Cqq1R",1,1,2,2)
                    +0.007066*ewgc("Cqq3R",0,0,2,2)+0.000416*ewgc("Cqq3R",1,1,2,2)+0.001101*ewgc("Cqu1R",0,0,2,2)
                    +0.000001*ewgc("Cqu1R",1,1,2,2)+0.000850*ewgc("Cqu1R",2,2,0,0)+0.000041*ewgc("Cqu1R",2,2,1,1)
                    +-0.000304*ewgc("Cud1R",2,2,0,0)+-0.000062*ewgc("Cud1R",2,2,1,1)+0.002439*ewgc("CuuR",0,0,2,2)
                    +0.000055*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_tta_bin_20_25/sigma_tta_bin_20_25_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
            
        }
    } else if(b_min == 25 && b_max == 30){
        
    
        double SM_sigma_tta_bin_25_30 = SM.getOptionalParameter("SM_sigma_tta_bin_25_30_ATLAS_emu");
        double sigma_tta_bin_25_30_madgraph = 0.187804; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tta_bin_25_30 +(-0.010872*ewgc("CG")+0.001566*ewgc("Cqd8R",2,2,0,0)
                    +0.000215*ewgc("Cqd8R",2,2,1,1)+0.026408*ewgc("Cqq1R",0,2,2,0)
                    +0.000870*ewgc("Cqq1R",1,2,2,1)+0.047162*ewgc("Cqq3R",0,2,2,0)
                    +0.003478*ewgc("Cqq3R",1,2,2,1)+0.000532*ewgc("Cqu8R",1,1,2,2)
                    +0.010073*ewgc("Cqu8R",2,2,0,0)+0.000333*ewgc("Cqu8R",2,2,1,1)
                    +-0.002372*ewgc("CuBR",2,2)+0.002749*ewgc("Cud8R",2,2,0,0)
                    +0.000356*ewgc("Cud8R",2,2,1,1)-0.049862*ewgc("CuGR",2,2)
                    +0.025851*ewgc("CuuR",0,2,2,0)+0.000867*ewgc("CuuR",1,2,2,1)
                    +-0.001656*ewgc("CuWR",2,2)
                    +-0.019822*ewgc("CHD")+-0.011396*ewgc("CHl3R",0,0)+-0.011396*ewgc("CHl3R",1,1)
                    +-0.042561*ewgc("CHWB")+0.011447*ewgc("CllR",0,1,1,0)+-0.000067*ewgc("Cqd1R",2,2,0,0)
                    +0.000071*ewgc("Cqd1R",2,2,1,1)+0.002047*ewgc("Cqq1R",0,0,2,2)
                    +-0.000040*ewgc("Cqq1R",1,1,2,2)+0.005286*ewgc("Cqq3R",0,0,2,2)
                    +0.000378*ewgc("Cqq3R",1,1,2,2)+0.000934*ewgc("Cqu1R",0,0,2,2)
                    +0.000053*ewgc("Cqu1R",1,1,2,2)+0.000629*ewgc("Cqu1R",2,2,0,0)
                    +0.000070*ewgc("Cqu1R",2,2,1,1)+-0.000164*ewgc("Cud1R",2,2,0,0)
                    +0.000048*ewgc("Cud1R",2,2,1,1)+0.001741*ewgc("CuuR",0,0,2,2)
                    +0.000125*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_tta_bin_25_30/sigma_tta_bin_25_30_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
    } else if(b_min == 30 && b_max == 35){
        
    
        double SM_sigma_tta_bin_30_35 = SM.getOptionalParameter("SM_sigma_tta_bin_30_35_ATLAS_emu");
        double sigma_tta_bin_30_35_madgraph = 0.148458; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tta_bin_30_35 +(-0.008634*ewgc("CG")+0.001254*ewgc("Cqd8R",2,2,0,0)
                    +0.000149*ewgc("Cqd8R",2,2,1,1)+0.020344*ewgc("Cqq1R",0,2,2,0)
                    +0.000623*ewgc("Cqq1R",1,2,2,1)+0.036940*ewgc("Cqq3R",0,2,2,0)
                    +0.002683*ewgc("Cqq3R",1,2,2,1)+0.009143*ewgc("Cqu8R",0,0,2,2)
                    +0.000393*ewgc("Cqu8R",1,1,2,2)+0.007956*ewgc("Cqu8R",2,2,0,0)
                    +0.000237*ewgc("Cqu8R",2,2,1,1)+-0.002334*ewgc("CuBR",2,2)
                    +0.002172*ewgc("Cud8R",2,2,0,0)+0.000281*ewgc("Cud8R",2,2,1,1)
                    -0.039223*ewgc("CuGR",2,2)+0.019888*ewgc("CuuR",0,2,2,0)
                    +0.000612*ewgc("CuuR",1,2,2,1)+-0.001536*ewgc("CuWR",2,2)
                    +-0.015702*ewgc("CHD")+-0.009016*ewgc("CHl3R",0,0)+-0.009016*ewgc("CHl3R",1,1)
                    +-0.033656*ewgc("CHWB")+0.009010*ewgc("CllR",0,1,1,0)+-0.000017*ewgc("Cqd1R",2,2,0,0)
                    +-0.000017*ewgc("Cqd1R",2,2,1,1)+0.001591*ewgc("Cqq1R",0,0,2,2)+-0.000074*ewgc("Cqq1R",1,1,2,2)
                    +0.004182*ewgc("Cqq3R",0,0,2,2)+0.000278*ewgc("Cqq3R",1,1,2,2)+0.000703*ewgc("Cqu1R",0,0,2,2)
                    +0.000024*ewgc("Cqu1R",1,1,2,2)+0.000526*ewgc("Cqu1R",2,2,0,0)+0.000020*ewgc("Cqu1R",2,2,1,1)
                    +-0.000162*ewgc("Cud1R",2,2,0,0)+-0.000027*ewgc("Cud1R",2,2,1,1)+0.001461*ewgc("CuuR",0,0,2,2)
                    +0.000046*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_tta_bin_30_35/sigma_tta_bin_30_35_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
    } else if(b_min == 35 && b_max == 40){
        
    
        double SM_sigma_tta_bin_35_40 = SM.getOptionalParameter("SM_sigma_tta_bin_35_40_ATLAS_emu");
        double sigma_tta_bin_35_40_madgraph = 0.120694; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            double total = SM_sigma_tta_bin_35_40 +(-0.007257*ewgc("CG")+0.000987*ewgc("Cqd8R",2,2,0,0)
                    +0.000102*ewgc("Cqd8R",2,2,1,1)+0.016185*ewgc("Cqq1R",0,2,2,0)
                    +0.000464*ewgc("Cqq1R",1,2,2,1)+0.029650*ewgc("Cqq3R",0,2,2,0)
                    +0.002137*ewgc("Cqq3R",1,2,2,1)+0.007430*ewgc("Cqu8R",0,0,2,2)
                    +0.000309*ewgc("Cqu8R",1,1,2,2)+0.006422*ewgc("Cqu8R",2,2,0,0)
                    +0.000183*ewgc("Cqu8R",2,2,1,1)+-0.002287*ewgc("CuBR",2,2)
                    +0.001776*ewgc("Cud8R",2,2,0,0)+0.000208*ewgc("Cud8R",2,2,1,1)
                    -0.031798*ewgc("CuGR",2,2)+0.015826*ewgc("CuuR",0,2,2,0)
                    +0.000435*ewgc("CuuR",1,2,2,1)+-0.001438*ewgc("CuWR",2,2)
                    +-0.012765*ewgc("CHD")+-0.007346*ewgc("CHl3R",0,0)+-0.007346*ewgc("CHl3R",1,1)
                    +-0.027384*ewgc("CHWB")+0.007318*ewgc("CllR",0,1,1,0)+-0.000040*ewgc("Cqd1R",2,2,0,0)
                    +0.000008*ewgc("Cqd1R",2,2,1,1)+0.001124*ewgc("Cqq1R",0,0,2,2)+-0.000087*ewgc("Cqq1R",1,1,2,2)
                    +0.003280*ewgc("Cqq3R",0,0,2,2)+0.000171*ewgc("Cqq3R",1,1,2,2)+0.000558*ewgc("Cqu1R",0,0,2,2)
                    +0.000005*ewgc("Cqu1R",1,1,2,2)+0.000385*ewgc("Cqu1R",2,2,0,0)+-0.000001*ewgc("Cqu1R",2,2,1,1)
                    +-0.000172*ewgc("Cud1R",2,2,0,0)+-0.000036*ewgc("Cud1R",2,2,1,1)+0.001123*ewgc("CuuR",0,0,2,2)
                    +0.000027*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_tta_bin_35_40/sigma_tta_bin_35_40_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
    } else if(b_min == 40 && b_max == 47){
        
    
        double SM_sigma_tta_bin_40_47 = SM.getOptionalParameter("SM_sigma_tta_bin_40_47_ATLAS_emu");
        double sigma_tta_bin_40_47_madgraph = 0.135826; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{

            double total = SM_sigma_tta_bin_40_47 +(-0.008293*ewgc("CG")+0.001098*ewgc("Cqd8R",2,2,0,0)
                    +0.000158*ewgc("Cqd8R",2,2,1,1)+0.017832*ewgc("Cqq1R",0,2,2,0)
                    +0.000514*ewgc("Cqq1R",1,2,2,1)+0.033274*ewgc("Cqq3R",0,2,2,0)
                    +0.002414*ewgc("Cqq3R",1,2,2,1)+0.008356*ewgc("Cqu8R",0,0,2,2)
                    +0.000370*ewgc("Cqu8R",1,1,2,2)+0.007249*ewgc("Cqu8R",2,2,0,0)
                    +0.000232*ewgc("Cqu8R",2,2,1,1)-0.00307051*ewgc("CuBR",2,2)
                    +0.002034*ewgc("Cud8R",2,2,0,0)+0.000278*ewgc("Cud8R",2,2,1,1)
                    -0.035596*ewgc("CuGR",2,2)+0.017419*ewgc("CuuR",0,2,2,0)
                    +0.000515*ewgc("CuuR",1,2,2,1)-0.001898*ewgc("CuWR",2,2)
                    +-0.014369*ewgc("CHD")+-0.008254*ewgc("CHl3R",0,0)+-0.008254*ewgc("CHl3R",1,1)
                    +-0.030793*ewgc("CHWB")+0.008271*ewgc("CllR",0,1,1,0)+-0.000044*ewgc("Cqd1R",2,2,0,0)
                    +0.000025*ewgc("Cqd1R",2,2,1,1)+0.001219*ewgc("Cqq1R",0,0,2,2)+-0.000075*ewgc("Cqq1R",1,1,2,2)
                    +0.003664*ewgc("Cqq3R",0,0,2,2)+0.000239*ewgc("Cqq3R",1,1,2,2)+0.000638*ewgc("Cqu1R",0,0,2,2)
                    +0.000039*ewgc("Cqu1R",1,1,2,2)+0.000460*ewgc("Cqu1R",2,2,0,0)+0.000023*ewgc("Cqu1R",2,2,1,1)
                    +-0.000134*ewgc("Cud1R",2,2,0,0)+-0.000025*ewgc("Cud1R",2,2,1,1)+0.001269*ewgc("CuuR",0,0,2,2)
                    +0.000027*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_tta_bin_40_47/sigma_tta_bin_40_47_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
    } else if(b_min == 47 && b_max == 55){
        
    
        double SM_sigma_tta_bin_47_55 = SM.getOptionalParameter("SM_sigma_tta_bin_47_55_ATLAS_emu");
        double sigma_tta_bin_47_55_madgraph = 0.121740; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tta_bin_47_55 +(-0.007575*ewgc("CG")+0.001001*ewgc("Cqd8R",2,2,0,0)
                    +0.000125*ewgc("Cqd8R",2,2,1,1)+0.015576*ewgc("Cqq1R",0,2,2,0)
                    +0.000477*ewgc("Cqq1R",1,2,2,1)+0.029798*ewgc("Cqq3R",0,2,2,0)
                    +0.002128*ewgc("Cqq3R",1,2,2,1)+0.007469*ewgc("Cqu8R",0,0,2,2)
                    +0.000318*ewgc("Cqu8R",1,1,2,2)+0.006480*ewgc("Cqu8R",2,2,0,0)
                    +0.000197*ewgc("Cqu8R",2,2,1,1)+-0.003278*ewgc("CuBR",2,2)
                    +0.001854*ewgc("Cud8R",2,2,0,0)+0.000228*ewgc("Cud8R",2,2,1,1)
                    -0.031752*ewgc("CuGR",2,2)+0.015239*ewgc("CuuR",0,2,2,0)
                    +0.000455*ewgc("CuuR",1,2,2,1)+-0.001958*ewgc("CuWR",2,2)
                    +-0.012849*ewgc("CHD")+-0.007408*ewgc("CHl3R",0,0)+-0.007408*ewgc("CHl3R",1,1)
                    +-0.027615*ewgc("CHWB")+0.007413*ewgc("CllR",0,1,1,0)+-0.000032*ewgc("Cqd1R",2,2,0,0)
                    +-0.000007*ewgc("Cqd1R",2,2,1,1)+0.001099*ewgc("Cqq1R",0,0,2,2)+-0.000060*ewgc("Cqq1R",1,1,2,2)
                    +0.003327*ewgc("Cqq3R",0,0,2,2)+0.000190*ewgc("Cqq3R",1,1,2,2)+0.000557*ewgc("Cqu1R",0,0,2,2)
                    +0.000035*ewgc("Cqu1R",1,1,2,2)+0.000418*ewgc("Cqu1R",2,2,0,0)+0.000028*ewgc("Cqu1R",2,2,1,1)
                    +-0.000150*ewgc("Cud1R",2,2,0,0)+-0.000001*ewgc("Cud1R",2,2,1,1)+0.001152*ewgc("CuuR",0,0,2,2)
                    +0.000044*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_tta_bin_47_55/sigma_tta_bin_47_55_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
    } else if(b_min == 55 && b_max == 70){
        
    
        double SM_sigma_tta_bin_55_70 = SM.getOptionalParameter("SM_sigma_tta_bin_55_70_ATLAS_emu");
        double sigma_tta_bin_55_70_madgraph = 0.165952; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tta_bin_55_70 +(+-0.010297*ewgc("CG")+0.001358*ewgc("Cqd8R",2,2,0,0)
                    +0.000129*ewgc("Cqd8R",2,2,1,1)+0.020636*ewgc("Cqq1R",0,2,2,0)
                    +0.000570*ewgc("Cqq1R",1,2,2,1)+0.039909*ewgc("Cqq3R",0,2,2,0)
                    +0.002860*ewgc("Cqq3R",1,2,2,1)+0.010263*ewgc("Cqu8R",0,0,2,2)
                    +0.000378*ewgc("Cqu8R",1,1,2,2)+0.008906*ewgc("Cqu8R",2,2,0,0)
                    +0.000208*ewgc("Cqu8R",2,2,1,1)+-0.005566*ewgc("CuBR",2,2)
                    +0.002574*ewgc("Cud8R",2,2,0,0)+0.000267*ewgc("Cud8R",2,2,1,1)
                    -0.042950*ewgc("CuGR",2,2)+0.020254*ewgc("CuuR",0,2,2,0)
                    +0.000551*ewgc("CuuR",1,2,2,1)+-0.003266*ewgc("CuWR",2,2)
                    +-0.017597*ewgc("CHD")+-0.010136*ewgc("CHl3R",0,0)+-0.010136*ewgc("CHl3R",1,1)
                    +-0.037685*ewgc("CHWB")+0.010050*ewgc("CllR",0,1,1,0)+-0.000022*ewgc("Cqd1R",2,2,0,0)
                    +-0.000027*ewgc("Cqd1R",2,2,1,1)+0.001286*ewgc("Cqq1R",0,0,2,2)+-0.000082*ewgc("Cqq1R",1,1,2,2)
                    +0.004444*ewgc("Cqq3R",0,0,2,2)+0.000294*ewgc("Cqq3R",1,1,2,2)+0.000794*ewgc("Cqu1R",0,0,2,2)
                    +-0.000025*ewgc("Cqu1R",1,1,2,2)+0.000629*ewgc("Cqu1R",2,2,0,0)+0.000005*ewgc("Cqu1R",2,2,1,1)
                    +-0.000205*ewgc("Cud1R",2,2,0,0)+-0.000035*ewgc("Cud1R",2,2,1,1)+0.001458*ewgc("CuuR",0,0,2,2)
                    +0.000014*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_tta_bin_55_70/sigma_tta_bin_55_70_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
        }
    } else if(b_min == 70 && b_max == 85){
        
    
        double SM_sigma_tta_bin_70_85 = SM.getOptionalParameter("SM_sigma_tta_bin_70_85_ATLAS_emu");
        double sigma_tta_bin_70_85_madgraph = 0.115196; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tta_bin_70_85 +(-0.007459*ewgc("CG")+0.000912*ewgc("Cqd8R",2,2,0,0)
                    +0.000085*ewgc("Cqd8R",2,2,1,1)+0.013773*ewgc("Cqq1R",0,2,2,0)
                    +0.000343*ewgc("Cqq1R",1,2,2,1)+0.027638*ewgc("Cqq3R",0,2,2,0)
                    +0.001935*ewgc("Cqq3R",1,2,2,1)+0.007129*ewgc("Cqu8R",0,0,2,2)
                    +0.000241*ewgc("Cqu8R",1,1,2,2)+0.006201*ewgc("Cqu8R",2,2,0,0)
                    +0.000146*ewgc("Cqu8R",2,2,1,1)+-0.004693*ewgc("CuBR",2,2)
                    +0.001789*ewgc("Cud8R",2,2,0,0)+0.000183*ewgc("Cud8R",2,2,1,1)
                    -0.029503*ewgc("CuGR",2,2)+0.013479*ewgc("CuuR",0,2,2,0)
                    +0.000332*ewgc("CuuR",1,2,2,1)+-0.002681*ewgc("CuWR",2,2)
                    +-0.012165*ewgc("CHD")+-0.007010*ewgc("CHl3R",0,0)+-0.007010*ewgc("CHl3R",1,1)
                    +-0.026130*ewgc("CHWB")+0.006989*ewgc("CllR",0,1,1,0)+-0.000044*ewgc("Cqd1R",2,2,0,0)
                    +-0.000025*ewgc("Cqd1R",2,2,1,1)+0.000816*ewgc("Cqq1R",0,0,2,2)+-0.000107*ewgc("Cqq1R",1,1,2,2)
                    +0.003080*ewgc("Cqq3R",0,0,2,2)+0.000154*ewgc("Cqq3R",1,1,2,2)+0.000491*ewgc("Cqu1R",0,0,2,2)
                    +-0.000015*ewgc("Cqu1R",1,1,2,2)+0.000390*ewgc("Cqu1R",2,2,0,0)+-0.000031*ewgc("Cqu1R",2,2,1,1)
                    +-0.000158*ewgc("Cud1R",2,2,0,0)+-0.000038*ewgc("Cud1R",2,2,1,1)+0.000956*ewgc("CuuR",0,0,2,2)
                    +-0.000004*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_tta_bin_70_85/sigma_tta_bin_70_85_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
            
        }
    } else if(b_min == 85 && b_max == 132){
        
    
        double SM_sigma_tta_bin_85_132 = SM.getOptionalParameter("SM_sigma_tta_bin_85_132_ATLAS_emu");
        double sigma_tta_bin_85_132_madgraph = 0.200083; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tta_bin_85_132 +(-0.014099*ewgc("CG")+0.001631*ewgc("Cqd8R",2,2,0,0)
                    +0.000171*ewgc("Cqd8R",2,2,1,1)+0.023131*ewgc("Cqq1R",0,2,2,0)
                    +0.000580*ewgc("Cqq1R",1,2,2,1)+0.048084*ewgc("Cqq3R",0,2,2,0)
                    +0.003229*ewgc("Cqq3R",1,2,2,1)+0.012662*ewgc("Cqu8R",0,0,2,2)
                    +0.000456*ewgc("Cqu8R",1,1,2,2)+0.011006*ewgc("Cqu8R",2,2,0,0)
                    +0.000281*ewgc("Cqu8R",2,2,1,1)+-0.010083*ewgc("CuBR",2,2)
                    +0.003292*ewgc("Cud8R",2,2,0,0)+0.000342*ewgc("Cud8R",2,2,1,1)
                    -0.050663*ewgc("CuGR",2,2)+0.022549*ewgc("CuuR",0,2,2,0)
                    +0.000559*ewgc("CuuR",1,2,2,1)+-0.005730*ewgc("CuWR",2,2)
                    +-0.021149*ewgc("CHD")+-0.012134*ewgc("CHl3R",0,0)+-0.012134*ewgc("CHl3R",1,1)
                    +-0.045370*ewgc("CHWB")+0.012130*ewgc("CllR",0,1,1,0)+-0.000086*ewgc("Cqd1R",2,2,0,0)
                    +-0.000004*ewgc("Cqd1R",2,2,1,1)+0.001348*ewgc("Cqq1R",0,0,2,2)+-0.000144*ewgc("Cqq1R",1,1,2,2)
                    +0.005438*ewgc("Cqq3R",0,0,2,2)+0.000306*ewgc("Cqq3R",1,1,2,2)+0.000845*ewgc("Cqu1R",0,0,2,2)
                    +0.000045*ewgc("Cqu1R",1,1,2,2)+0.000691*ewgc("Cqu1R",2,2,0,0)+0.000022*ewgc("Cqu1R",2,2,1,1)
                    +-0.000280*ewgc("Cud1R",2,2,0,0)+-0.000023*ewgc("Cud1R",2,2,1,1)+0.001599*ewgc("CuuR",0,0,2,2)
                    +0.000033*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_tta_bin_85_132/sigma_tta_bin_85_132_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
            
        }
    } else if(b_min == 132 && b_max == 180){
        
    
        double SM_sigma_tta_bin_132_180 = SM.getOptionalParameter("SM_sigma_tta_bin_132_180_ATLAS_emu");
        double sigma_tta_bin_132_180_madgraph = 0.091916; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tta_bin_132_180 +(-0.007853*ewgc("CG")+0.000817*ewgc("Cqd8R",2,2,0,0)
                    +0.000091*ewgc("Cqd8R",2,2,1,1)+0.010609*ewgc("Cqq1R",0,2,2,0)
                    +0.000257*ewgc("Cqq1R",1,2,2,1)+0.023194*ewgc("Cqq3R",0,2,2,0)
                    +0.001494*ewgc("Cqq3R",1,2,2,1)+0.006217*ewgc("Cqu8R",0,0,2,2)
                    +0.000215*ewgc("Cqu8R",1,1,2,2)+0.005413*ewgc("Cqu8R",2,2,0,0)
                    +0.000142*ewgc("Cqu8R",2,2,1,1)+-0.005438*ewgc("CuBR",2,2)
                    +0.001678*ewgc("Cud8R",2,2,0,0)+0.000183*ewgc("Cud8R",2,2,1,1)
                    -0.023098*ewgc("CuGR",2,2)+0.010347*ewgc("CuuR",0,2,2,0)
                    +0.000251*ewgc("CuuR",1,2,2,1)+-0.003075*ewgc("CuWR",2,2)
                    +-0.009700*ewgc("CHD")+-0.005572*ewgc("CHl3R",0,0)+-0.005572*ewgc("CHl3R",1,1)
                    +-0.020833*ewgc("CHWB")+0.005601*ewgc("CllR",0,1,1,0)+-0.000015*ewgc("Cqd1R",2,2,0,0)
                    +0.000014*ewgc("Cqd1R",2,2,1,1)+0.000597*ewgc("Cqq1R",0,0,2,2)+-0.000047*ewgc("Cqq1R",1,1,2,2)
                    +0.002666*ewgc("Cqq3R",0,0,2,2)+0.000156*ewgc("Cqq3R",1,1,2,2)+0.000449*ewgc("Cqu1R",0,0,2,2)
                    +0.000026*ewgc("Cqu1R",1,1,2,2)+0.000336*ewgc("Cqu1R",2,2,0,0)+0.000024*ewgc("Cqu1R",2,2,1,1)
                    +-0.000121*ewgc("Cud1R",2,2,0,0)+0.000002*ewgc("Cud1R",2,2,1,1)+0.000861*ewgc("CuuR",0,0,2,2)
                    +0.000036*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_tta_bin_132_180/sigma_tta_bin_132_180_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
    } else if(b_min == 180 && b_max == 300){
        
    
        double SM_sigma_tta_bin_180_300 = SM.getOptionalParameter("SM_sigma_tta_bin_180_300_ATLAS_emu");
        double sigma_tta_bin_180_300_madgraph = 0.077596; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tta_bin_180_300 +(-0.009519*ewgc("CG")+0.000800*ewgc("Cqd8R",2,2,0,0)
                    +0.000069*ewgc("Cqd8R",2,2,1,1)+0.009803*ewgc("Cqq1R",0,2,2,0)
                    +0.000190*ewgc("Cqq1R",1,2,2,1)+0.022629*ewgc("Cqq3R",0,2,2,0)
                    +0.001298*ewgc("Cqq3R",1,2,2,1)+0.006131*ewgc("Cqu8R",0,0,2,2)
                    +0.000175*ewgc("Cqu8R",1,1,2,2)+0.005339*ewgc("Cqu8R",2,2,0,0)
                    +0.000108*ewgc("Cqu8R",2,2,1,1)+-0.004601*ewgc("CuBR",2,2)
                    +0.001701*ewgc("Cud8R",2,2,0,0)+0.000146*ewgc("Cud8R",2,2,1,1)
                    -0.019949*ewgc("CuGR",2,2)+0.009534*ewgc("CuuR",0,2,2,0)
                    +0.000184*ewgc("CuuR",1,2,2,1)+-0.002603*ewgc("CuWR",2,2)
                    +-0.008197*ewgc("CHD")+-0.004712*ewgc("CHl3R",0,0)+-0.004712*ewgc("CHl3R",1,1)
                    +-0.017577*ewgc("CHWB")+0.004709*ewgc("CllR",0,1,1,0)+-0.000013*ewgc("Cqd1R",2,2,0,0)
                    +0.000001*ewgc("Cqd1R",2,2,1,1)+0.000424*ewgc("Cqq1R",0,0,2,2)+-0.000066*ewgc("Cqq1R",1,1,2,2)
                    +0.002602*ewgc("Cqq3R",0,0,2,2)+0.000120*ewgc("Cqq3R",1,1,2,2)+0.000448*ewgc("Cqu1R",0,0,2,2)
                    +0.000011*ewgc("Cqu1R",1,1,2,2)+0.000285*ewgc("Cqu1R",2,2,0,0)+0.000004*ewgc("Cqu1R",2,2,1,1)
                    +-0.000137*ewgc("Cud1R",2,2,0,0)+-0.000016*ewgc("Cud1R",2,2,1,1)+0.000755*ewgc("CuuR",0,0,2,2)
                    +0.000020*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_tta_bin_180_300/sigma_tta_bin_180_300_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct binning range for sigma_tta_diff_LO_ATLAS_emu.\n");
    }
    
    

}



sigma_tta_diff_LO_CMS_semileptonic_210701508::sigma_tta_diff_LO_CMS_semileptonic_210701508(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{

    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tta_bin_20_35_CMS_semileptonic"
            << "SM_sigma_tta_bin_35_50_CMS_semileptonic" << "SM_sigma_tta_bin_50_80_CMS_semileptonic"
            << "SM_sigma_tta_bin_80_120_CMS_semileptonic" << "SM_sigma_tta_bin_120_160_CMS_semileptonic" 
            << "SM_sigma_tta_bin_160_200_CMS_semileptonic" << "SM_sigma_tta_bin_200_260_CMS_semileptonic"
            << "SM_sigma_tta_bin_260_320_CMS_semileptonic" << "SM_sigma_tta_bin_320_380_CMS_semileptonic");

    
}


double sigma_tta_diff_LO_CMS_semileptonic_210701508::computeThValue()
{

    b_min = getBinMin();
    b_max = getBinMax();
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    

    
    
    
    if(b_min == 20 && b_max == 35){
        
        double SM_sigma_tta_bin_20_35 = SM.getOptionalParameter("SM_sigma_tta_bin_20_35_CMS_semileptonic");
        double sigma_tta_bin_20_35_madgraph = 0.583675; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
                    
            double total = SM_sigma_tta_bin_20_35 +(-0.034106*ewgc("CG")+0.004794*ewgc("Cqd8R",2,2,0,0)
                    +0.000605*ewgc("Cqd8R",2,2,1,1)+0.082493*ewgc("Cqq1R",0,2,2,0)
                    +0.002529*ewgc("Cqq1R",1,2,2,1)+0.146282*ewgc("Cqq3R",0,2,2,0)+0.010782*ewgc("Cqq3R",1,2,2,1)
                    +0.036042*ewgc("Cqu8R",0,0,2,2)+0.001578*ewgc("Cqu8R",1,1,2,2)+0.031190*ewgc("Cqu8R",2,2,0,0)
                    +0.000906*ewgc("Cqu8R",2,2,1,1)+-0.007204*ewgc("CuBR",2,2)+0.008516*ewgc("Cud8R",2,2,0,0)
                    +0.001057*ewgc("Cud8R",2,2,1,1)+-0.154829*ewgc("CuGR",2,2)
                    +0.080748*ewgc("CuuR",0,2,2,0)+0.002498*ewgc("CuuR",1,2,2,1)+-0.005006*ewgc("CuWR",2,2)
                    +-0.061811*ewgc("CHD")+-0.035571*ewgc("CHl3R",0,0)+-0.035571*ewgc("CHl3R",1,1)+-0.132541*ewgc("CHWB")
                    +0.035506*ewgc("CllR",0,1,1,0)+-0.000081*ewgc("Cqd1R",2,2,0,0)+-0.000027*ewgc("Cqd1R",2,2,1,1)
                    +0.006399*ewgc("Cqq1R",0,0,2,2)+-0.000210*ewgc("Cqq1R",1,1,2,2)+0.016394*ewgc("Cqq3R",0,0,2,2)
                    +0.000986*ewgc("Cqq3R",1,1,2,2)+0.002978*ewgc("Cqu1R",0,0,2,2)+0.000069*ewgc("Cqu1R",1,1,2,2)
                    +0.002181*ewgc("Cqu1R",2,2,0,0)+0.000072*ewgc("Cqu1R",2,2,1,1)+-0.000689*ewgc("Cud1R",2,2,0,0)
                    +-0.000057*ewgc("Cud1R",2,2,1,1)+0.005823*ewgc("CuuR",0,0,2,2)+0.000189*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_tta_bin_20_35/sigma_tta_bin_20_35_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
    } else if(b_min == 35 && b_max == 50){
        
    
        double SM_sigma_tta_bin_35_50 = SM.getOptionalParameter("SM_sigma_tta_bin_35_50_CMS_semileptonic");
        double sigma_tta_bin_35_50_madgraph = 0.305601; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tta_bin_35_50 +(-0.018555*ewgc("CG")+0.002547*ewgc("Cqd8R",2,2,0,0)
                    +0.000295*ewgc("Cqd8R",2,2,1,1)+0.040351*ewgc("Cqq1R",0,2,2,0)
                    +0.001256*ewgc("Cqq1R",1,2,2,1)+0.075386*ewgc("Cqq3R",0,2,2,0)+0.005401*ewgc("Cqq3R",1,2,2,1)
                    +0.018788*ewgc("Cqu8R",0,0,2,2)+0.000807*ewgc("Cqu8R",1,1,2,2)+0.016262*ewgc("Cqu8R",2,2,0,0)
                    +0.000470*ewgc("Cqu8R",2,2,1,1)+-0.006535*ewgc("CuBR",2,2)+0.004622*ewgc("Cud8R",2,2,0,0)
                    +0.000557*ewgc("Cud8R",2,2,1,1)-0.080445*ewgc("CuGR",2,2)
                    +0.039496*ewgc("CuuR",0,2,2,0)+0.001217*ewgc("CuuR",1,2,2,1)+-0.004069*ewgc("CuWR",2,2)
                    +-0.032332*ewgc("CHD")+-0.018575*ewgc("CHl3R",0,0)+-0.018575*ewgc("CHl3R",1,1)+-0.069327*ewgc("CHWB")
                    +0.018597*ewgc("CllR",0,1,1,0)+-0.000057*ewgc("Cqd1R",2,2,0,0)+-0.000002*ewgc("Cqd1R",2,2,1,1)
                    +0.003004*ewgc("Cqq1R",0,0,2,2)+-0.000132*ewgc("Cqq1R",1,1,2,2)+0.008326*ewgc("Cqq3R",0,0,2,2)
                    +0.000523*ewgc("Cqq3R",1,1,2,2)+0.001386*ewgc("Cqu1R",0,0,2,2)+0.000057*ewgc("Cqu1R",1,1,2,2)
                    +0.001019*ewgc("Cqu1R",2,2,0,0)+0.000055*ewgc("Cqu1R",2,2,1,1)+-0.000289*ewgc("Cud1R",2,2,0,0)
                    +-0.000021*ewgc("Cud1R",2,2,1,1)+0.002900*ewgc("CuuR",0,0,2,2)+0.000133*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_tta_bin_35_50/sigma_tta_bin_35_50_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
    } else if(b_min == 50 && b_max == 80){
        
        double SM_sigma_tta_bin_50_80 = SM.getOptionalParameter("SM_sigma_tta_bin_50_80_CMS_semileptonic");
        double sigma_tta_bin_50_80_madgraph = 0.319505; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
                    
            double total = SM_sigma_tta_bin_50_80 +(-0.020248*ewgc("CG")+0.002577*ewgc("Cqd8R",2,2,0,0)
                    +0.000287*ewgc("Cqd8R",2,2,1,1)+0.039600*ewgc("Cqq1R",0,2,2,0)
                    +0.001092*ewgc("Cqq1R",1,2,2,1)+0.077011*ewgc("Cqq3R",0,2,2,0)
                    +0.005369*ewgc("Cqq3R",1,2,2,1)+0.019650*ewgc("Cqu8R",0,0,2,2)
                    +0.000769*ewgc("Cqu8R",1,1,2,2)+0.017049*ewgc("Cqu8R",2,2,0,0)
                    +0.000467*ewgc("Cqu8R",2,2,1,1)+-0.010718*ewgc("CuBR",2,2)
                    +0.004924*ewgc("Cud8R",2,2,0,0)+0.000555*ewgc("Cud8R",2,2,1,1)
                    -0.082799*ewgc("CuGR",2,2)+0.038642*ewgc("CuuR",0,2,2,0)
                    +0.001064*ewgc("CuuR",1,2,2,1)+-0.006289*ewgc("CuWR",2,2)
                    +-0.033790*ewgc("CHD")+-0.019417*ewgc("CHl3R",0,0)+-0.019417*ewgc("CHl3R",1,1)
                    +-0.072465*ewgc("CHWB")+0.019405*ewgc("CllR",0,1,1,0)+-0.000134*ewgc("Cqd1R",2,2,0,0)
                    +-0.000033*ewgc("Cqd1R",2,2,1,1)+0.002483*ewgc("Cqq1R",0,0,2,2)+-0.000227*ewgc("Cqq1R",1,1,2,2)
                    +0.008481*ewgc("Cqq3R",0,0,2,2)+0.000442*ewgc("Cqq3R",1,1,2,2)+0.001386*ewgc("Cqu1R",0,0,2,2)
                    +0.000011*ewgc("Cqu1R",1,1,2,2)+0.000965*ewgc("Cqu1R",2,2,0,0)+0.000015*ewgc("Cqu1R",2,2,1,1)
                    +-0.000494*ewgc("Cud1R",2,2,0,0)+-0.000044*ewgc("Cud1R",2,2,1,1)+0.002770*ewgc("CuuR",0,0,2,2)
                    +0.000056*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_tta_bin_50_80/sigma_tta_bin_50_80_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
            
        }
    } else if(b_min == 80 && b_max == 120){
        
    
        double SM_sigma_tta_bin_80_120 = SM.getOptionalParameter("SM_sigma_tta_bin_80_120_CMS_semileptonic");
        double sigma_tta_bin_80_120_madgraph = 0.198070; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tta_bin_80_120 +(-0.013527*ewgc("CG")+0.001628*ewgc("Cqd8R",2,2,0,0)
                    +0.000164*ewgc("Cqd8R",2,2,1,1)+0.023160*ewgc("Cqq1R",0,2,2,0)
                    +0.000587*ewgc("Cqq1R",1,2,2,1)+0.047833*ewgc("Cqq3R",0,2,2,0)
                    +0.003240*ewgc("Cqq3R",1,2,2,1)+0.012405*ewgc("Cqu8R",0,0,2,2)
                    +0.000452*ewgc("Cqu8R",1,1,2,2)+0.010800*ewgc("Cqu8R",2,2,0,0)
                    +0.000273*ewgc("Cqu8R",2,2,1,1)+-0.009514*ewgc("CuBR",2,2)
                    +0.003221*ewgc("Cud8R",2,2,0,0)+0.000351*ewgc("Cud8R",2,2,1,1)
                    -0.050365*ewgc("CuGR",2,2)+0.022607*ewgc("CuuR",0,2,2,0)
                    +0.000576*ewgc("CuuR",1,2,2,1)+-0.005448*ewgc("CuWR",2,2)
                    +-0.020954*ewgc("CHD")+-0.012043*ewgc("CHl3R",0,0)+-0.012043*ewgc("CHl3R",1,1)
                    +-0.044910*ewgc("CHWB")+0.012023*ewgc("CllR",0,1,1,0)+-0.000046*ewgc("Cqd1R",2,2,0,0)
                    +-0.000006*ewgc("Cqd1R",2,2,1,1)+0.001338*ewgc("Cqq1R",0,0,2,2)+-0.000141*ewgc("Cqq1R",1,1,2,2)
                    +0.005279*ewgc("Cqq3R",0,0,2,2)+0.000310*ewgc("Cqq3R",1,1,2,2)+0.000893*ewgc("Cqu1R",0,0,2,2)
                    +0.000039*ewgc("Cqu1R",1,1,2,2)+0.000660*ewgc("Cqu1R",2,2,0,0)+0.000008*ewgc("Cqu1R",2,2,1,1)
                    +-0.000280*ewgc("Cud1R",2,2,0,0)+-0.000026*ewgc("Cud1R",2,2,1,1)+0.001709*ewgc("CuuR",0,0,2,2)
                    +0.000051*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_tta_bin_80_120/sigma_tta_bin_80_120_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
    } else if(b_min == 120 && b_max == 160){
        
    
        double SM_sigma_tta_bin_120_160 = SM.getOptionalParameter("SM_sigma_tta_bin_120_160_CMS_semileptonic");
        double sigma_tta_bin_120_160_madgraph = 0.097725; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tta_bin_120_160 +(-0.007834*ewgc("CG")+0.000808*ewgc("Cqd8R",2,2,0,0)
                    +0.000071*ewgc("Cqd8R",2,2,1,1)+0.011182*ewgc("Cqq1R",0,2,2,0)
                    +0.000245*ewgc("Cqq1R",1,2,2,1)+0.024260*ewgc("Cqq3R",0,2,2,0)
                    +0.001549*ewgc("Cqq3R",1,2,2,1)+0.006409*ewgc("Cqu8R",0,0,2,2)
                    +0.000204*ewgc("Cqu8R",1,1,2,2)+0.005568*ewgc("Cqu8R",2,2,0,0)
                    +0.000122*ewgc("Cqu8R",2,2,1,1)+-0.005632*ewgc("CuBR",2,2)
                    +0.001708*ewgc("Cud8R",2,2,0,0)+0.000162*ewgc("Cud8R",2,2,1,1)
                    -0.024540*ewgc("CuGR",2,2)+0.010926*ewgc("CuuR",0,2,2,0)
                    +0.000235*ewgc("CuuR",1,2,2,1)+-0.003185*ewgc("CuWR",2,2)
                    +-0.010353*ewgc("CHD")+-0.005952*ewgc("CHl3R",0,0)+-0.005952*ewgc("CHl3R",1,1)+-0.022160*ewgc("CHWB")
                    +0.005918*ewgc("CllR",0,1,1,0)+-0.000019*ewgc("Cqd1R",2,2,0,0)+-0.000016*ewgc("Cqd1R",2,2,1,1)
                    +0.000566*ewgc("Cqq1R",0,0,2,2)+-0.000079*ewgc("Cqq1R",1,1,2,2)+0.002746*ewgc("Cqq3R",0,0,2,2)
                    +0.000133*ewgc("Cqq3R",1,1,2,2)+0.000456*ewgc("Cqu1R",0,0,2,2)+-0.000005*ewgc("Cqu1R",1,1,2,2)
                    +0.000346*ewgc("Cqu1R",2,2,0,0)+-0.000009*ewgc("Cqu1R",2,2,1,1)+-0.000154*ewgc("Cud1R",2,2,0,0)
                    +-0.000030*ewgc("Cud1R",2,2,1,1)+0.000819*ewgc("CuuR",0,0,2,2)+0.000004*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_tta_bin_120_160/sigma_tta_bin_120_160_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
    } else if(b_min == 160 && b_max == 200){
        
    
        double SM_sigma_tta_bin_160_200 = SM.getOptionalParameter("SM_sigma_tta_bin_160_200_CMS_semileptonic");
        double sigma_tta_bin_160_200_madgraph = 0.053162; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tta_bin_160_200 +(-0.005116*ewgc("CG")+0.000482*ewgc("Cqd8R",2,2,0,0)
                    +0.000047*ewgc("Cqd8R",2,2,1,1)+0.006230*ewgc("Cqq1R",0,2,2,0)
                    +0.000137*ewgc("Cqq1R",1,2,2,1)+0.013887*ewgc("Cqq3R",0,2,2,0)
                    +0.000853*ewgc("Cqq3R",1,2,2,1)+0.003756*ewgc("Cqu8R",0,0,2,2)
                    +0.000123*ewgc("Cqu8R",1,1,2,2)+0.003272*ewgc("Cqu8R",2,2,0,0)
                    +0.000076*ewgc("Cqu8R",2,2,1,1)+-0.003224*ewgc("CuBR",2,2)
                    +0.001020*ewgc("Cud8R",2,2,0,0)+0.000101*ewgc("Cud8R",2,2,1,1)
                    -0.013391*ewgc("CuGR",2,2)+0.006063*ewgc("CuuR",0,2,2,0)
                    +0.000133*ewgc("CuuR",1,2,2,1)+-0.001792*ewgc("CuWR",2,2)
                    +-0.005617*ewgc("CHD")+-0.003223*ewgc("CHl3R",0,0)+-0.003223*ewgc("CHl3R",1,1)
                    +-0.012030*ewgc("CHWB")+0.003235*ewgc("CllR",0,1,1,0)+-0.000003*ewgc("Cqd1R",2,2,0,0)
                    +0.000004*ewgc("Cqd1R",2,2,1,1)+0.000283*ewgc("Cqq1R",0,0,2,2)+-0.000035*ewgc("Cqq1R",1,1,2,2)
                    +0.001579*ewgc("Cqq3R",0,0,2,2)+0.000085*ewgc("Cqq3R",1,1,2,2)+0.000238*ewgc("Cqu1R",0,0,2,2)
                    +0.000009*ewgc("Cqu1R",1,1,2,2)+0.000208*ewgc("Cqu1R",2,2,0,0)+0.000010*ewgc("Cqu1R",2,2,1,1)
                    +-0.000083*ewgc("Cud1R",2,2,0,0)+-0.000005*ewgc("Cud1R",2,2,1,1)+0.000471*ewgc("CuuR",0,0,2,2)
                    +0.000014*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_tta_bin_160_200/sigma_tta_bin_160_200_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
            
        }
    } else if(b_min == 200 && b_max == 260){
        
    
        double SM_sigma_tta_bin_200_260 = SM.getOptionalParameter("SM_sigma_tta_bin_200_260_CMS_semileptonic");
        double sigma_tta_bin_200_260_madgraph = 0.040676; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tta_bin_200_260 +(-0.004958*ewgc("CG")+0.000418*ewgc("Cqd8R",2,2,0,0)
                    +0.000036*ewgc("Cqd8R",2,2,1,1)+0.005096*ewgc("Cqq1R",0,2,2,0)
                    +0.000101*ewgc("Cqq1R",1,2,2,1)+0.011849*ewgc("Cqq3R",0,2,2,0)
                    +0.000680*ewgc("Cqq3R",1,2,2,1)+0.003215*ewgc("Cqu8R",0,0,2,2)
                    +0.000093*ewgc("Cqu8R",1,1,2,2)+0.002804*ewgc("Cqu8R",2,2,0,0)
                    +0.000057*ewgc("Cqu8R",2,2,1,1)+-0.002405*ewgc("CuBR",2,2)
                    +0.000885*ewgc("Cud8R",2,2,0,0)+0.000078*ewgc("Cud8R",2,2,1,1)
                    -0.010460*ewgc("CuGR",2,2)+0.004970*ewgc("CuuR",0,2,2,0)
                    +0.000098*ewgc("CuuR",1,2,2,1)+-0.001356*ewgc("CuWR",2,2)
                    +-0.004304*ewgc("CHD")+-0.002471*ewgc("CHl3R",0,0)+-0.002471*ewgc("CHl3R",1,1)
                    +-0.009223*ewgc("CHWB")+0.002470*ewgc("CllR",0,1,1,0)+-0.000011*ewgc("Cqd1R",2,2,0,0)
                    +0.000001*ewgc("Cqd1R",2,2,1,1)+0.000251*ewgc("Cqq1R",0,0,2,2)+-0.000039*ewgc("Cqq1R",1,1,2,2)
                    +0.001359*ewgc("Cqq3R",0,0,2,2)+0.000063*ewgc("Cqq3R",1,1,2,2)+0.000204*ewgc("Cqu1R",0,0,2,2)
                    +0.000004*ewgc("Cqu1R",1,1,2,2)+0.000151*ewgc("Cqu1R",2,2,0,0)+0.000004*ewgc("Cqu1R",2,2,1,1)
                    +-0.000074*ewgc("Cud1R",2,2,0,0)+-0.000005*ewgc("Cud1R",2,2,1,1)+0.000416*ewgc("CuuR",0,0,2,2)
                    +0.000007*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_tta_bin_200_260/sigma_tta_bin_200_260_madgraph);
            
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
    } else if(b_min == 260 && b_max == 320){
        
    
        double SM_sigma_tta_bin_260_320 = SM.getOptionalParameter("SM_sigma_tta_bin_260_320_CMS_semileptonic");
        double sigma_tta_bin_260_320_madgraph = 0.018994; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tta_bin_260_320 +(-0.003098*ewgc("CG")+0.000227*ewgc("Cqd8R",2,2,0,0)
                    +0.000013*ewgc("Cqd8R",2,2,1,1)+0.002712*ewgc("Cqq1R",0,2,2,0)
                    +0.000044*ewgc("Cqq1R",1,2,2,1)+0.006445*ewgc("Cqq3R",0,2,2,0)
                    +0.000337*ewgc("Cqq3R",1,2,2,1)+0.001768*ewgc("Cqu8R",0,0,2,2)
                    +0.000042*ewgc("Cqu8R",1,1,2,2)+0.001534*ewgc("Cqu8R",2,2,0,0)
                    +0.000025*ewgc("Cqu8R",2,2,1,1)+-0.001038*ewgc("CuBR",2,2)
                    +0.000492*ewgc("Cud8R",2,2,0,0)+0.000035*ewgc("Cud8R",2,2,1,1)
                    -0.005064*ewgc("CuGR",2,2)+0.002652*ewgc("CuuR",0,2,2,0)
                    +0.000042*ewgc("CuuR",1,2,2,1)+-0.000592*ewgc("CuWR",2,2)
                    +-0.002010*ewgc("CHD")+-0.001156*ewgc("CHl3R",0,0)+-0.001156*ewgc("CHl3R",1,1)
                    +-0.004307*ewgc("CHWB")+0.001150*ewgc("CllR",0,1,1,0)+-0.000011*ewgc("Cqd1R",2,2,0,0)
                    +-0.000005*ewgc("Cqd1R",2,2,1,1)+0.000094*ewgc("Cqq1R",0,0,2,2)+-0.000021*ewgc("Cqq1R",1,1,2,2)
                    +0.000730*ewgc("Cqq3R",0,0,2,2)+0.000028*ewgc("Cqq3R",1,1,2,2)+0.000112*ewgc("Cqu1R",0,0,2,2)
                    +-0.000003*ewgc("Cqu1R",1,1,2,2)+0.000081*ewgc("Cqu1R",2,2,0,0)+-0.000003*ewgc("Cqu1R",2,2,1,1)
                    +-0.000049*ewgc("Cud1R",2,2,0,0)+-0.000008*ewgc("Cud1R",2,2,1,1)+0.000214*ewgc("CuuR",0,0,2,2)
                    +-0.000002*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_tta_bin_260_320/sigma_tta_bin_260_320_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
            
        }
    } else if(b_min == 320 && b_max == 380){
        
    
        double SM_sigma_tta_bin_320_380 = SM.getOptionalParameter("SM_sigma_tta_bin_320_380_CMS_semileptonic");
        double sigma_tta_bin_320_380_madgraph = 0.009391; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tta_bin_320_380 +(-0.001921*ewgc("CG")+0.000136*ewgc("Cqd8R",2,2,0,0)
                    +0.000007*ewgc("Cqd8R",2,2,1,1)+0.001575*ewgc("Cqq1R",0,2,2,0)
                    +0.000023*ewgc("Cqq1R",1,2,2,1)+0.003819*ewgc("Cqq3R",0,2,2,0)
                    +0.000186*ewgc("Cqq3R",1,2,2,1)+0.001048*ewgc("Cqu8R",0,0,2,2)
                    +0.000022*ewgc("Cqu8R",1,1,2,2)+0.000910*ewgc("Cqu8R",2,2,0,0)
                    +0.000013*ewgc("Cqu8R",2,2,1,1)+-0.000451*ewgc("CuBR",2,2)
                    +0.000293*ewgc("Cud8R",2,2,0,0)+0.000019*ewgc("Cud8R",2,2,1,1)
                    -0.002605*ewgc("CuGR",2,2)+0.001534*ewgc("CuuR",0,2,2,0)
                    +0.000023*ewgc("CuuR",1,2,2,1)+-0.000262*ewgc("CuWR",2,2)
                    +-0.000996*ewgc("CHD")+-0.000574*ewgc("CHl3R",0,0)+-0.000574*ewgc("CHl3R",1,1)
                    +-0.002131*ewgc("CHWB")+0.000567*ewgc("CllR",0,1,1,0)+-0.000004*ewgc("Cqd1R",2,2,0,0)
                    +-0.000003*ewgc("Cqd1R",2,2,1,1)+0.000057*ewgc("Cqq1R",0,0,2,2)+-0.000012*ewgc("Cqq1R",1,1,2,2)
                    +0.000442*ewgc("Cqq3R",0,0,2,2)+0.000016*ewgc("Cqq3R",1,1,2,2)+0.000069*ewgc("Cqu1R",0,0,2,2)
                    +-0.000001*ewgc("Cqu1R",1,1,2,2)+0.000046*ewgc("Cqu1R",2,2,0,0)+-0.000001*ewgc("Cqu1R",2,2,1,1)
                    +-0.000029*ewgc("Cud1R",2,2,0,0)+-0.000005*ewgc("Cud1R",2,2,1,1)+0.000125*ewgc("CuuR",0,0,2,2)
                    +-0.0000005*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_tta_bin_320_380/sigma_tta_bin_320_380_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
    }  else{
            throw std::runtime_error("\nERROR: Please specify a correct binning range for sigma_tta_diff_LO_CMS_semileptonic.\n");
    }



}





sigma_tta_diff_LO_CMS_dilepton_220107301::sigma_tta_diff_LO_CMS_dilepton_220107301(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{

    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tta_bin_20_35_CMS_dilepton" << "SM_sigma_tta_bin_35_50_CMS_dilepton"
            << "SM_sigma_tta_bin_50_70_CMS_dilepton" << "SM_sigma_tta_bin_70_130_CMS_dilepton" 
            << "SM_sigma_tta_bin_130_200_CMS_dilepton" << "SM_sigma_tta_bin_200_300_CMS_dilepton");

    
}


double sigma_tta_diff_LO_CMS_dilepton_220107301::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    

    if(b_min == 20 && b_max == 35){
        
        double SM_sigma_tta_bin_20_35 = SM.getOptionalParameter("SM_sigma_tta_bin_20_35_CMS_dilepton");
        double sigma_tta_bin_20_35_madgraph = 0.583675; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            //There is no dependence on ewgc("Cqd8R",2,2,2,2) neither on ewgc("Cud8R",2,2,2,2) with the precision considered
        
            double total = SM_sigma_tta_bin_20_35 +(-0.034106*ewgc("CG")+0.004794*ewgc("Cqd8R",2,2,0,0)
                    +0.000605*ewgc("Cqd8R",2,2,1,1)+0.082493*ewgc("Cqq1R",0,2,2,0)
                    +0.002529*ewgc("Cqq1R",1,2,2,1)+0.146282*ewgc("Cqq3R",0,2,2,0)+0.010782*ewgc("Cqq3R",1,2,2,1)
                    +0.036042*ewgc("Cqu8R",0,0,2,2)+0.001578*ewgc("Cqu8R",1,1,2,2)+0.031190*ewgc("Cqu8R",2,2,0,0)
                    +0.000906*ewgc("Cqu8R",2,2,1,1)+-0.007204*ewgc("CuBR",2,2)+0.008516*ewgc("Cud8R",2,2,0,0)
                    +0.001057*ewgc("Cud8R",2,2,1,1)+-0.154829*ewgc("CuGR",2,2)
                    +0.080748*ewgc("CuuR",0,2,2,0)+0.002498*ewgc("CuuR",1,2,2,1)+-0.005006*ewgc("CuWR",2,2)
                    +-0.061811*ewgc("CHD")+-0.035571*ewgc("CHl3R",0,0)+-0.035571*ewgc("CHl3R",1,1)+-0.132541*ewgc("CHWB")
                    +0.035506*ewgc("CllR",0,1,1,0)+-0.000081*ewgc("Cqd1R",2,2,0,0)+-0.000027*ewgc("Cqd1R",2,2,1,1)
                    +0.006399*ewgc("Cqq1R",0,0,2,2)+-0.000210*ewgc("Cqq1R",1,1,2,2)+0.016394*ewgc("Cqq3R",0,0,2,2)
                    +0.000986*ewgc("Cqq3R",1,1,2,2)+0.002978*ewgc("Cqu1R",0,0,2,2)+0.000069*ewgc("Cqu1R",1,1,2,2)
                    +0.002181*ewgc("Cqu1R",2,2,0,0)+0.000072*ewgc("Cqu1R",2,2,1,1)+-0.000689*ewgc("Cud1R",2,2,0,0)
                    +-0.000057*ewgc("Cud1R",2,2,1,1)+0.005823*ewgc("CuuR",0,0,2,2)+0.000189*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_tta_bin_20_35/sigma_tta_bin_20_35_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
    } else if(b_min == 35 && b_max == 50){
        
    
        double SM_sigma_tta_bin_35_50 = SM.getOptionalParameter("SM_sigma_tta_bin_35_50_CMS_dilepton");
        double sigma_tta_bin_35_50_madgraph = 0.305601; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tta_bin_35_50 +(-0.018555*ewgc("CG")+0.002547*ewgc("Cqd8R",2,2,0,0)
                    +0.000295*ewgc("Cqd8R",2,2,1,1)+0.040351*ewgc("Cqq1R",0,2,2,0)
                    +0.001256*ewgc("Cqq1R",1,2,2,1)+0.075386*ewgc("Cqq3R",0,2,2,0)+0.005401*ewgc("Cqq3R",1,2,2,1)
                    +0.018788*ewgc("Cqu8R",0,0,2,2)+0.000807*ewgc("Cqu8R",1,1,2,2)+0.016262*ewgc("Cqu8R",2,2,0,0)
                    +0.000470*ewgc("Cqu8R",2,2,1,1)+-0.006535*ewgc("CuBR",2,2)+0.004622*ewgc("Cud8R",2,2,0,0)
                    +0.000557*ewgc("Cud8R",2,2,1,1)-0.080445*ewgc("CuGR",2,2)
                    +0.039496*ewgc("CuuR",0,2,2,0)+0.001217*ewgc("CuuR",1,2,2,1)+-0.004069*ewgc("CuWR",2,2)
                    +-0.032332*ewgc("CHD")+-0.018575*ewgc("CHl3R",0,0)+-0.018575*ewgc("CHl3R",1,1)+-0.069327*ewgc("CHWB")
                    +0.018597*ewgc("CllR",0,1,1,0)+-0.000057*ewgc("Cqd1R",2,2,0,0)+-0.000002*ewgc("Cqd1R",2,2,1,1)
                    +0.003004*ewgc("Cqq1R",0,0,2,2)+-0.000132*ewgc("Cqq1R",1,1,2,2)+0.008326*ewgc("Cqq3R",0,0,2,2)
                    +0.000523*ewgc("Cqq3R",1,1,2,2)+0.001386*ewgc("Cqu1R",0,0,2,2)+0.000057*ewgc("Cqu1R",1,1,2,2)
                    +0.001019*ewgc("Cqu1R",2,2,0,0)+0.000055*ewgc("Cqu1R",2,2,1,1)+-0.000289*ewgc("Cud1R",2,2,0,0)
                    +-0.000021*ewgc("Cud1R",2,2,1,1)+0.002900*ewgc("CuuR",0,0,2,2)+0.000133*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_tta_bin_35_50/sigma_tta_bin_35_50_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
            
        }
    } else if(b_min == 50 && b_max == 70){
        
    
        double SM_sigma_tta_bin_50_70 = SM.getOptionalParameter("SM_sigma_tta_bin_50_70_CMS_dilepton");
        double sigma_tta_bin_50_70_madgraph = 0.238500; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_tta_bin_50_70 +(-0.014859*ewgc("CG")+0.001934*ewgc("Cqd8R",2,2,0,0)
                    +0.000192*ewgc("Cqd8R",2,2,1,1)+
                    +0.029858*ewgc("Cqq1R",0,2,2,0)+0.000825*ewgc("Cqq1R",1,2,2,1)
                    +0.057705*ewgc("Cqq3R",0,2,2,0)+0.004114*ewgc("Cqq3R",1,2,2,1)
                    +0.014644*ewgc("Cqu8R",0,0,2,2)+0.000602*ewgc("Cqu8R",1,1,2,2)
                    +0.012747*ewgc("Cqu8R",2,2,0,0)+-0.007537*ewgc("CuBR",2,2)
                    +0.003692*ewgc("Cud8R",2,2,0,0)+0.000409*ewgc("Cud8R",2,2,1,1)
                    -0.061847*ewgc("CuGR",2,2)
                    +0.029209*ewgc("CuuR",0,2,2,0)+0.000813*ewgc("CuuR",1,2,2,1)
                    +-0.004460*ewgc("CuWR",2,2)
                    +-0.025253*ewgc("CHD")+-0.014510*ewgc("CHl3R",0,0)+-0.014510*ewgc("CHl3R",1,1)+-0.054091*ewgc("CHWB")
                    +0.014521*ewgc("CllR",0,1,1,0)+-0.000031*ewgc("Cqd1R",2,2,0,0)+-0.000008*ewgc("Cqd1R",2,2,1,1)
                    +0.001996*ewgc("Cqq1R",0,0,2,2)+-0.000150*ewgc("Cqq1R",1,1,2,2)+0.006394*ewgc("Cqq3R",0,0,2,2)
                    +0.000416*ewgc("Cqq3R",1,1,2,2)+0.001142*ewgc("Cqu1R",0,0,2,2)+0.000027*ewgc("Cqu1R",1,1,2,2)
                    +0.000807*ewgc("Cqu1R",2,2,0,0)+0.000037*ewgc("Cqu1R",2,2,1,1)+-0.000320*ewgc("Cud1R",2,2,0,0)
                    +-0.000011*ewgc("Cud1R",2,2,1,1)+0.002163*ewgc("CuuR",0,0,2,2)+0.000065*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_tta_bin_50_70/sigma_tta_bin_50_70_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
    } else if(b_min == 70 && b_max == 130){
        
    
        double SM_sigma_tta_bin_70_130 = SM.getOptionalParameter("SM_sigma_tta_bin_70_130_CMS_dilepton");
        double sigma_tta_bin_70_100_madgraph = 0.198831; //pb
        double sigma_tta_bin_100_130_madgraph =0.110866; //pb
        double sigma_tta_bin_70_130_madgraph =sigma_tta_bin_70_100_madgraph+sigma_tta_bin_100_130_madgraph; //pb
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double sigma_tta_bin_70_100_NP = (-0.013285*ewgc("CG")+0.001630*ewgc("Cqd8R",2,2,0,0)
                    +0.000250*ewgc("Cqd8R",2,2,1,1)
                    +0.023767*ewgc("Cqq1R",0,2,2,0)+0.000676*ewgc("Cqq1R",1,2,2,1)
                    +0.047671*ewgc("Cqq3R",0,2,2,0)+0.003336*ewgc("Cqq3R",1,2,2,1)
                    +0.012379*ewgc("Cqu8R",0,0,2,2)+0.000535*ewgc("Cqu8R",1,1,2,2)
                    +0.010765*ewgc("Cqu8R",2,2,0,0)+0.000350*ewgc("Cqu8R",2,2,1,1)
                    -0.008505*ewgc("CuBR",2,2)+0.003205*ewgc("Cud8R",2,2,0,0)
                    +0.000431*ewgc("Cud8R",2,2,1,1)
                    -0.050838*ewgc("CuGR",2,2)+0.023185*ewgc("CuuR",0,2,2,0)
                    +0.000655*ewgc("CuuR",1,2,2,1)-0.004930*ewgc("CuWR",2,2)
                    +-0.020982*ewgc("CHD")+-0.012029*ewgc("CHl3R",0,0)+-0.012029*ewgc("CHl3R",1,1)+-0.045039*ewgc("CHWB")
                    +0.012146*ewgc("CllR",0,1,1,0)+-0.000043*ewgc("Cqd1R",2,2,0,0)+0.000040*ewgc("Cqd1R",2,2,1,1)
                    +0.001365*ewgc("Cqq1R",0,0,2,2)+-0.000101*ewgc("Cqq1R",1,1,2,2)+0.005271*ewgc("Cqq3R",0,0,2,2)
                    +0.000313*ewgc("Cqq3R",1,1,2,2)+0.000917*ewgc("Cqu1R",0,0,2,2)+0.000089*ewgc("Cqu1R",1,1,2,2)
                    +0.000787*ewgc("Cqu1R",2,2,0,0)+0.000041*ewgc("Cqu1R",2,2,1,1)+-0.000264*ewgc("Cud1R",2,2,0,0)
                    +0.000019*ewgc("Cud1R",2,2,1,1)+0.001650*ewgc("CuuR",0,0,2,2)+0.000103*ewgc("CuuR",1,1,2,2));
             
            
            double sigma_tta_bin_100_130_NP = (-0.008092*ewgc("CG")+0.000916*ewgc("Cqd8R",2,2,0,0)+0.000088*ewgc("Cqd8R",2,2,1,1)
                    +0.012744*ewgc("Cqq1R",0,2,2,0)+0.000316*ewgc("Cqq1R",1,2,2,1)
                    +0.026925*ewgc("Cqq3R",0,2,2,0)+0.001792*ewgc("Cqq3R",1,2,2,1)+0.007060*ewgc("Cqu8R",0,0,2,2)
                    +0.000247*ewgc("Cqu8R",1,1,2,2)+0.006144*ewgc("Cqu8R",2,2,0,0)+0.000162*ewgc("Cqu8R",2,2,1,1)
                    +-0.005920*ewgc("CuBR",2,2)+0.001849*ewgc("Cud8R",2,2,0,0)+0.000200*ewgc("Cud8R",2,2,1,1)
                    +-0.027969*ewgc("CuGR",2,2)+0.012428*ewgc("CuuR",0,2,2,0)
                    +0.000308*ewgc("CuuR",1,2,2,1)+-0.003364*ewgc("CuWR",2,2)
                    +-0.011717*ewgc("CHD")+-0.006725*ewgc("CHl3R",0,0)+-0.006725*ewgc("CHl3R",1,1)+-0.025133*ewgc("CHWB")
                    +0.006730*ewgc("CllR",0,1,1,0)+-0.000029*ewgc("Cqd1R",2,2,0,0)+-0.000010*ewgc("Cqd1R",2,2,1,1)
                    +0.000698*ewgc("Cqq1R",0,0,2,2)+-0.000082*ewgc("Cqq1R",1,1,2,2)+0.002982*ewgc("Cqq3R",0,0,2,2)
                    +0.000170*ewgc("Cqq3R",1,1,2,2)+0.000509*ewgc("Cqu1R",0,0,2,2)+0.000010*ewgc("Cqu1R",1,1,2,2)
                    +0.000380*ewgc("Cqu1R",2,2,0,0)+0.000006*ewgc("Cqu1R",2,2,1,1)+-0.000174*ewgc("Cud1R",2,2,0,0)
                    +-0.000022*ewgc("Cud1R",2,2,1,1)+0.000917*ewgc("CuuR",0,0,2,2)+0.000026*ewgc("CuuR",1,1,2,2));
            
            double total = SM_sigma_tta_bin_70_130 + (sigma_tta_bin_70_100_NP + sigma_tta_bin_100_130_NP)*(SM_sigma_tta_bin_70_130/sigma_tta_bin_70_130_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
    } else if(b_min == 130 && b_max == 200){
        
    
        double SM_sigma_tta_bin_130_200 = SM.getOptionalParameter("SM_sigma_tta_bin_130_200_CMS_dilepton");
        double sigma_tta_bin_130_165_madgraph = 0.075543; //pb
        double sigma_tta_bin_165_200_madgraph = 0.044687; //pb
        double sigma_tta_bin_130_200_madgraph = sigma_tta_bin_130_165_madgraph+sigma_tta_bin_165_200_madgraph;
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double sigma_tta_bin_130_165_NP = (-0.006183*ewgc("CG")+0.000634*ewgc("Cqd8R",2,2,0,0)+0.000043*ewgc("Cqd8R",2,2,1,1)
                    +0.008637*ewgc("Cqq1R",0,2,2,0)+0.000173*ewgc("Cqq1R",1,2,2,1)
                    +0.018848*ewgc("Cqq3R",0,2,2,0)+0.001189*ewgc("Cqq3R",1,2,2,1)+0.005006*ewgc("Cqu8R",0,0,2,2)
                    +0.000145*ewgc("Cqu8R",1,1,2,2)+0.004346*ewgc("Cqu8R",2,2,0,0)+0.000081*ewgc("Cqu8R",2,2,1,1)
                    +-0.004429*ewgc("CuBR",2,2)+0.001330*ewgc("Cud8R",2,2,0,0)+0.000119*ewgc("Cud8R",2,2,1,1)
                    +-0.018971*ewgc("CuGR",2,2)+0.008402*ewgc("CuuR",0,2,2,0)
                    +0.000171*ewgc("CuuR",1,2,2,1)+-0.002499*ewgc("CuWR",2,2)
                    +-0.008004*ewgc("CHD")+-0.004609*ewgc("CHl3R",0,0)+-0.004609*ewgc("CHl3R",1,1)+-0.017140*ewgc("CHWB")
                    +0.004568*ewgc("CllR",0,1,1,0)+-0.000040*ewgc("Cqd1R",2,2,0,0)+-0.000023*ewgc("Cqd1R",2,2,1,1)
                    +0.000407*ewgc("Cqq1R",0,0,2,2)+-0.000077*ewgc("Cqq1R",1,1,2,2)+0.002110*ewgc("Cqq3R",0,0,2,2)
                    +0.000087*ewgc("Cqq3R",1,1,2,2)+0.000341*ewgc("Cqu1R",0,0,2,2)+-0.000013*ewgc("Cqu1R",1,1,2,2)
                    +0.000255*ewgc("Cqu1R",2,2,0,0)+-0.000015*ewgc("Cqu1R",2,2,1,1)+-0.000134*ewgc("Cud1R",2,2,0,0)
                    +-0.000034*ewgc("Cud1R",2,2,1,1)+0.000617*ewgc("CuuR",0,0,2,2)+-0.000013*ewgc("CuuR",1,1,2,2));
            
            double sigma_tta_bin_165_200_NP = (-0.004285*ewgc("CG")+0.000401*ewgc("Cqd8R",2,2,0,0)+0.000029*ewgc("Cqd8R",2,2,1,1)
                    +0.005252*ewgc("Cqq1R",0,2,2,0)+0.000104*ewgc("Cqq1R",1,2,2,1)
                    +0.011756*ewgc("Cqq3R",0,2,2,0)+0.000720*ewgc("Cqq3R",1,2,2,1)+0.003157*ewgc("Cqu8R",0,0,2,2)
                    +0.000089*ewgc("Cqu8R",1,1,2,2)+0.002752*ewgc("Cqu8R",2,2,0,0)+0.000051*ewgc("Cqu8R",2,2,1,1)
                    +-0.002716*ewgc("CuBR",2,2)+0.000855*ewgc("Cud8R",2,2,0,0)+0.000072*ewgc("Cud8R",2,2,1,1)
                    +-0.011315*ewgc("CuGR",2,2)+0.005106*ewgc("CuuR",0,2,2,0)
                    +0.000098*ewgc("CuuR",1,2,2,1)+-0.001515*ewgc("CuWR",2,2)
                    +-0.004733*ewgc("CHD")+-0.002719*ewgc("CHl3R",0,0)+-0.002719*ewgc("CHl3R",1,1)+-0.010136*ewgc("CHWB")
                    +0.002705*ewgc("CllR",0,1,1,0)+-0.000020*ewgc("Cqd1R",2,2,0,0)+-0.000010*ewgc("Cqd1R",2,2,1,1)
                    +0.000221*ewgc("Cqq1R",0,0,2,2)+-0.000045*ewgc("Cqq1R",1,1,2,2)+0.001322*ewgc("Cqq3R",0,0,2,2)
                    +0.000061*ewgc("Cqq3R",1,1,2,2)+0.000204*ewgc("Cqu1R",0,0,2,2)+-0.000003*ewgc("Cqu1R",1,1,2,2)
                    +0.000161*ewgc("Cqu1R",2,2,0,0)+-0.000007*ewgc("Cqu1R",2,2,1,1)+-0.000081*ewgc("Cud1R",2,2,0,0)
                    +-0.000014*ewgc("Cud1R",2,2,1,1)+0.000411*ewgc("CuuR",0,0,2,2)+-0.000002*ewgc("CuuR",1,1,2,2));
            
            double total = SM_sigma_tta_bin_130_200 + (sigma_tta_bin_130_165_NP+sigma_tta_bin_165_200_NP)*(SM_sigma_tta_bin_130_200/sigma_tta_bin_130_200_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
    } else if(b_min == 200 && b_max == 300){
        
    
        double SM_sigma_tta_bin_200_300 = SM.getOptionalParameter("SM_sigma_tta_bin_200_300_CMS_dilepton");
        double sigma_tta_bin_200_250_madgraph = 0.035911; //pb
        double sigma_tta_bin_250_300_madgraph = 0.018900; //pb
        double sigma_tta_bin_200_300_madgraph = sigma_tta_bin_200_250_madgraph+sigma_tta_bin_250_300_madgraph; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double sigma_tta_bin_200_250_NP = (-0.004232*ewgc("CG")+0.000362*ewgc("Cqd8R",2,2,0,0)+0.000028*ewgc("Cqd8R",2,2,1,1)
                    +0.004454*ewgc("Cqq1R",0,2,2,0)+0.000087*ewgc("Cqq1R",1,2,2,1)
                    +0.010395*ewgc("Cqq3R",0,2,2,0)+0.000590*ewgc("Cqq3R",1,2,2,1)+0.002802*ewgc("Cqu8R",0,0,2,2)
                    +0.000078*ewgc("Cqu8R",1,1,2,2)+0.002436*ewgc("Cqu8R",2,2,0,0)+0.000047*ewgc("Cqu8R",2,2,1,1)
                    -0.002135*ewgc("CuBR",2,2)+0.000773*ewgc("Cud8R",2,2,0,0)+0.000064*ewgc("Cud8R",2,2,1,1)
                    -0.009208*ewgc("CuGR",2,2)+0.004354*ewgc("CuuR",0,2,2,0)
                    +0.000086*ewgc("CuuR",1,2,2,1)+-0.001220*ewgc("CuWR",2,2)
                    +-0.003802*ewgc("CHD")+-0.002185*ewgc("CHl3R",0,0)+-0.002185*ewgc("CHl3R",1,1)+-0.008144*ewgc("CHWB")
                    +0.002181*ewgc("CllR",0,1,1,0)+-0.000016*ewgc("Cqd1R",2,2,0,0)+-0.000004*ewgc("Cqd1R",2,2,1,1)
                    +0.000176*ewgc("Cqq1R",0,0,2,2)+-0.000031*ewgc("Cqq1R",1,1,2,2)+0.001156*ewgc("Cqq3R",0,0,2,2)
                    +0.000055*ewgc("Cqq3R",1,1,2,2)+0.000179*ewgc("Cqu1R",0,0,2,2)+0.000003*ewgc("Cqu1R",1,1,2,2)
                    +0.000131*ewgc("Cqu1R",2,2,0,0)+0.000000*ewgc("Cqu1R",2,2,1,1)+-0.000073*ewgc("Cud1R",2,2,0,0)
                    +-0.000011*ewgc("Cud1R",2,2,1,1)+0.000320*ewgc("CuuR",0,0,2,2)+0.000005*ewgc("CuuR",1,1,2,2));
            
            double sigma_tta_bin_250_300_NP = (-0.002904*ewgc("CG")+0.000216*ewgc("Cqd8R",2,2,0,0)+0.000014*ewgc("Cqd8R",2,2,1,1)
                    +0.002603*ewgc("Cqq1R",0,2,2,0)+0.000045*ewgc("Cqq1R",1,2,2,1)
                    +0.006131*ewgc("Cqq3R",0,2,2,0)+0.000331*ewgc("Cqq3R",1,2,2,1)+0.001689*ewgc("Cqu8R",0,0,2,2)
                    +0.000042*ewgc("Cqu8R",1,1,2,2)+0.001470*ewgc("Cqu8R",2,2,0,0)+0.000026*ewgc("Cqu8R",2,2,1,1)
                    -0.001064*ewgc("CuBR",2,2)+0.000469*ewgc("Cud8R",2,2,0,0)+0.000035*ewgc("Cud8R",2,2,1,1)
                    +-0.004994*ewgc("CuGR",2,2)+0.002538*ewgc("CuuR",0,2,2,0)
                    +0.000044*ewgc("CuuR",1,2,2,1)+-0.000597*ewgc("CuWR",2,2)
                    +-0.001999*ewgc("CHD")+-0.001149*ewgc("CHl3R",0,0)+-0.001149*ewgc("CHl3R",1,1)+-0.004281*ewgc("CHWB")
                    +0.001146*ewgc("CllR",0,1,1,0)+-0.000008*ewgc("Cqd1R",2,2,0,0)+-0.000003*ewgc("Cqd1R",2,2,1,1)
                    +0.000087*ewgc("Cqq1R",0,0,2,2)+-0.000021*ewgc("Cqq1R",1,1,2,2)+0.000688*ewgc("Cqq3R",0,0,2,2)
                    +0.000029*ewgc("Cqq3R",1,1,2,2)+0.000105*ewgc("Cqu1R",0,0,2,2)+-0.000003*ewgc("Cqu1R",1,1,2,2)
                    +0.000073*ewgc("Cqu1R",2,2,0,0)+-0.000001*ewgc("Cqu1R",2,2,1,1)+-0.000043*ewgc("Cud1R",2,2,0,0)
                    +-0.000008*ewgc("Cud1R",2,2,1,1)+0.000197*ewgc("CuuR",0,0,2,2)+0.000001*ewgc("CuuR",1,1,2,2));
            
            double total = SM_sigma_tta_bin_200_300 + (sigma_tta_bin_200_250_NP + sigma_tta_bin_250_300_NP)*(SM_sigma_tta_bin_200_300/sigma_tta_bin_200_300_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct binning range for sigma_tta_diff_LO_CMS_dilepton.\n");
    }

}


//// ttz differential cross section //////


sigma_ttz_diff_LO_CMS_190711270::sigma_ttz_diff_LO_CMS_190711270(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{

    setParametersForObservable(make_vector<std::string>() << "SM_sigma_ttz_bin_0_75_CMS_190711270"
            << "SM_sigma_ttz_bin_75_150_CMS_190711270" << "SM_sigma_ttz_bin_150_250_CMS_190711270"
            << "SM_sigma_ttz_bin_250_500_CMS_190711270");

    
}


double sigma_ttz_diff_LO_CMS_190711270::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    bool   flag_Quadratic= false; //Needs to be properly defined
   
    if(b_min == 0 && b_max == 75){
        
        double SM_sigma_ttz_bin_0_75 = SM.getOptionalParameter("SM_sigma_ttz_bin_0_75_CMS_190711270");
        double sigma_ttz_bin_0_75_madgraph = 0.177589; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_ttz_bin_0_75 + (-0.031991*ewgc("CG")+-0.016710*ewgc("CHq1R",2,2)+0.016573*ewgc("CHq3R",2,2)
                    +0.009051*ewgc("CHuR",2,2)+0.003299*ewgc("Cqd8R",2,2,0,0)+0.000362*ewgc("Cqd8R",2,2,1,1)
                    +0.024362*ewgc("Cqq1R",0,2,2,0)+0.000624*ewgc("Cqq1R",1,2,2,1)+0.119583*ewgc("Cqq3R",0,2,2,0)
                    +0.012655*ewgc("Cqq3R",1,2,2,1)+0.015687*ewgc("Cqu8R",0,0,2,2)+0.001422*ewgc("Cqu8R",1,1,2,2)
                    +0.002967*ewgc("Cqu8R",2,2,0,0)+0.000052*ewgc("Cqu8R",2,2,1,1)+-0.000318*ewgc("CuBR",2,2)
                    +0.001835*ewgc("Cud8R",2,2,0,0)+0.000209*ewgc("Cud8R",2,2,1,1)-0.052919*ewgc("CuGR",2,2)
                    +0.006972*ewgc("CuuR",0,2,2,0)+0.000145*ewgc("CuuR",1,2,2,1)+-0.000317*ewgc("CuWR",2,2)
                    +0.000764*ewgc("CHD")+-0.010859*ewgc("CHl3R",0,0)+-0.010859*ewgc("CHl3R",1,1)+-0.000041*ewgc("CHW")
                    +0.006519*ewgc("CHWB")+0.010767*ewgc("CllR",0,1,1,0)+-0.000148*ewgc("Cqd1R",2,2,0,0)+-0.000064*ewgc("Cqd1R",2,2,1,1)
                    +-0.002423*ewgc("Cqq1R",0,0,2,2)+-0.000764*ewgc("Cqq1R",1,1,2,2)+0.011282*ewgc("Cqq3R",0,0,2,2)
                    +0.000883*ewgc("Cqq3R",1,1,2,2)+-0.000077*ewgc("Cqu1R",0,0,2,2)+-0.000116*ewgc("Cqu1R",1,1,2,2)
                    +0.000149*ewgc("Cqu1R",2,2,0,0)+-0.000054*ewgc("Cqu1R",2,2,1,1)+-0.000248*ewgc("Cud1R",2,2,0,0)
                    +-0.000072*ewgc("Cud1R",2,2,1,1)+0.000417*ewgc("CuuR",0,0,2,2)+-0.000049*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_ttz_bin_0_75/sigma_ttz_bin_0_75_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 75 && b_max == 150){
        
        double SM_sigma_ttz_bin_75_150 = SM.getOptionalParameter("SM_sigma_ttz_bin_75_150_CMS_190711270");
        double sigma_ttz_bin_75_150_madgraph = 0.209588; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
                    
            double total = SM_sigma_ttz_bin_75_150 + (-0.051622*ewgc("CG")+-0.021757*ewgc("CHq1R",2,2)+0.021813*ewgc("CHq3R",2,2)
                    +0.013632*ewgc("CHuR",2,2)+0.004262*ewgc("Cqd8R",2,2,0,0)+0.000521*ewgc("Cqd8R",2,2,1,1)
                    +0.030335*ewgc("Cqq1R",0,2,2,0)+0.000796*ewgc("Cqq1R",1,2,2,1)+0.129627*ewgc("Cqq3R",0,2,2,0)
                    +0.012107*ewgc("Cqq3R",1,2,2,1)+0.016607*ewgc("Cqu8R",0,0,2,2)+0.001371*ewgc("Cqu8R",1,1,2,2)
                    +0.004665*ewgc("Cqu8R",2,2,0,0)+0.000163*ewgc("Cqu8R",2,2,1,1)+-0.000026*ewgc("CuBR",2,2)
                    +0.002538*ewgc("Cud8R",2,2,0,0)+0.000335*ewgc("Cud8R",2,2,1,1)-0.067116*ewgc("CuGR",2,2)
                    +0.011624*ewgc("CuuR",0,2,2,0)+0.000337*ewgc("CuuR",1,2,2,1)+-0.000735*ewgc("CuWR",2,2)
                    +-0.000153*ewgc("CHD")+-0.012754*ewgc("CHl3R",0,0)+-0.012754*ewgc("CHl3R",1,1)+0.000067*ewgc("CHW")
                    +0.006566*ewgc("CHWB")+0.012822*ewgc("CllR",0,1,1,0)+-0.000189*ewgc("Cqd1R",2,2,0,0)+0.000022*ewgc("Cqd1R",2,2,1,1)
                    +-0.002576*ewgc("Cqq1R",0,0,2,2)+-0.000592*ewgc("Cqq1R",1,1,2,2)+0.012033*ewgc("Cqq3R",0,0,2,2)
                    +0.000970*ewgc("Cqq3R",1,1,2,2)+0.000170*ewgc("Cqu1R",0,0,2,2)+-0.000031*ewgc("Cqu1R",1,1,2,2)
                    +0.000301*ewgc("Cqu1R",2,2,0,0)+0.000043*ewgc("Cqu1R",2,2,1,1)+-0.000125*ewgc("Cud1R",2,2,0,0)
                    +0.000002*ewgc("Cud1R",2,2,1,1)+0.000863*ewgc("CuuR",0,0,2,2)+0.000061*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_ttz_bin_75_150/sigma_ttz_bin_75_150_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 150 && b_max == 250){
        
        double SM_sigma_ttz_bin_150_250 = SM.getOptionalParameter("SM_sigma_ttz_bin_150_250_CMS_190711270");
        double sigma_ttz_bin_150_250_madgraph = 0.132915; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
                    
            double total = SM_sigma_ttz_bin_150_250 + (-0.044763*ewgc("CG")-0.014558*ewgc("CHq1R",2,2)+0.014567*ewgc("CHq3R",2,2)
                    +0.010078*ewgc("CHuR",2,2)+0.003298*ewgc("Cqd8R",2,2,0,0)+0.000314*ewgc("Cqd8R",2,2,1,1)
                    +0.022847*ewgc("Cqq1R",0,2,2,0)+0.000511*ewgc("Cqq1R",1,2,2,1)+0.084393*ewgc("Cqq3R",0,2,2,0)
                    +0.006512*ewgc("Cqq3R",1,2,2,1)+0.010724*ewgc("Cqu8R",0,0,2,2)+0.000696*ewgc("Cqu8R",1,1,2,2)
                    +0.004171*ewgc("Cqu8R",2,2,0,0)+0.000098*ewgc("Cqu8R",2,2,1,1)+0.000118*ewgc("CuBR",2,2)
                    +0.002011*ewgc("Cud8R",2,2,0,0)+0.000196*ewgc("Cud8R",2,2,1,1)-0.045686*ewgc("CuGR",2,2)
                    +0.011095*ewgc("CuuR",0,2,2,0)+0.000255*ewgc("CuuR",1,2,2,1)+-0.000811*ewgc("CuWR",2,2)
                    +-0.000798*ewgc("CHD")+-0.008089*ewgc("CHl3R",0,0)+-0.008089*ewgc("CHl3R",1,1)+0.000024*ewgc("CHW")
                    +0.003455*ewgc("CHWB")+0.008087*ewgc("CllR",0,1,1,0)+-0.000114*ewgc("Cqd1R",2,2,0,0)+-0.000029*ewgc("Cqd1R",2,2,1,1)
                    +-0.001109*ewgc("Cqq1R",0,0,2,2)+-0.000367*ewgc("Cqq1R",1,1,2,2)+0.008379*ewgc("Cqq3R",0,0,2,2)
                    +0.000510*ewgc("Cqq3R",1,1,2,2)+0.000203*ewgc("Cqu1R",0,0,2,2)+-0.000032*ewgc("Cqu1R",1,1,2,2)
                    +0.000246*ewgc("Cqu1R",2,2,0,0)+-0.000005*ewgc("Cqu1R",2,2,1,1)+-0.000148*ewgc("Cud1R",2,2,0,0)
                    +-0.000031*ewgc("Cud1R",2,2,1,1)+0.000849*ewgc("CuuR",0,0,2,2)+0.000014*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_ttz_bin_150_250/sigma_ttz_bin_150_250_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 250 && b_max == 500){
        
        double SM_sigma_ttz_bin_250_500 = SM.getOptionalParameter("SM_sigma_ttz_bin_250_500_CMS_190711270");
        double sigma_ttz_bin_250_500_madgraph = 0.067464; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
                    
            double total = SM_sigma_ttz_bin_250_500 + (-0.031862*ewgc("CG")+-0.007411*ewgc("CHq1R",2,2)+0.007439*ewgc("CHq3R",2,2)
                    +0.005451*ewgc("CHuR",2,2)+0.002779*ewgc("Cqd8R",2,2,0,0)+0.000235*ewgc("Cqd8R",2,2,1,1)
                    +0.019482*ewgc("Cqq1R",0,2,2,0)+0.000381*ewgc("Cqq1R",1,2,2,1)+0.060750*ewgc("Cqq3R",0,2,2,0)
                    +0.003621*ewgc("Cqq3R",1,2,2,1)+0.007854*ewgc("Cqu8R",0,0,2,2)+0.000404*ewgc("Cqu8R",1,1,2,2)
                    +0.003938*ewgc("Cqu8R",2,2,0,0)+0.000089*ewgc("Cqu8R",2,2,1,1)+0.000042*ewgc("CuBR",2,2)
                    +0.001798*ewgc("Cud8R",2,2,0,0)+0.000162*ewgc("Cud8R",2,2,1,1)-0.026610*ewgc("CuGR",2,2)
                    +0.011353*ewgc("CuuR",0,2,2,0)+0.000235*ewgc("CuuR",1,2,2,1)+-0.000517*ewgc("CuWR",2,2)
                    +-0.000594*ewgc("CHD")+-0.004099*ewgc("CHl3R",0,0)+-0.004099*ewgc("CHl3R",1,1)+0.000046*ewgc("CHW")
                    +0.001552*ewgc("CHWB")+0.004130*ewgc("CllR",0,1,1,0)+-0.000078*ewgc("Cqd1R",2,2,0,0)+0.000009*ewgc("Cqd1R",2,2,1,1)
                    +-0.000036*ewgc("Cqq1R",0,0,2,2)+-0.000177*ewgc("Cqq1R",1,1,2,2)+0.006452*ewgc("Cqq3R",0,0,2,2)
                    +0.000324*ewgc("Cqq3R",1,1,2,2)+0.000333*ewgc("Cqu1R",0,0,2,2)+0.000016*ewgc("Cqu1R",1,1,2,2)
                    +0.000275*ewgc("Cqu1R",2,2,0,0)+0.000020*ewgc("Cqu1R",2,2,1,1)+-0.000138*ewgc("Cud1R",2,2,0,0)
                    +0.000004*ewgc("Cud1R",2,2,1,1)+0.000886*ewgc("CuuR",0,0,2,2)+0.000029*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_ttz_bin_250_500/sigma_ttz_bin_250_500_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    }  else throw std::runtime_error("wrong bin choice in sigma_ttz_diff_LO_CMS_190711270");
    
}





sigma_ttz_diff_LO_ATLAS_210312603::sigma_ttz_diff_LO_ATLAS_210312603(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{

    setParametersForObservable(make_vector<std::string>() << "SM_sigma_ttz_bin_0_40_ATLAS_210312603" << "SM_sigma_ttz_bin_40_70_ATLAS_210312603"
            << "SM_sigma_ttz_bin_70_110_ATLAS_210312603" << "SM_sigma_ttz_bin_110_160_ATLAS_210312603" << "SM_sigma_ttz_bin_160_220_ATLAS_210312603" 
            << "SM_sigma_ttz_bin_220_290_ATLAS_210312603" << "SM_sigma_ttz_bin_290_400_ATLAS_210312603");

    
}


double sigma_ttz_diff_LO_ATLAS_210312603::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    bool   flag_Quadratic= false; //Needs to be properly defined
   
    if(b_min == 0 && b_max == 40){
        
        double SM_sigma_ttz_bin_0_40 = SM.getOptionalParameter("SM_sigma_ttz_bin_0_40_ATLAS_210312603");
        double sigma_ttz_bin_0_40_madgraph = 0.064177; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = (sigma_ttz_bin_0_40_madgraph -0.010478*ewgc("CG")-0.005804*ewgc("CHq1R",2,2)
                    +0.005759*ewgc("CHq3R",2,2)+0.002967*ewgc("CHuR",2,2)+0.001126*ewgc("Cqd8R",2,2,0,0)
                    +0.000122*ewgc("Cqd8R",2,2,1,1)
                    +0.008742*ewgc("Cqq1R",0,2,2,0)+0.000211*ewgc("Cqq1R",1,2,2,1)
                    +0.045418*ewgc("Cqq3R",0,2,2,0)+0.004850*ewgc("Cqq3R",1,2,2,1)
                    +0.005860*ewgc("Cqu8R",0,0,2,2)+0.000532*ewgc("Cqu8R",1,1,2,2)
                    +0.000976*ewgc("Cqu8R",2,2,0,0)+0.000001*ewgc("Cqu8R",2,2,1,1)
                    -0.000170*ewgc("CuBR",2,2)+0.000666*ewgc("Cud8R",2,2,0,0)+0.000058*ewgc("Cud8R",2,2,1,1)
                    -0.018690*ewgc("CuGR",2,2)+0.002251*ewgc("CuuR",0,2,2,0)
                    +0.000036*ewgc("CuuR",1,2,2,1)-0.000098*ewgc("CuWR",2,2)
                    +0.000365*ewgc("CHD")+-0.003957*ewgc("CHl3R",0,0)+-0.003957*ewgc("CHl3R",1,1)+-0.000032*ewgc("CHW")
                    +0.002458*ewgc("CHWB")+0.003876*ewgc("CllR",0,1,1,0)+-0.000058*ewgc("Cqd1R",2,2,0,0)+-0.000030*ewgc("Cqd1R",2,2,1,1)
                    +-0.001074*ewgc("Cqq1R",0,0,2,2)+-0.000307*ewgc("Cqq1R",1,1,2,2)+0.004072*ewgc("Cqq3R",0,0,2,2)
                    +0.000317*ewgc("Cqq3R",1,1,2,2)+-0.000065*ewgc("Cqu1R",0,0,2,2)+-0.000069*ewgc("Cqu1R",1,1,2,2)
                    +0.000058*ewgc("Cqu1R",2,2,0,0)+-0.000030*ewgc("Cqu1R",2,2,1,1)+-0.000072*ewgc("Cud1R",2,2,0,0)
                    +-0.000031*ewgc("Cud1R",2,2,1,1)+0.000196*ewgc("CuuR",0,0,2,2)+-0.000029*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_ttz_bin_0_40/sigma_ttz_bin_0_40_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
        }
        
    } else if(b_min == 40 && b_max == 70){
        
        double SM_sigma_ttz_bin_40_70 = SM.getOptionalParameter("SM_sigma_ttz_bin_40_70_ATLAS_210312603");
        double sigma_ttz_bin_40_70_madgraph = 0.096642; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
                    
            double total = (sigma_ttz_bin_40_70_madgraph -0.018049*ewgc("CG")-0.009212*ewgc("CHq1R",2,2)
                    +0.009177*ewgc("CHq3R",2,2)+0.005078*ewgc("CHuR",2,2)+0.001840*ewgc("Cqd8R",2,2,0,0)
                    +0.000234*ewgc("Cqd8R",2,2,1,1)
                    +0.013449*ewgc("Cqq1R",0,2,2,0)+0.000372*ewgc("Cqq1R",1,2,2,1)
                    +0.065355*ewgc("Cqq3R",0,2,2,0)+0.006633*ewgc("Cqq3R",1,2,2,1)
                    +0.008494*ewgc("Cqu8R",0,0,2,2)+0.000775*ewgc("Cqu8R",1,1,2,2)
                    +0.001722*ewgc("Cqu8R",2,2,0,0)+0.000055*ewgc("Cqu8R",2,2,1,1)
                    -0.000190*ewgc("CuBR",2,2)+0.001097*ewgc("Cud8R",2,2,0,0)
                    +0.000148*ewgc("Cud8R",2,2,1,1)+-0.029042*ewgc("CuGR",2,2)
                    +0.003988*ewgc("CuuR",0,2,2,0)+0.000126*ewgc("CuuR",1,2,2,1)
                    -0.000106*ewgc("CuWR",2,2)
                    +0.000406*ewgc("CHD")+-0.005906*ewgc("CHl3R",0,0)+-0.005906*ewgc("CHl3R",1,1)+0.000019*ewgc("CHW")
                    +0.003564*ewgc("CHWB")+0.005878*ewgc("CllR",0,1,1,0)+-0.000048*ewgc("Cqd1R",2,2,0,0)+-0.000002*ewgc("Cqd1R",2,2,1,1)
                    +-0.001543*ewgc("Cqq1R",0,0,2,2)+-0.000396*ewgc("Cqq1R",1,1,2,2)+0.005957*ewgc("Cqq3R",0,0,2,2)
                    +0.000477*ewgc("Cqq3R",1,1,2,2)+0.000037*ewgc("Cqu1R",0,0,2,2)+-0.000046*ewgc("Cqu1R",1,1,2,2)
                    +0.000093*ewgc("Cqu1R",2,2,0,0)+0.000007*ewgc("Cqu1R",2,2,1,1)+-0.000040*ewgc("Cud1R",2,2,0,0)
                    +0.000012*ewgc("Cud1R",2,2,1,1)+0.000312*ewgc("CuuR",0,0,2,2)+0.000027*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_ttz_bin_40_70/sigma_ttz_bin_40_70_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 70 && b_max == 110){
        
        double SM_sigma_ttz_bin_70_110 = SM.getOptionalParameter("SM_sigma_ttz_bin_70_110_ATLAS_210312603");
        double sigma_ttz_bin_70_110_madgraph = 0.127419; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
                    
            double total = (sigma_ttz_bin_70_110_madgraph-0.028383*ewgc("CG")+-0.012875*ewgc("CHq1R",2,2)+0.012884*ewgc("CHq3R",2,2)
                    +0.007778*ewgc("CHuR",2,2)+0.002468*ewgc("Cqd8R",2,2,0,0)+0.000298*ewgc("Cqd8R",2,2,1,1)
                    +0.017859*ewgc("Cqq1R",0,2,2,0)+0.000464*ewgc("Cqq1R",1,2,2,1)
                    +0.080694*ewgc("Cqq3R",0,2,2,0)+0.007774*ewgc("Cqq3R",1,2,2,1)+0.010371*ewgc("Cqu8R",0,0,2,2)
                    +0.000893*ewgc("Cqu8R",1,1,2,2)+0.002610*ewgc("Cqu8R",2,2,0,0)+0.000063*ewgc("Cqu8R",2,2,1,1)
                    -0.000135*ewgc("CuBR",2,2)+0.001481*ewgc("Cud8R",2,2,0,0)+0.000183*ewgc("Cud8R",2,2,1,1)
                    -0.039573*ewgc("CuGR",2,2)+0.006317*ewgc("CuuR",0,2,2,0)+0.000171*ewgc("CuuR",1,2,2,1)
                    -0.000281*ewgc("CuWR",2,2)
                    +0.000076*ewgc("CHD")+-0.007774*ewgc("CHl3R",0,0)+-0.007774*ewgc("CHl3R",1,1)+0.000013*ewgc("CHW")
                    +0.004172*ewgc("CHWB")+0.007759*ewgc("CllR",0,1,1,0)+-0.000123*ewgc("Cqd1R",2,2,0,0)+-0.000006*ewgc("Cqd1R",2,2,1,1)
                    +-0.001720*ewgc("Cqq1R",0,0,2,2)+-0.000435*ewgc("Cqq1R",1,1,2,2)+0.007293*ewgc("Cqq3R",0,0,2,2)
                    +0.000589*ewgc("Cqq3R",1,1,2,2)+-0.000047*ewgc("Cqu1R",0,0,2,2)+-0.000039*ewgc("Cqu1R",1,1,2,2)
                    +0.000144*ewgc("Cqu1R",2,2,0,0)+-0.000001*ewgc("Cqu1R",2,2,1,1)+-0.000117*ewgc("Cud1R",2,2,0,0)
                    +-0.000012*ewgc("Cud1R",2,2,1,1)+0.000456*ewgc("CuuR",0,0,2,2)+0.000004*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_ttz_bin_70_110/sigma_ttz_bin_70_110_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 110 && b_max == 160){
        
        double SM_sigma_ttz_bin_110_160 = SM.getOptionalParameter("SM_sigma_ttz_bin_110_160_ATLAS_210312603");
        double sigma_ttz_bin_110_160_madgraph = 0.118917; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
                    
            double total = (sigma_ttz_bin_110_160_madgraph -0.032594*ewgc("CG")-0.012662*ewgc("CHq1R",2,2)
                    +0.012643*ewgc("CHq3R",2,2)+0.008236*ewgc("CHuR",2,2)+0.002540*ewgc("Cqd8R",2,2,0,0)
                    +0.000253*ewgc("Cqd8R",2,2,1,1)+0.017566*ewgc("Cqq1R",0,2,2,0)+0.000417*ewgc("Cqq1R",1,2,2,1)
                    +0.071863*ewgc("Cqq3R",0,2,2,0)+0.006279*ewgc("Cqq3R",1,2,2,1)+0.009224*ewgc("Cqu8R",0,0,2,2)
                    +0.000682*ewgc("Cqu8R",1,1,2,2)+0.002931*ewgc("Cqu8R",2,2,0,0)+0.000049*ewgc("Cqu8R",2,2,1,1)
                    +-0.000017*ewgc("CuBR",2,2)+0.001474*ewgc("Cud8R",2,2,0,0)+0.000150*ewgc("Cud8R",2,2,1,1)
                    -0.038866*ewgc("CuGR",2,2)+0.007314*ewgc("CuuR",0,2,2,0)
                    +0.000171*ewgc("CuuR",1,2,2,1)+-0.000599*ewgc("CuWR",2,2)
                    +-0.000370*ewgc("CHD")+-0.007261*ewgc("CHl3R",0,0)+-0.007261*ewgc("CHl3R",1,1)+-0.000003*ewgc("CHW")
                    +0.003479*ewgc("CHWB")+0.007219*ewgc("CllR",0,1,1,0)+-0.000093*ewgc("Cqd1R",2,2,0,0)+-0.000028*ewgc("Cqd1R",2,2,1,1)
                    +-0.001164*ewgc("Cqq1R",0,0,2,2)+-0.000361*ewgc("Cqq1R",1,1,2,2)+0.006786*ewgc("Cqq3R",0,0,2,2)
                    +0.000480*ewgc("Cqq3R",1,1,2,2)+0.000040*ewgc("Cqu1R",0,0,2,2)+-0.000053*ewgc("Cqu1R",1,1,2,2)
                    +0.000170*ewgc("Cqu1R",2,2,0,0)+-0.000016*ewgc("Cqu1R",2,2,1,1)+-0.000090*ewgc("Cud1R",2,2,0,0)
                    +-0.000038*ewgc("Cud1R",2,2,1,1)+0.000486*ewgc("CuuR",0,0,2,2)+-0.000013*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_ttz_bin_110_160/sigma_ttz_bin_110_160_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 160 && b_max == 220){
        
        double SM_sigma_ttz_bin_160_220 = SM.getOptionalParameter("SM_sigma_ttz_bin_160_220_ATLAS_210312603");
        double sigma_ttz_bin_160_220_madgraph = 0.086169; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
                    
            double total = (sigma_ttz_bin_160_220_madgraph -0.028727*ewgc("CG")+-0.009451*ewgc("CHq1R",2,2)
                    +0.009451*ewgc("CHq3R",2,2)+0.006505*ewgc("CHuR",2,2)+0.002109*ewgc("Cqd8R",2,2,0,0)
                    +0.000210*ewgc("Cqd8R",2,2,1,1)+0.014508*ewgc("Cqq1R",0,2,2,0)
                    +0.000323*ewgc("Cqq1R",1,2,2,1)+0.054025*ewgc("Cqq3R",0,2,2,0)+0.004204*ewgc("Cqq3R",1,2,2,1)
                    +0.006853*ewgc("Cqu8R",0,0,2,2)+0.000457*ewgc("Cqu8R",1,1,2,2)+0.002647*ewgc("Cqu8R",2,2,0,0)
                    +0.000061*ewgc("Cqu8R",2,2,1,1)+0.000070*ewgc("CuBR",2,2)+0.001274*ewgc("Cud8R",2,2,0,0)
                    +0.000130*ewgc("Cud8R",2,2,1,1)+-0.029532*ewgc("CuGR",2,2)
                    +0.006980*ewgc("CuuR",0,2,2,0)+0.000163*ewgc("CuuR",1,2,2,1)+-0.000537*ewgc("CuWR",2,2)
                    +-0.000506*ewgc("CHD")+-0.005252*ewgc("CHl3R",0,0)+-0.005252*ewgc("CHl3R",1,1)+0.000016*ewgc("CHW")
                    +0.002257*ewgc("CHWB")+0.005257*ewgc("CllR",0,1,1,0)+-0.000081*ewgc("Cqd1R",2,2,0,0)+-0.000017*ewgc("Cqd1R",2,2,1,1)
                    +-0.000799*ewgc("Cqq1R",0,0,2,2)+-0.000244*ewgc("Cqq1R",1,1,2,2)+0.005197*ewgc("Cqq3R",0,0,2,2)
                    +0.000326*ewgc("Cqq3R",1,1,2,2)+0.000178*ewgc("Cqu1R",0,0,2,2)+-0.000014*ewgc("Cqu1R",1,1,2,2)
                    +0.000180*ewgc("Cqu1R",2,2,0,0)+0.000002*ewgc("Cqu1R",2,2,1,1)+-0.000112*ewgc("Cud1R",2,2,0,0)
                    +-0.000009*ewgc("Cud1R",2,2,1,1)+0.000539*ewgc("CuuR",0,0,2,2)+0.000010*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_ttz_bin_160_220/sigma_ttz_bin_160_220_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 220 && b_max == 290){
        
        double SM_sigma_ttz_bin_220_290 = SM.getOptionalParameter("SM_sigma_ttz_bin_220_290_ATLAS_210312603");
        double sigma_ttz_bin_220_290_madgraph = 0.051619; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
                    
            double total = (sigma_ttz_bin_220_290_madgraph -0.020633*ewgc("CG")+-0.005708*ewgc("CHq1R",2,2)
                    +0.005730*ewgc("CHq3R",2,2)+0.004111*ewgc("CHuR",2,2)+0.001568*ewgc("Cqd8R",2,2,0,0)
                    +0.000156*ewgc("Cqd8R",2,2,1,1)+0.010800*ewgc("Cqq1R",0,2,2,0)
                    +0.000238*ewgc("Cqq1R",1,2,2,1)+0.036604*ewgc("Cqq3R",0,2,2,0)+0.002537*ewgc("Cqq3R",1,2,2,1)
                    +0.004688*ewgc("Cqu8R",0,0,2,2)+0.000287*ewgc("Cqu8R",1,1,2,2)+0.002098*ewgc("Cqu8R",2,2,0,0)
                    +0.000060*ewgc("Cqu8R",2,2,1,1)+0.000066*ewgc("CuBR",2,2)+0.000986*ewgc("Cud8R",2,2,0,0)
                    +0.000107*ewgc("Cud8R",2,2,1,1)-0.018706*ewgc("CuGR",2,2)
                    +0.005776*ewgc("CuuR",0,2,2,0)+0.000142*ewgc("CuuR",1,2,2,1)-0.000356*ewgc("CuWR",2,2)
                    +-0.000396*ewgc("CHD")+-0.003125*ewgc("CHl3R",0,0)+-0.003125*ewgc("CHl3R",1,1)+0.000030*ewgc("CHW")
                    +0.001241*ewgc("CHWB")+0.003151*ewgc("CllR",0,1,1,0)+-0.000024*ewgc("Cqd1R",2,2,0,0)+0.000006*ewgc("Cqd1R",2,2,1,1)
                    +-0.000280*ewgc("Cqq1R",0,0,2,2)+-0.000122*ewgc("Cqq1R",1,1,2,2)+0.003698*ewgc("Cqq3R",0,0,2,2)
                    +0.000221*ewgc("Cqq3R",1,1,2,2)+0.000151*ewgc("Cqu1R",0,0,2,2)+0.000009*ewgc("Cqu1R",1,1,2,2)
                    +0.000165*ewgc("Cqu1R",2,2,0,0)+0.000020*ewgc("Cqu1R",2,2,1,1)+-0.000062*ewgc("Cud1R",2,2,0,0)
                    +0.000012*ewgc("Cud1R",2,2,1,1)+0.000461*ewgc("CuuR",0,0,2,2)+0.000027*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_ttz_bin_220_290/sigma_ttz_bin_220_290_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 290 && b_max == 400){
        
        double SM_sigma_ttz_bin_290_400 = SM.getOptionalParameter("SM_sigma_ttz_bin_290_400_ATLAS_210312603");
        double sigma_ttz_bin_290_400_madgraph = 0.032507; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
                    
            double total = (sigma_ttz_bin_290_400_madgraph -0.015551*ewgc("CG")+-0.003576*ewgc("CHq1R",2,2)
                    +0.003588*ewgc("CHq3R",2,2)+0.002640*ewgc("CHuR",2,2)+0.001358*ewgc("Cqd8R",2,2,0,0)
                    +0.000110*ewgc("Cqd8R",2,2,1,1)+0.009422*ewgc("Cqq1R",0,2,2,0)
                    +0.000177*ewgc("Cqq1R",1,2,2,1)+0.029289*ewgc("Cqq3R",0,2,2,0)+0.001743*ewgc("Cqq3R",1,2,2,1)
                    +0.003812*ewgc("Cqu8R",0,0,2,2)+0.000192*ewgc("Cqu8R",1,1,2,2)+0.001922*ewgc("Cqu8R",2,2,0,0)
                    +0.000038*ewgc("Cqu8R",2,2,1,1)+0.000016*ewgc("CuBR",2,2)+0.000882*ewgc("Cud8R",2,2,0,0)
                    +0.000076*ewgc("Cud8R",2,2,1,1)-0.012764*ewgc("CuGR",2,2)
                    +0.005527*ewgc("CuuR",0,2,2,0)+0.000108*ewgc("CuuR",1,2,2,1)+-0.000264*ewgc("CuWR",2,2)
                    +-0.000295*ewgc("CHD")+-0.001978*ewgc("CHl3R",0,0)+-0.001978*ewgc("CHl3R",1,1)+0.000020*ewgc("CHW")
                    +0.000740*ewgc("CHWB")+0.001985*ewgc("CllR",0,1,1,0)+-0.000037*ewgc("Cqd1R",2,2,0,0)+-0.000001*ewgc("Cqd1R",2,2,1,1)
                    +-0.000103*ewgc("Cqq1R",0,0,2,2)+-0.000084*ewgc("Cqq1R",1,1,2,2)+0.003015*ewgc("Cqq3R",0,0,2,2)
                    +0.000159*ewgc("Cqq3R",1,1,2,2)+0.000151*ewgc("Cqu1R",0,0,2,2)+-0.000000*ewgc("Cqu1R",1,1,2,2)
                    +0.000122*ewgc("Cqu1R",2,2,0,0)+0.000006*ewgc("Cqu1R",2,2,1,1)+-0.000060*ewgc("Cud1R",2,2,0,0)
                    +-0.000001*ewgc("Cud1R",2,2,1,1)+0.000432*ewgc("CuuR",0,0,2,2)+0.000012*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_ttz_bin_290_400/sigma_ttz_bin_290_400_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    }
    else throw std::runtime_error("wrong bin choice in sigma_ttz_diff_LO_ATLAS_210312603");
   
}







sigma_ttz_diff_LO_ATLAS_231204450::sigma_ttz_diff_LO_ATLAS_231204450(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{

    setParametersForObservable(make_vector<std::string>() << "SM_sigma_ttz_bin_0_60_ATLAS_231204450" << "SM_sigma_ttz_bin_60_100_ATLAS_231204450"
            << "SM_sigma_ttz_bin_100_140_ATLAS_231204450" << "SM_sigma_ttz_bin_140_180_ATLAS_231204450" << "SM_sigma_ttz_bin_180_230_ATLAS_231204450"
            << "SM_sigma_ttz_bin_230_280_ATLAS_231204450" << "SM_sigma_ttz_bin_280_350_ATLAS_231204450" << "SM_sigma_ttz_bin_350_1000_ATLAS_231204450");

    
}


double sigma_ttz_diff_LO_ATLAS_231204450::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    bool   flag_Quadratic= false; //Needs to be properly defined
   
    if(b_min == 0 && b_max == 60){
        
        double SM_sigma_ttz_bin_0_60 = SM.getOptionalParameter("SM_sigma_ttz_bin_0_60_ATLAS_231204450");
        double sigma_ttz_bin_0_60_madgraph = 0.127174; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double total = SM_sigma_ttz_bin_0_60 + (+-0.021934*ewgc("CG")+-0.011795*ewgc("CHq1R",2,2)
            +0.011666*ewgc("CHq3R",2,2)+0.006189*ewgc("CHuR",2,2)+0.002304*ewgc("Cqd8R",2,2,0,0)
            +0.000241*ewgc("Cqd8R",2,2,1,1)+0.017168*ewgc("Cqq1R",0,2,2,0)+0.000433*ewgc("Cqq1R",1,2,2,1)
            +0.087805*ewgc("Cqq3R",0,2,2,0)+0.009151*ewgc("Cqq3R",1,2,2,1)+0.011486*ewgc("Cqu8R",0,0,2,2)
            +0.001023*ewgc("Cqu8R",1,1,2,2)+0.002106*ewgc("Cqu8R",2,2,0,0)+-0.000002*ewgc("Cqu8R",2,2,1,1)
            +-0.000289*ewgc("CuBR",2,2)+0.001302*ewgc("Cud8R",2,2,0,0)+0.000116*ewgc("Cud8R",2,2,1,1)
            +-0.037524*ewgc("CuGR",2,2)+0.004779*ewgc("CuuR",0,2,2,0)+0.000086*ewgc("CuuR",1,2,2,1)
            +-0.000247*ewgc("CuWR",2,2)
            +0.000577*ewgc("CHD")+-0.007807*ewgc("CHl3R",0,0)+-0.007807*ewgc("CHl3R",1,1)+0.004690*ewgc("CHWB")
            +0.007677*ewgc("CllR",0,1,1,0)+-0.000092*ewgc("Cqd1R",2,2,0,0)+-0.000094*ewgc("Cqd1R",2,2,1,1)
            +-0.001974*ewgc("Cqq1R",0,0,2,2)+-0.000549*ewgc("Cqq1R",1,1,2,2)+0.008083*ewgc("Cqq3R",0,0,2,2)
            +0.000640*ewgc("Cqq3R",1,1,2,2)+-0.000030*ewgc("Cqu1R",0,0,2,2)+-0.000129*ewgc("Cqu1R",1,1,2,2)
            +0.000094*ewgc("Cqu1R",2,2,0,0)+-0.000192*ewgc("Cud1R",2,2,0,0)+-0.000078*ewgc("Cud1R",2,2,1,1)
            +0.000339*ewgc("CuuR",0,0,2,2))
                    *(SM_sigma_ttz_bin_0_60/sigma_ttz_bin_0_60_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
        }
        
    } else if(b_min == 60 && b_max == 100){
        
        double SM_sigma_ttz_bin_60_100 = SM.getOptionalParameter("SM_sigma_ttz_bin_60_100_ATLAS_231204450");
        double sigma_ttz_bin_60_100_madgraph = 0.131287; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
                    
            double total = SM_sigma_ttz_bin_60_100 + (+-0.027680*ewgc("CG")+-0.013016*ewgc("CHq1R",2,2)
            +0.013157*ewgc("CHq3R",2,2)+0.007747*ewgc("CHuR",2,2)+0.002620*ewgc("Cqd8R",2,2,0,0)
            +0.000380*ewgc("Cqd8R",2,2,1,1)+0.018507*ewgc("Cqq1R",0,2,2,0)+0.000573*ewgc("Cqq1R",1,2,2,1)
            +0.084073*ewgc("Cqq3R",0,2,2,0)+0.008402*ewgc("Cqq3R",1,2,2,1)+0.010907*ewgc("Cqu8R",0,0,2,2)
            +0.001021*ewgc("Cqu8R",1,1,2,2)+0.002662*ewgc("Cqu8R",2,2,0,0)+0.000145*ewgc("Cqu8R",2,2,1,1)
            +-0.000059*ewgc("CuBR",2,2)+0.001529*ewgc("Cud8R",2,2,0,0)+0.000263*ewgc("Cud8R",2,2,1,1)
            +-0.040409*ewgc("CuGR",2,2)+0.006251*ewgc("CuuR",0,2,2,0)+0.000249*ewgc("CuuR",1,2,2,1)
            +-0.000289*ewgc("CuWR",2,2)
            +0.000262*ewgc("CHD")+-0.007922*ewgc("CHl3R",0,0)+-0.007922*ewgc("CHl3R",1,1)
            +0.000082*ewgc("CHW")+0.004519*ewgc("CHWB")+0.008054*ewgc("CllR",0,1,1,0)
            +-0.000051*ewgc("Cqd1R",2,2,0,0)+0.000055*ewgc("Cqd1R",2,2,1,1)+-0.001905*ewgc("Cqq1R",0,0,2,2)
            +-0.000372*ewgc("Cqq1R",1,1,2,2)+0.007480*ewgc("Cqq3R",0,0,2,2)+0.000691*ewgc("Cqq3R",1,1,2,2)
            +-0.000060*ewgc("Cqu1R",0,0,2,2)+0.000030*ewgc("Cqu1R",1,1,2,2)+0.000193*ewgc("Cqu1R",2,2,0,0)
            +0.000077*ewgc("Cqu1R",2,2,1,1)+-0.000077*ewgc("Cud1R",2,2,0,0)+0.000063*ewgc("Cud1R",2,2,1,1)
            +0.000561*ewgc("CuuR",0,0,2,2)+0.000089*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_ttz_bin_60_100/sigma_ttz_bin_60_100_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 100 && b_max == 140){
        
        double SM_sigma_ttz_bin_100_140 = SM.getOptionalParameter("SM_sigma_ttz_bin_100_140_ATLAS_231204450");
        double sigma_ttz_bin_100_140_madgraph = 0.106891; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
                    
            double total = SM_sigma_ttz_bin_100_140 + (+-0.027409*ewgc("CG")+-0.011223*ewgc("CHq1R",2,2)
            +0.011239*ewgc("CHq3R",2,2)+0.007180*ewgc("CHuR",2,2)+0.002215*ewgc("Cqd8R",2,2,0,0)
            +0.000267*ewgc("Cqd8R",2,2,1,1)+0.015523*ewgc("Cqq1R",0,2,2,0)+0.000417*ewgc("Cqq1R",1,2,2,1)
            +0.064904*ewgc("Cqq3R",0,2,2,0)+0.005949*ewgc("Cqq3R",1,2,2,1)+0.008353*ewgc("Cqu8R",0,0,2,2)
            +0.000678*ewgc("Cqu8R",1,1,2,2)+0.002497*ewgc("Cqu8R",2,2,0,0)+0.000075*ewgc("Cqu8R",2,2,1,1)
            +0.000008*ewgc("CuBR",2,2)+0.001307*ewgc("Cud8R",2,2,0,0)+0.000171*ewgc("Cud8R",2,2,1,1)
            -0.034466*ewgc("CuGR",2,2)+0.006148*ewgc("CuuR",0,2,2,0)+0.000180*ewgc("CuuR",1,2,2,1)
            +-0.000467*ewgc("CuWR",2,2)
            +-0.000192*ewgc("CHD")+-0.006497*ewgc("CHl3R",0,0)+-0.006497*ewgc("CHl3R",1,1)
            +0.000025*ewgc("CHW")+0.003234*ewgc("CHWB")+0.006517*ewgc("CllR",0,1,1,0)
            +-0.000091*ewgc("Cqd1R",2,2,0,0)+0.000008*ewgc("Cqd1R",2,2,1,1)+-0.001327*ewgc("Cqq1R",0,0,2,2)
            +-0.000321*ewgc("Cqq1R",1,1,2,2)+0.006050*ewgc("Cqq3R",0,0,2,2)+0.000473*ewgc("Cqq3R",1,1,2,2)
            +0.000067*ewgc("Cqu1R",0,0,2,2)+-0.000014*ewgc("Cqu1R",1,1,2,2)+0.000157*ewgc("Cqu1R",2,2,0,0)
            +0.000014*ewgc("Cqu1R",2,2,1,1)+-0.000087*ewgc("Cud1R",2,2,0,0)+0.000004*ewgc("Cud1R",2,2,1,1)
            +0.000476*ewgc("CuuR",0,0,2,2)+0.000028*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_ttz_bin_100_140/sigma_ttz_bin_100_140_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 140 && b_max == 180){
        
        double SM_sigma_ttz_bin_140_180 = SM.getOptionalParameter("SM_sigma_ttz_bin_140_180_ATLAS_231204450");
        double sigma_ttz_bin_140_180_madgraph = 0.076251; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
                    
            double total = SM_sigma_ttz_bin_140_180 + (+-0.023096*ewgc("CG")+-0.008249*ewgc("CHq1R",2,2)
            +0.008298*ewgc("CHq3R",2,2)+0.005563*ewgc("CHuR",2,2)+0.001755*ewgc("Cqd8R",2,2,0,0)
            +0.000205*ewgc("Cqd8R",2,2,1,1)+0.011891*ewgc("Cqq1R",0,2,2,0)+0.000314*ewgc("Cqq1R",1,2,2,1)
            +0.046442*ewgc("Cqq3R",0,2,2,0)+0.003888*ewgc("Cqq3R",1,2,2,1)+0.005906*ewgc("Cqu8R",0,0,2,2)
            +0.000446*ewgc("Cqu8R",1,1,2,2)+0.002075*ewgc("Cqu8R",2,2,0,0)+0.000079*ewgc("Cqu8R",2,2,1,1)
            +0.000070*ewgc("CuBR",2,2)+0.001065*ewgc("Cud8R",2,2,0,0)+0.000137*ewgc("Cud8R",2,2,1,1)
            +-0.025432*ewgc("CuGR",2,2)+0.005331*ewgc("CuuR",0,2,2,0)+0.000174*ewgc("CuuR",1,2,2,1)
            +-0.000417*ewgc("CuWR",2,2)
            +-0.000310*ewgc("CHD")+-0.004620*ewgc("CHl3R",0,0)+-0.004620*ewgc("CHl3R",1,1)+0.000023*ewgc("CHW")
            +0.002128*ewgc("CHWB")+0.004681*ewgc("CllR",0,1,1,0)+-0.000035*ewgc("Cqd1R",2,2,0,0)
            +0.000012*ewgc("Cqd1R",2,2,1,1)+-0.000730*ewgc("Cqq1R",0,0,2,2)+-0.000201*ewgc("Cqq1R",1,1,2,2)
            +0.004458*ewgc("Cqq3R",0,0,2,2)+0.000321*ewgc("Cqq3R",1,1,2,2)+0.000069*ewgc("Cqu1R",0,0,2,2)
            +0.000006*ewgc("Cqu1R",1,1,2,2)+0.000130*ewgc("Cqu1R",2,2,0,0)+0.000029*ewgc("Cqu1R",2,2,1,1)
            +-0.000038*ewgc("Cud1R",2,2,0,0)+0.000018*ewgc("Cud1R",2,2,1,1)+0.000376*ewgc("CuuR",0,0,2,2)
            +0.000038*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_ttz_bin_140_180/sigma_ttz_bin_140_180_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 180 && b_max == 230){
        
        double SM_sigma_ttz_bin_180_230 = SM.getOptionalParameter("SM_sigma_ttz_bin_180_230_ATLAS_231204450");
        double sigma_ttz_bin_180_230_madgraph = 0.061424; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
                    
            double total = SM_sigma_ttz_bin_180_230 + ( +-0.021559*ewgc("CG")+-0.006782*ewgc("CHq1R",2,2)
            +0.006767*ewgc("CHq3R",2,2)+0.004725*ewgc("CHuR",2,2)+0.001582*ewgc("Cqd8R",2,2,0,0)
            +0.000159*ewgc("Cqd8R",2,2,1,1)+0.010831*ewgc("Cqq1R",0,2,2,0)+0.000242*ewgc("Cqq1R",1,2,2,1)
            +0.039427*ewgc("Cqq3R",0,2,2,0)+0.002990*ewgc("Cqq3R",1,2,2,1)+0.005006*ewgc("Cqu8R",0,0,2,2)
            +0.000329*ewgc("Cqu8R",1,1,2,2)+0.002008*ewgc("Cqu8R",2,2,0,0)+0.000054*ewgc("Cqu8R",2,2,1,1)
            +0.000065*ewgc("CuBR",2,2)+0.000959*ewgc("Cud8R",2,2,0,0)+0.000107*ewgc("Cud8R",2,2,1,1)
            +-0.021340*ewgc("CuGR",2,2)+0.005369*ewgc("CuuR",0,2,2,0)+0.000130*ewgc("CuuR",1,2,2,1)
            +-0.000388*ewgc("CuWR",2,2)
            +-0.000403*ewgc("CHD")+-0.003740*ewgc("CHl3R",0,0)+-0.003740*ewgc("CHl3R",1,1)+0.000012*ewgc("CHW")
            +0.001550*ewgc("CHWB")+0.003740*ewgc("CllR",0,1,1,0)+-0.000027*ewgc("Cqd1R",2,2,0,0)
            +-0.000002*ewgc("Cqd1R",2,2,1,1)+-0.000533*ewgc("Cqq1R",0,0,2,2)+-0.000157*ewgc("Cqq1R",1,1,2,2)
            +0.003812*ewgc("Cqq3R",0,0,2,2)+0.000251*ewgc("Cqq3R",1,1,2,2)+0.000127*ewgc("Cqu1R",0,0,2,2)
            +-0.000004*ewgc("Cqu1R",1,1,2,2)+0.000135*ewgc("Cqu1R",2,2,0,0)+0.000011*ewgc("Cqu1R",2,2,1,1)
            +-0.000067*ewgc("Cud1R",2,2,0,0)+-0.000003*ewgc("Cud1R",2,2,1,1)+0.000420*ewgc("CuuR",0,0,2,2)
            +0.000013*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_ttz_bin_180_230/sigma_ttz_bin_180_230_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 230 && b_max == 280){
        
        double SM_sigma_ttz_bin_230_280 = SM.getOptionalParameter("SM_sigma_ttz_bin_230_280_ATLAS_231204450");
        double sigma_ttz_bin_230_280_madgraph = 0.036489; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
                    
            double total = SM_sigma_ttz_bin_230_280 + ( +-0.014603*ewgc("CG")+-0.004048*ewgc("CHq1R",2,2)
            +0.004054*ewgc("CHq3R",2,2)+0.002909*ewgc("CHuR",2,2)+0.001116*ewgc("Cqd8R",2,2,0,0)
            +0.000112*ewgc("Cqd8R",2,2,1,1)+0.007664*ewgc("Cqq1R",0,2,2,0)+0.000171*ewgc("Cqq1R",1,2,2,1)
            +0.026023*ewgc("Cqq3R",0,2,2,0)+0.001787*ewgc("Cqq3R",1,2,2,1)+0.003296*ewgc("Cqu8R",0,0,2,2)
            +0.000204*ewgc("Cqu8R",1,1,2,2)+0.001490*ewgc("Cqu8R",2,2,0,0)+0.000045*ewgc("Cqu8R",2,2,1,1)
            +0.000044*ewgc("CuBR",2,2)+0.000700*ewgc("Cud8R",2,2,0,0)+0.000078*ewgc("Cud8R",2,2,1,1)
            +-0.013239*ewgc("CuGR",2,2)+0.004109*ewgc("CuuR",0,2,2,0)+0.000100*ewgc("CuuR",1,2,2,1)
            +-0.000270*ewgc("CuWR",2,2)
            +-0.000293*ewgc("CHD")+-0.002210*ewgc("CHl3R",0,0)+-0.002210*ewgc("CHl3R",1,1)
            +0.000023*ewgc("CHW")+0.000870*ewgc("CHWB")+0.002232*ewgc("CllR",0,1,1,0)
            +-0.000029*ewgc("Cqd1R",2,2,0,0)+0.000008*ewgc("Cqd1R",2,2,1,1)+-0.000197*ewgc("Cqq1R",0,0,2,2)
            +-0.000086*ewgc("Cqq1R",1,1,2,2)+0.002618*ewgc("Cqq3R",0,0,2,2)+0.000159*ewgc("Cqq3R",1,1,2,2)
            +0.000100*ewgc("Cqu1R",0,0,2,2)+0.000005*ewgc("Cqu1R",1,1,2,2)+0.000102*ewgc("Cqu1R",2,2,0,0)
            +0.000014*ewgc("Cqu1R",2,2,1,1)+-0.000045*ewgc("Cud1R",2,2,0,0)+0.000005*ewgc("Cud1R",2,2,1,1)
            +0.000336*ewgc("CuuR",0,0,2,2)+0.000017*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_ttz_bin_230_280/sigma_ttz_bin_230_280_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 280 && b_max == 350){
        
        double SM_sigma_ttz_bin_280_350 = SM.getOptionalParameter("SM_sigma_ttz_bin_280_350_ATLAS_231204450");
        double sigma_ttz_bin_280_350_madgraph = 0.027390; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
                    
            double total = SM_sigma_ttz_bin_280_350 + (-0.012486*ewgc("CG")+-0.003025*ewgc("CHq1R",2,2)
            +0.003026*ewgc("CHq3R",2,2)+0.002208*ewgc("CHuR",2,2)+0.001030*ewgc("Cqd8R",2,2,0,0)
            +0.000085*ewgc("Cqd8R",2,2,1,1)+0.007220*ewgc("Cqq1R",0,2,2,0)+0.000133*ewgc("Cqq1R",1,2,2,1)
            +0.022915*ewgc("Cqq3R",0,2,2,0)+0.001400*ewgc("Cqq3R",1,2,2,1)+0.002949*ewgc("Cqu8R",0,0,2,2)
            +0.000152*ewgc("Cqu8R",1,1,2,2)+0.001441*ewgc("Cqu8R",2,2,0,0)+0.000028*ewgc("Cqu8R",2,2,1,1)
            +0.000023*ewgc("CuBR",2,2)+0.000662*ewgc("Cud8R",2,2,0,0)+0.000056*ewgc("Cud8R",2,2,1,1)
            +-0.010474*ewgc("CuGR",2,2)+0.004120*ewgc("CuuR",0,2,2,0)+0.000080*ewgc("CuuR",1,2,2,1)
            +-0.000219*ewgc("CuWR",2,2) 
            +-0.000243*ewgc("CHD")+-0.001670*ewgc("CHl3R",0,0)+-0.001670*ewgc("CHl3R",1,1)+0.000010*ewgc("CHW")
            +0.000626*ewgc("CHWB")+0.001667*ewgc("CllR",0,1,1,0)+-0.000035*ewgc("Cqd1R",2,2,0,0)
            +-0.000006*ewgc("Cqd1R",2,2,1,1)+-0.000058*ewgc("Cqq1R",0,0,2,2)+-0.000080*ewgc("Cqq1R",1,1,2,2)
            +0.002386*ewgc("Cqq3R",0,0,2,2)+0.000110*ewgc("Cqq3R",1,1,2,2)+0.000104*ewgc("Cqu1R",0,0,2,2)
            +-0.000005*ewgc("Cqu1R",1,1,2,2)+0.000086*ewgc("Cqu1R",2,2,0,0)+0.000001*ewgc("Cqu1R",2,2,1,1)
            +-0.000052*ewgc("Cud1R",2,2,0,0)+-0.000007*ewgc("Cud1R",2,2,1,1)+0.000324*ewgc("CuuR",0,0,2,2)
            +0.000004*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_ttz_bin_280_350/sigma_ttz_bin_280_350_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    } else if(b_min == 350 && b_max == 1000){
        
        double SM_sigma_ttz_bin_350_1000 = SM.getOptionalParameter("SM_sigma_ttz_bin_350_1000_ATLAS_231204450");
        double sigma_ttz_bin_350_1000_madgraph = 0.027540; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
                    
            double total = SM_sigma_ttz_bin_350_1000 + (-0.016044*ewgc("CG")+-0.002948*ewgc("CHq1R",2,2)
            +0.002943*ewgc("CHq3R",2,2)+0.002203*ewgc("CHuR",2,2)+0.001890*ewgc("Cqd8R",2,2,0,0)
            +0.000122*ewgc("Cqd8R",2,2,1,1)+0.013680*ewgc("Cqq1R",0,2,2,0)+0.000221*ewgc("Cqq1R",1,2,2,1)
            +0.038385*ewgc("Cqq3R",0,2,2,0)+0.001838*ewgc("Cqq3R",1,2,2,1)+0.005142*ewgc("Cqu8R",0,0,2,2)
            +0.000202*ewgc("Cqu8R",1,1,2,2)+0.002906*ewgc("Cqu8R",2,2,0,0)+0.000047*ewgc("Cqu8R",2,2,1,1)
            +-0.000016*ewgc("CuBR",2,2)+0.001290*ewgc("Cud8R",2,2,0,0)+0.000086*ewgc("Cud8R",2,2,1,1)
            +-0.012966*ewgc("CuGR",2,2)+0.008817*ewgc("CuuR",0,2,2,0)+0.000142*ewgc("CuuR",1,2,2,1)
            +-0.000262*ewgc("CuWR",2,2)
            +-0.000266*ewgc("CHD")+-0.001691*ewgc("CHl3R",0,0)+-0.001691*ewgc("CHl3R",1,1)+0.000018*ewgc("CHW")
            +0.000610*ewgc("CHWB")+0.001686*ewgc("CllR",0,1,1,0)+-0.000064*ewgc("Cqd1R",2,2,0,0)
            +-0.000007*ewgc("Cqd1R",2,2,1,1)+0.000124*ewgc("Cqq1R",0,0,2,2)+-0.000097*ewgc("Cqq1R",1,1,2,2)
            +0.004104*ewgc("Cqq3R",0,0,2,2)+0.000163*ewgc("Cqq3R",1,1,2,2)+0.000228*ewgc("Cqu1R",0,0,2,2)
            +-0.000002*ewgc("Cqu1R",1,1,2,2)+0.000200*ewgc("Cqu1R",2,2,0,0)+0.000000*ewgc("Cqu1R",2,2,1,1)
            +-0.000105*ewgc("Cud1R",2,2,0,0)+-0.000011*ewgc("Cud1R",2,2,1,1)+0.000717*ewgc("CuuR",0,0,2,2)
            +0.000011*ewgc("CuuR",1,1,2,2))
                    *(SM_sigma_ttz_bin_350_1000/sigma_ttz_bin_350_1000_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        
    }
    
    

}




//// s-channel 8 TeV ////

sigma_tb_8_LO::sigma_tb_8_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tb_8");
}

double sigma_tb_8_LO::computeThValue()
{
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    double sigma_tb_8_LO_madgraph = 4.142840; //pb?
    double sigma_tb_8_SM = SM.getOptionalParameter("SM_sigma_tb_8"); 

    
    
    if(flag_Quadratic){
        return  0.;
        }
    else{
         double total = sigma_tb_8_SM + (0.502537*ewgc("CHq3R",2,2)-1.915120*ewgc("CuWR",2,2)
                 +-0.504792*ewgc("CHl3R",0,0)+-0.504792*ewgc("CHl3R",1,1)+-0.014073*ewgc("CHudR",2,2)+0.505185*ewgc("CllR",0,1,1,0)
                 +2.585980*ewgc("Cqq1R",0,2,2,0)+0.179945*ewgc("Cqq1R",1,2,2,1)+15.532200*ewgc("Cqq3R",0,0,2,2)
                 +-2.597130*ewgc("Cqq3R",0,2,2,0)+1.080470*ewgc("Cqq3R",1,1,2,2)+-0.179831*ewgc("Cqq3R",1,2,2,1)
                 )*(sigma_tb_8_SM/sigma_tb_8_LO_madgraph);
         
         if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
         
    }

}


//// s-channel 13 TeV ////

sigma_tb_13_LO::sigma_tb_13_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tb_13");
}

double sigma_tb_13_LO::computeThValue()
{
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    double sigma_tb_13_LO_madgraph = 8.005990; //pb?
    double sigma_tb_13_SM = SM.getOptionalParameter("SM_sigma_tb_13"); 

    
    
    if(flag_Quadratic){
        return  0.;
        }
        else{
            double total = sigma_tb_13_SM + (0.970643*ewgc("CHq3R",2,2)-3.718540*ewgc("CuWR",2,2) 
                    +-0.975848*ewgc("CHl3R",0,0)+-0.975848*ewgc("CHl3R",1,1)+-0.026592*ewgc("CHudR",2,2)+0.975691*ewgc("CllR",0,1,1,0)
                    +5.362790*ewgc("Cqq1R",0,2,2,0)+0.520897*ewgc("Cqq1R",1,2,2,1)+32.193200*ewgc("Cqq3R",0,0,2,2)+-5.395460*ewgc("Cqq3R",0,2,2,0)
                    +3.133830*ewgc("Cqq3R",1,1,2,2)+-0.523430*ewgc("Cqq3R",1,2,2,1)
                    )*(sigma_tb_13_SM/sigma_tb_13_LO_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }

}





//// t-channel 7 TeV ////

sigma_tq_7_LO::sigma_tq_7_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tq_7");
}

double sigma_tq_7_LO::computeThValue()
{
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    double sigma_tq_7_LO_madgraph = 53.962300; //pb?
    double sigma_tq_7_SM = SM.getOptionalParameter("SM_sigma_tq_7"); 

    
    
    if(flag_Quadratic){
        return  0.;
        }
        else{
        
            double total = sigma_tq_7_SM +(6.512720*ewgc("CHq3R",2,2)-1.642790*ewgc("CuWR",2,2)
            -6.576310*ewgc("CHl3R",0,0)+-6.576310*ewgc("CHl3R",1,1)+0.008127*ewgc("CHudR",2,2)
            +6.580220*ewgc("CllR",0,1,1,0))*(sigma_tq_7_SM/sigma_tq_7_LO_madgraph);
            
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }

}






//// t-channel 8 TeV ////

sigma_tq_8_LO::sigma_tq_8_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tq_8");
}

double sigma_tq_8_LO::computeThValue()
{
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    double sigma_tq_8_LO_madgraph = 71.914000; //pb?
    double sigma_tq_8_SM = SM.getOptionalParameter("SM_sigma_tq_8"); 

    
    
    if(flag_Quadratic){
        return  0.;
        }
        else{
        
            double total = sigma_tq_8_SM +(8.681900*ewgc("CHq3R",2,2)-2.086100*ewgc("CuWR",2,2)
            -8.761960*ewgc("CHl3R",0,0)+-8.761960*ewgc("CHl3R",1,1)+8.771790*ewgc("CllR",0,1,1,0))*(sigma_tq_8_SM/sigma_tq_8_LO_madgraph);
            
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }

}









//// t-channel 13 TeV ////

sigma_tq_13_LO::sigma_tq_13_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tq_13");
}

double sigma_tq_13_LO::computeThValue()
{
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    double sigma_tq_13_LO_madgraph = 186.321000; //pb?
    double sigma_tq_13_SM = SM.getOptionalParameter("SM_sigma_tq_13"); 

    
    
    if(flag_Quadratic){
        return  0.;
        }
        else{
        
            double total = sigma_tq_13_SM +(22.491500*ewgc("CHq3R",2,2)-4.598280*ewgc("CuWR",2,2) 
            -22.710600*ewgc("CHl3R",0,0)+-22.710600*ewgc("CHl3R",1,1)+0.018216*ewgc("CHudR",2,2)
            +22.710300*ewgc("CllR",0,1,1,0))*(sigma_tq_13_SM/sigma_tq_13_LO_madgraph);
            
            
            
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }

}



//// taq ////

sigma_taq_LO_CMS::sigma_taq_LO_CMS(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_taq_CMS");
}

double sigma_taq_LO_CMS::computeThValue()
{
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    double sigma_taq_LO_madgraph = 2.176180; //pb?
    double sigma_taq_SM_CMS = SM.getOptionalParameter("SM_sigma_taq_CMS"); 

    
    
    if(flag_Quadratic){
        return  0.;
        }
        else{
           
            double total = sigma_taq_SM_CMS +(0.262660*ewgc("CHq3R",2,2)-0.003262*ewgc("CuBR",2,2)-0.030453*ewgc("CuWR",2,2)
                    +-0.401883*ewgc("CHl3R",0,0)+-0.401883*ewgc("CHl3R",1,1)+0.000348*ewgc("CHudR",2,2)
                    )*(sigma_taq_SM_CMS/sigma_taq_LO_madgraph);
            
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }

}


sigma_taq_LO_ATLAS::sigma_taq_LO_ATLAS(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_taq_ATLAS");
}

double sigma_taq_LO_ATLAS::computeThValue()
{
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    double sigma_taq_LO_madgraph = 2.176180; //pb?
    double sigma_taq_SM_ATLAS = SM.getOptionalParameter("SM_sigma_taq_ATLAS"); 

    
    
    if(flag_Quadratic){
        return  0.;
        }
        else{
        
            double total = sigma_taq_SM_ATLAS +(0.262660*ewgc("CHq3R",2,2)-0.003262*ewgc("CuBR",2,2)-0.030453*ewgc("CuWR",2,2)
                    +-0.401883*ewgc("CHl3R",0,0)+-0.401883*ewgc("CHl3R",1,1)+0.000348*ewgc("CHudR",2,2)
                    )*(sigma_taq_SM_ATLAS/sigma_taq_LO_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }

}



//// tzq ////

sigma_tzq_LO::sigma_tzq_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tzq");
}

double sigma_tzq_LO::computeThValue()
{
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    double sigma_tzq_LO_madgraph = 0.725915; //pb?
    double sigma_tzq_SM = SM.getOptionalParameter("SM_sigma_tzq"); 

    
    
    if(flag_Quadratic){
        return  0.;
        }
        else{
        
            double total = sigma_tzq_SM +(0.015938*ewgc("CHq1R",2,2)+0.141467*ewgc("CHq3R",2,2)
                    +0.004390*ewgc("CHuR",2,2)+-0.002539*ewgc("CuBR",2,2)-0.001949*ewgc("CuWR",2,2)
                    +0.016858*ewgc("CHD")+-0.000058*ewgc("CHdR",2,2)+-0.132067*ewgc("CHl3R",0,0)
                    +-0.132067*ewgc("CHl3R",1,1)+-0.000189*ewgc("CHudR",2,2)+0.060599*ewgc("CHWB")
                    +0.131970*ewgc("CllR",0,1,1,0))
                    *(sigma_tzq_SM/sigma_tzq_LO_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }

}






//// tw 7 TeV ////

sigma_tw_7_LO::sigma_tw_7_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tw_7");
}

double sigma_tw_7_LO::computeThValue()
{
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    double sigma_tw_7_LO_madgraph = 10.784300; //pb?
    double sigma_tw_7_SM = SM.getOptionalParameter("SM_sigma_tw_7"); 
      
    if(flag_Quadratic){
        return  0.;
        }
        else{
        
            double total = sigma_tw_7_SM + (1.299580*ewgc("CHq3R",2,2)+0.931471*ewgc("CuWR",2,2)
                    +-0.653786*ewgc("CHl3R",0,0)+-0.653786*ewgc("CHl3R",1,1)+0.010123*ewgc("CHudR",2,2)+0.653771*ewgc("CllR",0,1,1,0)
                    )*(sigma_tw_7_SM/sigma_tw_7_LO_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }

}







//// tw 8 TeV ////

sigma_tw_8_LO::sigma_tw_8_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tw_8");
}

double sigma_tw_8_LO::computeThValue()
{
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    double sigma_tw_8_LO_madgraph = 15.495600; //pb?
    double sigma_tw_8_SM = SM.getOptionalParameter("SM_sigma_tw_8"); 
      
    if(flag_Quadratic){
        return  0.;
        }
        else{
        
            double total = sigma_tw_8_SM + (1.867690*ewgc("CHq3R",2,2)+1.335140*ewgc("CuWR",2,2)
                    +-0.939459*ewgc("CHl3R",0,0)+-0.939459*ewgc("CHl3R",1,1)+0.014721*ewgc("CHudR",2,2)
                    +0.939357*ewgc("CllR",0,1,1,0)
                    )*(sigma_tw_8_SM/sigma_tw_8_LO_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }

}








//// tw 13 TeV ////

sigma_tw_13_LO::sigma_tw_13_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tw_13");
}

double sigma_tw_13_LO::computeThValue()
{
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    double sigma_tw_13_LO_madgraph = 51.021400; //pb?
    double sigma_tw_13_SM = SM.getOptionalParameter("SM_sigma_tw_13"); 
      
    if(flag_Quadratic){
        return  0.;
        }
        else{
        
            double total = sigma_tw_13_SM + (6.152940*ewgc("CHq3R",2,2)+4.361620*ewgc("CuWR",2,2)
                    +-3.094000*ewgc("CHl3R",0,0)+-3.094000*ewgc("CHl3R",1,1)+0.048982*ewgc("CHudR",2,2)
                    +3.091850*ewgc("CllR",0,1,1,0))
                    *(sigma_tw_13_SM/sigma_tw_13_LO_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }

}




//// ttw  ////

sigma_ttw_LO::sigma_ttw_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_ttw");
}

double sigma_ttw_LO::computeThValue()
{
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    double sigma_ttw_LO_madgraph = 0.516319; //pb?
    double sigma_ttw_SM = SM.getOptionalParameter("SM_sigma_ttw"); 

    
    
    if(flag_Quadratic){
        return  0.;
        }
        else{
        
            double total = sigma_ttw_SM +(-0.000325*ewgc("CHq1R",2,2)+0.000604*ewgc("CHq3R",2,2)-0.000146*ewgc("CHuR",2,2)
                    +0.184757*ewgc("Cqq1R",0,2,2,0)+0.010580*ewgc("Cqq1R",1,2,2,1)+0.872452*ewgc("Cqq3R",0,2,2,0)
                    +0.054539*ewgc("Cqq3R",1,2,2,1)+0.133233*ewgc("Cqu8R",0,0,2,2)+0.008296*ewgc("Cqu8R",1,1,2,2)
                    +0.000098*ewgc("CuBR",2,2)-0.125985*ewgc("CuGR",2,2)-0.005842*ewgc("CuWR",2,2)
                    +-0.031757*ewgc("CHl3R",0,0)+-0.031757*ewgc("CHl3R",1,1)+0.000243*ewgc("CHW")+-0.000192*ewgc("CHWB")
                    +0.031794*ewgc("CllR",0,1,1,0)+-0.012797*ewgc("Cqq1R",0,0,2,2)+-0.001191*ewgc("Cqq1R",1,1,2,2)
                    +0.098677*ewgc("Cqq3R",0,0,2,2)+0.004610*ewgc("Cqq3R",1,1,2,2)+0.003612*ewgc("Cqu1R",0,0,2,2)
                    +0.000106*ewgc("Cqu1R",1,1,2,2))
                    *(sigma_ttw_SM/sigma_ttw_LO_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }

}




//// R_ttw  ////

R_ttw_LO::R_ttw_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_R_ttw" << "SM_sigma_ttw");

}

double R_ttw_LO::computeThValue()
{
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    double sigma_ttwp_LO_madgraph = 0.349605; //pb?
    double sigma_ttwm_LO_madgraph = 0.166990; //pb?
    double sigma_ttw_SM = SM.getOptionalParameter("SM_sigma_ttw"); 
    double R_ttw_SM = SM.getOptionalParameter("SM_R_ttw"); 

    
    double sigma_ttwp_SM = sigma_ttw_SM/(1+1/R_ttw_SM);
    double sigma_ttwm_SM = sigma_ttw_SM/(1+R_ttw_SM);
    
   
    
    
    if(flag_Quadratic){
        return  0.;
    }
    else{
            
        double sigma_ttwp_NP = -0.000229*ewgc("CHq1R",2,2)+0.000409*ewgc("CHq3R",2,2)+-0.000105*ewgc("CHuR",2,2)
        +0.131627*ewgc("Cqq1R",0,2,2,0)+0.004173*ewgc("Cqq1R",1,2,2,1)+0.622300*ewgc("Cqq3R",0,2,2,0)
        +0.029216*ewgc("Cqq3R",1,2,2,1)+0.095061*ewgc("Cqu8R",0,0,2,2)+0.004239*ewgc("Cqu8R",1,1,2,2)
        +0.000066*ewgc("CuBR",2,2)-0.085559*ewgc("CuGR",2,2)-0.004044*ewgc("CuWR",2,2)
        +-0.000229*ewgc("CHq1R",2,2)+0.000409*ewgc("CHq3R",2,2)+-0.000105*ewgc("CHuR",2,2)+0.131627*ewgc("Cqq1R",0,2,2,0)
        +0.004173*ewgc("Cqq1R",1,2,2,1)+0.622300*ewgc("Cqq3R",0,2,2,0)+0.029216*ewgc("Cqq3R",1,2,2,1)
        +0.095061*ewgc("Cqu8R",0,0,2,2)+0.004239*ewgc("Cqu8R",1,1,2,2)+0.000066*ewgc("CuBR",2,2)
        +-0.085559*ewgc("CuGR",2,2)+-0.004044*ewgc("CuWR",2,2)
        +-0.021513*ewgc("CHl3R",0,0)+-0.021513*ewgc("CHl3R",1,1)+0.000171*ewgc("CHW")+-0.000141*ewgc("CHWB")
        +0.021524*ewgc("CllR",0,1,1,0)+-0.007924*ewgc("Cqq1R",0,0,2,2)+-0.000761*ewgc("Cqq1R",1,1,2,2)
        +0.068930*ewgc("Cqq3R",0,0,2,2)+0.002428*ewgc("Cqq3R",1,1,2,2)+0.003035*ewgc("Cqu1R",0,0,2,2)
        +-0.000034*ewgc("Cqu1R",1,1,2,2);
        double sigma_ttwm_NP = -0.000104*ewgc("CHq1R",2,2)+0.000188*ewgc("CHq3R",2,2)+-0.000050*ewgc("CHuR",2,2)
        +0.051853*ewgc("Cqq1R",0,2,2,0)+0.006487*ewgc("Cqq1R",1,2,2,1)+0.251097*ewgc("Cqq3R",0,2,2,0)
        +0.025657*ewgc("Cqq3R",1,2,2,1)+0.038346*ewgc("Cqu8R",0,0,2,2)+0.004035*ewgc("Cqu8R",1,1,2,2)
        +0.000030*ewgc("CuBR",2,2)-0.040504*ewgc("CuGR",2,2)-0.001812*ewgc("CuWR",2,2)
        +-0.000104*ewgc("CHq1R",2,2)+0.000188*ewgc("CHq3R",2,2)+-0.000050*ewgc("CHuR",2,2)+0.051853*ewgc("Cqq1R",0,2,2,0)
        +0.006487*ewgc("Cqq1R",1,2,2,1)+0.251097*ewgc("Cqq3R",0,2,2,0)+0.025657*ewgc("Cqq3R",1,2,2,1)
        +0.038346*ewgc("Cqu8R",0,0,2,2)+0.004035*ewgc("Cqu8R",1,1,2,2)+0.000030*ewgc("CuBR",2,2)
        +-0.040504*ewgc("CuGR",2,2)+-0.001812*ewgc("CuWR",2,2)
        +-0.010266*ewgc("CHl3R",0,0)+-0.010266*ewgc("CHl3R",1,1)+0.000075*ewgc("CHW")+-0.000066*ewgc("CHWB")
        +0.010267*ewgc("CllR",0,1,1,0)+-0.003590*ewgc("Cqq1R",0,0,2,2)+-0.000254*ewgc("Cqq1R",1,1,2,2)
        +0.025859*ewgc("Cqq3R",0,0,2,2)+0.002285*ewgc("Cqq3R",1,1,2,2)+0.000816*ewgc("Cqu1R",0,0,2,2)
        +0.000098*ewgc("Cqu1R",1,1,2,2);
        
        double sigma_ttwp_NP_Corrected = sigma_ttwp_NP*(sigma_ttwp_SM/sigma_ttwp_LO_madgraph);
        double sigma_ttwm_NP_Corrected = sigma_ttwm_NP*(sigma_ttwm_SM/sigma_ttwm_LO_madgraph);
            
        double total = R_ttw_SM + (sigma_ttwp_NP_Corrected - sigma_ttwm_NP_Corrected*sigma_ttwp_SM/sigma_ttwm_SM)/sigma_ttwm_SM ;
        
        if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
        return total;
        
    }

}



/////// Helicities  /////////////////////////////////////////////////////////////////////////////////////



F0_LO::F0_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_F0");
}

double F0_LO::computeThValue()
{
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    
    double width_total_madgraph = 1.4724;
    double width_long_madgraph = 1.0252;
    
    double F0_madgraph = width_long_madgraph/width_total_madgraph;
    double F0_SM = SM.getOptionalParameter("SM_F0"); 

    
    
    if(flag_Quadratic){
        return  0.;
        }
        else{
        
            //The dependence with CHq3R vanish when dividing by the total width so we have not included it in the final result
            //The dependences with CdWR and CHudR vanish in the limit mb->0 
            //double width_total_NP = 0.009300*ewgc("CdWR",2,2)+0.178518*ewgc("CHq3R",2,2)+0.002502*ewgc("CHudR",2,2)+-0.246222*ewgc("CuWR",2,2);
            //double width_long_NP = +0.003102*ewgc("CdWR",2,2)+0.124307*ewgc("CHq3R",2,2)+0.000840*ewgc("CHudR",2,2)+-0.082068*ewgc("CuWR",2,2);
            
            double width_total_NP = 0.009300*ewgc("CdWR",2,2)+0.002502*ewgc("CHudR",2,2)-0.246222*ewgc("CuWR",2,2)+-0.089460*ewgc("CHl3R",0,0)+-0.089460*ewgc("CHl3R",1,1)+0.089460*ewgc("CllR",0,1,1,0);
            double width_long_NP = +0.003102*ewgc("CdWR",2,2)+0.000840*ewgc("CHudR",2,2)-0.082068*ewgc("CuWR",2,2)+-0.062327*ewgc("CHl3R",0,0)+-0.062327*ewgc("CHl3R",1,1)+0.062327*ewgc("CllR",0,1,1,0);
            
            double F0_NP = (width_long_NP - width_total_NP*width_long_madgraph/width_total_madgraph)/width_total_madgraph;
            
            double total = F0_SM + (F0_NP)*(F0_SM/F0_madgraph);
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
}



FL_LO::FL_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_FL");
}

double FL_LO::computeThValue()
{
    
    
    
    bool   flag_Quadratic= false;//Needs to be properly defined
    
    double width_total_madgraph = 1.4724;
    double width_transleft_madgraph = 0.446800;
    double FL_madgraph = width_transleft_madgraph/width_total_madgraph;

    double FL_SM = SM.getOptionalParameter("SM_FL");
    
    

    if(flag_Quadratic){
        return  0.;
        }
    else{
        
        
        //The dependence with CHq3R vanish when dividing by the total width so we have not included it in the final result
        //The dependences with CdWR and CHudR vanish in the limit mb->0 
        //double width_total_NP = 0.009300*ewgc("CdWR",2,2)+0.178518*ewgc("CHq3R",2,2)+0.002502*ewgc("CHudR",2,2)+-0.246222*ewgc("CuWR",2,2);
        //double width_transleft_NP = +0.001108*ewgc("CdWR",2,2)+0.054172*ewgc("CHq3R",2,2)+0.000840*ewgc("CHudR",2,2)+-0.164114*ewgc("CuWR",2,2);
            
        double width_total_NP = 0.009300*ewgc("CdWR",2,2)+0.002502*ewgc("CHudR",2,2)-0.246222*ewgc("CuWR",2,2)
        +-0.089460*ewgc("CHl3R",0,0)+-0.089460*ewgc("CHl3R",1,1)+0.089460*ewgc("CllR",0,1,1,0);
        double width_transleft_NP = +0.001108*ewgc("CdWR",2,2)+0.000840*ewgc("CHudR",2,2)+-0.164114*ewgc("CuWR",2,2)
        -0.027124*ewgc("CHl3R",0,0)+-0.027124*ewgc("CHl3R",1,1)+0.027124*ewgc("CllR",0,1,1,0);
            
        double FL_NP = (width_transleft_NP - width_total_NP*width_transleft_madgraph/width_total_madgraph)/width_total_madgraph;
            
        double total = FL_SM + (FL_NP)*(FL_SM/FL_madgraph);
        
        if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
        
        
    }

    
}





FR_LO::FR_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_FR");
}

double FR_LO::computeThValue()
{
    
    
    
    bool   flag_Quadratic= false;//Needs to be properly defined
    
    double width_total_madgraph = 1.4724;
    double width_transright_madgraph = 0.000430;
    double FR_madgraph = width_transright_madgraph/width_total_madgraph;

    double FR_SM = SM.getOptionalParameter("SM_FR");
    
    

    if(flag_Quadratic){
        return  0.;
        }
    else{
        
        
        //The dependence with CHq3R vanish when dividing by the total width so we have not included it in the final result
        //The dependences with CdWR and CHudR vanish in the limit mb->0 
        //double width_total_NP = 0.009300*ewgc("CdWR",2,2)+0.178518*ewgc("CHq3R",2,2)+0.002502*ewgc("CHudR",2,2)+-0.246222*ewgc("CuWR",2,2);
        //double width_transright_NP = 0.005088*ewgc("CdWR",2,2)+0.000052*ewgc("CHq3R",2,2)+0.000840*ewgc("CHudR",2,2)+-0.000034*ewgc("CuWR",2,2);
            
        double width_total_NP = 0.009300*ewgc("CdWR",2,2)+0.002502*ewgc("CHudR",2,2)-0.246222*ewgc("CuWR",2,2)
        +-0.089460*ewgc("CHl3R",0,0)+-0.089460*ewgc("CHl3R",1,1)+0.089460*ewgc("CllR",0,1,1,0);
        double width_transright_NP = 0.005088*ewgc("CdWR",2,2)+0.000840*ewgc("CHudR",2,2)+-0.000034*ewgc("CuWR",2,2);
            
        double FR_NP = (width_transright_NP - width_total_NP*width_transright_madgraph/width_total_madgraph)/width_total_madgraph;
            
        double total = FR_SM + (FR_NP)*(FR_SM/FR_madgraph);
        
        if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
        
    }

    
}



Test_direct::Test_direct(const StandardModel& SM_i) 
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{ 

    setParametersForObservable(make_vector<std::string>() << "Test_direct_in");
}

double Test_direct::computeThValue()
{
    double Test_direct_in = SM.getOptionalParameter("Test_direct_in");
    
    return Test_direct_in;

    
}



