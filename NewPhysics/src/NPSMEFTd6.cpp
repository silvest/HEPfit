/*
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NPSMEFTd6.h"
#include <limits>
#include <gsl/gsl_sf.h>
#include <boost/bind.hpp>
#include "gslpp_function_adapter.h"

const std::string NPSMEFTd6::NPSMEFTd6Vars[NNPSMEFTd6Vars]
        = {"CG", "CW", "C2B", "C2W", "C2BS", "C2WS", "CHG", "CHW", "CHB", "CDHB", "CDHW", "CDB", "CDW", "CHWB", "CHD", "CT", "CHbox", "CH",
    "CHL1_11", "CHL1_12r", "CHL1_13r", "CHL1_22", "CHL1_23r", "CHL1_33",
    "CHL1_12i", "CHL1_13i", "CHL1_23i",
    "CHL3_11", "CHL3_12r", "CHL3_13r", "CHL3_22", "CHL3_23r", "CHL3_33",
    "CHL3_12i", "CHL3_13i", "CHL3_23i",
    "CHe_11", "CHe_12r", "CHe_13r", "CHe_22", "CHe_23r", "CHe_33",
    "CHe_12i", "CHe_13i", "CHe_23i",
    "CHQ1_11", "CHQ1_12r", "CHQ1_13r", "CHQ1_22", "CHQ1_23r", "CHQ1_33",
    "CHQ1_12i", "CHQ1_13i", "CHQ1_23i",
    "CHQ3_11", "CHQ3_12r", "CHQ3_13r", "CHQ3_22", "CHQ3_23r", "CHQ3_33",
    "CHQ3_12i", "CHQ3_13i", "CHQ3_23i",
    "CHu_11", "CHu_12r", "CHu_13r", "CHu_22", "CHu_23r", "CHu_33",
    "CHu_12i", "CHu_13i", "CHu_23i",
    "CHd_11", "CHd_12r", "CHd_13r", "CHd_22", "CHd_23r", "CHd_33",
    "CHd_12i", "CHd_13i", "CHd_23i",
    "CHud_11r", "CHud_12r", "CHud_13r", "CHud_22r", "CHud_23r", "CHud_33r",
    "CHud_11i", "CHud_12i", "CHud_13i", "CHud_22i", "CHud_23i", "CHud_33i",
    "CeH_11r", "CeH_12r", "CeH_13r", "CeH_22r", "CeH_23r", "CeH_33r",
    "CeH_11i", "CeH_12i", "CeH_13i", "CeH_22i", "CeH_23i", "CeH_33i",
    "CuH_11r", "CuH_12r", "CuH_13r", "CuH_22r", "CuH_23r", "CuH_33r",
    "CuH_11i", "CuH_12i", "CuH_13i", "CuH_22i", "CuH_23i", "CuH_33i",
    "CdH_11r", "CdH_12r", "CdH_13r", "CdH_22r", "CdH_23r", "CdH_33r",
    "CdH_11i", "CdH_12i", "CdH_13i", "CdH_22i", "CdH_23i", "CdH_33i",
    "CuG_11r", "CuG_12r", "CuG_13r", "CuG_22r", "CuG_23r", "CuG_33r",
    "CuG_11i", "CuG_12i", "CuG_13i", "CuG_22i", "CuG_23i", "CuG_33i",
    "CuW_11r", "CuW_12r", "CuW_13r", "CuW_22r", "CuW_23r", "CuW_33r",
    "CuW_11i", "CuW_12i", "CuW_13i", "CuW_22i", "CuW_23i", "CuW_33i",
    "CuB_11r", "CuB_12r", "CuB_13r", "CuB_22r", "CuB_23r", "CuB_33r",
    "CuB_11i", "CuB_12i", "CuB_13i", "CuB_22i", "CuB_23i", "CuB_33i",
    "CLL_1111","CLL_1221","CLL_1122",
    "CLL_1133","CLL_1331",
    "CLQ1_1111","CLQ1_1122","CLQ1_2211","CLQ1_1221","CLQ1_2112",
    "CLQ1_1133","CLQ1_3311","CLQ1_1331","CLQ1_3113",
    "CLQ1_1123","CLQ1_2223","CLQ1_3323",
    "CLQ1_1132","CLQ1_2232","CLQ1_3332",        
    "CLQ3_1111","CLQ3_1122","CLQ3_2211","CLQ3_1221","CLQ3_2112",
    "CLQ3_1133","CLQ3_3311","CLQ3_1331","CLQ3_3113",
    "CLQ3_1123","CLQ3_2223","CLQ3_3323",
    "CLQ3_1132","CLQ3_2232","CLQ3_3332",
    "Cee_1111","Cee_1122","Cee_1133",
    "Ceu_1111","Ceu_1122","Ceu_2211","Ceu_1133","Ceu_2233","Ceu_3311",
    "Ced_1111","Ced_1122","Ced_2211","Ced_1133","Ced_3311",
    "Ced_1123","Ced_2223","Ced_3323",
    "Ced_1132","Ced_2232","Ced_3332",
    "CLe_1111","CLe_1122","CLe_2211","CLe_1133","CLe_3311",
    "CLu_1111","CLu_1122","CLu_2211","CLu_1133","CLu_2233","CLu_3311",
    "CLd_1111","CLd_1122","CLd_2211","CLd_1133","CLd_3311",
    "CLd_1123","CLd_2223","CLd_3323",
    "CLd_1132","CLd_2232","CLd_3332",
    "CQe_1111","CQe_1122","CQe_2211","CQe_1133","CQe_3311", 
    "CQe_2311","CQe_2322","CQe_2333",
    "CQe_3211","CQe_3222","CQe_3233",
    "CLedQ_11","CLedQ_22","CpLedQ_11","CpLedQ_22",
    "Lambda_NP",
    "BrHinv","BrHexo",
    "dg1Z","dKappaga","lambZ",
    "eggFint","eggFpar","ettHint","ettHpar",
    "eVBFint","eVBFpar","eWHint","eWHpar","eZHint","eZHpar",
    "eeeWBFint","eeeWBFpar","eeeZHint","eeeZHpar","eeettHint","eeettHpar",
    "eepWBFint","eepWBFpar","eepZBFint","eepZBFpar",
    "eHggint","eHggpar","eHWWint","eHWWpar","eHZZint","eHZZpar","eHZgaint","eHZgapar",
    "eHgagaint","eHgagapar","eHmumuint","eHmumupar","eHtautauint","eHtautaupar",
    "eHccint","eHccpar","eHbbint","eHbbpar",
    "eggFHgaga","eggFHZga","eggFHZZ","eggFHWW","eggFHtautau","eggFHbb","eggFHmumu",   
    "eVBFHgaga","eVBFHZga","eVBFHZZ","eVBFHWW","eVBFHtautau","eVBFHbb","eVBFHmumu",   
    "eWHgaga","eWHZga","eWHZZ","eWHWW","eWHtautau","eWHbb","eWHmumu",    
    "eZHgaga","eZHZga","eZHZZ","eZHWW","eZHtautau","eZHbb","eZHmumu",
    "ettHgaga","ettHZga","ettHZZ","ettHWW","ettHtautau","ettHbb","ettHmumu",
    "eVBFHinv","eVHinv",
    "eVBF_2_Hbox", "eVBF_2_HQ1_11", "eVBF_2_Hu_11", "eVBF_2_Hd_11", "eVBF_2_HQ3_11",
    "eVBF_2_HD", "eVBF_2_HB", "eVBF_2_HW", "eVBF_2_HWB", "eVBF_2_HG", "eVBF_2_DHB",
    "eVBF_2_DHW", "eVBF_2_DeltaGF",
    "eVBF_78_Hbox", "eVBF_78_HQ1_11", "eVBF_78_Hu_11", "eVBF_78_Hd_11", "eVBF_78_HQ3_11",
    "eVBF_78_HD", "eVBF_78_HB", "eVBF_78_HW", "eVBF_78_HWB", "eVBF_78_HG", "eVBF_78_DHB",
    "eVBF_78_DHW", "eVBF_78_DeltaGF",
    "eVBF_1314_Hbox", "eVBF_1314_HQ1_11", "eVBF_1314_Hu_11", "eVBF_1314_Hd_11", "eVBF_1314_HQ3_11",
    "eVBF_1314_HD", "eVBF_1314_HB", "eVBF_1314_HW", "eVBF_1314_HWB", "eVBF_1314_HG", "eVBF_1314_DHB",
    "eVBF_1314_DHW", "eVBF_1314_DeltaGF",
    "eWH_2_Hbox", "eWH_2_HQ3_11", "eWH_2_HD", "eWH_2_HW", "eWH_2_HWB", "eWH_2_DHW", "eWH_2_DeltaGF",
    "eWH_78_Hbox", "eWH_78_HQ3_11", "eWH_78_HD", "eWH_78_HW", "eWH_78_HWB", "eWH_78_DHW", "eWH_78_DeltaGF",
    "eWH_1314_Hbox", "eWH_1314_HQ3_11", "eWH_1314_HD", "eWH_1314_HW", "eWH_1314_HWB", "eWH_1314_DHW", "eWH_1314_DeltaGF",
    "eZH_2_Hbox", "eZH_2_HQ1_11", "eZH_2_Hu_11", "eZH_2_Hd_11", "eZH_2_HQ3_11", "eZH_2_HD", "eZH_2_HB", "eZH_2_HW", "eZH_2_HWB", "eZH_2_DHB", "eZH_2_DHW", "eZH_2_DeltaGF",
    "eZH_78_Hbox", "eZH_78_HQ1_11", "eZH_78_Hu_11", "eZH_78_Hd_11", "eZH_78_HQ3_11", "eZH_78_HD", "eZH_78_HB", "eZH_78_HW", "eZH_78_HWB", "eZH_78_DHB", "eZH_78_DHW", "eZH_78_DeltaGF",
    "eZH_1314_Hbox", "eZH_1314_HQ1_11", "eZH_1314_Hu_11", "eZH_1314_Hd_11", "eZH_1314_HQ3_11", "eZH_1314_HD", "eZH_1314_HB", "eZH_1314_HW", "eZH_1314_HWB", "eZH_1314_DHB", "eZH_1314_DHW", "eZH_1314_DeltaGF",
    "ettH_2_HG", "ettH_2_G", "ettH_2_uG_33r", "ettH_2_DeltagHt",
    "ettH_78_HG", "ettH_78_G", "ettH_78_uG_33r", "ettH_78_DeltagHt",
    "ettH_1314_HG", "ettH_1314_G", "ettH_1314_uG_33r", "ettH_1314_DeltagHt"};

const std::string NPSMEFTd6::NPSMEFTd6VarsRot[NNPSMEFTd6Vars]
        = {"CG", "CW", "C2B", "C2W", "C2BS", "C2WS", "CHG", "CHWHB_gaga", "CHWHB_gagaorth", "CDHB", "CDHW", "CDB", "CDW", "CHWB", "CHD", "CT", "CHbox", "CH",
    "CHL1_11", "CHL1_12r", "CHL1_13r", "CHL1_22", "CHL1_23r", "CHL1_33",
    "CHL1_12i", "CHL1_13i", "CHL1_23i",
    "CHL3_11", "CHL3_12r", "CHL3_13r", "CHL3_22", "CHL3_23r", "CHL3_33",
    "CHL3_12i", "CHL3_13i", "CHL3_23i",
    "CHe_11", "CHe_12r", "CHe_13r", "CHe_22", "CHe_23r", "CHe_33",
    "CHe_12i", "CHe_13i", "CHe_23i",
    "CHQ1_11", "CHQ1_12r", "CHQ1_13r", "CHQ1_22", "CHQ1_23r", "CHQ1_33",
    "CHQ1_12i", "CHQ1_13i", "CHQ1_23i",
    "CHQ3_11", "CHQ3_12r", "CHQ3_13r", "CHQ3_22", "CHQ3_23r", "CHQ3_33",
    "CHQ3_12i", "CHQ3_13i", "CHQ3_23i",
    "CHu_11", "CHu_12r", "CHu_13r", "CHu_22", "CHu_23r", "CHu_33",
    "CHu_12i", "CHu_13i", "CHu_23i",
    "CHd_11", "CHd_12r", "CHd_13r", "CHd_22", "CHd_23r", "CHd_33",
    "CHd_12i", "CHd_13i", "CHd_23i",
    "CHud_11r", "CHud_12r", "CHud_13r", "CHud_22r", "CHud_23r", "CHud_33r",
    "CHud_11i", "CHud_12i", "CHud_13i", "CHud_22i", "CHud_23i", "CHud_33i",
    "CeH_11r", "CeH_12r", "CeH_13r", "CeH_22r", "CeH_23r", "CeH_33r",
    "CeH_11i", "CeH_12i", "CeH_13i", "CeH_22i", "CeH_23i", "CeH_33i",
    "CuH_11r", "CuH_12r", "CuH_13r", "CuH_22r", "CuH_23r", "CuH_33r",
    "CuH_11i", "CuH_12i", "CuH_13i", "CuH_22i", "CuH_23i", "CuH_33i",
    "CdH_11r", "CdH_12r", "CdH_13r", "CdH_22r", "CdH_23r", "CdH_33r",
    "CdH_11i", "CdH_12i", "CdH_13i", "CdH_22i", "CdH_23i", "CdH_33i",
    "CuG_11r", "CuG_12r", "CuG_13r", "CuG_22r", "CuG_23r", "CuG_33r",
    "CuG_11i", "CuG_12i", "CuG_13i", "CuG_22i", "CuG_23i", "CuG_33i",
    "CuW_11r", "CuW_12r", "CuW_13r", "CuW_22r", "CuW_23r", "CuW_33r",
    "CuW_11i", "CuW_12i", "CuW_13i", "CuW_22i", "CuW_23i", "CuW_33i",
    "CuB_11r", "CuB_12r", "CuB_13r", "CuB_22r", "CuB_23r", "CuB_33r",
    "CuB_11i", "CuB_12i", "CuB_13i", "CuB_22i", "CuB_23i", "CuB_33i",
    "CLL_1111","CLL_1221","CLL_1122",
    "CLL_1133","CLL_1331",
    "CLQ1_1111","CLQ1_1122","CLQ1_2211","CLQ1_1221","CLQ1_2112",
    "CLQ1_1133","CLQ1_3311","CLQ1_1331","CLQ1_3113",
    "CLQ1_1123","CLQ1_2223","CLQ1_3323",
    "CLQ1_1132","CLQ1_2232","CLQ1_3332",
    "CLQ3_1111","CLQ3_1122","CLQ3_2211","CLQ3_1221","CLQ3_2112",
    "CLQ3_1133","CLQ3_3311","CLQ3_1331","CLQ3_3113",
    "CLQ3_1123","CLQ3_2223","CLQ3_3323",
    "CLQ3_1132","CLQ3_2232","CLQ3_3332",
    "Cee_1111","Cee_1122","Cee_1133",
    "Ceu_1111","Ceu_1122","Ceu_2211","Ceu_1133","Ceu_2233","Ceu_3311",
    "Ced_1111","Ced_1122","Ced_2211","Ced_1133","Ced_3311",
    "Ced_1123","Ced_2223","Ced_3323",
    "Ced_1132","Ced_2232","Ced_3332",
    "CLe_1111","CLe_1122","CLe_2211","CLe_1133","CLe_3311",
    "CLu_1111","CLu_1122","CLu_2211","CLu_1133","CLu_2233","CLu_3311",
    "CLd_1111","CLd_1122","CLd_2211","CLd_1133","CLd_3311",
    "CLd_1123","CLd_2223","CLd_3323",
    "CLd_1132","CLd_2232","CLd_3332",
    "CQe_1111","CQe_1122","CQe_2211","CQe_1133","CQe_3311",
    "CQe_2311","CQe_2322","CQe_2333",
    "CQe_3211","CQe_3222","CQe_3233",
    "CLedQ_11","CLedQ_22","CpLedQ_11","CpLedQ_22",
    "Lambda_NP",
    "BrHinv","BrHexo",
    "dg1Z","dKappaga","lambZ",
    "eggFint","eggFpar","ettHint","ettHpar",
    "eVBFint","eVBFpar","eWHint","eWHpar","eZHint","eZHpar",
    "eeeWBFint","eeeWBFpar","eeeZHint","eeeZHpar","eeettHint","eeettHpar",
    "eepWBFint","eepWBFpar","eepZBFint","eepZBFpar",
    "eHggint","eHggpar","eHWWint","eHWWpar","eHZZint","eHZZpar","eHZgaint","eHZgapar",
    "eHgagaint","eHgagapar","eHmumuint","eHmumupar","eHtautauint","eHtautaupar",
    "eHccint","eHccpar","eHbbint","eHbbpar",
    "eggFHgaga","eggFHZga","eggFHZZ","eggFHWW","eggFHtautau","eggFHbb","eggFHmumu",   
    "eVBFHgaga","eVBFHZga","eVBFHZZ","eVBFHWW","eVBFHtautau","eVBFHbb","eVBFHmumu",   
    "eWHgaga","eWHZga","eWHZZ","eWHWW","eWHtautau","eWHbb","eWHmumu",    
    "eZHgaga","eZHZga","eZHZZ","eZHWW","eZHtautau","eZHbb","eZHmumu",
    "ettHgaga","ettHZga","ettHZZ","ettHWW","ettHtautau","ettHbb","ettHmumu",
    "eVBFHinv","eVHinv",
    "eVBF_2_Hbox", "eVBF_2_HQ1_11", "eVBF_2_Hu_11", "eVBF_2_Hd_11", "eVBF_2_HQ3_11",
    "eVBF_2_HD", "eVBF_2_HB", "eVBF_2_HW", "eVBF_2_HWB", "eVBF_2_HG", "eVBF_2_DHB",
    "eVBF_2_DHW", "eVBF_2_DeltaGF",
    "eVBF_78_Hbox", "eVBF_78_HQ1_11", "eVBF_78_Hu_11", "eVBF_78_Hd_11", "eVBF_78_HQ3_11",
    "eVBF_78_HD", "eVBF_78_HB", "eVBF_78_HW", "eVBF_78_HWB", "eVBF_78_HG", "eVBF_78_DHB",
    "eVBF_78_DHW", "eVBF_78_DeltaGF",
    "eVBF_1314_Hbox", "eVBF_1314_HQ1_11", "eVBF_1314_Hu_11", "eVBF_1314_Hd_11", "eVBF_1314_HQ3_11",
    "eVBF_1314_HD", "eVBF_1314_HB", "eVBF_1314_HW", "eVBF_1314_HWB", "eVBF_1314_HG", "eVBF_1314_DHB",
    "eVBF_1314_DHW", "eVBF_1314_DeltaGF",
    "eWH_2_Hbox", "eWH_2_HQ3_11", "eWH_2_HD", "eWH_2_HW", "eWH_2_HWB", "eWH_2_DHW", "eWH_2_DeltaGF",
    "eWH_78_Hbox", "eWH_78_HQ3_11", "eWH_78_HD", "eWH_78_HW", "eWH_78_HWB", "eWH_78_DHW", "eWH_78_DeltaGF",
    "eWH_1314_Hbox", "eWH_1314_HQ3_11", "eWH_1314_HD", "eWH_1314_HW", "eWH_1314_HWB", "eWH_1314_DHW", "eWH_1314_DeltaGF",
    "eZH_2_Hbox", "eZH_2_HQ1_11", "eZH_2_Hu_11", "eZH_2_Hd_11", "eZH_2_HQ3_11", "eZH_2_HD", "eZH_2_HB", "eZH_2_HW", "eZH_2_HWB", "eZH_2_DHB", "eZH_2_DHW", "eZH_2_DeltaGF",
    "eZH_78_Hbox", "eZH_78_HQ1_11", "eZH_78_Hu_11", "eZH_78_Hd_11", "eZH_78_HQ3_11", "eZH_78_HD", "eZH_78_HB", "eZH_78_HW", "eZH_78_HWB", "eZH_78_DHB", "eZH_78_DHW", "eZH_78_DeltaGF",
    "eZH_1314_Hbox", "eZH_1314_HQ1_11", "eZH_1314_Hu_11", "eZH_1314_Hd_11", "eZH_1314_HQ3_11", "eZH_1314_HD", "eZH_1314_HB", "eZH_1314_HW", "eZH_1314_HWB", "eZH_1314_DHB", "eZH_1314_DHW", "eZH_1314_DeltaGF",
    "ettH_2_HG", "ettH_2_G", "ettH_2_uG_33r", "ettH_2_DeltagHt",
    "ettH_78_HG", "ettH_78_G", "ettH_78_uG_33r", "ettH_78_DeltagHt",
    "ettH_1314_HG", "ettH_1314_G", "ettH_1314_uG_33r", "ettH_1314_DeltagHt"};

const std::string NPSMEFTd6::NPSMEFTd6Vars_LFU_QFU[NNPSMEFTd6Vars_LFU_QFU]
        = {"CG", "CW", "C2B", "C2W", "C2BS", "C2WS", "CHG", "CHW", "CHB", "CDHB", "CDHW", "CDB", "CDW", "CHWB", "CHD", "CT", "CHbox", "CH",
    "CHL1", "CHL3", "CHe", "CHQ1", "CHQ3", "CHu", "CHd", "CHud_r", "CHud_i",
    "CeH_11r", "CeH_22r", "CeH_33r", "CeH_11i", "CeH_22i", "CeH_33i", 
    "CuH_11r", "CuH_22r", "CuH_33r", "CuH_11i", "CuH_22i", "CuH_33i", 
    "CdH_11r", "CdH_22r", "CdH_33r", "CdH_11i", "CdH_22i", "CdH_33i",
    "CuG_r", "CuG_i", "CuW_r", "CuW_i", "CuB_r", "CuB_i",
    "CLL", "CLQ1", "CLQ3",
    "Cee", "Ceu", "Ced", "CLe", "CLu", "CLd", "CQe",
    "Lambda_NP",
    "BrHinv","BrHexo",
    "dg1Z","dKappaga","lambZ",
    "eggFint","eggFpar","ettHint","ettHpar",
    "eVBFint","eVBFpar","eWHint","eWHpar","eZHint","eZHpar",
    "eeeWBFint","eeeWBFpar","eeeZHint","eeeZHpar","eeettHint","eeettHpar",
    "eepWBFint","eepWBFpar","eepZBFint","eepZBFpar",
    "eHggint","eHggpar","eHWWint","eHWWpar","eHZZint","eHZZpar","eHZgaint","eHZgapar",
    "eHgagaint","eHgagapar","eHmumuint","eHmumupar","eHtautauint","eHtautaupar",
    "eHccint","eHccpar","eHbbint","eHbbpar",
    "eggFHgaga","eggFHZga","eggFHZZ","eggFHWW","eggFHtautau","eggFHbb","eggFHmumu",   
    "eVBFHgaga","eVBFHZga","eVBFHZZ","eVBFHWW","eVBFHtautau","eVBFHbb","eVBFHmumu",   
    "eWHgaga","eWHZga","eWHZZ","eWHWW","eWHtautau","eWHbb","eWHmumu",    
    "eZHgaga","eZHZga","eZHZZ","eZHWW","eZHtautau","eZHbb","eZHmumu",
    "ettHgaga","ettHZga","ettHZZ","ettHWW","ettHtautau","ettHbb","ettHmumu",
    "eVBFHinv","eVHinv",
    "eVBF_2_Hbox", "eVBF_2_HQ1_11", "eVBF_2_Hu_11", "eVBF_2_Hd_11", "eVBF_2_HQ3_11",
    "eVBF_2_HD", "eVBF_2_HB", "eVBF_2_HW", "eVBF_2_HWB", "eVBF_2_HG", "eVBF_2_DHB",
    "eVBF_2_DHW", "eVBF_2_DeltaGF",
    "eVBF_78_Hbox", "eVBF_78_HQ1_11", "eVBF_78_Hu_11", "eVBF_78_Hd_11", "eVBF_78_HQ3_11",
    "eVBF_78_HD", "eVBF_78_HB", "eVBF_78_HW", "eVBF_78_HWB", "eVBF_78_HG", "eVBF_78_DHB",
    "eVBF_78_DHW", "eVBF_78_DeltaGF",
    "eVBF_1314_Hbox", "eVBF_1314_HQ1_11", "eVBF_1314_Hu_11", "eVBF_1314_Hd_11", "eVBF_1314_HQ3_11",
    "eVBF_1314_HD", "eVBF_1314_HB", "eVBF_1314_HW", "eVBF_1314_HWB", "eVBF_1314_HG", "eVBF_1314_DHB",
    "eVBF_1314_DHW", "eVBF_1314_DeltaGF",
    "eWH_2_Hbox", "eWH_2_HQ3_11", "eWH_2_HD", "eWH_2_HW", "eWH_2_HWB", "eWH_2_DHW", "eWH_2_DeltaGF",
    "eWH_78_Hbox", "eWH_78_HQ3_11", "eWH_78_HD", "eWH_78_HW", "eWH_78_HWB", "eWH_78_DHW", "eWH_78_DeltaGF",
    "eWH_1314_Hbox", "eWH_1314_HQ3_11", "eWH_1314_HD", "eWH_1314_HW", "eWH_1314_HWB", "eWH_1314_DHW", "eWH_1314_DeltaGF",
    "eZH_2_Hbox", "eZH_2_HQ1_11", "eZH_2_Hu_11", "eZH_2_Hd_11", "eZH_2_HQ3_11", "eZH_2_HD", "eZH_2_HB", "eZH_2_HW", "eZH_2_HWB", "eZH_2_DHB", "eZH_2_DHW", "eZH_2_DeltaGF",
    "eZH_78_Hbox", "eZH_78_HQ1_11", "eZH_78_Hu_11", "eZH_78_Hd_11", "eZH_78_HQ3_11", "eZH_78_HD", "eZH_78_HB", "eZH_78_HW", "eZH_78_HWB", "eZH_78_DHB", "eZH_78_DHW", "eZH_78_DeltaGF",
    "eZH_1314_Hbox", "eZH_1314_HQ1_11", "eZH_1314_Hu_11", "eZH_1314_Hd_11", "eZH_1314_HQ3_11", "eZH_1314_HD", "eZH_1314_HB", "eZH_1314_HW", "eZH_1314_HWB", "eZH_1314_DHB", "eZH_1314_DHW", "eZH_1314_DeltaGF",
    "ettH_2_HG", "ettH_2_G", "ettH_2_uG_33r", "ettH_2_DeltagHt",
    "ettH_78_HG", "ettH_78_G", "ettH_78_uG_33r", "ettH_78_DeltagHt",
    "ettH_1314_HG", "ettH_1314_G", "ettH_1314_uG_33r", "ettH_1314_DeltagHt"};

const std::string NPSMEFTd6::NPSMEFTd6VarsRot_LFU_QFU[NNPSMEFTd6Vars_LFU_QFU]
        = {"CG", "CW", "C2B", "C2W", "C2BS", "C2WS", "CHG", "CHWHB_gaga", "CHWHB_gagaorth", "CDHB", "CDHW", "CDB", "CDW", "CHWB", "CHD", "CT", "CHbox", "CH",
    "CHL1", "CHL3", "CHe", "CHQ1", "CHQ3", "CHu", "CHd", "CHud_r", "CHud_i",
    "CeH_11r", "CeH_22r", "CeH_33r", "CeH_11i", "CeH_22i", "CeH_33i", 
    "CuH_11r", "CuH_22r", "CuH_33r", "CuH_11i", "CuH_22i", "CuH_33i", 
    "CdH_11r", "CdH_22r", "CdH_33r", "CdH_11i", "CdH_22i", "CdH_33i",
    "CuG_r", "CuG_i", "CuW_r", "CuW_i", "CuB_r", "CuB_i",
    "CLL", "CLQ1", "CLQ3",
    "Cee", "Ceu", "Ced", "CLe", "CLu", "CLd", "CQe",
    "Lambda_NP",
    "BrHinv","BrHexo",
    "dg1Z","dKappaga","lambZ",
    "eggFint","eggFpar","ettHint","ettHpar",
    "eVBFint","eVBFpar","eWHint","eWHpar","eZHint","eZHpar",
    "eeeWBFint","eeeWBFpar","eeeZHint","eeeZHpar","eeettHint","eeettHpar",
    "eepWBFint","eepWBFpar","eepZBFint","eepZBFpar",
    "eHggint","eHggpar","eHWWint","eHWWpar","eHZZint","eHZZpar","eHZgaint","eHZgapar",
    "eHgagaint","eHgagapar","eHmumuint","eHmumupar","eHtautauint","eHtautaupar",
    "eHccint","eHccpar","eHbbint","eHbbpar",
    "eggFHgaga","eggFHZga","eggFHZZ","eggFHWW","eggFHtautau","eggFHbb","eggFHmumu",   
    "eVBFHgaga","eVBFHZga","eVBFHZZ","eVBFHWW","eVBFHtautau","eVBFHbb","eVBFHmumu",   
    "eWHgaga","eWHZga","eWHZZ","eWHWW","eWHtautau","eWHbb","eWHmumu",    
    "eZHgaga","eZHZga","eZHZZ","eZHWW","eZHtautau","eZHbb","eZHmumu",
    "ettHgaga","ettHZga","ettHZZ","ettHWW","ettHtautau","ettHbb","ettHmumu",
    "eVBFHinv","eVHinv",
    "eVBF_2_Hbox", "eVBF_2_HQ1_11", "eVBF_2_Hu_11", "eVBF_2_Hd_11", "eVBF_2_HQ3_11",
    "eVBF_2_HD", "eVBF_2_HB", "eVBF_2_HW", "eVBF_2_HWB", "eVBF_2_HG", "eVBF_2_DHB",
    "eVBF_2_DHW", "eVBF_2_DeltaGF",
    "eVBF_78_Hbox", "eVBF_78_HQ1_11", "eVBF_78_Hu_11", "eVBF_78_Hd_11", "eVBF_78_HQ3_11",
    "eVBF_78_HD", "eVBF_78_HB", "eVBF_78_HW", "eVBF_78_HWB", "eVBF_78_HG", "eVBF_78_DHB",
    "eVBF_78_DHW", "eVBF_78_DeltaGF",
    "eVBF_1314_Hbox", "eVBF_1314_HQ1_11", "eVBF_1314_Hu_11", "eVBF_1314_Hd_11", "eVBF_1314_HQ3_11",
    "eVBF_1314_HD", "eVBF_1314_HB", "eVBF_1314_HW", "eVBF_1314_HWB", "eVBF_1314_HG", "eVBF_1314_DHB",
    "eVBF_1314_DHW", "eVBF_1314_DeltaGF",
    "eWH_2_Hbox", "eWH_2_HQ3_11", "eWH_2_HD", "eWH_2_HW", "eWH_2_HWB", "eWH_2_DHW", "eWH_2_DeltaGF",
    "eWH_78_Hbox", "eWH_78_HQ3_11", "eWH_78_HD", "eWH_78_HW", "eWH_78_HWB", "eWH_78_DHW", "eWH_78_DeltaGF",
    "eWH_1314_Hbox", "eWH_1314_HQ3_11", "eWH_1314_HD", "eWH_1314_HW", "eWH_1314_HWB", "eWH_1314_DHW", "eWH_1314_DeltaGF",
    "eZH_2_Hbox", "eZH_2_HQ1_11", "eZH_2_Hu_11", "eZH_2_Hd_11", "eZH_2_HQ3_11", "eZH_2_HD", "eZH_2_HB", "eZH_2_HW", "eZH_2_HWB", "eZH_2_DHB", "eZH_2_DHW", "eZH_2_DeltaGF",
    "eZH_78_Hbox", "eZH_78_HQ1_11", "eZH_78_Hu_11", "eZH_78_Hd_11", "eZH_78_HQ3_11", "eZH_78_HD", "eZH_78_HB", "eZH_78_HW", "eZH_78_HWB", "eZH_78_DHB", "eZH_78_DHW", "eZH_78_DeltaGF",
    "eZH_1314_Hbox", "eZH_1314_HQ1_11", "eZH_1314_Hu_11", "eZH_1314_Hd_11", "eZH_1314_HQ3_11", "eZH_1314_HD", "eZH_1314_HB", "eZH_1314_HW", "eZH_1314_HWB", "eZH_1314_DHB", "eZH_1314_DHW", "eZH_1314_DeltaGF",
    "ettH_2_HG", "ettH_2_G", "ettH_2_uG_33r", "ettH_2_DeltagHt",
    "ettH_78_HG", "ettH_78_G", "ettH_78_uG_33r", "ettH_78_DeltagHt",
    "ettH_1314_HG", "ettH_1314_G", "ettH_1314_uG_33r", "ettH_1314_DeltagHt"};

NPSMEFTd6::NPSMEFTd6(const bool FlagLeptonUniversal_in, const bool FlagQuarkUniversal_in)
: NPbase(), NPSMEFTd6M(*this), FlagLeptonUniversal(FlagLeptonUniversal_in), FlagQuarkUniversal(FlagQuarkUniversal_in)
{
    if ((!FlagLeptonUniversal && FlagQuarkUniversal)
            || (FlagLeptonUniversal && !FlagQuarkUniversal))
        throw std::runtime_error("Invalid arguments for NPSMEFTd6::NPSMEFTd6()");

    FlagQuadraticTerms = false;
    FlagRotateCHWCHB = false;
    FlagPartialQFU = false;
    FlagFlavU3OfX = false;
    FlagUnivOfX = false;
    FlagHiggsSM = false;
    FlagLoopHd6 = false;
    FlagLoopH3d6Quad = false;
    setModelLinearized();
    
    w_WW = gsl_integration_cquad_workspace_alloc(100);
    
    SMM.setObj((StandardModelMatching&) NPSMEFTd6M.getObj());

    ModelParamMap.insert(std::make_pair("CG", std::cref(CG)));
    ModelParamMap.insert(std::make_pair("CW", std::cref(CW)));
    ModelParamMap.insert(std::make_pair("C2B", std::cref(C2B)));
    ModelParamMap.insert(std::make_pair("C2W", std::cref(C2W)));
    ModelParamMap.insert(std::make_pair("C2BS", std::cref(C2BS)));
    ModelParamMap.insert(std::make_pair("C2WS", std::cref(C2WS)));
    ModelParamMap.insert(std::make_pair("CHG", std::cref(CHG)));
    ModelParamMap.insert(std::make_pair("CHW", std::cref(CHW)));
    ModelParamMap.insert(std::make_pair("CHB", std::cref(CHB)));
    ModelParamMap.insert(std::make_pair("CHWHB_gaga", std::cref(CHWHB_gaga)));
    ModelParamMap.insert(std::make_pair("CHWHB_gagaorth", std::cref(CHWHB_gagaorth)));        
    ModelParamMap.insert(std::make_pair("CDHB", std::cref(CDHB)));
    ModelParamMap.insert(std::make_pair("CDHW", std::cref(CDHW)));
    ModelParamMap.insert(std::make_pair("CDB", std::cref(CDB)));
    ModelParamMap.insert(std::make_pair("CDW", std::cref(CDW)));
    ModelParamMap.insert(std::make_pair("CHWB", std::cref(CHWB)));
    ModelParamMap.insert(std::make_pair("CHD", std::cref(CHD)));
    ModelParamMap.insert(std::make_pair("CT", std::cref(CT)));
    ModelParamMap.insert(std::make_pair("CHbox", std::cref(CHbox)));
    ModelParamMap.insert(std::make_pair("CH", std::cref(CH)));
    if (FlagLeptonUniversal) {
        ModelParamMap.insert(std::make_pair("CHL1", std::cref(CHL1_11)));
        ModelParamMap.insert(std::make_pair("CHL3", std::cref(CHL3_11)));
        ModelParamMap.insert(std::make_pair("CHe", std::cref(CHe_11)));
        ModelParamMap.insert(std::make_pair("CeH_11r", std::cref(CeH_11r)));
        ModelParamMap.insert(std::make_pair("CeH_22r", std::cref(CeH_22r)));
        ModelParamMap.insert(std::make_pair("CeH_33r", std::cref(CeH_33r)));
        ModelParamMap.insert(std::make_pair("CeH_11i", std::cref(CeH_11i)));
        ModelParamMap.insert(std::make_pair("CeH_22i", std::cref(CeH_22i)));
        ModelParamMap.insert(std::make_pair("CeH_33i", std::cref(CeH_33i)));
        ModelParamMap.insert(std::make_pair("CLL", std::cref(CLL_1221)));
        ModelParamMap.insert(std::make_pair("Cee", std::cref(Cee_1111)));
        ModelParamMap.insert(std::make_pair("CLe", std::cref(CLe_1111)));
    } else {
        ModelParamMap.insert(std::make_pair("CHL1_11", std::cref(CHL1_11)));
        ModelParamMap.insert(std::make_pair("CHL1_12r", std::cref(CHL1_12r)));
        ModelParamMap.insert(std::make_pair("CHL1_13r", std::cref(CHL1_13r)));
        ModelParamMap.insert(std::make_pair("CHL1_22", std::cref(CHL1_22)));
        ModelParamMap.insert(std::make_pair("CHL1_23r", std::cref(CHL1_23r)));
        ModelParamMap.insert(std::make_pair("CHL1_33", std::cref(CHL1_33)));
        ModelParamMap.insert(std::make_pair("CHL1_12i", std::cref(CHL1_12i)));
        ModelParamMap.insert(std::make_pair("CHL1_13i", std::cref(CHL1_13i)));
        ModelParamMap.insert(std::make_pair("CHL1_23i", std::cref(CHL1_23i)));
        ModelParamMap.insert(std::make_pair("CHL3_11", std::cref(CHL3_11)));
        ModelParamMap.insert(std::make_pair("CHL3_12r", std::cref(CHL3_12r)));
        ModelParamMap.insert(std::make_pair("CHL3_13r", std::cref(CHL3_13r)));
        ModelParamMap.insert(std::make_pair("CHL3_22", std::cref(CHL3_22)));
        ModelParamMap.insert(std::make_pair("CHL3_23r", std::cref(CHL3_23r)));
        ModelParamMap.insert(std::make_pair("CHL3_33", std::cref(CHL3_33)));
        ModelParamMap.insert(std::make_pair("CHL3_12i", std::cref(CHL3_12i)));
        ModelParamMap.insert(std::make_pair("CHL3_13i", std::cref(CHL3_13i)));
        ModelParamMap.insert(std::make_pair("CHL3_23i", std::cref(CHL3_23i)));
        ModelParamMap.insert(std::make_pair("CHe_11", std::cref(CHe_11)));
        ModelParamMap.insert(std::make_pair("CHe_12r", std::cref(CHe_12r)));
        ModelParamMap.insert(std::make_pair("CHe_13r", std::cref(CHe_13r)));
        ModelParamMap.insert(std::make_pair("CHe_22", std::cref(CHe_22)));
        ModelParamMap.insert(std::make_pair("CHe_23r", std::cref(CHe_23r)));
        ModelParamMap.insert(std::make_pair("CHe_33", std::cref(CHe_33)));
        ModelParamMap.insert(std::make_pair("CHe_12i", std::cref(CHe_12i)));
        ModelParamMap.insert(std::make_pair("CHe_13i", std::cref(CHe_13i)));
        ModelParamMap.insert(std::make_pair("CHe_23i", std::cref(CHe_23i)));
        ModelParamMap.insert(std::make_pair("CeH_11r", std::cref(CeH_11r)));
        ModelParamMap.insert(std::make_pair("CeH_12r", std::cref(CeH_12r)));
        ModelParamMap.insert(std::make_pair("CeH_13r", std::cref(CeH_13r)));
        ModelParamMap.insert(std::make_pair("CeH_22r", std::cref(CeH_22r)));
        ModelParamMap.insert(std::make_pair("CeH_23r", std::cref(CeH_23r)));
        ModelParamMap.insert(std::make_pair("CeH_33r", std::cref(CeH_33r)));
        ModelParamMap.insert(std::make_pair("CeH_11i", std::cref(CeH_11i)));
        ModelParamMap.insert(std::make_pair("CeH_12i", std::cref(CeH_12i)));
        ModelParamMap.insert(std::make_pair("CeH_13i", std::cref(CeH_13i)));
        ModelParamMap.insert(std::make_pair("CeH_22i", std::cref(CeH_22i)));
        ModelParamMap.insert(std::make_pair("CeH_23i", std::cref(CeH_23i)));
        ModelParamMap.insert(std::make_pair("CeH_33i", std::cref(CeH_33i)));
        ModelParamMap.insert(std::make_pair("CLL_1111", std::cref(CLL_1111)));
        ModelParamMap.insert(std::make_pair("CLL_1221", std::cref(CLL_1221)));
        ModelParamMap.insert(std::make_pair("CLL_1122", std::cref(CLL_1122)));
        ModelParamMap.insert(std::make_pair("CLL_1331", std::cref(CLL_1331)));
        ModelParamMap.insert(std::make_pair("CLL_1133", std::cref(CLL_1133)));
        ModelParamMap.insert(std::make_pair("Cee_1111", std::cref(Cee_1111)));
        ModelParamMap.insert(std::make_pair("Cee_1122", std::cref(Cee_1122)));
        ModelParamMap.insert(std::make_pair("Cee_1133", std::cref(Cee_1133)));
        ModelParamMap.insert(std::make_pair("CLe_1111", std::cref(CLe_1111)));
        ModelParamMap.insert(std::make_pair("CLe_1122", std::cref(CLe_1122)));
        ModelParamMap.insert(std::make_pair("CLe_2211", std::cref(CLe_2211)));
        ModelParamMap.insert(std::make_pair("CLe_1133", std::cref(CLe_1133)));
        ModelParamMap.insert(std::make_pair("CLe_3311", std::cref(CLe_3311)));       
    }
    if (FlagQuarkUniversal) {
        ModelParamMap.insert(std::make_pair("CHQ1", std::cref(CHQ1_11)));
        ModelParamMap.insert(std::make_pair("CHQ3", std::cref(CHQ3_11)));
        ModelParamMap.insert(std::make_pair("CHu", std::cref(CHu_11)));
        ModelParamMap.insert(std::make_pair("CHd", std::cref(CHd_11)));
        ModelParamMap.insert(std::make_pair("CHud_r", std::cref(CHud_11r)));
        ModelParamMap.insert(std::make_pair("CHud_i", std::cref(CHud_11i)));
        ModelParamMap.insert(std::make_pair("CuH_11r", std::cref(CuH_11r)));
        ModelParamMap.insert(std::make_pair("CuH_22r", std::cref(CuH_22r)));
        ModelParamMap.insert(std::make_pair("CuH_33r", std::cref(CuH_33r)));
        ModelParamMap.insert(std::make_pair("CuH_11i", std::cref(CuH_11i)));
        ModelParamMap.insert(std::make_pair("CuH_22i", std::cref(CuH_22i)));
        ModelParamMap.insert(std::make_pair("CuH_33i", std::cref(CuH_33i)));        
        ModelParamMap.insert(std::make_pair("CdH_11r", std::cref(CdH_11r)));
        ModelParamMap.insert(std::make_pair("CdH_22r", std::cref(CdH_22r)));
        ModelParamMap.insert(std::make_pair("CdH_33r", std::cref(CdH_33r)));        
        ModelParamMap.insert(std::make_pair("CdH_11i", std::cref(CdH_11i)));
        ModelParamMap.insert(std::make_pair("CdH_22i", std::cref(CdH_22i)));
        ModelParamMap.insert(std::make_pair("CdH_33i", std::cref(CdH_33i)));        
        ModelParamMap.insert(std::make_pair("CuG_r", std::cref(CuG_11r)));
        ModelParamMap.insert(std::make_pair("CuG_i", std::cref(CuG_11i)));
        ModelParamMap.insert(std::make_pair("CuW_r", std::cref(CuW_11r)));
        ModelParamMap.insert(std::make_pair("CuW_i", std::cref(CuW_11i)));
        ModelParamMap.insert(std::make_pair("CuB_r", std::cref(CuB_11r)));
        ModelParamMap.insert(std::make_pair("CuB_i", std::cref(CuB_11i)));
    } else {
        ModelParamMap.insert(std::make_pair("CHQ1_11", std::cref(CHQ1_11)));
        ModelParamMap.insert(std::make_pair("CHQ1_12r", std::cref(CHQ1_12r)));
        ModelParamMap.insert(std::make_pair("CHQ1_13r", std::cref(CHQ1_13r)));
        ModelParamMap.insert(std::make_pair("CHQ1_22", std::cref(CHQ1_22)));
        ModelParamMap.insert(std::make_pair("CHQ1_23r", std::cref(CHQ1_23r)));
        ModelParamMap.insert(std::make_pair("CHQ1_33", std::cref(CHQ1_33)));
        ModelParamMap.insert(std::make_pair("CHQ1_12i", std::cref(CHQ1_12i)));
        ModelParamMap.insert(std::make_pair("CHQ1_13i", std::cref(CHQ1_13i)));
        ModelParamMap.insert(std::make_pair("CHQ1_23i", std::cref(CHQ1_23i)));
        ModelParamMap.insert(std::make_pair("CHQ3_11", std::cref(CHQ3_11)));
        ModelParamMap.insert(std::make_pair("CHQ3_12r", std::cref(CHQ3_12r)));
        ModelParamMap.insert(std::make_pair("CHQ3_13r", std::cref(CHQ3_13r)));
        ModelParamMap.insert(std::make_pair("CHQ3_22", std::cref(CHQ3_22)));
        ModelParamMap.insert(std::make_pair("CHQ3_23r", std::cref(CHQ3_23r)));
        ModelParamMap.insert(std::make_pair("CHQ3_33", std::cref(CHQ3_33)));
        ModelParamMap.insert(std::make_pair("CHQ3_12i", std::cref(CHQ3_12i)));
        ModelParamMap.insert(std::make_pair("CHQ3_13i", std::cref(CHQ3_13i)));
        ModelParamMap.insert(std::make_pair("CHQ3_23i", std::cref(CHQ3_23i)));
        ModelParamMap.insert(std::make_pair("CHu_11", std::cref(CHu_11)));
        ModelParamMap.insert(std::make_pair("CHu_12r", std::cref(CHu_12r)));
        ModelParamMap.insert(std::make_pair("CHu_13r", std::cref(CHu_13r)));
        ModelParamMap.insert(std::make_pair("CHu_22", std::cref(CHu_22)));
        ModelParamMap.insert(std::make_pair("CHu_23r", std::cref(CHu_23r)));
        ModelParamMap.insert(std::make_pair("CHu_33", std::cref(CHu_33)));
        ModelParamMap.insert(std::make_pair("CHu_12i", std::cref(CHu_12i)));
        ModelParamMap.insert(std::make_pair("CHu_13i", std::cref(CHu_13i)));
        ModelParamMap.insert(std::make_pair("CHu_23i", std::cref(CHu_23i)));
        ModelParamMap.insert(std::make_pair("CHd_11", std::cref(CHd_11)));
        ModelParamMap.insert(std::make_pair("CHd_12r", std::cref(CHd_12r)));
        ModelParamMap.insert(std::make_pair("CHd_13r", std::cref(CHd_13r)));
        ModelParamMap.insert(std::make_pair("CHd_22", std::cref(CHd_22)));
        ModelParamMap.insert(std::make_pair("CHd_23r", std::cref(CHd_23r)));
        ModelParamMap.insert(std::make_pair("CHd_33", std::cref(CHd_33)));
        ModelParamMap.insert(std::make_pair("CHd_12i", std::cref(CHd_12i)));
        ModelParamMap.insert(std::make_pair("CHd_13i", std::cref(CHd_13i)));
        ModelParamMap.insert(std::make_pair("CHd_23i", std::cref(CHd_23i)));
        ModelParamMap.insert(std::make_pair("CHud_11r", std::cref(CHud_11r)));
        ModelParamMap.insert(std::make_pair("CHud_12r", std::cref(CHud_12r)));
        ModelParamMap.insert(std::make_pair("CHud_13r", std::cref(CHud_13r)));
        ModelParamMap.insert(std::make_pair("CHud_22r", std::cref(CHud_22r)));
        ModelParamMap.insert(std::make_pair("CHud_23r", std::cref(CHud_23r)));
        ModelParamMap.insert(std::make_pair("CHud_33r", std::cref(CHud_33r)));
        ModelParamMap.insert(std::make_pair("CHud_11i", std::cref(CHud_11i)));
        ModelParamMap.insert(std::make_pair("CHud_12i", std::cref(CHud_12i)));
        ModelParamMap.insert(std::make_pair("CHud_13i", std::cref(CHud_13i)));
        ModelParamMap.insert(std::make_pair("CHud_22i", std::cref(CHud_22i)));
        ModelParamMap.insert(std::make_pair("CHud_23i", std::cref(CHud_23i)));
        ModelParamMap.insert(std::make_pair("CHud_33i", std::cref(CHud_33i)));
        ModelParamMap.insert(std::make_pair("CuH_11r", std::cref(CuH_11r)));
        ModelParamMap.insert(std::make_pair("CuH_12r", std::cref(CuH_12r)));
        ModelParamMap.insert(std::make_pair("CuH_13r", std::cref(CuH_13r)));
        ModelParamMap.insert(std::make_pair("CuH_22r", std::cref(CuH_22r)));
        ModelParamMap.insert(std::make_pair("CuH_23r", std::cref(CuH_23r)));
        ModelParamMap.insert(std::make_pair("CuH_33r", std::cref(CuH_33r)));
        ModelParamMap.insert(std::make_pair("CuH_11i", std::cref(CuH_11i)));
        ModelParamMap.insert(std::make_pair("CuH_12i", std::cref(CuH_12i)));
        ModelParamMap.insert(std::make_pair("CuH_13i", std::cref(CuH_13i)));
        ModelParamMap.insert(std::make_pair("CuH_22i", std::cref(CuH_22i)));
        ModelParamMap.insert(std::make_pair("CuH_23i", std::cref(CuH_23i)));
        ModelParamMap.insert(std::make_pair("CuH_33i", std::cref(CuH_33i)));
        ModelParamMap.insert(std::make_pair("CdH_11r", std::cref(CdH_11r)));
        ModelParamMap.insert(std::make_pair("CdH_12r", std::cref(CdH_12r)));
        ModelParamMap.insert(std::make_pair("CdH_13r", std::cref(CdH_13r)));
        ModelParamMap.insert(std::make_pair("CdH_22r", std::cref(CdH_22r)));
        ModelParamMap.insert(std::make_pair("CdH_23r", std::cref(CdH_23r)));
        ModelParamMap.insert(std::make_pair("CdH_33r", std::cref(CdH_33r)));
        ModelParamMap.insert(std::make_pair("CdH_11i", std::cref(CdH_11i)));
        ModelParamMap.insert(std::make_pair("CdH_12i", std::cref(CdH_12i)));
        ModelParamMap.insert(std::make_pair("CdH_13i", std::cref(CdH_13i)));
        ModelParamMap.insert(std::make_pair("CdH_22i", std::cref(CdH_22i)));
        ModelParamMap.insert(std::make_pair("CdH_23i", std::cref(CdH_23i)));
        ModelParamMap.insert(std::make_pair("CdH_33i", std::cref(CdH_33i)));
        ModelParamMap.insert(std::make_pair("CuG_11r", std::cref(CuG_11r)));
        ModelParamMap.insert(std::make_pair("CuG_12r", std::cref(CuG_12r)));
        ModelParamMap.insert(std::make_pair("CuG_13r", std::cref(CuG_13r)));
        ModelParamMap.insert(std::make_pair("CuG_22r", std::cref(CuG_22r)));
        ModelParamMap.insert(std::make_pair("CuG_23r", std::cref(CuG_23r)));
        ModelParamMap.insert(std::make_pair("CuG_33r", std::cref(CuG_33r)));
        ModelParamMap.insert(std::make_pair("CuG_11i", std::cref(CuG_11i)));
        ModelParamMap.insert(std::make_pair("CuG_12i", std::cref(CuG_12i)));
        ModelParamMap.insert(std::make_pair("CuG_13i", std::cref(CuG_13i)));
        ModelParamMap.insert(std::make_pair("CuG_22i", std::cref(CuG_22i)));
        ModelParamMap.insert(std::make_pair("CuG_23i", std::cref(CuG_23i)));
        ModelParamMap.insert(std::make_pair("CuG_33i", std::cref(CuG_33i)));
        ModelParamMap.insert(std::make_pair("CuW_11r", std::cref(CuW_11r)));
        ModelParamMap.insert(std::make_pair("CuW_12r", std::cref(CuW_12r)));
        ModelParamMap.insert(std::make_pair("CuW_13r", std::cref(CuW_13r)));
        ModelParamMap.insert(std::make_pair("CuW_22r", std::cref(CuW_22r)));
        ModelParamMap.insert(std::make_pair("CuW_23r", std::cref(CuW_23r)));
        ModelParamMap.insert(std::make_pair("CuW_33r", std::cref(CuW_33r)));
        ModelParamMap.insert(std::make_pair("CuW_11i", std::cref(CuW_11i)));
        ModelParamMap.insert(std::make_pair("CuW_12i", std::cref(CuW_12i)));
        ModelParamMap.insert(std::make_pair("CuW_13i", std::cref(CuW_13i)));
        ModelParamMap.insert(std::make_pair("CuW_22i", std::cref(CuW_22i)));
        ModelParamMap.insert(std::make_pair("CuW_23i", std::cref(CuW_23i)));
        ModelParamMap.insert(std::make_pair("CuW_33i", std::cref(CuW_33i)));
        ModelParamMap.insert(std::make_pair("CuB_11r", std::cref(CuB_11r)));
        ModelParamMap.insert(std::make_pair("CuB_12r", std::cref(CuB_12r)));
        ModelParamMap.insert(std::make_pair("CuB_13r", std::cref(CuB_13r)));
        ModelParamMap.insert(std::make_pair("CuB_22r", std::cref(CuB_22r)));
        ModelParamMap.insert(std::make_pair("CuB_23r", std::cref(CuB_23r)));
        ModelParamMap.insert(std::make_pair("CuB_33r", std::cref(CuB_33r)));
        ModelParamMap.insert(std::make_pair("CuB_11i", std::cref(CuB_11i)));
        ModelParamMap.insert(std::make_pair("CuB_12i", std::cref(CuB_12i)));
        ModelParamMap.insert(std::make_pair("CuB_13i", std::cref(CuB_13i)));
        ModelParamMap.insert(std::make_pair("CuB_22i", std::cref(CuB_22i)));
        ModelParamMap.insert(std::make_pair("CuB_23i", std::cref(CuB_23i)));
        ModelParamMap.insert(std::make_pair("CuB_33i", std::cref(CuB_33i)));
    }
    if(FlagLeptonUniversal && FlagQuarkUniversal){
        ModelParamMap.insert(std::make_pair("CLQ1", std::cref(CLQ1_1111)));  
        ModelParamMap.insert(std::make_pair("CLQ3", std::cref(CLQ3_1111))); 
        ModelParamMap.insert(std::make_pair("Ceu", std::cref(Ceu_1111)));
        ModelParamMap.insert(std::make_pair("Ced", std::cref(Ced_1111)));
        ModelParamMap.insert(std::make_pair("CLu", std::cref(CLu_1111)));
        ModelParamMap.insert(std::make_pair("CLd", std::cref(CLd_1111)));
        ModelParamMap.insert(std::make_pair("CQe", std::cref(CQe_1111)));
    } else {
        ModelParamMap.insert(std::make_pair("CLQ1_1111", std::cref(CLQ1_1111)));
        ModelParamMap.insert(std::make_pair("CLQ1_1122", std::cref(CLQ1_1122)));
        ModelParamMap.insert(std::make_pair("CLQ1_2211", std::cref(CLQ1_2211)));
        ModelParamMap.insert(std::make_pair("CLQ1_1221", std::cref(CLQ1_1221)));
        ModelParamMap.insert(std::make_pair("CLQ1_2112", std::cref(CLQ1_2112)));
        ModelParamMap.insert(std::make_pair("CLQ1_1133", std::cref(CLQ1_1133)));
        ModelParamMap.insert(std::make_pair("CLQ1_3311", std::cref(CLQ1_3311)));
        ModelParamMap.insert(std::make_pair("CLQ1_1331", std::cref(CLQ1_1331)));
        ModelParamMap.insert(std::make_pair("CLQ1_3113", std::cref(CLQ1_3113)));
        ModelParamMap.insert(std::make_pair("CLQ1_1123", std::cref(CLQ1_1123)));
        ModelParamMap.insert(std::make_pair("CLQ1_2223", std::cref(CLQ1_2223)));
        ModelParamMap.insert(std::make_pair("CLQ1_3323", std::cref(CLQ1_3323)));
        ModelParamMap.insert(std::make_pair("CLQ1_1132", std::cref(CLQ1_1132)));
        ModelParamMap.insert(std::make_pair("CLQ1_2232", std::cref(CLQ1_2232)));       
        ModelParamMap.insert(std::make_pair("CLQ1_3332", std::cref(CLQ1_3332)));
        ModelParamMap.insert(std::make_pair("CLQ3_1111", std::cref(CLQ3_1111)));
        ModelParamMap.insert(std::make_pair("CLQ3_1122", std::cref(CLQ3_1122)));
        ModelParamMap.insert(std::make_pair("CLQ3_2211", std::cref(CLQ3_2211)));
        ModelParamMap.insert(std::make_pair("CLQ3_1221", std::cref(CLQ3_1221)));
        ModelParamMap.insert(std::make_pair("CLQ3_2112", std::cref(CLQ3_2112)));
        ModelParamMap.insert(std::make_pair("CLQ3_1133", std::cref(CLQ3_1133)));
        ModelParamMap.insert(std::make_pair("CLQ3_3311", std::cref(CLQ3_3311)));
        ModelParamMap.insert(std::make_pair("CLQ3_1331", std::cref(CLQ3_1331)));
        ModelParamMap.insert(std::make_pair("CLQ3_3113", std::cref(CLQ3_3113)));
        ModelParamMap.insert(std::make_pair("CLQ3_1123", std::cref(CLQ3_1123)));
        ModelParamMap.insert(std::make_pair("CLQ3_2223", std::cref(CLQ3_2223)));
        ModelParamMap.insert(std::make_pair("CLQ3_3323", std::cref(CLQ3_3323)));
        ModelParamMap.insert(std::make_pair("CLQ3_1132", std::cref(CLQ3_1132)));
        ModelParamMap.insert(std::make_pair("CLQ3_2232", std::cref(CLQ3_2232)));       
        ModelParamMap.insert(std::make_pair("CLQ3_3332", std::cref(CLQ3_3332)));
        ModelParamMap.insert(std::make_pair("Ceu_1111", std::cref(Ceu_1111)));
        ModelParamMap.insert(std::make_pair("Ceu_1122", std::cref(Ceu_1122)));
        ModelParamMap.insert(std::make_pair("Ceu_2211", std::cref(Ceu_2211)));
        ModelParamMap.insert(std::make_pair("Ceu_1133", std::cref(Ceu_1133)));
        ModelParamMap.insert(std::make_pair("Ceu_2233", std::cref(Ceu_2233)));
        ModelParamMap.insert(std::make_pair("Ceu_3311", std::cref(Ceu_3311)));
        ModelParamMap.insert(std::make_pair("Ced_1111", std::cref(Ced_1111)));
        ModelParamMap.insert(std::make_pair("Ced_1122", std::cref(Ced_1122)));
        ModelParamMap.insert(std::make_pair("Ced_2211", std::cref(Ced_2211)));
        ModelParamMap.insert(std::make_pair("Ced_1133", std::cref(Ced_1133)));
        ModelParamMap.insert(std::make_pair("Ced_3311", std::cref(Ced_3311)));
        ModelParamMap.insert(std::make_pair("Ced_1123", std::cref(Ced_1123)));
        ModelParamMap.insert(std::make_pair("Ced_2223", std::cref(Ced_2223)));
        ModelParamMap.insert(std::make_pair("Ced_3323", std::cref(Ced_3323)));
        ModelParamMap.insert(std::make_pair("Ced_1132", std::cref(Ced_1132)));
        ModelParamMap.insert(std::make_pair("Ced_2232", std::cref(Ced_2232)));       
        ModelParamMap.insert(std::make_pair("Ced_3332", std::cref(Ced_3332)));
        ModelParamMap.insert(std::make_pair("CLu_1111", std::cref(CLu_1111)));
        ModelParamMap.insert(std::make_pair("CLu_1122", std::cref(CLu_1122)));
        ModelParamMap.insert(std::make_pair("CLu_2211", std::cref(CLu_2211)));
        ModelParamMap.insert(std::make_pair("CLu_1133", std::cref(CLu_1133)));
        ModelParamMap.insert(std::make_pair("CLu_2233", std::cref(CLu_2233)));
        ModelParamMap.insert(std::make_pair("CLu_3311", std::cref(CLu_3311)));
        ModelParamMap.insert(std::make_pair("CLd_1111", std::cref(CLd_1111)));
        ModelParamMap.insert(std::make_pair("CLd_1122", std::cref(CLd_1122)));
        ModelParamMap.insert(std::make_pair("CLd_2211", std::cref(CLd_2211)));
        ModelParamMap.insert(std::make_pair("CLd_1133", std::cref(CLd_1133)));
        ModelParamMap.insert(std::make_pair("CLd_3311", std::cref(CLd_3311)));
        ModelParamMap.insert(std::make_pair("CLd_1123", std::cref(CLd_1123)));
        ModelParamMap.insert(std::make_pair("CLd_2223", std::cref(CLd_2223)));
        ModelParamMap.insert(std::make_pair("CLd_3323", std::cref(CLd_3323)));
        ModelParamMap.insert(std::make_pair("CLd_1132", std::cref(CLd_1132)));
        ModelParamMap.insert(std::make_pair("CLd_2232", std::cref(CLd_2232)));       
        ModelParamMap.insert(std::make_pair("CLd_3332", std::cref(CLd_3332)));
        ModelParamMap.insert(std::make_pair("CQe_1111", std::cref(CQe_1111)));
        ModelParamMap.insert(std::make_pair("CQe_1122", std::cref(CQe_1122)));
        ModelParamMap.insert(std::make_pair("CQe_2211", std::cref(CQe_2211)));
        ModelParamMap.insert(std::make_pair("CQe_1133", std::cref(CQe_1133)));
        ModelParamMap.insert(std::make_pair("CQe_3311", std::cref(CQe_3311)));
        ModelParamMap.insert(std::make_pair("CQe_2311", std::cref(CQe_2311)));
        ModelParamMap.insert(std::make_pair("CQe_2322", std::cref(CQe_2322)));
        ModelParamMap.insert(std::make_pair("CQe_2333", std::cref(CQe_2333)));
        ModelParamMap.insert(std::make_pair("CQe_3211", std::cref(CQe_3211)));
        ModelParamMap.insert(std::make_pair("CQe_3222", std::cref(CQe_3222)));
        ModelParamMap.insert(std::make_pair("CQe_3233", std::cref(CQe_3233)));
        ModelParamMap.insert(std::make_pair("CLedQ_11", std::cref(CLedQ_11)));
        ModelParamMap.insert(std::make_pair("CLedQ_22", std::cref(CLedQ_22)));
        ModelParamMap.insert(std::make_pair("CpLedQ_11", std::cref(CpLedQ_11)));
        ModelParamMap.insert(std::make_pair("CpLedQ_22", std::cref(CpLedQ_22)));
    }
    ModelParamMap.insert(std::make_pair("Lambda_NP", std::cref(Lambda_NP)));
    ModelParamMap.insert(std::make_pair("BrHinv", std::cref(BrHinv)));
    ModelParamMap.insert(std::make_pair("BrHexo", std::cref(BrHexo)));
    ModelParamMap.insert(std::make_pair("dg1Z", std::cref(dg1Z)));
    ModelParamMap.insert(std::make_pair("dKappaga", std::cref(dKappaga)));
    ModelParamMap.insert(std::make_pair("lambZ", std::cref(lambZ)));
    ModelParamMap.insert(std::make_pair("eggFint", std::cref(eggFint)));
    ModelParamMap.insert(std::make_pair("eggFpar", std::cref(eggFpar)));
    ModelParamMap.insert(std::make_pair("ettHint", std::cref(ettHint)));
    ModelParamMap.insert(std::make_pair("ettHpar", std::cref(ettHpar)));
    ModelParamMap.insert(std::make_pair("eVBFint", std::cref(eVBFint)));
    ModelParamMap.insert(std::make_pair("eVBFpar", std::cref(eVBFpar)));
    ModelParamMap.insert(std::make_pair("eWHint", std::cref(eWHint)));
    ModelParamMap.insert(std::make_pair("eWHpar", std::cref(eWHpar)));
    ModelParamMap.insert(std::make_pair("eZHint", std::cref(eZHint)));
    ModelParamMap.insert(std::make_pair("eZHpar", std::cref(eZHpar)));
    ModelParamMap.insert(std::make_pair("eeeWBFint", std::cref(eeeWBFint)));
    ModelParamMap.insert(std::make_pair("eeeWBFpar", std::cref(eeeWBFpar)));
    ModelParamMap.insert(std::make_pair("eeeZHint", std::cref(eeeZHint)));
    ModelParamMap.insert(std::make_pair("eeeZHpar", std::cref(eeeZHpar)));
    ModelParamMap.insert(std::make_pair("eeettHint", std::cref(eeettHint)));
    ModelParamMap.insert(std::make_pair("eeettHpar", std::cref(eeettHpar)));    
    ModelParamMap.insert(std::make_pair("eepWBFint", std::cref(eepWBFint)));
    ModelParamMap.insert(std::make_pair("eepWBFpar", std::cref(eepWBFpar)));
    ModelParamMap.insert(std::make_pair("eepZBFint", std::cref(eepZBFint)));
    ModelParamMap.insert(std::make_pair("eepZBFpar", std::cref(eepZBFpar)));        
    ModelParamMap.insert(std::make_pair("eHggint", std::cref(eHggint)));
    ModelParamMap.insert(std::make_pair("eHggpar", std::cref(eHggpar)));
    ModelParamMap.insert(std::make_pair("eHWWint", std::cref(eHWWint)));
    ModelParamMap.insert(std::make_pair("eHWWpar", std::cref(eHWWpar)));
    ModelParamMap.insert(std::make_pair("eHZZint", std::cref(eHZZint)));
    ModelParamMap.insert(std::make_pair("eHZZpar", std::cref(eHZZpar)));
    ModelParamMap.insert(std::make_pair("eHZgaint", std::cref(eHZgaint)));
    ModelParamMap.insert(std::make_pair("eHZgapar", std::cref(eHZgapar)));
    ModelParamMap.insert(std::make_pair("eHgagaint", std::cref(eHgagaint)));
    ModelParamMap.insert(std::make_pair("eHgagapar", std::cref(eHgagapar)));
    ModelParamMap.insert(std::make_pair("eHmumuint", std::cref(eHmumuint)));
    ModelParamMap.insert(std::make_pair("eHmumupar", std::cref(eHmumupar)));
    ModelParamMap.insert(std::make_pair("eHtautauint", std::cref(eHtautauint)));
    ModelParamMap.insert(std::make_pair("eHtautaupar", std::cref(eHtautaupar)));
    ModelParamMap.insert(std::make_pair("eHccint", std::cref(eHccint)));
    ModelParamMap.insert(std::make_pair("eHccpar", std::cref(eHccpar)));
    ModelParamMap.insert(std::make_pair("eHbbint", std::cref(eHbbint)));
    ModelParamMap.insert(std::make_pair("eHbbpar", std::cref(eHbbpar)));
    ModelParamMap.insert(std::make_pair("eggFHgaga", std::cref(eggFHgaga)));
    ModelParamMap.insert(std::make_pair("eggFHZga", std::cref(eggFHZga)));
    ModelParamMap.insert(std::make_pair("eggFHZZ", std::cref(eggFHZZ)));
    ModelParamMap.insert(std::make_pair("eggFHWW", std::cref(eggFHWW)));
    ModelParamMap.insert(std::make_pair("eggFHtautau", std::cref(eggFHtautau)));
    ModelParamMap.insert(std::make_pair("eggFHbb", std::cref(eggFHbb)));
    ModelParamMap.insert(std::make_pair("eggFHmumu", std::cref(eggFHmumu)));   
    ModelParamMap.insert(std::make_pair("eVBFHgaga", std::cref(eVBFHgaga)));
    ModelParamMap.insert(std::make_pair("eVBFHZga", std::cref(eVBFHZga)));
    ModelParamMap.insert(std::make_pair("eVBFHZZ", std::cref(eVBFHZZ)));
    ModelParamMap.insert(std::make_pair("eVBFHWW", std::cref(eVBFHWW)));
    ModelParamMap.insert(std::make_pair("eVBFHtautau", std::cref(eVBFHtautau)));
    ModelParamMap.insert(std::make_pair("eVBFHbb", std::cref(eVBFHbb)));
    ModelParamMap.insert(std::make_pair("eVBFHmumu", std::cref(eVBFHmumu)));   
    ModelParamMap.insert(std::make_pair("eWHgaga", std::cref(eWHgaga)));
    ModelParamMap.insert(std::make_pair("eWHZga", std::cref(eWHZga)));
    ModelParamMap.insert(std::make_pair("eWHZZ", std::cref(eWHZZ)));
    ModelParamMap.insert(std::make_pair("eWHWW", std::cref(eWHWW)));
    ModelParamMap.insert(std::make_pair("eWHtautau", std::cref(eWHtautau)));
    ModelParamMap.insert(std::make_pair("eWHbb", std::cref(eWHbb)));
    ModelParamMap.insert(std::make_pair("eWHmumu", std::cref(eWHmumu)));    
    ModelParamMap.insert(std::make_pair("eZHgaga", std::cref(eZHgaga)));
    ModelParamMap.insert(std::make_pair("eZHZga", std::cref(eZHZga)));
    ModelParamMap.insert(std::make_pair("eZHZZ", std::cref(eZHZZ)));
    ModelParamMap.insert(std::make_pair("eZHWW", std::cref(eZHWW)));
    ModelParamMap.insert(std::make_pair("eZHtautau", std::cref(eZHtautau)));
    ModelParamMap.insert(std::make_pair("eZHbb", std::cref(eZHbb)));
    ModelParamMap.insert(std::make_pair("eZHmumu", std::cref(eZHmumu)));
    ModelParamMap.insert(std::make_pair("ettHgaga", std::cref(ettHgaga)));
    ModelParamMap.insert(std::make_pair("ettHZga", std::cref(ettHZga)));
    ModelParamMap.insert(std::make_pair("ettHZZ", std::cref(ettHZZ)));
    ModelParamMap.insert(std::make_pair("ettHWW", std::cref(ettHWW)));
    ModelParamMap.insert(std::make_pair("ettHtautau", std::cref(ettHtautau)));
    ModelParamMap.insert(std::make_pair("ettHbb", std::cref(ettHbb)));
    ModelParamMap.insert(std::make_pair("ettHmumu", std::cref(ettHmumu)));
    ModelParamMap.insert(std::make_pair("eVBFHinv", std::cref(eVBFHinv)));
    ModelParamMap.insert(std::make_pair("eVHinv", std::cref(eVHinv)));
    ModelParamMap.insert(std::make_pair("eVBF_2_Hbox", std::cref(eVBF_2_Hbox)));
    ModelParamMap.insert(std::make_pair("eVBF_2_HQ1_11", std::cref(eVBF_2_HQ1_11)));
    ModelParamMap.insert(std::make_pair("eVBF_2_Hu_11", std::cref(eVBF_2_Hu_11)));
    ModelParamMap.insert(std::make_pair("eVBF_2_Hd_11", std::cref(eVBF_2_Hd_11)));
    ModelParamMap.insert(std::make_pair("eVBF_2_HQ3_11", std::cref(eVBF_2_HQ3_11)));
    ModelParamMap.insert(std::make_pair("eVBF_2_HD", std::cref(eVBF_2_HD)));
    ModelParamMap.insert(std::make_pair("eVBF_2_HB", std::cref(eVBF_2_HB)));
    ModelParamMap.insert(std::make_pair("eVBF_2_HW", std::cref(eVBF_2_HW)));
    ModelParamMap.insert(std::make_pair("eVBF_2_HWB", std::cref(eVBF_2_HWB)));
    ModelParamMap.insert(std::make_pair("eVBF_2_HG", std::cref(eVBF_2_HG)));
    ModelParamMap.insert(std::make_pair("eVBF_2_DHB", std::cref(eVBF_2_DHB)));
    ModelParamMap.insert(std::make_pair("eVBF_2_DHW", std::cref(eVBF_2_DHW)));
    ModelParamMap.insert(std::make_pair("eVBF_2_DeltaGF", std::cref(eVBF_2_DeltaGF)));
    ModelParamMap.insert(std::make_pair("eVBF_78_Hbox", std::cref(eVBF_78_Hbox)));
    ModelParamMap.insert(std::make_pair("eVBF_78_HQ1_11", std::cref(eVBF_78_HQ1_11)));
    ModelParamMap.insert(std::make_pair("eVBF_78_Hu_11", std::cref(eVBF_78_Hu_11)));
    ModelParamMap.insert(std::make_pair("eVBF_78_Hd_11", std::cref(eVBF_78_Hd_11)));
    ModelParamMap.insert(std::make_pair("eVBF_78_HQ3_11", std::cref(eVBF_78_HQ3_11)));
    ModelParamMap.insert(std::make_pair("eVBF_78_HD", std::cref(eVBF_78_HD)));
    ModelParamMap.insert(std::make_pair("eVBF_78_HB", std::cref(eVBF_78_HB)));
    ModelParamMap.insert(std::make_pair("eVBF_78_HW", std::cref(eVBF_78_HW)));
    ModelParamMap.insert(std::make_pair("eVBF_78_HWB", std::cref(eVBF_78_HWB)));
    ModelParamMap.insert(std::make_pair("eVBF_78_HG", std::cref(eVBF_78_HG)));
    ModelParamMap.insert(std::make_pair("eVBF_78_DHB", std::cref(eVBF_78_DHB)));
    ModelParamMap.insert(std::make_pair("eVBF_78_DHW", std::cref(eVBF_78_DHW)));
    ModelParamMap.insert(std::make_pair("eVBF_78_DeltaGF", std::cref(eVBF_78_DeltaGF)));
    ModelParamMap.insert(std::make_pair("eVBF_1314_Hbox", std::cref(eVBF_1314_Hbox)));
    ModelParamMap.insert(std::make_pair("eVBF_1314_HQ1_11", std::cref(eVBF_1314_HQ1_11)));
    ModelParamMap.insert(std::make_pair("eVBF_1314_Hu_11", std::cref(eVBF_1314_Hu_11)));
    ModelParamMap.insert(std::make_pair("eVBF_1314_Hd_11", std::cref(eVBF_1314_Hd_11)));
    ModelParamMap.insert(std::make_pair("eVBF_1314_HQ3_11", std::cref(eVBF_1314_HQ3_11)));
    ModelParamMap.insert(std::make_pair("eVBF_1314_HD", std::cref(eVBF_1314_HD)));
    ModelParamMap.insert(std::make_pair("eVBF_1314_HB", std::cref(eVBF_1314_HB)));
    ModelParamMap.insert(std::make_pair("eVBF_1314_HW", std::cref(eVBF_1314_HW)));
    ModelParamMap.insert(std::make_pair("eVBF_1314_HWB", std::cref(eVBF_1314_HWB)));
    ModelParamMap.insert(std::make_pair("eVBF_1314_HG", std::cref(eVBF_1314_HG)));
    ModelParamMap.insert(std::make_pair("eVBF_1314_DHB", std::cref(eVBF_1314_DHB)));
    ModelParamMap.insert(std::make_pair("eVBF_1314_DHW", std::cref(eVBF_1314_DHW)));
    ModelParamMap.insert(std::make_pair("eVBF_1314_DeltaGF", std::cref(eVBF_1314_DeltaGF)));
    ModelParamMap.insert(std::make_pair("eWH_2_Hbox", std::cref(eWH_2_Hbox)));
    ModelParamMap.insert(std::make_pair("eWH_2_HQ3_11", std::cref(eWH_2_HQ3_11)));
    ModelParamMap.insert(std::make_pair("eWH_2_HD", std::cref(eWH_2_HD)));
    ModelParamMap.insert(std::make_pair("eWH_2_HW", std::cref(eWH_2_HW)));
    ModelParamMap.insert(std::make_pair("eWH_2_HWB", std::cref(eWH_2_HWB)));
    ModelParamMap.insert(std::make_pair("eWH_2_DHW", std::cref(eWH_2_DHW)));
    ModelParamMap.insert(std::make_pair("eWH_2_DeltaGF", std::cref(eWH_2_DeltaGF)));
    ModelParamMap.insert(std::make_pair("eWH_78_Hbox", std::cref(eWH_78_Hbox)));
    ModelParamMap.insert(std::make_pair("eWH_78_HQ3_11", std::cref(eWH_78_HQ3_11)));
    ModelParamMap.insert(std::make_pair("eWH_78_HD", std::cref(eWH_78_HD)));
    ModelParamMap.insert(std::make_pair("eWH_78_HW", std::cref(eWH_78_HW)));
    ModelParamMap.insert(std::make_pair("eWH_78_HWB", std::cref(eWH_78_HWB)));
    ModelParamMap.insert(std::make_pair("eWH_78_DHW", std::cref(eWH_78_DHW)));
    ModelParamMap.insert(std::make_pair("eWH_78_DeltaGF", std::cref(eWH_78_DeltaGF)));
    ModelParamMap.insert(std::make_pair("eWH_1314_Hbox", std::cref(eWH_1314_Hbox)));
    ModelParamMap.insert(std::make_pair("eWH_1314_HQ3_11", std::cref(eWH_1314_HQ3_11)));
    ModelParamMap.insert(std::make_pair("eWH_1314_HD", std::cref(eWH_1314_HD)));
    ModelParamMap.insert(std::make_pair("eWH_1314_HW", std::cref(eWH_1314_HW)));
    ModelParamMap.insert(std::make_pair("eWH_1314_HWB", std::cref(eWH_1314_HWB)));
    ModelParamMap.insert(std::make_pair("eWH_1314_DHW", std::cref(eWH_1314_DHW)));
    ModelParamMap.insert(std::make_pair("eWH_1314_DeltaGF", std::cref(eWH_1314_DeltaGF)));
    ModelParamMap.insert(std::make_pair("eZH_2_Hbox", std::cref(eZH_2_Hbox)));
    ModelParamMap.insert(std::make_pair("eZH_2_HQ1_11", std::cref(eZH_2_HQ1_11)));
    ModelParamMap.insert(std::make_pair("eZH_2_Hu_11", std::cref(eZH_2_Hu_11)));
    ModelParamMap.insert(std::make_pair("eZH_2_Hd_11", std::cref(eZH_2_Hd_11)));
    ModelParamMap.insert(std::make_pair("eZH_2_HQ3_11", std::cref(eZH_2_HQ3_11)));
    ModelParamMap.insert(std::make_pair("eZH_2_HD", std::cref(eZH_2_HD)));
    ModelParamMap.insert(std::make_pair("eZH_2_HB", std::cref(eZH_2_HB)));
    ModelParamMap.insert(std::make_pair("eZH_2_HW", std::cref(eZH_2_HW)));
    ModelParamMap.insert(std::make_pair("eZH_2_HWB", std::cref(eZH_2_HWB)));
    ModelParamMap.insert(std::make_pair("eZH_2_DHB", std::cref(eZH_2_DHB)));
    ModelParamMap.insert(std::make_pair("eZH_2_DHW", std::cref(eZH_2_DHW)));
    ModelParamMap.insert(std::make_pair("eZH_2_DeltaGF", std::cref(eZH_2_DeltaGF)));
    ModelParamMap.insert(std::make_pair("eZH_78_Hbox", std::cref(eZH_78_Hbox)));
    ModelParamMap.insert(std::make_pair("eZH_78_HQ1_11", std::cref(eZH_78_HQ1_11)));
    ModelParamMap.insert(std::make_pair("eZH_78_Hu_11", std::cref(eZH_78_Hu_11)));
    ModelParamMap.insert(std::make_pair("eZH_78_Hd_11", std::cref(eZH_78_Hd_11)));
    ModelParamMap.insert(std::make_pair("eZH_78_HQ3_11", std::cref(eZH_78_HQ3_11)));
    ModelParamMap.insert(std::make_pair("eZH_78_HD", std::cref(eZH_78_HD)));
    ModelParamMap.insert(std::make_pair("eZH_78_HB", std::cref(eZH_78_HB)));
    ModelParamMap.insert(std::make_pair("eZH_78_HW", std::cref(eZH_78_HW)));
    ModelParamMap.insert(std::make_pair("eZH_78_HWB", std::cref(eZH_78_HWB)));
    ModelParamMap.insert(std::make_pair("eZH_78_DHB", std::cref(eZH_78_DHB)));
    ModelParamMap.insert(std::make_pair("eZH_78_DHW", std::cref(eZH_78_DHW)));
    ModelParamMap.insert(std::make_pair("eZH_78_DeltaGF", std::cref(eZH_78_DeltaGF)));
    ModelParamMap.insert(std::make_pair("eZH_1314_Hbox", std::cref(eZH_1314_Hbox)));
    ModelParamMap.insert(std::make_pair("eZH_1314_HQ1_11", std::cref(eZH_1314_HQ1_11)));
    ModelParamMap.insert(std::make_pair("eZH_1314_Hu_11", std::cref(eZH_1314_Hu_11)));
    ModelParamMap.insert(std::make_pair("eZH_1314_Hd_11", std::cref(eZH_1314_Hd_11)));
    ModelParamMap.insert(std::make_pair("eZH_1314_HQ3_11", std::cref(eZH_1314_HQ3_11)));
    ModelParamMap.insert(std::make_pair("eZH_1314_HD", std::cref(eZH_1314_HD)));
    ModelParamMap.insert(std::make_pair("eZH_1314_HB", std::cref(eZH_1314_HB)));
    ModelParamMap.insert(std::make_pair("eZH_1314_HW", std::cref(eZH_1314_HW)));
    ModelParamMap.insert(std::make_pair("eZH_1314_HWB", std::cref(eZH_1314_HWB)));
    ModelParamMap.insert(std::make_pair("eZH_1314_DHB", std::cref(eZH_1314_DHB)));
    ModelParamMap.insert(std::make_pair("eZH_1314_DHW", std::cref(eZH_1314_DHW)));
    ModelParamMap.insert(std::make_pair("eZH_1314_DeltaGF", std::cref(eZH_1314_DeltaGF)));
    ModelParamMap.insert(std::make_pair("ettH_2_HG", std::cref(ettH_2_HG)));
    ModelParamMap.insert(std::make_pair("ettH_2_G", std::cref(ettH_2_G)));
    ModelParamMap.insert(std::make_pair("ettH_2_uG_33r", std::cref(ettH_2_uG_33r)));
    ModelParamMap.insert(std::make_pair("ettH_2_DeltagHt", std::cref(ettH_2_DeltagHt)));
    ModelParamMap.insert(std::make_pair("ettH_78_HG", std::cref(ettH_78_HG)));
    ModelParamMap.insert(std::make_pair("ettH_78_G", std::cref(ettH_78_G)));
    ModelParamMap.insert(std::make_pair("ettH_78_uG_33r", std::cref(ettH_78_uG_33r)));
    ModelParamMap.insert(std::make_pair("ettH_78_DeltagHt", std::cref(ettH_78_DeltagHt)));
    ModelParamMap.insert(std::make_pair("ettH_1314_HG", std::cref(ettH_1314_HG)));
    ModelParamMap.insert(std::make_pair("ettH_1314_G", std::cref(ettH_1314_G)));
    ModelParamMap.insert(std::make_pair("ettH_1314_uG_33r", std::cref(ettH_1314_uG_33r)));
    ModelParamMap.insert(std::make_pair("ettH_1314_DeltagHt", std::cref(ettH_1314_DeltagHt)));
    
    if (FlagLeptonUniversal) {
        CeH_12r = 0.0;
        CeH_13r = 0.0;
        CeH_23r = 0.0;
        CeH_12i = 0.0;
        CeH_13i = 0.0;
        CeH_23i = 0.0;  
        
//  bsll/sbll entries only interesting (for the moment) if non-lepton universal. Set to 0 otherwise
        CLQ1_1123 = 0.0;
        CLQ1_2223 = 0.0;
        CLQ1_3323 = 0.0;
        CLQ1_1132 = 0.0;
        CLQ1_2232 = 0.0;
        CLQ1_3332 = 0.0;
        
        CLQ3_1123 = 0.0;
        CLQ3_2223 = 0.0;
        CLQ3_3323 = 0.0;
        CLQ3_1132 = 0.0;
        CLQ3_2232 = 0.0;
        CLQ3_3332 = 0.0;
        
        Ced_1123 = 0.0;
        Ced_2223 = 0.0;
        Ced_3323 = 0.0;
        Ced_1132 = 0.0;
        Ced_2232 = 0.0;
        Ced_3332 = 0.0;  
        
        CLd_1123 = 0.0;
        CLd_2223 = 0.0;
        CLd_3323 = 0.0;
        CLd_1132 = 0.0;
        CLd_2232 = 0.0;
        CLd_3332 = 0.0;
        
        CQe_2311 = 0.0;
        CQe_2322 = 0.0;
        CQe_2333 = 0.0;
        CQe_3211 = 0.0;
        CQe_3222 = 0.0;
        CQe_3233 = 0.0;
    }
    if (FlagQuarkUniversal) {
        CuH_12r = 0.0;
        CuH_13r = 0.0;
        CuH_23r = 0.0;
        CuH_12i = 0.0;
        CuH_13i = 0.0;
        CuH_23i = 0.0;

        CdH_12r = 0.0;
        CdH_13r = 0.0;
        CdH_23r = 0.0;
        CdH_12i = 0.0;
        CdH_13i = 0.0;
        CdH_23i = 0.0;        
    }
}

bool NPSMEFTd6::PostUpdate()
{
    if (!NPbase::PostUpdate()) return (false);

//  0) Post-update operations not involving SM nor NP parameters    
    if (!FlagHiggsSM) {
        cHSM = 0.0;
    } else {
        cHSM = 1.0;
    }
    
    if (!FlagLoopHd6) {
        cLHd6 = 0.0;
    } else {
        cLHd6 = 1.0;
    }
    
    if (FlagLoopH3d6Quad || FlagQuadraticTerms) {
        cLH3d62 = 1.0;
    } else {
        cLH3d62 = 0.0;
    }
   
//  1) Post-update operations involving SM parameters only
    LambdaNP2 = Lambda_NP * Lambda_NP;
    v2 = v() * v();
    v2_over_LambdaNP2 = v2 / LambdaNP2;
    aleMz = alphaMz();
    eeMz = sqrt( 4.0 * M_PI * aleMz );
    eeMz2 = eeMz*eeMz;
    cW_tree = Mw_tree() / Mz;
    cW2_tree = cW_tree * cW_tree;
    sW2_tree = 1.0 - cW2_tree;
    sW_tree = sqrt(sW2_tree);
    
    g1_tree = eeMz/cW_tree;
    g2_tree = eeMz/sW_tree;
    g3_tree = sqrt( 4.0 * M_PI * AlsMz );
    
    lambdaH_tree = mHl*mHl/2.0/v2;
    
    gZvL = (leptons[NEUTRINO_1].getIsospin());
    gZlL = (leptons[ELECTRON].getIsospin()) - (leptons[ELECTRON].getCharge())*sW2_tree;
    gZlR = - (leptons[ELECTRON].getCharge()) * sW2_tree;
    gZuL = (quarks[UP].getIsospin()) - (quarks[UP].getCharge())*sW2_tree;
    gZuR = - (quarks[UP].getCharge()) * sW2_tree;
    gZdL = (quarks[DOWN].getIsospin()) - (quarks[DOWN].getCharge())*sW2_tree;
    gZdR = - (quarks[DOWN].getCharge()) * sW2_tree;
    
    UevL = 1.0; // Neglect PMNS effects
    VudL = 1.0; // Neglect CKM effects    

    Yuke = sqrt(2.) * (leptons[ELECTRON].getMass()) / v();
    Yukmu = sqrt(2.) * (leptons[MU].getMass()) / v();
    Yuktau = sqrt(2.) * (leptons[TAU].getMass()) / v();
    Yuku = sqrt(2.) * (quarks[UP].getMass()) / v();
    Yukc = sqrt(2.) * (quarks[CHARM].getMass()) / v();
    Yukt = sqrt(2.) * mtpole / v();
    Yukd = sqrt(2.) * (quarks[DOWN].getMass()) / v();
    Yuks = sqrt(2.) * (quarks[STRANGE].getMass()) / v();
    Yukb = sqrt(2.) * (quarks[BOTTOM].getMass()) / v();
    
    dZH = -(9.0/16.0)*( GF*mHl*mHl/sqrt(2.0)/M_PI/M_PI )*( 2.0*M_PI/3.0/sqrt(3.0) - 1.0 );
        
//  2) Post-update operations related to assumptions in the form of the dimension-6 operators 
  
//  Rotated CHW and CHB parameters: Here I need to overwrite the model parameters (There are always 2 on/2 off but need the values of both in output)
    if (FlagRotateCHWCHB) {
        CHW = sW2_tree * CHWHB_gaga - cW2_tree * CHWHB_gagaorth;
        CHB = cW2_tree * CHWHB_gaga + sW2_tree * CHWHB_gagaorth;
    } else {
        CHWHB_gaga = sW2_tree * CHW + cW2_tree * CHB;
        CHWHB_gagaorth = - cW2_tree * CHW + sW2_tree * CHB;
    }    
    
//  Flavour universality assumptions    
    
//  Initialize the internal Wilson coeffs of the form CfH and CfV from the model parameters
    CieH_11r = CeH_11r;
    CieH_22r = CeH_22r;
    CieH_33r = CeH_33r;
    
    CiuH_11r = CuH_11r;
    CiuH_22r = CuH_22r;
    CiuH_33r = CuH_33r;
    
    CidH_11r = CdH_11r;
    CidH_22r = CdH_22r;
    CidH_33r = CdH_33r;
    
    CiuG_11r = CuG_11r;  
    CiuG_22r = CuG_22r;
    CiuG_33r = CuG_33r;

    CiuW_11r = CuW_11r;
    CiuW_22r = CuW_22r;
    CiuW_33r = CuW_33r;

    CiuB_11r = CuB_11r;
    CiuB_22r = CuB_22r;                   
    CiuB_33r = CuB_33r;

//  and depending on the flavour assumptions rewrite the values (but never rewritting the values model parameters)    
    
    if (FlagFlavU3OfX || FlagUnivOfX) {
        
        if (FlagUnivOfX) {      
//  All equal to uH_33r
            CieH_11r = CuH_33r;
            CieH_22r = CuH_33r;        
            CieH_33r = CuH_33r;
                
            CiuH_11r = CuH_33r;            
            CiuH_22r = CuH_33r; 
        // CiuH_33r = CuH_33r;

            CidH_11r = CuH_33r;
            CidH_22r = CuH_33r;            
            CidH_33r = CuH_33r;
        
        // Currently OfV are only implemented for u quarks so nothing else is needed to apply universality.
        }  
        
//  Proportional to Yukawa interactions Wilson coeff in Warsaw - C=y c - Wilson coeff in model par

        CieH_11r = Yuke * CeH_11r;
        CieH_22r = Yukmu * CeH_22r;        
        CieH_33r = Yuktau * CeH_33r;
                
        CiuH_11r = Yuku * CuH_11r;            
        CiuH_22r = Yukc * CuH_22r;     
        CiuH_33r = Yukt * CuH_33r;
        
        CidH_11r = Yukd * CdH_11r;
        CidH_22r = Yuks * CdH_22r;            
        CidH_33r = Yukb * CdH_33r;
        
        CiuG_11r = Yuku * CuG_11r;  
        CiuG_22r = Yukc * CuG_22r;
        CiuG_33r = Yukt * CuG_33r;

        CiuW_11r = Yuku * CuW_11r;
        CiuW_22r = Yukc * CuW_22r;
        CiuW_33r = Yukt * CuW_33r;

        CiuB_11r = Yuku * CuB_11r;
        CiuB_22r = Yukc * CuB_22r;                   
        CiuB_33r = Yukt * CuB_33r;
    }    
        
//  C2B, C2W, C2WS, C2BS, CDB, CDW, CT are incorporated by change of basis transformation:
//  Write here, before working with the dim 6 interactions,
//  the contributions from O2W and O2B to the other operators.
//  WARNING: Ignoring contributions to 4 fermion-processes for the moment. IMPORTANT FOR LEP2
    
//  WARNING (OBSOLETE MESSAGE?): if some of the parameters below, e.g. CHL1_11, are not floating in the fit this will
//  create a problem since the value generated below CHL1_11 will propagate to the next iteration
//  generating an uncontrolled value of the parameter.
//  (This is so because SetParameters is not called for non-floating parameters.)
//  Possible fix: Not modify model parameters but save everything into internal replicas
//  of each model relevant model par. Those then have to be used in the calculations.
//  Comment out the following lines until this is resolved
    
//  Contributionsfrom C2W, C2B, C2WS, C2BS, CT
    CiHL1_11 = CHL1_11 - (g1_tree*g1_tree/2.0) * (C2B + 0.5 * C2BS);
    CiHL1_22 = CHL1_22 - (g1_tree*g1_tree/2.0) * (C2B + 0.5 * C2BS);
    CiHL1_33 = CHL1_33 - (g1_tree*g1_tree/2.0) * (C2B + 0.5 * C2BS);
    CiHL3_11 = CHL3_11 + (g2_tree*g2_tree/2.0) * (C2W + 0.5 * C2WS);
    CiHL3_22 = CHL3_22 + (g2_tree*g2_tree/2.0) * (C2W + 0.5 * C2WS);
    CiHL3_33 = CHL3_33 + (g2_tree*g2_tree/2.0) * (C2W + 0.5 * C2WS);
    
    CiHQ1_11 = CHQ1_11 + (g1_tree*g1_tree/6.0) * (C2B + 0.5 * C2BS);
    CiHQ1_22 = CHQ1_22 + (g1_tree*g1_tree/6.0) * (C2B + 0.5 * C2BS);
    CiHQ1_33 = CHQ1_33 + (g1_tree*g1_tree/6.0) * (C2B + 0.5 * C2BS);
    CiHQ3_11 = CHQ3_11 + (g2_tree*g2_tree/2.0) * (C2W + 0.5 * C2WS);
    CiHQ3_22 = CHQ3_22 + (g2_tree*g2_tree/2.0) * (C2W + 0.5 * C2WS);
    CiHQ3_33 = CHQ3_33 + (g2_tree*g2_tree/2.0) * (C2W + 0.5 * C2WS);
    
    CiHe_11 = CHe_11 - (g1_tree*g1_tree) * (C2B + 0.5 * C2BS);
    CiHe_22 = CHe_22 - (g1_tree*g1_tree) * (C2B + 0.5 * C2BS);
    CiHe_33 = CHe_33 - (g1_tree*g1_tree) * (C2B + 0.5 * C2BS);
    
    CiHu_11 = CHu_11 + (2.0*g1_tree*g1_tree/3.0) * (C2B + 0.5 * C2BS);
    CiHu_22 = CHu_22 + (2.0*g1_tree*g1_tree/3.0) * (C2B + 0.5 * C2BS);
    CiHu_33 = CHu_33 + (2.0*g1_tree*g1_tree/3.0) * (C2B + 0.5 * C2BS);
    
    CiHd_11 = CHd_11 - (g1_tree*g1_tree/3.0) * (C2B + 0.5 * C2BS);
    CiHd_22 = CHd_22 - (g1_tree*g1_tree/3.0) * (C2B + 0.5 * C2BS);
    CiHd_33 = CHd_33 - (g1_tree*g1_tree/3.0) * (C2B + 0.5 * C2BS);
    
    CiW = CW + g2_tree * C2W;
        
    CiHbox = CHbox - 0.5 * CT + (g1_tree*g1_tree/4.0) * (C2B + 0.5 * C2BS) + (3.0*g2_tree*g2_tree/4.0) * (C2W + 0.5 * C2WS);
    CiHD = CHD - 2.0 * CT + (g1_tree*g1_tree/4.0) * (C2B + 0.5 * C2BS);
    CiH = CH + (2.0*g2_tree*g2_tree*lambdaH_tree) * (C2W + 0.5 * C2WS);
    
//  For the CfH I must use CifH = CifH + ... to account for previous operations.  
    
    CieH_11r = CieH_11r + (g2_tree*g2_tree*Yuke) * (C2W + 0.5 * C2WS);
    CieH_22r = CieH_22r + (g2_tree*g2_tree*Yukmu) * (C2W + 0.5 * C2WS);
    CieH_33r = CieH_33r + (g2_tree*g2_tree*Yuktau) * (C2W + 0.5 * C2WS);
    
    CiuH_11r = CiuH_11r + (g2_tree*g2_tree*Yuku) * (C2W + 0.5 * C2WS);
    CiuH_22r = CiuH_22r + (g2_tree*g2_tree*Yukc) * (C2W + 0.5 * C2WS);
    CiuH_33r = CiuH_33r + (g2_tree*g2_tree*Yukt) * (C2W + 0.5 * C2WS);
    
    CidH_11r = CidH_11r + (g2_tree*g2_tree*Yukd) * (C2W + 0.5 * C2WS);
    CidH_22r = CidH_22r + (g2_tree*g2_tree*Yuks) * (C2W + 0.5 * C2WS);
    CidH_33r = CidH_33r + (g2_tree*g2_tree*Yukb) * (C2W + 0.5 * C2WS);
    
    CiLL_1221 = CLL_1221 + (g2_tree*g2_tree/2.0) * (C2W + 0.5 * C2WS);
    CiLL_2112 = CiLL_1221;
    
//  Contributionsfrom CDW, DB    
    CiHB = CHB + (g1_tree/4.0) * CDB;    
    CiHW = CHW + (g2_tree/4.0) * CDW; 
//    CiHWHB_gaga = CHWHB_gaga; 
//    CiHWHB_gagaorth = CHWHB_gagaorth; 
    CiDHB = CDHB + CDB; 
    CiDHW = CDHW + CDW;
    CiHWB = CHWB + (1.0/4.0) * ( g1_tree * CDW + g2_tree * CDB );
    
//  3) Post-update operations working directly with the dimension six operators  
    
    delta_ZZ = (cW2_tree * CiHW + sW2_tree * CiHB + sW_tree * cW_tree * CiHWB) * v2_over_LambdaNP2;
    delta_AA = (sW2_tree * CiHW + cW2_tree * CiHB - sW_tree * cW_tree * CiHWB) * v2_over_LambdaNP2;
    delta_AZ = 2.0 * sW_tree * cW_tree * (CiHW - CiHB) * v2_over_LambdaNP2
            - (cW2_tree - sW2_tree) * CiHWB * v2_over_LambdaNP2;
    delta_h = (-CiHD / 4.0 + CiHbox) * v2_over_LambdaNP2;
    
//  Calculation of some quantities repeteadly used in the code
    
//  NP corrections to Total Higgs width
    dGammaHTotR1 = deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        dGammaHTotR2 = deltaGammaTotalRatio2();
    } else {
        dGammaHTotR2 = 0.0;  
    }
    
//  Total: to be used in BR functions to check positivity
    GammaHTotR = 1.0 + dGammaHTotR1 + dGammaHTotR2;
    
    // The total theory error in the H width: set to 0.0 for the moment
    eHwidth = deltaGammaTotalRatio1() - deltaGammaTotalRatio1noError();
      
//  Dimension-6 coefficients used in the STXS parameterization
    aiG = 16.0 * M_PI * M_PI * CHG * Mw_tree() * Mw_tree() / g3_tree / g3_tree / LambdaNP2;
    ai3G = CG * Mw_tree() * Mw_tree() / g3_tree / g3_tree / g3_tree / LambdaNP2;
    ai2G =0.0; // Add
    aiT = 2.0 * CiHD * v2_over_LambdaNP2;
    aiH = - 2.0 * CiHbox * v2_over_LambdaNP2;
    aiWW = 0.0; // Add
    aiB = 0.0; // Add
    aiHW = CiDHW * Mw_tree() * Mw_tree() / 2.0 / g2_tree / LambdaNP2;
    aiHB = CiDHB * Mw_tree() * Mw_tree() / 2.0 / g1_tree / LambdaNP2;
    aiA = CiHB * Mw_tree() * Mw_tree() / g1_tree / g1_tree / LambdaNP2;
    aiHQ = CiHQ1_11 * v2_over_LambdaNP2; // Valid only for flavour universal NP
    aipHQ = CiHQ3_11 * v2_over_LambdaNP2; // Valid only for flavour universal NP
    aiHL = CiHL1_11 * v2_over_LambdaNP2; // Valid only for flavour universal NP
    aipHL = CiHL3_11 * v2_over_LambdaNP2; // Valid only for flavour universal NP. From HEL Lagrangian. Not in original note
    aiHu = CiHu_11 * v2_over_LambdaNP2; // Valid only for flavour universal NP
    aiHd = CiHd_11 * v2_over_LambdaNP2; // Valid only for flavour universal NP
    aiHe = CiHe_11 * v2_over_LambdaNP2; // Valid only for flavour universal NP
    aiu = - CiuH_33r * v2_over_LambdaNP2 / Yukt;
    aiuG = CiuG_33r * Mw_tree() * Mw_tree() / g3_tree / LambdaNP2 / Yukt / 4.0; // From HEL.fr Lagrangian. Not in original note. Valid only for flavour universal NP
    
   
//  Dim 6 SMEFT matching
    
    NPSMEFTd6M.getObj().updateNPSMEFTd6Parameters();

    return (true);
}

void NPSMEFTd6::setParameter(const std::string name, const double& value)
{
    if (name.compare("CG") == 0)
        CG = value;
    else if (name.compare("CW") == 0)
        CW = value;
    else if (name.compare("C2B") == 0)
        C2B = value;
    else if (name.compare("C2W") == 0)
        C2W = value;
    else if (name.compare("C2BS") == 0)
        C2BS = value;
    else if (name.compare("C2WS") == 0)
        C2WS = value;
    else if (name.compare("CHG") == 0)
        CHG = value;
    else if (name.compare("CHW") == 0)
        CHW = value;
    else if (name.compare("CHB") == 0)
        CHB = value;
    else if (name.compare("CHWHB_gaga") == 0)
        CHWHB_gaga = value;
    else if (name.compare("CHWHB_gagaorth") == 0)
        CHWHB_gagaorth = value;
    else if (name.compare("CDHB") == 0)
        CDHB = value;
    else if (name.compare("CDHW") == 0)
        CDHW = value;
    else if (name.compare("CDB") == 0)
        CDB = value;
    else if (name.compare("CDW") == 0)
        CDW = value;
    else if (name.compare("CHWB") == 0)
        CHWB = value;
    else if (name.compare("CHD") == 0)
        CHD = value;
    else if (name.compare("CT") == 0)
        CT = value;
    else if (name.compare("CHbox") == 0)
        CHbox = value;
    else if (name.compare("CH") == 0)
        CH = value;
    else if (name.compare("CHL1_11") == 0)
        CHL1_11 = value;
    else if (name.compare("CHL1_12r") == 0)
        CHL1_12r = value;
    else if (name.compare("CHL1_13r") == 0)
        CHL1_13r = value;
    else if (name.compare("CHL1_22") == 0)
        CHL1_22 = value;
    else if (name.compare("CHL1_23r") == 0)
        CHL1_23r = value;
    else if (name.compare("CHL1_33") == 0)
        CHL1_33 = value;
    else if (name.compare("CHL1_12i") == 0)
        CHL1_12i = value;
    else if (name.compare("CHL1_13i") == 0)
        CHL1_13i = value;
    else if (name.compare("CHL1_23i") == 0)
        CHL1_23i = value;
    else if (name.compare("CHL1") == 0) {
        CHL1_11 = value;
        CHL1_12r = 0.0;
        CHL1_13r = 0.0;
        CHL1_22 = value;
        CHL1_23r = 0.0;
        CHL1_33 = value;
        CHL1_12i = 0.0;
        CHL1_13i = 0.0;
        CHL1_23i = 0.0;
    } else if (name.compare("CHL3_11") == 0)
        CHL3_11 = value;
    else if (name.compare("CHL3_12r") == 0)
        CHL3_12r = value;
    else if (name.compare("CHL3_13r") == 0)
        CHL3_13r = value;
    else if (name.compare("CHL3_22") == 0)
        CHL3_22 = value;
    else if (name.compare("CHL3_23r") == 0)
        CHL3_23r = value;
    else if (name.compare("CHL3_33") == 0)
        CHL3_33 = value;
    else if (name.compare("CHL3_12i") == 0)
        CHL3_12i = value;
    else if (name.compare("CHL3_13i") == 0)
        CHL3_13i = value;
    else if (name.compare("CHL3_23i") == 0)
        CHL3_23i = value;
    else if (name.compare("CHL3") == 0) {
        CHL3_11 = value;
        CHL3_12r = 0.0;
        CHL3_13r = 0.0;
        CHL3_22 = value;
        CHL3_23r = 0.0;
        CHL3_33 = value;
        CHL3_12i = 0.0;
        CHL3_13i = 0.0;
        CHL3_23i = 0.0;
    } else if (name.compare("CHe_11") == 0)
        CHe_11 = value;
    else if (name.compare("CHe_12r") == 0)
        CHe_12r = value;
    else if (name.compare("CHe_13r") == 0)
        CHe_13r = value;
    else if (name.compare("CHe_22") == 0)
        CHe_22 = value;
    else if (name.compare("CHe_23r") == 0)
        CHe_23r = value;
    else if (name.compare("CHe_33") == 0)
        CHe_33 = value;
    else if (name.compare("CHe_12i") == 0)
        CHe_12i = value;
    else if (name.compare("CHe_13i") == 0)
        CHe_13i = value;
    else if (name.compare("CHe_23i") == 0)
        CHe_23i = value;
    else if (name.compare("CHe") == 0) {
        CHe_11 = value;
        CHe_12r = 0.0;
        CHe_13r = 0.0;
        CHe_22 = value;
        CHe_23r = 0.0;
        CHe_33 = value;
        CHe_12i = 0.0;
        CHe_13i = 0.0;
        CHe_23i = 0.0;
    } else if (name.compare("CHQ1_11") == 0) {
        CHQ1_11 = value;
        if (FlagPartialQFU){
            CHQ1_22 = value;
        }
    } else if (name.compare("CHQ1_12r") == 0)
        CHQ1_12r = value;
    else if (name.compare("CHQ1_13r") == 0)
        CHQ1_13r = value;
    else if (name.compare("CHQ1_22") == 0) {
        if (!FlagPartialQFU){
            CHQ1_22 = value;
        }
    } else if (name.compare("CHQ1_23r") == 0)
        CHQ1_23r = value;
    else if (name.compare("CHQ1_33") == 0)
        CHQ1_33 = value;
    else if (name.compare("CHQ1_12i") == 0)
        CHQ1_12i = value;
    else if (name.compare("CHQ1_13i") == 0)
        CHQ1_13i = value;
    else if (name.compare("CHQ1_23i") == 0)
        CHQ1_23i = value;
    else if (name.compare("CHQ1") == 0) {
        CHQ1_11 = value;
        CHQ1_12r = 0.0;
        CHQ1_13r = 0.0;
        CHQ1_22 = value;
        CHQ1_23r = 0.0;
        CHQ1_33 = value;
        CHQ1_12i = 0.0;
        CHQ1_13i = 0.0;
        CHQ1_23i = 0.0;
    } else if (name.compare("CHQ3_11") == 0){
        CHQ3_11 = value;
        if (FlagPartialQFU){
            CHQ3_22 = value;
        }
    } else if (name.compare("CHQ3_12r") == 0)
        CHQ3_12r = value;
    else if (name.compare("CHQ3_13r") == 0)
        CHQ3_13r = value;
    else if (name.compare("CHQ3_22") == 0){
        if (!FlagPartialQFU){
            CHQ3_22 = value;
        }
    } else if (name.compare("CHQ3_23r") == 0)
        CHQ3_23r = value;
    else if (name.compare("CHQ3_33") == 0)
        CHQ3_33 = value;
    else if (name.compare("CHQ3_12i") == 0)
        CHQ3_12i = value;
    else if (name.compare("CHQ3_13i") == 0)
        CHQ3_13i = value;
    else if (name.compare("CHQ3_23i") == 0)
        CHQ3_23i = value;
    else if (name.compare("CHQ3") == 0) {
        CHQ3_11 = value;
        CHQ3_12r = 0.0;
        CHQ3_13r = 0.0;
        CHQ3_22 = value;
        CHQ3_23r = 0.0;
        CHQ3_33 = value;
        CHQ3_12i = 0.0;
        CHQ3_13i = 0.0;
        CHQ3_23i = 0.0;
    } else if (name.compare("CHu_11") == 0){
        CHu_11 = value;
        if (FlagPartialQFU){
            CHu_22 = value;
        }
    } else if (name.compare("CHu_12r") == 0)
        CHu_12r = value;
    else if (name.compare("CHu_13r") == 0)
        CHu_13r = value;
    else if (name.compare("CHu_22") == 0){
        if (!FlagPartialQFU){
            CHu_22 = value;
        }
    } else if (name.compare("CHu_23r") == 0)
        CHu_23r = value;
    else if (name.compare("CHu_33") == 0)
        CHu_33 = value;
    else if (name.compare("CHu_12i") == 0)
        CHu_12i = value;
    else if (name.compare("CHu_13i") == 0)
        CHu_13i = value;
    else if (name.compare("CHu_23i") == 0)
        CHu_23i = value;
    else if (name.compare("CHu") == 0) {
        CHu_11 = value;
        CHu_12r = 0.0;
        CHu_13r = 0.0;
        CHu_22 = value;
        CHu_23r = 0.0;
        CHu_33 = value;
        CHu_12i = 0.0;
        CHu_13i = 0.0;
        CHu_23i = 0.0;
    } else if (name.compare("CHd_11") == 0){
        CHd_11 = value;
        if (FlagPartialQFU){
            CHd_22 = value;
        }
    } else if (name.compare("CHd_12r") == 0)
        CHd_12r = value;
    else if (name.compare("CHd_13r") == 0)
        CHd_13r = value;
    else if (name.compare("CHd_22") == 0){
        if (!FlagPartialQFU){
            CHd_22 = value;
        }
    } else if (name.compare("CHd_23r") == 0)
        CHd_23r = value;
    else if (name.compare("CHd_33") == 0)
        CHd_33 = value;
    else if (name.compare("CHd_12i") == 0)
        CHd_12i = value;
    else if (name.compare("CHd_13i") == 0)
        CHd_13i = value;
    else if (name.compare("CHd_23i") == 0)
        CHd_23i = value;
    else if (name.compare("CHd") == 0) {
        CHd_11 = value;
        CHd_12r = 0.0;
        CHd_13r = 0.0;
        CHd_22 = value;
        CHd_23r = 0.0;
        CHd_33 = value;
        CHd_12i = 0.0;
        CHd_13i = 0.0;
        CHd_23i = 0.0;
    } else if (name.compare("CHud_11r") == 0){
        CHud_11r = value;
        if (FlagPartialQFU){
            CHud_22r = value;
        }
    } else if (name.compare("CHud_12r") == 0)
        CHud_12r = value;
    else if (name.compare("CHud_13r") == 0)
        CHud_13r = value;
    else if (name.compare("CHud_22r") == 0){
        if (!FlagPartialQFU){
            CHud_22r = value;
        }
    } else if (name.compare("CHud_23r") == 0)
        CHud_23r = value;
    else if (name.compare("CHud_33r") == 0)
        CHud_33r = value;
    else if (name.compare("CHud_r") == 0) {
        CHud_11r = value;
        CHud_12r = 0.0;
        CHud_13r = 0.0;
        CHud_22r = value;
        CHud_23r = 0.0;
        CHud_33r = value;
    } else if (name.compare("CHud_11i") == 0){
        CHud_11i = value;
        if (FlagPartialQFU){
            CHud_22i = value;
        }
    } else if (name.compare("CHud_12i") == 0)
        CHud_12i = value;
    else if (name.compare("CHud_13i") == 0)
        CHud_13i = value;
    else if (name.compare("CHud_22i") == 0){
        if (!FlagPartialQFU){
            CHud_22i = value;
        }
    } else if (name.compare("CHud_23i") == 0)
        CHud_23i = value;
    else if (name.compare("CHud_33i") == 0)
        CHud_33i = value;
    else if (name.compare("CHud_i") == 0) {
        CHud_11i = value;
        CHud_12i = 0.0;
        CHud_13i = 0.0;
        CHud_22i = value;
        CHud_23i = 0.0;
        CHud_33i = value;
    } else if (name.compare("CeH_11r") == 0){
        if (!FlagFlavU3OfX){
            CeH_11r = value;
        }
    } else if (name.compare("CeH_12r") == 0)
        CeH_12r = value;
    else if (name.compare("CeH_13r") == 0)
        CeH_13r = value;
    else if (name.compare("CeH_22r") == 0){
        if (!FlagFlavU3OfX){
            CeH_22r = value;
        }
    } else if (name.compare("CeH_23r") == 0)
        CeH_23r = value;
    else if (name.compare("CeH_33r") == 0){
        CeH_33r = value;
        if (FlagFlavU3OfX){
            CeH_11r = value;
            CeH_22r = value;
        }
    } else if (name.compare("CeH_11i") == 0)
        CeH_11i = value;
    else if (name.compare("CeH_12i") == 0)
        CeH_12i = value;
    else if (name.compare("CeH_13i") == 0)
        CeH_13i = value;
    else if (name.compare("CeH_22i") == 0)
        CeH_22i = value;
    else if (name.compare("CeH_23i") == 0)
        CeH_23i = value;
    else if (name.compare("CeH_33i") == 0)
        CeH_33i = value;
    else if (name.compare("CuH_11r") == 0){
        if (!FlagFlavU3OfX){
            CuH_11r = value;
        }
    } else if (name.compare("CuH_12r") == 0)
        CuH_12r = value;
    else if (name.compare("CuH_13r") == 0)
        CuH_13r = value;
    else if (name.compare("CuH_22r") == 0){
        if (!FlagFlavU3OfX){
            CuH_22r = value;
        }
    } else if (name.compare("CuH_23r") == 0)
        CuH_23r = value;
    else if (name.compare("CuH_33r") == 0){
        CuH_33r = value;
        if (FlagFlavU3OfX){
            CuH_11r = value;            
            CuH_22r = value;            
        }
    } else if (name.compare("CuH_11i") == 0)
        CuH_11i = value;
    else if (name.compare("CuH_12i") == 0)
        CuH_12i = value;
    else if (name.compare("CuH_13i") == 0)
        CuH_13i = value;
    else if (name.compare("CuH_22i") == 0)
        CuH_22i = value;
    else if (name.compare("CuH_23i") == 0)
        CuH_23i = value;
    else if (name.compare("CuH_33i") == 0)
        CuH_33i = value;
    else if (name.compare("CdH_11r") == 0){
        if (!FlagFlavU3OfX){
            CdH_11r = value;
        }
    } else if (name.compare("CdH_12r") == 0)
        CdH_12r = value;
    else if (name.compare("CdH_13r") == 0)
        CdH_13r = value;
    else if (name.compare("CdH_22r") == 0){
        if (!FlagFlavU3OfX){
            CdH_22r = value;
        }
    } else if (name.compare("CdH_23r") == 0)
        CdH_23r = value;
    else if (name.compare("CdH_33r") == 0){
        CdH_33r = value;
        if (FlagFlavU3OfX){
            CdH_11r = value;
            CdH_22r = value;            
        }
    } else if (name.compare("CdH_11i") == 0)
        CdH_11i = value;
    else if (name.compare("CdH_12i") == 0)
        CdH_12i = value;
    else if (name.compare("CdH_13i") == 0)
        CdH_13i = value;
    else if (name.compare("CdH_22i") == 0)
        CdH_22i = value;
    else if (name.compare("CdH_23i") == 0)
        CdH_23i = value;
    else if (name.compare("CdH_33i") == 0)
        CdH_33i = value;
    else if (name.compare("CuG_11r") == 0){
        if (!FlagFlavU3OfX){
            CuG_11r = value;
        }
    } else if (name.compare("CuG_12r") == 0)
        CuG_12r = value;
    else if (name.compare("CuG_13r") == 0)
        CuG_13r = value;
    else if (name.compare("CuG_22r") == 0){
        if (!FlagFlavU3OfX){
            CuG_22r = value;
        }
    } else if (name.compare("CuG_23r") == 0)
        CuG_23r = value;
    else if (name.compare("CuG_33r") == 0){
        CuG_33r = value;
        if (FlagFlavU3OfX){
            CuG_11r = value;  
            CuG_22r = value;
        }
    } else if (name.compare("CuG_r") == 0) {
        CuG_11r = value;
        CuG_12r = 0.0;
        CuG_13r = 0.0;
        CuG_22r = value;
        CuG_23r = 0.0;
        CuG_33r = value;
    } else if (name.compare("CuG_11i") == 0)
        CuG_11i = value;
    else if (name.compare("CuG_12i") == 0)
        CuG_12i = value;
    else if (name.compare("CuG_13i") == 0)
        CuG_13i = value;
    else if (name.compare("CuG_22i") == 0)
        CuG_22i = value;
    else if (name.compare("CuG_23i") == 0)
        CuG_23i = value;
    else if (name.compare("CuG_33i") == 0)
        CuG_33i = value;
    else if (name.compare("CuG_i") == 0) {
        CuG_11i = value;
        CuG_12i = 0.0;
        CuG_13i = 0.0;
        CuG_22i = value;
        CuG_23i = 0.0;
        CuG_33i = value;
    } else if (name.compare("CuW_11r") == 0){
        if (!FlagFlavU3OfX){
            CuW_11r = value;
        }
    } else if (name.compare("CuW_12r") == 0)
        CuW_12r = value;
    else if (name.compare("CuW_13r") == 0)
        CuW_13r = value;
    else if (name.compare("CuW_22r") == 0){
        if (!FlagFlavU3OfX){
            CuW_22r = value;
        }
    } else if (name.compare("CuW_23r") == 0)
        CuW_23r = value;
    else if (name.compare("CuW_33r") == 0){
        CuW_33r = value;
        if (FlagFlavU3OfX){
            CuW_11r = value;
            CuW_22r = value;
        }
    } else if (name.compare("CuW_r") == 0) {
        CuW_11r = value;
        CuW_12r = 0.0;
        CuW_13r = 0.0;
        CuW_22r = value;
        CuW_23r = 0.0;
        CuW_33r = value;
    } else if (name.compare("CuW_11i") == 0)
        CuW_11i = value;
    else if (name.compare("CuW_12i") == 0)
        CuW_12i = value;
    else if (name.compare("CuW_13i") == 0)
        CuW_13i = value;
    else if (name.compare("CuW_22i") == 0)
        CuW_22i = value;
    else if (name.compare("CuW_23i") == 0)
        CuW_23i = value;
    else if (name.compare("CuW_33i") == 0)
        CuW_33i = value;
    else if (name.compare("CuW_i") == 0) {
        CuW_11i = value;
        CuW_12i = 0.0;
        CuW_13i = 0.0;
        CuW_22i = value;
        CuW_23i = 0.0;
        CuW_33i = value;
    } else if (name.compare("CuB_11r") == 0){
        if (!FlagFlavU3OfX){
            CuB_11r = value;
        }
    } else if (name.compare("CuB_12r") == 0)
        CuB_12r = value;
    else if (name.compare("CuB_13r") == 0)
        CuB_13r = value;
    else if (name.compare("CuB_22r") == 0){
        if (!FlagFlavU3OfX){
            CuB_22r = value;
        }
    } else if (name.compare("CuB_23r") == 0)
        CuB_23r = value;
    else if (name.compare("CuB_33r") == 0){
        CuB_33r = value;
        if (FlagFlavU3OfX){
            CuB_11r = value;
            CuB_22r = value;            
        }
    } else if (name.compare("CuB_r") == 0) {
        CuB_11r = value;
        CuB_12r = 0.0;
        CuB_13r = 0.0;
        CuB_22r = value;
        CuB_23r = 0.0;
        CuB_33r = value;
    } else if (name.compare("CuB_11i") == 0)
        CuB_11i = value;
    else if (name.compare("CuB_12i") == 0)
        CuB_12i = value;
    else if (name.compare("CuB_13i") == 0)
        CuB_13i = value;
    else if (name.compare("CuB_22i") == 0)
        CuB_22i = value;
    else if (name.compare("CuB_23i") == 0)
        CuB_23i = value;
    else if (name.compare("CuB_33i") == 0)
        CuB_33i = value;
    else if (name.compare("CuB_i") == 0) {
        CuB_11i = value;
        CuB_12i = 0.0;
        CuB_13i = 0.0;
        CuB_22i = value;
        CuB_23i = 0.0;
        CuB_33i = value;
//  Several redundancies for the 4-fermionn operators below
    } else if (name.compare("CLL_1111") == 0) {
        CLL_1111 = value;
    } else if (name.compare("CLL_1122") == 0) {
        CLL_1122 = value;
        CLL_2211 = value;
    } else if (name.compare("CLL_1133") == 0) {
        CLL_1133 = value;
        CLL_3311 = value;
    } else if (name.compare("CLL_1221") == 0) {
        CLL_1221 = value;
        CLL_2112 = value;
    } else if (name.compare("CLL_1331") == 0) {
        CLL_1331 = value;
        CLL_3113 = value;
    } else if (name.compare("CLL") == 0) {
        CLL_1111 = value;
        CLL_1221 = value;
        CLL_2112 = value;
        CLL_2211 = value;
        CLL_1122 = value;
        CLL_3311 = value;
        CLL_1133 = value;
        CLL_1331 = value;
        CLL_3113 = value;
    } else if (name.compare("CLQ1_1111") == 0) {
        CLQ1_1111 = value;
    } else if (name.compare("CLQ1_1122") == 0) {
        CLQ1_1122 = value;
    } else if (name.compare("CLQ1_2211") == 0) {
        CLQ1_2211 = value;
    } else if (name.compare("CLQ1_2112") == 0) {
        CLQ1_2112 = value;
    } else if (name.compare("CLQ1_1221") == 0) {
        CLQ1_1221 = value;
    } else if (name.compare("CLQ1_1133") == 0) {
        CLQ1_1133 = value;
    } else if (name.compare("CLQ1_3311") == 0) {
        CLQ1_3311 = value;
    } else if (name.compare("CLQ1_3113") == 0) {
        CLQ1_3113 = value;
    } else if (name.compare("CLQ1_1331") == 0) {
        CLQ1_1331 = value;
    } else if (name.compare("CLQ1_1123") == 0) {
        CLQ1_1123 = value;
    } else if (name.compare("CLQ1_2223") == 0) {
        CLQ1_2223 = value;
    } else if (name.compare("CLQ1_3323") == 0) {
        CLQ1_3323 = value;
    } else if (name.compare("CLQ1_1132") == 0) {
        CLQ1_1132 = value;
    } else if (name.compare("CLQ1_2232") == 0) {
        CLQ1_2232 = value;
    } else if (name.compare("CLQ1_3332") == 0) {
        CLQ1_3332 = value;
    } else if (name.compare("CLQ1") == 0) {
        CLQ1_1111 = value;
        CLQ1_1122 = value;
        CLQ1_2211 = value;
        CLQ1_1221 = value;
        CLQ1_2112 = value;
        CLQ1_1133 = value;
        CLQ1_3311 = value;
        CLQ1_1331 = value;
        CLQ1_3113 = value;
    } else if (name.compare("CLQ3_1111") == 0) {
        CLQ3_1111 = value;
    } else if (name.compare("CLQ3_1122") == 0) {
        CLQ3_1122 = value;
    } else if (name.compare("CLQ3_2211") == 0) {
        CLQ3_2211 = value;
    } else if (name.compare("CLQ3_2112") == 0) {
        CLQ3_2112 = value;
    } else if (name.compare("CLQ3_1221") == 0) {
        CLQ3_1221 = value;
    } else if (name.compare("CLQ3_1133") == 0) {
        CLQ3_1133 = value;
    } else if (name.compare("CLQ3_3311") == 0) {
        CLQ3_3311 = value;
    } else if (name.compare("CLQ3_3113") == 0) {
        CLQ3_3113 = value;
    } else if (name.compare("CLQ3_1331") == 0) {
        CLQ3_1331 = value;
    } else if (name.compare("CLQ3_1123") == 0) {
        CLQ3_1123 = value;
    } else if (name.compare("CLQ3_2223") == 0) {
        CLQ3_2223 = value;
    } else if (name.compare("CLQ3_3323") == 0) {
        CLQ3_3323 = value;
    } else if (name.compare("CLQ3_1132") == 0) {
        CLQ3_1132 = value;
    } else if (name.compare("CLQ3_2232") == 0) {
        CLQ3_2232 = value;
    } else if (name.compare("CLQ3_3332") == 0) {
        CLQ3_3332 = value;
    } else if (name.compare("CLQ3") == 0) {
        CLQ3_1111 = value;
        CLQ3_1122 = value;
        CLQ3_2211 = value;
        CLQ3_1221 = value;
        CLQ3_2112 = value;
        CLQ3_1133 = value;
        CLQ3_3311 = value;
        CLQ3_1331 = value;
        CLQ3_3113 = value;
    } else if (name.compare("Cee") == 0) {
        Cee_1111 = value;
        Cee_1122 = value;
        Cee_2211 = value;
        Cee_1133 = value;
        Cee_3311 = value;
    } else if (name.compare("Cee_1111") == 0) {
        Cee_1111 = value;
    } else if (name.compare("Cee_1122") == 0) {
        Cee_1122 = value;
        Cee_2211 = value;
    } else if (name.compare("Cee_1133") == 0) {
        Cee_1133 = value;
        Cee_3311 = value;
    } else if (name.compare("Ceu") == 0) {
        Ceu_1111 = value;
        Ceu_1122 = value;
        Ceu_2211 = value;
        Ceu_1133 = value;
        Ceu_2233 = value;
        Ceu_3311 = value;
    } else if (name.compare("Ceu_1111") == 0) {
        Ceu_1111 = value;
    } else if (name.compare("Ceu_1122") == 0) {
        Ceu_1122 = value;
    } else if (name.compare("Ceu_2211") == 0) {
        Ceu_2211 = value;
    } else if (name.compare("Ceu_1133") == 0) {
        Ceu_1133 = value;
    } else if (name.compare("Ceu_2233") == 0) {
        Ceu_2233 = value;
    } else if (name.compare("Ceu_3311") == 0) {
        Ceu_3311 = value;
    } else if (name.compare("Ced") == 0) {
        Ced_1111 = value;
        Ced_1122 = value;
        Ced_2211 = value;
        Ced_1133 = value;
        Ced_3311 = value;
    } else if (name.compare("Ced_1111") == 0) {
        Ced_1111 = value;
    } else if (name.compare("Ced_1122") == 0) {
        Ced_1122 = value;
    } else if (name.compare("Ced_2211") == 0) {
        Ced_2211 = value;
    } else if (name.compare("Ced_1133") == 0) {
        Ced_1133 = value;
    } else if (name.compare("Ced_3311") == 0) {
        Ced_3311 = value;
    } else if (name.compare("Ced_1123") == 0) {
        Ced_1123 = value;
    } else if (name.compare("Ced_2223") == 0) {
        Ced_2223 = value;
    } else if (name.compare("Ced_3323") == 0) {
        Ced_3323 = value;
    } else if (name.compare("Ced_1132") == 0) {
        Ced_1132 = value;
    } else if (name.compare("Ced_2232") == 0) {
        Ced_2232 = value;
    } else if (name.compare("Ced_3332") == 0) {
        Ced_3332 = value;
    } else if (name.compare("CLe") == 0) {
        CLe_1111 = value;
        CLe_1122 = value;
        CLe_2211 = value;
        CLe_1133 = value;
        CLe_3311 = value;
    } else if (name.compare("CLe_1111") == 0) {
        CLe_1111 = value;
    } else if (name.compare("CLe_1122") == 0) {
        CLe_1122 = value;
    } else if (name.compare("CLe_2211") == 0) {
        CLe_2211 = value;
    } else if (name.compare("CLe_1133") == 0) {
        CLe_1133 = value;
    } else if (name.compare("CLe_3311") == 0) {
        CLe_3311 = value;
    } else if (name.compare("CLu") == 0) {
        CLu_1111 = value;
        CLu_1122 = value;
        CLu_2211 = value;
        CLu_1133 = value;
        CLu_2233 = value;
        CLu_3311 = value;
    } else if (name.compare("CLu_1111") == 0) {
        CLu_1111 = value;
    } else if (name.compare("CLu_1122") == 0) {
        CLu_1122 = value;
    } else if (name.compare("CLu_2211") == 0) {
        CLu_2211 = value;
    } else if (name.compare("CLu_1133") == 0) {
        CLu_1133 = value;
    } else if (name.compare("CLu_2233") == 0) {
        CLu_2233 = value;
    } else if (name.compare("CLu_3311") == 0) {
        CLu_3311 = value;
    } else if (name.compare("CLd") == 0) {
        CLd_1111 = value;
        CLd_1122 = value;
        CLd_2211 = value;
        CLd_1133 = value;
        CLd_3311 = value;
    } else if (name.compare("CLd_1111") == 0) {
        CLd_1111 = value;
    } else if (name.compare("CLd_1122") == 0) {
        CLd_1122 = value;
    } else if (name.compare("CLd_2211") == 0) {
        CLd_2211 = value;
    } else if (name.compare("CLd_1133") == 0) {
        CLd_1133 = value;
    } else if (name.compare("CLd_3311") == 0) {
        CLd_3311 = value;
    } else if (name.compare("CLd_1123") == 0) {
        CLd_1123 = value;
    } else if (name.compare("CLd_2223") == 0) {
        CLd_2223 = value;
    } else if (name.compare("CLd_3323") == 0) {
        CLd_3323 = value;
    } else if (name.compare("CLd_1132") == 0) {
        CLd_1132 = value;
    } else if (name.compare("CLd_2232") == 0) {
        CLd_2232 = value;
    } else if (name.compare("CLd_3332") == 0) {
        CLd_3332 = value;
    } else if (name.compare("CQe") == 0) {
        CQe_1111 = value;
        CQe_1122 = value;
        CQe_2211 = value;
        CQe_1133 = value;
        CQe_3311 = value;
    } else if (name.compare("CQe_1111") == 0) {
        CQe_1111 = value;
    } else if (name.compare("CQe_1122") == 0) {
        CQe_1122 = value;
    } else if (name.compare("CQe_2211") == 0) {
        CQe_2211 = value;
    } else if (name.compare("CQe_1133") == 0) {
        CQe_1133 = value;
    } else if (name.compare("CQe_3311") == 0) {
        CQe_3311 = value;
    } else if (name.compare("CQe_2311") == 0) {
        CQe_2311 = value;
    } else if (name.compare("CQe_2322") == 0) {
        CQe_2322 = value;
    } else if (name.compare("CQe_2333") == 0) {
        CQe_2333 = value;
    } else if (name.compare("CQe_3211") == 0) {
        CQe_3211 = value;
    } else if (name.compare("CQe_3222") == 0) {
        CQe_3222 = value;
    } else if (name.compare("CLedQ_11") == 0) {
        CLedQ_11 = value;
    } else if (name.compare("CLedQ_22") == 0) {
        CLedQ_22 = value; 
    } else if (name.compare("CpLedQ_11") == 0) {
        CpLedQ_11 = value;
    } else if (name.compare("CpLedQ_22") == 0) {
        CpLedQ_22 = value; 
    } else if (name.compare("CQe_3233") == 0) {
        CQe_3233 = value;
    } else if (name.compare("Lambda_NP") == 0) {
        Lambda_NP = value;        
    } else if (name.compare("BrHinv") == 0) {
//  Always positive
        BrHinv = fabs(value);       
    } else if (name.compare("BrHexo") == 0) {
//  Always positive
        BrHexo = fabs(value);      
    } else if (name.compare("dg1Z") == 0) {
        dg1Z = value;
    } else if (name.compare("dKappaga") == 0) {
        dKappaga = value;
    } else if (name.compare("lambZ") == 0) {
        lambZ = value;        
    } else if (name.compare("eggFint") == 0) {
        eggFint = value;
    } else if (name.compare("eggFpar") == 0) {
        eggFpar = value;
    } else if (name.compare("ettHint") == 0) {
        ettHint = value;
    } else if (name.compare("ettHpar") == 0) {
        ettHpar = value;
    } else if (name.compare("eVBFint") == 0) {
        eVBFint = value;
    } else if (name.compare("eVBFpar") == 0) {
        eVBFpar = value;
    } else if (name.compare("eWHint") == 0) {
        eWHint = value;
    } else if (name.compare("eWHpar") == 0) {
        eWHpar = value;
    } else if (name.compare("eZHint") == 0) {
        eZHint = value;
    } else if (name.compare("eZHpar") == 0) {
        eZHpar = value;
    } else if (name.compare("eeeWBFint") == 0) {
        eeeWBFint = value;
    } else if (name.compare("eeeWBFpar") == 0) {
        eeeWBFpar = value;
    } else if (name.compare("eeeZHint") == 0) {
        eeeZHint = value;
    } else if (name.compare("eeeZHpar") == 0) {
        eeeZHpar = value;
    } else if (name.compare("eeettHint") == 0) {
        eeettHint = value;
    } else if (name.compare("eeettHpar") == 0) {
        eeettHpar = value;
    } else if (name.compare("eepWBFint") == 0) {
        eepWBFint = value;
    } else if (name.compare("eepWBFpar") == 0) {
        eepWBFpar = value;
    } else if (name.compare("eepZBFint") == 0) {
        eepZBFint = value;
    } else if (name.compare("eepZBFpar") == 0) {
        eepZBFpar = value;
    } else if (name.compare("eHggint") == 0) {
        eHggint = value;
    } else if (name.compare("eHggpar") == 0) {
        eHggpar = value;
    } else if (name.compare("eHWWint") == 0) {
        eHWWint = value;
    } else if (name.compare("eHWWpar") == 0) {
        eHWWpar = value;
    } else if (name.compare("eHZZint") == 0) {
        eHZZint = value;
    } else if (name.compare("eHZZpar") == 0) {
        eHZZpar = value;
    } else if (name.compare("eHZgaint") == 0) {
        eHZgaint = value;
    } else if (name.compare("eHZgapar") == 0) {
        eHZgapar = value;
    } else if (name.compare("eHgagaint") == 0) {
        eHgagaint = value;
    } else if (name.compare("eHgagapar") == 0) {
        eHgagapar = value;
    } else if (name.compare("eHmumuint") == 0) {
        eHmumuint = value;
    } else if (name.compare("eHmumupar") == 0) {
        eHmumupar = value;
    } else if (name.compare("eHtautauint") == 0) {
        eHtautauint = value;
    } else if (name.compare("eHtautaupar") == 0) {
        eHtautaupar = value;
    } else if (name.compare("eHccint") == 0) {
        eHccint = value;
    } else if (name.compare("eHccpar") == 0) {
        eHccpar = value;
    } else if (name.compare("eHbbint") == 0) {
        eHbbint = value;
    } else if (name.compare("eHbbpar") == 0) {
        eHbbpar = value;
    } else if (name.compare("eggFHgaga") == 0) {
        eggFHgaga = value;
    } else if (name.compare("eggFHZga") == 0) {
        eggFHZga = value;
    } else if (name.compare("eggFHZZ") == 0) {
        eggFHZZ = value;
    } else if (name.compare("eggFHWW") == 0) {
        eggFHWW = value;
    } else if (name.compare("eggFHtautau") == 0) {
        eggFHtautau = value;
    } else if (name.compare("eggFHbb") == 0) {
        eggFHbb = value;
    } else if (name.compare("eggFHmumu") == 0) {
        eggFHmumu = value;     
    } else if (name.compare("eVBFHgaga") == 0) {
        eVBFHgaga = value;
    } else if (name.compare("eVBFHZga") == 0) {
        eVBFHZga = value;
    } else if (name.compare("eVBFHZZ") == 0) {
        eVBFHZZ = value;
    } else if (name.compare("eVBFHWW") == 0) {
        eVBFHWW = value;
    } else if (name.compare("eVBFHtautau") == 0) {
        eVBFHtautau = value;
    } else if (name.compare("eVBFHbb") == 0) {
        eVBFHbb = value;
    } else if (name.compare("eVBFHmumu") == 0) {
        eVBFHmumu = value;  
    } else if (name.compare("eWHgaga") == 0) {
        eWHgaga = value;
    } else if (name.compare("eWHZga") == 0) {
        eWHZga = value;
    } else if (name.compare("eWHZZ") == 0) {
        eWHZZ = value;
    } else if (name.compare("eWHWW") == 0) {
        eWHWW = value;
    } else if (name.compare("eWHtautau") == 0) {
        eWHtautau = value;
    } else if (name.compare("eWHbb") == 0) {
        eWHbb = value;
    } else if (name.compare("eWHmumu") == 0) {
        eWHmumu = value; 
    } else if (name.compare("eZHgaga") == 0) {
        eZHgaga = value;
    } else if (name.compare("eZHZga") == 0) {
        eZHZga = value;
    } else if (name.compare("eZHZZ") == 0) {
        eZHZZ = value;
    } else if (name.compare("eZHWW") == 0) {
        eZHWW = value;
    } else if (name.compare("eZHtautau") == 0) {
        eZHtautau = value;
    } else if (name.compare("eZHbb") == 0) {
        eZHbb = value;
    } else if (name.compare("eZHmumu") == 0) {
        eZHmumu = value; 
    } else if (name.compare("ettHgaga") == 0) {
        ettHgaga = value;
    } else if (name.compare("ettHZga") == 0) {
        ettHZga = value;
    } else if (name.compare("ettHZZ") == 0) {
        ettHZZ = value;
    } else if (name.compare("ettHWW") == 0) {
        ettHWW = value;
    } else if (name.compare("ettHtautau") == 0) {
        ettHtautau = value;
    } else if (name.compare("ettHbb") == 0) {
        ettHbb = value;
    } else if (name.compare("ettHmumu") == 0) {
        ettHmumu = value;
    } else if (name.compare("eVBFHinv") == 0) {
        eVBFHinv = value;
    } else if (name.compare("eVHinv") == 0) {
        eVHinv = value;
    } else if (name.compare("eVBF_2_Hbox") == 0) {
         eVBF_2_Hbox = value;
    } else if (name.compare("eVBF_2_HQ1_11") == 0) {
         eVBF_2_HQ1_11 = value;
    } else if (name.compare("eVBF_2_Hu_11") == 0) {
         eVBF_2_Hu_11 = value;
    } else if (name.compare("eVBF_2_Hd_11") == 0) {
         eVBF_2_Hd_11 = value;
    } else if (name.compare("eVBF_2_HQ3_11") == 0) {
         eVBF_2_HQ3_11 = value;
    } else if (name.compare("eVBF_2_HD") == 0) {
         eVBF_2_HD = value;
    } else if (name.compare("eVBF_2_HB") == 0) {
         eVBF_2_HB = value;
    } else if (name.compare("eVBF_2_HW") == 0) {
         eVBF_2_HW = value;
    } else if (name.compare("eVBF_2_HWB") == 0) {
         eVBF_2_HWB = value;
    } else if (name.compare("eVBF_2_HG") == 0) {
         eVBF_2_HG = value;
    } else if (name.compare("eVBF_2_DHB") == 0) {
         eVBF_2_DHB = value;
    } else if (name.compare("eVBF_2_DHW") == 0) {
         eVBF_2_DHW = value;
    } else if (name.compare("eVBF_2_DeltaGF") == 0) {
         eVBF_2_DeltaGF = value;
    } else if (name.compare("eVBF_78_Hbox") == 0) {
         eVBF_78_Hbox = value;
    } else if (name.compare("eVBF_78_HQ1_11") == 0) {
         eVBF_78_HQ1_11 = value;
    } else if (name.compare("eVBF_78_Hu_11") == 0) {
         eVBF_78_Hu_11 = value;
    } else if (name.compare("eVBF_78_Hd_11") == 0) {
         eVBF_78_Hd_11 = value;
    } else if (name.compare("eVBF_78_HQ3_11") == 0) {
         eVBF_78_HQ3_11 = value;
    } else if (name.compare("eVBF_78_HD") == 0) {
         eVBF_78_HD = value;
    } else if (name.compare("eVBF_78_HB") == 0) {
         eVBF_78_HB = value;
    } else if (name.compare("eVBF_78_HW") == 0) {
         eVBF_78_HW = value;
    } else if (name.compare("eVBF_78_HWB") == 0) {
         eVBF_78_HWB = value;
    } else if (name.compare("eVBF_78_HG") == 0) {
         eVBF_78_HG = value;
    } else if (name.compare("eVBF_78_DHB") == 0) {
         eVBF_78_DHB = value;
    } else if (name.compare("eVBF_78_DHW") == 0) {
         eVBF_78_DHW = value;
    } else if (name.compare("eVBF_78_DeltaGF") == 0) {
         eVBF_78_DeltaGF = value;
    } else if (name.compare("eVBF_1314_Hbox") == 0) {
         eVBF_1314_Hbox = value;
    } else if (name.compare("eVBF_1314_HQ1_11") == 0) {
         eVBF_1314_HQ1_11 = value;
    } else if (name.compare("eVBF_1314_Hu_11") == 0) {
         eVBF_1314_Hu_11 = value;
    } else if (name.compare("eVBF_1314_Hd_11") == 0) {
         eVBF_1314_Hd_11 = value;
    } else if (name.compare("eVBF_1314_HQ3_11") == 0) {
         eVBF_1314_HQ3_11 = value;
    } else if (name.compare("eVBF_1314_HD") == 0) {
         eVBF_1314_HD = value;
    } else if (name.compare("eVBF_1314_HB") == 0) {
         eVBF_1314_HB = value;
    } else if (name.compare("eVBF_1314_HW") == 0) {
         eVBF_1314_HW = value;
    } else if (name.compare("eVBF_1314_HWB") == 0) {
         eVBF_1314_HWB = value;
    } else if (name.compare("eVBF_1314_HG") == 0) {
         eVBF_1314_HG = value;
    } else if (name.compare("eVBF_1314_DHB") == 0) {
         eVBF_1314_DHB = value;
    } else if (name.compare("eVBF_1314_DHW") == 0) {
         eVBF_1314_DHW = value;
    } else if (name.compare("eVBF_1314_DeltaGF") == 0) {
         eVBF_1314_DeltaGF = value;
    } else if (name.compare("eWH_2_Hbox") == 0) {
         eWH_2_Hbox = value;
    } else if (name.compare("eWH_2_HQ3_11") == 0) {
         eWH_2_HQ3_11 = value;
    } else if (name.compare("eWH_2_HD") == 0) {
         eWH_2_HD = value;
    } else if (name.compare("eWH_2_HW") == 0) {
         eWH_2_HW = value;
    } else if (name.compare("eWH_2_HWB") == 0) {
         eWH_2_HWB = value;
    } else if (name.compare("eWH_2_DHW") == 0) {
         eWH_2_DHW = value;
    } else if (name.compare("eWH_2_DeltaGF") == 0) {
         eWH_2_DeltaGF = value;
    } else if (name.compare("eWH_78_Hbox") == 0) {
         eWH_78_Hbox = value;
    } else if (name.compare("eWH_78_HQ3_11") == 0) {
         eWH_78_HQ3_11 = value;
    } else if (name.compare("eWH_78_HD") == 0) {
         eWH_78_HD = value;
    } else if (name.compare("eWH_78_HW") == 0) {
         eWH_78_HW = value;
    } else if (name.compare("eWH_78_HWB") == 0) {
         eWH_78_HWB = value;
    } else if (name.compare("eWH_78_DHW") == 0) {
         eWH_78_DHW = value;
    } else if (name.compare("eWH_78_DeltaGF") == 0) {
         eWH_78_DeltaGF = value;
    } else if (name.compare("eWH_1314_Hbox") == 0) {
         eWH_1314_Hbox = value;
    } else if (name.compare("eWH_1314_HQ3_11") == 0) {
         eWH_1314_HQ3_11 = value;
    } else if (name.compare("eWH_1314_HD") == 0) {
         eWH_1314_HD = value;
    } else if (name.compare("eWH_1314_HW") == 0) {
         eWH_1314_HW = value;
    } else if (name.compare("eWH_1314_HWB") == 0) {
         eWH_1314_HWB = value;
    } else if (name.compare("eWH_1314_DHW") == 0) {
         eWH_1314_DHW = value;
    } else if (name.compare("eWH_1314_DeltaGF") == 0) {
         eWH_1314_DeltaGF = value;
    } else if (name.compare("eZH_2_Hbox") == 0) {
         eZH_2_Hbox = value;
    } else if (name.compare("eZH_2_HQ1_11") == 0) {
         eZH_2_HQ1_11 = value;
    } else if (name.compare("eZH_2_Hu_11") == 0) {
         eZH_2_Hu_11 = value;
    } else if (name.compare("eZH_2_Hd_11") == 0) {
         eZH_2_Hd_11 = value;
    } else if (name.compare("eZH_2_HQ3_11") == 0) {
         eZH_2_HQ3_11 = value;
    } else if (name.compare("eZH_2_HD") == 0) {
         eZH_2_HD = value;
    } else if (name.compare("eZH_2_HB") == 0) {
         eZH_2_HB = value;
    } else if (name.compare("eZH_2_HW") == 0) {
         eZH_2_HW = value;
    } else if (name.compare("eZH_2_HWB") == 0) {
         eZH_2_HWB = value;
    } else if (name.compare("eZH_2_DHB") == 0) {
         eZH_2_DHB = value;
    } else if (name.compare("eZH_2_DHW") == 0) {
         eZH_2_DHW = value;
    } else if (name.compare("eZH_2_DeltaGF") == 0) {
         eZH_2_DeltaGF = value;
    } else if (name.compare("eZH_78_Hbox") == 0) {
         eZH_78_Hbox = value;
    } else if (name.compare("eZH_78_HQ1_11") == 0) {
         eZH_78_HQ1_11 = value;
    } else if (name.compare("eZH_78_Hu_11") == 0) {
         eZH_78_Hu_11 = value;
    } else if (name.compare("eZH_78_Hd_11") == 0) {
         eZH_78_Hd_11 = value;
    } else if (name.compare("eZH_78_HQ3_11") == 0) {
         eZH_78_HQ3_11 = value;
    } else if (name.compare("eZH_78_HD") == 0) {
         eZH_78_HD = value;
    } else if (name.compare("eZH_78_HB") == 0) {
         eZH_78_HB = value;
    } else if (name.compare("eZH_78_HW") == 0) {
         eZH_78_HW = value;
    } else if (name.compare("eZH_78_HWB") == 0) {
         eZH_78_HWB = value;
    } else if (name.compare("eZH_78_DHB") == 0) {
         eZH_78_DHB = value;
    } else if (name.compare("eZH_78_DHW") == 0) {
         eZH_78_DHW = value;
    } else if (name.compare("eZH_78_DeltaGF") == 0) {
         eZH_78_DeltaGF = value;
    } else if (name.compare("eZH_1314_Hbox") == 0) {
         eZH_1314_Hbox = value;
    } else if (name.compare("eZH_1314_HQ1_11") == 0) {
         eZH_1314_HQ1_11 = value;
    } else if (name.compare("eZH_1314_Hu_11") == 0) {
         eZH_1314_Hu_11 = value;
    } else if (name.compare("eZH_1314_Hd_11") == 0) {
         eZH_1314_Hd_11 = value;
    } else if (name.compare("eZH_1314_HQ3_11") == 0) {
         eZH_1314_HQ3_11 = value;
    } else if (name.compare("eZH_1314_HD") == 0) {
         eZH_1314_HD = value;
    } else if (name.compare("eZH_1314_HB") == 0) {
         eZH_1314_HB = value;
    } else if (name.compare("eZH_1314_HW") == 0) {
         eZH_1314_HW = value;
    } else if (name.compare("eZH_1314_HWB") == 0) {
         eZH_1314_HWB = value;
    } else if (name.compare("eZH_1314_DHB") == 0) {
         eZH_1314_DHB = value;
    } else if (name.compare("eZH_1314_DHW") == 0) {
         eZH_1314_DHW = value;
    } else if (name.compare("eZH_1314_DeltaGF") == 0) {
         eZH_1314_DeltaGF = value;
    } else if (name.compare("ettH_2_HG") == 0) {
         ettH_2_HG = value;
    } else if (name.compare("ettH_2_G") == 0) {
         ettH_2_G = value;
    } else if (name.compare("ettH_2_uG_33r") == 0) {
         ettH_2_uG_33r = value;
    } else if (name.compare("ettH_2_DeltagHt") == 0) {
         ettH_2_DeltagHt = value;
    } else if (name.compare("ettH_78_HG") == 0) {
         ettH_78_HG = value;
    } else if (name.compare("ettH_78_G") == 0) {
         ettH_78_G = value;
    } else if (name.compare("ettH_78_uG_33r") == 0) {
         ettH_78_uG_33r = value;
    } else if (name.compare("ettH_78_DeltagHt") == 0) {
         ettH_78_DeltagHt = value;
    } else if (name.compare("ettH_1314_HG") == 0) {
         ettH_1314_HG = value;
    } else if (name.compare("ettH_1314_G") == 0) {
         ettH_1314_G = value;
    } else if (name.compare("ettH_1314_uG_33r") == 0) {
         ettH_1314_uG_33r = value;
    } else if (name.compare("ettH_1314_DeltagHt") == 0) {
         ettH_1314_DeltagHt = value;
    } else
        NPbase::setParameter(name, value);
}

bool NPSMEFTd6::CheckParameters(const std::map<std::string, double>& DPars)
{
    if (FlagLeptonUniversal && FlagQuarkUniversal) {
        if (FlagRotateCHWCHB) {
              for (int i = 0; i < NNPSMEFTd6Vars_LFU_QFU; i++) {
                    if (DPars.find(NPSMEFTd6VarsRot_LFU_QFU[i]) == DPars.end()) {
                          std::cout << "ERROR: Missing mandatory NPSMEFTd6_LFU_QFU parameter "
                          << NPSMEFTd6VarsRot_LFU_QFU[i] << std::endl;
                          raiseMissingModelParameterCount();
                          addMissingModelParameter(NPSMEFTd6VarsRot_LFU_QFU[i]);
                    }
              }
        } else {
              for (int i = 0; i < NNPSMEFTd6Vars_LFU_QFU; i++) {
                    if (DPars.find(NPSMEFTd6Vars_LFU_QFU[i]) == DPars.end()) {
                          std::cout << "ERROR: Missing mandatory NPSMEFTd6_LFU_QFU parameter "
                          << NPSMEFTd6Vars_LFU_QFU[i] << std::endl;
                          raiseMissingModelParameterCount();
                          addMissingModelParameter(NPSMEFTd6Vars_LFU_QFU[i]);
                    }
              }
        }
    } else if (!FlagLeptonUniversal && !FlagQuarkUniversal) {
        if (FlagRotateCHWCHB) {
              for (int i = 0; i < NNPSMEFTd6Vars; i++) {
                    if (DPars.find(NPSMEFTd6VarsRot[i]) == DPars.end()) {
                          std::cout << "ERROR: Missing mandatory NPSMEFTd6 parameter "
                          << NPSMEFTd6VarsRot[i] << std::endl;
                          raiseMissingModelParameterCount();
                          addMissingModelParameter(NPSMEFTd6VarsRot[i]);
                    }
              }
        } else {
              for (int i = 0; i < NNPSMEFTd6Vars; i++) {
                    if (DPars.find(NPSMEFTd6Vars[i]) == DPars.end()) {
                          std::cout << "ERROR: Missing mandatory NPSMEFTd6 parameter "
                          << NPSMEFTd6Vars[i] << std::endl;
                          raiseMissingModelParameterCount();
                          addMissingModelParameter(NPSMEFTd6Vars[i]);
                    }
              }
        }
          
    } else
        throw std::runtime_error("Error in NPSMEFTd6::CheckParameters()");

    return (NPbase::CheckParameters(DPars));
}

bool NPSMEFTd6::setFlag(const std::string name, const bool value)
{
    bool res = false;
    if (name.compare("QuadraticTerms") == 0) {
        FlagQuadraticTerms = value;
        if(value) setModelLinearized(false);
        res = true;
    } else if (name.compare("RotateCHWCHB") == 0) {
        FlagRotateCHWCHB = value;
        res = true;
    } else if (name.compare("PartialQFU") == 0) {        
        FlagPartialQFU = value;
        res = true; 
    } else if (name.compare("FlavU3OfX") == 0) {        
        FlagFlavU3OfX = value;
        res = true; 
    } else if (name.compare("UnivOfX") == 0) {        
        FlagUnivOfX = value;
        res = true; 
    } else if (name.compare("HiggsSM") == 0) {        
        FlagHiggsSM = value;
        res = true;
    } else if (name.compare("LoopHd6") == 0) {        
        FlagLoopHd6 = value;
        res = true;
    } else if (name.compare("LoopH3d6Quad") == 0) {        
        FlagLoopH3d6Quad = value;
        res = true;
    } else
        res = NPbase::setFlag(name, value);

    return (res);
}


////////////////////////////////////////////////////////////////////////

double NPSMEFTd6::CHF1_diag(const Particle F) const
{
    if (F.is("NEUTRINO_1") || F.is("ELECTRON"))
        return CiHL1_11;
    else if (F.is("NEUTRINO_2") || F.is("MU"))
        return CiHL1_22;
    else if (F.is("NEUTRINO_3") || F.is("TAU"))
        return CiHL1_33;
    else if (F.is("UP") || F.is("DOWN"))
        return CiHQ1_11;
    else if (F.is("CHARM") || F.is("STRANGE"))
        return CiHQ1_22;
    else if (F.is("TOP") || F.is("BOTTOM"))
        return CiHQ1_33;
    else
        throw std::runtime_error("NPSMEFTd6::CHF1_diag(): wrong argument");
}

double NPSMEFTd6::CHF3_diag(const Particle F) const
{
    if (F.is("NEUTRINO_1") || F.is("ELECTRON"))
        return CiHL3_11;
    else if (F.is("NEUTRINO_2") || F.is("MU"))
        return CiHL3_22;
    else if (F.is("NEUTRINO_3") || F.is("TAU"))
        return CiHL3_33;
    else if (F.is("UP") || F.is("DOWN"))
        return CiHQ3_11;
    else if (F.is("CHARM") || F.is("STRANGE"))
        return CiHQ3_22;
    else if (F.is("TOP") || F.is("BOTTOM"))
        return CiHQ3_33;
    else
        throw std::runtime_error("NPSMEFTd6::CHF3_diag(): wrong argument");
}

double NPSMEFTd6::CHf_diag(const Particle f) const
{
    if (f.is("NEUTRINO_1") || f.is("NEUTRINO_2") || f.is("NEUTRINO_3"))
        return 0.0;
    else if (f.is("ELECTRON"))
        return CiHe_11;
    else if (f.is("MU"))
        return CiHe_22;
    else if (f.is("TAU"))
        return CiHe_33;
    else if (f.is("UP"))
        return CiHu_11;
    else if (f.is("CHARM"))
        return CiHu_22;
    else if (f.is("TOP"))
        return CiHu_33;
    else if (f.is("DOWN"))
        return CiHd_11;
    else if (f.is("STRANGE"))
        return CiHd_22;
    else if (f.is("BOTTOM"))
        return CiHd_33;
    else
        throw std::runtime_error("NPSMEFTd6::CHf_diag(): wrong argument");
}

gslpp::complex NPSMEFTd6::CHud_diag(const Particle u) const
{
    if (!u.is("QUARK") || u.getIndex() % 2 != 0)
        throw std::runtime_error("NPSMEFTd6::CHud_diag(): wrong argument");

    if (u.is("UP"))
        return gslpp::complex(CHud_11r, CHud_11i, false);
    else if (u.is("CHARM"))
        return gslpp::complex(CHud_22r, CHud_22i, false);
    else if (u.is("TOP"))
        return gslpp::complex(CHud_22r, CHud_33i, false);
    else
        throw std::runtime_error("NPSMEFTd6::CHud_diag(): wrong argument");
}

gslpp::complex NPSMEFTd6::CfH_diag(const Particle f) const
{
    if (f.is("NEUTRINO_1") || f.is("NEUTRINO_2") || f.is("NEUTRINO_3"))
        return 0.0;
    else if (f.is("ELECTRON"))
        return gslpp::complex(CieH_11r, CeH_11i, false);
    else if (f.is("MU"))
        return gslpp::complex(CieH_22r, CeH_22i, false);
    else if (f.is("TAU"))
        return gslpp::complex(CieH_33r, CeH_33i, false);
    else if (f.is("UP"))
        return gslpp::complex(CiuH_11r, CuH_11i, false);
    else if (f.is("CHARM"))
        return gslpp::complex(CiuH_22r, CuH_22i, false);
    else if (f.is("TOP"))
        return gslpp::complex(CiuH_33r, CuH_33i, false);
    else if (f.is("DOWN"))
        return gslpp::complex(CidH_11r, CdH_11i, false);
    else if (f.is("STRANGE"))
        return gslpp::complex(CidH_22r, CdH_22i, false);
    else if (f.is("BOTTOM"))
        return gslpp::complex(CidH_33r, CdH_33i, false);
    else
        throw std::runtime_error("NPSMEFTd6::CfH_diag(): wrong argument");
}
      
gslpp::complex NPSMEFTd6::CfG_diag(const Particle f) const
{
    if (f.is("NEUTRINO_1") || f.is("NEUTRINO_2") || f.is("NEUTRINO_3"))
        return 0.0;
    else if (f.is("ELECTRON"))
        return 0.0;
    else if (f.is("MU"))
        return 0.0;
    else if (f.is("TAU"))
        return 0.0;
    else if (f.is("UP"))
        return gslpp::complex(CiuG_11r, CuG_11i, false);
    else if (f.is("CHARM"))
        return gslpp::complex(CiuG_22r, CuG_22i, false);
    else if (f.is("TOP"))
        return gslpp::complex(CiuG_33r, CuG_33i, false);
    else if (f.is("DOWN"))
        return 0.0;
    else if (f.is("STRANGE"))
        return 0.0;
    else if (f.is("BOTTOM"))
        return 0.0;
    else
        throw std::runtime_error("NPSMEFTd6::CfG_diag(): wrong argument");
}
      
gslpp::complex NPSMEFTd6::CfW_diag(const Particle f) const
{
    if (f.is("NEUTRINO_1") || f.is("NEUTRINO_2") || f.is("NEUTRINO_3"))
        return 0.0;
    else if (f.is("ELECTRON"))
        return 0.0;
    else if (f.is("MU"))
        return 0.0;
    else if (f.is("TAU"))
        return 0.0;
    else if (f.is("UP"))
        return gslpp::complex(CiuW_11r, CuW_11i, false);
    else if (f.is("CHARM"))
        return gslpp::complex(CiuW_22r, CuW_22i, false);
    else if (f.is("TOP"))
        return gslpp::complex(CiuW_33r, CuW_33i, false);
    else if (f.is("DOWN"))
        return 0.0;
    else if (f.is("STRANGE"))
        return 0.0;
    else if (f.is("BOTTOM"))
        return 0.0;
    else
        throw std::runtime_error("NPSMEFTd6::CfW_diag(): wrong argument");
}
      
gslpp::complex NPSMEFTd6::CfB_diag(const Particle f) const
{
    if (f.is("NEUTRINO_1") || f.is("NEUTRINO_2") || f.is("NEUTRINO_3"))
        return 0.0;
    else if (f.is("ELECTRON"))
        return 0.0;
    else if (f.is("MU"))
        return 0.0;
    else if (f.is("TAU"))
        return 0.0;
    else if (f.is("UP"))
        return gslpp::complex(CiuB_11r, CuB_11i, false);
    else if (f.is("CHARM"))
        return gslpp::complex(CiuB_22r, CuB_22i, false);
    else if (f.is("TOP"))
        return gslpp::complex(CiuB_33r, CuB_33i, false);
    else if (f.is("DOWN"))
        return 0.0;
    else if (f.is("STRANGE"))
        return 0.0;
    else if (f.is("BOTTOM"))
        return 0.0;
    else
        throw std::runtime_error("NPSMEFTd6::CfB_diag(): wrong argument");
}


////////////////////////////////////////////////////////////////////////

double NPSMEFTd6::DeltaGF() const
{
    return ((CiHL3_11 + CiHL3_22 - 0.5 * (CiLL_1221 + CiLL_2112)) * v2_over_LambdaNP2);
}

double NPSMEFTd6::obliqueS() const
{
    return (4.0 * sW_tree * cW_tree * CiHWB / alphaMz() * v2_over_LambdaNP2);
}

double NPSMEFTd6::obliqueT() const
{
    return (-CiHD / 2.0 / alphaMz() * v2_over_LambdaNP2);
}

double NPSMEFTd6::obliqueU() const
{
    return 0.0;
}

double NPSMEFTd6::obliqueW() const
{
    return (- g2_tree * g2_tree * (C2W + 0.5 * C2WS) * v2_over_LambdaNP2 / 2.0);
}

double NPSMEFTd6::obliqueY() const
{
    return (- g2_tree * g2_tree * (C2B + 0.5 * C2BS) * v2_over_LambdaNP2 / 2.0);
}

/////////////////////////////// Deviations in the experimental values of the SM input parameters /////////////////////////////////////////

double NPSMEFTd6::deltaMz() const
{
    //  Ref. value from SM EW fit 2018
    return ( (Mz - 91.1882) / 91.1882 );
}
    
double NPSMEFTd6::deltaMz2() const
{
    return ( 0.0 );
}
        
double NPSMEFTd6::deltaMh() const
{
    //  Ref. value from SM EW fit 2018
    return ( (mHl - 125.1) / 125.1 );
}
    
double NPSMEFTd6::deltaMh2() const
{
    return ( 0.0 );
}
    
double NPSMEFTd6::deltamt() const
{
    //  Ref. value from SM EW fit 2018
    return ( (mtpole - 173.2) / 173.2 );
}
    
double NPSMEFTd6::deltamt2() const
{
    return ( 0.0 );
}

double NPSMEFTd6::deltamb() const
{
    //  Ref. value fixed in SM EW fit 2018: from PDG 2018
    return ( ((quarks[BOTTOM].getMass()) - 4.18) / 4.18 );
}
    
double NPSMEFTd6::deltamb2() const
{
    return ( 0.0 );
}
    
double NPSMEFTd6::deltamc() const
{
    //  Ref. value fixed in SM EW fit 2018: from PDG 2018
    return ( ((quarks[CHARM].getMass()) - 1.275) / 1.275 );
}
    
double NPSMEFTd6::deltamc2() const
{
    return ( 0.0 );
}
   
double NPSMEFTd6::deltamtau() const
{
    //  Ref. value fixed in SM EW fit 2018: from PDG 2018
    return ( ((leptons[TAU].getMass()) - 1.77686) / 1.77686 );
}
    
double NPSMEFTd6::deltamtau2() const
{
    return ( 0.0 );
}
 
double NPSMEFTd6::deltaGmu() const
{
    //  Ref. value fixed in SM EW fit 2018: from PDG 2018
    return ( (GF - 1.1663787/100000.0 ) / (1.1663787/100000.0) );
}
    
double NPSMEFTd6::deltaGmu2() const
{
    return ( 0.0 );
}
    
double NPSMEFTd6::deltaaMZ() const
{
    //  Ref. value from SM EW fit 2018
    return ( (aleMz - 0.007754941997887603) / 0.007754941997887603 );
}
    
double NPSMEFTd6::deltaaMZ2() const
{
    return ( 0.0 );
}

double NPSMEFTd6::deltaa0() const
{
    //  Ref. value fixed in SM EW fit 2018: from PDG 2018
    return ( (aleMz - 0.0072973525664) / 0.0072973525664 );
}
    
double NPSMEFTd6::deltaa02() const
{
    return ( 0.0 );
}
    
double NPSMEFTd6::deltaaSMZ() const
{
    //  Ref. value from SM EW fit 2018
    return ( (AlsMz - 0.1180) / 0.1180 );
}
    
double NPSMEFTd6::deltaaSMZ2() const
{
    return ( 0.0 );
}
    
    
////////////////////////////////////////////////////////////////////////

double NPSMEFTd6::Mw() const
{
    return (trueSM.Mw() - Mw_tree() / 4.0 / (cW2_tree - sW2_tree)
            *(4.0 * sW_tree * cW_tree * CiHWB * v2_over_LambdaNP2
            + cW2_tree * CiHD * v2_over_LambdaNP2
            + 2.0 * sW2_tree * DeltaGF()));
}

double NPSMEFTd6::deltaMwd6() const
{
    return (- 1.0 / 4.0 / (cW2_tree - sW2_tree)
            *(4.0 * sW_tree * cW_tree * CiHWB * v2_over_LambdaNP2
            + cW2_tree * CiHD * v2_over_LambdaNP2
            + 2.0 * sW2_tree * DeltaGF()));
}

double NPSMEFTd6::deltaMwd62() const
{
    double dMW = 0.0;
    
    return (dMW*dMW);
}

double NPSMEFTd6::deltaGamma_Wff(const Particle fi, const Particle fj) const
{
    double G0 = GF * pow(Mz*cW_tree, 3.0) / 6.0 / sqrt(2.0) / M_PI;
    double deltaGamma_Wij;
    double GammaW_tree;
    double CHF3ij;
    
    if (fj.getIndex() - fi.getIndex() == 1)
        CHF3ij = CHF3_diag(fi);
    else
        CHF3ij = 0.;
    
    if (fi.is("QUARK")) {        
        GammaW_tree = Nc * G0;
    } else {
        GammaW_tree = G0;        
    }
    
    deltaGamma_Wij = - 3.0 * GammaW_tree / 4.0 / (cW2_tree - sW2_tree)
            *(4.0 * sW_tree * cW_tree * CiHWB * v2_over_LambdaNP2
            + cW2_tree * CiHD * v2_over_LambdaNP2
            + 2.0 * (1.0 + cW2_tree) / 3.0 * DeltaGF());
    
    deltaGamma_Wij = deltaGamma_Wij + 2.0 * GammaW_tree * CHF3ij * v2_over_LambdaNP2;

    return deltaGamma_Wij;    
}


double NPSMEFTd6::GammaW(const Particle fi, const Particle fj) const
{
    return ( trueSM.GammaW(fi, fj) + deltaGamma_Wff(fi, fj) );          
}

double NPSMEFTd6::deltaGamma_W() const
{
    double G0 = GF * pow(Mz*cW_tree, 3.0) / 6.0 / sqrt(2.0) / M_PI;
    double GammaW_tree = (3.0 + 2.0 * Nc) * G0;

    return (- 3.0 * GammaW_tree / 4.0 / (cW2_tree - sW2_tree)
            *(4.0 * sW_tree * cW_tree * CiHWB * v2_over_LambdaNP2
            + cW2_tree * CiHD * v2_over_LambdaNP2
            + 2.0 * (1.0 + cW2_tree) / 3.0 * DeltaGF())
            + 2.0 * G0 * (CiHL3_11 + CiHL3_22 + CiHL3_33 + Nc*(CiHQ3_11 + CiHQ3_22)) * v2_over_LambdaNP2);          
//            + 2.0 * GammaW_tree / 3.0 * (CiHL3_11 + CiHQ3_11 + CiHQ3_22) * v2_over_LambdaNP2);    
}

double NPSMEFTd6::GammaW() const
{
    return ( trueSM.GammaW() + deltaGamma_W() );
}

double NPSMEFTd6::deltaGwd6() const
{
    return ( deltaGamma_W() / trueSM.GammaW() );
}
    
double NPSMEFTd6::deltaGwd62() const
{
    double dWW = 0.0;
    
    return (dWW*dWW);
}
    
double NPSMEFTd6::deltaGzd6() const
{
    return ( deltaGamma_Z() / trueSM.Gamma_Z() );
}
    
double NPSMEFTd6::deltaGzd62() const
{
    double dWZ = 0.0;
    
    return (dWZ*dWZ);
}

double NPSMEFTd6::deltaGV_f(const Particle p) const
{
    return (deltaGL_f(p) + deltaGR_f(p));
}

double NPSMEFTd6::deltaGA_f(const Particle p) const
{
    return (deltaGL_f(p) - deltaGR_f(p));
}

double NPSMEFTd6::deltaGL_f(const Particle p) const
{
    double I3p = p.getIsospin(), Qp = p.getCharge();
    double CHF1 = CHF1_diag(p);
    double CHF3 = CHF3_diag(p);
    double NPindirect;

    NPindirect = -I3p / 4.0 * (CiHD * v2_over_LambdaNP2 + 2.0 * DeltaGF())
                - Qp * sW2_tree / 4.0 / (cW2_tree - sW2_tree)
                *((4.0 * cW_tree / sW_tree * CiHWB + CiHD) * v2_over_LambdaNP2 + 2.0 * DeltaGF());

    double NPdirect = -0.5 * (CHF1 - 2.0 * I3p * CHF3) * v2_over_LambdaNP2;
    return (NPindirect + NPdirect);
}

double NPSMEFTd6::deltaGR_f(const Particle p) const
{
    double Qp = p.getCharge();
    double CHf = CHf_diag(p);
    double NPindirect;

    NPindirect = -Qp * sW2_tree / 4.0 / (cW2_tree - sW2_tree)
                *((4.0 * cW_tree / sW_tree * CiHWB + CiHD) * v2_over_LambdaNP2 + 2.0 * DeltaGF());

    double NPdirect = -0.5 * CHf*v2_over_LambdaNP2;
    return (NPindirect + NPdirect);
}


////////////////////////////////////////////////////////////////////////

gslpp::complex NPSMEFTd6::deltaGL_Wff(const Particle pbar, const Particle p) const
{
    if (pbar.getIndex() + 1 != p.getIndex() || pbar.getIndex() % 2 != 0)
        throw std::runtime_error("NPSMEFTd6::deltaGL_Wff(): Not implemented");

    double CHF3 = CHF3_diag(pbar);
    double NPindirect;

    NPindirect = -cW2_tree / 4.0 / (cW2_tree - sW2_tree)
                * ((4.0 * sW_tree / cW_tree * CiHWB + CiHD) * v2_over_LambdaNP2 + 2.0 * DeltaGF());

    double NPdirect = CHF3 * v2_over_LambdaNP2;
    return (NPindirect + NPdirect);
}

gslpp::complex NPSMEFTd6::deltaGR_Wff(const Particle pbar, const Particle p) const
{
    if (pbar.getIndex() + 1 != p.getIndex() || pbar.getIndex() % 2 != 0)
        throw std::runtime_error("NPSMEFTd6::deltaGR_Wff(): Not implemented");

    gslpp::complex CHud = CHud_diag(pbar);
    return (0.5 * CHud * v2_over_LambdaNP2);
}

double NPSMEFTd6::deltaG_hgg() const
{
    return (CHG * v2_over_LambdaNP2 / v());
}

double NPSMEFTd6::deltaG_hggRatio() const
{    
    double m_t = mtpole;
    double m_b = quarks[BOTTOM].getMass();
    double m_c = quarks[CHARM].getMass();
    double tau_t = 4.0 * m_t * m_t / mHl / mHl;
    double tau_b = 4.0 * m_b * m_b / mHl / mHl;
    double tau_c = 4.0 * m_c * m_c / mHl / mHl;
    double aSPiv = AlsMz / 16.0 / M_PI / v();
    gslpp::complex gSM, dg;
    gslpp::complex dKappa_t = cLHd6 * deltaG_hff(quarks[TOP]) / (-m_t / v());
    gslpp::complex dKappa_b = cLHd6 * deltaG_hff(quarks[BOTTOM]) / (-m_b / v());
    gslpp::complex dKappa_c = cLHd6 * deltaG_hff(quarks[CHARM]) / (-m_c / v());
    double deltaloc = deltaG_hgg();
    
    gSM = aSPiv * (AH_f(tau_t) + AH_f(tau_b) + AH_f(tau_c));
    
    dg = deltaloc/gSM + (aSPiv/gSM) * (dKappa_t*AH_f(tau_t) + dKappa_b*AH_f(tau_b) + dKappa_c*AH_f(tau_c));

    return dg.real();
}

double NPSMEFTd6::deltaG1_hWW() const
{
    return (( 2.0 * CiHW - sqrt( M_PI * aleMz ) * CiDHW / sW_tree ) * v2_over_LambdaNP2 / v());
}

double NPSMEFTd6::deltaG2_hWW() const
{
    return ( - sqrt( M_PI * aleMz ) * ( CiDHW / sW_tree ) * v2_over_LambdaNP2 / v());
}

double NPSMEFTd6::deltaG3_hWW() const
{
    double NPindirect;

    NPindirect = 2.0 * cW2_tree * Mz * Mz / v()
                * (delta_h - 1.0 / 2.0 / (cW2_tree - sW2_tree)
                * ((4.0 * sW_tree * cW_tree * CiHWB + cW2_tree * CiHD) * v2_over_LambdaNP2 + DeltaGF()));

    return NPindirect;
}

double NPSMEFTd6::deltaG1_hZZ() const
{
    return ( (delta_ZZ - 0.5 * sqrt( M_PI * aleMz ) * (CiDHB / cW_tree + CiDHW / sW_tree) * v2_over_LambdaNP2 )/ v());
}

double NPSMEFTd6::deltaG2_hZZ() const
{
    return ( - sqrt( M_PI * aleMz ) * ( CiDHB / cW_tree + CiDHW / sW_tree ) * v2_over_LambdaNP2 / v());
}

double NPSMEFTd6::deltaG3_hZZ() const
{
    double NPindirect = Mz * Mz / v() * (-0.5 * CiHD * v2_over_LambdaNP2 + delta_h - 0.5 * DeltaGF());
    double NPdirect = Mz * Mz / v() * CiHD * v2_over_LambdaNP2;
    return (NPindirect + NPdirect);
}

double NPSMEFTd6::deltaG1_hZA() const
{
    return ( (delta_AZ + 0.5 * sqrt( M_PI * aleMz ) * (CiDHB / sW_tree - CiDHW / cW_tree) * v2_over_LambdaNP2 )/ v());
}

double NPSMEFTd6::deltaG1_hZARatio() const
{    
    double m_t = mtpole;
    double m_b = quarks[BOTTOM].getMass();
    double m_c = quarks[CHARM].getMass();
    double m_tau = leptons[TAU].getMass();
    double m_mu = leptons[MU].getMass();
    
    double M_w_2 = (trueSM.Mw())*(trueSM.Mw());
    
    double Qt = quarks[TOP].getCharge();
    double Qb = quarks[BOTTOM].getCharge();
    double Qc = quarks[CHARM].getCharge();
    double Qtau = leptons[TAU].getCharge();
    double Qmu = leptons[MU].getCharge();

    double tau_t = 4.0 * m_t * m_t / mHl / mHl;
    double tau_b = 4.0 * m_b * m_b / mHl / mHl;
    double tau_c = 4.0 * m_c * m_c / mHl / mHl;
    double tau_tau = 4.0 * m_tau * m_tau / mHl / mHl;
    double tau_mu = 4.0 * m_mu * m_mu / mHl / mHl;
    double tau_W = 4.0 * M_w_2 / mHl / mHl;    
    
    double lambda_t = 4.0 * m_t * m_t / Mz / Mz;
    double lambda_b = 4.0 * m_b * m_b / Mz / Mz;
    double lambda_c = 4.0 * m_c * m_c / Mz / Mz;
    double lambda_tau = 4.0 * m_tau * m_tau / Mz / Mz;
    double lambda_mu = 4.0 * m_mu * m_mu / Mz / Mz;
    double lambda_W = 4.0 * M_w_2 / Mz / Mz; 
    double alpha2 = sqrt(2.0)*GF*M_w_2 / M_PI;
    double aPiv = sqrt(ale*alpha2) / 4.0 / M_PI / v();
    
//  mod. of Higgs couplings
    gslpp::complex gSM, dg;
    gslpp::complex dKappa_t = cLHd6 * deltaG_hff(quarks[TOP]) / (-m_t / v());
    gslpp::complex dKappa_b = cLHd6 * deltaG_hff(quarks[BOTTOM]) / (-m_b / v());
    gslpp::complex dKappa_c = cLHd6 * deltaG_hff(quarks[CHARM]) / (-m_c / v());
    gslpp::complex dKappa_tau = cLHd6 * deltaG_hff(leptons[TAU]) / (-m_tau / v());
    gslpp::complex dKappa_mu = cLHd6 * deltaG_hff(leptons[MU]) / (-m_mu / v());
    double dKappa_W = cLHd6 * (0.5 * v() / M_w_2)*deltaG3_hWW();
    
//  mod of EW vector couplings vf =2 gvf    
    double vSMt = 2.0*(quarks[TOP].getIsospin()) - 4.0 * Qt * sW2_tree;
    double vSMb = 2.0*(quarks[BOTTOM].getIsospin()) - 4.0 * Qb * sW2_tree;
    double vSMc = 2.0*(quarks[CHARM].getIsospin()) - 4.0 * Qc * sW2_tree;
    double vSMtau = 2.0*(leptons[TAU].getIsospin()) - 4.0 * Qtau * sW2_tree;
    double vSMmu = 2.0*(leptons[MU].getIsospin()) - 4.0 * Qmu * sW2_tree;
    
    double dvSMt = cLHd6 * 2.0*deltaGV_f(quarks[TOP]);
    double dvSMb = cLHd6 * 2.0*deltaGV_f(quarks[BOTTOM]);
    double dvSMc = cLHd6 * 2.0*deltaGV_f(quarks[CHARM]);
    double dvSMtau = cLHd6 * 2.0*deltaGV_f(leptons[TAU]);
    double dvSMmu = cLHd6 * 2.0*deltaGV_f(leptons[MU]);
    
    double deltaloc = deltaG1_hZA();
    
    gSM = -aPiv * ((3.0*vSMt*Qt*AHZga_f(tau_t,lambda_t) + 
            3.0*vSMb*Qb*AHZga_f(tau_b,lambda_b) + 
            3.0*vSMc*Qc*AHZga_f(tau_c,lambda_c) +
            vSMtau*Qtau*AHZga_f(tau_tau,lambda_tau)+
            vSMmu*Qmu*AHZga_f(tau_mu,lambda_mu))/cW_tree + 
            AHZga_W(tau_W,lambda_W));
    
    dg = deltaloc/gSM - (aPiv/gSM) * (
            (3.0*vSMt*dKappa_t*Qt*AHZga_f(tau_t,lambda_t) + 
            3.0*vSMb*dKappa_b*Qb*AHZga_f(tau_b,lambda_b) + 
            3.0*vSMc*dKappa_c*Qc*AHZga_f(tau_c,lambda_c)+ 
            dKappa_tau*vSMtau*Qtau*AHZga_f(tau_tau,lambda_tau)+ 
            dKappa_mu*vSMmu*Qmu*AHZga_f(tau_mu,lambda_mu))/cW_tree + 
            dKappa_W*AHZga_W(tau_W,lambda_W) + 
            (3.0*dvSMt*Qt*AHZga_f(tau_t,lambda_t) + 
            3.0*dvSMb*Qb*AHZga_f(tau_b,lambda_b) + 
            3.0*dvSMc*Qc*AHZga_f(tau_c,lambda_c)+ 
            dvSMtau*Qtau*AHZga_f(tau_tau,lambda_tau)+ 
            dvSMmu*Qmu*AHZga_f(tau_mu,lambda_mu))/cW_tree
            );

    return dg.real();
}

double NPSMEFTd6::deltaG2_hZA() const
{
    return ( sqrt( M_PI * aleMz ) * ( CiDHB / sW_tree - CiDHW / cW_tree ) * v2_over_LambdaNP2 / v());
}

double NPSMEFTd6::deltaG_hAA() const
{
    return (delta_AA / v());
}

double NPSMEFTd6::deltaG_hAARatio() const
{    
    double m_t = mtpole;
    double m_b = quarks[BOTTOM].getMass();
    double m_c = quarks[CHARM].getMass();
    double m_tau = leptons[TAU].getMass();
    double m_mu = leptons[MU].getMass();
    
    double M_w_2 = (trueSM.Mw())*(trueSM.Mw());
    
    double Qt = quarks[TOP].getCharge();
    double Qb = quarks[BOTTOM].getCharge();
    double Qc = quarks[CHARM].getCharge();
    double Qtau = leptons[TAU].getCharge();
    double Qmu = leptons[MU].getCharge();

    double tau_t = 4.0 * m_t * m_t / mHl / mHl;
    double tau_b = 4.0 * m_b * m_b / mHl / mHl;
    double tau_c = 4.0 * m_c * m_c / mHl / mHl;
    double tau_tau = 4.0 * m_tau * m_tau / mHl / mHl;
    double tau_mu = 4.0 * m_mu * m_mu / mHl / mHl;
    double tau_W = 4.0 * M_w_2 / mHl / mHl;    

    double aPiv = ale / 8.0 / M_PI / v();
    gslpp::complex gSM, dg;
    gslpp::complex dKappa_t = cLHd6 * deltaG_hff(quarks[TOP]) / (-m_t / v());
    gslpp::complex dKappa_b = cLHd6 * deltaG_hff(quarks[BOTTOM]) / (-m_b / v());
    gslpp::complex dKappa_c = cLHd6 * deltaG_hff(quarks[CHARM]) / (-m_c / v());
    gslpp::complex dKappa_tau = cLHd6 * deltaG_hff(leptons[TAU]) / (-m_tau / v());
    gslpp::complex dKappa_mu = cLHd6 * deltaG_hff(leptons[MU]) / (-m_mu / v());
    double dKappa_W = cLHd6 * (0.5 * v() / M_w_2)*deltaG3_hWW();
    
    double deltaloc = deltaG_hAA();
    
    gSM = aPiv * (3.0*Qt*Qt*AH_f(tau_t) + 
            3.0*Qb*Qb*AH_f(tau_b) + 
            3.0*Qc*Qc*AH_f(tau_c) + 
            Qtau*Qtau*AH_f(tau_tau) + 
            Qmu*Qmu*AH_f(tau_mu) + 
            AH_W(tau_W));
    
    dg = deltaloc/gSM + (aPiv/gSM) * (
            3.0*Qt*Qt*dKappa_t*AH_f(tau_t) + 
            3.0*Qb*Qb*dKappa_b*AH_f(tau_b) + 
            3.0*Qc*Qc*dKappa_c*AH_f(tau_c) +
            dKappa_tau*Qtau*Qtau*AH_f(tau_tau) +
            dKappa_mu*Qmu*Qmu*AH_f(tau_mu) +
            dKappa_W*AH_W(tau_W)
            );

    return dg.real();
}

gslpp::complex NPSMEFTd6::deltaG_hff(const Particle p) const
{
    /* The effects of the RG running are neglected. */
    double mf;
    if (p.is("TOP"))
        //mf = p.getMass(); // m_t(m_t)
        mf = mtpole; // pole mass
    else
        mf = p.getMass();
    gslpp::complex CfH = CfH_diag(p);
    return (-mf / v() * (delta_h - 0.5 * DeltaGF())
            + CfH * v2_over_LambdaNP2 / sqrt(2.0));
}

double NPSMEFTd6::deltaG_hhhRatio() const
{    
    double dg;
    
    dg = -0.5 * DeltaGF() + 3.0 * delta_h - 2.0 * CiH * v2_over_LambdaNP2 * v2/mHl/mHl;

    return dg;
}

gslpp::complex NPSMEFTd6::deltaGL_Wffh(const Particle pbar, const Particle p) const
{
    if (pbar.getIndex() + 1 != p.getIndex() || pbar.getIndex() % 2 != 0)
        throw std::runtime_error("NPSMEFTd6::deltaGL_Wffh(): Not implemented");

    double CHF3 = CHF3_diag(pbar);
    return (2.0 * sqrt(2.0) * Mz * cW_tree / v() / v() * CHF3 * v2_over_LambdaNP2);
}

gslpp::complex NPSMEFTd6::deltaGR_Wffh(const Particle pbar, const Particle p) const
{
    if (pbar.getIndex() + 1 != p.getIndex() || pbar.getIndex() % 2 != 0)
        throw std::runtime_error("NPSMEFTd6::deltaGR_Wffh(): Not implemented");

    gslpp::complex CHud = CHud_diag(pbar);
    return (sqrt(2.0) * Mz * cW_tree / v() / v() * CHud * v2_over_LambdaNP2);
}

double NPSMEFTd6::deltaGL_Zffh(const Particle p) const
{
    double I3p = p.getIsospin();
    double CHF1 = CHF1_diag(p);
    double CHF3 = CHF3_diag(p);
    return (-2.0 * Mz / v() / v() * (CHF1 - 2.0 * I3p * CHF3) * v2_over_LambdaNP2);
}

double NPSMEFTd6::deltaGR_Zffh(const Particle p) const
{
    double CHf = CHf_diag(p);
    return (-2.0 * Mz / v() / v() * CHf * v2_over_LambdaNP2);
}
      
gslpp::complex NPSMEFTd6::deltaG_hGff(const Particle p) const
{
    /* Set to 0. for the moment */

    return 0.;
}
      
gslpp::complex NPSMEFTd6::deltaG_hZff(const Particle p) const
{
    /* Set to 0. for the moment */

    return 0.;
}
      
gslpp::complex NPSMEFTd6::deltaG_hAff(const Particle p) const
{
    /* Set to 0. for the moment */

    return 0.;
}
      
gslpp::complex NPSMEFTd6::deltaG_Gff(const Particle p) const
{
    /* Set to 0. for the moment */

    return 0.;
}
      
gslpp::complex NPSMEFTd6::deltaG_Zff(const Particle p) const
{
    /* Set to 0. for the moment */

    return 0.;
}
      
gslpp::complex NPSMEFTd6::deltaG_Aff(const Particle p) const
{
    /* Set to 0. for the moment */

    return 0.;
}
      
double NPSMEFTd6::deltag3G() const
{
    /* Set to 0. for the moment */

    return 0.;
}


////////////////////////////////////////////////////////////////////////

gslpp::complex NPSMEFTd6::f_triangle(const double tau) const
{
    gslpp::complex tmp;
    if (tau >= 1.0) {
        tmp = asin(1.0 / sqrt(tau));
        return (tmp * tmp);
    } else {
        tmp = log((1.0 + sqrt(1.0 - tau)) / (1.0 - sqrt(1.0 - tau))) - M_PI * gslpp::complex::i();
        return (-0.25 * tmp * tmp);
    }
}

gslpp::complex NPSMEFTd6::g_triangle(const double tau) const
{
    gslpp::complex tmp;
    if (tau >= 1.0) {
        tmp = sqrt(tau -1.0) * asin(1.0 / sqrt(tau));
        return tmp;
    } else {
        tmp = sqrt(1.0 - tau) * ( log((1.0 + sqrt(1.0 - tau)) / (1.0 - sqrt(1.0 - tau))) - M_PI * gslpp::complex::i() );
        return 0.5 * tmp;
    }
}

gslpp::complex NPSMEFTd6::I_triangle_1(const double tau, const double lambda) const
{
    gslpp::complex tmp;
    
    tmp = ( tau*lambda * (f_triangle(tau)- f_triangle(lambda)) + 2.0 * tau * (g_triangle(tau)- g_triangle(lambda)) ) / (tau-lambda);
    
    tmp = tau*lambda * ( 1.0 + tmp ) / (2.0*(tau-lambda));
    
    return tmp;
}

gslpp::complex NPSMEFTd6::I_triangle_2(const double tau, const double lambda) const
{
    gslpp::complex tmp;
    
    tmp = - 0.5 * tau*lambda * (f_triangle(tau)- f_triangle(lambda)) / (tau-lambda);
    
    return tmp;
}

gslpp::complex NPSMEFTd6::AH_f(const double tau) const
{
    return (2.0 * tau * (1.0 + (1.0 - tau) * f_triangle(tau)));
}

gslpp::complex NPSMEFTd6::AH_W(const double tau) const
{
    return -( 2.0 + 3.0 * tau + 3.0 * tau * (2.0 - tau) * f_triangle(tau) );
}

gslpp::complex NPSMEFTd6::AHZga_f(const double tau, const double lambda) const
{
    return I_triangle_1(tau,lambda) - I_triangle_2(tau,lambda);
}

gslpp::complex NPSMEFTd6::AHZga_W(const double tau, const double lambda) const
{
    gslpp::complex tmp;
    
    double tan2w = trueSM.sW2() / trueSM.cW2();
    
    tmp = 4.0 * (3.0 - tan2w ) * I_triangle_2(tau,lambda);
    
    tmp = tmp + ((1.0 +2.0 / tau)* tan2w - (5.0 + 2.0/tau)) * I_triangle_1(tau,lambda);

    return sqrt(trueSM.cW2()) * tmp;
}

double NPSMEFTd6::muggH(const double sqrt_s) const
{
    
    double C1 = 0.0066; //It seems to be independent of energy 
    
    double m_t = mtpole;
    //doulbe m_t = quarks[TOP].getMass();
    double m_b = quarks[BOTTOM].getMass();
    double m_c = quarks[CHARM].getMass();

    /* L_eff_SM = (G_eff_t_SM + G_eff_b_SM)*hGG */
    gslpp::complex G_eff_t_SM = AlsMz / 16.0 / M_PI / v() * AH_f(4.0 * m_t * m_t / mHl / mHl);
    gslpp::complex G_eff_b_SM = AlsMz / 16.0 / M_PI / v() * AH_f(4.0 * m_b * m_b / mHl / mHl);
    gslpp::complex G_eff_c_SM = AlsMz / 16.0 / M_PI / v() * AH_f(4.0 * m_c * m_c / mHl / mHl);
    gslpp::complex G_eff_SM = G_eff_t_SM + G_eff_b_SM + G_eff_c_SM;

    //double sigma_tt_SM = trueSM.computeSigmaggH_tt(sqrt_s);
    //double sigma_bb_SM = trueSM.computeSigmaggH_bb(sqrt_s);
    //double sigma_tb_SM = trueSM.computeSigmaggH_tb(sqrt_s);
    //gslpp::complex tmp = (2.0 * dKappa_t * sigma_tt_SM
    //        + 2.0 * dKappa_b * sigma_bb_SM
    //        + (dKappa_t + dKappa_b) * sigma_tb_SM)
    //        / (sigma_tt_SM + sigma_bb_SM + sigma_tb_SM);
    
    gslpp::complex dKappa_t = cLHd6 * deltaG_hff(quarks[TOP]) / (-m_t / v());
    gslpp::complex dKappa_b = cLHd6 * deltaG_hff(quarks[BOTTOM]) / (-m_b / v());
    gslpp::complex dKappa_c = cLHd6 * deltaG_hff(quarks[CHARM]) / (-m_c / v());

    gslpp::complex tmpHG = CHG / v() * v2_over_LambdaNP2 / G_eff_SM;
    gslpp::complex tmpt = G_eff_t_SM * dKappa_t / G_eff_SM;
    gslpp::complex tmpb = G_eff_b_SM * dKappa_b / G_eff_SM;
    gslpp::complex tmpc = G_eff_c_SM * dKappa_c / G_eff_SM;
    
    double mu = (1.0 + 2.0 * ( tmpt.real() + tmpb.real() + tmpc.real() + tmpHG.real() ) );
    
//  Linear contribution from Higgs self-coupling
    mu = mu + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    mu = mu + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        gslpp::complex tmp2 = tmpt +tmpb +tmpc + tmpHG;
        
        mu += tmp2.abs2();
        
    }
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += eggFint + eggFpar;
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::muggHH(const double sqrt_s) const
{
    double mu = 1.0;
    double A1HH = 0.0, A2HH = 0.0, A3HH = 0.0, A4HH = 0.0, A5HH = 0.0;
    double A6HH = 0.0, A7HH = 0.0, A8HH = 0.0, A9HH = 0.0, A10HH = 0.0;
    double A11HH = 0.0, A12HH = 0.0, A13HH = 0.0, A14HH = 0.0, A15HH = 0.0;    
    double ct,c2t,c3,cg,c2g;
    
    if (sqrt_s == 14.0) {
        
        // From the cut-based analysis. Table IV

        A1HH = 1.70;
        A2HH = 10.7;
        A3HH = 0.117;
        A4HH = 6.11;
        A5HH = 217.0;
        A6HH = -7.56;
        A7HH = -0.819;
        A8HH = 1.95;
        A9HH = 10.90;
        A10HH = 51.6;
        A11HH = -3.86;
        A12HH = -12.5;
        A13HH = 1.46;
        A14HH = 5.49;
        A15HH = 58.4;
        
    } else if (sqrt_s == 100.0) {
        
        // From the cut-based analysis. Table IV

        A1HH = 1.59;
        A2HH = 12.8;
        A3HH = 0.090;
        A4HH = 5.2;
        A5HH = 358.0;
        A6HH = -7.66;
        A7HH = -0.681;
        A8HH = 1.83;
        A9HH = 9.25;
        A10HH = 51.2;
        A11HH = -2.61;
        A12HH = -7.35;
        A13HH = 1.03;
        A14HH = 4.65;
        A15HH = 65.5;
        
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::muggHH()");
    
    ct= 1.0 - 0.5 * DeltaGF() + delta_h - v() * CiuH_33r * v2_over_LambdaNP2 / sqrt(2.0)/ mtpole;
    c2t = delta_h - 3.0 *v() * CiuH_33r * v2_over_LambdaNP2 / 2.0 /sqrt(2.0)/ mtpole;
    c3 = 1.0 + deltaG_hhhRatio();
    cg = M_PI * CHG * v2_over_LambdaNP2 / AlsMz; 
    c2g = cg;

// In the SM the Eq. returns 0.999. Fix that small offset by adding 0.0010    
    mu = 0.0010 + A1HH*ct*ct*ct*ct +
            A2HH*c2t*c2t +
            A3HH*ct*ct*c3*c3 +
            A4HH*cg*cg*c3*c3 +
            A5HH*c2g*c2g +
            A6HH*c2t*ct*ct +
            A7HH*ct*ct*ct*c3 +
            A8HH*c2t*ct*c3 +
            A9HH*c2t*cg*c3 +
            A10HH*c2t*c2g +
            A11HH*ct*ct*cg*c3 +
            A12HH*ct*ct*c2g +
            A13HH*ct*c3*c3*cg +
            A14HH*ct*c3*c2g +
            A15HH*cg*c3*c2g;

    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::muVBF(const double sqrt_s) const
{
    double mu = 1.0;
    
    double C1 = 0.0; 
    
    if (sqrt_s == 1.96) {
        
        C1 = 0.0; // N.A.

        mu += 
                +120936. * (1. + eVBF_2_Hbox ) * CiHbox / LambdaNP2
                -9422.68 * (1. + eVBF_2_HQ1_11 ) * CiHQ1_11 / LambdaNP2
                -10683.8 * (1. + eVBF_2_Hu_11 ) * CiHu_11 / LambdaNP2
                +4055.59 * (1. + eVBF_2_Hd_11 ) * CiHd_11 / LambdaNP2
                -229691. * (1. + eVBF_2_HQ3_11 ) * CiHQ3_11 / LambdaNP2
                -170093. * (1. + eVBF_2_HD ) * CiHD / LambdaNP2
                +8971.22 * (1. + eVBF_2_HB ) * CiHB / LambdaNP2
                -65827.6 * (1. + eVBF_2_HW ) * CiHW / LambdaNP2
                -323514. * (1. + eVBF_2_HWB ) * CiHWB / LambdaNP2
                +481332. * (1. + eVBF_2_HG ) * CHG / LambdaNP2
                +1255.16 * (1. + eVBF_2_DHB ) * CiDHB / LambdaNP2
                -34956.7 * (1. + eVBF_2_DHW ) * CiDHW / LambdaNP2
                -4.511 * (1. + eVBF_2_DeltaGF ) * DeltaGF()
                -3.481 * deltaMwd6()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

        mu +=  0.0;

        }
        
    } else if (sqrt_s == 7.0) {
        
        C1 = 0.0065;

        mu += 
                +121582. * (1. + eVBF_78_Hbox ) * CiHbox / LambdaNP2
                +13546.6 * (1. + eVBF_78_HQ1_11 ) * CiHQ1_11 / LambdaNP2
                -27657.6 * (1. + eVBF_78_Hu_11 ) * CiHu_11 / LambdaNP2
                +8892.12 * (1. + eVBF_78_Hd_11 ) * CiHd_11 / LambdaNP2
                -411400. * (1. + eVBF_78_HQ3_11 ) * CiHQ3_11 / LambdaNP2
                -164286. * (1. + eVBF_78_HD ) * CiHD / LambdaNP2
                -423.123 * (1. + eVBF_78_HB ) * CiHB / LambdaNP2
                -89854. * (1. + eVBF_78_HW ) * CiHW / LambdaNP2
                -312617. * (1. + eVBF_78_HWB ) * CiHWB / LambdaNP2
                -82956.8 * (1. + eVBF_78_HG ) * CHG / LambdaNP2
                -279.08 * (1. + eVBF_78_DHB ) * CiDHB / LambdaNP2
                -54861. * (1. + eVBF_78_DHW ) * CiDHW / LambdaNP2
                -4.479 * (1. + eVBF_78_DeltaGF ) * DeltaGF()
                -3.22 * deltaMwd6()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

        mu +=  0.0;

        }
        
    } else if (sqrt_s == 8.0) {
        
        C1 = 0.0065;

        mu += 
                +121042. * (1. + eVBF_78_Hbox ) * CiHbox / LambdaNP2
                +12739.3 * (1. + eVBF_78_HQ1_11 ) * CiHQ1_11 / LambdaNP2
                -28367.7 * (1. + eVBF_78_Hu_11 ) * CiHu_11 / LambdaNP2
                +9134.21 * (1. + eVBF_78_Hd_11 ) * CiHd_11 / LambdaNP2
                -423704. * (1. + eVBF_78_HQ3_11 ) * CiHQ3_11 / LambdaNP2
                -165182. * (1. + eVBF_78_HD ) * CiHD / LambdaNP2
                -349.242 * (1. + eVBF_78_HB ) * CiHB / LambdaNP2
                -87279.4 * (1. + eVBF_78_HW ) * CiHW / LambdaNP2
                -313449. * (1. + eVBF_78_HWB ) * CiHWB / LambdaNP2
                -69421.9 * (1. + eVBF_78_HG ) * CHG / LambdaNP2
                -373.338 * (1. + eVBF_78_DHB ) * CiDHB / LambdaNP2
                -57028.1 * (1. + eVBF_78_DHW ) * CiDHW / LambdaNP2
                -4.472 * (1. + eVBF_78_DeltaGF ) * DeltaGF()
                -3.138 * deltaMwd6()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

        mu +=  0.0;

        }
    } else if (sqrt_s == 13.0) {
        
        C1 = 0.0064;

        mu += 
                +121798. * (1. + eVBF_1314_Hbox ) * CiHbox / LambdaNP2
                +10339.7 * (1. + eVBF_1314_HQ1_11 ) * CiHQ1_11 / LambdaNP2
                -30827.2 * (1. + eVBF_1314_Hu_11 ) * CiHu_11 / LambdaNP2
                +10564.3 * (1. + eVBF_1314_Hd_11 ) * CiHd_11 / LambdaNP2
                -466270. * (1. + eVBF_1314_HQ3_11 ) * CiHQ3_11 / LambdaNP2
                -164119. * (1. + eVBF_1314_HD ) * CiHD / LambdaNP2
                -61.471 * (1. + eVBF_1314_HB ) * CiHB / LambdaNP2
                -82985.3 * (1. + eVBF_1314_HW ) * CiHW / LambdaNP2
                -313815. * (1. + eVBF_1314_HWB ) * CiHWB / LambdaNP2
                -36554. * (1. + eVBF_1314_HG ) * CHG / LambdaNP2
                -725.694 * (1. + eVBF_1314_DHB ) * CiDHB / LambdaNP2
                -65253.4 * (1. + eVBF_1314_DHW ) * CiDHW / LambdaNP2
                -4.474 * (1. + eVBF_1314_DeltaGF ) * DeltaGF()
                -3.109 * deltaMwd6()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 14.0) {
        
        C1 = 0.0064;

        mu += 
                +120948. * (1. + eVBF_1314_Hbox ) * CiHbox / LambdaNP2
                +9896.36 * (1. + eVBF_1314_HQ1_11 ) * CiHQ1_11 / LambdaNP2
                -31371. * (1. + eVBF_1314_Hu_11 ) * CiHu_11 / LambdaNP2
                +10716.4 * (1. + eVBF_1314_Hd_11 ) * CiHd_11 / LambdaNP2
                -473497. * (1. + eVBF_1314_HQ3_11 ) * CiHQ3_11 / LambdaNP2
                -164672. * (1. + eVBF_1314_HD ) * CiHD / LambdaNP2
                -60.253 * (1. + eVBF_1314_HB ) * CiHB / LambdaNP2
                -83504.9 * (1. + eVBF_1314_HW ) * CiHW / LambdaNP2
                -314059. * (1. + eVBF_1314_HWB ) * CiHWB / LambdaNP2
                -33627.6 * (1. + eVBF_1314_HG ) * CHG / LambdaNP2
                -775.959 * (1. + eVBF_1314_DHB ) * CiDHB / LambdaNP2
                -66336.3 * (1. + eVBF_1314_DHW ) * CiDHW / LambdaNP2
                -4.474 * (1. + eVBF_1314_DeltaGF ) * DeltaGF()
                -3.193 * deltaMwd6()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
   
    } else if (sqrt_s == 27.0) {
        
        C1 = 0.0062; // From arXiv: 1902.00134

        mu += 
                +120777. * CiHbox / LambdaNP2
                +6664.27 * CiHQ1_11 / LambdaNP2
                -34230.7 * CiHu_11 / LambdaNP2
                +12917.3 * CiHd_11 / LambdaNP2
                -536216. * CiHQ3_11 / LambdaNP2
                -163493. * CiHD / LambdaNP2
                +58.33 * CiHB / LambdaNP2
                -81360.5 * CiHW / LambdaNP2
                -313026. * CiHWB / LambdaNP2
                -16430. * CHG / LambdaNP2
                -1314.45 * CiDHB / LambdaNP2
                -75884.6 * CiDHW / LambdaNP2
                -4.475 * DeltaGF()
                -2.99 * deltaMwd6()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
        
    } else if (sqrt_s == 100.0) {
        
        C1 = 0.0; // N.A.

        mu += 
                +121714. * CiHbox / LambdaNP2
                -2261.73 * CiHQ1_11 / LambdaNP2
                -42045.4 * CiHu_11 / LambdaNP2
                +17539.2 * CiHd_11 / LambdaNP2
                -674206. * CiHQ3_11 / LambdaNP2
                -163344. * CiHD / LambdaNP2
                +71.488 * CiHB / LambdaNP2
                -90808.2 * CiHW / LambdaNP2
                -312544. * CiHWB / LambdaNP2
                -8165.65 * CHG / LambdaNP2
                -2615.48 * CiDHB / LambdaNP2
                -96539.6 * CiDHW / LambdaNP2
                -4.452 * DeltaGF()
                -2.949 * deltaMwd6()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
        
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::muVBF()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += eVBFint + eVBFpar;
    
//  Linear contribution from Higgs self-coupling
    mu = mu + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    mu = mu + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();

    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}




double NPSMEFTd6::muVBFgamma(const double sqrt_s) const
{
    double mu = 1.0;
    
    double C1 = 0.0; //Use same values as VBF
    
    if (sqrt_s == 13.0) {
        
        C1 = 0.0064;

        mu += 
                +119630. * CiHbox / LambdaNP2
                -501300. * CiHQ3_11 / LambdaNP2
                -200890. * CiHD / LambdaNP2
                +11852.5 * CiHB / LambdaNP2
                -131586. * CiHW / LambdaNP2
                -361991. * CiHWB / LambdaNP2
                -18894.5 * CiDHB / LambdaNP2
                -69025.4 * CiDHW / LambdaNP2
                +23773.1 * CiW / LambdaNP2
                -4.629 * DeltaGF()
                -5.637 * deltaMwd6()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::muVBFgamma()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy. Use same as VBF.)
    mu += eVBFint + eVBFpar;
    
//  Linear contribution from Higgs self-coupling
    mu = mu + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    mu = mu + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();

    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::mueeWBF(const double sqrt_s) const
{
    double mu = 1.0;
    
    double C1 = 0.0;
    
    if (sqrt_s == 0.240) {
        
        C1 = 0.0064;

        mu += 
                +121120. * CiHbox / LambdaNP2
                -138682. * CiHL3_11 / LambdaNP2
                -203727. * CiHD / LambdaNP2
                -24699.7 * CiHW / LambdaNP2
                -379830. * CiHWB / LambdaNP2
                -18173.7 * CiDHW / LambdaNP2
                -4.716 * DeltaGF()
                -5.665 * deltaMwd6()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( 
                +3.307 * deltaMz()
                -3.995 * deltaMh()
                -0.486 * deltaaMZ()
                +3.507 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
          
    } else if (sqrt_s == 0.250) {
        
        C1 = 0.0064;

        mu += 
                +121142. * CiHbox / LambdaNP2
                -147357. * CiHL3_11 / LambdaNP2
                -203726. * CiHD / LambdaNP2
                -26559.2 * CiHW / LambdaNP2
                -379797. * CiHWB / LambdaNP2
                -19265.3 * CiDHW / LambdaNP2
                -4.717 * DeltaGF()
                -5.593 * deltaMwd6()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( 
                +3.413 * deltaMz()
                -3.644 * deltaMh()
                -0.502 * deltaaMZ()
                +3.523 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 0.350) {
        
        C1 = 0.0062;

        mu += 
                +121107. * CiHbox / LambdaNP2
                -219582. * CiHL3_11 / LambdaNP2
                -203717. * CiHD / LambdaNP2
                -39722.3 * CiHW / LambdaNP2
                -379795. * CiHWB / LambdaNP2
                -28864.2 * CiDHW / LambdaNP2
                -4.714 * DeltaGF()
                -5.13 * deltaMwd6()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( 
                +4.073 * deltaMz()
                -1.94 * deltaMh()
                -0.598 * deltaaMZ()
                +3.623 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }      
        
    } else if (sqrt_s == 0.365) {
        
        C1 = 0.0062; // Use the same as 350 GeV

        mu += 
                +121071. * CiHbox / LambdaNP2
                -228452. * CiHL3_11 / LambdaNP2
                -203725. * CiHD / LambdaNP2
                -40966.9 * CiHW / LambdaNP2
                -379798. * CiHWB / LambdaNP2
                -30110.4 * CiDHW / LambdaNP2
                -4.714 * DeltaGF()
                -5.08 * deltaMwd6()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( 
                +4.136 * deltaMz()
                -1.817 * deltaMh()
                -0.609 * deltaaMZ()
                +3.635 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 0.380) {
        
        C1 = 0.0062; // Use the same as 350 GeV

        mu += 
                +121001. * CiHbox / LambdaNP2
                -237126. * CiHL3_11 / LambdaNP2
                -203726. * CiHD / LambdaNP2
                -42070.9 * CiHW / LambdaNP2
                -379788. * CiHWB / LambdaNP2
                -31352.7 * CiDHW / LambdaNP2
                -4.714 * DeltaGF()
                -5.044 * deltaMwd6()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( 
                +4.192 * deltaMz()
                -1.711 * deltaMh()
                -0.618 * deltaaMZ()
                +3.64 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }

    } else if (sqrt_s == 0.500) {
        
        C1 = 0.0061;
        
        mu += 
                +121063. * CiHbox / LambdaNP2
                -295115. * CiHL3_11 / LambdaNP2
                -203679. * CiHD / LambdaNP2
                -47539.5 * CiHW / LambdaNP2
                -379773. * CiHWB / LambdaNP2
                -39825.1 * CiDHW / LambdaNP2
                -4.715 * DeltaGF()
                -4.817 * deltaMwd6()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( 
                +4.509 * deltaMz()
                -1.178 * deltaMh()
                -0.666 * deltaaMZ()
                +3.692 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 1.0) {
        
        C1 = 0.0059;

        mu += 
                +120960. * CiHbox / LambdaNP2
                -442647. * CiHL3_11 / LambdaNP2
                -203748. * CiHD / LambdaNP2
                -49375.4 * CiHW / LambdaNP2
                -379685. * CiHWB / LambdaNP2
                -63503.9 * CiDHW / LambdaNP2
                -4.712 * DeltaGF()
                -4.481 * deltaMwd6()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( 
                +4.99 * deltaMz()
                -0.582 * deltaMh()
                -0.734 * deltaaMZ()
                +3.765 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
  
    } else if (sqrt_s == 1.4) {
        
        C1 = 0.0058;

        mu += 
                +121118. * CiHbox / LambdaNP2
                -515189. * CiHL3_11 / LambdaNP2
                -203684. * CiHD / LambdaNP2
                -46619.5 * CiHW / LambdaNP2
                -379667. * CiHWB / LambdaNP2
                -75747.8 * CiDHW / LambdaNP2
                -4.714 * DeltaGF()
                -4.391 * deltaMwd6()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( 
                +5.13 * deltaMz()
                -0.446 * deltaMh()
                -0.754 * deltaaMZ()
                +3.784 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 1.5) {
        
        C1 = 0.0058;// Use the same as 1400 GeV

        mu += 
                +121200. * CiHbox / LambdaNP2
                -530152. * CiHL3_11 / LambdaNP2
                -203649. * CiHD / LambdaNP2
                -45921.3 * CiHW / LambdaNP2
                -379591. * CiHWB / LambdaNP2
                -78241.3 * CiDHW / LambdaNP2
                -4.715 * DeltaGF()
                -4.38 * deltaMwd6()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( 
                +5.154 * deltaMz()
                -0.424 * deltaMh()
                -0.757 * deltaaMZ()
                +3.786 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 3.0) {
        
        C1 = 0.0057;
        
        mu += 
                +121321. * CiHbox / LambdaNP2
                -684382. * CiHL3_11 / LambdaNP2
                -203585. * CiHD / LambdaNP2
                -38239. * CiHW / LambdaNP2
                -379518. * CiHWB / LambdaNP2
                -104465. * CiDHW / LambdaNP2
                -4.714 * DeltaGF()
                -4.258 * deltaMwd6()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( 
                +5.331 * deltaMz()
                -0.279 * deltaMh()
                -0.785 * deltaaMZ()
                +3.81 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }

    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::mueeWBF()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += eeeWBFint + eeeWBFpar;
    
//  Linear contribution from Higgs self-coupling
    mu = mu + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    mu = mu + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}


double NPSMEFTd6::mueeWBFPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{

//  Pure WBF, hence only initiated by LH fermions. No difference between polarizations at the linear level.
//  Expand like other functions when quadratic terms are included
    
    return mueeWBF(sqrt_s);
}

double NPSMEFTd6::mueeHvv(const double sqrt_s) const
{
    double mu = 1.0;
    
    double C1 = 0.0;
    
//  For the Higgs trilinear dependence assume the WBF mechanism dominates
    
    if (sqrt_s == 0.240) {
        
        C1 = 0.0064;

        mu += 
                +121539. * CiHbox / LambdaNP2
                +328845. * CiHL1_11 / LambdaNP2
                -37798.9 * CiHe_11 / LambdaNP2
                +279733. * CiHL3_11 / LambdaNP2
                -196039. * CiHD / LambdaNP2
                -70718.5 * CiHB / LambdaNP2
                +29671.9 * CiHW / LambdaNP2
                -401378. * CiHWB / LambdaNP2
                -23969.3 * CiDHB / LambdaNP2
                -1814.47 * CiDHW / LambdaNP2
                -4.698 * DeltaGF()
                -5.463 * deltaMwd6()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( 
                +4.842 * deltaMz()
                -2.535 * deltaMh()
                -0.528 * deltaaMZ()
                +3.46 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
          
    } else if (sqrt_s == 0.250) {
        
        C1 = 0.0064;

        mu += 
                +120627. * CiHbox / LambdaNP2
                +256825. * CiHL1_11 / LambdaNP2
                -38677.5 * CiHe_11 / LambdaNP2
                +175735. * CiHL3_11 / LambdaNP2
                -201059. * CiHD / LambdaNP2
                -57405. * CiHB / LambdaNP2
                -9860.82 * CiHW / LambdaNP2
                -403474. * CiHWB / LambdaNP2
                -20447.1 * CiDHB / LambdaNP2
                -9672.74 * CiDHW / LambdaNP2
                -4.656 * DeltaGF()
                -5.633 * deltaMwd6()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( 
                +4.194 * deltaMz()
                -2.783 * deltaMh()
                -0.477 * deltaaMZ()
                +3.414 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 0.350) {
        
        C1 = 0.0062;

        mu += 
                +120666. * CiHbox / LambdaNP2
                -19184.6 * CiHL1_11 / LambdaNP2
                -27432.1 * CiHe_11 / LambdaNP2
                -238244. * CiHL3_11 / LambdaNP2
                -204898. * CiHD / LambdaNP2
                +11833.5 * CiHB / LambdaNP2
                -94273.3 * CiHW / LambdaNP2
                -377703. * CiHWB / LambdaNP2
                +1111.63 * CiDHB / LambdaNP2
                -31735.2 * CiDHW / LambdaNP2
                -4.669 * DeltaGF()
                -5.329 * deltaMwd6()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( 
                +3.738 * deltaMz()
                -1.994 * deltaMh()
                -0.537 * deltaaMZ()
                +3.484 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }      
        
    } else if (sqrt_s == 0.365) {
        
        C1 = 0.0062; // Use the same as 350 GeV

        mu += 
                +120864. * CiHbox / LambdaNP2
                -24430. * CiHL1_11 / LambdaNP2
                -24398.7 * CiHe_11 / LambdaNP2
                -253414. * CiHL3_11 / LambdaNP2
                -204817. * CiHD / LambdaNP2
                +12826.5 * CiHB / LambdaNP2
                -93455. * CiHW / LambdaNP2
                -377489. * CiHWB / LambdaNP2
                +1693.48 * CiDHB / LambdaNP2
                -32834.7 * CiDHW / LambdaNP2
                -4.68 * DeltaGF()
                -5.265 * deltaMwd6()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( 
                +3.834 * deltaMz()
                -1.867 * deltaMh()
                -0.556 * deltaaMZ()
                +3.512 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 0.380) {
        
        C1 = 0.0062; // Use the same as 350 GeV

        mu += 
                +120775. * CiHbox / LambdaNP2
                -27548.7 * CiHL1_11 / LambdaNP2
                -22022.3 * CiHe_11 / LambdaNP2
                -266603. * CiHL3_11 / LambdaNP2
                -204782. * CiHD / LambdaNP2
                +13052.3 * CiHB / LambdaNP2
                -92560.2 * CiHW / LambdaNP2
                -377461. * CiHWB / LambdaNP2
                +1916.19 * CiDHB / LambdaNP2
                -33824.9 * CiDHW / LambdaNP2
                -4.684 * DeltaGF()
                -5.221 * deltaMwd6()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( 
                +3.931 * deltaMz()
                -1.75 * deltaMh()
                -0.574 * deltaaMZ()
                +3.532 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }

    } else if (sqrt_s == 0.500) {
        
        C1 = 0.0061;
        
        mu += 
                +120683. * CiHbox / LambdaNP2
                -26906.2 * CiHL1_11 / LambdaNP2
                -11055.8 * CiHe_11 / LambdaNP2
                -326940. * CiHL3_11 / LambdaNP2
                -204335. * CiHD / LambdaNP2
                +10505.8 * CiHB / LambdaNP2
                -82453.1 * CiHW / LambdaNP2
                -378407. * CiHWB / LambdaNP2
                +1889.64 * CiDHB / LambdaNP2
                -41332.3 * CiDHW / LambdaNP2
                -4.705 * DeltaGF()
                -4.943 * deltaMwd6()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( 
                +4.412 * deltaMz()
                -1.191 * deltaMh()
                -0.659 * deltaaMZ()
                +3.633 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 1.0) {
        
        C1 = 0.0059;

        mu += 
                +120462. * CiHbox / LambdaNP2
                -9025.99 * CiHL1_11 / LambdaNP2
                -3124.38 * CiHe_11 / LambdaNP2
                -454282. * CiHL3_11 / LambdaNP2
                -204077. * CiHD / LambdaNP2
                +3421.94 * CiHB / LambdaNP2
                -61892.5 * CiHW / LambdaNP2
                -379786. * CiHWB / LambdaNP2
                +396.747 * CiDHB / LambdaNP2
                -63826.6 * CiDHW / LambdaNP2
                -4.711 * DeltaGF()
                -4.587 * deltaMwd6()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( 
                +4.969 * deltaMz()
                -0.583 * deltaMh()
                -0.745 * deltaaMZ()
                +3.729 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
  
    } else if (sqrt_s == 1.4) {
        
        C1 = 0.0058;

        mu += 
                +120512. * CiHbox / LambdaNP2
                -4746.27 * CiHL1_11 / LambdaNP2
                -2212.55 * CiHe_11 / LambdaNP2
                -521829. * CiHL3_11 / LambdaNP2
                -204054. * CiHD / LambdaNP2
                +1891.37 * CiHB / LambdaNP2
                -54492.9 * CiHW / LambdaNP2
                -379916. * CiHWB / LambdaNP2
                +142.745 * CiDHB / LambdaNP2
                -75976. * CiDHW / LambdaNP2
                -4.712 * DeltaGF()
                -4.486 * deltaMwd6()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( 
                +5.108 * deltaMz()
                -0.447 * deltaMh()
                -0.767 * deltaaMZ()
                +3.751 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 1.5) {
        
        C1 = 0.0058;// Use the same as 1400 GeV

        mu += 
                +120512. * CiHbox / LambdaNP2
                -4105.67 * CiHL1_11 / LambdaNP2
                -2086.49 * CiHe_11 / LambdaNP2
                -536150. * CiHL3_11 / LambdaNP2
                -204072. * CiHD / LambdaNP2
                +1682.65 * CiHB / LambdaNP2
                -53138.1 * CiHW / LambdaNP2
                -379943. * CiHWB / LambdaNP2
                +134.612 * CiDHB / LambdaNP2
                -78546.2 * CiDHW / LambdaNP2
                -4.711 * DeltaGF()
                -4.469 * deltaMwd6()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( 
                +5.132 * deltaMz()
                -0.424 * deltaMh()
                -0.773 * deltaaMZ()
                +3.757 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 3.0) {
        
        C1 = 0.0057;
        
        mu += 
                +120404. * CiHbox / LambdaNP2
                -1215.14 * CiHL1_11 / LambdaNP2
                -1382.75 * CiHe_11 / LambdaNP2
                -686451. * CiHL3_11 / LambdaNP2
                -204039. * CiHD / LambdaNP2
                +293.31 * CiHB / LambdaNP2
                -41440.6 * CiHW / LambdaNP2
                -380130. * CiHWB / LambdaNP2
                -272.36 * CiDHB / LambdaNP2
                -104900. * CiDHW / LambdaNP2
                -4.706 * DeltaGF()
                -4.343 * deltaMwd6()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( 
                +5.307 * deltaMz()
                -0.283 * deltaMh()
                -0.802 * deltaaMZ()
                +3.789 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }

    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::mueeHvv()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += eeeWBFint + eeeWBFpar;
    
//  Linear contribution from Higgs self-coupling
    mu = mu + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    mu = mu + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}


double NPSMEFTd6::mueeHvvPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    double mu = 1.0;
    
    double C1 = 0.0;    
    
//  For the Higgs trilinear dependence assume the WBF mechanism dominates

    if (sqrt_s == 0.240) {
        
        C1 = 0.0064;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121180. * CiHbox / LambdaNP2
                +221479. * CiHL1_11 / LambdaNP2
                -508958. * CiHe_11 / LambdaNP2
                +220003. * CiHL3_11 / LambdaNP2
                -149238. * CiHD / LambdaNP2
                +24268.3 * CiHB / LambdaNP2
                -32411.5 * CiHW / LambdaNP2
                -194663. * CiHWB / LambdaNP2
                +29267.1 * CiDHB / LambdaNP2
                -11610.1 * CiDHW / LambdaNP2
                -3.633 * DeltaGF()
                -4.394 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +2.975 * deltaMz()
                -2.624 * deltaMh()
                +0.379 * deltaaMZ()
                +2.282 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121456. * CiHbox / LambdaNP2
                +337881. * CiHL1_11 / LambdaNP2
                +931.718 * CiHe_11 / LambdaNP2
                +283908. * CiHL3_11 / LambdaNP2
                -199920. * CiHD / LambdaNP2
                -78796.8 * CiHB / LambdaNP2
                +34606.7 * CiHW / LambdaNP2
                -418335. * CiHWB / LambdaNP2
                -28484. * CiDHB / LambdaNP2
                -1197.92 * CiDHW / LambdaNP2
                -4.781 * DeltaGF()
                -5.537 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +5.005 * deltaMz()
                -2.529 * deltaMh()
                -0.603 * deltaaMZ()
                +3.57 * deltaGmu() );
    
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121483. * CiHbox / LambdaNP2
                +266382. * CiHL1_11 / LambdaNP2
                -313151. * CiHe_11 / LambdaNP2
                +245682. * CiHL3_11 / LambdaNP2
                -168446. * CiHD / LambdaNP2
                -15072.1 * CiHB / LambdaNP2
                -6209.98 * CiHW / LambdaNP2
                -281195. * CiHWB / LambdaNP2
                +6468.72 * CiDHB / LambdaNP2
                -7633.09 * CiDHW / LambdaNP2
                -4.079 * DeltaGF()
                -4.832 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +3.758 * deltaMz()
                -2.579 * deltaMh()
                +0.009 * deltaaMZ()
                +2.778 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121500. * CiHbox / LambdaNP2
                +337280. * CiHL1_11 / LambdaNP2
                -1209.82 * CiHe_11 / LambdaNP2
                +283754. * CiHL3_11 / LambdaNP2
                -199723. * CiHD / LambdaNP2
                -78465.3 * CiHB / LambdaNP2
                +34393.4 * CiHW / LambdaNP2
                -417413. * CiHWB / LambdaNP2
                -28344.3 * CiDHB / LambdaNP2
                -1296.23 * CiDHW / LambdaNP2
                -4.777 * DeltaGF()
                -5.539 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.99 * deltaMz()
                -2.528 * deltaMh()
                -0.6 * deltaaMZ()
                +3.56 * deltaGmu() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeHvvPol()");
        }
        
    } else if (sqrt_s == 0.250) {
        
        C1 = 0.0064;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +120626. * CiHbox / LambdaNP2
                +172936. * CiHL1_11 / LambdaNP2
                -516799. * CiHe_11 / LambdaNP2
                +146366. * CiHL3_11 / LambdaNP2
                -156275. * CiHD / LambdaNP2
                +30993.1 * CiHB / LambdaNP2
                -62277.2 * CiHW / LambdaNP2
                -213096. * CiHWB / LambdaNP2
                +32593.7 * CiDHB / LambdaNP2
                -18479.4 * CiDHW / LambdaNP2
                -3.678 * DeltaGF()
                -4.598 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +2.739 * deltaMz()
                -2.661 * deltaMh()
                +0.356 * deltaaMZ()
                +2.343 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +120567. * CiHbox / LambdaNP2
                +263666. * CiHL1_11 / LambdaNP2
                -351.165 * CiHe_11 / LambdaNP2
                -396055. * CiHL3_11 / LambdaNP2
                -204612. * CiHD / LambdaNP2
                -64672.8 * CiHB / LambdaNP2
                -5618.64 * CiHW / LambdaNP2
                -418629. * CiHWB / LambdaNP2
                -24815.6 * CiDHB / LambdaNP2
                -9013.23 * CiDHW / LambdaNP2
                +286902. * CiLL_1221 / LambdaNP2
                -5.706 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.313 * deltaMz()
                -2.793 * deltaMh()
                -0.544 * deltaaMZ()
                +3.494 * deltaGmu() );
    
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +120240. * CiHbox / LambdaNP2
                +208124. * CiHL1_11 / LambdaNP2
                -315248. * CiHe_11 / LambdaNP2
                +158895. * CiHL3_11 / LambdaNP2
                -175074. * CiHD / LambdaNP2
                -6529.15 * CiHB / LambdaNP2
                -40099.4 * CiHW / LambdaNP2
                -293696. * CiHWB / LambdaNP2
                +10284.9 * CiDHB / LambdaNP2
                -15311.7 * CiDHW / LambdaNP2
                -4.092 * DeltaGF()
                -5.01 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +3.351 * deltaMz()
                -2.698 * deltaMh()
                -0.006 * deltaaMZ()
                +2.791 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +120459. * CiHbox / LambdaNP2
                +263262. * CiHL1_11 / LambdaNP2
                -2507.98 * CiHe_11 / LambdaNP2
                +177390. * CiHL3_11 / LambdaNP2
                -204514. * CiHD / LambdaNP2
                -64371.5 * CiHB / LambdaNP2
                -5927.95 * CiHW / LambdaNP2
                -417860. * CiHWB / LambdaNP2
                -24699.8 * CiDHB / LambdaNP2
                -9119.93 * CiDHW / LambdaNP2
                -4.726 * DeltaGF()
                -5.715 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.305 * deltaMz()
                -2.793 * deltaMh()
                -0.54 * deltaaMZ()
                +3.492 * deltaGmu() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeHvvPol()");
        }
        
    } else if (sqrt_s == 0.350) {
        
        C1 = 0.0062;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +120937. * CiHbox / LambdaNP2
                -41080.7 * CiHL1_11 / LambdaNP2
                -416801. * CiHe_11 / LambdaNP2
                -192794. * CiHL3_11 / LambdaNP2
                -182281. * CiHD / LambdaNP2
                +102909. * CiHB / LambdaNP2
                -87947.8 * CiHW / LambdaNP2
                -228111. * CiHWB / LambdaNP2
                +40181.7 * CiDHB / LambdaNP2
                -37530.5 * CiDHW / LambdaNP2
                -4.236 * DeltaGF()
                -4.832 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +3.177 * deltaMz()
                -1.894 * deltaMh()
                -0.171 * deltaaMZ()
                +3.022 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +120796. * CiHbox / LambdaNP2
                -17710.6 * CiHL1_11 / LambdaNP2
                -1357.61 * CiHe_11 / LambdaNP2
                -241114. * CiHL3_11 / LambdaNP2
                -206464. * CiHD / LambdaNP2
                +5738.97 * CiHB / LambdaNP2
                -94600.4 * CiHW / LambdaNP2
                -387581. * CiHWB / LambdaNP2
                -1403.89 * CiDHB / LambdaNP2
                -31363.8 * CiDHW / LambdaNP2
                -4.699 * DeltaGF()
                -5.361 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +3.768 * deltaMz()
                -2. * deltaMh()
                -0.556 * deltaaMZ()
                +3.512 * deltaGmu() );
    
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121065. * CiHbox / LambdaNP2
                -30567.4 * CiHL1_11 / LambdaNP2
                -235832. * CiHe_11 / LambdaNP2
                -213581. * CiHL3_11 / LambdaNP2
                -192620. * CiHD / LambdaNP2
                +60320.1 * CiHB / LambdaNP2
                -90446.2 * CiHW / LambdaNP2
                -297833. * CiHWB / LambdaNP2
                +22132.1 * CiDHB / LambdaNP2
                -34844.4 * CiDHW / LambdaNP2
                -4.439 * DeltaGF()
                -5.054 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +3.437 * deltaMz()
                -1.943 * deltaMh()
                -0.343 * deltaaMZ()
                +3.237 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +120725. * CiHbox / LambdaNP2
                -17741.9 * CiHL1_11 / LambdaNP2
                -2786.58 * CiHe_11 / LambdaNP2
                -241197. * CiHL3_11 / LambdaNP2
                -206387. * CiHD / LambdaNP2
                +6134.48 * CiHB / LambdaNP2
                -94603.3 * CiHW / LambdaNP2
                -387053. * CiHWB / LambdaNP2
                -1323.12 * CiDHB / LambdaNP2
                -31434.2 * CiDHW / LambdaNP2
                -4.696 * DeltaGF()
                -5.365 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +3.764 * deltaMz()
                -2. * deltaMh()
                -0.556 * deltaaMZ()
                +3.517 * deltaGmu() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeHvvPol()");
        } 
    
    } else if (sqrt_s == 0.365) {
        
        C1 = 0.0062; // Use the same as 350 GeV

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121120. * CiHbox / LambdaNP2
                -43274.8 * CiHL1_11 / LambdaNP2
                -379332. * CiHe_11 / LambdaNP2
                -213151. * CiHL3_11 / LambdaNP2
                -185704. * CiHD / LambdaNP2
                +95027.9 * CiHB / LambdaNP2
                -87042.2 * CiHW / LambdaNP2
                -246839. * CiHWB / LambdaNP2
                +37834.6 * CiDHB / LambdaNP2
                -38594.2 * CiDHW / LambdaNP2
                -4.314 * DeltaGF()
                -4.867 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +3.356 * deltaMz()
                -1.787 * deltaMh()
                -0.246 * deltaaMZ()
                +3.12 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +120708. * CiHbox / LambdaNP2
                -23163.4 * CiHL1_11 / LambdaNP2
                -1266.64 * CiHe_11 / LambdaNP2
                -256145. * CiHL3_11 / LambdaNP2
                -206112. * CiHD / LambdaNP2
                +7209.08 * CiHB / LambdaNP2
                -94095.3 * CiHW / LambdaNP2
                -386056. * CiHWB / LambdaNP2
                -673.745 * CiDHB / LambdaNP2
                -32528.4 * CiDHW / LambdaNP2
                -4.703 * DeltaGF()
                -5.297 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +3.865 * deltaMz()
                -1.869 * deltaMh()
                -0.577 * deltaaMZ()
                +3.533 * deltaGmu() );
    
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +120872. * CiHbox / LambdaNP2
                -34492.1 * CiHL1_11 / LambdaNP2
                -212361. * CiHe_11 / LambdaNP2
                -232050. * CiHL3_11 / LambdaNP2
                -194801. * CiHD / LambdaNP2
                +56353. * CiHB / LambdaNP2
                -90080.9 * CiHW / LambdaNP2
                -308151. * CiHWB / LambdaNP2
                +20707.2 * CiDHB / LambdaNP2
                -35840.6 * CiDHW / LambdaNP2
                -4.485 * DeltaGF()
                -5.033 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +3.586 * deltaMz()
                -1.817 * deltaMh()
                -0.393 * deltaaMZ()
                +3.287 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +120806. * CiHbox / LambdaNP2
                -23082.3 * CiHL1_11 / LambdaNP2
                -2521.89 * CiHe_11 / LambdaNP2
                -255807. * CiHL3_11 / LambdaNP2
                -205972. * CiHD / LambdaNP2
                +7600.7 * CiHB / LambdaNP2
                -94080.6 * CiHW / LambdaNP2
                -385587. * CiHWB / LambdaNP2
                -525.394 * CiDHB / LambdaNP2
                -32486.9 * CiDHW / LambdaNP2
                -4.703 * DeltaGF()
                -5.294 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +3.87 * deltaMz()
                -1.873 * deltaMh()
                -0.577 * deltaaMZ()
                +3.533 * deltaGmu() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeHvvPol()");
        } 
    
    } else if (sqrt_s == 0.380) {
        
        C1 = 0.0062; // Use the same as 350 GeV

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +120907. * CiHbox / LambdaNP2
                -43917.7 * CiHL1_11 / LambdaNP2
                -344628. * CiHe_11 / LambdaNP2
                -230932. * CiHL3_11 / LambdaNP2
                -188656. * CiHD / LambdaNP2
                +86802.5 * CiHB / LambdaNP2
                -86378.3 * CiHW / LambdaNP2
                -262732. * CiHWB / LambdaNP2
                +35211.7 * CiDHB / LambdaNP2
                -39122. * CiDHW / LambdaNP2
                -4.375 * DeltaGF()
                -4.833 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +3.526 * deltaMz()
                -1.675 * deltaMh()
                -0.322 * deltaaMZ()
                +3.202 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +120826. * CiHbox / LambdaNP2
                -26397.1 * CiHL1_11 / LambdaNP2
                -1156.51 * CiHe_11 / LambdaNP2
                -268680. * CiHL3_11 / LambdaNP2
                -205752. * CiHD / LambdaNP2
                +8226.72 * CiHB / LambdaNP2
                -92973.9 * CiHW / LambdaNP2
                -384868. * CiHWB / LambdaNP2
                -154.996 * CiDHB / LambdaNP2
                -33479.2 * CiDHW / LambdaNP2
                -4.706 * DeltaGF()
                -5.24 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +3.957 * deltaMz()
                -1.756 * deltaMh()
                -0.592 * deltaaMZ()
                +3.551 * deltaGmu() );
    
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121123. * CiHbox / LambdaNP2
                -35934.5 * CiHL1_11 / LambdaNP2
                -191922. * CiHe_11 / LambdaNP2
                -247636. * CiHL3_11 / LambdaNP2
                -196255. * CiHD / LambdaNP2
                +52143.1 * CiHB / LambdaNP2
                -89227.7 * CiHW / LambdaNP2
                -317018. * CiHWB / LambdaNP2
                +19725.8 * CiDHB / LambdaNP2
                -36723.5 * CiDHW / LambdaNP2
                -4.524 * DeltaGF()
                -5.007 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +3.729 * deltaMz()
                -1.706 * deltaMh()
                -0.439 * deltaaMZ()
                +3.366 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +120839. * CiHbox / LambdaNP2
                -26545. * CiHL1_11 / LambdaNP2
                -2293.44 * CiHe_11 / LambdaNP2
                -268673. * CiHL3_11 / LambdaNP2
                -205696. * CiHD / LambdaNP2
                +8476.41 * CiHB / LambdaNP2
                -92899.6 * CiHW / LambdaNP2
                -384414. * CiHWB / LambdaNP2
                +15.496 * CiDHB / LambdaNP2
                -33502.8 * CiDHW / LambdaNP2
                -4.704 * DeltaGF()
                -5.232 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +3.958 * deltaMz()
                -1.755 * deltaMh()
                -0.59 * deltaaMZ()
                +3.555 * deltaGmu() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeHvvPol()");
        } 
    
    } else if (sqrt_s == 0.500) {
        
        C1 = 0.0061;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +120734. * CiHbox / LambdaNP2
                -33626. * CiHL1_11 / LambdaNP2
                -177471. * CiHe_11 / LambdaNP2
                -312922. * CiHL3_11 / LambdaNP2
                -199388. * CiHD / LambdaNP2
                +44288.8 * CiHB / LambdaNP2
                -78960.3 * CiHW / LambdaNP2
                -332501. * CiHWB / LambdaNP2
                +20615.5 * CiDHB / LambdaNP2
                -43923.9 * CiDHW / LambdaNP2
                -4.614 * DeltaGF()
                -4.84 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.296 * deltaMz()
                -1.178 * deltaMh()
                -0.582 * deltaaMZ()
                +3.535 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +120746. * CiHbox / LambdaNP2
                -26369.8 * CiHL1_11 / LambdaNP2
                -905.141 * CiHe_11 / LambdaNP2
                -327709. * CiHL3_11 / LambdaNP2
                -204622. * CiHD / LambdaNP2
                +8508.33 * CiHB / LambdaNP2
                -82669.6 * CiHW / LambdaNP2
                -381185. * CiHWB / LambdaNP2
                +784.456 * CiDHB / LambdaNP2
                -41153.8 * CiDHW / LambdaNP2
                -4.711 * DeltaGF()
                -4.948 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.417 * deltaMz()
                -1.196 * deltaMh()
                -0.664 * deltaaMZ()
                +3.639 * deltaGmu() );
    
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +120667. * CiHbox / LambdaNP2
                -30480.6 * CiHL1_11 / LambdaNP2
                -96672.9 * CiHe_11 / LambdaNP2
                -320011. * CiHL3_11 / LambdaNP2
                -201855. * CiHD / LambdaNP2
                +27690.6 * CiHB / LambdaNP2
                -80770. * CiHW / LambdaNP2
                -355060. * CiHWB / LambdaNP2
                +11299.4 * CiDHB / LambdaNP2
                -42756.5 * CiDHW / LambdaNP2
                -4.656 * DeltaGF()
                -4.875 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.345 * deltaMz()
                -1.186 * deltaMh()
                -0.621 * deltaaMZ()
                +3.589 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +120715. * CiHbox / LambdaNP2
                -26433.4 * CiHL1_11 / LambdaNP2
                -1490.31 * CiHe_11 / LambdaNP2
                -327665. * CiHL3_11 / LambdaNP2
                -204644. * CiHD / LambdaNP2
                +8471.25 * CiHB / LambdaNP2
                -82673.2 * CiHW / LambdaNP2
                -381049. * CiHWB / LambdaNP2
                +862.813 * CiDHB / LambdaNP2
                -41179.7 * CiDHW / LambdaNP2
                -4.711 * DeltaGF()
                -4.942 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.416 * deltaMz()
                -1.194 * deltaMh()
                -0.664 * deltaaMZ()
                +3.64 * deltaGmu() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeHvvPol()");
        } 
    
    } else if (sqrt_s == 1.0) {
        
        C1 = 0.0059;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +120494. * CiHbox / LambdaNP2
                -9728.66 * CiHL1_11 / LambdaNP2
                -46166.9 * CiHe_11 / LambdaNP2
                -452752. * CiHL3_11 / LambdaNP2
                -203700. * CiHD / LambdaNP2
                +8561.22 * CiHB / LambdaNP2
                -61449.7 * CiHW / LambdaNP2
                -374076. * CiHWB / LambdaNP2
                +6473.98 * CiDHB / LambdaNP2
                -64032.3 * CiDHW / LambdaNP2
                -4.706 * DeltaGF()
                -4.581 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.956 * deltaMz()
                -0.583 * deltaMh()
                -0.739 * deltaaMZ()
                +3.723 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +120522. * CiHbox / LambdaNP2
                -8881.26 * CiHL1_11 / LambdaNP2
                -529.908 * CiHe_11 / LambdaNP2
                -454326. * CiHL3_11 / LambdaNP2
                -204057. * CiHD / LambdaNP2
                +3158.25 * CiHB / LambdaNP2
                -61850.9 * CiHW / LambdaNP2
                -380114. * CiHWB / LambdaNP2
                +63.589 * CiDHB / LambdaNP2
                -63800.9 * CiDHW / LambdaNP2
                -4.712 * DeltaGF()
                -4.587 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.967 * deltaMz()
                -0.582 * deltaMh()
                -0.746 * deltaaMZ()
                +3.731 * deltaGmu() );
    
        } else if (Pol_em == 80. && Pol_ep == -20.){
            mu += 
                +120541. * CiHbox / LambdaNP2
                -9598.71 * CiHL1_11 / LambdaNP2
                -37435. * CiHe_11 / LambdaNP2
                -453118. * CiHL3_11 / LambdaNP2
                -203771. * CiHD / LambdaNP2
                +7555.11 * CiHB / LambdaNP2
                -61524.6 * CiHW / LambdaNP2
                -375155. * CiHWB / LambdaNP2
                +5263.81 * CiDHB / LambdaNP2
                -64001.7 * CiDHW / LambdaNP2
                -4.706 * DeltaGF()
                -4.589 * deltaMwd6()
                ;     
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.959 * deltaMz()
                -0.583 * deltaMh()
                -0.741 * deltaaMZ()
                +3.726 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 20.){
            mu += 
                +120482. * CiHbox / LambdaNP2
                -8932.26 * CiHL1_11 / LambdaNP2
                -597.015 * CiHe_11 / LambdaNP2
                -454406. * CiHL3_11 / LambdaNP2
                -204110. * CiHD / LambdaNP2
                +3145.81 * CiHB / LambdaNP2
                -61837. * CiHW / LambdaNP2
                -380115. * CiHWB / LambdaNP2
                +45.924 * CiDHB / LambdaNP2
                -63834.7 * CiDHW / LambdaNP2
                -4.711 * DeltaGF()
                -4.588 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.968 * deltaMz()
                -0.582 * deltaMh()
                -0.746 * deltaaMZ()
                +3.73 * deltaGmu() );
    
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +120509. * CiHbox / LambdaNP2
                -9342.32 * CiHL1_11 / LambdaNP2
                -25028.5 * CiHe_11 / LambdaNP2
                -453487. * CiHL3_11 / LambdaNP2
                -203871. * CiHD / LambdaNP2
                +6021.71 * CiHB / LambdaNP2
                -61580. * CiHW / LambdaNP2
                -376790. * CiHWB / LambdaNP2
                +3494.08 * CiDHB / LambdaNP2
                -63959. * CiDHW / LambdaNP2
                -4.708 * DeltaGF()
                -4.589 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.962 * deltaMz()
                -0.582 * deltaMh()
                -0.742 * deltaaMZ()
                +3.726 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +120526. * CiHbox / LambdaNP2
                -8927.83 * CiHL1_11 / LambdaNP2
                -633.766 * CiHe_11 / LambdaNP2
                -454337. * CiHL3_11 / LambdaNP2
                -204073. * CiHD / LambdaNP2
                +3196.39 * CiHB / LambdaNP2
                -61833.5 * CiHW / LambdaNP2
                -380094. * CiHWB / LambdaNP2
                +82.665 * CiDHB / LambdaNP2
                -63817.5 * CiDHW / LambdaNP2
                -4.712 * DeltaGF()
                -4.588 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.967 * deltaMz()
                -0.582 * deltaMh()
                -0.746 * deltaaMZ()
                +3.731 * deltaGmu() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeHvvPol()");
        } 
    
    } else if (sqrt_s == 1.4) {
        
        C1 = 0.0058;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +120516. * CiHbox / LambdaNP2
                -5019.36 * CiHL1_11 / LambdaNP2
                -29937.8 * CiHe_11 / LambdaNP2
                -521211. * CiHL3_11 / LambdaNP2
                -203908. * CiHD / LambdaNP2
                +4153.08 * CiHB / LambdaNP2
                -54219.3 * CiHW / LambdaNP2
                -377548. * CiHWB / LambdaNP2
                +4509.78 * CiDHB / LambdaNP2
                -76054.8 * CiDHW / LambdaNP2
                -4.71 * DeltaGF()
                -4.484 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +5.105 * deltaMz()
                -0.447 * deltaMh()
                -0.765 * deltaaMZ()
                +3.747 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +120530. * CiHbox / LambdaNP2
                -4727.84 * CiHL1_11 / LambdaNP2
                -488.036 * CiHe_11 / LambdaNP2
                -521821. * CiHL3_11 / LambdaNP2
                -204045. * CiHD / LambdaNP2
                +1784.38 * CiHB / LambdaNP2
                -54507.5 * CiHW / LambdaNP2
                -380042. * CiHWB / LambdaNP2
                -122.009 * CiDHB / LambdaNP2
                -75950.5 * CiDHW / LambdaNP2
                -4.712 * DeltaGF()
                -4.487 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +5.108 * deltaMz()
                -0.447 * deltaMh()
                -0.768 * deltaaMZ()
                +3.749 * deltaGmu() );
    
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +120542. * CiHbox / LambdaNP2
                -4870.22 * CiHL1_11 / LambdaNP2
                -16376.8 * CiHe_11 / LambdaNP2
                -521472. * CiHL3_11 / LambdaNP2
                -203960. * CiHD / LambdaNP2
                +3068.42 * CiHB / LambdaNP2
                -54375.2 * CiHW / LambdaNP2
                -378699. * CiHWB / LambdaNP2
                +2390.51 * CiDHB / LambdaNP2
                -75996.8 * CiDHW / LambdaNP2
                -4.711 * DeltaGF()
                -4.485 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +5.107 * deltaMz()
                -0.448 * deltaMh()
                -0.766 * deltaaMZ()
                +3.749 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +120504. * CiHbox / LambdaNP2
                -4718.66 * CiHL1_11 / LambdaNP2
                -574.963 * CiHe_11 / LambdaNP2
                -521805. * CiHL3_11 / LambdaNP2
                -204053. * CiHD / LambdaNP2
                +1784.37 * CiHB / LambdaNP2
                -54482.7 * CiHW / LambdaNP2
                -380051. * CiHWB / LambdaNP2
                -99.132 * CiDHB / LambdaNP2
                -75974.5 * CiDHW / LambdaNP2
                -4.712 * DeltaGF()
                -4.487 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +5.107 * deltaMz()
                -0.447 * deltaMh()
                -0.767 * deltaaMZ()
                +3.749 * deltaGmu() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeHvvPol()");
        } 
        
    } else if (sqrt_s == 1.5) {
        
        C1 = 0.0058;// Use the same as 1400 GeV

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +120531. * CiHbox / LambdaNP2
                -4421.38 * CiHL1_11 / LambdaNP2
                -28114.2 * CiHe_11 / LambdaNP2
                -535633. * CiHL3_11 / LambdaNP2
                -203960. * CiHD / LambdaNP2
                +3556.32 * CiHB / LambdaNP2
                -52816.2 * CiHW / LambdaNP2
                -377932. * CiHWB / LambdaNP2
                +4253.17 * CiDHB / LambdaNP2
                -78599.6 * CiDHW / LambdaNP2
                -4.71 * DeltaGF()
                -4.465 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +5.128 * deltaMz()
                -0.424 * deltaMh()
                -0.772 * deltaaMZ()
                +3.755 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +120491. * CiHbox / LambdaNP2
                -4113.21 * CiHL1_11 / LambdaNP2
                -517.747 * CiHe_11 / LambdaNP2
                -536169. * CiHL3_11 / LambdaNP2
                -204050. * CiHD / LambdaNP2
                +1553.24 * CiHB / LambdaNP2
                -53097.9 * CiHW / LambdaNP2
                -380055. * CiHWB / LambdaNP2
                -129.437 * CiDHB / LambdaNP2
                -78539.4 * CiDHW / LambdaNP2
                -4.711 * DeltaGF()
                -4.468 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +5.131 * deltaMz()
                -0.424 * deltaMh()
                -0.773 * deltaaMZ()
                +3.755 * deltaGmu() );
    
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +120525. * CiHbox / LambdaNP2
                -4256.39 * CiHL1_11 / LambdaNP2
                -15376.9 * CiHe_11 / LambdaNP2
                -535845. * CiHL3_11 / LambdaNP2
                -203987. * CiHD / LambdaNP2
                +2641.32 * CiHB / LambdaNP2
                -53045.1 * CiHW / LambdaNP2
                -378920. * CiHWB / LambdaNP2
                +2237.55 * CiDHB / LambdaNP2
                -78549.8 * CiDHW / LambdaNP2
                -4.711 * DeltaGF()
                -4.468 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +5.129 * deltaMz()
                -0.424 * deltaMh()
                -0.772 * deltaaMZ()
                +3.753 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +120499. * CiHbox / LambdaNP2
                -4113.23 * CiHL1_11 / LambdaNP2
                -616.984 * CiHe_11 / LambdaNP2
                -536155. * CiHL3_11 / LambdaNP2
                -204035. * CiHD / LambdaNP2
                +1570.5 * CiHB / LambdaNP2
                -53079.3 * CiHW / LambdaNP2
                -380043. * CiHWB / LambdaNP2
                -112.179 * CiDHB / LambdaNP2
                -78543.9 * CiDHW / LambdaNP2
                -4.711 * DeltaGF()
                -4.468 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +5.13 * deltaMz()
                -0.424 * deltaMh()
                -0.773 * deltaaMZ()
                +3.755 * deltaGmu() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeHvvPol()");
        }
    
    } else if (sqrt_s == 3.0) {
        
        C1 = 0.0057;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +120384. * CiHbox / LambdaNP2
                -1301.85 * CiHL1_11 / LambdaNP2
                -16370.4 * CiHe_11 / LambdaNP2
                -686389. * CiHL3_11 / LambdaNP2
                -204031. * CiHD / LambdaNP2
                +628.479 * CiHB / LambdaNP2
                -41464.7 * CiHW / LambdaNP2
                -379766. * CiHWB / LambdaNP2
                +2259.53 * CiDHB / LambdaNP2
                -104941. * CiDHW / LambdaNP2
                -4.706 * DeltaGF()
                -4.342 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +5.306 * deltaMz()
                -0.283 * deltaMh()
                -0.802 * deltaaMZ()
                +3.787 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +120423. * CiHbox / LambdaNP2
                -1253.47 * CiHL1_11 / LambdaNP2
                -537.201 * CiHe_11 / LambdaNP2
                -686427. * CiHL3_11 / LambdaNP2
                -204047. * CiHD / LambdaNP2
                +268.601 * CiHB / LambdaNP2
                -41454. * CiHW / LambdaNP2
                -380141. * CiHWB / LambdaNP2
                -447.668 * CiDHB / LambdaNP2
                -104906. * CiDHW / LambdaNP2
                -4.707 * DeltaGF()
                -4.342 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +5.305 * deltaMz()
                -0.284 * deltaMh()
                -0.802 * deltaaMZ()
                +3.787 * deltaGmu() );
    
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +120399. * CiHbox / LambdaNP2
                -1267.47 * CiHL1_11 / LambdaNP2
                -9008.44 * CiHe_11 / LambdaNP2
                -686485. * CiHL3_11 / LambdaNP2
                -204052. * CiHD / LambdaNP2
                +439.947 * CiHB / LambdaNP2
                -41459.8 * CiHW / LambdaNP2
                -379947. * CiHWB / LambdaNP2
                +1005.59 * CiDHB / LambdaNP2
                -104927. * CiDHW / LambdaNP2
                -4.706 * DeltaGF()
                -4.342 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +5.303 * deltaMz()
                -0.283 * deltaMh()
                -0.802 * deltaaMZ()
                +3.789 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +120385. * CiHbox / LambdaNP2
                -1245.4 * CiHL1_11 / LambdaNP2
                -535.407 * CiHe_11 / LambdaNP2
                -686461. * CiHL3_11 / LambdaNP2
                -204048. * CiHD / LambdaNP2
                +244.425 * CiHB / LambdaNP2
                -41447.5 * CiHW / LambdaNP2
                -380150. * CiHWB / LambdaNP2
                -430.653 * CiDHB / LambdaNP2
                -104905. * CiDHW / LambdaNP2
                -4.706 * DeltaGF()
                -4.343 * deltaMwd6()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +5.307 * deltaMz()
                -0.283 * deltaMh()
                -0.802 * deltaaMZ()
                +3.789 * deltaGmu() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeHvvPol()");
        } 
    
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::mueeHvvPol()");   
    
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += eeeWBFint + eeeWBFpar;
    
//  Linear contribution from Higgs self-coupling
    mu = mu + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    mu = mu + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::mueeZBF(const double sqrt_s) const
{
    double mu = 1.0;
    
    double C1 = 0.0;
    
    if (sqrt_s == 0.240) {
        
    C1 = 0.0070;

        mu += 
                +121661. * CiHbox / LambdaNP2
                +489617. * CiHL1_11 / LambdaNP2
                -357163. * CiHe_11 / LambdaNP2
                +489617. * CiHL3_11 / LambdaNP2
                -39217.8 * CiHD / LambdaNP2
                +1525468. * CiHB / LambdaNP2
                +378019. * CiHW / LambdaNP2
                +215983. * CiHWB / LambdaNP2
                -6554.11 * CiDHB / LambdaNP2
                +1175.47 * CiDHW / LambdaNP2
                -3.161 * DeltaGF()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( +0.908 * deltaMz()
                -5.799 * deltaMh()
                -0.248 * deltaaMZ()
                +3.158 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 0.250) {
        
    C1 = 0.0070;

        mu += 
                +122144. * CiHbox / LambdaNP2
                +444406. * CiHL1_11 / LambdaNP2
                -315727. * CiHe_11 / LambdaNP2
                +444406. * CiHL3_11 / LambdaNP2
                -41440.8 * CiHD / LambdaNP2
                +1186855. * CiHB / LambdaNP2
                +301913. * CiHW / LambdaNP2
                +98540.5 * CiHWB / LambdaNP2
                -5766.35 * CiDHB / LambdaNP2
                +294.724 * CiDHW / LambdaNP2
                -3.279 * DeltaGF()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( +2.044 * deltaMz()
                -4.578 * deltaMh()
                -0.341 * deltaaMZ()
                +3.283 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 0.350) {
        
    C1 = 0.0069;

        mu += 
                +121556. * CiHbox / LambdaNP2
                +46354.9 * CiHL1_11 / LambdaNP2
                -251.929 * CiHe_11 / LambdaNP2
                +46354.9 * CiHL3_11 / LambdaNP2
                -43426.2 * CiHD / LambdaNP2
                +450512. * CiHB / LambdaNP2
                +166493. * CiHW / LambdaNP2
                -198898. * CiHWB / LambdaNP2
                -4408.76 * CiDHB / LambdaNP2
                -17005.2 * CiDHW / LambdaNP2
                -3.427 * DeltaGF()
                ; 
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( +3.845 * deltaMz()
                -1.857 * deltaMh()
                -0.423 * deltaaMZ()
                +3.407 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 0.365) {
        
    C1 = 0.0069; // use same as 350 GeV

        mu += 
                +121067. * CiHbox / LambdaNP2
                +9887.64 * CiHL1_11 / LambdaNP2
                +27809. * CiHe_11 / LambdaNP2
                +9887.64 * CiHL3_11 / LambdaNP2
                -43174.2 * CiHD / LambdaNP2
                +417865. * CiHB / LambdaNP2
                +154270. * CiHW / LambdaNP2
                -201517. * CiHWB / LambdaNP2
                -4943.82 * CiDHB / LambdaNP2
                -19213.5 * CiDHW / LambdaNP2
                -3.423 * DeltaGF()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( +3.861 * deltaMz()
                -1.736 * deltaMh()
                -0.426 * deltaaMZ()
                +3.375 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 0.380) {
        
    C1 = 0.0069; // use same as 350 GeV

        mu += 
                +121214. * CiHbox / LambdaNP2
                -22289.7 * CiHL1_11 / LambdaNP2
                +52903.2 * CiHe_11 / LambdaNP2
                -22289.7 * CiHL3_11 / LambdaNP2
                -43137.3 * CiHD / LambdaNP2
                +388336. * CiHB / LambdaNP2
                +140923. * CiHW / LambdaNP2
                -202884. * CiHWB / LambdaNP2
                -5363.69 * CiDHB / LambdaNP2
                -21404.2 * CiDHW / LambdaNP2
                -3.418 * DeltaGF()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( +3.887 * deltaMz()
                -1.633 * deltaMh()
                -0.419 * deltaaMZ()
                +3.393 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 0.500) {
        
    C1 = 0.0067;

        mu += 
                +121453. * CiHbox / LambdaNP2
                -185326. * CiHL1_11 / LambdaNP2
                +178925. * CiHe_11 / LambdaNP2
                -185326. * CiHL3_11 / LambdaNP2
                -42051.6 * CiHD / LambdaNP2
                +236945. * CiHB / LambdaNP2
                +67833.5 * CiHW / LambdaNP2
                -178623. * CiHWB / LambdaNP2
                -8004.61 * CiDHB / LambdaNP2
                -33567.3 * CiDHW / LambdaNP2
                -3.416 * DeltaGF()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( +3.963 * deltaMz()
                -1.143 * deltaMh()
                -0.408 * deltaaMZ()
                +3.383 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 1.0) {
        
    C1 = 0.0065;

        mu += 
                +121062. * CiHbox / LambdaNP2
                -409543. * CiHL1_11 / LambdaNP2
                +356730. * CiHe_11 / LambdaNP2
                -409543. * CiHL3_11 / LambdaNP2
                -42133.9 * CiHD / LambdaNP2
                +69851. * CiHB / LambdaNP2
                -14416.8 * CiHW / LambdaNP2
                -113198. * CiHWB / LambdaNP2
                -18688.4 * CiDHB / LambdaNP2
                -61696. * CiDHW / LambdaNP2
                -3.405 * DeltaGF()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( +4.216 * deltaMz()
                -0.546 * deltaMh()
                -0.407 * deltaaMZ()
                +3.393 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 1.4) {
        
    C1 = 0.0065;

        mu += 
                +120749. * CiHbox / LambdaNP2
                -493617. * CiHL1_11 / LambdaNP2
                +426669. * CiHe_11 / LambdaNP2
                -493617. * CiHL3_11 / LambdaNP2
                -42486.9 * CiHD / LambdaNP2
                +34633.1 * CiHB / LambdaNP2
                -27609.6 * CiHW / LambdaNP2
                -97014.2 * CiHWB / LambdaNP2
                -23942.2 * CiDHB / LambdaNP2
                -74940.3 * CiDHW / LambdaNP2
                -3.405 * DeltaGF()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( +4.309 * deltaMz()
                -0.422 * deltaMh()
                -0.402 * deltaaMZ()
                +3.379 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 1.5) {
        
        C1 = 0.0065;// Use the same as 1400 GeV

        mu += 
                +120587. * CiHbox / LambdaNP2
                -510290. * CiHL1_11 / LambdaNP2
                +440504. * CiHe_11 / LambdaNP2
                -510290. * CiHL3_11 / LambdaNP2
                -42529.6 * CiHD / LambdaNP2
                +30448.1 * CiHB / LambdaNP2
                -30741.2 * CiHW / LambdaNP2
                -95903.3 * CiHWB / LambdaNP2
                -25074.9 * CiDHB / LambdaNP2
                -77634.5 * CiDHW / LambdaNP2
                -3.401 * DeltaGF()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( +4.326 * deltaMz()
                -0.4 * deltaMh()
                -0.403 * deltaaMZ()
                +3.37 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 3.0) {
        
    C1 = 0.0063;
        
        mu += 
                +120474. * CiHbox / LambdaNP2
                -677185. * CiHL1_11 / LambdaNP2
                +582037. * CiHe_11 / LambdaNP2
                -677185. * CiHL3_11 / LambdaNP2
                -42541.3 * CiHD / LambdaNP2
                +6810.6 * CiHB / LambdaNP2
                -32994.5 * CiHW / LambdaNP2
                -78012.3 * CiHWB / LambdaNP2
                -36250. * CiDHB / LambdaNP2
                -105734. * CiDHW / LambdaNP2
                -3.405 * DeltaGF()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( +4.463 * deltaMz()
                -0.265 * deltaMh()
                -0.405 * deltaaMZ()
                +3.351 * deltaGmu() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }

    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZBF()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    //(Assume similar to WBF.)
    mu += eeeWBFint + eeeWBFpar;
    
//  Linear contribution from Higgs self-coupling
    mu = mu + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    mu = mu + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}


double NPSMEFTd6::mueeZBFPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    double mu = 1.0;
    
    double C1 = 0.0;

    if (sqrt_s == 0.240) {
        
        C1 = 0.0070;
        
        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121531. * CiHbox / LambdaNP2
                +58943.5 * CiHL1_11 / LambdaNP2
                -939512. * CiHe_11 / LambdaNP2
                +58943.5 * CiHL3_11 / LambdaNP2
                +77442.6 * CiHD / LambdaNP2
                +2082256. * CiHB / LambdaNP2
                +108043. * CiHW / LambdaNP2
                +1362693. * CiHWB / LambdaNP2
                +40385. * CiDHB / LambdaNP2
                -21886. * CiDHW / LambdaNP2
                +0.563 * DeltaGF()
                ;            

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -6.582 * deltaMz()
                -5.732 * deltaMh()
                +3.573 * deltaaMZ()
                -0.708 * deltaGmu() );
            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +122065. * CiHbox / LambdaNP2
                +905327. * CiHL1_11 / LambdaNP2
                -55689. * CiHe_11 / LambdaNP2
                +905327. * CiHL3_11 / LambdaNP2
                -124548. * CiHD / LambdaNP2
                +905057. * CiHB / LambdaNP2
                +540185. * CiHW / LambdaNP2
                -329708. * CiHWB / LambdaNP2
                -37296.9 * CiDHB / LambdaNP2
                +20497.1 * CiDHW / LambdaNP2
                -5.854 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +6.473 * deltaMz()
                -5.971 * deltaMh()
                -3.019 * deltaaMZ()
                +5.959 * deltaGmu() );
            
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121947. * CiHbox / LambdaNP2
                +88774.4 * CiHL1_11 / LambdaNP2
                -753269. * CiHe_11 / LambdaNP2
                +88774.4 * CiHL3_11 / LambdaNP2
                +54593.2 * CiHD / LambdaNP2
                +2101955. * CiHB / LambdaNP2
                +182237. * CiHW / LambdaNP2
                +972861. * CiHWB / LambdaNP2
                +29346.2 * CiDHB / LambdaNP2
                -18562.1 * CiDHW / LambdaNP2
                -0.206 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -5.131 * deltaMz()
                -5.658 * deltaMh()
                +2.794 * deltaaMZ()
                +0.082 * deltaGmu() );
            
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +122265. * CiHbox / LambdaNP2
                +785643. * CiHL1_11 / LambdaNP2
                -66907.6 * CiHe_11 / LambdaNP2
                +785643. * CiHL3_11 / LambdaNP2
                -107673. * CiHD / LambdaNP2
                +1115316. * CiHB / LambdaNP2
                +521873. * CiHW / LambdaNP2
                -331727. * CiHWB / LambdaNP2
                -32442.4 * CiDHB / LambdaNP2
                +15348.7 * CiDHW / LambdaNP2
                -5.334 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +5.367 * deltaMz()
                -5.87 * deltaMh()
                -2.491 * deltaaMZ()
                +5.409 * deltaGmu() );
            
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZBFPol()");
        }
        
    } else if (sqrt_s == 0.250) {
        
        C1 = 0.0070;
        
        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121054. * CiHbox / LambdaNP2
                +51113. * CiHL1_11 / LambdaNP2
                -851357. * CiHe_11 / LambdaNP2
                +51113. * CiHL3_11 / LambdaNP2
                +76762.9 * CiHD / LambdaNP2
                +1629614. * CiHB / LambdaNP2
                +72741.6 * CiHW / LambdaNP2
                +1130834. * CiHWB / LambdaNP2
                +34381.7 * CiDHB / LambdaNP2
                -19876.5 * CiDHW / LambdaNP2
                +0.563 * DeltaGF()
                ;     

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -5.658 * deltaMz()
                -4.485 * deltaMh()
                +3.577 * deltaaMZ()
                -0.638 * deltaGmu() );
            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121471. * CiHbox / LambdaNP2
                +824294. * CiHL1_11 / LambdaNP2
                -45066.5 * CiHe_11 / LambdaNP2
                +824294. * CiHL3_11 / LambdaNP2
                -128864. * CiHD / LambdaNP2
                +644513. * CiHB / LambdaNP2
                +425051. * CiHW / LambdaNP2
                -383720. * CiHWB / LambdaNP2
                -32434.3 * CiDHB / LambdaNP2
                +15329.4 * CiDHW / LambdaNP2
                -6.022 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +7.852 * deltaMz()
                -4.536 * deltaMh()
                -3.165 * deltaaMZ()
                +6.136 * deltaGmu() );
            
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121494. * CiHbox / LambdaNP2
                +77372.1 * CiHL1_11 / LambdaNP2
                -676199. * CiHe_11 / LambdaNP2
                +77372.1 * CiHL3_11 / LambdaNP2
                +53294.7 * CiHD / LambdaNP2
                +1668830. * CiHB / LambdaNP2
                +145010. * CiHW / LambdaNP2
                +772902. * CiHWB / LambdaNP2
                +23910.6 * CiDHB / LambdaNP2
                -16890.6 * CiDHW / LambdaNP2
                -0.226 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -4.183 * deltaMz()
                -4.557 * deltaMh()
                +2.773 * deltaaMZ()
                +0.148 * deltaGmu() );
            
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121947. * CiHbox / LambdaNP2
                +713174. * CiHL1_11 / LambdaNP2
                -53393.3 * CiHe_11 / LambdaNP2
                +713174. * CiHL3_11 / LambdaNP2
                -111120. * CiHD / LambdaNP2
                +843388. * CiHB / LambdaNP2
                +417838. * CiHW / LambdaNP2
                -386753. * CiHWB / LambdaNP2
                -27915.7 * CiDHB / LambdaNP2
                +11946.5 * CiDHW / LambdaNP2
                -5.496 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +6.641 * deltaMz()
                -4.576 * deltaMh()
                -2.605 * deltaaMZ()
                +5.56 * deltaGmu() );
            
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZBFPol()");
        }
        
    } else if (sqrt_s == 0.350) {
        
        C1 = 0.0069;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121674. * CiHbox / LambdaNP2
                -47420.2 * CiHL1_11 / LambdaNP2
                -172088. * CiHe_11 / LambdaNP2
                -47420.2 * CiHL3_11 / LambdaNP2
                +59728. * CiHD / LambdaNP2
                +544205. * CiHB / LambdaNP2
                +83604.4 * CiHW / LambdaNP2
                +435393. * CiHWB / LambdaNP2
                -24800.4 * CiDHB / LambdaNP2
                -4583.09 * CiDHW / LambdaNP2
                -0.05 * DeltaGF()
                ;           

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.905 * deltaMz()
                -1.842 * deltaMh()
                +2.966 * deltaaMZ()
                +0.009 * deltaGmu() );
            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121541. * CiHbox / LambdaNP2
                +197618. * CiHL1_11 / LambdaNP2
                +42238.9 * CiHe_11 / LambdaNP2
                +197618. * CiHL3_11 / LambdaNP2
                -124376. * CiHD / LambdaNP2
                +181154. * CiHB / LambdaNP2
                +195329. * CiHW / LambdaNP2
                -505800. * CiHWB / LambdaNP2
                +13082.6 * CiDHB / LambdaNP2
                -26607.4 * CiDHW / LambdaNP2
                -6.096 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +9.303 * deltaMz()
                -1.82 * deltaMh()
                -3.105 * deltaaMZ()
                +6.071 * deltaGmu() );
            
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121760. * CiHbox / LambdaNP2
                -62853. * CiHL1_11 / LambdaNP2
                -83019.6 * CiHe_11 / LambdaNP2
                -62853. * CiHL3_11 / LambdaNP2
                +34395.4 * CiHD / LambdaNP2
                +623389. * CiHB / LambdaNP2
                +123932. * CiHW / LambdaNP2
                +181789. * CiHWB / LambdaNP2
                -20420. * CiDHB / LambdaNP2
                -7820.42 * CiDHW / LambdaNP2
                -0.875 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -1.322 * deltaMz()
                -1.873 * deltaMh()
                +2.14 * deltaaMZ()
                +0.844 * deltaGmu() );
            
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121557. * CiHbox / LambdaNP2
                +131443. * CiHL1_11 / LambdaNP2
                +63326.7 * CiHe_11 / LambdaNP2
                +131443. * CiHL3_11 / LambdaNP2
                -103038. * CiHD / LambdaNP2
                +323596. * CiHB / LambdaNP2
                +201676. * CiHW / LambdaNP2
                -491019. * CiHWB / LambdaNP2
                +7992.43 * CiDHB / LambdaNP2
                -24283.6 * CiDHW / LambdaNP2
                -5.391 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +7.818 * deltaMz()
                -1.846 * deltaMh()
                -2.402 * deltaaMZ()
                +5.358 * deltaGmu() );
            
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZBFPol()");
        } 
    
    } else if (sqrt_s == 0.365) {
        
        C1 = 0.0069; // Use same as 350 GeV

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121458. * CiHbox / LambdaNP2
                -58695.1 * CiHL1_11 / LambdaNP2
                -109686. * CiHe_11 / LambdaNP2
                -58695.1 * CiHL3_11 / LambdaNP2
                +58496.7 * CiHD / LambdaNP2
                +489137. * CiHB / LambdaNP2
                +80751.3 * CiHW / LambdaNP2
                +410304. * CiHWB / LambdaNP2
                -30918.3 * CiDHB / LambdaNP2
                -3571.31 * CiDHW / LambdaNP2
                -0.085 * DeltaGF()
                ;            

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.809 * deltaMz()
                -1.721 * deltaMh()
                +2.93 * deltaaMZ()
                +0.026 * deltaGmu() );
            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121152. * CiHbox / LambdaNP2
                +136019. * CiHL1_11 / LambdaNP2
                +50762. * CiHe_11 / LambdaNP2
                +136019. * CiHL3_11 / LambdaNP2
                -123859. * CiHD / LambdaNP2
                +165799. * CiHB / LambdaNP2
                +176652. * CiHW / LambdaNP2
                -504889. * CiHWB / LambdaNP2
                +16920.7 * CiDHB / LambdaNP2
                -31414.1 * CiDHW / LambdaNP2
                -6.076 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +9.271 * deltaMz()
                -1.7 * deltaMh()
                -3.092 * deltaaMZ()
                +6.031 * deltaGmu() );
            
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121193. * CiHbox / LambdaNP2
                -76905.7 * CiHL1_11 / LambdaNP2
                -32264.3 * CiHe_11 / LambdaNP2
                -76905.7 * CiHL3_11 / LambdaNP2
                +33650.3 * CiHD / LambdaNP2
                +573505. * CiHB / LambdaNP2
                +117937. * CiHW / LambdaNP2
                +166382. * CiHWB / LambdaNP2
                -25012.1 * CiDHB / LambdaNP2
                -7703.47 * CiDHW / LambdaNP2
                -0.911 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -1.233 * deltaMz()
                -1.746 * deltaMh()
                +2.101 * deltaaMZ()
                +0.861 * deltaGmu() );
            
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121177. * CiHbox / LambdaNP2
                +77981.5 * CiHL1_11 / LambdaNP2
                +74274.1 * CiHe_11 / LambdaNP2
                +77981.5 * CiHL3_11 / LambdaNP2
                -102068. * CiHD / LambdaNP2
                +305730. * CiHB / LambdaNP2
                +183682. * CiHW / LambdaNP2
                -487770. * CiHWB / LambdaNP2
                +10624.8 * CiDHB / LambdaNP2
                -28092.3 * CiDHW / LambdaNP2
                -5.366 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +7.791 * deltaMz()
                -1.726 * deltaMh()
                -2.377 * deltaaMZ()
                +5.325 * deltaGmu() );
            
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZBFPol()");
        } 
    
    } else if (sqrt_s == 0.380) {
        
        C1 = 0.0069; // Use same as 350 GeV

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121392. * CiHbox / LambdaNP2
                -68799.8 * CiHL1_11 / LambdaNP2
                -54383.2 * CiHe_11 / LambdaNP2
                -68799.8 * CiHL3_11 / LambdaNP2
                +57427.7 * CiHD / LambdaNP2
                +439155. * CiHB / LambdaNP2
                +76978.2 * CiHW / LambdaNP2
                +392293. * CiHWB / LambdaNP2
                -36175.9 * CiDHB / LambdaNP2
                -3193.74 * CiDHW / LambdaNP2
                -0.11 * DeltaGF()
                ;           

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.74 * deltaMz()
                -1.62 * deltaMh()
                +2.907 * deltaaMZ()
                +0.079 * deltaGmu() );
            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121306. * CiHbox / LambdaNP2
                +80159.7 * CiHL1_11 / LambdaNP2
                +58002.2 * CiHe_11 / LambdaNP2
                +80159.7 * CiHL3_11 / LambdaNP2
                -123524. * CiHD / LambdaNP2
                +151617. * CiHB / LambdaNP2
                +154342. * CiHW / LambdaNP2
                -500961. * CiHWB / LambdaNP2
                +20509.9 * CiDHB / LambdaNP2
                -35718.1 * CiDHW / LambdaNP2
                -6.064 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +9.254 * deltaMz()
                -1.608 * deltaMh()
                -3.07 * deltaaMZ()
                +6.04 * deltaGmu() );
            
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121171. * CiHbox / LambdaNP2
                -89494.3 * CiHL1_11 / LambdaNP2
                +11882.3 * CiHe_11 / LambdaNP2
                -89494.3 * CiHL3_11 / LambdaNP2
                +32430.1 * CiHD / LambdaNP2
                +524620. * CiHB / LambdaNP2
                +111520. * CiHW / LambdaNP2
                +156122. * CiHWB / LambdaNP2
                -29271.1 * CiDHB / LambdaNP2
                -8056.8 * CiDHW / LambdaNP2
                -0.928 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -1.145 * deltaMz()
                -1.643 * deltaMh()
                +2.077 * deltaaMZ()
                +0.898 * deltaGmu() );
            
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121286. * CiHbox / LambdaNP2
                +30046.7 * CiHL1_11 / LambdaNP2
                +84014. * CiHe_11 / LambdaNP2
                +30046.7 * CiHL3_11 / LambdaNP2
                -101539. * CiHD / LambdaNP2
                +286981. * CiHB / LambdaNP2
                +164662. * CiHW / LambdaNP2
                -480410. * CiHWB / LambdaNP2
                +13149.6 * CiDHB / LambdaNP2
                -31886.7 * CiDHW / LambdaNP2
                -5.346 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +7.766 * deltaMz()
                -1.629 * deltaMh()
                -2.353 * deltaaMZ()
                +5.316 * deltaGmu() );
            
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZBFPol()");
        } 
    
    } else if (sqrt_s == 0.500) {
        
        C1 = 0.0067;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121372. * CiHbox / LambdaNP2
                -121062. * CiHL1_11 / LambdaNP2
                +224754. * CiHe_11 / LambdaNP2
                -121062. * CiHL3_11 / LambdaNP2
                +55161.7 * CiHD / LambdaNP2
                +201238. * CiHB / LambdaNP2
                +52456.6 * CiHW / LambdaNP2
                +335517. * CiHWB / LambdaNP2
                -63733.4 * CiDHB / LambdaNP2
                -2379.21 * CiDHW / LambdaNP2
                -0.207 * DeltaGF()
                ;            

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.453 * deltaMz()
                -1.136 * deltaMh()
                +2.81 * deltaaMZ()
                +0.175 * deltaGmu() );
            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121399. * CiHbox / LambdaNP2
                -200849. * CiHL1_11 / LambdaNP2
                +96427.7 * CiHe_11 / LambdaNP2
                -200849. * CiHL3_11 / LambdaNP2
                -121178. * CiHD / LambdaNP2
                +83220.9 * CiHB / LambdaNP2
                +42832.2 * CiHW / LambdaNP2
                -464173. * CiHWB / LambdaNP2
                +37654.2 * CiDHB / LambdaNP2
                -59029.6 * CiDHW / LambdaNP2
                -6.025 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +9.205 * deltaMz()
                -1.133 * deltaMh()
                -3.019 * deltaaMZ()
                +5.99 * deltaGmu() );
            
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121435. * CiHbox / LambdaNP2
                -154953. * CiHL1_11 / LambdaNP2
                +235326. * CiHe_11 / LambdaNP2
                -154953. * CiHL3_11 / LambdaNP2
                +30472. * CiHD / LambdaNP2
                +298145. * CiHB / LambdaNP2
                +75047.6 * CiHW / LambdaNP2
                +137304. * CiHWB / LambdaNP2
                -49636.1 * CiDHB / LambdaNP2
                -10277.1 * CiDHW / LambdaNP2
                -1.027 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -0.829 * deltaMz()
                -1.142 * deltaMh()
                +1.988 * deltaaMZ()
                +0.989 * deltaGmu() );
            
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121468. * CiHbox / LambdaNP2
                -208577. * CiHL1_11 / LambdaNP2
                +134790. * CiHe_11 / LambdaNP2
                -208577. * CiHL3_11 / LambdaNP2
                -98708.1 * CiHD / LambdaNP2
                +190310. * CiHB / LambdaNP2
                +62321.4 * CiHW / LambdaNP2
                -429412. * CiHWB / LambdaNP2
                +24628.2 * CiDHB / LambdaNP2
                -51722.9 * CiDHW / LambdaNP2
                -5.287 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +7.714 * deltaMz()
                -1.14 * deltaMh()
                -2.279 * deltaaMZ()
                +5.251 * deltaGmu() );
            
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZBFPol()");
        } 
    
    } else if (sqrt_s == 1.0) {
        
        C1 = 0.0065;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121044. * CiHbox / LambdaNP2
                -206156. * CiHL1_11 / LambdaNP2
                +586357. * CiHe_11 / LambdaNP2
                -206156. * CiHL3_11 / LambdaNP2
                +54157.3 * CiHD / LambdaNP2
                -30839.6 * CiHB / LambdaNP2
                +18110.3 * CiHW / LambdaNP2
                +345253. * CiHWB / LambdaNP2
                -108488. * CiDHB / LambdaNP2
                -12324.2 * CiDHW / LambdaNP2
                -0.229 * DeltaGF()
                ;         

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.141 * deltaMz()
                -0.544 * deltaMh()
                +2.775 * deltaaMZ()
                +0.211 * deltaGmu() );
            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121085. * CiHbox / LambdaNP2
                -565700. * CiHL1_11 / LambdaNP2
                +157498. * CiHe_11 / LambdaNP2
                -565700. * CiHL3_11 / LambdaNP2
                -120795. * CiHD / LambdaNP2
                +7953.6 * CiHB / LambdaNP2
                -79908.9 * CiHW / LambdaNP2
                -402278. * CiHWB / LambdaNP2
                +54805.3 * CiDHB / LambdaNP2
                -101988. * CiDHW / LambdaNP2
                -6.001 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +9.412 * deltaMz()
                -0.546 * deltaMh()
                -3.005 * deltaaMZ()
                +5.986 * deltaGmu() );
            
        } else if (Pol_em == 80. && Pol_ep == -20.){
            mu += 
                +121091. * CiHbox / LambdaNP2
                -225779. * CiHL1_11 / LambdaNP2
                +568149. * CiHe_11 / LambdaNP2
                -225779. * CiHL3_11 / LambdaNP2
                +45736.7 * CiHD / LambdaNP2
                +2164.38 * CiHB / LambdaNP2
                +20504.6 * CiHW / LambdaNP2
                +290141. * CiHWB / LambdaNP2
                -100416. * CiDHB / LambdaNP2
                -16574.6 * CiDHW / LambdaNP2
                -0.51 * DeltaGF()
                ;     
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -1.569 * deltaMz()
                -0.555 * deltaMh()
                +2.507 * deltaaMZ()
                +0.493 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 20.){
            mu += 
                +121091. * CiHbox / LambdaNP2
                -552286. * CiHL1_11 / LambdaNP2
                +177286. * CiHe_11 / LambdaNP2
                -552286. * CiHL3_11 / LambdaNP2
                -113484. * CiHD / LambdaNP2
                +29757.9 * CiHB / LambdaNP2
                -69897.4 * CiHW / LambdaNP2
                -385087. * CiHWB / LambdaNP2
                +47999.3 * CiDHB / LambdaNP2
                -98310.4 * CiDHW / LambdaNP2
                -5.76 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +8.942 * deltaMz()
                -0.556 * deltaMh()
                -2.75 * deltaaMZ()
                +5.748 * deltaGmu() );
    
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +120996. * CiHbox / LambdaNP2
                -263143. * CiHL1_11 / LambdaNP2
                +533190. * CiHe_11 / LambdaNP2
                -263143. * CiHL3_11 / LambdaNP2
                +29434.5 * CiHD / LambdaNP2
                +63176.5 * CiHB / LambdaNP2
                +26728.5 * CiHW / LambdaNP2
                +184228. * CiHWB / LambdaNP2
                -85487.1 * CiDHB / LambdaNP2
                -24906.1 * CiDHW / LambdaNP2
                -1.044 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -0.508 * deltaMz()
                -0.545 * deltaMh()
                +1.958 * deltaaMZ()
                +1.027 * deltaGmu() );
            
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121114. * CiHbox / LambdaNP2
                -524119. * CiHL1_11 / LambdaNP2
                +218758. * CiHe_11 / LambdaNP2
                -524119. * CiHL3_11 / LambdaNP2
                -98164. * CiHD / LambdaNP2
                +74694.7 * CiHB / LambdaNP2
                -49060.4 * CiHW / LambdaNP2
                -348619. * CiHWB / LambdaNP2
                +33861.6 * CiDHB / LambdaNP2
                -90369.8 * CiDHW / LambdaNP2
                -5.256 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +7.922 * deltaMz()
                -0.546 * deltaMh()
                -2.261 * deltaaMZ()
                +5.242 * deltaGmu() );
            
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZBFPol()");
        } 
    
    } else if (sqrt_s == 1.4) {
        
        C1 = 0.0065;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +120762. * CiHbox / LambdaNP2
                -242720. * CiHL1_11 / LambdaNP2
                +714345. * CiHe_11 / LambdaNP2
                -242720. * CiHL3_11 / LambdaNP2
                +53823.3 * CiHD / LambdaNP2
                -64876.7 * CiHB / LambdaNP2
                +9362.37 * CiHW / LambdaNP2
                +355440. * CiHWB / LambdaNP2
                -127361. * CiDHB / LambdaNP2
                -18147.3 * CiDHW / LambdaNP2
                -0.228 * DeltaGF()
                ;         

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.05 * deltaMz()
                -0.422 * deltaMh()
                +2.78 * deltaaMZ()
                +0.2 * deltaGmu() );
            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +120818. * CiHbox / LambdaNP2
                -692905. * CiHL1_11 / LambdaNP2
                +184416. * CiHe_11 / LambdaNP2
                -692905. * CiHL3_11 / LambdaNP2
                -121143. * CiHD / LambdaNP2
                -4989.81 * CiHB / LambdaNP2
                -93241.6 * CiHW / LambdaNP2
                -392394. * CiHWB / LambdaNP2
                +60556.9 * CiDHB / LambdaNP2
                -121409. * CiDHW / LambdaNP2
                -6.003 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +9.501 * deltaMz()
                -0.422 * deltaMh()
                -2.999 * deltaaMZ()
                +5.972 * deltaGmu() );
            
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +120773. * CiHbox / LambdaNP2
                -309806. * CiHL1_11 / LambdaNP2
                +643900. * CiHe_11 / LambdaNP2
                -309806. * CiHL3_11 / LambdaNP2
                +29091.1 * CiHD / LambdaNP2
                +22438.3 * CiHB / LambdaNP2
                +16021.7 * CiHW / LambdaNP2
                +202496. * CiHWB / LambdaNP2
                -100775. * CiDHB / LambdaNP2
                -32830.8 * CiDHW / LambdaNP2
                -1.043 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -0.415 * deltaMz()
                -0.422 * deltaMh()
                +1.961 * deltaaMZ()
                +1.014 * deltaGmu() );
            
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +120795. * CiHbox / LambdaNP2
                -637584. * CiHL1_11 / LambdaNP2
                +256188. * CiHe_11 / LambdaNP2
                -637584. * CiHL3_11 / LambdaNP2
                -98543.3 * CiHD / LambdaNP2
                +49040.2 * CiHB / LambdaNP2
                -63051.7 * CiHW / LambdaNP2
                -332850. * CiHWB / LambdaNP2
                +36510.1 * CiDHB / LambdaNP2
                -108018. * CiDHW / LambdaNP2
                -5.256 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +8.01 * deltaMz()
                -0.423 * deltaMh()
                -2.255 * deltaaMZ()
                +5.227 * deltaGmu() );
            
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZBFPol()");
        } 
        
    } else if (sqrt_s == 1.5) {
        
        C1 = 0.0065;// Use the same as 1400 GeV

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +120570. * CiHbox / LambdaNP2
                -250340. * CiHL1_11 / LambdaNP2
                +739684. * CiHe_11 / LambdaNP2
                -250340. * CiHL3_11 / LambdaNP2
                +53685.8 * CiHD / LambdaNP2
                -71192.9 * CiHB / LambdaNP2
                +9743.41 * CiHW / LambdaNP2
                +357556. * CiHWB / LambdaNP2
                -131206. * CiDHB / LambdaNP2
                -19448. * CiDHW / LambdaNP2
                -0.224 * DeltaGF()
                ;          

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.032 * deltaMz()
                -0.4 * deltaMh()
                +2.778 * deltaaMZ()
                +0.194 * deltaGmu() );
            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +120602. * CiHbox / LambdaNP2
                -718001. * CiHL1_11 / LambdaNP2
                +189852. * CiHe_11 / LambdaNP2
                -718001. * CiHL3_11 / LambdaNP2
                -121214. * CiHD / LambdaNP2
                -6057.91 * CiHB / LambdaNP2
                -95148.1 * CiHW / LambdaNP2
                -390958. * CiHWB / LambdaNP2
                +61690.7 * CiDHB / LambdaNP2
                -125382. * CiDHW / LambdaNP2
                -5.997 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +9.519 * deltaMz()
                -0.399 * deltaMh()
                -3.001 * deltaaMZ()
                +5.965 * deltaGmu() );
            
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +120563. * CiHbox / LambdaNP2
                -319378. * CiHL1_11 / LambdaNP2
                +665765. * CiHe_11 / LambdaNP2
                -319378. * CiHL3_11 / LambdaNP2
                +29010.7 * CiHD / LambdaNP2
                +14190.4 * CiHB / LambdaNP2
                +16080. * CiHW / LambdaNP2
                +205187. * CiHWB / LambdaNP2
                -103927. * CiDHB / LambdaNP2
                -34420.2 * CiDHW / LambdaNP2
                -1.04 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -0.398 * deltaMz()
                -0.4 * deltaMh()
                +1.96 * deltaaMZ()
                +1.01 * deltaGmu() );
            
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +120607. * CiHbox / LambdaNP2
                -659879. * CiHL1_11 / LambdaNP2
                +263841. * CiHe_11 / LambdaNP2
                -659879. * CiHL3_11 / LambdaNP2
                -98617.3 * CiHD / LambdaNP2
                +46418.4 * CiHB / LambdaNP2
                -64166.6 * CiHW / LambdaNP2
                -330855. * CiHWB / LambdaNP2
                +36774.5 * CiDHB / LambdaNP2
                -111573. * CiDHW / LambdaNP2
                -5.253 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +8.03 * deltaMz()
                -0.4 * deltaMh()
                -2.257 * deltaaMZ()
                +5.221 * deltaGmu() );
            
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZBFPol()");
        }
    
    } else if (sqrt_s == 3.0) {
        
        C1 = 0.0063;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +120539. * CiHbox / LambdaNP2
                -327096. * CiHL1_11 / LambdaNP2
                +988310. * CiHe_11 / LambdaNP2
                -327096. * CiHL3_11 / LambdaNP2
                +53758.1 * CiHD / LambdaNP2
                -79161. * CiHB / LambdaNP2
                +3856.87 * CiHW / LambdaNP2
                +369878. * CiHWB / LambdaNP2
                -170059. * CiDHB / LambdaNP2
                -32235.8 * CiDHW / LambdaNP2
                -0.226 * DeltaGF()
                ;    

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -1.896 * deltaMz()
                -0.264 * deltaMh()
                +2.778 * deltaaMZ()
                +0.174 * deltaGmu() );
            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +120565. * CiHbox / LambdaNP2
                -961658. * CiHL1_11 / LambdaNP2
                +247947. * CiHe_11 / LambdaNP2
                -961658. * CiHL3_11 / LambdaNP2
                -121230. * CiHD / LambdaNP2
                -10752.9 * CiHB / LambdaNP2
                -92123.7 * CiHW / LambdaNP2
                -391807. * CiHWB / LambdaNP2
                +73242.2 * CiDHB / LambdaNP2
                -165690. * CiDHW / LambdaNP2
                -6.002 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +9.659 * deltaMz()
                -0.264 * deltaMh()
                -3.003 * deltaaMZ()
                +5.943 * deltaGmu() );
            
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +120534. * CiHbox / LambdaNP2
                -417962. * CiHL1_11 / LambdaNP2
                +884851. * CiHe_11 / LambdaNP2
                -417962. * CiHL3_11 / LambdaNP2
                +29065.5 * CiHD / LambdaNP2
                -10885.4 * CiHB / LambdaNP2
                +8249.25 * CiHW / LambdaNP2
                +228820. * CiHWB / LambdaNP2
                -135851. * CiDHB / LambdaNP2
                -51177.2 * CiDHW / LambdaNP2
                -1.04 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -0.262 * deltaMz()
                -0.264 * deltaMh()
                +1.959 * deltaaMZ()
                +0.987 * deltaGmu() );
            
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +120480. * CiHbox / LambdaNP2
                -880604. * CiHL1_11 / LambdaNP2
                +344657. * CiHe_11 / LambdaNP2
                -880604. * CiHL3_11 / LambdaNP2
                -98656.8 * CiHD / LambdaNP2
                +28681.4 * CiHB / LambdaNP2
                -66216.6 * CiHW / LambdaNP2
                -320715. * CiHWB / LambdaNP2
                +41721.6 * CiDHB / LambdaNP2
                -148698. * CiDHW / LambdaNP2
                -5.256 * DeltaGF()
                ; 

    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +8.169 * deltaMz()
                -0.264 * deltaMh()
                -2.259 * deltaaMZ()
                +5.202 * deltaGmu() );
            
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZBFPol()");
        } 
    
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZBFPol()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    //(Assume similar to WBF.)
    mu += eeeWBFint + eeeWBFpar;
    
//  Linear contribution from Higgs self-coupling
    mu = mu + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    mu = mu + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::muepWBF(const double sqrt_s) const
{
    double mu = 1.0;
    
    if (sqrt_s == 1.3) {
        
        mu += 
                +121790. * CiHbox / LambdaNP2
                -161604. * CiHL3_11 / LambdaNP2
                -161282. * CiHQ3_11 / LambdaNP2
                -203141. * CiHD / LambdaNP2
                -88171.6 * CiHW / LambdaNP2
                -377218. * CiHWB / LambdaNP2
                -37738.9 * CiDHW / LambdaNP2
                -4.676 * DeltaGF()
                -4.916 * deltaMwd6()
                ;
        
//        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

//        }
        
    } else if (sqrt_s == 1.8) {
        
        mu += 
                +121867. * CiHbox / LambdaNP2
                -182643. * CiHL3_11 / LambdaNP2
                -181961. * CiHQ3_11 / LambdaNP2
                -202400. * CiHD / LambdaNP2
                -78295.8 * CiHW / LambdaNP2
                -377193. * CiHWB / LambdaNP2
                -45757.3 * CiDHW / LambdaNP2
                -4.672 * DeltaGF()
                -4.637 * deltaMwd6()
                ;
        
//        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

//        }
        
    } else if (sqrt_s == 3.5) {

        mu += 
                +121250. * CiHbox / LambdaNP2
                -216885. * CiHL3_11 / LambdaNP2
                -218544. * CiHQ3_11 / LambdaNP2
                -202390. * CiHD / LambdaNP2
                -64783.2 * CiHW / LambdaNP2
                -377727. * CiHWB / LambdaNP2
                -60431.2 * CiDHW / LambdaNP2
                -4.688 * DeltaGF()
                -4.573 * deltaMwd6()
                ;
        
//        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

//        }
          
    } else if (sqrt_s == 5.0) {

        mu += 
                +119662. * CiHbox / LambdaNP2
                -237868. * CiHL3_11 / LambdaNP2
                -236470. * CiHQ3_11 / LambdaNP2
                -203294. * CiHD / LambdaNP2
                -60911. * CiHW / LambdaNP2
                -378045. * CiHWB / LambdaNP2
                -67483.7 * CiDHW / LambdaNP2
                -4.667 * DeltaGF()
                -4.437 * deltaMwd6()
                ;
        
//        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

//        }
        
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::muepWBF()");
        
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += eepWBFint + eepWBFpar;

    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::muepZBF(const double sqrt_s) const
{
    double mu = 1.0;

    if (sqrt_s == 1.3) {
        
        mu += 
                +121280. * CiHbox / LambdaNP2
                -152367. * CiHL1_11 / LambdaNP2
                +32200. * CiHQ1_11 / LambdaNP2
                +124934. * CiHe_11 / LambdaNP2
                -42209.5 * CiHu_11 / LambdaNP2
                +12445.7 * CiHd_11 / LambdaNP2
                -152367. * CiHL3_11 / LambdaNP2
                -165343. * CiHQ3_11 / LambdaNP2
                -173922. * CiHD / LambdaNP2
                -34636.2 * CiHB / LambdaNP2
                -121438. * CiHW / LambdaNP2
                -74939.1 * CiHWB / LambdaNP2
                -5454.93 * CiDHB / LambdaNP2
                -39349.6 * CiDHW / LambdaNP2
                -3.719 * DeltaGF()
                ;
        
//        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

//        }
        
    } else if (sqrt_s == 1.8) {
        
        mu += 
                +120218. * CiHbox / LambdaNP2
                -173566. * CiHL1_11 / LambdaNP2
                +26307.1 * CiHQ1_11 / LambdaNP2
                +142600. * CiHe_11 / LambdaNP2
                -47449. * CiHu_11 / LambdaNP2
                +14356.2 * CiHd_11 / LambdaNP2
                -173566. * CiHL3_11 / LambdaNP2
                -188606. * CiHQ3_11 / LambdaNP2
                -174301. * CiHD / LambdaNP2
                -19800. * CiHB / LambdaNP2
                -103254. * CiHW / LambdaNP2
                -89049.2 * CiHWB / LambdaNP2
                -8304.85 * CiDHB / LambdaNP2
                -48942.9 * CiDHW / LambdaNP2
                -3.714 * DeltaGF()
                ;
        
//        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

//        }
        
    } else if (sqrt_s == 3.5) {

        mu +=              
                +123119. * CiHbox / LambdaNP2
                -206981. * CiHL1_11 / LambdaNP2
                +18620.9 * CiHQ1_11 / LambdaNP2
                +177706. * CiHe_11 / LambdaNP2
                -53822. * CiHu_11 / LambdaNP2
                +20491.5 * CiHd_11 / LambdaNP2
                -206981. * CiHL3_11 / LambdaNP2
                -227549. * CiHQ3_11 / LambdaNP2
                -172298. * CiHD / LambdaNP2
                -6887.17 * CiHB / LambdaNP2
                -79245. * CiHW / LambdaNP2
                -103223. * CiHWB / LambdaNP2
                -9863.11 * CiDHB / LambdaNP2
                -61304.3 * CiDHW / LambdaNP2
                -3.721 * DeltaGF()
                ;
        
//        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

//        }
          
    } else if (sqrt_s == 5.0) {

        mu += 
                +121709. * CiHbox / LambdaNP2
                -225267. * CiHL1_11 / LambdaNP2
                +13471.8 * CiHQ1_11 / LambdaNP2
                +193542. * CiHe_11 / LambdaNP2
                -57640.9 * CiHu_11 / LambdaNP2
                +22573. * CiHd_11 / LambdaNP2
                -225267. * CiHL3_11 / LambdaNP2
                -247738. * CiHQ3_11 / LambdaNP2
                -172768. * CiHD / LambdaNP2
                -4524.89 * CiHB / LambdaNP2
                -71935.4 * CiHW / LambdaNP2
                -104998. * CiHWB / LambdaNP2
                -11877.8 * CiDHB / LambdaNP2
                -69467.3 * CiDHW / LambdaNP2
                -3.71 * DeltaGF()
                ;
        
//        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

//        }
        
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::muepZBF()");
    
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += eepZBFint + eepZBFpar;

    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::muWH(const double sqrt_s) const
{
    double mu = 1.0;
    
    double C1 = 0.0;
    
    if (sqrt_s == 1.96) {
        
        C1 = 0.0; // N.A.

        mu += 
                +121173. * (1. + eWH_2_Hbox ) * CiHbox / LambdaNP2
                +1566788. * (1. + eWH_2_HQ3_11 ) * CiHQ3_11 / LambdaNP2
                -160914. * (1. + eWH_2_HD ) * CiHD / LambdaNP2
                +860916. * (1. + eWH_2_HW ) * CiHW / LambdaNP2
                -286409. * (1. + eWH_2_HWB ) * CiHWB / LambdaNP2
                +134641. * (1. + eWH_2_DHW ) * CiDHW / LambdaNP2
                -3.31 * (1. + eWH_2_DeltaGF ) * DeltaGF()
                -2.199 * deltaMwd6()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
        
    } else if (sqrt_s == 7.0) {
        
        C1 = 0.0106;

        mu += 
                +121015. * (1. + eWH_78_Hbox ) * CiHbox / LambdaNP2
                +1792020. * (1. + eWH_78_HQ3_11 ) * CiHQ3_11 / LambdaNP2
                -159689. * (1. + eWH_78_HD ) * CiHD / LambdaNP2
                +881065. * (1. + eWH_78_HW ) * CiHW / LambdaNP2
                -283895. * (1. + eWH_78_HWB ) * CiHWB / LambdaNP2
                +168173. * (1. + eWH_78_DHW ) * CiDHW / LambdaNP2
                -3.273 * (1. + eWH_78_DeltaGF ) * DeltaGF()
                -2.143 * deltaMwd6()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
        
    } else if (sqrt_s == 8.0) {
        
        C1 = 0.0105;

        mu += 
                +121226. * (1. + eWH_78_Hbox ) * CiHbox / LambdaNP2
                +1830192. * (1. + eWH_78_HQ3_11 ) * CiHQ3_11 / LambdaNP2
                -159543. * (1. + eWH_78_HD ) * CiHD / LambdaNP2
                +884671. * (1. + eWH_78_HW ) * CiHW / LambdaNP2
                -283662. * (1. + eWH_78_HWB ) * CiHWB / LambdaNP2
                +174061. * (1. + eWH_78_DHW ) * CiDHW / LambdaNP2
                -3.278 * (1. + eWH_78_DeltaGF ) * DeltaGF()
                -2.147 * deltaMwd6()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }

    } else if (sqrt_s == 13.0) {
        
        C1 = 0.0103;

        mu += 
                +120439. * (1. + eWH_1314_Hbox ) * CiHbox / LambdaNP2
                +1953200. * (1. + eWH_1314_HQ3_11 ) * CiHQ3_11 / LambdaNP2
                -159847. * (1. + eWH_1314_HD ) * CiHD / LambdaNP2
                +892264. * (1. + eWH_1314_HW ) * CiHW / LambdaNP2
                -283830. * (1. + eWH_1314_HWB ) * CiHWB / LambdaNP2
                +192168. * (1. + eWH_1314_DHW ) * CiDHW / LambdaNP2
                -3.269 * (1. + eWH_1314_DeltaGF ) * DeltaGF()
                -2.101 * deltaMwd6()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
      
    } else if (sqrt_s == 14.0) {
        
        C1 = 0.0103;

        mu += 
                +120284. * (1. + eWH_1314_Hbox ) * CiHbox / LambdaNP2
                +1971011. * (1. + eWH_1314_HQ3_11 ) * CiHQ3_11 / LambdaNP2
                -159830. * (1. + eWH_1314_HD ) * CiHD / LambdaNP2
                +893216. * (1. + eWH_1314_HW ) * CiHW / LambdaNP2
                -283818. * (1. + eWH_1314_HWB ) * CiHWB / LambdaNP2
                +194877. * (1. + eWH_1314_DHW ) * CiDHW / LambdaNP2
                -3.272 * (1. + eWH_1314_DeltaGF ) * DeltaGF()
                -2.103 * deltaMwd6()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }

    } else if (sqrt_s == 27.0) {
        
        C1 = 0.0101; // From arXiv: 1902.00134

        mu += 
                +120696. * CiHbox / LambdaNP2
                +2105646. * CiHQ3_11 / LambdaNP2
                -159695. * CiHD / LambdaNP2
                +900162. * CiHW / LambdaNP2
                -283257. * CiHWB / LambdaNP2
                +215592. * CiDHW / LambdaNP2
                -3.256 * DeltaGF()
                -2.063 * deltaMwd6()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
    
    } else if (sqrt_s == 100.0) {
        
        C1 = 0.0; // N.A. 

        mu += 
                +121319. * CiHbox / LambdaNP2
                +2294991. * CiHQ3_11 / LambdaNP2
                -159242. * CiHD / LambdaNP2
                +908130. * CiHW / LambdaNP2
                -282574. * CiHWB / LambdaNP2
                +245406. * CiDHW / LambdaNP2
                -3.259 * DeltaGF()
                -2.047 * deltaMwd6()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
          
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::muWH()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += eWHint + eWHpar;
    
//  Linear contribution from Higgs self-coupling
    mu = mu + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    mu = mu + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::muZH(const double sqrt_s) const
{
    double mu = 1.0;
    
    double C1 = 0.0;
    
    if (sqrt_s == 1.96) {
        
        C1 = 0.0; // N.A.

        mu += 
                +121197. * (1. + eZH_2_Hbox ) * CiHbox / LambdaNP2
                -810445. * (1. + eZH_2_HQ1_11 ) * CiHQ1_11 / LambdaNP2
                +529340. * (1. + eZH_2_Hu_11 ) * CiHu_11 / LambdaNP2
                -69410.3 * (1. + eZH_2_Hd_11 ) * CiHd_11 / LambdaNP2
                +1567161. * (1. + eZH_2_HQ3_11 ) * CiHQ3_11 / LambdaNP2
                -16992.5 * (1. + eZH_2_HD ) * CiHD / LambdaNP2
                +79314.5 * (1. + eZH_2_HB ) * CiHB / LambdaNP2
                +711710. * (1. + eZH_2_HW ) * CiHW / LambdaNP2
                +189054. * (1. + eZH_2_HWB ) * CiHWB / LambdaNP2
                +9774.73 * (1. + eZH_2_DHB ) * CiDHB / LambdaNP2
                +130777. * (1. + eZH_2_DHW ) * CiDHW / LambdaNP2
                -2.535 * (1. + eZH_2_DeltaGF ) * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
        
    } else if (sqrt_s == 7.0) {
        
        C1 = 0.0123;

        mu += 
                +121069. * (1. + eZH_78_Hbox ) * CiHbox / LambdaNP2
                -182215. * (1. + eZH_78_HQ1_11 ) * CiHQ1_11 / LambdaNP2
                +421780. * (1. + eZH_78_Hu_11 ) * CiHu_11 / LambdaNP2
                -139169. * (1. + eZH_78_Hd_11 ) * CiHd_11 / LambdaNP2
                +1712111. * (1. + eZH_78_HQ3_11 ) * CiHQ3_11 / LambdaNP2
                -15395.4 * (1. + eZH_78_HD ) * CiHD / LambdaNP2
                +87094.9 * (1. + eZH_78_HB ) * CiHB / LambdaNP2
                +717388. * (1. + eZH_78_HW ) * CiHW / LambdaNP2
                +203105. * (1. + eZH_78_HWB ) * CiHWB / LambdaNP2
                +17532.4 * (1. + eZH_78_DHB ) * CiDHB / LambdaNP2
                +152950. * (1. + eZH_78_DHW ) * CiDHW / LambdaNP2
                -2.502 * (1. + eZH_78_DeltaGF ) * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
        
    } else if (sqrt_s == 8.0) {
        
        C1 = 0.0122;

        mu += 
                +121334. * (1. + eZH_78_Hbox ) * CiHbox / LambdaNP2
                -176804. * (1. + eZH_78_HQ1_11 ) * CiHQ1_11 / LambdaNP2
                +428587. * (1. + eZH_78_Hu_11 ) * CiHu_11 / LambdaNP2
                -142508. * (1. + eZH_78_Hd_11 ) * CiHd_11 / LambdaNP2
                +1747367. * (1. + eZH_78_HQ3_11 ) * CiHQ3_11 / LambdaNP2
                -15002.7 * (1. + eZH_78_HD ) * CiHD / LambdaNP2
                +87781.5 * (1. + eZH_78_HB ) * CiHB / LambdaNP2
                +721405. * (1. + eZH_78_HW ) * CiHW / LambdaNP2
                +204540. * (1. + eZH_78_HWB ) * CiHWB / LambdaNP2
                +18868.6 * (1. + eZH_78_DHB ) * CiDHB / LambdaNP2
                +158432. * (1. + eZH_78_DHW ) * CiDHW / LambdaNP2
                -2.507 * (1. + eZH_78_DeltaGF ) * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }

    } else if (sqrt_s == 13.0) {
        
        C1 = 0.0119;

        mu += 
                +121374. * (1. + eZH_1314_Hbox ) * CiHbox / LambdaNP2
                -152273. * (1. + eZH_1314_HQ1_11 ) * CiHQ1_11 / LambdaNP2
                +448168. * (1. + eZH_1314_Hu_11 ) * CiHu_11 / LambdaNP2
                -155999. * (1. + eZH_1314_Hd_11 ) * CiHd_11 / LambdaNP2
                +1862364. * (1. + eZH_1314_HQ3_11 ) * CiHQ3_11 / LambdaNP2
                -15185. * (1. + eZH_1314_HD ) * CiHD / LambdaNP2
                +88937.9 * (1. + eZH_1314_HB ) * CiHB / LambdaNP2
                +728207. * (1. + eZH_1314_HW ) * CiHW / LambdaNP2
                +207857. * (1. + eZH_1314_HWB ) * CiHWB / LambdaNP2
                +21647.4 * (1. + eZH_1314_DHB ) * CiDHB / LambdaNP2
                +175015. * (1. + eZH_1314_DHW ) * CiDHW / LambdaNP2
                -2.506 * (1. + eZH_1314_DeltaGF ) * DeltaGF()
                ;

        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
        
    } else if (sqrt_s == 14.0) {
        
        C1 = 0.0118;

        mu += 
                +121437. * (1. + eZH_1314_Hbox ) * CiHbox / LambdaNP2
                -147580. * (1. + eZH_1314_HQ1_11 ) * CiHQ1_11 / LambdaNP2
                +450628. * (1. + eZH_1314_Hu_11 ) * CiHu_11 / LambdaNP2
                -157625. * (1. + eZH_1314_Hd_11 ) * CiHd_11 / LambdaNP2
                +1878132. * (1. + eZH_1314_HQ3_11 ) * CiHQ3_11 / LambdaNP2
                -15299.4 * (1. + eZH_1314_HD ) * CiHD / LambdaNP2
                +88761.8 * (1. + eZH_1314_HB ) * CiHB / LambdaNP2
                +729243. * (1. + eZH_1314_HW ) * CiHW / LambdaNP2
                +207707. * (1. + eZH_1314_HWB ) * CiHWB / LambdaNP2
                +21958.9 * (1. + eZH_1314_DHB ) * CiDHB / LambdaNP2
                +177300. * (1. + eZH_1314_DHW ) * CiDHW / LambdaNP2
                -2.507 * (1. + eZH_1314_DeltaGF ) * DeltaGF()
                ;

        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }

    } else if (sqrt_s == 27.0) {
        
        C1 = 0.0116; // From arXiv: 1902.00134

        mu += 
                +121206. * CiHbox / LambdaNP2
                -101865. * CiHQ1_11 / LambdaNP2
                +468029. * CiHu_11 / LambdaNP2
                -173377. * CiHd_11 / LambdaNP2
                +2002478. * CiHQ3_11 / LambdaNP2
                -15486.3 * CiHD / LambdaNP2
                +89958. * CiHB / LambdaNP2
                +735013. * CiHW / LambdaNP2
                +211026. * CiHWB / LambdaNP2
                +25604. * CiDHB / LambdaNP2
                +196710. * CiDHW / LambdaNP2
                -2.505 * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
    
    } else if (sqrt_s == 100.0) {
        
        C1 = 0.0; // N.A.
 
        mu += 
                +121269. * CiHbox / LambdaNP2
                +90.68 * CiHQ1_11 / LambdaNP2
                +484275. * CiHu_11 / LambdaNP2
                -197878. * CiHd_11 / LambdaNP2
                +2175601. * CiHQ3_11 / LambdaNP2
                -14992.4 * CiHD / LambdaNP2
                +91707.3 * CiHB / LambdaNP2
                +741805. * CiHW / LambdaNP2
                +215319. * CiHWB / LambdaNP2
                +31435.6 * CiDHB / LambdaNP2
                +223843. * CiDHW / LambdaNP2
                -2.504 * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::muZH()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += eZHint + eZHpar;
    
//  Linear contribution from Higgs self-coupling
    mu = mu + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    mu = mu + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::mueeZH(const double sqrt_s) const
{
    double mu = 1.0;
    
    double C1 = 0.0;

    if (sqrt_s == 0.240) {
        
        C1 = 0.017;

        mu += 
                +121263. * CiHbox / LambdaNP2
                +898682. * CiHL1_11 / LambdaNP2
                -767820. * CiHe_11 / LambdaNP2
                +898682. * CiHL3_11 / LambdaNP2
                -6046.36 * CiHD / LambdaNP2
                +122439. * CiHB / LambdaNP2
                +540057. * CiHW / LambdaNP2
                +231063. * CiHWB / LambdaNP2
                +17593.2 * CiDHB / LambdaNP2
                +53409.5 * CiDHW / LambdaNP2
                -2.2 * DeltaGF()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( -0.2 * deltaaMZ()
                +2.2 * deltaGmu()
                +4.775 * deltaMz()
                -3.071 * deltaMh() ); 
        
        if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }

    } else if (sqrt_s == 0.250) {
        
        C1 = 0.015;

        mu += 
                +121263. * CiHbox / LambdaNP2
                +975101. * CiHL1_11 / LambdaNP2
                -833750. * CiHe_11 / LambdaNP2
                +975101. * CiHL3_11 / LambdaNP2
                -6046.36 * CiHD / LambdaNP2
                +128443. * CiHB / LambdaNP2
                +568273. * CiHW / LambdaNP2
                +244206. * CiHWB / LambdaNP2
                +19818.6 * CiDHB / LambdaNP2
                +60127.6 * CiDHW / LambdaNP2
                -2.2 * DeltaGF()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( -0.2 * deltaaMZ()
                +2.2 * deltaGmu()
                +5.219 * deltaMz()
                -2.27 * deltaMh() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 0.350) {
        
        C1 = 0.0057;

        mu += 
                +121283. * CiHbox / LambdaNP2
                +1911340. * CiHL1_11 / LambdaNP2
                -1640958. * CiHe_11 / LambdaNP2
                +1911340. * CiHL3_11 / LambdaNP2
                -6009.52 * CiHD / LambdaNP2
                +173183. * CiHB / LambdaNP2
                +785843. * CiHW / LambdaNP2
                +344494. * CiHWB / LambdaNP2
                +59158.7 * CiDHB / LambdaNP2
                +167954. * CiDHW / LambdaNP2
                -2.201 * DeltaGF()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( -0.2 * deltaaMZ()
                +2.2 * deltaGmu()
                +5.396 * deltaMz()
                -0.729 * deltaMh() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 0.365) {
        
        C1 = 0.0057; // Use same as 350 GeV

        mu += 
                +121243. * CiHbox / LambdaNP2
                +2078482. * CiHL1_11 / LambdaNP2
                -1785085. * CiHe_11 / LambdaNP2
                +2078482. * CiHL3_11 / LambdaNP2
                -6010.65 * CiHD / LambdaNP2
                +178173. * CiHB / LambdaNP2
                +809806. * CiHW / LambdaNP2
                +355487. * CiHWB / LambdaNP2
                +67662.7 * CiDHB / LambdaNP2
                +190194. * CiDHW / LambdaNP2
                -2.201 * DeltaGF()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( -0.2 * deltaaMZ()
                +2.2 * deltaGmu()
                +5.348 * deltaMz()
                -0.664 * deltaMh() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 0.380) {
        
        C1 = 0.0057; // Use same as 350 GeV

        mu += 
                +121281. * CiHbox / LambdaNP2
                +2253013. * CiHL1_11 / LambdaNP2
                -1934557. * CiHe_11 / LambdaNP2
                +2253013. * CiHL3_11 / LambdaNP2
                -6026.37 * CiHD / LambdaNP2
                +182674. * CiHB / LambdaNP2
                +832109. * CiHW / LambdaNP2
                +365819. * CiHWB / LambdaNP2
                +76742. * CiDHB / LambdaNP2
                +214030. * CiDHW / LambdaNP2
                -2.202 * DeltaGF()
                ;    
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( -0.2 * deltaaMZ()
                +2.2 * deltaGmu()
                +5.301 * deltaMz()
                -0.609 * deltaMh() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 0.500) {
        
        C1 = 0.00099;

        mu += 
                +121264. * CiHbox / LambdaNP2
                +3900384. * CiHL1_11 / LambdaNP2
                -3350136. * CiHe_11 / LambdaNP2
                +3900384. * CiHL3_11 / LambdaNP2
                -6019.22 * CiHD / LambdaNP2
                +209229. * CiHB / LambdaNP2
                +959942. * CiHW / LambdaNP2
                +425112. * CiHWB / LambdaNP2
                +169841. * CiDHB / LambdaNP2
                +455437. * CiDHW / LambdaNP2
                -2.202 * DeltaGF()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( -0.2 * deltaaMZ()
                +2.2 * deltaGmu()
                +5. * deltaMz()
                -0.351 * deltaMh() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 1.0) {
        
        C1 = -0.0012;

        mu += 
                +121274. * CiHbox / LambdaNP2
                +15601820. * CiHL1_11 / LambdaNP2
                -13395670. * CiHe_11 / LambdaNP2
                +15601820. * CiHL3_11 / LambdaNP2
                -6040.16 * CiHD / LambdaNP2
                +243960. * CiHB / LambdaNP2
                +1128805. * CiHW / LambdaNP2
                +503138. * CiHWB / LambdaNP2
                +899357. * CiDHB / LambdaNP2
                +2321619. * CiDHW / LambdaNP2
                -2.202 * DeltaGF()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( -0.2 * deltaaMZ()
                +2.2 * deltaGmu()
                +4.574 * deltaMz()
                -0.092 * deltaMh() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 1.4) {
        
        C1 = -0.0011;

        mu += 
                +121283. * CiHbox / LambdaNP2
                +30579278. * CiHL1_11 / LambdaNP2
                -26253064. * CiHe_11 / LambdaNP2
                +30579278. * CiHL3_11 / LambdaNP2
                -6010.77 * CiHD / LambdaNP2
                +250804. * CiHB / LambdaNP2
                +1161208. * CiHW / LambdaNP2
                +518040. * CiHWB / LambdaNP2
                +1848758. * CiDHB / LambdaNP2
                +4747422. * CiDHW / LambdaNP2
                -2.203 * DeltaGF()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( -0.2 * deltaaMZ()
                +2.2 * deltaGmu()
                +4.491 * deltaMz()
                -0.047 * deltaMh() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 1.5) {
        
        C1 = -0.0011;// Use the same as 1400 GeV

        mu += 
                +121262. * CiHbox / LambdaNP2
                +35102329. * CiHL1_11 / LambdaNP2
                -30135878. * CiHe_11 / LambdaNP2
                +35102329. * CiHL3_11 / LambdaNP2
                -6034.22 * CiHD / LambdaNP2
                +251576. * CiHB / LambdaNP2
                +1165634. * CiHW / LambdaNP2
                +519954. * CiHWB / LambdaNP2
                +2132554. * CiDHB / LambdaNP2
                +5481906. * CiDHW / LambdaNP2
                -2.203 * DeltaGF()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( -0.2 * deltaaMZ()
                +2.2 * deltaGmu()
                +4.479 * deltaMz()
                -0.041 * deltaMh() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 3.0) {
        
        C1 = -0.00054;
        
        mu += 
                +121279. * CiHbox / LambdaNP2
                +140413697. * CiHL1_11 / LambdaNP2
                -120540988. * CiHe_11 / LambdaNP2
                +140413697. * CiHL3_11 / LambdaNP2
                -6012.61 * CiHD / LambdaNP2
                +257222. * CiHB / LambdaNP2
                +1188444. * CiHW / LambdaNP2
                +530503. * CiHWB / LambdaNP2
                +8839419. * CiDHB / LambdaNP2
                +22583370. * CiDHW / LambdaNP2
                -2.202 * DeltaGF()
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( -0.2 * deltaaMZ()
                +2.2 * deltaGmu()
                +4.42 * deltaMz()
                -0.01 * deltaMh() );
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZH()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += eeeZHint + eeeZHpar;
    
//  Linear contribution from Higgs self-coupling
    mu = mu + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    mu = mu + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::mueeZllH(const double sqrt_s) const
{
    
//  The signal strength eeZH
    double mu = mueeZH(sqrt_s);
    
//  The (relative) linear correction to the Z>ll BR
    double deltaBRratio;

    deltaBRratio = deltaGamma_Zf(leptons[ELECTRON]) 
            + deltaGamma_Zf(leptons[MU]);
    
    deltaBRratio = deltaBRratio / 
            ( trueSM.GammaZ(leptons[ELECTRON]) + trueSM.GammaZ(leptons[MU]) );
    
    deltaBRratio = deltaBRratio - deltaGamma_Z() / trueSM.Gamma_Z();

    return mu + deltaBRratio;
}

double NPSMEFTd6::mueeZqqH(const double sqrt_s) const
{
    
//  The signal strength eeZH
    double mu = mueeZH(sqrt_s);
    
//  The (relative) linear correction to the Z>qq BR
    double deltaBRratio;

    deltaBRratio = deltaGamma_Zf(quarks[UP]) 
            + deltaGamma_Zf(quarks[DOWN])
            + deltaGamma_Zf(quarks[CHARM])
            + deltaGamma_Zf(quarks[STRANGE])
            + deltaGamma_Zf(quarks[BOTTOM]);
    
    deltaBRratio = deltaBRratio / 
            ( trueSM.GammaZ(quarks[UP]) + trueSM.GammaZ(quarks[DOWN])
            + trueSM.GammaZ(quarks[CHARM]) + trueSM.GammaZ(quarks[STRANGE])
            + trueSM.GammaZ(quarks[BOTTOM]));
    
    deltaBRratio = deltaBRratio - deltaGamma_Z() / trueSM.Gamma_Z();

    return mu + deltaBRratio;
}

double NPSMEFTd6::mueeZHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    double mu = 1.0;
    
    double C1 = 0.0;

    if (sqrt_s == 0.240) {
        
        C1 = 0.017;
        
        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121260. * CiHbox / LambdaNP2
                +117191. * CiHL1_11 / LambdaNP2
                -1681596. * CiHe_11 / LambdaNP2
                +117191. * CiHL3_11 / LambdaNP2
                +74555.1 * CiHD / LambdaNP2
                +528105. * CiHB / LambdaNP2
                +134403. * CiHW / LambdaNP2
                +872560. * CiHWB / LambdaNP2
                +137571. * CiDHB / LambdaNP2
                -12321.5 * CiDHW / LambdaNP2
                +0.459 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +2.46 * deltaaMZ()
                -0.46 * deltaGmu()
                -0.544 * deltaMz()
                -3.071 * deltaMh() );
    
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121254. * CiHbox / LambdaNP2
                +1495015. * CiHL1_11 / LambdaNP2
                -76567.2 * CiHe_11 / LambdaNP2
                +1495015. * CiHL3_11 / LambdaNP2
                -67582.1 * CiHD / LambdaNP2
                -187104. * CiHB / LambdaNP2
                +849552. * CiHW / LambdaNP2
                -258537. * CiHWB / LambdaNP2
                -73970.1 * CiDHB / LambdaNP2
                +103582. * CiDHW / LambdaNP2
                -4.23 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.23 * deltaaMZ()
                +4.23 * deltaGmu()
                +8.834 * deltaMz()
                -3.071 * deltaMh() );
    
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121256. * CiHbox / LambdaNP2
                +204529. * CiHL1_11 / LambdaNP2
                -1578998. * CiHe_11 / LambdaNP2
                +204529. * CiHL3_11 / LambdaNP2
                +65548.7 * CiHD / LambdaNP2
                +482729. * CiHB / LambdaNP2
                +179733. * CiHW / LambdaNP2
                +800870. * CiHWB / LambdaNP2
                +124170. * CiDHB / LambdaNP2
                -5016.48 * CiDHW / LambdaNP2
                +0.162 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +2.163 * deltaaMZ()
                -0.163 * deltaGmu()
                +0.05 * deltaMz()
                -3.071 * deltaMh() );
    
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121264. * CiHbox / LambdaNP2
                +1442776. * CiHL1_11 / LambdaNP2
                -137405. * CiHe_11 / LambdaNP2
                +1442776. * CiHL3_11 / LambdaNP2
                -62167.6 * CiHD / LambdaNP2
                -159988. * CiHB / LambdaNP2
                +822448. * CiHW / LambdaNP2
                -215639. * CiHWB / LambdaNP2
                -65950.1 * CiDHB / LambdaNP2
                +99206.1 * CiDHW / LambdaNP2
                -4.052 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.052 * deltaaMZ()
                +4.052 * deltaGmu()
                +8.479 * deltaMz()
                -3.071 * deltaMh() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZHPol()");
        }
        
    } else if (sqrt_s == 0.250) {
        
        C1 = 0.015;
        
        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121264. * CiHbox / LambdaNP2
                +127210. * CiHL1_11 / LambdaNP2
                -1824910. * CiHe_11 / LambdaNP2
                +127210. * CiHL3_11 / LambdaNP2
                +74597.1 * CiHD / LambdaNP2
                +560319. * CiHB / LambdaNP2
                +136129. * CiHW / LambdaNP2
                +902676. * CiHWB / LambdaNP2
                +154358. * CiDHB / LambdaNP2
                -13612.9 * CiDHW / LambdaNP2
                +0.459 * DeltaGF()
                ;      
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +2.46 * deltaaMZ()
                -0.46 * deltaGmu()
                -0.1 * deltaMz()
                -2.27 * deltaMh() );
    
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121257. * CiHbox / LambdaNP2
                +1622228. * CiHL1_11 / LambdaNP2
                -83107. * CiHe_11 / LambdaNP2
                +1622228. * CiHL3_11 / LambdaNP2
                -67554.3 * CiHD / LambdaNP2
                -201409. * CiHB / LambdaNP2
                +898116. * CiHW / LambdaNP2
                -258306. * CiHWB / LambdaNP2
                -82898. * CiDHB / LambdaNP2
                +116421. * CiDHW / LambdaNP2
                -4.23 * DeltaGF()
                ; 
                    
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.23 * deltaaMZ()
                +4.23 * deltaGmu()
                +9.279 * deltaMz()
                -2.27 * deltaMh() );
    
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121309. * CiHbox / LambdaNP2
                +221930. * CiHL1_11 / LambdaNP2
                -1714047. * CiHe_11 / LambdaNP2
                +221930. * CiHL3_11 / LambdaNP2
                +65599.6 * CiHD / LambdaNP2
                +512136. * CiHB / LambdaNP2
                +184424. * CiHW / LambdaNP2
                +829145. * CiHWB / LambdaNP2
                +139369. * CiDHB / LambdaNP2
                -5351.17 * CiDHW / LambdaNP2
                +0.162 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +2.163 * deltaaMZ()
                -0.163 * deltaGmu()
                +0.494 * deltaMz()
                -2.27 * deltaMh() );
    
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121269. * CiHbox / LambdaNP2
                +1565559. * CiHL1_11 / LambdaNP2
                -148908. * CiHe_11 / LambdaNP2
                +1565559. * CiHL3_11 / LambdaNP2
                -62170. * CiHD / LambdaNP2
                -172540. * CiHB / LambdaNP2
                +869218. * CiHW / LambdaNP2
                -214299. * CiHWB / LambdaNP2
                -73929.8 * CiDHB / LambdaNP2
                +111494. * CiDHW / LambdaNP2
                -4.053 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.052 * deltaaMZ()
                +4.052 * deltaGmu()
                +8.923 * deltaMz()
                -2.27 * deltaMh() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZHPol()");
        }
        
    } else if (sqrt_s == 0.350) {
        
        C1 = 0.0057;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121274. * CiHbox / LambdaNP2
                +249309. * CiHL1_11 / LambdaNP2
                -3576996. * CiHe_11 / LambdaNP2
                +249309. * CiHL3_11 / LambdaNP2
                +74596.5 * CiHD / LambdaNP2
                +812491. * CiHB / LambdaNP2
                +146212. * CiHW / LambdaNP2
                +1135161. * CiHWB / LambdaNP2
                +395085. * CiDHB / LambdaNP2
                -16140.8 * CiDHW / LambdaNP2
                +0.458 * DeltaGF()
                ;
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +2.46 * deltaaMZ()
                -0.46 * deltaGmu()
                +0.077 * deltaMz()
                -0.729 * deltaMh() );
                
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121289. * CiHbox / LambdaNP2
                +3179548. * CiHL1_11 / LambdaNP2
                -163347. * CiHe_11 / LambdaNP2
                +3179548. * CiHL3_11 / LambdaNP2
                -67524.8 * CiHD / LambdaNP2
                -314653. * CiHB / LambdaNP2
                +1273817. * CiHW / LambdaNP2
                -258947. * CiHWB / LambdaNP2
                -197137. * CiDHB / LambdaNP2
                +308384. * CiDHW / LambdaNP2
                -4.231 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.23 * deltaaMZ()
                +4.23 * deltaGmu()
                +9.456 * deltaMz()
                -0.729 * deltaMh() );
    
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121304. * CiHbox / LambdaNP2
                +434952. * CiHL1_11 / LambdaNP2
                -3360980. * CiHe_11 / LambdaNP2
                +434952. * CiHL3_11 / LambdaNP2
                +65624.7 * CiHD / LambdaNP2
                +741142. * CiHB / LambdaNP2
                +217654. * CiHW / LambdaNP2
                +1046799. * CiHWB / LambdaNP2
                +357606. * CiDHB / LambdaNP2
                +4440.1 * CiDHW / LambdaNP2
                +0.161 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +2.163 * deltaaMZ()
                -0.163 * deltaGmu()
                +0.671 * deltaMz()
                -0.729 * deltaMh() );
    
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121259. * CiHbox / LambdaNP2
                +3068356. * CiHL1_11 / LambdaNP2
                -292427. * CiHe_11 / LambdaNP2
                +3068356. * CiHL3_11 / LambdaNP2
                -62160.7 * CiHD / LambdaNP2
                -271962. * CiHB / LambdaNP2
                +1231171. * CiHW / LambdaNP2
                -206112. * CiHWB / LambdaNP2
                -174718. * CiDHB / LambdaNP2
                +296046. * CiDHW / LambdaNP2
                -4.053 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.052 * deltaaMZ()
                +4.052 * deltaGmu()
                +9.1 * deltaMz()
                -0.729 * deltaMh() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZHPol()");
        } 
    
    } else if (sqrt_s == 0.365) {
        
        C1 = 0.0057; // Use same as 350 GeV

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121270. * CiHbox / LambdaNP2
                +271098. * CiHL1_11 / LambdaNP2
                -3890169. * CiHe_11 / LambdaNP2
                +271098. * CiHL3_11 / LambdaNP2
                +74554. * CiHD / LambdaNP2
                +840573. * CiHB / LambdaNP2
                +147108. * CiHW / LambdaNP2
                +1160947. * CiHWB / LambdaNP2
                +442125. * CiDHB / LambdaNP2
                -15038.8 * CiDHW / LambdaNP2
                +0.459 * DeltaGF()
                ;     
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +2.46 * deltaaMZ()
                -0.46 * deltaGmu()
                +0.029 * deltaMz()
                -0.664 * deltaMh() );
    
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121238. * CiHbox / LambdaNP2
                +3457848. * CiHL1_11 / LambdaNP2
                -177584. * CiHe_11 / LambdaNP2
                +3457848. * CiHL3_11 / LambdaNP2
                -67578.3 * CiHD / LambdaNP2
                -327391. * CiHB / LambdaNP2
                +1315671. * CiHW / LambdaNP2
                -259142. * CiHWB / LambdaNP2
                -218241. * CiDHB / LambdaNP2
                +346804. * CiDHW / LambdaNP2
                -4.231 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.23 * deltaaMZ()
                +4.23 * deltaGmu()
                +9.408 * deltaMz()
                -0.664 * deltaMh() );
    
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121251. * CiHbox / LambdaNP2
                +472985. * CiHL1_11 / LambdaNP2
                -3655203. * CiHe_11 / LambdaNP2
                +472985. * CiHL3_11 / LambdaNP2
                +65559.4 * CiHD / LambdaNP2
                +766585. * CiHB / LambdaNP2
                +221202. * CiHW / LambdaNP2
                +1070933. * CiHWB / LambdaNP2
                +400293. * CiDHB / LambdaNP2
                +7914.02 * CiDHW / LambdaNP2
                +0.161 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +2.163 * deltaaMZ()
                -0.163 * deltaGmu()
                +0.623 * deltaMz()
                -0.664 * deltaMh() );
    
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121238. * CiHbox / LambdaNP2
                +3336984. * CiHL1_11 / LambdaNP2
                -317944. * CiHe_11 / LambdaNP2
                +3336984. * CiHL3_11 / LambdaNP2
                -62188.9 * CiHD / LambdaNP2
                -283174. * CiHB / LambdaNP2
                +1271272. * CiHW / LambdaNP2
                -205330. * CiHWB / LambdaNP2
                -193153. * CiDHB / LambdaNP2
                +333078. * CiDHW / LambdaNP2
                -4.053 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.052 * deltaaMZ()
                +4.052 * deltaGmu()
                +9.052 * deltaMz()
                -0.664 * deltaMh() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZHPol()");
        } 
    
    } else if (sqrt_s == 0.380) {
        
        C1 = 0.0057; // Use same as 350 GeV

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121228. * CiHbox / LambdaNP2
                +293860. * CiHL1_11 / LambdaNP2
                -4216491. * CiHe_11 / LambdaNP2
                +293860. * CiHL3_11 / LambdaNP2
                +74561.4 * CiHD / LambdaNP2
                +866754. * CiHB / LambdaNP2
                +147982. * CiHW / LambdaNP2
                +1184912. * CiHWB / LambdaNP2
                +492018. * CiDHB / LambdaNP2
                -13596.5 * CiDHW / LambdaNP2
                +0.459 * DeltaGF()
                ;           
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +2.46 * deltaaMZ()
                -0.46 * deltaGmu()
                -0.018 * deltaMz()
                -0.609 * deltaMh() );
    
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121226. * CiHbox / LambdaNP2
                +3747707. * CiHL1_11 / LambdaNP2
                -192650. * CiHe_11 / LambdaNP2
                +3747707. * CiHL3_11 / LambdaNP2
                -67608.3 * CiHD / LambdaNP2
                -339193. * CiHB / LambdaNP2
                +1354040. * CiHW / LambdaNP2
                -259321. * CiHWB / LambdaNP2
                -240311. * CiDHB / LambdaNP2
                +387710. * CiDHW / LambdaNP2
                -4.23 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.23 * deltaaMZ()
                +4.23 * deltaGmu()
                +9.361 * deltaMz()
                -0.609 * deltaMh() );
    
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121325. * CiHbox / LambdaNP2
                +512707. * CiHL1_11 / LambdaNP2
                -3961665. * CiHe_11 / LambdaNP2
                +512707. * CiHL3_11 / LambdaNP2
                +65601.7 * CiHD / LambdaNP2
                +790306. * CiHB / LambdaNP2
                +224394. * CiHW / LambdaNP2
                +1093297. * CiHWB / LambdaNP2
                +445530. * CiDHB / LambdaNP2
                +11860.4 * CiDHW / LambdaNP2
                +0.161 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +2.163 * deltaaMZ()
                -0.163 * deltaGmu()
                +0.576 * deltaMz()
                -0.609 * deltaMh() );
    
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121273. * CiHbox / LambdaNP2
                +3617032. * CiHL1_11 / LambdaNP2
                -344629. * CiHe_11 / LambdaNP2
                +3617032. * CiHL3_11 / LambdaNP2
                -62148.3 * CiHD / LambdaNP2
                -293491. * CiHB / LambdaNP2
                +1308558. * CiHW / LambdaNP2
                -204594. * CiHWB / LambdaNP2
                -212514. * CiDHB / LambdaNP2
                +372554. * CiDHW / LambdaNP2
                -4.053 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.052 * deltaaMZ()
                +4.052 * deltaGmu()
                +9.005 * deltaMz()
                -0.609 * deltaMh() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZHPol()");
        } 
    
    } else if (sqrt_s == 0.500) {
        
        C1 = 0.00099;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121268. * CiHbox / LambdaNP2
                +508715. * CiHL1_11 / LambdaNP2
                -7299333. * CiHe_11 / LambdaNP2
                +508715. * CiHL3_11 / LambdaNP2
                +74603.6 * CiHD / LambdaNP2
                +1018069. * CiHB / LambdaNP2
                +151257. * CiHW / LambdaNP2
                +1323862. * CiHWB / LambdaNP2
                +985604. * CiDHB / LambdaNP2
                +8362.16 * CiDHW / LambdaNP2
                +0.458 * DeltaGF()
                ;    
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +2.46 * deltaaMZ()
                -0.46 * deltaGmu()
                -0.319 * deltaMz()
                -0.351 * deltaMh() );
    
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121273. * CiHbox / LambdaNP2
                +6488707. * CiHL1_11 / LambdaNP2
                -332950. * CiHe_11 / LambdaNP2
                +6488707. * CiHL3_11 / LambdaNP2
                -67530.9 * CiHD / LambdaNP2
                -408101. * CiHB / LambdaNP2
                +1576859. * CiHW / LambdaNP2
                -260777. * CiHWB / LambdaNP2
                -452746. * CiDHB / LambdaNP2
                +796569. * CiDHW / LambdaNP2
                -4.231 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.23 * deltaaMZ()
                +4.23 * deltaGmu()
                +9.06 * deltaMz()
                -0.351 * deltaMh() );
    
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121280. * CiHbox / LambdaNP2
                +887632. * CiHL1_11 / LambdaNP2
                -6858533. * CiHe_11 / LambdaNP2
                +887632. * CiHL3_11 / LambdaNP2
                +65606.6 * CiHD / LambdaNP2
                +927745. * CiHB / LambdaNP2
                +241619. * CiHW / LambdaNP2
                +1223535. * CiHWB / LambdaNP2
                +894441. * CiDHB / LambdaNP2
                +58317. * CiDHW / LambdaNP2
                +0.161 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +2.163 * deltaaMZ()
                -0.163 * deltaGmu()
                +0.275 * deltaMz()
                -0.351 * deltaMh() );
    
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121268. * CiHbox / LambdaNP2
                +6262095. * CiHL1_11 / LambdaNP2
                -597046. * CiHe_11 / LambdaNP2
                +6262095. * CiHL3_11 / LambdaNP2
                -62148.8 * CiHD / LambdaNP2
                -353914. * CiHB / LambdaNP2
                +1522841. * CiHW / LambdaNP2
                -200684. * CiHWB / LambdaNP2
                -398214. * CiDHB / LambdaNP2
                +766821. * CiDHW / LambdaNP2
                -4.054 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.052 * deltaaMZ()
                +4.052 * deltaGmu()
                +8.704 * deltaMz()
                -0.351 * deltaMh() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZHPol()");
        } 
    
    } else if (sqrt_s == 1.0) {
        
        C1 = -0.0012;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121236. * CiHbox / LambdaNP2
                +2034785. * CiHL1_11 / LambdaNP2
                -29195703. * CiHe_11 / LambdaNP2
                +2034785. * CiHL3_11 / LambdaNP2
                +74612.7 * CiHD / LambdaNP2
                +1218284. * CiHB / LambdaNP2
                +154779. * CiHW / LambdaNP2
                +1507673. * CiHWB / LambdaNP2
                +4701988. * CiDHB / LambdaNP2
                +239404. * CiDHW / LambdaNP2
                +0.458 * DeltaGF()
                ;     
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +2.46 * deltaaMZ()
                -0.46 * deltaGmu()
                -0.745 * deltaMz()
                -0.092 * deltaMh() );
    
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121298. * CiHbox / LambdaNP2
                +25954994. * CiHL1_11 / LambdaNP2
                -1333713. * CiHe_11 / LambdaNP2
                +25954994. * CiHL3_11 / LambdaNP2
                -67536.7 * CiHD / LambdaNP2
                -499699. * CiHB / LambdaNP2
                +1872177. * CiHW / LambdaNP2
                -263454. * CiHWB / LambdaNP2
                -1999387. * CiDHB / LambdaNP2
                +3910434. * CiDHW / LambdaNP2
                -4.233 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.23 * deltaaMZ()
                +4.23 * deltaGmu()
                +8.633 * deltaMz()
                -0.092 * deltaMh() );
    
        } else if (Pol_em == 80. && Pol_ep == -20.){
            mu += 
                +121257. * CiHbox / LambdaNP2
                +2475072. * CiHL1_11 / LambdaNP2
                -28682974. * CiHe_11 / LambdaNP2
                +2475072. * CiHL3_11 / LambdaNP2
                +72023. * CiHD / LambdaNP2
                +1186280. * CiHB / LambdaNP2
                +186435. * CiHW / LambdaNP2
                +1475072. * CiHWB / LambdaNP2
                +4578518. * CiDHB / LambdaNP2
                +307070. * CiDHW / LambdaNP2
                +0.371 * DeltaGF()
                ;     
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -0.572 * deltaMz()
                -0.091 * deltaMh()
                +2.375 * deltaaMZ()
                -0.377 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 20.){
            mu += 
                +121306. * CiHbox / LambdaNP2
                +25696973. * CiHL1_11 / LambdaNP2
                -1634825. * CiHe_11 / LambdaNP2
                +25696973. * CiHL3_11 / LambdaNP2
                -65976.8 * CiHD / LambdaNP2
                -480973. * CiHB / LambdaNP2
                +1853631. * CiHW / LambdaNP2
                -244288. * CiHWB / LambdaNP2
                -1927204. * CiDHB / LambdaNP2
                +3870798. * CiDHW / LambdaNP2
                -4.182 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +8.536 * deltaMz()
                -0.09 * deltaMh()
                -2.178 * deltaaMZ()
                +4.178 * deltaGmu() );
    
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121307. * CiHbox / LambdaNP2
                +3550656. * CiHL1_11 / LambdaNP2
                -27432206. * CiHe_11 / LambdaNP2
                +3550656. * CiHL3_11 / LambdaNP2
                +65607.4 * CiHD / LambdaNP2
                +1109435. * CiHB / LambdaNP2
                +263679. * CiHW / LambdaNP2
                +1395519. * CiHWB / LambdaNP2
                +4277336. * CiDHB / LambdaNP2
                +472106. * CiDHW / LambdaNP2
                +0.159 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +2.163 * deltaaMZ()
                -0.163 * deltaGmu()
                -0.151 * deltaMz()
                -0.092 * deltaMh() );
    
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121327. * CiHbox / LambdaNP2
                +25048839. * CiHL1_11 / LambdaNP2
                -2390358. * CiHe_11 / LambdaNP2
                +25048839. * CiHL3_11 / LambdaNP2
                -62132.7 * CiHD / LambdaNP2
                -434824. * CiHB / LambdaNP2
                +1807095. * CiHW / LambdaNP2
                -196264. * CiHWB / LambdaNP2
                -1746222. * CiDHB / LambdaNP2
                +3771341. * CiDHW / LambdaNP2
                -4.056 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.052 * deltaaMZ()
                +4.052 * deltaGmu()
                +8.278 * deltaMz()
                -0.092 * deltaMh() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZHPol()");
        } 
    
    } else if (sqrt_s == 1.4) {
        
        C1 = -0.0011;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121277. * CiHbox / LambdaNP2
                +3988231. * CiHL1_11 / LambdaNP2
                -57226150. * CiHe_11 / LambdaNP2
                +3988231. * CiHL3_11 / LambdaNP2
                +74608.5 * CiHD / LambdaNP2
                +1256970. * CiHB / LambdaNP2
                +155358. * CiHW / LambdaNP2
                +1542655. * CiHWB / LambdaNP2
                +9506894. * CiDHB / LambdaNP2
                +553431. * CiDHW / LambdaNP2
                +0.457 * DeltaGF()
                ;        
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +2.46 * deltaaMZ()
                -0.46 * deltaGmu()
                -0.828 * deltaMz()
                -0.047 * deltaMh() );
    
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121314. * CiHbox / LambdaNP2
                +50871646. * CiHL1_11 / LambdaNP2
                -2614134. * CiHe_11 / LambdaNP2
                +50871646. * CiHL3_11 / LambdaNP2
                -67535.5 * CiHD / LambdaNP2
                -516385. * CiHB / LambdaNP2
                +1928805. * CiHW / LambdaNP2
                -264072. * CiHWB / LambdaNP2
                -3989947. * CiDHB / LambdaNP2
                +7948308. * CiDHW / LambdaNP2
                -4.233 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.23 * deltaaMZ()
                +4.23 * deltaGmu()
                +8.55 * deltaMz()
                -0.047 * deltaMh() );
    
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121250. * CiHbox / LambdaNP2
                +6958750. * CiHL1_11 / LambdaNP2
                -53762500. * CiHe_11 / LambdaNP2
                +6958750. * CiHL3_11 / LambdaNP2
                +65589.3 * CiHD / LambdaNP2
                +1144464. * CiHB / LambdaNP2
                +267732. * CiHW / LambdaNP2
                +1428214. * CiHWB / LambdaNP2
                +8650536. * CiDHB / LambdaNP2
                +1021964. * CiDHW / LambdaNP2
                +0.16 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +2.163 * deltaaMZ()
                -0.163 * deltaGmu()
                -0.234 * deltaMz()
                -0.047 * deltaMh() );
    
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121278. * CiHbox / LambdaNP2
                +49094486. * CiHL1_11 / LambdaNP2
                -4685522. * CiHe_11 / LambdaNP2
                +49094486. * CiHL3_11 / LambdaNP2
                -62150.9 * CiHD / LambdaNP2
                -450090. * CiHB / LambdaNP2
                +1861602. * CiHW / LambdaNP2
                -195621. * CiHWB / LambdaNP2
                -3478338. * CiDHB / LambdaNP2
                +7668095. * CiDHW / LambdaNP2
                -4.055 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.052 * deltaaMZ()
                +4.052 * deltaGmu()
                +8.195 * deltaMz()
                -0.047 * deltaMh() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZHPol()");
        } 
        
    } else if (sqrt_s == 1.5) {
        
        C1 = -0.0011;// Use the same as 1400 GeV

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121268. * CiHbox / LambdaNP2
                +4578315. * CiHL1_11 / LambdaNP2
                -65691823. * CiHe_11 / LambdaNP2
                +4578315. * CiHL3_11 / LambdaNP2
                +74595.2 * CiHD / LambdaNP2
                +1262261. * CiHB / LambdaNP2
                +155435. * CiHW / LambdaNP2
                +1547379. * CiHWB / LambdaNP2
                +10961322. * CiDHB / LambdaNP2
                +649157. * CiDHW / LambdaNP2
                +0.457 * DeltaGF()
                ;          
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +2.46 * deltaaMZ()
                -0.46 * deltaGmu()
                -0.84 * deltaMz()
                -0.041 * deltaMh() );
    
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121277. * CiHbox / LambdaNP2
                +58398883. * CiHL1_11 / LambdaNP2
                -3000385. * CiHe_11 / LambdaNP2
                +58398883. * CiHL3_11 / LambdaNP2
                -67535.8 * CiHD / LambdaNP2
                -518798. * CiHB / LambdaNP2
                +1936613. * CiHW / LambdaNP2
                -264171. * CiHWB / LambdaNP2
                -4590136. * CiDHB / LambdaNP2
                +9169803. * CiDHW / LambdaNP2
                -4.233 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.23 * deltaaMZ()
                +4.23 * deltaGmu()
                +8.539 * deltaMz()
                -0.041 * deltaMh() );
    
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121289. * CiHbox / LambdaNP2
                +7988570. * CiHL1_11 / LambdaNP2
                -61718691. * CiHe_11 / LambdaNP2
                +7988570. * CiHL3_11 / LambdaNP2
                +65599. * CiHD / LambdaNP2
                +1149083. * CiHB / LambdaNP2
                +268317. * CiHW / LambdaNP2
                +1432777. * CiHWB / LambdaNP2
                +9972576. * CiDHB / LambdaNP2
                +1188554. * CiDHW / LambdaNP2
                +0.16 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +2.163 * deltaaMZ()
                -0.163 * deltaGmu()
                -0.246 * deltaMz()
                -0.041 * deltaMh() );
    
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121259. * CiHbox / LambdaNP2
                +56356946. * CiHL1_11 / LambdaNP2
                -5378233. * CiHe_11 / LambdaNP2
                +56356946. * CiHL3_11 / LambdaNP2
                -62168.7 * CiHD / LambdaNP2
                -452149. * CiHB / LambdaNP2
                +1869136. * CiHW / LambdaNP2
                -195562. * CiHWB / LambdaNP2
                -4000306. * CiDHB / LambdaNP2
                +8846432. * CiDHW / LambdaNP2
                -4.055 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.052 * deltaaMZ()
                +4.052 * deltaGmu()
                +8.183 * deltaMz()
                -0.041 * deltaMh() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZHPol()");
        }
    
    } else if (sqrt_s == 3.0) {
        
        C1 = -0.00054;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121320. * CiHbox / LambdaNP2
                +18314161. * CiHL1_11 / LambdaNP2
                -262773345. * CiHe_11 / LambdaNP2
                +18314161. * CiHL3_11 / LambdaNP2
                +74663.6 * CiHD / LambdaNP2
                +1289569. * CiHB / LambdaNP2
                +155612. * CiHW / LambdaNP2
                +1572580. * CiHWB / LambdaNP2
                +44806408. * CiDHB / LambdaNP2
                +2877519. * CiDHW / LambdaNP2
                +0.456 * DeltaGF()
                ;          
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +2.46 * deltaaMZ()
                -0.46 * deltaGmu()
                -0.899 * deltaMz()
                -0.01 * deltaMh() );
    
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121305. * CiHbox / LambdaNP2
                +233598342. * CiHL1_11 / LambdaNP2
                -12002450. * CiHe_11 / LambdaNP2
                +233598342. * CiHL3_11 / LambdaNP2
                -67507.7 * CiHD / LambdaNP2
                -531387. * CiHB / LambdaNP2
                +1976750. * CiHW / LambdaNP2
                -264661. * CiHWB / LambdaNP2
                -18587969. * CiDHB / LambdaNP2
                +37618569. * CiDHW / LambdaNP2
                -4.233 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.23 * deltaaMZ()
                +4.23 * deltaGmu()
                +8.48 * deltaMz()
                -0.01 * deltaMh() );
    
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121225. * CiHbox / LambdaNP2
                +31953446. * CiHL1_11 / LambdaNP2
                -246870182. * CiHe_11 / LambdaNP2
                +31953446. * CiHL3_11 / LambdaNP2
                +65576.5 * CiHD / LambdaNP2
                +1173703. * CiHB / LambdaNP2
                +270983. * CiHW / LambdaNP2
                +1456032. * CiHWB / LambdaNP2
                +40783748. * CiDHB / LambdaNP2
                +5077924. * CiDHW / LambdaNP2
                +0.16 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +2.163 * deltaaMZ()
                -0.163 * deltaGmu()
                -0.305 * deltaMz()
                -0.01 * deltaMh() );
    
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121248. * CiHbox / LambdaNP2
                +225427310. * CiHL1_11 / LambdaNP2
                -21505526. * CiHe_11 / LambdaNP2
                +225427310. * CiHL3_11 / LambdaNP2
                -62193.4 * CiHD / LambdaNP2
                -463403. * CiHB / LambdaNP2
                +1907593. * CiHW / LambdaNP2
                -195017. * CiHWB / LambdaNP2
                -16188019. * CiDHB / LambdaNP2
                +36299719. * CiDHW / LambdaNP2
                -4.054 * DeltaGF()
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -2.052 * deltaaMZ()
                +4.052 * deltaGmu()
                +8.124 * deltaMz()
                -0.01 * deltaMh() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZHPol()");
        } 
    
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZHPol()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += eeeZHint + eeeZHpar;
    
//  Linear contribution from Higgs self-coupling
    mu = mu + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    mu = mu + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::mueeZllHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    
//  The signal strength eeZH
    double mu = mueeZHPol(sqrt_s, Pol_em, Pol_ep);
    
//  The (relative) linear correction to the Z>ll BR
    double deltaBRratio;

    deltaBRratio = deltaGamma_Zf(leptons[ELECTRON]) 
            + deltaGamma_Zf(leptons[MU]);
    
    deltaBRratio = deltaBRratio / 
            ( trueSM.GammaZ(leptons[ELECTRON]) + trueSM.GammaZ(leptons[MU]) );
    
    deltaBRratio = deltaBRratio - deltaGamma_Z() / trueSM.Gamma_Z();

    return mu + deltaBRratio;
}

double NPSMEFTd6::mueeZqqHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    
//  The signal strength eeZH
    double mu = mueeZHPol(sqrt_s, Pol_em, Pol_ep);
    
//  The (relative) linear correction to the Z>qq BR
    double deltaBRratio;

    deltaBRratio = deltaGamma_Zf(quarks[UP]) 
            + deltaGamma_Zf(quarks[DOWN])
            + deltaGamma_Zf(quarks[CHARM])
            + deltaGamma_Zf(quarks[STRANGE])
            + deltaGamma_Zf(quarks[BOTTOM]);
    
    deltaBRratio = deltaBRratio / 
            ( trueSM.GammaZ(quarks[UP]) + trueSM.GammaZ(quarks[DOWN])
            + trueSM.GammaZ(quarks[CHARM]) + trueSM.GammaZ(quarks[STRANGE])
            + trueSM.GammaZ(quarks[BOTTOM]));
    
    deltaBRratio = deltaBRratio - deltaGamma_Z() / trueSM.Gamma_Z();

    return mu + deltaBRratio;
}

double NPSMEFTd6::aPskPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    
    // Expression missing CLL contributions!
    
    double aL, aR, aPol;
    double sM = sqrt_s * sqrt_s;
    double Mz2 = Mz*Mz;
    double MH2 = mHl*mHl;
    double dMz = 0.0;
    double dMH = 0.0;
    double dv,dg,dgp,dgL,dgR;
    double kCM, kCM2, EZ, EZ2, kZ, kH;
    double EtaZ;
    double CHpsk, CTpsk,CHL,CHLp, CHE;
    double CWB, CBB, CWW;
    
    // Convention for dim 6 operators
    CWB = g2_tree*g2_tree/(8.0*g2_tree*g1_tree) * CiHWB * v2_over_LambdaNP2;
    CBB = 0.25 * (g2_tree*g2_tree/g1_tree/g1_tree) * CiHB * v2_over_LambdaNP2;
    CWW = 0.25 * CiHW * v2_over_LambdaNP2;

    CHpsk = ( -2.0 * CiHbox + 0.25 * CiHD ) * v2_over_LambdaNP2;    
    CTpsk = -0.5 * CiHD * v2_over_LambdaNP2;    
    CHL = CiHL1_11 * v2_over_LambdaNP2;     
    CHLp = CiHL3_11 * v2_over_LambdaNP2; 
    CHE = CiHe_11 * v2_over_LambdaNP2;
    
    //  Other parameters (1): Missing CLL!!!
    dv = 0.5 * ( CiHL3_11 + CiHL3_22 )* v2_over_LambdaNP2;
    
    // WFR
    EtaZ = -(1.0/2.0)*CHpsk + 2.0*dMz - dv - CTpsk;
    
    // Kinematics
    kCM = sqrt( (sM*sM + (MH2 - Mz2)*(MH2 - Mz2) - 2.0*sM*(MH2 + Mz2))/(4.0*sM) );   
    kCM2 = kCM*kCM;
    
    EZ = sqrt( Mz2 + kCM2);
    EZ2 = EZ*EZ;
    
    kZ = 2.0*Mz2/(sM - Mz2) + (EZ*Mz2)/(2*kCM2*sqrt_s) - Mz2/(2*kCM2) - (EZ2/Mz2)/(2.0 + EZ2/Mz2)*(1.0 - Mz2/(EZ*sqrt_s));
   
    kH = -((EZ*MH2)/(2*kCM2*sqrt_s)) - (EZ2/Mz2)/(2 + EZ2/Mz2)*MH2/(EZ*sqrt_s); 
    
    //  Other parameters (2): Missing CLL!!!
    dg = -(1.0/(g1_tree * ( cW2_tree*cW2_tree - sW2_tree*sW2_tree))) * ( dv * cW2_tree * g1_tree 
            - cW2_tree * dMz * g1_tree + 0.25 * CiHD * cW2_tree * g1_tree * v2_over_LambdaNP2 
            + CiHW * cW2_tree*cW2_tree * g1_tree * v2_over_LambdaNP2 + CiHWB * cW2_tree * g2_tree * sW2_tree * v2_over_LambdaNP2 
            - CiHW * g1_tree * sW2_tree*sW2_tree * v2_over_LambdaNP2 + CiHWB * g2_tree * sW2_tree*sW2_tree * v2_over_LambdaNP2);
    
    
    dgp = -(1.0/(cW2_tree * g1_tree * g1_tree * (-cW2_tree*cW2_tree + sW2_tree*sW2_tree))) * ( dv * cW2_tree * g1_tree * g1_tree * sW2_tree 
            - cW2_tree * dMz * g1_tree * g1_tree * sW2_tree + 0.25 * CiHD * cW2_tree * g1_tree*g1_tree * sW2_tree * v2_over_LambdaNP2 
            + CiHWB * cW2_tree * cW2_tree * g1_tree * g2_tree * sW2_tree * v2_over_LambdaNP2 
            - CiHB * cW2_tree * cW2_tree * g2_tree * g2_tree * sW2_tree * v2_over_LambdaNP2 
            + CiHWB * cW2_tree * g1_tree * g2_tree * sW2_tree * sW2_tree * v2_over_LambdaNP2 
            + CiHB * g2_tree* g2_tree * sW2_tree*sW2_tree*sW2_tree * v2_over_LambdaNP2 );
    
    dgL = (1.0/(0.5 - sW2_tree))*(cW2_tree*(0.5 + sW2_tree)*dg 
            - sW2_tree*(0.5 + cW2_tree)*dgp 
            + 0.5*(CHL + CHLp) 
            + 0.25*cW2_tree*(1.0 + 2.0*sW2_tree)*8.0*CWW 
            - 0.5*sW2_tree*(1.0 - 2.0*sW2_tree)*8.0*CWB 
            - 0.25*sW2_tree*sW2_tree/cW2_tree*(1.0 + 2.0*cW2_tree)*8.0*CBB);
    
    dgR = -cW2_tree*dg + (1.0 + cW2_tree)*dgp 
            - 1.0/(2.0*sW2_tree)*CHE - 0.5*cW2_tree*8*CWW 
            + cW2_tree*8.0*CWB + 0.5*sW2_tree/cW2_tree*(1.0 + cW2_tree)*8.0*CBB;
    
            
    //  LH and RH pars            
    
    aL = dgL + 2*dMz - dv + EtaZ + (sM - Mz2)/(2*Mz2)*(CHL + CHLp)/(0.5 - sW2_tree) + kZ*dMz + kH*dMH;
    aR = dgR + 2*dMz - dv + EtaZ - (sM - Mz2)/(2*Mz2)*CHE/sW2_tree + kZ*dMz + kH*dMH;
    
    //  Polarized a parameter
    aPol = 0.25 * ( (1.0 - Pol_em/100.0)*(1.0 + Pol_ep/100.0) * aL
            + (1.0 + Pol_em/100.0)*(1.0 - Pol_ep/100.0) * aR );
    
    return aPol;
}

double NPSMEFTd6::bPskPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    double bL, bR, bPol;
    double sM = sqrt_s * sqrt_s;
    double Mz2 = Mz*Mz;
    
    double ZetaZ, ZetaAZ;
    double CWB, CBB, CWW;
    
    // Convention for dim 6 operators
    CWB = g2_tree*g2_tree/(8.0*g2_tree*g1_tree) * CiHWB * v2_over_LambdaNP2;
    CBB = 0.25 * (g2_tree*g2_tree/g1_tree/g1_tree) * CiHB * v2_over_LambdaNP2;
    CWW = 0.25 * CiHW * v2_over_LambdaNP2;
 
    ZetaZ = cW2_tree*8.0*CWW + 2.0*sW2_tree*8*CWB + (sW2_tree*sW2_tree/cW2_tree)*8.0*CBB;
    ZetaAZ = sW_tree*cW_tree*(8.0*CWW - (1.0 - sW2_tree/cW2_tree)*8*CWB - (sW2_tree/cW2_tree)*8.0*CBB);
    
    //  LH and RH pars 
    bL = ZetaZ + (sW_tree*cW_tree)/(0.5 - sW2_tree)*(sM - Mz2)/sM*ZetaAZ;
    bR = ZetaZ - (cW_tree/sW_tree)*(sM - Mz2)/sM*ZetaAZ;
    
    //  Polarized b parameter
    bPol = 0.25 * ( (1.0 - Pol_em/100.0)*(1.0 + Pol_ep/100.0) * bL
            + (1.0 + Pol_em/100.0)*(1.0 - Pol_ep/100.0) * bR );
    
    return bPol;
}

double NPSMEFTd6::muVH(const double sqrt_s) const
{
    double sigmaWH_SM = computeSigmaWH(sqrt_s);
    double sigmaZH_SM = computeSigmaZH(sqrt_s);
    double sigmaWH = muWH(sqrt_s) * sigmaWH_SM;
    double sigmaZH = muZH(sqrt_s) * sigmaZH_SM;
    double mu = ((sigmaWH + sigmaZH) / (sigmaWH_SM + sigmaZH_SM));
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::muVBFpVH(const double sqrt_s) const
{
    double sigmaWH_SM = computeSigmaWH(sqrt_s);
    double sigmaZH_SM = computeSigmaZH(sqrt_s);
    double sigmaVBF_SM = computeSigmaVBF(sqrt_s);
    double sigmaWH = muWH(sqrt_s) * sigmaWH_SM;
    double sigmaZH = muZH(sqrt_s) * sigmaZH_SM;
    double sigmaVBF = muVBF(sqrt_s) * sigmaVBF_SM;
    double mu = ((sigmaWH + sigmaZH + sigmaVBF) / (sigmaWH_SM + sigmaZH_SM + sigmaVBF_SM));
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::muttH(const double sqrt_s) const
{
    double mu = 1.0;
    
     double C1 = 0.0;
    
    if (sqrt_s == 1.96) {
        
        C1 = 0.0; // N.A.

        mu += 
                +423420. * (1. + ettH_2_HG ) * CHG / LambdaNP2
                -4269.4 * (1. + ettH_2_G ) * CG / LambdaNP2
                +566792. * (1. + ettH_2_uG_33r ) * CiuG_33r / LambdaNP2
                -2.854 * (1. + ettH_2_DeltagHt ) * deltaG_hff(quarks[TOP]).real()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
        
    } else if (sqrt_s == 7.0) {
        
        C1 = 0.0387;

        mu += 
                +532200. * (1. + ettH_78_HG ) * CHG / LambdaNP2
                -85145.2 * (1. + ettH_78_G ) * CG / LambdaNP2
                +811678. * (1. + ettH_78_uG_33r ) * CiuG_33r / LambdaNP2
                -2.841 * (1. + ettH_78_DeltagHt ) * deltaG_hff(quarks[TOP]).real()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
        
    } else if (sqrt_s == 8.0) {
        
        C1 = 0.0378;

        mu += 
                +535632. * (1. + ettH_78_HG ) * CHG / LambdaNP2
                -86537.2 * (1. + ettH_78_G ) * CG / LambdaNP2
                +825379. * (1. + ettH_78_uG_33r ) * CiuG_33r / LambdaNP2
                -2.849 * (1. + ettH_78_DeltagHt ) * deltaG_hff(quarks[TOP]).real()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
        
    } else if (sqrt_s == 13.0) {
        
        C1 = 0.0351;

        mu += 
                +538764. * (1. + ettH_1314_HG ) * CHG / LambdaNP2
                -84648. * (1. + ettH_1314_G ) * CG / LambdaNP2
                +860470. * (1. + ettH_1314_uG_33r ) * CiuG_33r / LambdaNP2
                -2.834 * (1. + ettH_1314_DeltagHt ) * deltaG_hff(quarks[TOP]).real()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
        
    } else if (sqrt_s == 14.0) {
        
        C1 = 0.0347;

        mu += 
                +536600. * (1. + ettH_1314_HG ) * CHG / LambdaNP2
                -83824.6 * (1. + ettH_1314_G ) * CG / LambdaNP2
                +863670. * (1. + ettH_1314_uG_33r ) * CiuG_33r / LambdaNP2
                -2.839 * (1. + ettH_1314_DeltagHt ) * deltaG_hff(quarks[TOP]).real()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }

    } else if (sqrt_s == 27.0) {
        
        C1 = 0.0320; // From arXiv: 1902.00134

        mu += 
                +519682. * CHG / LambdaNP2
                -68463.1 * CG / LambdaNP2
                +884060. * CiuG_33r / LambdaNP2
                -2.849 * deltaG_hff(quarks[TOP]).real()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
            
    } else if (sqrt_s == 100.0) {
        
        C1 = 0.0; // N.A.

        mu += 
                +467438. * CHG / LambdaNP2
                -22519. * CG / LambdaNP2
                +880378. * CiuG_33r / LambdaNP2
                -2.837 * deltaG_hff(quarks[TOP]).real()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
        
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::muttH()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += ettHint + ettHpar;
    
//  Linear contribution from Higgs self-coupling
    mu = mu + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    mu = mu + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}


double NPSMEFTd6::mutHq(const double sqrt_s) const
{
    double mu = 1.0;
    
     double C1 = 0.0;
    
    if (sqrt_s == 7.0) {
        
        C1 = 0.0;

        mu +=  0.0;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
        
    } else if (sqrt_s == 8.0) {
        
        C1 = 0.0;

        mu +=  0.0;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
        
    } else if (sqrt_s == 13.0) {
        
        C1 = 0.0;

        mu +=  0.0;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
        
    } else if (sqrt_s == 14.0) {
        
        C1 = 0.0;

        mu +=  0.0;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }

    } else if (sqrt_s == 27.0) {
        
        C1 = 0.0;

        mu +=  0.0;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
            
    } else if (sqrt_s == 100.0) {
        
        C1 = 0.0;

        mu +=  0.0;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
        
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::mutHq()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    //mu += etHqint + etHqpar;
    
//  Linear contribution from Higgs self-coupling
    mu = mu + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    mu = mu + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}


double NPSMEFTd6::muggHpttH(const double sqrt_s) const
{
    double sigmaggH_SM = computeSigmaggH(sqrt_s);
    double sigmattH_SM = computeSigmattH(sqrt_s);
    double sigmaggH = muggH(sqrt_s) * sigmaggH_SM;
    double sigmattH = muttH(sqrt_s) * sigmattH_SM;

    double mu = ((sigmaggH + sigmattH) / (sigmaggH_SM + sigmattH_SM));
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
} 

double NPSMEFTd6::mueettH(const double sqrt_s) const
{
    double mu = 1.0;
    
    double C1 = 0.0;
    
    if (sqrt_s == 0.500) {
        
        C1 = 0.086;

        mu += 
                +121901. * CiHbox / LambdaNP2
                +84038.2 * CiHL1_11 / LambdaNP2
                +41671.2 * CiHe_11 / LambdaNP2
                -31418.2 * CiHu_11 / LambdaNP2
                +84038.2 * CiHL3_11 / LambdaNP2
                -121791. * CiuH_33r / LambdaNP2
                -59467.6 * CiHD / LambdaNP2
                +138929. * CiHB / LambdaNP2
                +130909. * CiHW / LambdaNP2
                -253030. * CiHWB / LambdaNP2
                -1757.66 * CiDHB / LambdaNP2
                +1501.34 * CiDHW / LambdaNP2
                +1386027. * CiuW_33r / LambdaNP2
                +1698012. * CiuB_33r / LambdaNP2
                -1.965 * DeltaGF()
                -1.187 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ;
                    
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( +1.932 * deltaMz()
                -9.827 * deltaMh()
                +1.04 * deltaaMZ()
                +1.992 * deltaGmu()
                -18.476 * deltamt() );
    
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 1.0) {
        
        C1 = 0.017;

        mu += 
                +122013. * CiHbox / LambdaNP2
                +889282. * CiHL1_11 / LambdaNP2
                -543424. * CiHe_11 / LambdaNP2
                -8240.83 * CiHu_11 / LambdaNP2
                +889282. * CiHL3_11 / LambdaNP2
                -116099. * CiuH_33r / LambdaNP2
                -60351.9 * CiHD / LambdaNP2
                +352804. * CiHB / LambdaNP2
                +361918. * CiHW / LambdaNP2
                -397547. * CiHWB / LambdaNP2
                +37326.1 * CiDHB / LambdaNP2
                +113772. * CiDHW / LambdaNP2
                +2758980. * CiuW_33r / LambdaNP2
                +3462941. * CiuB_33r / LambdaNP2
                -2.08 * DeltaGF()
                -2.575 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( +2.185 * deltaMz()
                -1.195 * deltaMh()
                +0.92 * deltaaMZ()
                +2.096 * deltaGmu()
                +2.136 * deltamt() );
    
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 1.4) {
        
        C1 = 0.0094;
       
        mu += 
                +122081. * CiHbox / LambdaNP2
                +2544832. * CiHL1_11 / LambdaNP2
                -1901938. * CiHe_11 / LambdaNP2
                +3241.73 * CiHu_11 / LambdaNP2
                +2544832. * CiHL3_11 / LambdaNP2
                -112208. * CiuH_33r / LambdaNP2
                -60340.4 * CiHD / LambdaNP2
                +464967. * CiHB / LambdaNP2
                +487659. * CiHW / LambdaNP2
                -471053. * CiHWB / LambdaNP2
                +134900. * CiDHB / LambdaNP2
                +371767. * CiDHW / LambdaNP2
                +3804096. * CiuW_33r / LambdaNP2
                +4800265. * CiuB_33r / LambdaNP2
                -2.139 * DeltaGF()
                -3.203 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( +2.309 * deltaMz()
                -0.898 * deltaMh()
                +0.872 * deltaaMZ()
                +2.157 * deltaGmu()
                +2.262 * deltamt() );
    
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 1.5) {
        
        C1 = 0.0094;// Use the same as 1400 GeV

        mu += 
                +122173. * CiHbox / LambdaNP2
                +3117293. * CiHL1_11 / LambdaNP2
                -2378233. * CiHe_11 / LambdaNP2
                +5531.15 * CiHu_11 / LambdaNP2
                +3117293. * CiHL3_11 / LambdaNP2
                -111274. * CiuH_33r / LambdaNP2
                -60192. * CiHD / LambdaNP2
                +487962. * CiHB / LambdaNP2
                +513503. * CiHW / LambdaNP2
                -485782. * CiHWB / LambdaNP2
                +170734. * CiDHB / LambdaNP2
                +462665. * CiDHW / LambdaNP2
                +4068326. * CiuW_33r / LambdaNP2
                +5138930. * CiuB_33r / LambdaNP2
                -2.149 * DeltaGF()
                -3.325 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( +2.322 * deltaMz()
                -0.858 * deltaMh()
                +0.866 * deltaaMZ()
                +2.164 * deltaGmu()
                +2.265 * deltamt() );
    
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 3.0) {
        
        C1 = 0.0037;
        
        mu += 
                +121915. * CiHbox / LambdaNP2
                +19529668. * CiHL1_11 / LambdaNP2
                -16356276. * CiHe_11 / LambdaNP2
                +23142.9 * CiHu_11 / LambdaNP2
                +19529668. * CiHL3_11 / LambdaNP2
                -104011. * CiuH_33r / LambdaNP2
                -58710.4 * CiHD / LambdaNP2
                +697868. * CiHB / LambdaNP2
                +751003. * CiHW / LambdaNP2
                -625171. * CiHWB / LambdaNP2
                +1204441. * CiDHB / LambdaNP2
                +3111413. * CiDHW / LambdaNP2
                +8604912. * CiuW_33r / LambdaNP2
                +10946841. * CiuB_33r / LambdaNP2
                -2.224 * DeltaGF()
                -4.279 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ;
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( +2.483 * deltaMz()
                -0.572 * deltaMh()
                +0.771 * deltaaMZ()
                +2.242 * deltaGmu()
                +2.182 * deltamt() );
    
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::mueettH()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += eeettHint + eeettHpar;
    
//  Linear contribution from Higgs self-coupling
    mu = mu + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    mu = mu + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::mueettHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    double mu = 1.0;
    
    double C1 = 0.0; 
    
    if (sqrt_s == 0.500) {
        
        C1 = 0.086;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121861. * CiHbox / LambdaNP2
                +14207.9 * CiHL1_11 / LambdaNP2
                +124191. * CiHe_11 / LambdaNP2
                +112591. * CiHu_11 / LambdaNP2
                +14207.9 * CiHL3_11 / LambdaNP2
                -123399. * CiuH_33r / LambdaNP2
                -12437.7 * CiHD / LambdaNP2
                +249779. * CiHB / LambdaNP2
                +18912.8 * CiHW / LambdaNP2
                -109936. * CiHWB / LambdaNP2
                -5170.73 * CiDHB / LambdaNP2
                +3167.65 * CiDHW / LambdaNP2
                +174267. * CiuW_33r / LambdaNP2
                +3032981. * CiuB_33r / LambdaNP2
                -0.388 * DeltaGF()
                +3.51 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ;  

            // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -1.319 * deltaMz()
                -9.866 * deltaMh()
                +2.617 * deltaaMZ()
                +0.421 * deltaGmu()
                -18.44 * deltamt() );
        
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121809. * CiHbox / LambdaNP2
                +116253. * CiHL1_11 / LambdaNP2
                +3415.4 * CiHe_11 / LambdaNP2
                -98311.8 * CiHu_11 / LambdaNP2
                +116253. * CiHL3_11 / LambdaNP2
                -121117. * CiuH_33r / LambdaNP2
                -81321.2 * CiHD / LambdaNP2
                +87352.2 * CiHB / LambdaNP2
                +182702. * CiHW / LambdaNP2
                -319427. * CiHWB / LambdaNP2
                -21.616 * CiDHB / LambdaNP2
                +799.81 * CiDHW / LambdaNP2
                +1948272. * CiuW_33r / LambdaNP2
                +1078489. * CiuB_33r / LambdaNP2
                -2.697 * DeltaGF()
                -3.37 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ;

            // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +3.441 * deltaMz()
                -9.806 * deltaMh()
                +0.308 * deltaaMZ()
                +2.725 * deltaGmu()
                -18.491 * deltamt() );
        
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121837. * CiHbox / LambdaNP2
                +24323.6 * CiHL1_11 / LambdaNP2
                +111998. * CiHe_11 / LambdaNP2
                +91391.1 * CiHu_11 / LambdaNP2
                +24323.6 * CiHL3_11 / LambdaNP2
                -123203. * CiuH_33r / LambdaNP2
                -19404.2 * CiHD / LambdaNP2
                +233452. * CiHB / LambdaNP2
                +35310.2 * CiHW / LambdaNP2
                -131019. * CiHWB / LambdaNP2
                -4810.06 * CiDHB / LambdaNP2
                +2842.31 * CiDHW / LambdaNP2
                +351790. * CiuW_33r / LambdaNP2
                +2837005. * CiuB_33r / LambdaNP2
                -0.617 * DeltaGF()
                +2.818 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ; 

            // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -0.843 * deltaMz()
                -9.86 * deltaMh()
                +2.385 * deltaaMZ()
                +0.645 * deltaGmu()
                -18.45 * deltamt() );
        
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121814. * CiHbox / LambdaNP2
                +113858. * CiHL1_11 / LambdaNP2
                +6221.44 * CiHe_11 / LambdaNP2
                -93321.6 * CiHu_11 / LambdaNP2
                +113858. * CiHL3_11 / LambdaNP2
                -121180. * CiuH_33r / LambdaNP2
                -79695. * CiHD / LambdaNP2
                +91201.9 * CiHB / LambdaNP2
                +178853. * CiHW / LambdaNP2
                -314513. * CiHWB / LambdaNP2
                -137.642 * CiDHB / LambdaNP2
                +853.383 * CiDHW / LambdaNP2
                +1906734. * CiuW_33r / LambdaNP2
                +1124181. * CiuB_33r / LambdaNP2
                -2.642 * DeltaGF()
                -3.21 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ; 

            // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +3.33 * deltaMz()
                -9.807 * deltaMh()
                +0.362 * deltaaMZ()
                +2.671 * deltaGmu()
                -18.489 * deltamt() );
        
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueettHPol()");
        }
          
    } else if (sqrt_s == 1.0) {
        
        C1 = 0.017;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +122269. * CiHbox / LambdaNP2
                +148925. * CiHL1_11 / LambdaNP2
                -1516295. * CiHe_11 / LambdaNP2
                +181376. * CiHu_11 / LambdaNP2
                +148925. * CiHL3_11 / LambdaNP2
                -115721. * CiuH_33r / LambdaNP2
                -9966.97 * CiHD / LambdaNP2
                +648027. * CiHB / LambdaNP2
                +58990.6 * CiHW / LambdaNP2
                -166947. * CiHWB / LambdaNP2
                +258446. * CiDHB / LambdaNP2
                +27641. * CiDHW / LambdaNP2
                +416063. * CiuW_33r / LambdaNP2
                +5771745. * CiuB_33r / LambdaNP2
                -0.426 * DeltaGF()
                +3.026 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ;       

            // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -1.159 * deltaMz()
                -1.211 * deltaMh()
                +2.586 * deltaaMZ()
                +0.445 * deltaGmu()
                +2.101 * deltamt() );
        
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +122212. * CiHbox / LambdaNP2
                +1266376. * CiHL1_11 / LambdaNP2
                -47326.8 * CiHe_11 / LambdaNP2
                -104685. * CiHu_11 / LambdaNP2
                +1266376. * CiHL3_11 / LambdaNP2
                -116193. * CiuH_33r / LambdaNP2
                -85861. * CiHD / LambdaNP2
                +202732. * CiHB / LambdaNP2
                +516612. * CiHW / LambdaNP2
                -514723. * CiHWB / LambdaNP2
                -75504.5 * CiDHB / LambdaNP2
                +158356. * CiDHW / LambdaNP2
                +3954267. * CiuW_33r / LambdaNP2
                +2288387. * CiuB_33r / LambdaNP2
                -2.929 * DeltaGF()
                -5.432 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ; 

            // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +3.902 * deltaMz()
                -1.192 * deltaMh()
                +0.075 * deltaaMZ()
                +2.94 * deltaGmu()
                +2.16 * deltamt() );
        
        } else if (Pol_em == 80. && Pol_ep == -20.){
            mu += 
                +122563. * CiHbox / LambdaNP2
                +179718. * CiHL1_11 / LambdaNP2
                -1476392. * CiHe_11 / LambdaNP2
                +173910. * CiHu_11 / LambdaNP2
                +179718. * CiHL3_11 / LambdaNP2
                -115349. * CiuH_33r / LambdaNP2
                -11797.8 * CiHD / LambdaNP2
                +636347. * CiHB / LambdaNP2
                +71703.6 * CiHW / LambdaNP2
                -176417. * CiHWB / LambdaNP2
                +249649. * CiDHB / LambdaNP2
                +31542.3 * CiDHW / LambdaNP2
                +513357. * CiuW_33r / LambdaNP2
                +5678281. * CiuB_33r / LambdaNP2
                -0.497 * DeltaGF()
                +2.823 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ;     
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -0.986 * deltaMz()
                -1.242 * deltaMh()
                +2.514 * deltaaMZ()
                +0.529 * deltaGmu()
                +2.133 * deltamt() );
    
        } else if (Pol_em == -80. && Pol_ep == 20.){
            mu += 
                +122316. * CiHbox / LambdaNP2
                +1258544. * CiHL1_11 / LambdaNP2
                -57807.1 * CiHe_11 / LambdaNP2
                -102560. * CiHu_11 / LambdaNP2
                +1258544. * CiHL3_11 / LambdaNP2
                -116091. * CiuH_33r / LambdaNP2
                -85249.7 * CiHD / LambdaNP2
                +206295. * CiHB / LambdaNP2
                +513404. * CiHW / LambdaNP2
                -512197. * CiHWB / LambdaNP2
                -72925.9 * CiDHB / LambdaNP2
                +157286. * CiDHW / LambdaNP2
                +3929488. * CiuW_33r / LambdaNP2
                +2314064. * CiuB_33r / LambdaNP2
                -2.911 * DeltaGF()
                -5.37 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +3.877 * deltaMz()
                -1.222 * deltaMh()
                +0.099 * deltaaMZ()
                +2.937 * deltaGmu()
                +2.184 * deltamt() );
    
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +122564. * CiHbox / LambdaNP2
                +252265. * CiHL1_11 / LambdaNP2
                -1381101. * CiHe_11 / LambdaNP2
                +155161. * CiHu_11 / LambdaNP2
                +252265. * CiHL3_11 / LambdaNP2
                -115358. * CiuH_33r / LambdaNP2
                -16813.1 * CiHD / LambdaNP2
                +607466. * CiHB / LambdaNP2
                +101359. * CiHW / LambdaNP2
                -198737. * CiHWB / LambdaNP2
                +227834. * CiDHB / LambdaNP2
                +39939.6 * CiDHW / LambdaNP2
                +742520. * CiuW_33r / LambdaNP2
                +5453267. * CiuB_33r / LambdaNP2
                -0.659 * DeltaGF()
                +2.273 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ; 

            // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -0.69 * deltaMz()
                -1.205 * deltaMh()
                +2.349 * deltaaMZ()
                +0.676 * deltaGmu()
                +2.105 * deltamt() );
        
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +122380. * CiHbox / LambdaNP2
                +1238124. * CiHL1_11 / LambdaNP2
                -84811.2 * CiHe_11 / LambdaNP2
                -97259.2 * CiHu_11 / LambdaNP2
                +1238124. * CiHL3_11 / LambdaNP2
                -116044. * CiuH_33r / LambdaNP2
                -83798.9 * CiHD / LambdaNP2
                +214128. * CiHB / LambdaNP2
                +505118. * CiHW / LambdaNP2
                -505830. * CiHWB / LambdaNP2
                -66814.1 * CiDHB / LambdaNP2
                +155075. * CiDHW / LambdaNP2
                +3863710. * CiuW_33r / LambdaNP2
                +2378351. * CiuB_33r / LambdaNP2
                -2.867 * DeltaGF()
                -5.212 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ; 

            // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +3.771 * deltaMz()
                -1.195 * deltaMh()
                +0.137 * deltaaMZ()
                +2.878 * deltaGmu()
                +2.166 * deltamt() );
        
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueettHPol()");
        } 
        
    } else if (sqrt_s == 1.4) {
        
        C1 = 0.0094;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121945. * CiHbox / LambdaNP2
                +416437. * CiHL1_11 / LambdaNP2
                -5198451. * CiHe_11 / LambdaNP2
                +211446. * CiHu_11 / LambdaNP2
                +416437. * CiHL3_11 / LambdaNP2
                -110413. * CiuH_33r / LambdaNP2
                -8089.5 * CiHD / LambdaNP2
                +852065. * CiHB / LambdaNP2
                +78915.7 * CiHW / LambdaNP2
                -191411. * CiHWB / LambdaNP2
                +881670. * CiDHB / LambdaNP2
                +72289.2 * CiDHW / LambdaNP2
                +588296. * CiuW_33r / LambdaNP2
                +7812392. * CiuB_33r / LambdaNP2
                -0.441 * DeltaGF()
                +2.819 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ;            

            // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -1.109 * deltaMz()
                -0.905 * deltaMh()
                +2.571 * deltaaMZ()
                +0.451 * deltaGmu()
                +2.225 * deltamt() );
        
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +122124. * CiHbox / LambdaNP2
                +3668482. * CiHL1_11 / LambdaNP2
                -164738. * CiHe_11 / LambdaNP2
                -106285. * CiHu_11 / LambdaNP2
                +3668482. * CiHL3_11 / LambdaNP2
                -112775. * CiuH_33r / LambdaNP2
                -87497.2 * CiHD / LambdaNP2
                +261266. * CiHB / LambdaNP2
                +703789. * CiHW / LambdaNP2
                -618584. * CiHWB / LambdaNP2
                -257636. * CiDHB / LambdaNP2
                +530202. * CiDHW / LambdaNP2
                +5501929. * CiuW_33r / LambdaNP2
                +3213842. * CiuB_33r / LambdaNP2
                -3.038 * DeltaGF()
                -6.378 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ; 

            // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.12 * deltaMz()
                -0.898 * deltaMh()
                -0.029 * deltaaMZ()
                +3.056 * deltaGmu()
                +2.28 * deltamt() );
        
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121843. * CiHbox / LambdaNP2
                +706068. * CiHL1_11 / LambdaNP2
                -4748505. * CiHe_11 / LambdaNP2
                +182964. * CiHu_11 / LambdaNP2
                +706068. * CiHL3_11 / LambdaNP2
                -110672. * CiuH_33r / LambdaNP2
                -15249.5 * CiHD / LambdaNP2
                +798771. * CiHB / LambdaNP2
                +134415. * CiHW / LambdaNP2
                -229663. * CiHWB / LambdaNP2
                +779863. * CiDHB / LambdaNP2
                +112951. * CiDHW / LambdaNP2
                +1026697. * CiuW_33r / LambdaNP2
                +7402171. * CiuB_33r / LambdaNP2
                -0.673 * DeltaGF()
                +1.996 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ; 

            // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -0.648 * deltaMz()
                -0.901 * deltaMh()
                +2.34 * deltaaMZ()
                +0.693 * deltaGmu()
                +2.232 * deltamt() );
        
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +122069. * CiHbox / LambdaNP2
                +3581543. * CiHL1_11 / LambdaNP2
                -298692. * CiHe_11 / LambdaNP2
                -97874.3 * CiHu_11 / LambdaNP2
                +3581543. * CiHL3_11 / LambdaNP2
                -112737. * CiuH_33r / LambdaNP2
                -85431.2 * CiHD / LambdaNP2
                +276629. * CiHB / LambdaNP2
                +687136. * CiHW / LambdaNP2
                -607155. * CiHWB / LambdaNP2
                -227375. * CiDHB / LambdaNP2
                +517945. * CiDHW / LambdaNP2
                +5370183. * CiuW_33r / LambdaNP2
                +3335906. * CiuB_33r / LambdaNP2
                -2.969 * DeltaGF()
                -6.138 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ; 

            // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +3.976 * deltaMz()
                -0.895 * deltaMh()
                +0.039 * deltaaMZ()
                +2.986 * deltaGmu()
                +2.271 * deltamt() );
        
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueettHPol()");
        } 
        
    } else if (sqrt_s == 1.5) {
        
        C1 = 0.0094;// Use the same as 1400 GeV

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121854. * CiHbox / LambdaNP2
                +507190. * CiHL1_11 / LambdaNP2
                -6475118. * CiHe_11 / LambdaNP2
                +216935. * CiHu_11 / LambdaNP2
                +507190. * CiHL3_11 / LambdaNP2
                -109820. * CiuH_33r / LambdaNP2
                -7568.59 * CiHD / LambdaNP2
                +893094. * CiHB / LambdaNP2
                +82781.5 * CiHW / LambdaNP2
                -196556. * CiHWB / LambdaNP2
                +1099527. * CiDHB / LambdaNP2
                +87228. * CiDHW / LambdaNP2
                +630747. * CiuW_33r / LambdaNP2
                +8328477. * CiuB_33r / LambdaNP2
                -0.442 * DeltaGF()
                +2.756 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ;            

            // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -1.104 * deltaMz()
                -0.856 * deltaMh()
                +2.568 * deltaaMZ()
                +0.455 * deltaGmu()
                +2.232 * deltamt() );
        
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121994. * CiHbox / LambdaNP2
                +4501280. * CiHL1_11 / LambdaNP2
                -206085. * CiHe_11 / LambdaNP2
                -106381. * CiHu_11 / LambdaNP2
                +4501280. * CiHL3_11 / LambdaNP2
                -112104. * CiuH_33r / LambdaNP2
                -87805.6 * CiHD / LambdaNP2
                +273106. * CiHB / LambdaNP2
                +741955. * CiHW / LambdaNP2
                -639545. * CiHWB / LambdaNP2
                -322155. * CiDHB / LambdaNP2
                +661931. * CiDHW / LambdaNP2
                +5892414. * CiuW_33r / LambdaNP2
                +3448015. * CiuB_33r / LambdaNP2
                -3.057 * DeltaGF()
                -6.552 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ; 

            // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.154 * deltaMz()
                -0.856 * deltaMh()
                -0.045 * deltaaMZ()
                +3.071 * deltaGmu()
                +2.287 * deltamt() );
        
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121793. * CiHbox / LambdaNP2
                +861242. * CiHL1_11 / LambdaNP2
                -5919951. * CiHe_11 / LambdaNP2
                +188249. * CiHu_11 / LambdaNP2
                +861242. * CiHL3_11 / LambdaNP2
                -109696. * CiuH_33r / LambdaNP2
                -14806.7 * CiHD / LambdaNP2
                +837632. * CiHB / LambdaNP2
                +141142. * CiHW / LambdaNP2
                -235907. * CiHWB / LambdaNP2
                +973107. * CiDHB / LambdaNP2
                +138331. * CiDHW / LambdaNP2
                +1097452. * CiuW_33r / LambdaNP2
                +7895510. * CiuB_33r / LambdaNP2
                -0.673 * DeltaGF()
                +1.935 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ; 

            // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -0.637 * deltaMz()
                -0.859 * deltaMh()
                +2.339 * deltaaMZ()
                +0.68 * deltaGmu()
                +2.236 * deltamt() );
        
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +122029. * CiHbox / LambdaNP2
                +4394189. * CiHL1_11 / LambdaNP2
                -373205. * CiHe_11 / LambdaNP2
                -97750.6 * CiHu_11 / LambdaNP2
                +4394189. * CiHL3_11 / LambdaNP2
                -112024. * CiuH_33r / LambdaNP2
                -85643.3 * CiHD / LambdaNP2
                +289620. * CiHB / LambdaNP2
                +724463. * CiHW / LambdaNP2
                -627885. * CiHWB / LambdaNP2
                -284076. * CiDHB / LambdaNP2
                +646658. * CiDHW / LambdaNP2
                +5753330. * CiuW_33r / LambdaNP2
                +3578793. * CiuB_33r / LambdaNP2
                -2.989 * DeltaGF()
                -6.311 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ; 

            // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.014 * deltaMz()
                -0.855 * deltaMh()
                +0.024 * deltaaMZ()
                +3.011 * deltaGmu()
                +2.286 * deltamt() );
        
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueettHPol()");
        }
             
    } else if (sqrt_s == 3.0) {
        
        C1 = 0.0037;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +122442. * CiHbox / LambdaNP2
                +3092340. * CiHL1_11 / LambdaNP2
                -43264264. * CiHe_11 / LambdaNP2
                +259622. * CiHu_11 / LambdaNP2
                +3092340. * CiHL3_11 / LambdaNP2
                -100510. * CiuH_33r / LambdaNP2
                -3230.01 * CiHD / LambdaNP2
                +1267548. * CiHB / LambdaNP2
                +118886. * CiHW / LambdaNP2
                -247164. * CiHWB / LambdaNP2
                +7397753. * CiDHB / LambdaNP2
                +510206. * CiDHW / LambdaNP2
                +1343630. * CiuW_33r / LambdaNP2
                +17234081. * CiuB_33r / LambdaNP2
                -0.459 * DeltaGF()
                +2.453 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ;            

            // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -1.07 * deltaMz()
                -0.576 * deltaMh()
                +2.542 * deltaaMZ()
                +0.468 * deltaGmu()
                +2.145 * deltamt() );
        
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +122230. * CiHbox / LambdaNP2
                +28686134. * CiHL1_11 / LambdaNP2
                -1435177. * CiHe_11 / LambdaNP2
                -108195. * CiHu_11 / LambdaNP2
                +28686134. * CiHL3_11 / LambdaNP2
                -105858. * CiuH_33r / LambdaNP2
                -89803.1 * CiHD / LambdaNP2
                +381886. * CiHB / LambdaNP2
                +1102843. * CiHW / LambdaNP2
                -834821. * CiHWB / LambdaNP2
                -2237555. * CiDHB / LambdaNP2
                +4557030. * CiDHW / LambdaNP2
                +12639913. * CiuW_33r / LambdaNP2
                +7455995. * CiuB_33r / LambdaNP2
                -3.212 * DeltaGF()
                -8.009 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ; 

            // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.469 * deltaMz()
                -0.595 * deltaMh()
                -0.222 * deltaaMZ()
                +3.22 * deltaGmu()
                +2.195 * deltamt() );
        
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +122688. * CiHbox / LambdaNP2
                +5271741. * CiHL1_11 / LambdaNP2
                -39707692. * CiHe_11 / LambdaNP2
                +228729. * CiHu_11 / LambdaNP2
                +5271741. * CiHL3_11 / LambdaNP2
                -100891. * CiuH_33r / LambdaNP2
                -10526.3 * CiHD / LambdaNP2
                +1192421. * CiHB / LambdaNP2
                +202915. * CiHW / LambdaNP2
                -296939. * CiHWB / LambdaNP2
                +6582510. * CiDHB / LambdaNP2
                +853895. * CiDHW / LambdaNP2
                +2303644. * CiuW_33r / LambdaNP2
                +16407287. * CiuB_33r / LambdaNP2
                -0.693 * DeltaGF()
                +1.565 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ; 

            // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( -0.597 * deltaMz()
                -0.565 * deltaMh()
                +2.305 * deltaaMZ()
                +0.708 * deltaGmu()
                +2.153 * deltamt() );
        
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121781. * CiHbox / LambdaNP2
                +27966374. * CiHL1_11 / LambdaNP2
                -2597153. * CiHe_11 / LambdaNP2
                -98089.4 * CiHu_11 / LambdaNP2
                +27966374. * CiHL3_11 / LambdaNP2
                -105885. * CiuH_33r / LambdaNP2
                -87600.3 * CiHD / LambdaNP2
                +406305. * CiHB / LambdaNP2
                +1075086. * CiHW / LambdaNP2
                -818808. * CiHWB / LambdaNP2
                -1967062. * CiDHB / LambdaNP2
                +4442109. * CiDHW / LambdaNP2
                +12322125. * CiuW_33r / LambdaNP2
                +7728315. * CiuB_33r / LambdaNP2
                -3.134 * DeltaGF()
                -7.724 * 0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2
                ; 

            // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.305 * deltaMz()
                -0.59 * deltaMh()
                -0.147 * deltaaMZ()
                +3.144 * deltaGmu()
                +2.192 * deltamt() );
        
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueettHPol()");
        } 
    
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::mueettHPol()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += eeettHint + eeettHpar;
    
//  Linear contribution from Higgs self-coupling
    mu = mu + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    mu = mu + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::mummH(const double sqrt_s) const
{
    double mu = 1.0;
    
    double dymu = deltaG_hff(leptons[MU]).real();
    double ymuSM = -(leptons[MU].getMass()) / v();
    
//  The ratio at all energies is given by a scaling of the muon Yukawa.
    mu = 1.0 + 2.0 * dymu/ymuSM ;
        
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu += dymu*dymu/ymuSM/ymuSM;
    }

    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::BrHggRatio() const
{
    double Br = 1.0;
    double dGHiR1=0.0, dGHiR2=0.0, GHiR=1.0;
    
    dGHiR1= deltaGammaHggRatio1();
    
    Br += dGHiR1 - dGammaHTotR1;
    
    if (FlagQuadraticTerms) {
        
        dGHiR2= deltaGammaHggRatio2();
        
        //Add contributions that are quadratic in the effective coefficients
        Br += - dGHiR1 * dGammaHTotR1
                + dGHiR2 - dGammaHTotR2
                + pow(dGammaHTotR1,2.0);            
        }
    
    GHiR += dGHiR1 + dGHiR2;    
    if ((Br < 0) || (GHiR < 0) || (GammaHTotR < 0)) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPSMEFTd6::BrHWWRatio() const
{
    
    return BrHWW4fRatio();

}

double NPSMEFTd6::BrHWlvRatio() const
{
    double Br = 1.0;
    double dGHiR1=0.0, dGHiR2=0.0, GHiR=1.0;
    
    dGHiR1= deltaGammaHWlvRatio1();
    
    Br += dGHiR1 - dGammaHTotR1;
    
    if (FlagQuadraticTerms) {
        
        dGHiR2= deltaGammaHWlvRatio2();
        
        //Add contributions that are quadratic in the effective coefficients
        Br += - dGHiR1 * dGammaHTotR1
                + dGHiR2 - dGammaHTotR2
                + pow(dGammaHTotR1,2.0);            
        }
    
    GHiR += dGHiR1 + dGHiR2;        
    if ((Br < 0) || (GHiR < 0) || (GammaHTotR < 0)) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;
}
    
double NPSMEFTd6::BrHWW2l2vRatio() const
{
    double Br = 1.0;
    double dGHiR1=0.0, dGHiR2=0.0, GHiR=1.0;
    
    dGHiR1= deltaGammaHWW2l2vRatio1();
    
    Br += dGHiR1 - dGammaHTotR1;
    
    if (FlagQuadraticTerms) {
        
        dGHiR2= deltaGammaHWW2l2vRatio2();
        
        //Add contributions that are quadratic in the effective coefficients
        Br += - dGHiR1 * dGammaHTotR1
                + dGHiR2 - dGammaHTotR2
                + pow(dGammaHTotR1,2.0);            
        }
    
    GHiR += dGHiR1 + dGHiR2;        
    if ((Br < 0) || (GHiR < 0) || (GammaHTotR < 0)) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;
}

double NPSMEFTd6::BrHWjjRatio() const
{
    double Br = 1.0;
    double dGHiR1=0.0, dGHiR2=0.0, GHiR=1.0;
    
    dGHiR1= deltaGammaHWjjRatio1();
    
    Br += dGHiR1 - dGammaHTotR1;
    
    if (FlagQuadraticTerms) {
        
        dGHiR2= deltaGammaHWjjRatio2();
        
        //Add contributions that are quadratic in the effective coefficients
        Br += - dGHiR1 * dGammaHTotR1
                + dGHiR2 - dGammaHTotR2
                + pow(dGammaHTotR1,2.0);            
        }
    
    GHiR += dGHiR1 + dGHiR2;        
    if ((Br < 0) || (GHiR < 0) || (GammaHTotR < 0)) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;
}

double NPSMEFTd6::BrHWW4jRatio() const
{
    double Br = 1.0;
    double dGHiR1=0.0, dGHiR2=0.0, GHiR=1.0;
    
    dGHiR1= deltaGammaHWW4jRatio1();
    
    Br += dGHiR1 - dGammaHTotR1;
    
    if (FlagQuadraticTerms) {
        
        dGHiR2= deltaGammaHWW4jRatio2();
        
        //Add contributions that are quadratic in the effective coefficients
        Br += - dGHiR1 * dGammaHTotR1
                + dGHiR2 - dGammaHTotR2
                + pow(dGammaHTotR1,2.0);            
        }
    
    GHiR += dGHiR1 + dGHiR2;        
    if ((Br < 0) || (GHiR < 0) || (GammaHTotR < 0)) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;
}

double NPSMEFTd6::BrHWffRatio() const
{
    double Br = 1.0;
    double dGHiR1=0.0, dGHiR2=0.0, GHiR=1.0;
    
    dGHiR1= deltaGammaHWffRatio1();
    
    Br += dGHiR1 - dGammaHTotR1;
    
    if (FlagQuadraticTerms) {
        
        dGHiR2= deltaGammaHWffRatio2();
        
        //Add contributions that are quadratic in the effective coefficients
        Br += - dGHiR1 * dGammaHTotR1
                + dGHiR2 - dGammaHTotR2
                + pow(dGammaHTotR1,2.0);            
        }
    
    GHiR += dGHiR1 + dGHiR2;        
    if ((Br < 0) || (GHiR < 0) || (GammaHTotR < 0)) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;
}


double NPSMEFTd6::BrHWW4fRatio() const
{    
    double Br = 1.0;
    double dGHiR1=0.0, dGHiR2=0.0, GHiR=1.0;
    
    dGHiR1= deltaGammaHWW4fRatio1();
    
    Br += dGHiR1 - dGammaHTotR1;
    
    if (FlagQuadraticTerms) {
        
        dGHiR2= deltaGammaHWW4fRatio2();
        
        //Add contributions that are quadratic in the effective coefficients
        Br += - dGHiR1 * dGammaHTotR1
                + dGHiR2 - dGammaHTotR2
                + pow(dGammaHTotR1,2.0);            
        }
    
    GHiR += dGHiR1 + dGHiR2;        
    if ((Br < 0) || (GHiR < 0) || (GammaHTotR < 0)) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;
}

double NPSMEFTd6::BrHZZRatio() const
{    
    return BrHZZ4fRatio();
}

double NPSMEFTd6::BrHZllRatio() const
{
    double Br = 1.0;
    double dGHiR1=0.0, dGHiR2=0.0, GHiR=1.0;
    
    dGHiR1= deltaGammaHZllRatio1();
    
    Br += dGHiR1 - dGammaHTotR1;
    
    if (FlagQuadraticTerms) {
        
        dGHiR2= deltaGammaHZllRatio2();
        
        //Add contributions that are quadratic in the effective coefficients
        Br += - dGHiR1 * dGammaHTotR1
                + dGHiR2 - dGammaHTotR2
                + pow(dGammaHTotR1,2.0);            
        }
    
    GHiR += dGHiR1 + dGHiR2;        
    if ((Br < 0) || (GHiR < 0) || (GammaHTotR < 0)) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;
}

double NPSMEFTd6::BrHZZ4lRatio() const
{
    double Br = 1.0;
    double dGHiR1=0.0, dGHiR2=0.0, GHiR=1.0;
    
    dGHiR1= deltaGammaHZZ4lRatio1();
    
    Br += dGHiR1 - dGammaHTotR1;
    
    if (FlagQuadraticTerms) {
        
        dGHiR2= deltaGammaHZZ4lRatio2();
        
        //Add contributions that are quadratic in the effective coefficients
        Br += - dGHiR1 * dGammaHTotR1
                + dGHiR2 - dGammaHTotR2
                + pow(dGammaHTotR1,2.0);            
        }
    
    GHiR += dGHiR1 + dGHiR2;        
    if ((Br < 0) || (GHiR < 0) || (GammaHTotR < 0)) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;
}

double NPSMEFTd6::BrHZZ4eRatio() const
{
    double Br = 1.0;
    double dGHiR1=0.0, dGHiR2=0.0, GHiR=1.0;
    
    dGHiR1= deltaGammaHZZ4eRatio1();
    
    Br += dGHiR1 - dGammaHTotR1;
    
    if (FlagQuadraticTerms) {
        
        dGHiR2= deltaGammaHZZ4eRatio2();
        
        //Add contributions that are quadratic in the effective coefficients
        Br += - dGHiR1 * dGammaHTotR1
                + dGHiR2 - dGammaHTotR2
                + pow(dGammaHTotR1,2.0);            
        }
    
    GHiR += dGHiR1 + dGHiR2;        
    if ((Br < 0) || (GHiR < 0) || (GammaHTotR < 0)) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;
}

double NPSMEFTd6::BrHZZ2e2muRatio() const
{
    double Br = 1.0;
    double dGHiR1=0.0, dGHiR2=0.0, GHiR=1.0;
    
    dGHiR1= deltaGammaHZZ2e2muRatio1();
    
    Br += dGHiR1 - dGammaHTotR1;
    
    if (FlagQuadraticTerms) {
        
        dGHiR2= deltaGammaHZZ2e2muRatio2();
        
        //Add contributions that are quadratic in the effective coefficients
        Br += - dGHiR1 * dGammaHTotR1
                + dGHiR2 - dGammaHTotR2
                + pow(dGammaHTotR1,2.0);            
        }
    
    GHiR += dGHiR1 + dGHiR2;        
    if ((Br < 0) || (GHiR < 0) || (GammaHTotR < 0)) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;
}

double NPSMEFTd6::BrHZZ4muRatio() const
{
    double Br = 1.0;
    double dGHiR1=0.0, dGHiR2=0.0, GHiR=1.0;
    
    dGHiR1= deltaGammaHZZ4muRatio1();
    
    Br += dGHiR1 - dGammaHTotR1;
    
    if (FlagQuadraticTerms) {
        
        dGHiR2= deltaGammaHZZ4muRatio2();
        
        //Add contributions that are quadratic in the effective coefficients
        Br += - dGHiR1 * dGammaHTotR1
                + dGHiR2 - dGammaHTotR2
                + pow(dGammaHTotR1,2.0);            
        }
    
    GHiR += dGHiR1 + dGHiR2;        
    if ((Br < 0) || (GHiR < 0) || (GammaHTotR < 0)) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;
}

double NPSMEFTd6::BrHZvvRatio() const
{
    double Br = 1.0;
    double dGHiR1=0.0, dGHiR2=0.0, GHiR=1.0;
    
    dGHiR1= deltaGammaHZvvRatio1();
    
    Br += dGHiR1 - dGammaHTotR1;
    
    if (FlagQuadraticTerms) {
        
        dGHiR2= deltaGammaHZvvRatio2();
        
        //Add contributions that are quadratic in the effective coefficients
        Br += - dGHiR1 * dGammaHTotR1
                + dGHiR2 - dGammaHTotR2
                + pow(dGammaHTotR1,2.0);            
        }
    
    GHiR += dGHiR1 + dGHiR2;        
    if ((Br < 0) || (GHiR < 0) || (GammaHTotR < 0)) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;
}

double NPSMEFTd6::BrHZZ4vRatio() const
{
    double Br = 1.0;
    double dGHiR1=0.0, dGHiR2=0.0, GHiR=1.0;
    
    dGHiR1= deltaGammaHZZ4vRatio1();
    
    Br += dGHiR1 - dGammaHTotR1;
    
    if (FlagQuadraticTerms) {
        
        dGHiR2= deltaGammaHZZ4vRatio2();
        
        //Add contributions that are quadratic in the effective coefficients
        Br += - dGHiR1 * dGammaHTotR1
                + dGHiR2 - dGammaHTotR2
                + pow(dGammaHTotR1,2.0);            
        }
    
    GHiR += dGHiR1 + dGHiR2;        
    if ((Br < 0) || (GHiR < 0) || (GammaHTotR < 0)) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;
}

double NPSMEFTd6::BrHZuuRatio() const
{
    double Br = 1.0;
    double dGHiR1=0.0, dGHiR2=0.0, GHiR=1.0;
    
    dGHiR1= deltaGammaHZuuRatio1();
    
    Br += dGHiR1 - dGammaHTotR1;
    
    if (FlagQuadraticTerms) {
        
        dGHiR2= deltaGammaHZuuRatio2();
        
        //Add contributions that are quadratic in the effective coefficients
        Br += - dGHiR1 * dGammaHTotR1
                + dGHiR2 - dGammaHTotR2
                + pow(dGammaHTotR1,2.0);            
        }
    
    GHiR += dGHiR1 + dGHiR2;        
    if ((Br < 0) || (GHiR < 0) || (GammaHTotR < 0)) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;
}

double NPSMEFTd6::BrHZZ4uRatio() const
{
    double deltaBRratio;
    
    deltaBRratio = deltaGamma_Zf(quarks[UP]) 
            + deltaGamma_Zf(quarks[CHARM]);
    
    deltaBRratio = deltaBRratio / 
            ( trueSM.GammaZ(quarks[UP]) + trueSM.GammaZ(quarks[CHARM]) );
    
    deltaBRratio = deltaBRratio - deltaGamma_Z() / trueSM.Gamma_Z();
    
    return ( BrHZuuRatio() + deltaBRratio );
}

double NPSMEFTd6::BrHZddRatio() const
{
    double Br = 1.0;
    double dGHiR1=0.0, dGHiR2=0.0, GHiR=1.0;
    
    dGHiR1= deltaGammaHZddRatio1();
    
    Br += dGHiR1 - dGammaHTotR1;
    
    if (FlagQuadraticTerms) {
        
        dGHiR2= deltaGammaHZddRatio2();
        
        //Add contributions that are quadratic in the effective coefficients
        Br += - dGHiR1 * dGammaHTotR1
                + dGHiR2 - dGammaHTotR2
                + pow(dGammaHTotR1,2.0);            
        }
    
    GHiR += dGHiR1 + dGHiR2;        
    if ((Br < 0) || (GHiR < 0) || (GammaHTotR < 0)) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;
}

double NPSMEFTd6::BrHZZ4dRatio() const
{
    double deltaBRratio;
    
    deltaBRratio = deltaGamma_Zf(quarks[DOWN]) 
            + deltaGamma_Zf(quarks[STRANGE]) 
            + deltaGamma_Zf(quarks[BOTTOM]);
    
    deltaBRratio = deltaBRratio / 
            ( trueSM.GammaZ(quarks[DOWN]) 
            + trueSM.GammaZ(quarks[STRANGE]) 
            + trueSM.GammaZ(quarks[BOTTOM]) );
    
    deltaBRratio = deltaBRratio - deltaGamma_Z() / trueSM.Gamma_Z();
    
    return ( BrHZddRatio() + deltaBRratio );
}

double NPSMEFTd6::BrHZffRatio() const
{
    double Br = 1.0;
    double dGHiR1=0.0, dGHiR2=0.0, GHiR=1.0;
    
    dGHiR1= deltaGammaHZffRatio1();
    
    Br += dGHiR1 - dGammaHTotR1;
    
    if (FlagQuadraticTerms) {
        
        dGHiR2= deltaGammaHZffRatio2();
        
        //Add contributions that are quadratic in the effective coefficients
        Br += - dGHiR1 * dGammaHTotR1
                + dGHiR2 - dGammaHTotR2
                + pow(dGammaHTotR1,2.0);            
        }
    
    GHiR += dGHiR1 + dGHiR2;        
    if ((Br < 0) || (GHiR < 0) || (GammaHTotR < 0)) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;
}

double NPSMEFTd6::BrHZZ4fRatio() const
{    
    double Br = 1.0;
    double dGHiR1=0.0, dGHiR2=0.0, GHiR=1.0;
    
    dGHiR1= deltaGammaHZZ4fRatio1();
    
    Br += dGHiR1 - dGammaHTotR1;
    
    if (FlagQuadraticTerms) {
        
        dGHiR2= deltaGammaHZZ4fRatio2();
        
        //Add contributions that are quadratic in the effective coefficients
        Br += - dGHiR1 * dGammaHTotR1
                + dGHiR2 - dGammaHTotR2
                + pow(dGammaHTotR1,2.0);            
        }
    
    GHiR += dGHiR1 + dGHiR2;        
    if ((Br < 0) || (GHiR < 0) || (GammaHTotR < 0)) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;
}

double NPSMEFTd6::BrHZgaRatio() const
{
    double Br = 1.0;
    double dGHiR1=0.0, dGHiR2=0.0, GHiR=1.0;
    
    dGHiR1= deltaGammaHZgaRatio1();
    
    Br += dGHiR1 - dGammaHTotR1;
    
    if (FlagQuadraticTerms) {
        
        dGHiR2= deltaGammaHZgaRatio2();
        
        //Add contributions that are quadratic in the effective coefficients
        Br += - dGHiR1 * dGammaHTotR1
                + dGHiR2 - dGammaHTotR2
                + pow(dGammaHTotR1,2.0);            
        }
    
    GHiR += dGHiR1 + dGHiR2;        
    if ((Br < 0) || (GHiR < 0) || (GammaHTotR < 0)) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPSMEFTd6::BrHZgallRatio() const
{
    double deltaBRratio;
    
    deltaBRratio = deltaGamma_Zf(leptons[ELECTRON]) 
            + deltaGamma_Zf(leptons[MU]);
    
    deltaBRratio = deltaBRratio / 
            ( trueSM.GammaZ(leptons[ELECTRON]) + trueSM.GammaZ(leptons[MU]) );
    
    deltaBRratio = deltaBRratio - deltaGamma_Z() / trueSM.Gamma_Z();
    
    return ( BrHZgaRatio() + deltaBRratio );
}

double NPSMEFTd6::BrHZgaeeRatio() const
{
    double deltaBRratio;
    
    deltaBRratio = deltaGamma_Zf(leptons[ELECTRON]) / (trueSM.GammaZ(leptons[ELECTRON]));
    
    deltaBRratio = deltaBRratio - deltaGamma_Z() / trueSM.Gamma_Z();
    
    return ( BrHZgaRatio() + deltaBRratio );
}

double NPSMEFTd6::BrHZgamumuRatio() const
{
    double deltaBRratio;
    
    deltaBRratio = deltaGamma_Zf(leptons[MU])/(trueSM.GammaZ(leptons[MU]));
    
    deltaBRratio = deltaBRratio - deltaGamma_Z() / trueSM.Gamma_Z();
    
    return ( BrHZgaRatio() + deltaBRratio );
}

double NPSMEFTd6::BrHgagaRatio() const
{
    double Br = 1.0;
    double dGHiR1=0.0, dGHiR2=0.0, GHiR=1.0;
    
    dGHiR1= deltaGammaHgagaRatio1();
    
    Br += dGHiR1 - dGammaHTotR1;
    
    if (FlagQuadraticTerms) {
        
        dGHiR2= deltaGammaHgagaRatio2();
        
        //Add contributions that are quadratic in the effective coefficients
        Br += - dGHiR1 * dGammaHTotR1
                + dGHiR2 - dGammaHTotR2
                + pow(dGammaHTotR1,2.0);            
        }
    
    GHiR += dGHiR1 + dGHiR2;        
    if ((Br < 0) || (GHiR < 0) || (GammaHTotR < 0)) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}
      
double NPSMEFTd6::BrHmumuRatio() const
{
    double Br = 1.0;
    double dGHiR1=0.0, dGHiR2=0.0, GHiR=1.0;
    
    dGHiR1= deltaGammaHmumuRatio1();
    
    Br += dGHiR1 - dGammaHTotR1;
    
    if (FlagQuadraticTerms) {
        
        dGHiR2= deltaGammaHmumuRatio2();
        
        //Add contributions that are quadratic in the effective coefficients
        Br += - dGHiR1 * dGammaHTotR1
                + dGHiR2 - dGammaHTotR2
                + pow(dGammaHTotR1,2.0);            
        }
    
    GHiR += dGHiR1 + dGHiR2;        
    if ((Br < 0) || (GHiR < 0) || (GammaHTotR < 0)) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPSMEFTd6::BrHtautauRatio() const
{
    double Br = 1.0;
    double dGHiR1=0.0, dGHiR2=0.0, GHiR=1.0;
    
    dGHiR1= deltaGammaHtautauRatio1();
    
    Br += dGHiR1 - dGammaHTotR1;
    
    if (FlagQuadraticTerms) {
        
        dGHiR2= deltaGammaHtautauRatio2();
        
        //Add contributions that are quadratic in the effective coefficients
        Br += - dGHiR1 * dGammaHTotR1
                + dGHiR2 - dGammaHTotR2
                + pow(dGammaHTotR1,2.0);            
        }
    
    GHiR += dGHiR1 + dGHiR2;        
    if ((Br < 0) || (GHiR < 0) || (GammaHTotR < 0)) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPSMEFTd6::BrHccRatio() const
{
    double Br = 1.0;
    double dGHiR1=0.0, dGHiR2=0.0, GHiR=1.0;
    
    dGHiR1= deltaGammaHccRatio1();
    
    Br += dGHiR1 - dGammaHTotR1;
    
    if (FlagQuadraticTerms) {
        
        dGHiR2= deltaGammaHccRatio2();
        
        //Add contributions that are quadratic in the effective coefficients
        Br += - dGHiR1 * dGammaHTotR1
                + dGHiR2 - dGammaHTotR2
                + pow(dGammaHTotR1,2.0);            
        }
    
    GHiR += dGHiR1 + dGHiR2;        
    if ((Br < 0) || (GHiR < 0) || (GammaHTotR < 0)) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPSMEFTd6::BrHbbRatio() const
{
    double Br = 1.0;
    double dGHiR1=0.0, dGHiR2=0.0, GHiR=1.0;
    
    dGHiR1= deltaGammaHbbRatio1();
    
    Br += dGHiR1 - dGammaHTotR1;
    
    if (FlagQuadraticTerms) {
        
        dGHiR2= deltaGammaHbbRatio2();
        
        //Add contributions that are quadratic in the effective coefficients
        Br += - dGHiR1 * dGammaHTotR1
                + dGHiR2 - dGammaHTotR2
                + pow(dGammaHTotR1,2.0);            
        }
    
    GHiR += dGHiR1 + dGHiR2;        
    if ((Br < 0) || (GHiR < 0) || (GammaHTotR < 0)) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPSMEFTd6::computeGammaTotalRatio() const
{
    double width = 1.0;

    width += dGammaHTotR1;
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += dGammaHTotR2;
        }
    
    if (width < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return width;

}

double NPSMEFTd6::deltaGammaTotalRatio1() const
{
    double deltaGammaRatio;
    
//  The change in the ratio asumming only SM decays
    deltaGammaRatio = ( trueSM.computeBrHtogg() * deltaGammaHggRatio1()
            + trueSM.computeBrHtoWW() * deltaGammaHWWRatio1()
            + trueSM.computeBrHtoZZ() * deltaGammaHZZRatio1()
            + trueSM.computeBrHtoZga() * deltaGammaHZgaRatio1()
            + trueSM.computeBrHtogaga() * deltaGammaHgagaRatio1()
            + trueSM.computeBrHtomumu() * deltaGammaHmumuRatio1()
            + trueSM.computeBrHtotautau() * deltaGammaHtautauRatio1()
            + trueSM.computeBrHtocc() * deltaGammaHccRatio1()
            + trueSM.computeBrHtobb() * deltaGammaHbbRatio1() );
    
//  Add the effect of the invisible and exotic BR. Include also here the
//  pure contribution from BrHinv and BrHexo even in case of no dim 6 contibutions    
    deltaGammaRatio = -1.0 + (1.0 + deltaGammaRatio) / (1.0 - BrHinv - BrHexo);
    
    return deltaGammaRatio;
}

double NPSMEFTd6::deltaGammaTotalRatio1noError() const
{
    double deltaGammaRatio;
    
//  The change in the ratio asumming only SM decays
    deltaGammaRatio = ( trueSM.computeBrHtogg() * (deltaGammaHggRatio1() - eHggint - eHggpar )
            + trueSM.computeBrHtoWW() * (deltaGammaHWWRatio1() - eHWWint - eHWWpar )
            + trueSM.computeBrHtoZZ() * (deltaGammaHZZRatio1() - eHZZint - eHZZpar )
            + trueSM.computeBrHtoZga() * (deltaGammaHZgaRatio1() - eHZgaint - eHZgapar )
            + trueSM.computeBrHtogaga() * (deltaGammaHgagaRatio1() - eHgagaint - eHgagapar )
            + trueSM.computeBrHtomumu() * (deltaGammaHmumuRatio1() - eHmumuint - eHmumupar )
            + trueSM.computeBrHtotautau() * (deltaGammaHtautauRatio1() - eHtautauint - eHtautaupar )
            + trueSM.computeBrHtocc() * (deltaGammaHccRatio1() - eHccint - eHccpar )
            + trueSM.computeBrHtobb() * (deltaGammaHbbRatio1() - eHbbint - eHbbpar ) );
    
//  Add the effect of the invisible and exotic BR. Include also here the
//  pure contribution from BrHinv and BrHexo even in case of no dim 6 contibutions    
    deltaGammaRatio = -1.0 + (1.0 + deltaGammaRatio) / (1.0 - BrHinv - BrHexo);
    
    return deltaGammaRatio;
}

double NPSMEFTd6::deltaGammaTotalRatio2() const
{
    double deltaGammaRatio;
    
//  The change in the ratio asumming only SM decays
    deltaGammaRatio = trueSM.computeBrHtogg() * deltaGammaHggRatio2()
            + trueSM.computeBrHtoWW() * deltaGammaHWWRatio2()
            + trueSM.computeBrHtoZZ() * deltaGammaHZZRatio2()
            + trueSM.computeBrHtoZga() * deltaGammaHZgaRatio2()
            + trueSM.computeBrHtogaga() * deltaGammaHgagaRatio2()
            + trueSM.computeBrHtomumu() * deltaGammaHmumuRatio2()
            + trueSM.computeBrHtotautau() * deltaGammaHtautauRatio2()
            + trueSM.computeBrHtocc() * deltaGammaHccRatio2()
            + trueSM.computeBrHtobb() * deltaGammaHbbRatio2();
    
//  Add the effect of the invisible and exotic BR and return     
    return (deltaGammaRatio / (1.0 - BrHinv - BrHexo));
}

double NPSMEFTd6::GammaHggRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;

    width += deltaGammaHggRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHggRatio2();  
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHggRatio1() const
{
    double dwidth = 0.0;
    
    double C1 = 0.0066;
    
    dwidth = ( +37526258. * CHG / LambdaNP2
            + cLHd6 * (
            +121248. * CiHbox / LambdaNP2
            +173353. * CiuH_22r / LambdaNP2
            -129155. * CiuH_33r / LambdaNP2
            +248530. * CidH_33r / LambdaNP2
            -30312.1 * CiHD / LambdaNP2
            -60624.1 * DeltaGF() / v() / v() ) 
            );
    
//  Linear contribution from Higgs self-coupling
    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // Add modifications due to small variations of the SM parameters    
    dwidth += cHSM * ( +1.003 * deltaGmu()
                +2.31 * deltaaSMZ()
                +3.276 * deltaMh()
                -0.134 * deltamt()
                -0.106 * deltamb()
                -0.03 * deltamc() );    

    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHggint + eHggpar;
    
    return dwidth;
}

double NPSMEFTd6::deltaGammaHggRatio2() const
{
    double dwidth = 0.0;
     

    //Contributions that are quadratic in the effective coefficients
    return ( dwidth );
    
}

double NPSMEFTd6::GammaHWWRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;
    
    width += deltaGammaHWWRatio1();
      
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHWWRatio2(); 
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHWWRatio1() const
{    
    double dwidth = 0.0;
    
//    double C1 = 0.0073;

    dwidth = deltaGammaHWW4fRatio1();
    
//  Linear contribution from Higgs self-coupling
//    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio(); 
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
//    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
//    dwidth += eHWWint + eHWWpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHWWRatio2() const
{    
    double dwidth = 0.0;
    
    //Contributions that are quadratic in the effective coefficients
    dwidth = deltaGammaHWW4fRatio2();
    
    
    return dwidth;
    
}

double NPSMEFTd6::GammaHWlvRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;

    width += deltaGammaHWlvRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHWlvRatio2(); 
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHWlvRatio1() const
{    
    double dwidth = 0.0;
    
    double C1 = 0.0073;

    dwidth = ( +121875. * CiHbox / LambdaNP2
                +18351.9 * (1.0/2.0) * ( CiHL3_11 + CiHL3_22 ) / LambdaNP2
                -159873. * CiHD / LambdaNP2
                -91288.7 * CiHW / LambdaNP2
                -284689. * CiHWB / LambdaNP2
                +37703.7 * CiDHW / LambdaNP2
                -3.292 * DeltaGF()
                -15.14 * deltaMwd6()  );
    
//  Linear contribution from Higgs self-coupling
    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // Add modifications due to small variations of the SM parameters    
    //dwidth += cHSM * ( 0.0 );
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    //dwidth += eHWWint + eHWWpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHWlvRatio2() const
{
    double dwidth = 0.0;
     

    //Contributions that are quadratic in the effective coefficients
    return ( dwidth );
    
}

double NPSMEFTd6::GammaHWW2l2vRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;
    
    width += deltaGammaHWW2l2vRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHWW2l2vRatio2();
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHWW2l2vRatio1() const
{    
    double dwidth = 0.0;
    
    double C1 = 0.0073;
        
    dwidth = ( +120742. * CiHbox / LambdaNP2
                +131582. * (1.0/2.0) * ( CiHL3_11 + CiHL3_22 ) / LambdaNP2
                -204043. * CiHD / LambdaNP2
                -91463.9 * CiHW / LambdaNP2
                -379529. * CiHWB / LambdaNP2
                +36848.2 * CiDHW / LambdaNP2
                -4.705 * DeltaGF()
                -13.735 * deltaMwd6()
                -0.965 * deltaGwd6()
                  );
    
//  Linear contribution from Higgs self-coupling
    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // Add modifications due to small variations of the SM parameters    
    dwidth += cHSM * ( -12.123 * deltaMz()
                +13.615 * deltaMh()
                +1.756 * deltaaMZ()
                +0.216 * deltaGmu() );
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHWWint + eHWWpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHWW2l2vRatio2() const
{
    double dwidth = 0.0;
    

    //Contributions that are quadratic in the effective coefficients
    return ( dwidth );
    
}

double NPSMEFTd6::GammaHWjjRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;

    width += deltaGammaHWjjRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHWjjRatio2(); 
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHWjjRatio1() const
{    
    double dwidth = 0.0;
    
    double C1 = 0.0073;

    dwidth = ( +121611. * CiHbox / LambdaNP2
                +17701.4 * (1.0/2.0) * ( CiHQ3_11 + CiHQ3_22 ) / LambdaNP2
                -159273. * CiHD / LambdaNP2
                -91021.6 * CiHW / LambdaNP2
                -282574. * CiHWB / LambdaNP2
                +37917.5 * CiDHW / LambdaNP2
                -3.259 * DeltaGF()
                -15.198 * deltaMwd6()  );
    
//  Linear contribution from Higgs self-coupling
    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // Add modifications due to small variations of the SM parameters    
    //dwidth += cHSM * ( 0.0 ); 
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    //dwidth += eHWWint + eHWWpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHWjjRatio2() const
{
    double dwidth = 0.0;
     

    //Contributions that are quadratic in the effective coefficients
    return ( dwidth );
    
}

double NPSMEFTd6::GammaHWW4jRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;

    width += deltaGammaHWW4jRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHWW4jRatio2(); 
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHWW4jRatio1() const
{    
    double dwidth = 0.0;
    
    double C1 = 0.0073;

    dwidth = ( +121936. * CiHbox / LambdaNP2
                +138860. * (1.0/2.0) * ( CiHQ3_11 + CiHQ3_22 ) / LambdaNP2
                -205023. * CiHD / LambdaNP2
                -89938.5 * CiHW / LambdaNP2
                -383944. * CiHWB / LambdaNP2
                +38244.6 * CiDHW / LambdaNP2
                -4.816 * DeltaGF()
                -13.647 * deltaMwd6()
                -0.959 * deltaGwd6()  );
    
//  Linear contribution from Higgs self-coupling
    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // Add modifications due to small variations of the SM parameters    
    dwidth += cHSM * ( -12.168 * deltaMz()
                +13.66 * deltaMh()
                +1.899 * deltaaMZ()
                +0.189 * deltaGmu() ); 
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHWWint + eHWWpar;

    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHWW4jRatio2() const
{
    double dwidth = 0.0;
     

    //Contributions that are quadratic in the effective coefficients
    return ( dwidth );
    
}

double NPSMEFTd6::GammaHWffRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;
    
    width += deltaGammaHWffRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHWffRatio2();
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHWffRatio1() const
{    
    double dwidth = 0.0;
    
    double C1 = 0.0073;

    dwidth = ( +121288. * CiHbox / LambdaNP2
                +5395.21 * (1.0/3.0) * ( CiHL3_11 + CiHL3_22 + CiHL3_33 ) / LambdaNP2
                +11680.9 * (1.0/2.0) * ( CiHQ3_11 + CiHQ3_22 ) / LambdaNP2
                -159787. * CiHD / LambdaNP2
                -91509.1 * CiHW / LambdaNP2
                -283092. * CiHWB / LambdaNP2
                +37845.1 * CiDHW / LambdaNP2
                -3.259 * DeltaGF()
                -15.196 * deltaMwd6()  );
    
//  Linear contribution from Higgs self-coupling
    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // Add modifications due to small variations of the SM parameters    
    //dwidth += cHSM * ( 0.0 ); 
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    //dwidth += eHWWint + eHWWpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHWffRatio2() const
{
    double dwidth = 0.0;
     

    //Contributions that are quadratic in the effective coefficients
    return ( dwidth );
    
}

double NPSMEFTd6::GammaHWW4fRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;

    width += deltaGammaHWW4fRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHWW4fRatio2();
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHWW4fRatio1() const
{    
    double dwidth = 0.0;
    
    double C1 = 0.0073;
    
    double CWff, sf;
    
    CWff = ( CiHL3_11 + CiHL3_22 + CiHL3_33 ) * v2_over_LambdaNP2 +
            Nc * ( CiHQ3_11 + CiHQ3_22 ) * v2_over_LambdaNP2;
    
    CWff = CWff/( 3.0 + 2.0*Nc );
    
    sf = 90362.5 * (1.0/2.0) * ( 3.0 + 2.0*Nc )/(Nc*v2) ; // Coefficient of the CWff term. From the CiHQ3_11 term in the ME.

    dwidth = ( +121886. * CiHbox / LambdaNP2
                + sf* CWff
                -204009. * CiHD / LambdaNP2
                -91455.7 * CiHW / LambdaNP2
                -382903. * CiHWB / LambdaNP2
                +38314.9 * CiDHW / LambdaNP2
                -4.757 * DeltaGF()
                -13.716 * deltaMwd6()
                -0.963 * deltaGwd6()               
                );
        
//  Linear contribution from Higgs self-coupling
    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // Add modifications due to small variations of the SM parameters    
    dwidth += cHSM * ( -12.271 * deltaMz()
                +13.665 * deltaMh()
                +1.85 * deltaaMZ()
                +0.224 * deltaGmu() ); 
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHWWint + eHWWpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHWW4fRatio2() const
{
    double dwidth = 0.0;
    

    //Contributions that are quadratic in the effective coefficients
    return ( dwidth );
    
}

double NPSMEFTd6::GammaHZZRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;
    
    width += deltaGammaHZZRatio1();
       
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHZZRatio2();   
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHZZRatio1() const
{
    double dwidth = 0.0;
    
//    double C1 = 0.0083;

    dwidth = deltaGammaHZZ4fRatio1();
    
//  Linear contribution from Higgs self-coupling
//    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
//    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
//    dwidth += eHZZint + eHZZpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHZZRatio2() const
{
    double dwidth = 0.0;
    
    //Contributions that are quadratic in the effective coefficients
    dwidth = deltaGammaHZZ4fRatio2();
    
    
    return dwidth;
    
}

double NPSMEFTd6::GammaHZllRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;

    width += deltaGammaHZllRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHZllRatio2();   
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHZllRatio1() const
{
    double dwidth = 0.0;
    
    double C1 = 0.0083;

    dwidth = ( +121715. * CiHbox / LambdaNP2
                +8726.9 * (1.0/2.0) * ( CiHL1_11 + CiHL1_22 ) / LambdaNP2
                -7315.2 * (1.0/2.0) * ( CiHe_11 + CiHe_22 ) / LambdaNP2
                +8726.9 * (1.0/2.0) * ( CiHL3_11 + CiHL3_22 ) / LambdaNP2
                -5544.15 * CiHD / LambdaNP2
                -13560.9 * CiHB / LambdaNP2
                -45585.2 * CiHW / LambdaNP2
                -53507.9 * CiHWB / LambdaNP2
                +16829.2 * CiDHB / LambdaNP2
                +30766.6 * CiDHW / LambdaNP2
                -2.204 * DeltaGF()  );
    
//  Linear contribution from Higgs self-coupling
    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // Add modifications due to small variations of the SM parameters    
    //dwidth += cHSM * ( 0.0 ); 
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    //dwidth += eHZZint + eHZZpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHZllRatio2() const
{
    double dwidth = 0.0;
     

    //Contributions that are quadratic in the effective coefficients
    return ( dwidth );
    
}

double NPSMEFTd6::GammaHZeeRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;

    width += deltaGammaHZeeRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHZeeRatio2();   
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHZeeRatio1() const
{
    double dwidth = 0.0;
    
    double C1 = 0.0083;
    
//  Derived from the HZll expression for l=e only

    dwidth = ( +121715. * CiHbox / LambdaNP2
                +8726.9 * CiHL1_11 / LambdaNP2
                -7315.2 * CiHe_11 / LambdaNP2
                +8726.9 * CiHL3_11 / LambdaNP2
                -5544.15 * CiHD / LambdaNP2
                -13560.9 * CiHB / LambdaNP2
                -45585.2 * CiHW / LambdaNP2
                -53507.9 * CiHWB / LambdaNP2
                +16829.2 * CiDHB / LambdaNP2
                +30766.6 * CiDHW / LambdaNP2
                -2.204 * DeltaGF()  );
    
//  Linear contribution from Higgs self-coupling
    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // Add modifications due to small variations of the SM parameters    
    //dwidth += cHSM * ( 0.0 ); 
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    //dwidth += eHZZint + eHZZpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHZeeRatio2() const
{
    double dwidth = 0.0;
    
 
    //Contributions that are quadratic in the effective coefficients
    return ( dwidth );
    
}

double NPSMEFTd6::GammaHZmumuRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;

    width += deltaGammaHZmumuRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHZmumuRatio2();   
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHZmumuRatio1() const
{
    double dwidth = 0.0;
    
    double C1 = 0.0083;
    
//  Derived from the HZll expression for l=mu only

    dwidth = ( +121715. * CiHbox / LambdaNP2
                +8726.9 * CiHL1_22 / LambdaNP2
                -7315.2 * CiHe_22 / LambdaNP2
                +8726.9 * CiHL3_22 / LambdaNP2
                -5544.15 * CiHD / LambdaNP2
                -13560.9 * CiHB / LambdaNP2
                -45585.2 * CiHW / LambdaNP2
                -53507.9 * CiHWB / LambdaNP2
                +16829.2 * CiDHB / LambdaNP2
                +30766.6 * CiDHW / LambdaNP2
                -2.204 * DeltaGF()  );
    
//  Linear contribution from Higgs self-coupling
    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // Add modifications due to small variations of the SM parameters    
    //dwidth += cHSM * ( 0.0 ); 
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    //dwidth += eHZZint + eHZZpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHZmumuRatio2() const
{
    double dwidth = 0.0;
     

    //Contributions that are quadratic in the effective coefficients
    return ( dwidth );
    
}

double NPSMEFTd6::GammaHZZ4lRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;

    width += deltaGammaHZZ4lRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHZZ4lRatio2();   
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHZZ4lRatio1() const
{
    double dwidth = 0.0;
    
    double C1 = 0.0083;
    
    double CZll, sf;
    
    CZll = gZlL*(-0.5 * (CiHL1_11 + CiHL1_22 + CiHL3_11 + CiHL3_22) * v2_over_LambdaNP2) +
            gZlR*(-0.5 * (CiHe_11 + CiHe_22) * v2_over_LambdaNP2);
    
    CZll = CZll/(2.0*(gZlL*gZlL + gZlR*gZlR));
    
    sf = 124479. * (1.0/2.0) * (2.0*(gZlL*gZlL + gZlR*gZlR))/(-0.5*gZlL*v2) ; // Coefficient of the CZll term. From the CiHL1_11 term in the ME.
    
    dwidth = ( +122273. * CiHbox / LambdaNP2
                + sf*CZll
                -44025.7 * CiHD / LambdaNP2
                -13602.6 * CiHB / LambdaNP2
                -45248.6 * CiHW / LambdaNP2
                -88372.1 * CiHWB / LambdaNP2
                +16088.6 * CiDHB / LambdaNP2
                +29210.1 * CiDHW / LambdaNP2
                -3.462 * DeltaGF()
                -0.808 * deltaGzd6() );    
    
//  Linear contribution from Higgs self-coupling
    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // Add modifications due to small variations of the SM parameters    
    dwidth += cHSM * ( -9.734 * deltaMz()
                +15.37 * deltaMh()
                -0.154 * deltaaMZ()
                +2.339 * deltaGmu() ); 
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHZZint + eHZZpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHZZ4lRatio2() const
{
    double dwidth = 0.0;
    

    //Contributions that are quadratic in the effective coefficients
    return ( dwidth );
    
}

double NPSMEFTd6::GammaHZZ4eRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;

    width += deltaGammaHZZ4eRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHZZ4eRatio2();   
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHZZ4eRatio1() const
{
    double dwidth = 0.0;
    
    double C1 = 0.0083;
    
    dwidth = ( +121386. * CiHbox / LambdaNP2
                +123413. * CiHL1_11 / LambdaNP2
                -103717. * CiHe_11 / LambdaNP2
                +123413. * CiHL3_11 / LambdaNP2
                -44056.9 * CiHD / LambdaNP2
                -13385.3 * CiHB / LambdaNP2
                -45127.7 * CiHW / LambdaNP2
                -91708.7 * CiHWB / LambdaNP2
                +16138.9 * CiDHB / LambdaNP2
                +28759.4 * CiDHW / LambdaNP2
                -3.462 * DeltaGF()
                -0.769 * deltaGzd6() );
    
//  Linear contribution from Higgs self-coupling
    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // Add modifications due to small variations of the SM parameters    
    dwidth += cHSM * ( -9.228 * deltaMz()
                +15.148 * deltaMh()
                -0.229 * deltaaMZ()
                +2.493 * deltaGmu() ); 
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHZZint + eHZZpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHZZ4eRatio2() const
{
    double dwidth = 0.0;
     

    //Contributions that are quadratic in the effective coefficients
    return ( dwidth );
    
}

double NPSMEFTd6::GammaHZZ2e2muRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;

    width += deltaGammaHZZ2e2muRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHZZ2e2muRatio2();   
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHZZ2e2muRatio1() const
{
    double dwidth = 0.0;
    
    double C1 = 0.0083;
    
    dwidth = ( +120836. * CiHbox / LambdaNP2
                +126374. * (1.0/2.0) * ( CiHL1_11 + CiHL1_22 ) / LambdaNP2
                -109064. * (1.0/2.0) * ( CiHe_11 + CiHe_22 ) / LambdaNP2
                +126374. * (1.0/2.0) * ( CiHL3_11 + CiHL3_22 ) / LambdaNP2
                -42370.4 * CiHD / LambdaNP2
                -14299. * CiHB / LambdaNP2
                -47298.2 * CiHW / LambdaNP2
                -83098.2 * CiHWB / LambdaNP2
                +16362.7 * CiDHB / LambdaNP2
                +29503.4 * CiDHW / LambdaNP2
                -3.378 * DeltaGF()
                -0.85 * deltaGzd6() );
    
//  Linear contribution from Higgs self-coupling
    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
        
    // Add modifications due to small variations of the SM parameters    
    dwidth += cHSM * ( -10.07 * deltaMz()
                +15.626 * deltaMh()
                -0.128 * deltaaMZ()
                +2.258 * deltaGmu() ); 
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHZZint + eHZZpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHZZ2e2muRatio2() const
{
    double dwidth = 0.0;

    //Contributions that are quadratic in the effective coefficients
    return ( dwidth );
    
}

double NPSMEFTd6::GammaHZZ4muRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;

    width += deltaGammaHZZ4muRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHZZ4muRatio2();   
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHZZ4muRatio1() const
{
    double dwidth = 0.0;
    
    double C1 = 0.0083;
        
    dwidth = ( +120688. * CiHbox / LambdaNP2
                +123059. * CiHL1_22 / LambdaNP2
                -103862. * CiHe_22 / LambdaNP2
                +123059. * CiHL3_22 / LambdaNP2
                -43977.1 * CiHD / LambdaNP2
                -13575.5 * CiHB / LambdaNP2
                -45200.8 * CiHW / LambdaNP2
                -91625.2 * CiHWB / LambdaNP2
                +15449.3 * CiDHB / LambdaNP2
                +28489.5 * CiDHW / LambdaNP2
                -3.471 * DeltaGF()
                -0.774 * deltaGzd6() );
    
//  Linear contribution from Higgs self-coupling
    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // Add modifications due to small variations of the SM parameters    
    dwidth += cHSM * ( -9.254 * deltaMz()
                +15.109 * deltaMh()
                -0.207 * deltaaMZ()
                +2.405 * deltaGmu() ); 
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHZZint + eHZZpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHZZ4muRatio2() const
{
    double dwidth = 0.0;
     

    //Contributions that are quadratic in the effective coefficients
    return ( dwidth );
    
}

double NPSMEFTd6::GammaHZvvRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;

    width += deltaGammaHZvvRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHZvvRatio2();   
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHZvvRatio1() const
{
    double dwidth = 0.0;

    double C1 = 0.0083;
    
    dwidth = ( +121530. * CiHbox / LambdaNP2
                -7943.34 * (1.0/3.0) * ( CiHL1_11 + CiHL1_22 + CiHL1_33 ) / LambdaNP2
                +7943.34 * (1.0/3.0) * ( CiHL3_11 + CiHL3_22 + CiHL3_33 ) / LambdaNP2
                -229.41 * CiHD / LambdaNP2
                -13535.2 * CiHB / LambdaNP2
                -45480.6 * CiHW / LambdaNP2
                -24891. * CiHWB / LambdaNP2
                +16833. * CiDHB / LambdaNP2
                +30597.6 * CiDHW / LambdaNP2
                -2. * DeltaGF()  );
    
//  Linear contribution from Higgs self-coupling
    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // Add modifications due to small variations of the SM parameters    
    //dwidth += cHSM * ( 0.0 ); 
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    //dwidth += eHZZint + eHZZpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHZvvRatio2() const
{
    double dwidth = 0.0;
     

    //Contributions that are quadratic in the effective coefficients
    return ( dwidth );
    
}

double NPSMEFTd6::GammaHZZ4vRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;

    width += deltaGammaHZZ4vRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHZZ4vRatio2();   
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHZZ4vRatio1() const
{
    double dwidth = 0.0;

    double C1 = 0.0083;
    
    dwidth = ( +120596. * CiHbox / LambdaNP2
                -115532. * (1.0/3.0) * ( CiHL1_11 + CiHL1_22 + CiHL1_33 ) / LambdaNP2
                +115532. * (1.0/3.0) * ( CiHL3_11 + CiHL3_22 + CiHL3_33 ) / LambdaNP2
                -28744.1 * CiHD / LambdaNP2
                -13816.7 * CiHB / LambdaNP2
                -44782.1 * CiHW / LambdaNP2
                -25256.6 * CiHWB / LambdaNP2
                +15982.5 * CiDHB / LambdaNP2
                +28910.7 * CiDHW / LambdaNP2
                -3.013 * DeltaGF()
                -0.787 * deltaGzd6()
                 );
    
//  Linear contribution from Higgs self-coupling
    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // Add modifications due to small variations of the SM parameters    
    dwidth += cHSM * ( -10.49 * deltaMz()
                +15.294 * deltaMh()
                +0.255 * deltaaMZ()
                +1.979 * deltaGmu() ); 
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHZZint + eHZZpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHZZ4vRatio2() const
{
    double dwidth = 0.0;
     

    //Contributions that are quadratic in the effective coefficients
    return ( dwidth );
    
}

double NPSMEFTd6::GammaHZuuRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;

    width += deltaGammaHZuuRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHZuuRatio2();   
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHZuuRatio1() const
{
    double dwidth = 0.0;
    
    double C1 = 0.0083;

    dwidth = ( +121512. * CiHbox / LambdaNP2
                -9648.28 * (1.0/2.0) * ( CiHQ1_11 + CiHQ1_22 ) / LambdaNP2
                +4218.6 * (1.0/2.0) * ( CiHu_11 + CiHu_22 ) / LambdaNP2
                +9648.28 * (1.0/2.0) * ( CiHQ3_11 + CiHQ3_22 ) / LambdaNP2
                -17762.5 * CiHD / LambdaNP2
                -13473.3 * CiHB / LambdaNP2
                -45667.9 * CiHW / LambdaNP2
                -110057. * CiHWB / LambdaNP2
                +16854.2 * CiDHB / LambdaNP2
                +30781.7 * CiDHW / LambdaNP2
                -2.6 * DeltaGF()  );
    
//  Linear contribution from Higgs self-coupling
    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // Add modifications due to small variations of the SM parameters    
    //dwidth += cHSM * ( 0.0 ); 
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    //dwidth += eHZZint + eHZZpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHZuuRatio2() const
{
    double dwidth = 0.0;
     

    //Contributions that are quadratic in the effective coefficients
    return ( dwidth );
    
}

double NPSMEFTd6::GammaHZddRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;

    width += deltaGammaHZddRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHZddRatio2();   
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHZddRatio1() const
{
    double dwidth = 0.0;
    
    double C1 = 0.0083;

    dwidth = ( +121756. * CiHbox / LambdaNP2
                +9252.73 * (1.0/3.0) * ( CiHQ1_11 + CiHQ1_22 + CiHQ1_33 ) / LambdaNP2
                -1471.03 * (1.0/3.0) * ( CiHd_11 + CiHd_22 + CiHd_33 ) / LambdaNP2
                +9252.73 * (1.0/3.0) * ( CiHQ3_11 + CiHQ3_22 + CiHQ3_33 ) / LambdaNP2
                -12714.3 * CiHD / LambdaNP2
                -13589.3 * CiHB / LambdaNP2
                -45689.4 * CiHW / LambdaNP2
                -85582.3 * CiHWB / LambdaNP2
                +17007.1 * CiDHB / LambdaNP2
                +30733.1 * CiDHW / LambdaNP2
                -2.427 * DeltaGF()  );
    
//  Linear contribution from Higgs self-coupling
    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // Add modifications due to small variations of the SM parameters    
    //dwidth += cHSM * ( 0.0 ); 
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    //dwidth += eHZZint + eHZZpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHZddRatio2() const
{
    double dwidth = 0.0;
     

    //Contributions that are quadratic in the effective coefficients
    return ( dwidth );
    
}

double NPSMEFTd6::GammaHZffRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;

    width += deltaGammaHZffRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHZffRatio2();   
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHZffRatio1() const
{
    double dwidth = 0.0;
    
    double C1 = 0.0083;

    dwidth = ( +121551. * CiHbox / LambdaNP2
                -824.482 * (1.0/3.0) * ( CiHL1_11 + CiHL1_22 + CiHL1_33 ) / LambdaNP2
                +1840.54 * (1.0/12.0) * ( 5.0 * CiHQ1_11 + 5.0 * CiHQ1_22 + 2.0 * CiHQ1_33 - CiHQ3_11 - CiHQ3_22 + 2.0 * CiHQ3_33 ) / LambdaNP2
                -795.383 * (1.0/3.0) * ( CiHe_11 + CiHe_22 + CiHe_33 ) / LambdaNP2
                +1069.4 * (1.0/2.0) * ( CiHu_11 + CiHu_22 ) / LambdaNP2
                -579.563 * (1.0/3.0) * ( CiHd_11 + CiHd_22 + CiHd_33 ) / LambdaNP2
                +3164.56 * (1.0/3.0) * ( CiHL3_11 + CiHL3_22 + CiHL3_33 ) / LambdaNP2
                +6413.99 * (-1.0/12.0) * ( CiHQ1_11 + CiHQ1_22 - 2.0 * CiHQ1_33 - 5.0 * CiHQ3_11 - 5.0 * CiHQ3_22 - 2.0 * CiHQ3_33) / LambdaNP2
                -10839.5 * CiHD / LambdaNP2
                -14222.3 * CiHB / LambdaNP2
                -45455.6 * CiHW / LambdaNP2
                -75343.1 * CiHWB / LambdaNP2
                +16804.9 * CiDHB / LambdaNP2
                +30421. * CiDHW / LambdaNP2
                -2.356 * DeltaGF()  );
    
//  Linear contribution from Higgs self-coupling
    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // Add modifications due to small variations of the SM parameters    
    //dwidth += cHSM * ( 0.0 ); 
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    //dwidth += eHZZint + eHZZpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHZffRatio2() const
{
    double dwidth = 0.0;
     

    //Contributions that are quadratic in the effective coefficients
    return ( dwidth );
    
}

double NPSMEFTd6::GammaHZZ4fRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;

    width += deltaGammaHZZ4fRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHZZ4fRatio2();   
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHZZ4fRatio1() const
{
    double dwidth = 0.0;
    
    double C1 = 0.0083;
    
    double CZff, sf;
    
    CZff = gZvL*(-0.5 * (CiHL1_11 + CiHL1_22 + CiHL1_33 - CiHL3_11 - CiHL3_22 - CiHL3_33) * v2_over_LambdaNP2) +
            gZlL*(-0.5 * (CiHL1_11 + CiHL1_22 + CiHL1_33 + CiHL3_11 + CiHL3_22 + CiHL3_33) * v2_over_LambdaNP2) +
            gZlR*(-0.5 * (CiHe_11 + CiHe_22 + CiHe_33) * v2_over_LambdaNP2) +
            Nc * (
            gZdL*(-0.5 * (CiHQ1_11 + CiHQ1_22 + CiHQ1_33 + CiHQ3_11 + CiHQ3_22 + CiHQ3_33) * v2_over_LambdaNP2) +          
            gZdR*(-0.5 * (CiHd_11 + CiHd_22 + CiHd_33) * v2_over_LambdaNP2) +
            gZuL*(-0.5 * (CiHQ1_11 + CiHQ1_22 - CiHQ3_11 - CiHQ3_22) * v2_over_LambdaNP2) +          
            gZuR*(-0.5 * (CiHu_11 + CiHu_22) * v2_over_LambdaNP2)
            );
    
    CZff = CZff/(            
            3.0*( gZvL*gZvL + gZlL*gZlL + gZlR*gZlR ) +
            Nc * ( 3.0*( gZdL*gZdL + gZdR*gZdR ) + 2.0*( gZuL*gZuL + gZuR*gZuR ) )           
            );
    
    sf = -11267.6 * (1.0/3.0) * (            
            3.0*( gZvL*gZvL + gZlL*gZlL + gZlR*gZlR ) +
            Nc * ( 3.0*( gZdL*gZdL + gZdR*gZdR ) + 2.0*( gZuL*gZuL + gZuR*gZuR ) )           
            );
            
    sf = sf/(-0.5*(gZlL + gZvL)*v2) ; // Coefficient of the CZff term. From the CiHL1_11 term in the ME.
    
    dwidth = ( +121373. * CiHbox / LambdaNP2
                + sf*CZff
                -50927.1 * CiHD / LambdaNP2
                -14137.9 * CiHB / LambdaNP2
                -46350.1 * CiHW / LambdaNP2
                -126336. * CiHWB / LambdaNP2
                +16558.7 * CiDHB / LambdaNP2
                +29628.7 * CiDHW / LambdaNP2
                -3.715 * DeltaGF()
                -0.834 * deltaGzd6()
                  );    
    
//  Linear contribution from Higgs self-coupling
    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // Add modifications due to small variations of the SM parameters    
    dwidth += cHSM * ( -9.548 * deltaMz()
                +15.799 * deltaMh()
                -0.412 * deltaaMZ()
                +2.569 * deltaGmu() ); 
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHZZint + eHZZpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHZZ4fRatio2() const
{
    double dwidth = 0.0;
    

    //Contributions that are quadratic in the effective coefficients
    return ( dwidth );
    
}

double NPSMEFTd6::GammaHZgaRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;

    width += deltaGammaHZgaRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHZgaRatio2();      
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHZgaRatio1() const
{
    double dwidth = 0.0;
    
    double C1 = 0.0;

//  It includes modifications of Zff vertices and MW, but not on the pure VVV and VVVV vertices
    
//  Write the tree-level contributions directly as a function 
//  of delta_ZA (or deltaG1_hZA()) to account for variations of sw2 and cw2

    dwidth = ( -71769.02 * deltaG1_hZA() 
//            +14894914. * CiHB / LambdaNP2
//            -14894913. * CiHW / LambdaNP2
//            +9508089. * CiHWB / LambdaNP2
//            -2869576. * CiDHB / LambdaNP2
//            +1572613. * CiDHW / LambdaNP2            
            + cLHd6 * (
            +120002. * CiHbox / LambdaNP2
            +50.12 * CiHL1_33 / LambdaNP2
            +17401. * CiHQ1_33 / LambdaNP2
            +50.12 * CiHe_33 / LambdaNP2
            +17188.7 * CiHu_33 / LambdaNP2
            +212.376 * CiHd_33 / LambdaNP2
            +50.12 * CiHL3_33 / LambdaNP2
            -16976.3 * CiHQ3_33 / LambdaNP2
            -373.856 * CieH_33r / LambdaNP2
            -2953.05 * CiuH_22r / LambdaNP2
            +6636.34 * CiuH_33r / LambdaNP2
            -6121.66 * CidH_33r / LambdaNP2
            -111254. * CiHD / LambdaNP2
            -162538. * CiHWB / LambdaNP2
            -96076.1 * DeltaGF() / v() / v()
            -0.123 * deltaMwd6() )
             );
    
//  Linear contribution from Higgs self-coupling
    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // Add modifications due to small variations of the SM parameters    
    dwidth += cHSM * ( +1. * deltaa0()
                -0.629 * deltaaMZ()
                +2.629 * deltaGmu()
                -4.926 * deltaMz()
                +0.004 * deltaaSMZ()
                +11.167 * deltaMh()
                +0.013 * deltamt()
                +0.004 * deltamb()
                +0.001 * deltamc()
                +0. * deltamtau() ); 
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHZgaint + eHZgapar;
    
    return dwidth;
}

double NPSMEFTd6::deltaGammaHZgaRatio2() const
{
    double dwidth = 0.0;
     

    //Contributions that are quadratic in the effective coefficients
    return ( dwidth );
    
}

double NPSMEFTd6::GammaHgagaRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;

    width += deltaGammaHgagaRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHgagaRatio2();  
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHgagaRatio1() const
{
    double dwidth = 0.0;
    
    double C1 = 0.0049;

//  It does not include modifications of MW
    
//  Write the tree-level contributions directly as a function 
//  of delta_AA (or deltaG_hAA) to account for variations of sw2 and cw2
    
    dwidth = ( -255156.97*deltaG_hAA()            
//            -48314158. * CiHB / LambdaNP2
//            -14510502. * CiHW / LambdaNP2
//            +26477588. * CiHWB / LambdaNP2 
            + cLHd6 * (
            +119766. * CiHbox / LambdaNP2
            -42565.7 * CieH_33r / LambdaNP2
            -48868.1 * CiuH_22r / LambdaNP2
            +32078.2 * CiuH_33r / LambdaNP2
            -18428.3 * CidH_33r / LambdaNP2
            -137452. * CiHD / LambdaNP2
            -235677. * CiHWB / LambdaNP2
            -124462. * DeltaGF() / v() / v() 
            -1.257 * deltaMwd6() )
            );
    
//  Linear contribution from Higgs self-coupling
    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // Add modifications due to small variations of the SM parameters    
    dwidth += cHSM * ( +2. * deltaa0()
                +0.27 * deltaaMZ()
                +0.736 * deltaGmu()
                -1.797 * deltaMz()
                +0.02 * deltaaSMZ()
                +4.195 * deltaMh()
                +0.047 * deltamt()
                +0.008 * deltamb()
                +0.009 * deltamc()
                +0.01 * deltamtau() ); 
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHgagaint + eHgagapar;
    
    return dwidth;
}

double NPSMEFTd6::deltaGammaHgagaRatio2() const
{
    double dwidth = 0.0;
     

    //Contributions that are quadratic in the effective coefficients
    return ( dwidth );
    
}
      
double NPSMEFTd6::GammaHmumuRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;

    width += deltaGammaHmumuRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHmumuRatio2();
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHmumuRatio1() const
{
    double dwidth = 0.0;
    
    double C1 = 0.0;
    
    dwidth = ( +121248. * CiHbox / LambdaNP2
            -199792511. * CieH_22r / LambdaNP2
            -30312.1 * CiHD / LambdaNP2
            -60624.1 * DeltaGF() / v() / v() ); 
    
//  Linear contribution from Higgs self-coupling
    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // Add modifications due to small variations of the SM parameters    
    dwidth += cHSM * ( +1. * deltaGmu()
                +1. * deltaMh() ); 
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHmumuint + eHmumupar;
    
    return dwidth;
}

double NPSMEFTd6::deltaGammaHmumuRatio2() const
{
    double dwidth = 0.0;
     

    //Contributions that are quadratic in the effective coefficients
    return ( dwidth );
    
}

double NPSMEFTd6::GammaHtautauRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;

    width += deltaGammaHtautauRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHtautauRatio2();    
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHtautauRatio1() const
{
    double dwidth = 0.0;
    
    double C1 = 0.0;

    dwidth = ( +121248. * CiHbox / LambdaNP2
            -11880369. * CieH_33r / LambdaNP2
            -30312.1 * CiHD / LambdaNP2
            -60624.1 * DeltaGF() / v() / v() );
    
//  Linear contribution from Higgs self-coupling
    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // Add modifications due to small variations of the SM parameters    
    dwidth += cHSM * ( +1. * deltaGmu()
                +1.002 * deltaMh()
                +1.998 * deltamtau() ); 
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHtautauint + eHtautaupar;
    
    return dwidth;
}

double NPSMEFTd6::deltaGammaHtautauRatio2() const
{
    double dwidth = 0.0;
     

    //Contributions that are quadratic in the effective coefficients
    return ( dwidth );
    
}

double NPSMEFTd6::GammaHccRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;

    width += deltaGammaHccRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHccRatio2();   
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHccRatio1() const
{
    double dwidth = 0.0;
    
    double C1 = 0.0;
    
    if (FlagLoopHd6) {
        
        dwidth = ( +121248. * CiHbox / LambdaNP2
            -16421890. * CiuH_22r / LambdaNP2
            -992.159 * CiuH_33r / LambdaNP2
            -30312.1 * CiHD / LambdaNP2
            -60624.1 * DeltaGF() / v() / v() );
        
    } else {
        
        dwidth = ( +121248. * CiHbox / LambdaNP2
            -16556668. * CiuH_22r / LambdaNP2
            -30312.1 * CiHD / LambdaNP2
            -60624.1 * DeltaGF() / v() / v() );  
    }
    
//  Linear contribution from Higgs self-coupling
    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // Add modifications due to small variations of the SM parameters    
    dwidth += cHSM * ( +1. * deltaGmu()
                -0.789 * deltaaSMZ()
                +1.004 * deltaMh()
                +0.001 * deltamt()
                +1.995 * deltamc() ); 
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHccint + eHccpar;
    
    return dwidth;
}

double NPSMEFTd6::deltaGammaHccRatio2() const
{
    double dwidth = 0.0;
     

    //Contributions that are quadratic in the effective coefficients
    return ( dwidth );
    
}

double NPSMEFTd6::GammaHbbRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;
    
    width += deltaGammaHbbRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHbbRatio2();    
        }
    
    return width;
}

double NPSMEFTd6::deltaGammaHbbRatio1() const
{
    double dwidth = 0.0;
    
    double C1 = 0.0;

    if (FlagLoopHd6) {
        
        dwidth = ( +121248. * CiHbox / LambdaNP2
            -558.186 * CiuH_33r / LambdaNP2
            -5027051. * CidH_33r / LambdaNP2
            -30312.1 * CiHD / LambdaNP2
            -60624.1 * DeltaGF() / v() / v() );
        
    } else {
        
        dwidth = ( +121248. * CiHbox / LambdaNP2
            -5050180. * CidH_33r / LambdaNP2
            -30312.1 * CiHD / LambdaNP2
            -60624.1 * DeltaGF() / v() / v() );  
    }

//  Linear contribution from Higgs self-coupling
    dwidth = dwidth + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    dwidth = dwidth + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    // Add modifications due to small variations of the SM parameters    
    dwidth += cHSM * ( +1. * deltaGmu()
                -0.23 * deltaaSMZ()
                +1.007 * deltaMh()
                +0.001 * deltamt()
                +1.992 * deltamb() ); 
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHbbint + eHbbpar;
    
    return dwidth;
}

double NPSMEFTd6::deltaGammaHbbRatio2() const
{
    double dwidth = 0.0;
     

    //Contributions that are quadratic in the effective coefficients
    return ( dwidth );        
    
}

double NPSMEFTd6::Br_H_exo() const
{    
    if (BrHexo < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return BrHexo;
}

double NPSMEFTd6::Br_H_inv() const
{    
//  Contributions from both modifications in H->ZZ->4v and the extra invisible decays
    double BR4v;
    
    BR4v = BrHZZ4vRatio()*(trueSM.computeBrHtoZZinv());
    
//  BR4v positivity is already checked inside BrHZZ4vRatio()
//  and will be nan if negative. Check here BrHinv, to make sure both are positive
    if (BrHinv < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return BR4v + BrHinv;
}


double NPSMEFTd6::Br_H_inv_NP() const
{    
    
//  BR4v positivity is already checked inside BrHZZ4vRatio()
//  and will be nan if negative. Check here BrHinv, to make sure both are positive
    if (BrHinv < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return BrHinv;
}


double NPSMEFTd6::BrHvisRatio() const
{    
    double Br = 1.0;
    double dvis1 = 0.0, dvis2 = 0.0, delta2SM;
    double GHvisR = 1.0;

//  Sum over decays of visible SM and exotic modes
    dvis1 = (trueSM.computeBrHtogg() * deltaGammaHggRatio1()
            + trueSM.computeBrHtoWW() * deltaGammaHWWRatio1()
            + trueSM.computeBrHtoZZ() * deltaGammaHZZRatio1()
            + trueSM.computeBrHtoZga() * deltaGammaHZgaRatio1()
            + trueSM.computeBrHtogaga() * deltaGammaHgagaRatio1()
            + trueSM.computeBrHtomumu() * deltaGammaHmumuRatio1()
            + trueSM.computeBrHtotautau() * deltaGammaHtautauRatio1()
            + trueSM.computeBrHtocc() * deltaGammaHccRatio1()
            + trueSM.computeBrHtobb() * deltaGammaHbbRatio1()
            + BrHexo);    
    
    Br += dvis1 - dGammaHTotR1;
    
    if (FlagQuadraticTerms) {

//  Sum over decays of visible SM and exotic modes
        delta2SM = trueSM.computeBrHtogg() * deltaGammaHggRatio2()
            + trueSM.computeBrHtoWW() * deltaGammaHWWRatio2()
            + trueSM.computeBrHtoZZ() * deltaGammaHZZRatio2()
            + trueSM.computeBrHtoZga() * deltaGammaHZgaRatio2()
            + trueSM.computeBrHtogaga() * deltaGammaHgagaRatio2()
            + trueSM.computeBrHtomumu() * deltaGammaHmumuRatio2()
            + trueSM.computeBrHtotautau() * deltaGammaHtautauRatio2()
            + trueSM.computeBrHtocc() * deltaGammaHccRatio2()
            + trueSM.computeBrHtobb() * deltaGammaHbbRatio2();

        dvis2 = delta2SM + (BrHexo)*(BrHexo + delta2SM);

        //Add contributions that are quadratic in the effective coefficients
        Br += - dvis1 * dGammaHTotR1
                + dvis2 - dGammaHTotR2
                + pow(dGammaHTotR1,2.0);            
        }
    
    GHvisR += dvis1 + dvis2; 
    if ((Br < 0) || (GHvisR < 0) || (GammaHTotR < 0)) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;
}


double NPSMEFTd6::BrHtoinvRatio() const
{        
    return (Br_H_inv()/(trueSM.computeBrHtoZZinv()));
}


///////////////////////SPECIAL OBSERVABLES/////////////////////////

double NPSMEFTd6::muttHZbbboost(const double sqrt_s) const
{
    /* Ratios of BR with the SM*/
    double BrHbbrat = BrHbbRatio();
    double BrZbbSM = (trueSM.GammaZ(quarks[BOTTOM]))/trueSM.Gamma_Z();
    double BrZbbrat = BR_Zf(quarks[BOTTOM])/BrZbbSM;

//    gslpp::complex dKappa_t = deltaG_hff(quarks[TOP]) / (-mtpole / v());    
//    double dkt = dKappa_t.real();
    
//    double dgV = deltaGV_f(quarks[TOP]);
//    double dgA = deltaGA_f(quarks[TOP]);
//    double gLSM = quarks[TOP].getIsospin() 
//    - (quarks[TOP].getCharge())*sW2_tree;
//    double gRSM = - (quarks[TOP].getCharge())*sW2_tree;
    
//    double dgL = 0.5*(dgV + dgA)/gLSM;
//    double dgR = 0.5*(dgV - dgA)/gRSM;

    double dsigmarat;

    /* VERY CRUDE APPROX. */    
    //dsigmarat = 1.0 + 
    //        2.0 * dkt -
    //        2.0 * (gLSM*gLSM*dgL + gRSM*gRSM*dgR)/(gLSM*gLSM + gRSM*gRSM);

    dsigmarat = 1.0;
//  ttH 100 TeV (from muttH func): NOT BOOSTED YET
    dsigmarat += +467438. * CHG / LambdaNP2
                -22519. * CG / LambdaNP2
                +880378. * CiuG_33r / LambdaNP2
                -2.837 * deltaG_hff(quarks[TOP]).real()
                ;
//  Divided (linearized) by ttZ 100 TeV
    dsigmarat = dsigmarat - (
                -40869.4 * CiHD / LambdaNP2
                -52607.9 * CiHWB / LambdaNP2
                -90424.9 * CHG / LambdaNP2
                +432089. * CG / LambdaNP2
                +326525. * CiuG_33r / LambdaNP2
                -2028.11 * CiuW_33r / LambdaNP2
                +1679.85 * CiuB_33r / LambdaNP2
                +1454.5 * CiHQ1_11 / LambdaNP2
                +1065.27 * CiHu_11 / LambdaNP2
                +82169.1 * CiHu_33 / LambdaNP2
                -1229.16 * CiHd_11 / LambdaNP2
                +6780.84 * CiHQ3_11 / LambdaNP2
                -1.374 * DeltaGF()
                +4.242 * -0.5 * (CiHQ1_33 - CiHQ3_33) * v2_over_LambdaNP2            
            );    
    
    return dsigmarat * (BrHbbrat / BrZbbrat);
    
}


double NPSMEFTd6::muggHgaga(const double sqrt_s) const
{   
    return muggH(sqrt_s) * BrHgagaRatio();
    
}

double NPSMEFTd6::muVBFHgaga(const double sqrt_s) const
{   
    return muVBF(sqrt_s) * BrHgagaRatio();
    
}

double NPSMEFTd6::muZHgaga(const double sqrt_s) const
{   
    return muZH(sqrt_s) * BrHgagaRatio();
    
}

double NPSMEFTd6::muWHgaga(const double sqrt_s) const
{   
    return muWH(sqrt_s) * BrHgagaRatio();
    
}

double NPSMEFTd6::muVHgaga(const double sqrt_s) const
{   
    return muVH(sqrt_s) * BrHgagaRatio();
    
}

double NPSMEFTd6::muttHgaga(const double sqrt_s) const
{   
    return muttH(sqrt_s) * BrHgagaRatio();
    
}

double NPSMEFTd6::muggHZga(const double sqrt_s) const
{   
    return muggH(sqrt_s) * BrHZgaRatio();
    
}

double NPSMEFTd6::muVBFHZga(const double sqrt_s) const
{   
    return muVBF(sqrt_s) * BrHZgaRatio();
    
}

double NPSMEFTd6::muZHZga(const double sqrt_s) const
{   
    return muZH(sqrt_s) * BrHZgaRatio();
    
}

double NPSMEFTd6::muWHZga(const double sqrt_s) const
{   
    return muWH(sqrt_s) * BrHZgaRatio();
    
}

double NPSMEFTd6::muVHZga(const double sqrt_s) const
{   
    return muVH(sqrt_s) * BrHZgaRatio();
    
}

double NPSMEFTd6::muttHZga(const double sqrt_s) const
{   
    return muttH(sqrt_s) * BrHZgaRatio();
    
}

double NPSMEFTd6::muggHZZ(const double sqrt_s) const
{   
    return muggH(sqrt_s) * BrHZZRatio();
    
}

double NPSMEFTd6::muVBFHZZ(const double sqrt_s) const
{   
    return muVBF(sqrt_s) * BrHZZRatio();
    
}

double NPSMEFTd6::muZHZZ(const double sqrt_s) const
{   
    return muZH(sqrt_s) * BrHZZRatio();
    
}

double NPSMEFTd6::muWHZZ(const double sqrt_s) const
{   
    return muWH(sqrt_s) * BrHZZRatio();
    
}

double NPSMEFTd6::muVHZZ(const double sqrt_s) const
{   
    return muVH(sqrt_s) * BrHZZRatio();
    
}

double NPSMEFTd6::muttHZZ(const double sqrt_s) const
{   
    return muttH(sqrt_s) * BrHZZRatio();
    
}

double NPSMEFTd6::muggHZZ4l(const double sqrt_s) const
{   
    return muggH(sqrt_s) * BrHZZ4lRatio();
    
}

double NPSMEFTd6::muVBFHZZ4l(const double sqrt_s) const
{   
    return muVBF(sqrt_s) * BrHZZ4lRatio();
    
}

double NPSMEFTd6::muZHZZ4l(const double sqrt_s) const
{   
    return muZH(sqrt_s) * BrHZZ4lRatio();
    
}

double NPSMEFTd6::muWHZZ4l(const double sqrt_s) const
{   
    return muWH(sqrt_s) * BrHZZ4lRatio();
    
}

double NPSMEFTd6::muVHZZ4l(const double sqrt_s) const
{   
    return muVH(sqrt_s) * BrHZZ4lRatio();
    
}

double NPSMEFTd6::muttHZZ4l(const double sqrt_s) const
{   
    return muttH(sqrt_s) * BrHZZ4lRatio();
    
}

double NPSMEFTd6::muggHWW(const double sqrt_s) const
{   
    return muggH(sqrt_s) * BrHWWRatio();
    
}

double NPSMEFTd6::muVBFHWW(const double sqrt_s) const
{   
    return muVBF(sqrt_s) * BrHWWRatio();
    
}

double NPSMEFTd6::muZHWW(const double sqrt_s) const
{   
    return muZH(sqrt_s) * BrHWWRatio();
    
}

double NPSMEFTd6::muWHWW(const double sqrt_s) const
{   
    return muWH(sqrt_s) * BrHWWRatio();
    
}

double NPSMEFTd6::muVHWW(const double sqrt_s) const
{   
    return muVH(sqrt_s) * BrHWWRatio();
    
}

double NPSMEFTd6::muttHWW(const double sqrt_s) const
{   
    return muttH(sqrt_s) * BrHWWRatio();
    
}

double NPSMEFTd6::muggHWW2l2v(const double sqrt_s) const
{   
    return muggH(sqrt_s) * BrHWW2l2vRatio();
    
}

double NPSMEFTd6::muVBFHWW2l2v(const double sqrt_s) const
{   
    return muVBF(sqrt_s) * BrHWW2l2vRatio();
    
}

double NPSMEFTd6::muZHWW2l2v(const double sqrt_s) const
{   
    return muZH(sqrt_s) * BrHWW2l2vRatio();
    
}

double NPSMEFTd6::muWHWW2l2v(const double sqrt_s) const
{   
    return muWH(sqrt_s) * BrHWW2l2vRatio();
    
}

double NPSMEFTd6::muVHWW2l2v(const double sqrt_s) const
{   
    return muVH(sqrt_s) * BrHWW2l2vRatio();
    
}

double NPSMEFTd6::muttHWW2l2v(const double sqrt_s) const
{   
    return muttH(sqrt_s) * BrHWW2l2vRatio();
    
}

double NPSMEFTd6::muggHmumu(const double sqrt_s) const
{   
    return muggH(sqrt_s) * BrHmumuRatio();
    
}

double NPSMEFTd6::muVBFHmumu(const double sqrt_s) const
{   
    return muVBF(sqrt_s) * BrHmumuRatio();
    
}

double NPSMEFTd6::muZHmumu(const double sqrt_s) const
{   
    return muZH(sqrt_s) * BrHmumuRatio();
    
}

double NPSMEFTd6::muWHmumu(const double sqrt_s) const
{   
    return muWH(sqrt_s) * BrHmumuRatio();
    
}

double NPSMEFTd6::muVHmumu(const double sqrt_s) const
{   
    return muVH(sqrt_s) * BrHmumuRatio();
    
}

double NPSMEFTd6::muttHmumu(const double sqrt_s) const
{   
    return muttH(sqrt_s) * BrHmumuRatio();
    
}

double NPSMEFTd6::muggHtautau(const double sqrt_s) const
{   
    return muggH(sqrt_s) * BrHtautauRatio();
    
}

double NPSMEFTd6::muVBFHtautau(const double sqrt_s) const
{   
    return muVBF(sqrt_s) * BrHtautauRatio();
    
}

double NPSMEFTd6::muZHtautau(const double sqrt_s) const
{   
    return muZH(sqrt_s) * BrHtautauRatio();
    
}

double NPSMEFTd6::muWHtautau(const double sqrt_s) const
{   
    return muWH(sqrt_s) * BrHtautauRatio();
    
}

double NPSMEFTd6::muVHtautau(const double sqrt_s) const
{   
    return muVH(sqrt_s) * BrHtautauRatio();
    
}

double NPSMEFTd6::muttHtautau(const double sqrt_s) const
{   
    return muttH(sqrt_s) * BrHtautauRatio();
    
}

double NPSMEFTd6::muggHbb(const double sqrt_s) const
{   
    return muggH(sqrt_s) * BrHbbRatio();
    
}

double NPSMEFTd6::muVBFHbb(const double sqrt_s) const
{   
    return muVBF(sqrt_s) * BrHbbRatio();
    
}

double NPSMEFTd6::muZHbb(const double sqrt_s) const
{   
    return muZH(sqrt_s) * BrHbbRatio();
    
}

double NPSMEFTd6::muWHbb(const double sqrt_s) const
{   
    return muWH(sqrt_s) * BrHbbRatio();
    
}

double NPSMEFTd6::muVHbb(const double sqrt_s) const
{   
    return muVH(sqrt_s) * BrHbbRatio();
    
}

double NPSMEFTd6::muttHbb(const double sqrt_s) const
{   
    return muttH(sqrt_s) * BrHbbRatio();
    
}

////////////////////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------------------
//-- Special Hadron collider signal strengths with separate full TH unc U(prod x decay) ---
//-----------------------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////////////////// 

double NPSMEFTd6::muTHUggHgaga(const double sqrt_s) const
{   
    if (FlagQuadraticTerms) {
        return ( muggH(sqrt_s)*BrHgagaRatio() * (1.0 + eggFHgaga ) * (1.0 + eHwidth)/(1.0 + eggFint + eggFpar)/(1.0 + eHgagaint + eHgagapar) );
    } else {
        return ( muggH(sqrt_s) + BrHgagaRatio() - 1.0 + eggFHgaga - eggFint - eggFpar - eHgagaint - eHgagapar + eHwidth );
    }    
}
   
double NPSMEFTd6::muTHUVBFHgaga(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muVBF(sqrt_s)*BrHgagaRatio() * (1.0 + eVBFHgaga ) * (1.0 + eHwidth)/(1.0 + eVBFint + eVBFpar)/(1.0 + eHgagaint + eHgagapar) );
    } else {
        return ( muVBF(sqrt_s) + BrHgagaRatio() - 1.0 + eVBFHgaga - eVBFint - eVBFpar - eHgagaint - eHgagapar + eHwidth );
    }   
}

double NPSMEFTd6::muTHUZHgaga(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muZH(sqrt_s)*BrHgagaRatio() * (1.0 + eZHgaga ) * (1.0 + eHwidth)/(1.0 + eZHint + eZHpar)/(1.0 + eHgagaint + eHgagapar) );
    } else {
        return ( muZH(sqrt_s) + BrHgagaRatio() - 1.0 + eZHgaga - eZHint - eZHpar - eHgagaint - eHgagapar + eHwidth );
    }   
}

double NPSMEFTd6::muTHUWHgaga(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muWH(sqrt_s)*BrHgagaRatio() * (1.0 + eWHgaga ) * (1.0 + eHwidth)/(1.0 + eWHint + eWHpar)/(1.0 + eHgagaint + eHgagapar) );
    } else {
        return ( muWH(sqrt_s) + BrHgagaRatio() - 1.0 + eWHgaga - eWHint - eWHpar - eHgagaint - eHgagapar + eHwidth );
    }   
}

double NPSMEFTd6::muTHUVHgaga(const double sqrt_s) const
{
    //  Theory uncertainty in VH production, from the WH and ZH ones
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);    
    double eVHtot,eVHgaga;
    
    eVHtot = ((eWHint + eWHpar) * sigmaWH_SM + (eZHint + eZHpar) * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    eVHgaga = (eWHgaga * sigmaWH_SM + eZHgaga * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    if (FlagQuadraticTerms) {
        return ( muVH(sqrt_s)*BrHgagaRatio() * (1.0 + eVHgaga ) * (1.0 + eHwidth)/(1.0 + eVHtot)/(1.0 + eHgagaint + eHgagapar) );
    } else {
        return ( muVH(sqrt_s) + BrHgagaRatio() - 1.0 + eVHgaga - eVHtot - eHgagaint - eHgagapar + eHwidth );
    }   
}

double NPSMEFTd6::muTHUttHgaga(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muttH(sqrt_s)*BrHgagaRatio() * (1.0 + ettHgaga ) * (1.0 + eHwidth)/(1.0 + eeettHint + eeettHpar)/(1.0 + eHgagaint + eHgagapar) );
    } else {
        return ( muttH(sqrt_s) + BrHgagaRatio() - 1.0 + ettHgaga - eeettHint - eeettHpar - eHgagaint - eHgagapar + eHwidth );
    }      
}

double NPSMEFTd6::muTHUggHZga(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muggH(sqrt_s)*BrHZgaRatio() * (1.0 + eggFHZga ) * (1.0 + eHwidth)/(1.0 + eggFint + eggFpar)/(1.0 + eHZgaint + eHZgapar) );
    } else {
        return ( muggH(sqrt_s) + BrHZgaRatio() - 1.0 + eggFHZga - eggFint - eggFpar - eHZgaint - eHZgapar + eHwidth );
    }     
}

double NPSMEFTd6::muTHUVBFHZga(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muVBF(sqrt_s)*BrHZgaRatio() * (1.0 + eVBFHZga ) * (1.0 + eHwidth)/(1.0 + eVBFint + eVBFpar)/(1.0 + eHZgaint + eHZgapar) );
    } else {
        return ( muVBF(sqrt_s) + BrHZgaRatio() - 1.0 + eVBFHZga - eVBFint - eVBFpar - eHZgaint - eHZgapar + eHwidth );
    }     
}

double NPSMEFTd6::muTHUZHZga(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muZH(sqrt_s)*BrHZgaRatio() * (1.0 + eZHZga ) * (1.0 + eHwidth)/(1.0 + eZHint + eZHpar)/(1.0 + eHZgaint + eHZgapar) );
    } else {
        return ( muZH(sqrt_s) + BrHZgaRatio() - 1.0 + eZHZga - eZHint - eZHpar - eHZgaint - eHZgapar + eHwidth );
    }     
}

double NPSMEFTd6::muTHUWHZga(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muWH(sqrt_s)*BrHZgaRatio() * (1.0 + eWHZga ) * (1.0 + eHwidth)/(1.0 + eWHint + eWHpar)/(1.0 + eHZgaint + eHZgapar) );
    } else {
        return ( muWH(sqrt_s) + BrHZgaRatio() - 1.0 + eWHZga - eWHint - eWHpar - eHZgaint - eHZgapar + eHwidth );
    }   
}

double NPSMEFTd6::muTHUVHZga(const double sqrt_s) const
{
    //  Theory uncertainty in VH production, from the WH and ZH ones
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);    
    double eVHtot,eVHZga;
    
    eVHtot = ((eWHint + eWHpar) * sigmaWH_SM + (eZHint + eZHpar) * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    eVHZga = (eWHZga * sigmaWH_SM + eZHZga * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    if (FlagQuadraticTerms) {
        return ( muVH(sqrt_s)*BrHZgaRatio() * (1.0 + eVHZga ) * (1.0 + eHwidth)/(1.0 + eVHtot)/(1.0 + eHZgaint + eHZgapar) );
    } else {
        return ( muVH(sqrt_s) + BrHZgaRatio() - 1.0 + eVHZga - eVHtot - eHZgaint - eHZgapar + eHwidth );
    }   
}

double NPSMEFTd6::muTHUttHZga(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muttH(sqrt_s)*BrHZgaRatio() * (1.0 + ettHZga ) * (1.0 + eHwidth)/(1.0 + eeettHint + eeettHpar)/(1.0 + eHZgaint + eHZgapar) );
    } else {
        return ( muttH(sqrt_s) + BrHZgaRatio() - 1.0 + ettHZga - eeettHint - eeettHpar - eHZgaint - eHZgapar + eHwidth );
    }    
}

double NPSMEFTd6::muTHUggHZZ(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muggH(sqrt_s)*BrHZZRatio() * (1.0 + eggFHZZ ) * (1.0 + eHwidth)/(1.0 + eggFint + eggFpar)/(1.0 + eHZZint + eHZZpar) );
    } else {
        return ( muggH(sqrt_s) + BrHZZRatio()  - 1.0 + eggFHZZ - eggFint - eggFpar - eHZZint - eHZZpar + eHwidth );
    }   
}

double NPSMEFTd6::muTHUVBFHZZ(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muVBF(sqrt_s)*BrHZZRatio() * (1.0 + eVBFHZZ ) * (1.0 + eHwidth)/(1.0 + eVBFint + eVBFpar)/(1.0 + eHZZint + eHZZpar) );
    } else {
        return ( muVBF(sqrt_s) + BrHZZRatio() - 1.0 + eVBFHZZ - eVBFint - eVBFpar - eHZZint - eHZZpar + eHwidth );
    }    
}

double NPSMEFTd6::muTHUZHZZ(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muZH(sqrt_s)*BrHZZRatio() * (1.0 + eZHZZ ) * (1.0 + eHwidth)/(1.0 + eZHint + eZHpar)/(1.0 + eHZZint + eHZZpar) );
    } else {
        return ( muZH(sqrt_s) + BrHZZRatio() - 1.0 + eZHZZ - eZHint - eZHpar - eHZZint - eHZZpar + eHwidth );
    }   
}

double NPSMEFTd6::muTHUWHZZ(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muWH(sqrt_s)*BrHZZRatio() * (1.0 + eWHZZ ) * (1.0 + eHwidth)/(1.0 + eWHint + eWHpar)/(1.0 + eHZZint + eHZZpar) );
    } else {
        return ( muWH(sqrt_s) + BrHZZRatio() - 1.0 + eWHZZ - eWHint - eWHpar - eHZZint - eHZZpar + eHwidth );
    }   
}

double NPSMEFTd6::muTHUVHZZ(const double sqrt_s) const
{
    //  Theory uncertainty in VH production, from the WH and ZH ones
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);    
    double eVHtot,eVHZZ;
    
    eVHtot = ((eWHint + eWHpar) * sigmaWH_SM + (eZHint + eZHpar) * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    eVHZZ = (eWHZZ * sigmaWH_SM + eZHZZ * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    if (FlagQuadraticTerms) {
        return ( muVH(sqrt_s)*BrHZZRatio() * (1.0 + eVHZZ ) * (1.0 + eHwidth)/(1.0 + eVHtot)/(1.0 + eHZZint + eHZZpar) );
    } else {
        return ( muVH(sqrt_s) + BrHZZRatio() - 1.0 + eVHZZ - eVHtot - eHZZint - eHZZpar + eHwidth );
    }       
}

double NPSMEFTd6::muTHUttHZZ(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muttH(sqrt_s)*BrHZZRatio() * (1.0 + ettHZZ ) * (1.0 + eHwidth)/(1.0 + eeettHint + eeettHpar)/(1.0 + eHZZint + eHZZpar) );
    } else {
        return ( muttH(sqrt_s) + BrHZZRatio() - 1.0 + ettHZZ - eeettHint - eeettHpar - eHZZint - eHZZpar + eHwidth );
    }    
}

double NPSMEFTd6::muTHUggHZZ4l(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muggH(sqrt_s)*BrHZZ4lRatio() * (1.0 + eggFHZZ ) * (1.0 + eHwidth)/(1.0 + eggFint + eggFpar)/(1.0 + eHZZint + eHZZpar) );
    } else {
        return ( muggH(sqrt_s) + BrHZZ4lRatio() - 1.0 + eggFHZZ - eggFint - eggFpar - eHZZint - eHZZpar + eHwidth );
    }  
}

double NPSMEFTd6::muTHUVBFHZZ4l(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muVBF(sqrt_s)*BrHZZ4lRatio() * (1.0 + eVBFHZZ ) * (1.0 + eHwidth)/(1.0 + eVBFint + eVBFpar)/(1.0 + eHZZint + eHZZpar) );
    } else {
        return ( muVBF(sqrt_s) + BrHZZ4lRatio() - 1.0 + eVBFHZZ - eVBFint - eVBFpar - eHZZint - eHZZpar + eHwidth );
    }  
}

double NPSMEFTd6::muTHUZHZZ4l(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muZH(sqrt_s)*BrHZZ4lRatio() * (1.0 + eZHZZ ) * (1.0 + eHwidth)/(1.0 + eZHint + eZHpar)/(1.0 + eHZZint + eHZZpar) );
    } else {
        return ( muZH(sqrt_s) + BrHZZ4lRatio() - 1.0 + eZHZZ - eZHint - eZHpar - eHZZint - eHZZpar + eHwidth );
    }   
}

double NPSMEFTd6::muTHUWHZZ4l(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muWH(sqrt_s)*BrHZZ4lRatio() * (1.0 + eWHZZ ) * (1.0 + eHwidth)/(1.0 + eWHint + eWHpar)/(1.0 + eHZZint + eHZZpar) );
    } else {
        return ( muWH(sqrt_s) + BrHZZ4lRatio() - 1.0 + eWHZZ - eWHint - eWHpar - eHZZint - eHZZpar + eHwidth );
    }  
}

double NPSMEFTd6::muTHUVHZZ4l(const double sqrt_s) const
{
    //  Theory uncertainty in VH production, from the WH and ZH ones
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);    
    double eVHtot,eVHZZ;
    
    eVHtot = ((eWHint + eWHpar) * sigmaWH_SM + (eZHint + eZHpar) * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    eVHZZ = (eWHZZ * sigmaWH_SM + eZHZZ * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    if (FlagQuadraticTerms) {
        return ( muVH(sqrt_s)*BrHZZ4lRatio() * (1.0 + eVHZZ ) * (1.0 + eHwidth)/(1.0 + eVHtot)/(1.0 + eHZZint + eHZZpar) );
    } else {
        return ( muVH(sqrt_s) + BrHZZ4lRatio() - 1.0 + eVHZZ - eVHtot - eHZZint - eHZZpar + eHwidth );
    }  
}

double NPSMEFTd6::muTHUttHZZ4l(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muttH(sqrt_s)*BrHZZ4lRatio() * (1.0 + ettHZZ ) * (1.0 + eHwidth)/(1.0 + eeettHint + eeettHpar)/(1.0 + eHZZint + eHZZpar) );
    } else {
        return ( muttH(sqrt_s) + BrHZZ4lRatio() - 1.0 + ettHZZ - eeettHint - eeettHpar - eHZZint - eHZZpar + eHwidth );
    }  
}

double NPSMEFTd6::muTHUggHWW(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muggH(sqrt_s)*BrHWWRatio() * (1.0 + eggFHWW ) * (1.0 + eHwidth)/(1.0 + eggFint + eggFpar)/(1.0 + eHWWint + eHWWpar) );
    } else {
        return ( muggH(sqrt_s) + BrHWWRatio() - 1.0 + eggFHWW - eggFint - eggFpar - eHWWint - eHWWpar + eHwidth );
    }   
}

double NPSMEFTd6::muTHUVBFHWW(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muVBF(sqrt_s)*BrHWWRatio() * (1.0 + eVBFHWW ) * (1.0 + eHwidth)/(1.0 + eVBFint + eVBFpar)/(1.0 + eHWWint + eHWWpar) );
    } else {
        return ( muVBF(sqrt_s) + BrHWWRatio() - 1.0 + eVBFHWW - eVBFint - eVBFpar - eHWWint - eHWWpar + eHwidth );
    }   
}

double NPSMEFTd6::muTHUZHWW(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muZH(sqrt_s)*BrHWWRatio() * (1.0 + eZHWW ) * (1.0 + eHwidth)/(1.0 + eZHint + eZHpar)/(1.0 + eHWWint + eHWWpar) );
    } else {
        return ( muZH(sqrt_s) + BrHWWRatio() - 1.0 + eZHWW - eZHint - eZHpar - eHWWint - eHWWpar + eHwidth );
    }   
}

double NPSMEFTd6::muTHUWHWW(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muWH(sqrt_s)*BrHWWRatio() * (1.0 + eWHWW ) * (1.0 + eHwidth)/(1.0 + eWHint + eWHpar)/(1.0 + eHWWint + eHWWpar) );
    } else {
        return ( muWH(sqrt_s) + BrHWWRatio() - 1.0 + eWHWW - eWHint - eWHpar - eHWWint - eHWWpar + eHwidth );
    }   
}

double NPSMEFTd6::muTHUVHWW(const double sqrt_s) const
{
    //  Theory uncertainty in VH production, from the WH and ZH ones
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);    
    double eVHtot,eVHWW;
    
    eVHtot = ((eWHint + eWHpar) * sigmaWH_SM + (eZHint + eZHpar) * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    eVHWW = (eWHWW * sigmaWH_SM + eZHWW * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    if (FlagQuadraticTerms) {
        return ( muVH(sqrt_s)*BrHWWRatio() * (1.0 + eVHWW ) * (1.0 + eHwidth)/(1.0 + eVHtot)/(1.0 + eHWWint + eHWWpar) );
    } else {
        return ( muVH(sqrt_s) + BrHWWRatio() - 1.0 + eVHWW - eVHtot - eHWWint - eHWWpar + eHwidth );
    }    
}

double NPSMEFTd6::muTHUttHWW(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muttH(sqrt_s)*BrHWWRatio() * (1.0 + ettHWW ) * (1.0 + eHwidth)/(1.0 + eeettHint + eeettHpar)/(1.0 + eHWWint + eHWWpar) );
    } else {
        return ( muttH(sqrt_s) + BrHWWRatio() - 1.0 + ettHWW - eeettHint - eeettHpar - eHWWint - eHWWpar + eHwidth );
    }   
}

double NPSMEFTd6::muTHUggHWW2l2v(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muggH(sqrt_s)*BrHWW2l2vRatio() * (1.0 + eggFHWW ) * (1.0 + eHwidth)/(1.0 + eggFint + eggFpar)/(1.0 + eHWWint + eHWWpar) );
    } else {
        return ( muggH(sqrt_s) + BrHWW2l2vRatio() - 1.0 + eggFHWW - eggFint - eggFpar - eHWWint - eHWWpar + eHwidth );
    }      
}

double NPSMEFTd6::muTHUVBFHWW2l2v(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muVBF(sqrt_s)*BrHWW2l2vRatio() * (1.0 + eVBFHWW ) * (1.0 + eHwidth)/(1.0 + eVBFint + eVBFpar)/(1.0 + eHWWint + eHWWpar) );
    } else {
        return ( muVBF(sqrt_s) + BrHWW2l2vRatio() - 1.0 + eVBFHWW - eVBFint - eVBFpar - eHWWint - eHWWpar + eHwidth );
    }   
}

double NPSMEFTd6::muTHUZHWW2l2v(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muZH(sqrt_s)*BrHWW2l2vRatio() * (1.0 + eZHWW ) * (1.0 + eHwidth)/(1.0 + eZHint + eZHpar)/(1.0 + eHWWint + eHWWpar) );
    } else {
        return ( muZH(sqrt_s) + BrHWW2l2vRatio() - 1.0 + eZHWW - eZHint - eZHpar - eHWWint - eHWWpar + eHwidth );
    }       
}

double NPSMEFTd6::muTHUWHWW2l2v(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muWH(sqrt_s)*BrHWW2l2vRatio() * (1.0 + eWHWW ) * (1.0 + eHwidth)/(1.0 + eWHint + eWHpar)/(1.0 + eHWWint + eHWWpar) );
    } else {
        return ( muWH(sqrt_s) + BrHWW2l2vRatio() - 1.0 + eWHWW - eWHint - eWHpar - eHWWint - eHWWpar + eHwidth );
    }        
}

double NPSMEFTd6::muTHUVHWW2l2v(const double sqrt_s) const
{
    //  Theory uncertainty in VH production, from the WH and ZH ones
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);    
    double eVHtot,eVHWW;
    
    eVHtot = ((eWHint + eWHpar) * sigmaWH_SM + (eZHint + eZHpar) * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    eVHWW = (eWHWW * sigmaWH_SM + eZHWW * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    if (FlagQuadraticTerms) {
        return ( muVH(sqrt_s)*BrHWW2l2vRatio() * (1.0 + eVHWW ) * (1.0 + eHwidth)/(1.0 + eVHtot)/(1.0 + eHWWint + eHWWpar) );
    } else {
        return ( muVH(sqrt_s) + BrHWW2l2vRatio() - 1.0 + eVHWW - eVHtot - eHWWint - eHWWpar + eHwidth );
    }       
}

double NPSMEFTd6::muTHUttHWW2l2v(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muttH(sqrt_s)*BrHWW2l2vRatio() * (1.0 + ettHWW ) * (1.0 + eHwidth)/(1.0 + eeettHint + eeettHpar)/(1.0 + eHWWint + eHWWpar) );
    } else {
        return ( muttH(sqrt_s) + BrHWW2l2vRatio() - 1.0 + ettHWW - eeettHint - eeettHpar - eHWWint - eHWWpar + eHwidth );
    }        
}

double NPSMEFTd6::muTHUggHmumu(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muggH(sqrt_s)*BrHmumuRatio() * (1.0 + eggFHmumu ) * (1.0 + eHwidth)/(1.0 + eggFint + eggFpar)/(1.0 + eHmumuint + eHmumupar) );
    } else {
        return ( muggH(sqrt_s) + BrHmumuRatio() - 1.0 + eggFHmumu - eggFint - eggFpar - eHmumuint - eHmumupar + eHwidth );
    }       
}

double NPSMEFTd6::muTHUVBFHmumu(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muVBF(sqrt_s)*BrHmumuRatio() * (1.0 + eVBFHmumu ) * (1.0 + eHwidth)/(1.0 + eVBFint + eVBFpar)/(1.0 + eHmumuint + eHmumupar) );
    } else {
        return ( muVBF(sqrt_s) + BrHmumuRatio() - 1.0 + eVBFHmumu - eVBFint - eVBFpar - eHmumuint - eHmumupar + eHwidth );
    }       
}

double NPSMEFTd6::muTHUZHmumu(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muZH(sqrt_s)*BrHmumuRatio() * (1.0 + eZHmumu ) * (1.0 + eHwidth)/(1.0 + eZHint + eZHpar)/(1.0 + eHmumuint + eHmumupar) );
    } else {
        return ( muZH(sqrt_s) + BrHmumuRatio() - 1.0 + eZHmumu - eZHint - eZHpar - eHmumuint - eHmumupar + eHwidth );
    }      
}

double NPSMEFTd6::muTHUWHmumu(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muWH(sqrt_s)*BrHmumuRatio() * (1.0 + eWHmumu ) * (1.0 + eHwidth)/(1.0 + eWHint + eWHpar)/(1.0 + eHmumuint + eHmumupar) );
    } else {
        return ( muWH(sqrt_s) + BrHmumuRatio() - 1.0 + eWHmumu - eWHint - eWHpar - eHmumuint - eHmumupar + eHwidth );
    }      
}

double NPSMEFTd6::muTHUVHmumu(const double sqrt_s) const
{
    //  Theory uncertainty in VH production, from the WH and ZH ones
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);    
    double eVHtot,eVHmumu;
    
    eVHtot = ((eWHint + eWHpar) * sigmaWH_SM + (eZHint + eZHpar) * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    eVHmumu = (eWHmumu * sigmaWH_SM + eZHmumu * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    if (FlagQuadraticTerms) {
        return ( muVH(sqrt_s)*BrHmumuRatio() * (1.0 + eVHmumu ) * (1.0 + eHwidth)/(1.0 + eVHtot)/(1.0 + eHmumuint + eHmumupar) );
    } else {
        return ( muVH(sqrt_s) + BrHmumuRatio() - 1.0 + eVHmumu - eVHtot - eHmumuint - eHmumupar + eHwidth );
    }       
}

double NPSMEFTd6::muTHUttHmumu(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muttH(sqrt_s)*BrHmumuRatio() * (1.0 + ettHmumu ) * (1.0 + eHwidth)/(1.0 + eeettHint + eeettHpar)/(1.0 + eHmumuint + eHmumupar) );
    } else {
        return ( muttH(sqrt_s) + BrHmumuRatio() - 1.0 + ettHmumu - eeettHint - eeettHpar - eHmumuint - eHmumupar + eHwidth );
    }       
}

double NPSMEFTd6::muTHUggHtautau(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muggH(sqrt_s)*BrHtautauRatio() * (1.0 + eggFHtautau ) * (1.0 + eHwidth)/(1.0 + eggFint + eggFpar)/(1.0 + eHtautauint + eHtautaupar) );
    } else {
        return ( muggH(sqrt_s) + BrHtautauRatio() - 1.0 + eggFHtautau - eggFint - eggFpar - eHtautauint - eHtautaupar + eHwidth );
    }      
}

double NPSMEFTd6::muTHUVBFHtautau(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muVBF(sqrt_s)*BrHtautauRatio() * (1.0 + eVBFHtautau ) * (1.0 + eHwidth)/(1.0 + eVBFint + eVBFpar)/(1.0 + eHtautauint + eHtautaupar) );
    } else {
        return ( muVBF(sqrt_s) + BrHtautauRatio() - 1.0 + eVBFHtautau - eVBFint - eVBFpar - eHtautauint - eHtautaupar + eHwidth );
    }        
}

double NPSMEFTd6::muTHUZHtautau(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muZH(sqrt_s)*BrHtautauRatio() * (1.0 + eZHtautau ) * (1.0 + eHwidth)/(1.0 + eZHint + eZHpar)/(1.0 + eHtautauint + eHtautaupar) );
    } else {
        return ( muZH(sqrt_s) + BrHtautauRatio() - 1.0 + eZHtautau - eZHint - eZHpar - eHtautauint - eHtautaupar + eHwidth );
    }       
}

double NPSMEFTd6::muTHUWHtautau(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muWH(sqrt_s)*BrHtautauRatio() * (1.0 + eWHtautau ) * (1.0 + eHwidth)/(1.0 + eWHint + eWHpar)/(1.0 + eHtautauint + eHtautaupar) );
    } else {
        return ( muWH(sqrt_s) + BrHtautauRatio() - 1.0 + eWHtautau - eWHint - eWHpar - eHtautauint - eHtautaupar + eHwidth );
    }       
}

double NPSMEFTd6::muTHUVHtautau(const double sqrt_s) const
{
    //  Theory uncertainty in VH production, from the WH and ZH ones
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);    
    double eVHtot,eVHtautau;
    
    eVHtot = ((eWHint + eWHpar) * sigmaWH_SM + (eZHint + eZHpar) * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    eVHtautau = (eWHtautau * sigmaWH_SM + eZHtautau * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    if (FlagQuadraticTerms) {
        return ( muVH(sqrt_s)*BrHtautauRatio() * (1.0 + eVHtautau ) * (1.0 + eHwidth)/(1.0 + eVHtot)/(1.0 + eHtautauint + eHtautaupar) );
    } else {
        return ( muVH(sqrt_s) + BrHtautauRatio() - 1.0 + eVHtautau - eVHtot - eHtautauint - eHtautaupar + eHwidth );
    }       
}

double NPSMEFTd6::muTHUttHtautau(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muttH(sqrt_s)*BrHtautauRatio() * (1.0 + ettHtautau ) * (1.0 + eHwidth)/(1.0 + eeettHint + eeettHpar)/(1.0 + eHtautauint + eHtautaupar) );
    } else {
        return ( muttH(sqrt_s) + BrHtautauRatio() - 1.0 + ettHtautau - eeettHint - eeettHpar - eHtautauint - eHtautaupar + eHwidth );
    }        
}

double NPSMEFTd6::muTHUggHbb(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muggH(sqrt_s)*BrHbbRatio() * (1.0 + eggFHbb ) * (1.0 + eHwidth)/(1.0 + eggFint + eggFpar)/(1.0 + eHbbint + eHbbpar) );
    } else {
        return ( muggH(sqrt_s) + BrHbbRatio() - 1.0 + eggFHbb - eggFint - eggFpar - eHbbint - eHbbpar + eHwidth );
    }        
}

double NPSMEFTd6::muTHUVBFHbb(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muVBF(sqrt_s)*BrHbbRatio() * (1.0 + eVBFHbb ) * (1.0 + eHwidth)/(1.0 + eVBFint + eVBFpar)/(1.0 + eHbbint + eHbbpar) );
    } else {
        return ( muVBF(sqrt_s) + BrHbbRatio() - 1.0 + eVBFHbb - eVBFint - eVBFpar - eHbbint - eHbbpar + eHwidth );
    }       
}

double NPSMEFTd6::muTHUZHbb(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muZH(sqrt_s)*BrHbbRatio() * (1.0 + eZHbb ) * (1.0 + eHwidth)/(1.0 + eZHint + eZHpar)/(1.0 + eHbbint + eHbbpar) );
    } else {
        return ( muZH(sqrt_s) + BrHbbRatio() - 1.0 + eZHbb - eZHint - eZHpar - eHbbint - eHbbpar + eHwidth );
    }       
}

double NPSMEFTd6::muTHUWHbb(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muWH(sqrt_s)*BrHbbRatio() * (1.0 + eWHbb ) * (1.0 + eHwidth)/(1.0 + eWHint + eWHpar)/(1.0 + eHbbint + eHbbpar) );
    } else {
        return ( muWH(sqrt_s) + BrHbbRatio() - 1.0 + eWHbb - eWHint - eWHpar - eHbbint - eHbbpar + eHwidth );
    }        
}

double NPSMEFTd6::muTHUVHbb(const double sqrt_s) const
{
    //  Theory uncertainty in VH production, from the WH and ZH ones
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);    
    double eVHtot,eVHbb;
    
    eVHtot = ((eWHint + eWHpar) * sigmaWH_SM + (eZHint + eZHpar) * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    eVHbb = (eWHbb * sigmaWH_SM + eZHbb * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    if (FlagQuadraticTerms) {
        return ( muVH(sqrt_s)*BrHbbRatio() * (1.0 + eVHbb ) * (1.0 + eHwidth)/(1.0 + eVHtot)/(1.0 + eHbbint + eHbbpar) );
    } else {
        return ( muVH(sqrt_s) + BrHbbRatio() - 1.0 + eVHbb - eVHtot - eHbbint - eHbbpar + eHwidth );
    }       
}

double NPSMEFTd6::muTHUttHbb(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muttH(sqrt_s)*BrHbbRatio() * (1.0 + ettHbb ) * (1.0 + eHwidth)/(1.0 + eeettHint + eeettHpar)/(1.0 + eHbbint + eHbbpar) );
    } else {
        return ( muttH(sqrt_s) + BrHbbRatio() - 1.0 + ettHbb - eeettHint - eeettHpar - eHbbint - eHbbpar + eHwidth );
    }        
}

double NPSMEFTd6::muTHUVBFBRinv(const double sqrt_s) const
{
    return ( muVBF(sqrt_s)*Br_H_inv() * (1.0 + eVBFHinv )/(1.0 + eVBFint + eVBFpar) );   
}

double NPSMEFTd6::muTHUVBFHinv(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muVBF(sqrt_s)*BrHtoinvRatio() * (1.0 + eVBFHinv )/(1.0 + eVBFint + eVBFpar) );
    } else {
        return ( muVBF(sqrt_s) + BrHtoinvRatio() - 1.0 + eVBFHinv - eVBFint - eVBFpar );
    }        
}

double NPSMEFTd6::muTHUVHBRinv(const double sqrt_s) const
{
    //  Theory uncertainty in VH production, from the WH and ZH ones
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);    
    double eVHtot;
    
    eVHtot = ((eWHint + eWHpar) * sigmaWH_SM + (eZHint + eZHpar) * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    return ( muVH(sqrt_s)*Br_H_inv() * (1.0 + eVHinv )/(1.0 + eVHtot) );   
}

double NPSMEFTd6::muTHUVHinv(const double sqrt_s) const
{
    //  Theory uncertainty in VH production, from the WH and ZH ones
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);    
    double eVHtot;
    
    eVHtot = ((eWHint + eWHpar) * sigmaWH_SM + (eZHint + eZHpar) * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    if (FlagQuadraticTerms) {
        return ( muVH(sqrt_s)*BrHtoinvRatio() * (1.0 + eVHinv )/(1.0 + eVHtot) );
    } else {
        return ( muVH(sqrt_s) + BrHtoinvRatio() - 1.0 + eVHinv - eVHtot );
    }        
}


double NPSMEFTd6::muTHUggHZZ4mu(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muggH(sqrt_s)*BrHZZ4muRatio() * (1.0 + eggFHZZ ) * (1.0 + eHwidth)/(1.0 + eggFint + eggFpar)/(1.0 + eHZZint + eHZZpar) );
    } else {
        return ( muggH(sqrt_s) + BrHZZ4muRatio() - 1.0 + eggFHZZ - eggFint - eggFpar - eHZZint - eHZZpar + eHwidth );
    } 
}

double NPSMEFTd6::muTHUggHZgamumu(const double sqrt_s) const
{
    if (FlagQuadraticTerms) {
        return ( muggH(sqrt_s)*BrHZgamumuRatio() * (1.0 + eggFHZga ) * (1.0 + eHwidth)/(1.0 + eggFint + eggFpar)/(1.0 + eHZgaint + eHZgapar) );
    } else {
        return ( muggH(sqrt_s) + BrHZgamumuRatio() - 1.0 + eggFHZga - eggFint - eggFpar - eHZgaint - eHZgapar + eHwidth );
    } 
}


///////////////////////////////////////////////////////////////////////////////

double NPSMEFTd6::deltag1ZNP() const
{
      double NPdirect, NPindirect;
      
      /*    From own calculations. Agrees with with LHCHXWG-INT-2015-001 for common interactions */
      NPdirect = sW_tree / sqrt( 4.0 * M_PI * aleMz );
      NPdirect = - NPdirect * (Mz * Mz / v () / v() ) * CiDHW * v2_over_LambdaNP2;
      
      NPindirect = - 1.0 / (cW2_tree-sW2_tree);
      
      NPindirect = NPindirect * (sW_tree * CiHWB / cW_tree 
              + 0.25 * CiHD ) * v2_over_LambdaNP2
              + 0.5 * NPindirect * DeltaGF() ;
      
      return NPdirect + NPindirect + dg1Z ;
}
      
double NPSMEFTd6::deltaKgammaNP() const
{
      double NPdirect;

      /*    Translate from LHCHXWG-INT-2015-001: Checked with own calculations  */
      NPdirect = sqrt( 4.0 * M_PI * aleMz ) / 4.0 / sW2_tree;
      
      NPdirect = NPdirect * ( (4.0 * sW_tree * cW_tree / sqrt( 4.0 * M_PI * aleMz ) ) * CiHWB 
              - sW_tree * CiDHW 
              - cW_tree * CiDHB ) * v2_over_LambdaNP2;
      
      return NPdirect + dKappaga ;
}
      
double NPSMEFTd6::lambdaZNP() const
{
      double NPdirect;

      /*    Translate from LHCHXWG-INT-2015-001: Checked with own calculations  */
      NPdirect = - (3.0 / 2.0) * (sqrt( 4.0 * M_PI * aleMz ) / sW_tree) * CiW * v2_over_LambdaNP2;

      return NPdirect + lambZ ;
}

///////////////////////////////////////////////////////////////////////////////

double NPSMEFTd6::deltag1ZNPEff() const
{
    /* From arXiv:1708.09079 [hep-ph]. In our case, delta_e=0 since it is taken as inputs and its effects propagated
     * everywhere else */
    double dgEff;
    
    dgEff = (1.0/ cW2_tree) * ( (cW2_tree - sW2_tree)*deltaGL_f(leptons[ELECTRON])/gZlL +
            sW2_tree * deltaGR_f(leptons[ELECTRON])/gZlR -
            2.0 * deltaGL_Wff(leptons[NEUTRINO_1], leptons[ELECTRON]).real()/UevL );
      
    return dgEff + deltag1ZNP() ;
}
      
double NPSMEFTd6::deltaKgammaNPEff() const
{     
    /* From arXiv:1708.09079 [hep-ph]. In our case, delta_e=0 since it is taken as inputs and its effects propagated
     * everywhere else */
    double dgEff;
    
    dgEff = (cW2_tree - sW2_tree)*( deltaGL_f(leptons[ELECTRON])/gZlL - deltaGR_f(leptons[ELECTRON])/gZlR )
            - 2.0 * deltaGL_Wff(leptons[NEUTRINO_1], leptons[ELECTRON]).real()/UevL;
    
    return dgEff + deltaKgammaNP() ;
}

///////////////////////////////////////////////////////////////////////////////

double NPSMEFTd6::dxseeWWdcos(const double sqrt_s, const double cos) const
{
    double sqrt_sGeV = 1000. * sqrt_s;
    double s = sqrt_sGeV * sqrt_sGeV;
    double cos2 = cos * cos;
    double sin2 = 1.0 - cos2;
    double sin = sqrt(sin2);
    
    double topb = 0.3894*1000000000.0;
    
//  NC and CC couplings
    double gLe, gRe;
    gslpp::complex Uenu;
    
    gLe = -0.5 + sW2_tree + deltaGL_f(leptons[ELECTRON]);
    gRe = sW2_tree + deltaGR_f(leptons[ELECTRON]);
    
    Uenu = deltaGL_Wff(leptons[NEUTRINO_1], leptons[ELECTRON]);
    Uenu = 1.0 + Uenu;
    
//  W mass
    double mw;
    
    mw = Mw();

//  Wigner functions
    double d1pp[2],d1mm[2],d1p0[2],d1m0[2],d10p[2],d10m[2],d100[2];
    
    d1pp[0]=sqrt((1.0 - cos2)/2.0);
    d1pp[1]=-sqrt((1.0 - cos2)/2.0);
    
    d1mm[0]=d1pp[0];
    d1mm[1]=d1pp[1];
    
    d1p0[0]=(1.0 - cos)/2.0;
    d1p0[1]=(1.0 + cos)/2.0;
    
    d1m0[0]=d1p0[1];
    d1m0[1]=d1p0[0];
    
    d10p[0]=d1p0[1];
    d10p[1]=d1p0[0];
    
    d10m[0]=d1p0[0];
    d10m[1]=d1p0[1];
    
    d100[0]=d1pp[0];
    d100[1]=d1pp[1];
    
    gslpp::matrix<double> d1LH(3, 3, 0.0);
    
    gslpp::matrix<double> d1RH(3, 3, 0.0);
    
    d1LH.assign(0,0, d1pp[0]);
    d1LH.assign(0,1, d1p0[0]);
    d1LH.assign(0,2, 0.0);
    
    d1LH.assign(1,0, d10p[0]);
    d1LH.assign(1,1, d100[0]);
    d1LH.assign(1,2, d10m[0]);
    
    d1LH.assign(2,0, 0.0);
    d1LH.assign(2,1, d1m0[0]);
    d1LH.assign(2,2, d1mm[0]);
    
    d1RH.assign(0,0, d1pp[1]);
    d1RH.assign(0,1, d1p0[1]);
    d1RH.assign(0,2, 0.0);
    
    d1RH.assign(1,0, d10p[1]);
    d1RH.assign(1,1, d100[1]);
    d1RH.assign(1,2, d10m[1]);
    
    d1RH.assign(2,0, 0.0);
    d1RH.assign(2,1, d1m0[1]);
    d1RH.assign(2,2, d1mm[1]);
    
//  TGC parameterization
    double g1Z,g1ga,kZ,kga,lambdaZ,lambdaga,g4Z,g4ga,g5Z,g5ga,ktZ,ktga,lambdatZ,lambdatga;

//  TGC present in the SM     
    g1Z=1.0 + deltag1ZNP();
    g1ga=1.0;
    kZ=1.0 + deltag1ZNP() - (sW2_tree/cW2_tree) * deltaKgammaNP();
    kga=1.0 + deltaKgammaNP();
//  TGC not present in the SM
    lambdaZ=lambdaZNP(); //Check normalization
    lambdaga=lambdaZ;
    g4Z=0.0;
    g4ga=0.0;
    g5Z=0.0;
    g5ga=0.0;
    ktZ=0.0;
    ktga=0.0;
    lambdatZ=0.0;
    lambdatga=0.0;
    
    double f3Z, f3ga;
    
    f3Z = g1Z + kZ + lambdaZ;    
    f3ga = g1ga + kga + lambdaga;
    
 // Kinematic factors
    double beta, gamma, gamma2;
    
    beta = sqrt(1.0 - 4.0 * mw * mw / s);
    gamma = sqrt_sGeV/(2.0 * mw);
    gamma2= gamma*gamma;
    
//  J=1 Subamplitudes: Z
    gslpp::complex AZpp, AZmm, AZp0, AZm0, AZ0p, AZ0m, AZ00;
    
    AZpp = gslpp::complex(g1Z + 2.0* gamma2* lambdaZ, (ktZ + lambdatZ - 2.0*lambdatZ)/beta , false);
    AZmm = gslpp::complex(g1Z + 2.0* gamma2* lambdaZ, -(ktZ + lambdatZ - 2.0*lambdatZ)/beta , false);
    AZp0 = gslpp::complex(f3Z + beta * g5Z , -g4Z + (ktZ-lambdatZ)/beta , false);
    AZp0 = gamma * AZp0;
    AZm0 = gslpp::complex(f3Z - beta * g5Z , -g4Z - (ktZ-lambdatZ)/beta , false);
    AZm0 = gamma * AZm0;
    AZ0p = gslpp::complex(f3Z - beta * g5Z , g4Z + (ktZ-lambdatZ)/beta , false);
    AZ0p = gamma * AZ0p;
    AZ0m = gslpp::complex(f3Z + beta * g5Z , g4Z - (ktZ-lambdatZ)/beta , false);
    AZ0m = gamma * AZ0m;
    AZ00 = gslpp::complex( g1Z + 2.0*gamma2*kZ, 0.0 , false);
    
//  Collect in matrices and separate LH and RH
    gslpp::matrix<gslpp::complex> AmpZLH(3, 3, 0.0);
    gslpp::matrix<gslpp::complex> AmpZRH(3, 3, 0.0);
    
    AmpZLH.assign(0,0, AZpp *  d1LH(0,0) );
    AmpZLH.assign(0,1, AZp0 *  d1LH(0,1));
    AmpZLH.assign(0,2, 0.0);
    
    AmpZLH.assign(1,0, AZ0p *  d1LH(1,0));
    AmpZLH.assign(1,1, AZ00 *  d1LH(1,1));
    AmpZLH.assign(1,2, AZ0m *  d1LH(1,2));
    
    AmpZLH.assign(2,0, 0.0);
    AmpZLH.assign(2,1, AZm0 *  d1LH(2,1));
    AmpZLH.assign(2,2, AZmm *  d1LH(2,2));
    
    AmpZLH = AmpZLH * beta * s/(s-Mz*Mz);
    
//  Add the correct Zff coupling
    AmpZLH = AmpZLH * gLe / sW2_tree;
    
    AmpZRH.assign(0,0, AZpp *  d1RH(0,0) );
    AmpZRH.assign(0,1, AZp0 *  d1RH(0,1));
    AmpZRH.assign(0,2, 0.0);
    
    AmpZRH.assign(1,0, AZ0p *  d1RH(1,0));
    AmpZRH.assign(1,1, AZ00 *  d1RH(1,1));
    AmpZRH.assign(1,2, AZ0m *  d1RH(1,2));
    
    AmpZRH.assign(2,0, 0.0);
    AmpZRH.assign(2,1, AZm0 *  d1RH(2,1));
    AmpZRH.assign(2,2, AZmm *  d1RH(2,2));
    
    AmpZRH = AmpZRH * beta * s/(s-Mz*Mz);

//  Add the correct Zff coupling    
    AmpZRH = AmpZRH * gRe / sW2_tree;

//  J=1 Subamplitudes: gamma
    gslpp::complex Agapp, Agamm, Agap0, Agam0, Aga0p, Aga0m, Aga00;
    
    Agapp = gslpp::complex(g1ga + 2.0* gamma2* lambdaga, (ktga + lambdatga - 2.0*lambdatga)/beta , false);
    Agamm = gslpp::complex(g1ga + 2.0* gamma2* lambdaga, -(ktga + lambdatga - 2.0*lambdatga)/beta , false);
    Agap0 = gslpp::complex(f3ga + beta * g5ga , -g4ga + (ktga-lambdatga)/beta , false);
    Agap0 = gamma * Agap0;
    Agam0 = gslpp::complex(f3ga - beta * g5ga , -g4ga - (ktga-lambdatga)/beta , false);
    Agam0 = gamma * Agam0;
    Aga0p = gslpp::complex(f3ga - beta * g5ga , g4ga + (ktga-lambdatga)/beta , false);
    Aga0p = gamma * Aga0p;
    Aga0m = gslpp::complex(f3ga + beta * g5ga , g4ga - (ktga-lambdatga)/beta , false);
    Aga0m = gamma * Aga0m;
    Aga00 = gslpp::complex( g1ga + 2.0*gamma2*kga, 0.0 , false);

//  Collect in matrices. Here LH = RH, except for the Wigner functions  
    gslpp::matrix<gslpp::complex> AmpgaLH(3, 3, 0.0);
    gslpp::matrix<gslpp::complex> AmpgaRH(3, 3, 0.0);
    
    AmpgaLH.assign(0,0, Agapp * d1LH(0,0));
    AmpgaLH.assign(0,1, Agap0 * d1LH(0,1));
    AmpgaLH.assign(0,2, 0.0);
    
    AmpgaLH.assign(1,0, Aga0p * d1LH(1,0));
    AmpgaLH.assign(1,1, Aga00 * d1LH(1,1));
    AmpgaLH.assign(1,2, Aga0m * d1LH(1,2));
    
    AmpgaLH.assign(2,0, 0.0);
    AmpgaLH.assign(2,1, Agam0 * d1LH(2,1));
    AmpgaLH.assign(2,2, Agamm * d1LH(2,2));
    
    AmpgaRH.assign(0,0, Agapp * d1RH(0,0));
    AmpgaRH.assign(0,1, Agap0 * d1RH(0,1));
    AmpgaRH.assign(0,2, 0.0);
    
    AmpgaRH.assign(1,0, Aga0p * d1RH(1,0));
    AmpgaRH.assign(1,1, Aga00 * d1RH(1,1));
    AmpgaRH.assign(1,2, Aga0m * d1RH(1,2));
    
    AmpgaRH.assign(2,0, 0.0);
    AmpgaRH.assign(2,1, Agam0 * d1RH(2,1));
    AmpgaRH.assign(2,2, Agamm * d1RH(2,2));
    
    AmpgaLH = -beta * AmpgaLH;
    AmpgaRH = -beta * AmpgaRH;
    
//  J=1 Subamplitudes: neutrino
    gslpp::complex Bpp, Bmm, Bp0, Bm0, B0p, B0m, B00;
    gslpp::complex Cpp, Cmm, Cp0, Cm0, C0p, C0m, C00;
    
    Bpp = gslpp::complex(1.0 , 0.0 , false);
    Bmm = Bpp;
    Bp0 = gslpp::complex( 2.0 * gamma, 0.0 , false);
    Bm0 = Bp0;
    B0p = Bp0;
    B0m = Bp0;
    B00 = gslpp::complex( 2.0 * gamma2, 0.0 , false);
    
    Cpp = gslpp::complex(1.0/gamma2 , 0.0 , false);
    Cmm = Cpp;
    Cp0 = gslpp::complex( 2.0 * (1.0 + beta)/gamma, 0.0 , false);
    Cm0 = gslpp::complex( 2.0 * (1.0 - beta)/gamma, 0.0 , false);
    C0p = Cm0;
    C0m = Cp0;
    C00 = gslpp::complex( 2.0 / gamma2, 0.0 , false);
    
//  Collect in matrices. Here LH = RH    
    gslpp::matrix<gslpp::complex> Bnu(3, 3, 0.0);
    gslpp::matrix<gslpp::complex> Cnu(3, 3, 0.0);
    
    Bnu.assign(0,0, Bpp * d1LH(0,0));
    Bnu.assign(0,1, Bp0 * d1LH(0,1));
    Bnu.assign(0,2, 0.0);
    
    Bnu.assign(1,0, B0p * d1LH(1,0));
    Bnu.assign(1,1, B00 * d1LH(1,1));
    Bnu.assign(1,2, B0m * d1LH(1,2));
    
    Bnu.assign(2,0, 0.0);
    Bnu.assign(2,1, Bm0 * d1LH(2,1));
    Bnu.assign(2,2, Bmm * d1LH(2,2));
    
    Cnu.assign(0,0, Cpp * d1LH(0,0));
    Cnu.assign(0,1, Cp0 * d1LH(0,1));
    Cnu.assign(0,2, 0.0);
    
    Cnu.assign(1,0, C0p * d1LH(1,0));
    Cnu.assign(1,1, C00 * d1LH(1,1));
    Cnu.assign(1,2, C0m * d1LH(1,2));
    
    Cnu.assign(2,0, 0.0);
    Cnu.assign(2,1, Cm0 * d1LH(2,1));
    Cnu.assign(2,2, Cmm * d1LH(2,2));
 
//  The matrix with the total J=1 neutrino amplitude (only LH neutrinos)  
    gslpp::matrix<gslpp::complex> Ampnu1(3, 3, 0.0);
    
    Ampnu1 = Bnu - Cnu/(1.0 + beta*beta - 2.0 * beta * cos);
    
    Ampnu1 = Uenu * Uenu.conjugate() * Ampnu1 / (2.0 * beta * sW2_tree);
    
    gslpp::matrix<gslpp::complex> Ampnu2(3, 3, 0.0);

    Ampnu2.assign(0,2, (1.0 - cos)/2.0 );
    Ampnu2.assign(1,1, 0.0);
    Ampnu2.assign(2,0, -(1.0 + cos)/2.0);
    
    Ampnu2 = (8.0 * M_PI * aleMz / sW2_tree)* Uenu * Uenu.conjugate() * Ampnu2 * sin / (1.0 + beta*beta - 2.0*beta*cos);
    
//  Total amplitudes 
    gslpp::matrix<gslpp::complex> MRH(3, 3, 0.0);
    gslpp::matrix<gslpp::complex> MLH(3, 3, 0.0);

    MRH = sqrt(2.0) * 4.0 * M_PI * aleMz * (AmpZRH + AmpgaRH);
    MLH = - sqrt(2.0) * 4.0 * M_PI * aleMz * (AmpZLH + AmpgaLH + Ampnu1) + Ampnu2;
    
//  Total amplitude squared and differential cross section (in pb)
    gslpp::matrix<double> M2(3, 3, 0.0);
    double dxsdcos;
    
    dxsdcos = 0.0;
    
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            M2.assign(i,j, (MRH(i,j)* (MRH(i,j).conjugate()) 
                    + MLH(i,j)* (MLH(i,j).conjugate())).real() );
            
            dxsdcos = dxsdcos + M2(i,j);
        }
    }
    
//  Differential cross section in pb
    dxsdcos = (topb * beta / 32.0 / M_PI / s) * dxsdcos;

    return dxsdcos;
}

double NPSMEFTd6::dxseeWWdcosBin(const double sqrt_s, const double cos1, const double cos2) const
{    
    double xsWWbin;/**< Gsl integral variable */
    double errWW;/**< Gsl integral variable */
    
    gsl_function FR;/**< Gsl integral variable */
    
    FR = convertToGslFunction(boost::bind(&NPSMEFTd6::dxseeWWdcos,&(*this), sqrt_s, _1));
    
    gsl_integration_cquad(&FR, cos1, cos2, 1.e-5, 1.e-4, w_WW, &xsWWbin, &errWW, NULL);
    
//  Simple integration for testing
//    double cosx;
    
//    xsWWbin = 0.0;
    
//    for (int i=1; i<100; i++){
//        cosx = cos1 +  i*(cos2-cos1)/100;
//        xsWWbin = xsWWbin + dxseeWWdcos(sqrt_s, cosx);
//    }
    
//    xsWWbin = xsWWbin + 0.5 * (dxseeWWdcos(sqrt_s, cos1) + dxseeWWdcos(sqrt_s, cos2));
    
//    xsWWbin = xsWWbin * (cos2-cos1)/100;
    
//  Compute the BR into e nu, mu nu for one W and into jets for the other
    double BRlv, BRjj;
    
    BRlv = GammaW(leptons[NEUTRINO_1], leptons[ELECTRON]) + 
            GammaW(leptons[NEUTRINO_2], leptons[MU]) + 
            GammaW(leptons[NEUTRINO_3], leptons[TAU]);
    
    BRjj = GammaW() - BRlv;
    
    BRlv = BRlv - GammaW(leptons[NEUTRINO_3], leptons[TAU]);
    
    BRlv =BRlv / GammaW();
    
    BRjj =BRjj / GammaW();
    
    
    
    return xsWWbin * BRlv * BRjj;
}

double NPSMEFTd6::xseeWW(const double sqrt_s) const
{    
    return dxseeWWdcosBin(sqrt_s, -1.0, 1.0);
}


double NPSMEFTd6::mueeWW(const double sqrt_s) const
{
    double mu = 1.0;

     if (sqrt_s == 0.161) {
        
        mu += 
                -127.685 * CiHL1_11 / LambdaNP2
                -175.567 * CiHe_11 / LambdaNP2
                +242506. * CiHL3_11 / LambdaNP2
                -86570.7 * CiHD / LambdaNP2
                -189772. * CiHWB / LambdaNP2
                +12.769 * CiDHB / LambdaNP2
                +6.384 * CiDHW / LambdaNP2
                +0. * CiW / LambdaNP2
                -2.858 * DeltaGF()
                -70.01 * deltaMwd6();
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( -13.134 * deltaMz()
                +0. * deltaaMZ()
                +18.795 * deltaGmu() ); 
        
        if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }

    } else if (sqrt_s == 0.240) {
        
        mu += 
                -26882.4 * CiHL1_11 / LambdaNP2
                -17485.4 * CiHe_11 / LambdaNP2
                +267456. * CiHL3_11 / LambdaNP2
                -83799.2 * CiHD / LambdaNP2
                -168074. * CiHWB / LambdaNP2
                +3199.72 * CiDHB / LambdaNP2
                +3401.93 * CiDHW / LambdaNP2
                +6649.22 * CiW / LambdaNP2
                -2.812 * DeltaGF()
                -0.993 * deltaMwd6();
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( +4.101 * deltaMz()
                -0.584 * deltaaMZ()
                +2.688 * deltaGmu() ); 
        
        if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }

    } else if (sqrt_s == 0.250) {

        mu += 
                -29442.7 * CiHL1_11 / LambdaNP2
                -18494.5 * CiHe_11 / LambdaNP2
                +269747. * CiHL3_11 / LambdaNP2
                -83750.9 * CiHD / LambdaNP2
                -167811. * CiHWB / LambdaNP2
                +3401.99 * CiDHB / LambdaNP2
                +3624.67 * CiDHW / LambdaNP2
                +7249.33 * CiW / LambdaNP2
                -2.812 * DeltaGF()
                -0.959 * deltaMwd6();
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( +4.184 * deltaMz()
                -0.585 * deltaaMZ()
                +2.709 * deltaGmu() ); 
        
        if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 0.350) {

        mu += 
                -47552.4 * CiHL1_11 / LambdaNP2
                -23798.8 * CiHe_11 / LambdaNP2
                +289379. * CiHL3_11 / LambdaNP2
                -83905.3 * CiHD / LambdaNP2
                -168326. * CiHWB / LambdaNP2
                +5979.05 * CiDHB / LambdaNP2
                +6520.95 * CiDHW / LambdaNP2
                +10476.9 * CiW / LambdaNP2
                -2.832 * DeltaGF()
                -0.781 * deltaMwd6();
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( +4.516 * deltaMz()
                -0.659 * deltaaMZ()
                +2.768 * deltaGmu()); 
        
        if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 0.365) {

        mu += 
                -49800.4 * CiHL1_11 / LambdaNP2
                -24520.1 * CiHe_11 / LambdaNP2
                +290743. * CiHL3_11 / LambdaNP2
                -84033.5 * CiHD / LambdaNP2
                -168466. * CiHWB / LambdaNP2
                +6310.59 * CiDHB / LambdaNP2
                +6842.81 * CiDHW / LambdaNP2
                +10606.3 * CiW / LambdaNP2
                -2.828 * DeltaGF()
                -0.775 * deltaMwd6();
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( +4.533 * deltaMz()
                -0.661 * deltaaMZ()
                +2.789 * deltaGmu() ); 
        
        if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 0.500) {

        mu += 
                -68234.1 * CiHL1_11 / LambdaNP2
                -31290. * CiHe_11 / LambdaNP2
                +309504. * CiHL3_11 / LambdaNP2
                -84926.8 * CiHD / LambdaNP2
                -171658. * CiHWB / LambdaNP2
                +9325.19 * CiDHB / LambdaNP2
                +10009.9 * CiDHW / LambdaNP2
                +10896.4 * CiW / LambdaNP2
                -2.84 * DeltaGF()
                -0.705 * deltaMwd6();
        
    // Add modifications due to small variations of the SM parameters    
        mu += cHSM * ( +4.7 * deltaMz()
                -0.683 * deltaaMZ()
                +2.799 * deltaGmu() ); 
        
        if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::mueeWW()");
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}


double NPSMEFTd6::mueeWWPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    double mu = 1.0;

    if (sqrt_s == 0.240) {
                
        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                -23395. * CiHL1_11 / LambdaNP2
                -261092. * CiHe_11 / LambdaNP2
                +231526. * CiHL3_11 / LambdaNP2
                -72645.8 * CiHD / LambdaNP2
                -25084.5 * CiHWB / LambdaNP2
                +27060.4 * CiDHB / LambdaNP2
                -7822.83 * CiDHW / LambdaNP2
                -587.63 * CiW / LambdaNP2
                -2.437 * DeltaGF()
                -1.554 * deltaMwd6(); 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +3.226 * deltaMz()
                -0.083 * deltaaMZ()
                +2.189 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                -27334.5 * CiHL1_11 / LambdaNP2
                -564.392 * CiHe_11 / LambdaNP2
                +269600. * CiHL3_11 / LambdaNP2
                -84684.5 * CiHD / LambdaNP2
                -178168. * CiHWB / LambdaNP2
                +1539.25 * CiDHB / LambdaNP2
                +4130.32 * CiDHW / LambdaNP2
                +7121.6 * CiW / LambdaNP2
                -2.838 * DeltaGF()
                -0.949 * deltaMwd6(); 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.156 * deltaMz()
                -0.607 * deltaaMZ()
                +2.724 * deltaGmu() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeWWPol()");
        }
        
    } else if (sqrt_s == 0.250) {
        
        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                -25554.9 * CiHL1_11 / LambdaNP2
                -274633. * CiHe_11 / LambdaNP2
                +234621. * CiHL3_11 / LambdaNP2
                -72498.3 * CiHD / LambdaNP2
                -23308.5 * CiHWB / LambdaNP2
                +29321.9 * CiDHB / LambdaNP2
                -7518.62 * CiDHW / LambdaNP2
                +314.876 * CiW / LambdaNP2
                -2.444 * DeltaGF()
                -1.448 * deltaMwd6(); 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +3.37 * deltaMz()
                -0.119 * deltaaMZ()
                +2.223 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                -29714.6 * CiHL1_11 / LambdaNP2
                -693.518 * CiHe_11 / LambdaNP2
                +271032. * CiHL3_11 / LambdaNP2
                -84929.3 * CiHD / LambdaNP2
                -177727. * CiHWB / LambdaNP2
                +1648.44 * CiDHB / LambdaNP2
                +4443.85 * CiDHW / LambdaNP2
                +7778.07 * CiW / LambdaNP2
                -2.829 * DeltaGF()
                -0.914 * deltaMwd6(); 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.233 * deltaMz()
                -0.62 * deltaaMZ()
                +2.73 * deltaGmu() );
    
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                -27418.7 * CiHL1_11 / LambdaNP2
                -157891. * CiHe_11 / LambdaNP2
                +250086. * CiHL3_11 / LambdaNP2
                -77904.2 * CiHD / LambdaNP2
                -89451.9 * CiHWB / LambdaNP2
                +17499.7 * CiDHB / LambdaNP2
                -2499.14 * CiDHW / LambdaNP2
                +3435.6 * CiW / LambdaNP2
                -2.607 * DeltaGF()
                -1.242 * deltaMwd6(); 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +3.759 * deltaMz()
                -0.343 * deltaaMZ()
                +2.459 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                -29686. * CiHL1_11 / LambdaNP2
                -1698.32 * CiHe_11 / LambdaNP2
                +271004. * CiHL3_11 / LambdaNP2
                -84881.5 * CiHD / LambdaNP2
                -177249. * CiHWB / LambdaNP2
                +1732.98 * CiDHB / LambdaNP2
                +4380.98 * CiDHW / LambdaNP2
                +7742.96 * CiW / LambdaNP2
                -2.828 * DeltaGF()
                -0.915 * deltaMwd6(); 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.244 * deltaMz()
                -0.624 * deltaaMZ()
                +2.729 * deltaGmu() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeWWPol()");
        }
        
    } else if (sqrt_s == 0.350) {
        
        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                -43312.4 * CiHL1_11 / LambdaNP2
                -370403. * CiHe_11 / LambdaNP2
                +262809. * CiHL3_11 / LambdaNP2
                -76119.5 * CiHD / LambdaNP2
                -35565.5 * CiHWB / LambdaNP2
                +48488.8 * CiDHB / LambdaNP2
                -4519.05 * CiDHW / LambdaNP2
                +6279.71 * CiW / LambdaNP2
                -2.571 * DeltaGF()
                -1.059 * deltaMwd6(); 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.035 * deltaMz()
                -0.336 * deltaaMZ()
                +2.471 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                -47925. * CiHL1_11 / LambdaNP2
                -912.302 * CiHe_11 / LambdaNP2
                +290384. * CiHL3_11 / LambdaNP2
                -84475.3 * CiHD / LambdaNP2
                -177142. * CiHWB / LambdaNP2
                +3105.71 * CiDHB / LambdaNP2
                +7205.25 * CiDHW / LambdaNP2
                +10660.4 * CiW / LambdaNP2
                -2.841 * DeltaGF()
                -0.773 * deltaMwd6(); 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.542 * deltaMz()
                -0.672 * deltaaMZ()
                +2.797 * deltaGmu() );
    
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                -45448.7 * CiHL1_11 / LambdaNP2
                -208484. * CiHe_11 / LambdaNP2
                +274583. * CiHL3_11 / LambdaNP2
                -80024.1 * CiHD / LambdaNP2
                -97902.7 * CiHWB / LambdaNP2
                +28562.8 * CiDHB / LambdaNP2
                +575.898 * CiDHW / LambdaNP2
                +8122.74 * CiW / LambdaNP2
                -2.687 * DeltaGF()
                -0.928 * deltaMwd6(); 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.257 * deltaMz()
                -0.496 * deltaaMZ()
                +2.607 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                -47903.7 * CiHL1_11 / LambdaNP2
                -2144.19 * CiHe_11 / LambdaNP2
                +290349. * CiHL3_11 / LambdaNP2
                -84405.4 * CiHD / LambdaNP2
                -176530. * CiHWB / LambdaNP2
                +3309.62 * CiDHB / LambdaNP2
                +7174.21 * CiDHW / LambdaNP2
                +10675.5 * CiW / LambdaNP2
                -2.84 * DeltaGF()
                -0.777 * deltaMwd6(); 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.543 * deltaMz()
                -0.674 * deltaaMZ()
                +2.798 * deltaGmu() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeWWPol()");
        }
        
    } else if (sqrt_s == 0.365) {
        
        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                -45618.2 * CiHL1_11 / LambdaNP2
                -382668. * CiHe_11 / LambdaNP2
                +265703. * CiHL3_11 / LambdaNP2
                -77085.4 * CiHD / LambdaNP2
                -38791. * CiHWB / LambdaNP2
                +51079.9 * CiDHB / LambdaNP2
                -3972.2 * CiDHW / LambdaNP2
                +6727.91 * CiW / LambdaNP2
                -2.582 * DeltaGF()
                -1.04 * deltaMwd6(); 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.09 * deltaMz()
                -0.349 * deltaaMZ()
                +2.483 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                -50230.7 * CiHL1_11 / LambdaNP2
                -1000.53 * CiHe_11 / LambdaNP2
                +291951. * CiHL3_11 / LambdaNP2
                -84657.2 * CiHD / LambdaNP2
                -177196. * CiHWB / LambdaNP2
                +3348.72 * CiDHB / LambdaNP2
                +7579.53 * CiDHW / LambdaNP2
                +10879.2 * CiW / LambdaNP2
                -2.84 * DeltaGF()
                -0.753 * deltaMwd6(); 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.576 * deltaMz()
                -0.681 * deltaaMZ()
                +2.795 * deltaGmu() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeWWPol()");
        }
    
    } else if (sqrt_s == 0.380) {
        
        if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                -49806.5 * CiHL1_11 / LambdaNP2
                -221155. * CiHe_11 / LambdaNP2
                +280445. * CiHL3_11 / LambdaNP2
                -80550.4 * CiHD / LambdaNP2
                -101476. * CiHWB / LambdaNP2
                +31723.3 * CiDHB / LambdaNP2
                +1672.16 * CiDHW / LambdaNP2
                +8838.57 * CiW / LambdaNP2
                -2.707 * DeltaGF()
                -0.891 * deltaMwd6(); 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.331 * deltaMz()
                -0.503 * deltaaMZ()
                +2.64 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                -52386.5 * CiHL1_11 / LambdaNP2
                -2537.08 * CiHe_11 / LambdaNP2
                +294134. * CiHL3_11 / LambdaNP2
                -84922.5 * CiHD / LambdaNP2
                -176871. * CiHWB / LambdaNP2
                +3635.55 * CiDHB / LambdaNP2
                +7973.68 * CiDHW / LambdaNP2
                +10984.7 * CiW / LambdaNP2
                -2.838 * DeltaGF()
                -0.753 * deltaMwd6(); 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.589 * deltaMz()
                -0.68 * deltaaMZ()
                +2.81 * deltaGmu() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeWWPol()");
        }
    
    } else if (sqrt_s == 0.500) {
        
        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                -64264.6 * CiHL1_11 / LambdaNP2
                -495727. * CiHe_11 / LambdaNP2
                +289682. * CiHL3_11 / LambdaNP2
                -80108.8 * CiHD / LambdaNP2
                -61678. * CiHWB / LambdaNP2
                +75403.3 * CiDHB / LambdaNP2
                +458.146 * CiDHW / LambdaNP2
                +8723.87 * CiW / LambdaNP2
                -2.664 * DeltaGF()
                -0.849 * deltaMwd6(); 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.362 * deltaMz()
                -0.496 * deltaaMZ()
                +2.591 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                -68310.7 * CiHL1_11 / LambdaNP2
                -1341.22 * CiHe_11 / LambdaNP2
                +311528. * CiHL3_11 / LambdaNP2
                -84984.5 * CiHD / LambdaNP2
                -178260. * CiHWB / LambdaNP2
                +5206.37 * CiDHB / LambdaNP2
                +10705.4 * CiDHW / LambdaNP2
                +11071.1 * CiW / LambdaNP2
                -2.855 * DeltaGF()
                -0.671 * deltaMwd6(); 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.728 * deltaMz()
                -0.698 * deltaaMZ()
                +2.817 * deltaGmu() );
    
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                -66178. * CiHL1_11 / LambdaNP2
                -274919. * CiHe_11 / LambdaNP2
                +299745. * CiHL3_11 / LambdaNP2
                -82524.6 * CiHD / LambdaNP2
                -113979. * CiHWB / LambdaNP2
                +43898.4 * CiDHB / LambdaNP2
                +5024.43 * CiDHW / LambdaNP2
                +9759.79 * CiW / LambdaNP2
                -2.752 * DeltaGF()
                -0.778 * deltaMwd6(); 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.515 * deltaMz()
                -0.602 * deltaaMZ()
                +2.695 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                -68435.6 * CiHL1_11 / LambdaNP2
                -3089.11 * CiHe_11 / LambdaNP2
                +310020. * CiHL3_11 / LambdaNP2
                -85227.7 * CiHD / LambdaNP2
                -178139. * CiHWB / LambdaNP2
                +5322.77 * CiDHB / LambdaNP2
                +10598. * CiDHW / LambdaNP2
                +11009.9 * CiW / LambdaNP2
                -2.846 * DeltaGF()
                -0.681 * deltaMwd6(); 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.725 * deltaMz()
                -0.695 * deltaaMZ()
                +2.828 * deltaGmu() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeWWPol()");
        }
        
    } else if (sqrt_s == 1.0) {
        
        if (Pol_em == 80. && Pol_ep == -20.){
            mu += 
                -145951. * CiHL1_11 / LambdaNP2
                -885593. * CiHe_11 / LambdaNP2
                +383080. * CiHL3_11 / LambdaNP2
                -83628.6 * CiHD / LambdaNP2
                -114732. * CiHWB / LambdaNP2
                +159832. * CiDHB / LambdaNP2
                +17735.5 * CiDHW / LambdaNP2
                +8916.37 * CiW / LambdaNP2
                -2.787 * DeltaGF()
                -0.57 * deltaMwd6() ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.793 * deltaMz()
                -0.653 * deltaaMZ()
                +2.677 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 20.){
            mu += 
                -150086. * CiHL1_11 / LambdaNP2
                -4395.1 * CiHe_11 / LambdaNP2
                +394641. * CiHL3_11 / LambdaNP2
                -85925.1 * CiHD / LambdaNP2
                -181046. * CiHWB / LambdaNP2
                +13333.6 * CiDHB / LambdaNP2
                +23871.2 * CiDHW / LambdaNP2
                +9450.35 * CiW / LambdaNP2
                -2.871 * DeltaGF()
                -0.492 * deltaMwd6() ; 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +5.001 * deltaMz()
                -0.752 * deltaaMZ()
                +2.79 * deltaGmu() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeWWPol()");
        }
        
    } else if (sqrt_s == 1.5) {
        
        if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                -261040. * CiHL1_11 / LambdaNP2
                -1059495. * CiHe_11 / LambdaNP2
                +500666. * CiHL3_11 / LambdaNP2
                -84992.3 * CiHD / LambdaNP2
                -144925. * CiHWB / LambdaNP2
                +205215. * CiDHB / LambdaNP2
                +38777.5 * CiDHW / LambdaNP2
                +7857.84 * CiW / LambdaNP2
                -2.817 * DeltaGF()
                -0.471 * deltaMwd6(); 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +4.975 * deltaMz()
                -0.718 * deltaaMZ()
                +2.688 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                -265008. * CiHL1_11 / LambdaNP2
                -13002.4 * CiHe_11 / LambdaNP2
                +507924. * CiHL3_11 / LambdaNP2
                -86313.9 * CiHD / LambdaNP2
                -182113. * CiHWB / LambdaNP2
                +24953.6 * CiDHB / LambdaNP2
                +42429.8 * CiDHW / LambdaNP2
                +8014.86 * CiW / LambdaNP2
                -2.857 * DeltaGF()
                -0.429 * deltaMwd6(); 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +5.094 * deltaMz()
                -0.768 * deltaaMZ()
                +2.739 * deltaGmu() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeWWPol()");
        }
    
    } else if (sqrt_s == 3.0) {
        
        if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                -776767. * CiHL1_11 / LambdaNP2
                -3168410. * CiHe_11 / LambdaNP2
                +1016120. * CiHL3_11 / LambdaNP2
                -85414.3 * CiHD / LambdaNP2
                -155729. * CiHWB / LambdaNP2
                +628130. * CiDHB / LambdaNP2
                +123368. * CiDHW / LambdaNP2
                +6454.34 * CiW / LambdaNP2
                -2.831 * DeltaGF()
                -0.352 * deltaMwd6(); 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +5.165 * deltaMz()
                -0.755 * deltaaMZ()
                +2.77 * deltaGmu() );
    
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                -785359. * CiHL1_11 / LambdaNP2
                -39533. * CiHe_11 / LambdaNP2
                +1027322. * CiHL3_11 / LambdaNP2
                -86621.7 * CiHD / LambdaNP2
                -184516. * CiHWB / LambdaNP2
                +75975.5 * CiDHB / LambdaNP2
                +127086. * CiDHW / LambdaNP2
                +6519.78 * CiW / LambdaNP2
                -2.86 * DeltaGF()
                -0.328 * deltaMwd6(); 
            
    // Add modifications due to small variations of the SM parameters    
            mu += cHSM * ( +5.246 * deltaMz()
                -0.79 * deltaaMZ()
                +2.81 * deltaGmu() );
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeWWPol()");
        }
    
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::mueeWWPol()");
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

////////////////////////////////////////////////////////////////////////
    
    //----- High Energy diboson observables at hadron colliders


double NPSMEFTd6::ppZHprobe(const double sqrt_s) const
{    

    double gpZ=0.0;
    
    double ghZuL,ghZdL,ghZuR,ghZdR;
    
    // In the Warsaw basis the contact interactions are generated only by CHF ops but
    // in the modified basis ODHB, ODHW also contribute
    
    ghZuL = -(eeMz/sW_tree/cW_tree)*(CiHQ1_11 - CiHQ3_11 + g1_tree * (1.0/12.0) * CiDHB - (g2_tree/4.0) * CiDHW) * v2_over_LambdaNP2;
    ghZdL = -(eeMz/sW_tree/cW_tree)*(CiHQ1_11 + CiHQ3_11 + g1_tree * (1.0/12.0) * CiDHB + (g2_tree/4.0) * CiDHW) * v2_over_LambdaNP2;
    ghZuR = -(eeMz/sW_tree/cW_tree)*(CiHu_11 + g1_tree * (1.0/3.0) * CiDHB) * v2_over_LambdaNP2;
    ghZdR = -(eeMz/sW_tree/cW_tree)*(CiHd_11 - g1_tree * (1.0/6.0) * CiDHB) * v2_over_LambdaNP2;   
    
    if (sqrt_s == 14.0) {
         
        gpZ = ghZuL - 0.76 * ghZdL - 0.45 * ghZuR + 0.14 * ghZdR;
        
    } else if (sqrt_s == 27.0) {
        // Use the same as for 14 TeV for the moment 
        
        gpZ = ghZuL - 0.76 * ghZdL - 0.45 * ghZuR + 0.14 * ghZdR;
        
    } else if (sqrt_s == 100.0) {
         
        gpZ = ghZuL - 0.90 * ghZdL - 0.45 * ghZuR + 0.17 * ghZdR;
        
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::ppZHprobe()");
    

    return gpZ;
    
}
    
double NPSMEFTd6::mupTVppWZ(const double sqrt_s, const double pTV1, const double pTV2) const
{    
    double mu = 1.0;
    
    double cHWp = 0.0;
    
    // In the Warsaw basis the contact interactions are generated only by CiHQ3 but
    // in the modified basis ODHW also contribute
    // Master Equations below are for cHWp = Ci/Lambda^2 in units of TeV^{-2}, 
    // but LambdaNP is in GeV. Add conversion factor.
    
    cHWp = 4.0 * (sW2_tree/eeMz2) * (CiHQ3_11 + (g2_tree/4.0) * CiDHW) * 1000000.0 / LambdaNP2;

//  Bin dependences assuming cutoff of the EFT at 5 TeV
//  Normalize to the total number of events to remove the dependence on Lumi
//  (Numbers correspond to 3/ab)
    if (sqrt_s == 14.0) {
        
        if (pTV1 == 100.){
            mu += (558.0 * cHWp + 56.8 * cHWp * cHWp) / 3450.0;          
    
        } else if (pTV1 == 150.){
            mu += (410.0 * cHWp + 17.64 * cHWp * cHWp) / 2690.0;           
    
        } else if (pTV1 == 220.){
            mu += (266.0 * cHWp + 45.6 * cHWp * cHWp) / 925.0;
    
        } else if (pTV1 == 300.){
            mu += (304.0 * cHWp + 108.0 * cHWp * cHWp) / 563.0;

        } else if (pTV1 == 500.){
            mu += (114.40 * cHWp + 96.8 * cHWp * cHWp) / 85.1 ;
    
        } else if (pTV1 == 750.){
            mu += (46.20 * cHWp + 86.8 * cHWp * cHWp) / 14.9;
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mupTVppWZ()");
        } 
        
    } else if (sqrt_s == 27.0) {
        
        if (pTV1 == 150.){
            mu += (824.0 * cHWp + 71.6 * cHWp * cHWp) / 5370.0;          
    
        } else if (pTV1 == 220.){
            mu += (510.0 * cHWp + 75.2 * cHWp * cHWp) / 2210.0;           
    
        } else if (pTV1 == 300.){
            mu += (808.0 * cHWp + 268.4 * cHWp * cHWp) / 1610.0;
    
        } else if (pTV1 == 500.){
            mu += (374.0 * cHWp + 308.0 * cHWp * cHWp) / 331.0;

        } else if (pTV1 == 750.){
            mu += (216.0 * cHWp + 420.0 * cHWp * cHWp) / 85.9;
    
        } else if (pTV1 == 1200.){
            mu += (78.2 * cHWp + 325.2 * cHWp * cHWp) / 10.0;
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mupTVppWZ()");
        } 
        
    } else if (sqrt_s == 100.0) {

        if (pTV1 == 220.){
            mu += (2000.0 * cHWp + 368.4 * cHWp * cHWp) / 8030.0;          
    
        } else if (pTV1 == 300.){
            mu += (2780.0 * cHWp + 1000.0 * cHWp * cHWp) / 7270.0;           
    
        } else if (pTV1 == 500.){
            mu += (1544.0 * cHWp + 1428.0 * cHWp * cHWp) / 2000.0;
    
        } else if (pTV1 == 750.){
            mu += (1256.0 * cHWp + 2668.0 * cHWp * cHWp) / 717.0;

        } else if (pTV1 == 1200.){
            mu += (678.0 * cHWp + 3400.0 * cHWp * cHWp) / 142.0;
    
        } else if (pTV1 == 1800.){
            mu += (234.0 * cHWp + 2540.0 * cHWp * cHWp) / 27.5;
    
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mupTVppWZ()");
        } 
    
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::mupTVppWZ()");    
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;

}



    ////////////////////////////////////////////////////////////////////////
    
    //----- Simplified Template Cross Sections Bins
    // NOTE: Not our own calculations. From https://twiki.cern.ch/twiki/bin/view/LHCPhysics/STXStoEFT

double NPSMEFTd6::STXS_ggH_VBFtopo_j3v(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 + 56.6*aiG + 5.5*ai3G + 4.36*ai2G;
    
    return STXSb;
}

double NPSMEFTd6::STXS_ggH_VBFtopo_j3(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 + 55.9*aiG + 9.04*ai3G + 8.1*ai2G;
    
    return STXSb;
}
    

double NPSMEFTd6::STXS_ggH0j(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 + 55.2*aiG + 0.362*ai3G + 0.276*ai2G;
    
    return STXSb;
}
    

double NPSMEFTd6::STXS_ggH1j_pTH_0_60(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 + 56.0*aiG + 1.52*ai3G + 1.19*ai2G;
    
    return STXSb;
}
    

double NPSMEFTd6::STXS_ggH1j_pTH_60_120(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 + 55.5*aiG + 4.12*ai3G + 2.76*ai2G;
    
    return STXSb;
}
    

double NPSMEFTd6::STXS_ggH1j_pTH_120_200(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 + 56.5*aiG + 17.8*ai3G + 11.2*ai2G;
    
    return STXSb;
}
    

double NPSMEFTd6::STXS_ggH1j_pTH_200(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 + 55.0*aiG + 52.0*ai3G + 34.0*ai2G;
    
    return STXSb;
}
    

double NPSMEFTd6::STXS_ggH2j_pTH_0_200(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    return STXSb;
}


double NPSMEFTd6::STXS_ggH2j_pTH_0_60(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 + 55.6*aiG + 3.66*ai3G + 4.23*ai2G;
    
    return STXSb;
}

double NPSMEFTd6::STXS_ggH2j_pTH_60_120(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 + 56.1*aiG + 7.73*ai3G + 6.81*ai2G;
    
    return STXSb;
}

double NPSMEFTd6::STXS_ggH2j_pTH_120_200(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 + 55.8*aiG + 23.0*ai3G + 17.5*ai2G;
    
    return STXSb;
}

double NPSMEFTd6::STXS_ggH2j_pTH_200(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 + 56.0*aiG + 89.8*ai3G + 68.1*ai2G;
    
    return STXSb;
}
    

double NPSMEFTd6::STXS_qqHqq_VBFtopo_Rest(const double sqrt_s) const{
    
    return STXS_qqHqq_Rest(sqrt_s);
}


double NPSMEFTd6::STXS_qqHqq_VBFtopo_j3v(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 + 1.256*aiWW - 0.02319*aiB - 4.31*aiHW - 0.2907*aiHB;
    
    return STXSb;
}

double NPSMEFTd6::STXS_qqHqq_VBFtopo_j3(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 + 1.204*aiWW - 0.02692*aiB - 5.76*aiHW - 0.4058*aiHB;
    
    return STXSb;
}
    

double NPSMEFTd6::STXS_qqHqq_VHtopo(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 + 1.389*aiWW - 0.0284*aiB - 6.23*aiHW - 0.417*aiHB;
    
    return STXSb;
}


double NPSMEFTd6::STXS_qqHqq_Rest(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 + 1.546*aiWW - 0.02509*aiB - 3.631*aiHW - 0.2361*aiHB;
    
    return STXSb;
}
    

double NPSMEFTd6::STXS_qqHqq_pTj_200(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 + 7.82*aiWW - 0.1868*aiB - 30.65*aiHW - 2.371*aiHB;
    
    return STXSb;
}
    

double NPSMEFTd6::STXS_qqHlv_pTV_0_250(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    return STXSb;
}


double NPSMEFTd6::STXS_qqHlv_pTV_0_150(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 - 1.001*aiH + 33.63*aiWW + 11.49*aiHW + 23.62*aipHQ + 2.013*aipHL;
    
    return STXSb;
}


double NPSMEFTd6::STXS_qqHlv_pTV_150_250_0j(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 - 0.998*aiH + 76.3*aiWW + 50.7*aiHW + 66.5*aipHQ + 2.03*aipHL;
    
    return STXSb;
}


double NPSMEFTd6::STXS_qqHlv_pTV_150_250_1j(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 - 1.006*aiH + 70.9*aiWW + 45.5*aiHW + 60.8*aipHQ + 2.04*aipHL;
    
    return STXSb;
}
    

double NPSMEFTd6::STXS_qqHlv_pTV_250(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 - 1.001*aiH + 196.5*aiWW + 169.4*aiHW + 186.3*aipHQ + 2.03*aipHL;
    
    return STXSb;
}
    

double NPSMEFTd6::STXS_qqHll_pTV_0_150(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 - 1.0*aiH - 4.001*aiT + 29.82*aiWW + 8.43*aiB + 8.5*aiHW 
            + 2.545*aiHB + 0.0315*aiA - 1.89*aiHQ + 22.84*aipHQ + 5.247*aiHu 
            - 2.0*aiHd - 0.963*aiHL + 2.042*aipHL - 0.2307*aiHe;
    
    return STXSb;
}
    

double NPSMEFTd6::STXS_qqHll_pTV_150_250(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    return STXSb;
}


double NPSMEFTd6::STXS_qqHll_pTV_150_250_0j(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 - 0.993*aiH - 4.0*aiT + 62.4*aiWW + 18.08*aiB + 37.6*aiHW 
            + 11.22*aiHB - 5.03*aiHQ + 61.0*aipHQ + 14.39*aiHu - 5.17*aiHd 
            - 0.977*aiHL + 2.08*aipHL - 0.234*aiHe;
    
    return STXSb;
}


double NPSMEFTd6::STXS_qqHll_pTV_150_250_1j(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 - 1.002*aiH - 4.01*aiT + 57.9*aiWW + 16.78*aiB + 32.8*aiHW 
            + 9.86*aiHB - 4.58*aiHQ + 55.6*aipHQ + 13.54*aiHu - 4.56*aiHd 
            - 0.989*aiHL + 2.09*aipHL - 0.235*aiHe;
    
    return STXSb;
}
    

double NPSMEFTd6::STXS_qqHll_pTV_250(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 - 0.998*aiH - 4.0*aiT + 153.1*aiWW + 45.6*aiB + 126.4*aiHW 
            + 37.9*aiHB - 13.85*aiHQ + 168.6*aipHQ + 41.7*aiHu - 13.48*aiHd 
            - 0.977*aiHL + 2.09*aipHL - 0.238*aiHe;
    
    return STXSb;
}
    
    
double NPSMEFTd6::STXS_ttHtH(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 - 0.983*aiH + 2.949*aiu + 0.928*aiG + 313.6*aiuG 
            + 27.48*ai3G - 13.09*ai2G;
    
    return STXSb;
}

double NPSMEFTd6::STXS_WHqqHqq_VBFtopo_j3v(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 - 0.94*aiH + 39.5*aiWW + 13.8*aiHW + 32.1*aipHQ;
    
    return STXSb;
}
    
double NPSMEFTd6::STXS_WHqqHqq_VBFtopo_j3(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 - 1.04*aiH + 44.9*aiWW + 20.3*aiHW + 36.8*aipHQ;
    
    return STXSb;
}
    
double NPSMEFTd6::STXS_WHqqHqq_VH2j(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 - 0.996*aiH + 45.57*aiWW + 23.66*aiHW + 37.55*aipHQ;
    
    return STXSb;
}
    
double NPSMEFTd6::STXS_WHqqHqq_Rest(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 - 1.002*aiH + 34.29*aiWW + 11.56*aiHW + 26.27*aipHQ;
    
    return STXSb;
}
    
double NPSMEFTd6::STXS_WHqqHqq_pTj1_200(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 - 1.003*aiH + 181.2*aiWW + 152.3*aiHW + 173.7*aipHQ;
    
    return STXSb;
}
    
double NPSMEFTd6::STXS_ZHqqHqq_VBFtopo_j3v(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 - 0.94*aiH - 4.0*aiT + 34.8*aiWW + 10.0*aiB + 9.9*aiHW 
            + 3.04*aiHB - 2.14*aiHQ + 31.1*aipHQ + 7.6*aiHu - 2.59*aiHd;
    
    return STXSb;
}
    
double NPSMEFTd6::STXS_ZHqqHqq_VBFtopo_j3(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 - 0.97*aiH - 3.98*aiT + 38.1*aiWW + 10.5*aiB + 14.2*aiHW 
            + 4.15*aiHB - 2.36*aiHQ + 34.5*aipHQ + 8.4*aiHu - 2.79*aiHd;
    
    return STXSb;
}
    
double NPSMEFTd6::STXS_ZHqqHqq_VH2j(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 - 0.998*aiH - 4.002*aiT + 37.99*aiWW + 10.47*aiB + 16.45*aiHW 
            + 4.927*aiHB - 2.401*aiHQ + 34.45*aipHQ + 7.94*aiHu - 2.993*aiHd;
    
    return STXSb;
}
    
double NPSMEFTd6::STXS_ZHqqHqq_Rest(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 - 1.001*aiH - 3.998*aiT + 30.89*aiWW + 8.35*aiB + 8.71*aiHW 
            + 2.616*aiHB - 1.782*aiHQ + 26.1*aipHQ + 5.942*aiHu - 2.305*aiHd;
    
    return STXSb;
}
    
double NPSMEFTd6::STXS_ZHqqHqq_pTj1_200(const double sqrt_s) const{
    
    double STXSb = 1.0;
    
    STXSb = 1.0 - 1.003*aiH - 4.03*aiT + 141.5*aiWW + 41.6*aiB + 112.5*aiHW 
            + 33.6*aiHB - 11.52*aiHQ + 156.2*aipHQ + 38.9*aiHu - 12.53*aiHd;
    
    return STXSb;
}
    

///////////////////////////////////////////////////////////////////////////////

double NPSMEFTd6::kappamueff() const
{
      return sqrt(GammaHmumuRatio());
}

double NPSMEFTd6::kappataueff() const
{
      return sqrt(GammaHtautauRatio());
}

double NPSMEFTd6::kappaceff() const
{
      return sqrt(GammaHccRatio());
}

double NPSMEFTd6::kappabeff() const
{
      return sqrt(GammaHbbRatio());
}

double NPSMEFTd6::kappaGeff() const
{
      return sqrt(GammaHggRatio());
}

double NPSMEFTd6::kappaZeff() const
{
      return sqrt(GammaHZZRatio());
}

double NPSMEFTd6::kappaWeff() const
{
      return sqrt(GammaHWWRatio());
}

double NPSMEFTd6::kappaAeff() const
{
      return sqrt(GammaHgagaRatio());
}

double NPSMEFTd6::kappaZAeff() const
{
      return sqrt(GammaHZgaRatio());
}


/////////////Basic interactions of the so-called Higgs basis////////////////
    

double NPSMEFTd6::deltayt_HB() const
{
    double mf= mtpole;
    double ciHB;
    
    ciHB = - (v()/mf/sqrt(2.0))*CiuH_33r*v2_over_LambdaNP2 + delta_h - 0.5*DeltaGF();
    
    return ciHB;
}
    

double NPSMEFTd6::deltayb_HB() const
{
    double mf= (quarks[BOTTOM].getMass());
    double ciHB;
    
    ciHB = - (v()/mf/sqrt(2.0))*CidH_33r*v2_over_LambdaNP2 + delta_h - 0.5*DeltaGF();
    
    return ciHB;
}
    

double NPSMEFTd6::deltaytau_HB() const
{
    double mf= (leptons[TAU].getMass());
    double ciHB;
    
    ciHB = - (v()/mf/sqrt(2.0))*CieH_33r*v2_over_LambdaNP2 + delta_h - 0.5*DeltaGF();
    
    return ciHB;
}
    

double NPSMEFTd6::deltayc_HB() const
{
    double mf= (quarks[CHARM].getMass());
    double ciHB;
    
    ciHB = - (v()/mf/sqrt(2.0))*CiuH_22r*v2_over_LambdaNP2 + delta_h - 0.5*DeltaGF();
    
    return ciHB;
}
    

double NPSMEFTd6::deltaymu_HB() const
{
    double mf= (leptons[MU].getMass());
    double ciHB;
    
    ciHB = - (v()/mf/sqrt(2.0))*CieH_22r*v2_over_LambdaNP2 + delta_h - 0.5*DeltaGF();
    
    return ciHB;
}
    

double NPSMEFTd6::deltacZ_HB() const
{
    double ciHB;
    
    ciHB = delta_h - (3.0/2.0)*DeltaGF();
    
    return ciHB;
}
    

double NPSMEFTd6::cZBox_HB() const
{
    double ciHB;
    
    ciHB = (sW2_tree/eeMz2)*( DeltaGF() + 0.5*CiHD*v2_over_LambdaNP2 );
    
    ciHB = ciHB + 0.5*(sW2_tree/eeMz)*(CiDHB / cW_tree + CiDHW / sW_tree)*v2_over_LambdaNP2; // Extra, not in Warsaw basis
    
    return ciHB;
}
    

double NPSMEFTd6::cZZ_HB() const
{
    double ciHB;
    
    ciHB = (4.0*sW2_tree*cW2_tree/eeMz2)*( cW2_tree*CiHW + sW2_tree*CiHB + sW_tree*cW_tree*CiHWB )*v2_over_LambdaNP2;
    
    ciHB = ciHB - (sW2_tree*cW2_tree/eeMz)*(CiDHB / cW_tree + CiDHW / sW_tree)*v2_over_LambdaNP2; // Extra, not in Warsaw basis
    
    return ciHB;
}
    

double NPSMEFTd6::cZga_HB() const
{
    double ciHB;
    
    ciHB = (sW2_tree*cW2_tree/eeMz2)*( 4.0*CiHW - 4.0*CiHB - (2.0*(cW2_tree-sW2_tree)/sW_tree/cW_tree)*CiHWB )*v2_over_LambdaNP2;
    
    ciHB = ciHB + 0.5*(sW_tree*cW_tree/eeMz)*(CiDHB / sW_tree - CiDHW / cW_tree)*v2_over_LambdaNP2; // Extra, not in Warsaw basis
    
    return ciHB;
}
    

double NPSMEFTd6::cgaga_HB() const
{
    double ciHB;
    
    ciHB = (4.0/eeMz2)*( sW2_tree*CiHW + cW2_tree*CiHB - sW_tree*cW_tree*CiHWB )*v2_over_LambdaNP2;
    
    return ciHB;
}
    

double NPSMEFTd6::cgg_HB() const
{
    double ciHB;
    
    ciHB = (1.0/(M_PI * AlsMz))*CHG*v2_over_LambdaNP2;
    
    return ciHB;
}

double NPSMEFTd6::cggEff_HB() const
{
    double ciHB;
    
    double m_t = mtpole;
    //doulbe m_t = quarks[TOP].getMass();
    double m_b = quarks[BOTTOM].getMass();
    double m_c = quarks[CHARM].getMass();
    
    double At = deltayt_HB() * AH_f(4.0 * m_t * m_t / mHl / mHl).real();
    double Ab = deltayb_HB() * AH_f(4.0 * m_b * m_b / mHl / mHl).real();
    double Ac = deltayc_HB() * AH_f(4.0 * m_c * m_c / mHl / mHl).real();
    
    ciHB = cgg_HB() + (1.0/16.0/M_PI/M_PI) * (At + Ab + Ac) ;
    
    return ciHB;
}
    

double NPSMEFTd6::lambz_HB() const
{
    double ciHB;
    
    ciHB = -(3.0/2.0)*(eeMz/sW_tree)*CiW*v2_over_LambdaNP2;
    
    return ciHB;
}

////////////////////////////Auxiliary observables//////////////////////////////

double NPSMEFTd6::AuxObs_NP1() const
{
    // To be used for some temporary observable
    
    // WY analysis at 13 TeV for HL-LHC 3/ab
    double Wpar, Ypar, Wpar2, Ypar2;
    double Chi2NC13, Chi2CC13, Chi2Tot;
    
    Wpar = 10000.0 * obliqueW();
    Ypar = 10000.0 * obliqueY();
    
    Wpar2 = Wpar*Wpar;
    Ypar2 = Ypar*Ypar;
    
    Chi2CC13 = Wpar2 * (18.365037149441695 + 2.422904241798858 * Wpar + 0.12120594308623695 * Wpar2);
    
    Chi2NC13 = 0.032772034538390675 * Wpar2*Wpar2 + 2.815243944990361 * Ypar2 - 0.36522061776278516 * Ypar2*Ypar 
            + 0.017375258924241194 * Ypar2*Ypar2 + Wpar2*Wpar * (-0.7059117582389635 + 0.006816297425306027 * Ypar) 
            + Wpar * Ypar * (7.988302197022343 + Ypar * (-0.5450119819316416 + 0.0050292149953719766 * Ypar)) 
            + Wpar2 * (5.68581760491364 + Ypar * (-0.5794111075840261 + 0.048026245835369625 * Ypar));
    
    Chi2Tot = Chi2CC13 + Chi2NC13;
    
    // To be used as Gaussian observable with mean=0, var=1 I must return the sqrt.
    return sqrt(Chi2Tot);
}

double NPSMEFTd6::AuxObs_NP2() const
{
    // To be used for some temporary observable
    
    // WY analysis at 13 TeV for HL-LHC 3/ab for the CC
    // WY analysis at 27 TeV for HE-LHC 15/ab for the NC. 5% systematics (corr and uncorr)
    double Wpar, Ypar, Wpar2, Ypar2;
    double Chi2NC27, Chi2CC13, Chi2Tot;
    
    Wpar = 10000.0 * obliqueW();
    Ypar = 10000.0 * obliqueY();
    
    Wpar2 = Wpar*Wpar;
    Ypar2 = Ypar*Ypar;
    
    Chi2CC13 = Wpar2 * (18.365037149441695 + 2.422904241798858 * Wpar + 0.12120594308623695 * Wpar2);
    
    Chi2NC27 = 21.139285368181907 * Wpar2*Wpar2 + Wpar2*Wpar * (-89.16828370317616 + 7.182929295852857 * Ypar) 
            + Wpar * Ypar * (208.8092257396059 + Ypar * (-81.00102926445666 + 6.203591096144735 * Ypar)) 
            + Ypar2 * (81.01075991905888 + Ypar * (-58.822719932531164 + 14.670206406369107 * Ypar)) 
            + Wpar2 * (136.70787790194357 + Ypar * (-86.48485007990255 + 35.67671393730628 * Ypar));
    
    Chi2Tot = Chi2CC13 + Chi2NC27;
    
    // To be used as Gaussian observable with mean=0, var=1 I must return the sqrt.
    return sqrt(Chi2Tot);
}

double NPSMEFTd6::AuxObs_NP3() const
{
    // To be used for some temporary observable
    
    // WY analysis at 13 TeV for HL-LHC 3/ab for the CC
    // WY analysis at 27 TeV for HE-LHC 15/ab for the NC. 1% systematics (corr and uncorr)
    double Wpar, Ypar, Wpar2, Ypar2;
    double Chi2NC27, Chi2CC13, Chi2Tot;
    
    Wpar = 10000.0 * obliqueW();
    Ypar = 10000.0 * obliqueY();
    
    Wpar2 = Wpar*Wpar;
    Ypar2 = Ypar*Ypar;
    
    Chi2CC13 = Wpar2 * (18.365037149441695 + 2.422904241798858 * Wpar + 0.12120594308623695 * Wpar2);
    
    Chi2NC27 = 25.148424251427552 * Wpar2*Wpar2 + Wpar2*Wpar * (-105.31753344410277 + 8.01723084630248 * Ypar) 
            + Wpar * Ypar * (253.11721255992683 + Ypar * (-93.18990615818014 + 6.8250043104055816 * Ypar)) 
            + Ypar2 * (97.52107126224298 + Ypar * (-67.961770347904945 + 16.80046890875678 * Ypar)) 
            + Wpar2 * (166.84179829911304 + Ypar * (-100.88118582829852 + 41.55424691040131 * Ypar));
        
    Chi2Tot = Chi2CC13 + Chi2NC27;
    
    // To be used as Gaussian observable with mean=0, var=1 I must return the sqrt.
    return sqrt(Chi2Tot);
}

double NPSMEFTd6::AuxObs_NP4() const
{
    // WH distribution at 14 TeV: From 1704.01953 + hvqq terms
    
    double Bin1 = 1.0, Bin2 = 1.0, Bin3 = 1.0, Bin4 = 1.0, Bin5 = 1.0;
    
    double dVud = 0.0, dVcs = 0.0;
    double dcZ = 0.0, cZBox = 0.0, cZZ = 0.0, cZA = 0.0, cAA = 0.0;
    
    double C11 = 0.0178, C12 = 0.0144, C13 = 0.0102, C14 = 0.0052, C15 = 0.0006;
    
    double dchi2;
    
//  Production in each bin (signal strength)
    
    Bin1 += 12.8 * dVud + 1.75 * dVcs
            + 2.00 * dcZ + 5.01 * cZBox + 2.72 * cZZ - 0.0267 * cZA - 0.0217 * cAA;
    
//  Linear contribution from Higgs self-coupling
    Bin1 = Bin1 + cLHd6*(C11 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    Bin1 = Bin1 + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    Bin2 += 15.3 * dVud + 1.91 * dVcs
            + 2.00 * dcZ + 5.81 * cZBox + 3.10 * cZZ - 0.0337 * cZA - 0.0255 * cAA;
    
//  Linear contribution from Higgs self-coupling
    Bin2 = Bin2 + cLHd6*(C12 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    Bin2 = Bin2 + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    Bin3 += 20.7 * dVud + 2.49 * dVcs
            + 2.01 * dcZ + 7.44 * cZBox + 3.76 * cZZ - 0.0535 * cZA - 0.0340 * cAA;
    
//  Linear contribution from Higgs self-coupling
    Bin3 = Bin3 + cLHd6*(C13 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    Bin3 = Bin3 + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    Bin4 += 35.1 * dVud + 3.63 * dVcs
            + 1.98 * dcZ + 11.8 * cZBox + 5.40 * cZZ - 0.112 * cZA - 0.0572 * cAA;
    
//  Linear contribution from Higgs self-coupling
    Bin4 = Bin4 + cLHd6*(C14 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    Bin4 = Bin4 + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();

    Bin5 += 67.7 * dVud + 5.41 * dVcs
            + 2.03 * dcZ + 22.6 * cZBox + 9.05 * cZZ - 0.276 * cZA - 0.117 * cAA;
    
//  Linear contribution from Higgs self-coupling
    Bin5 = Bin5 + cLHd6*(C15 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    Bin5 = Bin5 + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
//  Compute Chi square using only the last bin and the diphoton, ZZ and bb channels
    dchi2 = ( Bin5 * BrHZZ4lRatio() - 1.0 ) * ( Bin5 * BrHZZ4lRatio() - 1.0 )/(0.07*0.07 + 0.48*0.48)
            + ( Bin5 * BrHgagaRatio() - 1.0 ) * ( Bin5 * BrHgagaRatio() - 1.0 )/(0.08*0.08 + 0.54*0.54)
            + ( Bin5 * BrHbbRatio() - 1.0 ) * ( Bin5 * BrHbbRatio() - 1.0 )/(0.33*0.33 + 0.61*0.61);

    // To be used as Gaussian observable with mean=0, var=1 I must return the sqrt.
    return sqrt(dchi2);
}

double NPSMEFTd6::AuxObs_NP5() const
{
    // ZH distribution at 14 TeV: From 1704.01953 + hvqq terms
    
    double Bin1 = 1.0, Bin2 = 1.0, Bin3 = 1.0, Bin4 = 1.0, Bin5 = 1.0;
    
    double dgLZuu = 0.0, dgRZuu = 0.0, dgLZcc = 0.0, dgRZcc = 0.0;
    double dgLZdd = 0.0, dgRZdd = 0.0, dgLZss = 0.0, dgRZss = 0.0;
    
    double dcZ = 0.0, cZBox = 0.0, cZZ = 0.0, cZA = 0.0, cAA = 0.0;
    
    double C11 = 0.0208, C12 = 0.0164, C13 = 0.0112, C14 = 0.0051, C15 = 0.0021;
    
    double dchi2;
    
//  Production in each bin (signal strength)
    
    Bin1 += 14.6 * dgLZuu - 6.74 * dgRZuu - 11.6 * dgLZdd + 2.28 * dgRZdd
            + 1.35 * dgLZcc - 0.589 * dgRZcc - 2.35 * dgLZss + 0.431 * dgRZss
            + 2.01 * dcZ + 4.14 * cZBox + 2.12 * cZZ - 0.0237 * cZA - 0.0126 * cAA;
    
//  Linear contribution from Higgs self-coupling
    Bin1 = Bin1 + cLHd6*(C11 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    Bin1 = Bin1 + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    Bin2 += 16.2 * dgLZuu - 7.77 * dgRZuu - 13.4 * dgLZdd + 2.63 * dgRZdd
            + 1.44 * dgLZcc - 0.668 * dgRZcc - 2.52 * dgLZss + 0.462 * dgRZss
            + 2.01 * dcZ + 4.86* cZBox + 2.49 * cZZ - 0.0284 * cZA - 0.0156 * cAA;
    
//  Linear contribution from Higgs self-coupling
    Bin2 = Bin2 + cLHd6*(C12 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    Bin2 = Bin2 + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    Bin3 += 23.0* dgLZuu - 10.8 * dgRZuu - 19.0 * dgLZdd + 3.64 * dgRZdd
            + 1.88 * dgLZcc - 0.891 * dgRZcc - 3.19 * dgLZss + 0.591 * dgRZss
            + 2.00 * dcZ + 6.35 * cZBox + 3.02 * cZZ - 0.0448 * cZA - 0.0221 * cAA;
    
//  Linear contribution from Higgs self-coupling
    Bin3 = Bin3 + cLHd6*(C13 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    Bin3 = Bin3 + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
    Bin4 += 39.2 * dgLZuu - 18.4 * dgRZuu - 31.4 * dgLZdd + 5.88 * dgRZdd
            + 2.78 * dgLZcc - 1.36 * dgRZcc - 4.64 * dgLZss + 0.919 * dgRZss
            + 1.98 * dcZ + 10.5 * cZBox + 4.44 * cZZ - 0.0873 * cZA - 0.0396 * cAA;
    
//  Linear contribution from Higgs self-coupling
    Bin4 = Bin4 + cLHd6*(C14 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    Bin4 = Bin4 + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();

    Bin5 += 73.4 * dgLZuu - 35.5 * dgRZuu - 58.5 * dgLZdd + 11.2 * dgRZdd
            + 4.13 * dgLZcc - 1.95 * dgRZcc - 6.97 * dgLZss + 1.41 * dgRZss
            + 1.96 * dcZ + 20.3 * cZBox + 7.27 * cZZ - 0.193 * cZA - 0.0800 * cAA;
    
//  Linear contribution from Higgs self-coupling
    Bin5 = Bin5 + cLHd6*(C15 + 2.0*dZH)*deltaG_hhhRatio();
//  Quadratic contribution from Higgs self-coupling: add separately from FlagQuadraticTerms
    Bin5 = Bin5 + cLHd6*cLH3d62*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
    
//  Compute Chi square using only the last bin and the diphoton, ZZ and bb channels
    dchi2 = ( Bin5 * BrHZZ4lRatio() - 1.0 ) * ( Bin5 * BrHZZ4lRatio() - 1.0 )/(0.09*0.09 + 0.65*0.65)
            + ( Bin5 * BrHgagaRatio() - 1.0 ) * ( Bin5 * BrHgagaRatio() - 1.0 )/(0.03*0.03 + 0.99*0.99)
            + ( Bin5 * BrHbbRatio() - 1.0 ) * ( Bin5 * BrHbbRatio() - 1.0 )/(0.10*0.10 + 0.34*0.34);

    // To be used as Gaussian observable with mean=0, var=1 I must return the sqrt.
    return sqrt(dchi2);
}

double NPSMEFTd6::AuxObs_NP6() const
{
    // To be used for some temporary observable
    
    // HL-LHC DiHiggs invariant mass distribution: 14 TeV 3/ab
    
    double Chi2Tot;

//  NP in decays
    double dGH2,dGgaga,dGbb, dBRTot;
    
//  Contributions from the different bins    
    double Bin1,Bin2,Bin3,Bin4,Bin5,Bin6;
    double LLBin1, LLBin2, LLBin3, LLBin4, LLBin5, LLBin6;

//  Higgs basis parameters
    double dcZHB,cZboxHB,cZZHB,cZgaHB,cgagaHB,cggHB;
    double dytHB,dybHB,dytauHB;
    double dKlambda;
    
    dcZHB = deltacZ_HB();
    cZboxHB = cZBox_HB();
    cZZHB = cZZ_HB();
    
// In the paper it seems they use diff. norm but in the chi 2.nb
//  they translate into that convention, so I assume their calculation
//  is directly in the HB for the following 3 couplings 
    cZgaHB = cZga_HB();
    cgagaHB = cgaga_HB();
    cggHB = cgg_HB();

    dytHB = deltayt_HB();
    dybHB = deltayb_HB();
    dytauHB = deltaytau_HB();
    
    dKlambda = deltaG_hhhRatio();
        
//  Corrections to the different Higgs widths
    dGH2 = 1. + 0.010512791990056657 * cZboxHB 
            -  0.003819752423722165 * cZZHB + 0.0016024991450954641 * cZgaHB 
            -  0.0005968238492400916 * (2.8975474398595105 * cZboxHB 
            +  1.8975474398595107 * cZZHB - cZgaHB - 0.3426378481886507 * cgagaHB) 
            +  0.0990750425382019 * (1.4487737199297552 * cZboxHB + 0.44877371992975534 * cZZHB 
            -  0.2365019764475461 * cZgaHB - 0.08103452830235015 * cgagaHB) 
            -  0.0330404571742506 * (cZZHB + 0.4730039528950922 * cZgaHB + 0.055933184863595636 * cgagaHB) 
            -  0.00033171593951211893 * cgagaHB + 0.48287726036165796 * dcZHB 
            +  1.1541846695471276 * dybHB + 0.12642022723635785 * dytauHB 
            +  0.1704272683629381 * (0. + 118.68284969347252 * cggHB 
            -  0.031082871395970327 * dybHB + 1.034601498835783 * dytHB) 
            +  0.004560729716754681 * (0. - 12.079950077697095 * cgagaHB 
            +  1.2739859351743013 * dcZHB + 0.0022136399615102554 * dybHB 
            -  0.28081416399029446 * dytHB + 0.0036305606562964158 * dytauHB) 
            +  0.003080492878860618 * (0. - 17.021015025105033 * cZgaHB 
            +  1.0557935963831278 * dcZHB + 0.0006235357344154619 * dybHB 
            -  0.05644023795399054 * dytHB + 0.000023105836447458856 * dytauHB);
    
    dGH2 = dGH2 * dGH2;
    
    dGgaga = 1.0 + 2.0 * (0. - 12.079950077697095 * cgagaHB 
            + 1.2739859351743013 * dcZHB + 0.0022136399615102554 * dybHB 
            - 0.28081416399029446 * dytHB + 0.0036305606562964158 * dytauHB);
    
    dGbb = 1.0 + 2.0 * dybHB;
    
    dBRTot = dGbb * dGgaga / dGH2;
    
    // Bin 1
    Bin1 = 0.17*(1.0 + 3.9863794294589585 * cggHB 
            +    21.333394807321064 * cggHB*cggHB + 3.9527789724382836 * dcZHB 
            +    0.5566823785534646 * cggHB*dcZHB + 9.077153576669469 * dcZHB*dcZHB 
            -    7.713285621354339 * dytHB + 6.573887966178747 * cggHB*dytHB 
            -    45.88983201032187 * dcZHB*dytHB + 62.42156375416841 * dytHB*dytHB 
            +    4.257555672380181 * cggHB*dytHB*dytHB +   4.620310477256665 * dcZHB*dytHB*dytHB 
            -    9.403185493195476 * dytHB*dytHB*dytHB +   1.1563473213070041 * dytHB*dytHB*dytHB*dytHB 
            -    0.14505129596051047 * dKlambda -    0.1418831193390564 * cggHB*dKlambda 
            +    1.3502693869386464 * cggHB*cggHB*dKlambda -   0.6675315048183816 * dcZHB*dKlambda 
            -    0.002999558395846163 * cggHB*dcZHB*dKlambda 
            +    1.5448485758806263 * dytHB * dKlambda 
            -    0.005002986050963205 * cggHB*dytHB*dKlambda 
            -    0.6675315048183816 * dcZHB*dytHB * dKlambda 
            +    1.5222565251876392 * dytHB*dytHB * dKlambda 
            +    0.1278814581005547 * cggHB*dytHB*dytHB * dKlambda 
            -    0.1676433466534976 * dytHB*dytHB*dytHB * dKlambda 
            +    0.011296025346493552 * dKlambda*dKlambda 
            +    0.0014116654816114353 * cggHB*dKlambda*dKlambda 
            +    0.022260157195710357 * cggHB*cggHB*dKlambda*dKlambda 
            +    0.022592050692987104 * dytHB * dKlambda*dKlambda 
            +    0.0014116654816114353 * cggHB*dytHB*dKlambda*dKlambda 
            +    0.011296025346493552 * dytHB*dytHB * dKlambda*dKlambda);
    
    Bin1 = 0.67944 + Bin1 * dBRTot;
    
    // Exclude points with negative values of BinX
    if ( Bin1 < 0 ) return std::numeric_limits<double>::quiet_NaN();
    
    // Delta chi2 = -2*LL for the bin
    // Add an abs in the denominator of the log, 
    // even if events with negative BinX are not supposed to reach here.
    LLBin1 = 2.0 * (Bin1 - 0.84944 + 0.84944 * log( 0.84944 / fabs(Bin1) ) );
    
    // Bin 2
    Bin2 = 0.33*(1.0 + 1.8019627645351037 * cggHB 
            +    7.953163597932105 * cggHB*cggHB + 3.735123481549394 * dcZHB 
            -    2.654186900737259 * cggHB*dcZHB + 6.403420811368324 * dcZHB*dcZHB 
            -    6.991501690350679 * dytHB + 11.425848100026737 * cggHB*dytHB 
            -    30.219763494155394 * dcZHB*dytHB + 39.692409895713936 * dytHB*dytHB 
            +    1.661324633279857 * cggHB*dytHB*dytHB +    4.46563789250516 * dcZHB*dytHB*dytHB 
            -    8.710706509282613 * dytHB*dytHB*dytHB +    1.2361692069676826 * dytHB*dytHB*dytHB*dytHB 
            -    0.21386875429750188 * dKlambda +    0.2363972133088796 * cggHB*dKlambda 
            +    0.8549707073528667 * cggHB*cggHB*dKlambda -    0.7305144109557659 * dcZHB*dKlambda 
            -    0.14136602060890807 * cggHB*dcZHB*dKlambda +    1.50533606463443 * dytHB * dKlambda 
            +    0.747017712869579 * cggHB*dytHB*dKlambda -    0.7305144109557659 * dcZHB*dytHB * dKlambda 
            +    1.4607351592940678 * dytHB*dytHB * dKlambda 
            +    0.08652243773397514 * cggHB*dytHB*dytHB * dKlambda 
            -    0.25846965963786395 * dytHB*dytHB*dytHB * dKlambda 
            +    0.022300452670181038 * dKlambda*dKlambda +    0.009236644319657653 * cggHB*dKlambda*dKlambda 
            +    0.023125582948149842 * cggHB*cggHB*dKlambda*dKlambda 
            +    0.044600905340362075 * dytHB * dKlambda*dKlambda 
            +    0.009236644319657653 * cggHB*dytHB*dKlambda*dKlambda 
            +    0.022300452670181038 * dytHB*dytHB * dKlambda*dKlambda) ;
    
    Bin2 = 1.4312 + Bin2 * dBRTot;
    
    // Exclude points with negative values of BinX
    if ( Bin2 < 0 ) return std::numeric_limits<double>::quiet_NaN();
    
    // Delta chi2 = -2*LL for the bin
    // Add an abs in the denominator of the log, 
    // even if events with negative BinX are not supposed to reach here.
    LLBin2 = 2.0 * (Bin2 - 1.7612 + 1.7612 * log( 1.7612 / fabs(Bin2) ) );
    
    // Bin 3
    Bin3 = 0.99*(1.0 + 0.6707152151845268 * cggHB 
            +    4.113022405261353 * cggHB*cggHB + 3.4241906309399726 * dcZHB 
            -    2.9926046286644703 * cggHB*dcZHB + 4.72026565086762 * dcZHB*dcZHB 
            -    5.98522416048399 * dytHB + 10.012680455917307 * cggHB*dytHB 
            -    20.69102310585157 * dcZHB*dytHB + 26.4871108999121 * dytHB*dytHB 
            +    0.36415135473936855 * cggHB*dytHB*dytHB 
            +    4.206380168414172 * dcZHB*dytHB*dytHB - 7.688318821918381 * dytHB*dytHB*dytHB 
            +    1.3217369754941033 * dytHB*dytHB*dytHB*dytHB - 0.2873477323359291 * dKlambda 
            +    0.35631144357921507 * cggHB*dKlambda 
            +    0.6197019283831009 * cggHB*cggHB*dKlambda 
            -    0.7821895374741993 * dcZHB*dKlambda 
            -    0.23172596419155064 * cggHB*dcZHB*dKlambda 
            +    1.415746929098462 * dytHB * dKlambda 
            +    1.0816714186441074 * cggHB*dytHB*dKlambda 
            -    0.7821895374741993 * dcZHB*dytHB * dKlambda 
            +    1.3469684427821131 * dytHB*dytHB * dKlambda 
            +    0.030182082490240562 * cggHB*dytHB*dytHB * dKlambda 
            -    0.35612621865227795 * dytHB*dytHB*dytHB * dKlambda 
            +    0.03438924315817444 * dKlambda*dKlambda 
            +    0.019565500643816278 * cggHB*dKlambda*dKlambda 
            +    0.02382411268034237 * cggHB*cggHB*dKlambda*dKlambda 
            +    0.06877848631634888 * dytHB * dKlambda*dKlambda 
            +    0.019565500643816278 * cggHB*dytHB*dKlambda*dKlambda 
            +    0.03438924315817444 * dytHB*dytHB * dKlambda*dKlambda);
    
    Bin3 = 1.9764 + Bin3 * dBRTot;
    
    // Exclude points with negative values of BinX
    if ( Bin3 < 0 ) return std::numeric_limits<double>::quiet_NaN();
    
    // Delta chi2 = -2*LL for the bin
    // Add an abs in the denominator of the log, 
    // even if events with negative BinX are not supposed to reach here.
    LLBin3 = 2.0 * (Bin3 - 2.9664 + 2.9664 * log( 2.9664 / fabs(Bin3) ) );
    
    // Bin 4
    Bin4 = 2.86*(1.0 - 0.27406342847042814 * cggHB 
            +    1.9597360046161074 * cggHB*cggHB + 3.0113078755334115 * dcZHB 
            -    2.776019265892887 * cggHB*dcZHB + 3.1917709639679823 * dcZHB*dcZHB 
            -    4.6362529563760955 * dytHB + 7.377234185667426 * cggHB*dytHB 
            -    12.294598143269557 * dcZHB*dytHB + 15.407456380301479 * dytHB*dytHB 
            -    0.6767601835408067 * cggHB*dytHB*dytHB 
            +    3.844719765004924 * dcZHB*dytHB*dytHB 
            -    6.227970053277897 * dytHB*dytHB*dytHB +    1.4542592857563688 * dytHB*dytHB*dytHB*dytHB 
            -    0.39767067022413716 * dKlambda +    0.3661464075997459 * cggHB*dKlambda 
            +    0.4464409042746693 * cggHB*cggHB*dKlambda 
            -    0.8334118894715125 * dcZHB*dKlambda 
            -    0.3263197431214281 * cggHB*dcZHB*dKlambda 
            +    1.1940464266776625 * dytHB * dKlambda 
            +    1.2643073873631234 * cggHB*dytHB*dKlambda 
            -    0.8334118894715125 * dcZHB*dytHB * dKlambda 
            +    1.0808691956131988 * dytHB*dytHB * dKlambda 
            -    0.0807982496009068 * cggHB*dytHB*dytHB * dKlambda 
            -    0.5108479012886007 * dytHB*dytHB*dytHB * dKlambda 
            +    0.05658861553223176 * dKlambda*dKlambda 
            +    0.04424790213027415 * cggHB*dKlambda*dKlambda 
            +    0.02585578262020257 * cggHB*cggHB*dKlambda*dKlambda 
            +    0.11317723106446352 * dytHB * dKlambda*dKlambda 
            +    0.04424790213027415 * cggHB*dytHB*dKlambda*dKlambda 
            +    0.05658861553223176 * dytHB*dytHB * dKlambda*dKlambda);
    
    Bin4 = 5.167 + Bin4 * dBRTot;
    
    // Exclude points with negative values of BinX
    if ( Bin4 < 0 ) return std::numeric_limits<double>::quiet_NaN();
    
    // Delta chi2 = -2*LL for the bin
    // Add an abs in the denominator of the log, 
    // even if events with negative BinX are not supposed to reach here.
    LLBin4 = 2.0 * (Bin4 - 8.027 + 8.027 * log( 8.027 / fabs(Bin4) ) );
    
    // Bin 5
    Bin5 = 6.34* (1.0 - 1.094329254675176 * cggHB 
            +    1.0393648302909912 * cggHB*cggHB + 2.6000916816530903 * dcZHB 
            -    2.4448264513323226 * cggHB*dcZHB + 2.073935963891534 * dcZHB*dcZHB 
            -    3.192332240205929 * dytHB + 4.5914586198385 * cggHB*dytHB 
            -    6.2871857258718595 * dcZHB*dytHB + 8.134770266934664 * dytHB*dytHB 
            -    1.648691479483292 * cggHB*dytHB*dytHB +    3.5563383758242524 * dcZHB*dytHB*dytHB 
            -    4.615570013047001 * dytHB*dytHB*dytHB +    1.7227511548362076 * dytHB*dytHB*dytHB*dytHB 
            -    0.6079428047533413 * dKlambda +    0.33825211279194234 * cggHB*dKlambda 
            +    0.3879052211526028 * cggHB*cggHB*dKlambda -    0.956246694171162 * dcZHB*dKlambda 
            -    0.4572431444456198 * cggHB*dcZHB*dKlambda +    0.8152949680877302 * dytHB * dKlambda 
            +    1.3814632626914451 * cggHB*dytHB*dKlambda 
            -    0.956246694171162 * dcZHB*dytHB * dKlambda +    0.5856782679219981 * dytHB*dytHB * dKlambda 
            -    0.3285182834373566 * cggHB*dytHB*dytHB * dKlambda 
            -    0.8375595049190734 * dytHB*dytHB*dytHB * dKlambda +    0.11480835008286604 * dKlambda*dKlambda 
            +    0.11240817142118299 * cggHB*dKlambda*dKlambda +    0.03688252014841459 * cggHB*cggHB*dKlambda*dKlambda 
            +    0.22961670016573207 * dytHB * dKlambda*dKlambda 
            +    0.11240817142118299 * cggHB*dytHB*dKlambda*dKlambda 
            +    0.11480835008286604 * dytHB*dytHB * dKlambda*dKlambda);
    
    Bin5 = 15.93 + Bin5 * dBRTot;
    
    // Exclude points with negative values of BinX
    if ( Bin5 < 0 ) return std::numeric_limits<double>::quiet_NaN();
    
    // Delta chi2 = -2*LL for the bin
    // Add an abs in the denominator of the log, 
    // even if events with negative BinX are not supposed to reach here.
    LLBin5 = 2.0 * (Bin5 - 22.27 + 22.27 * log( 22.27 / fabs(Bin5) ) );
    
    // Bin 6
    Bin6 = 2.14*(1.0 - 2.007855065799201 * cggHB + 1.1994575008850934 * cggHB*cggHB 
            +    2.5987763498382352 * dcZHB - 2.908713303420072 * cggHB*dcZHB 
            +    1.804645897901265 * dcZHB*dcZHB - 2.806900956988577 * dytHB 
            +    3.5621616844486415 * cggHB*dytHB - 4.250685020965587 * dcZHB*dytHB 
            +    5.7468374752045515 * dytHB*dytHB - 3.1561231600123736 * cggHB*dytHB*dytHB 
            +    3.9784140166037667 * dcZHB*dytHB*dytHB - 4.4303353405513395 * dytHB*dytHB*dytHB 
            +    2.257739308366916 * dytHB*dytHB*dytHB*dytHB - 0.9894280925261291 * dKlambda 
            +    0.589956279744333 * cggHB*dKlambda +   0.6687315933211253 * cggHB*cggHB*dKlambda 
            -    1.3796376667655315 * dcZHB*dKlambda -    0.8069993678124955 * cggHB*dcZHB*dKlambda 
            +    0.6340062910366335 * dytHB * dKlambda +    2.127573647123277 * cggHB*dytHB*dKlambda 
            -    1.3796376667655315 * dcZHB*dytHB * dKlambda +   0.09738385935505989 * dytHB*dytHB * dKlambda 
            -    0.8833807360585424 * cggHB*dytHB*dytHB * dKlambda -    1.5260505242077027 * dytHB*dytHB*dytHB * dKlambda 
            +    0.2683112158407868 * dKlambda*dKlambda +   0.32506892158970235 * cggHB*dKlambda*dKlambda 
            +    0.09418943796384227 * cggHB*cggHB*dKlambda*dKlambda +    0.5366224316815736 * dytHB * dKlambda*dKlambda 
            +    0.32506892158970235 * cggHB*dytHB*dKlambda*dKlambda 
            +    0.2683112158407868 * dytHB*dytHB * dKlambda*dKlambda);
    
    Bin6 = 12.01 + Bin6 * dBRTot;
    
    // Exclude points with negative values of BinX
    if ( Bin6 < 0 ) return std::numeric_limits<double>::quiet_NaN();
    
    // Delta chi2 = -2*LL for the bin
    // Add an abs in the denominator of the log, 
    // even if events with negative BinX are not supposed to reach here.
    LLBin6 = 2.0 * (Bin6 - 14.15 + 14.15 * log( 14.15 / fabs(Bin6) ) );
    
    // The total contributions to the log-likelihood/chi-square
    Chi2Tot = LLBin1 + LLBin2 + LLBin3 + LLBin4 + LLBin5 + LLBin6;
        
    // To be used as Gaussian observable with mean=0, var=1 I must return the sqrt.
    return sqrt(Chi2Tot);
}

double NPSMEFTd6::AuxObs_NP7() const
{
    // To be used for some temporary observable
    
    // CLIC STWY using difermion production at all energies: 380, 1500 and 3000 GeV
    double Spar, Tpar, Wpar, Ypar, Spar2, Tpar2, Wpar2, Ypar2;
    double Chi2Tot;

    Spar = obliqueS();
    Tpar = obliqueT();
    Wpar = 10000.0 * obliqueW();
    Ypar = 10000.0 * obliqueY();

    Spar2 = Spar*Spar;
    Tpar2 = Tpar*Tpar;    
    Wpar2 = Wpar*Wpar;
    Ypar2 = Ypar*Ypar;
    
    Chi2Tot = 442.84977653097394 * Spar2 
            - 728.5215604181935 * Spar * Tpar 
            + 404.15957807101813 * Tpar2 
            + 400.03987723904224 * Spar * Wpar 
            - 639.6154242400826 * Tpar * Wpar 
            + 4337.791457515823 * Wpar2  
            - 106.87313892453362 * Spar * Ypar 
            - 72.94355609762007 * Tpar * Ypar 
            + 3002.848116515672 * Wpar * Ypar 
            + 3040.1630882458923 * Ypar2;
    
    // To be used as Gaussian observable with mean=0, var=1 I must return the sqrt.
    return sqrt(Chi2Tot);
}

double NPSMEFTd6::AuxObs_NP8() const
{
    // To be used for some temporary observable
    
    // CLIC DiHiggs: exclusive analysis. Full CLIC run   
    double Chi2Tot;

//  Higgs basis parameters
    double dKlambda;
        
    dKlambda = deltaG_hhhRatio();
    
    Chi2Tot = dKlambda * dKlambda * (50.04473972806045  
            - 104.47283225861888 * dKlambda 
            + 84.48333683635175 * dKlambda*dKlambda );
    
    // To be used as Gaussian observable with mean=0, var=1 I must return the sqrt.
    return sqrt(Chi2Tot);
}

double NPSMEFTd6::AuxObs_NP9() const
{
    // To be used for some temporary observable
    
    // ILC DiHiggs at 500 GeV: 2/ab per polarization (+-80,-+30) 
    
    double Chi2p80m30, Chi2m80p30, Chi2Tot;

//  Higgs basis parameters
    double dcZHB,cZboxHB,cZZHB,cZgaHB,cgagaHB;
    double dKlambda;
    
    dcZHB = deltacZ_HB();
    cZboxHB = cZBox_HB();
    cZZHB = cZZ_HB(); 
    cZgaHB = cZga_HB();
    cgagaHB = cgaga_HB();
    
    dKlambda = deltaG_hhhRatio();

//  The signal strength -1    
    Chi2p80m30 = 13.6982 * cZZHB 
            - 7.58943 * cZgaHB 
            + 14.6843 * cZboxHB 
            - 1.51882 * cgagaHB 
            + 5.46836 * dcZHB 
            + 0.565585 * dKlambda 
            + 0.000631004 * cZZHB * dKlambda 
            - 0.195079 * cZgaHB * dKlambda 
            + 0.064441 * cZboxHB * dKlambda 
            + 0.440061 * cgagaHB * dKlambda 
            + 2.13192 * dcZHB * dKlambda 
            + 0.0968208 * dKlambda * dKlambda;
    
//  ILC report (1903.01629) gives total cross section a 4/ab: 16.8%. 
// Assume the precision for each polarization is the same as they do for single Higgs in ZH...    
    Chi2p80m30 = Chi2p80m30 * Chi2p80m30 / 0.168 / 0.168 / 2.0;

//  The signal strength -1 
    Chi2m80p30 = - 2.57112 * cZZHB 
            + 6.97966 * cZgaHB 
            - 10.2626 * cZboxHB 
            + 1.39647 * cgagaHB 
            + 5.4684 * dcZHB 
            + 0.565577 * dKlambda 
            + 4.71916 * cZZHB * dKlambda 
            + 0.179045 * cZgaHB * dKlambda 
            + 7.28766 * cZboxHB * dKlambda 
            - 0.405166 * cgagaHB * dKlambda 
            + 2.13189 * dcZHB * dKlambda 
            + 0.0968201 * dKlambda * dKlambda;

//  ILC report (1903.01629) gives total cross section a 4/ab: 16.8%. 
// Assume the precision for each polarization is the same as they do for single Higgs in ZH...        
    Chi2m80p30 = Chi2m80p30 * Chi2m80p30 / 0.168 / 0.168 / 2.0;
    
    Chi2Tot = Chi2p80m30 + Chi2m80p30;

    // To be used as Gaussian observable with mean=0, var=1 I must return the sqrt.
    return sqrt(Chi2Tot);
}

double NPSMEFTd6::AuxObs_NP10() const
{
    // CLIC STWY using difermion production at all energies: 380 and 1500 GeV
    double Spar, Tpar, Wpar, Ypar, Spar2, Tpar2, Wpar2, Ypar2;
    double Chi2Tot;

    Spar = obliqueS();
    Tpar = obliqueT();
    Wpar = 10000.0 * obliqueW();
    Ypar = 10000.0 * obliqueY();

    Spar2 = Spar*Spar;
    Tpar2 = Tpar*Tpar;    
    Wpar2 = Wpar*Wpar;
    Ypar2 = Ypar*Ypar;
    
    Chi2Tot = 375.63808963031073 * Spar2 
            - 617.8864704052573 * Spar * Tpar 
            + 353.1650032169891 * Tpar2 
            + 215.96605851087603 * Spar * Wpar 
            - 309.3469843690006 * Tpar * Wpar 
            + 518.10263970583244 * Wpar2 
            - 45.972763923203014 * Spar * Ypar 
            - 40.670385844305705 * Tpar * Ypar 
            + 340.56677318671185 * Wpar * Ypar 
            + 364.5290176991845 * Ypar2;
    
    // To be used as Gaussian observable with mean=0, var=1 I must return the sqrt.
    return sqrt(Chi2Tot);
}

double NPSMEFTd6::AuxObs_NP11() const
{
    // CLIC STWY using difermion production at all energies: 380 GeV
    double Spar, Tpar, Wpar, Ypar, Spar2, Tpar2, Wpar2, Ypar2;
    double Chi2Tot;

    Spar = obliqueS();
    Tpar = obliqueT();
    Wpar = 10000.0 * obliqueW();
    Ypar = 10000.0 * obliqueY();

    Spar2 = Spar*Spar;
    Tpar2 = Tpar*Tpar;    
    Wpar2 = Wpar*Wpar;
    Ypar2 = Ypar*Ypar;
    
    Chi2Tot = 282.9842573293628 * Spar2 
            - 462.32090035841725 * Spar * Tpar 
            + 276.2496928300019 * Tpar2 
            + 66.08702076419566 * Spar * Wpar 
            - 87.95794393624075 * Tpar * Wpar 
            + 9.5435699879102 * Wpar2 
            - 26.170009941328716 * Spar * Ypar 
            - 9.695238064023518 * Tpar * Ypar 
            + 6.519573295893438 * Wpar * Ypar 
            + 12.858593910798793 * Ypar2;
    
    // To be used as Gaussian observable with mean=0, var=1 I must return the sqrt.
    return sqrt(Chi2Tot);
}

double NPSMEFTd6::AuxObs_NP12() const
{
    // CLIC dim6 Top fit 1500 GeV: only for SVF operators
    double CHqminus, CHt;
    double Chi2Tot;
    
    // The chi2 is given assuming C/Lambda^2 is in units of TeV^-2
    CHqminus= 0.5 * (CiHQ1_33 - CiHQ3_33) * (1000000.0 / LambdaNP2);
    CHt= 0.5 * CiHu_33 * (1000000.0 / LambdaNP2);
    
    Chi2Tot= 1203.58 * CHqminus * CHqminus + 1661.59 * CHqminus * CHt + 1257.83 * CHt * CHt;
    
    // To be used as Gaussian observable with mean=0, var=1 I must return the sqrt.
    return sqrt(Chi2Tot);
}

double NPSMEFTd6::AuxObs_NP13() const
{
    // CLIC dim6 Top fit 3000 GeV: only for SVF operators
    double CHqminus, CHt;
    double Chi2Tot;
    
    // The chi2 is given assuming C/Lambda^2 is in units of TeV^-2
    CHqminus= 0.5 * (CiHQ1_33 - CiHQ3_33) * (1000000.0 / LambdaNP2);
    CHt= 0.5 * CiHu_33 * (1000000.0 / LambdaNP2);
    
    Chi2Tot= 5756.01 * CHqminus * CHqminus + 8013.79 * CHqminus * CHt + 3380.7 * CHt * CHt;
    
    // To be used as Gaussian observable with mean=0, var=1 I must return the sqrt.
    return sqrt(Chi2Tot);
}

double NPSMEFTd6::AuxObs_NP14() const
{
    // Test chi2 for HH production at 100 TeV: only the first two bins in 1704.01953 are included,
    // with the same coefficients (including ratios of cross sections in each bin) its table 4.  The EFT parameterization of Higgs decays are not included.
    double Chi2Tot;
    
//  Higgs basis parameters
    double dcZHB,cggHB;
    double dytHB;
    double dKlambda;
    
    dcZHB = deltacZ_HB();
    cggHB = cgg_HB();
    dytHB = deltayt_HB();
    dKlambda = deltaG_hhhRatio();
    
    double dcZHB2,dcZHB3,dcZHB4;
    double cggHB2,cggHB3,cggHB4;
    double dytHB2,dytHB3,dytHB4,dytHB5,dytHB6,dytHB7,dytHB8;
    double dKlambda2,dKlambda3,dKlambda4;
    
    dcZHB2 = dcZHB * dcZHB;
    dcZHB3 = dcZHB2 * dcZHB;
    dcZHB4 = dcZHB3 * dcZHB;
    
    cggHB2 = cggHB * cggHB;
    cggHB3 = cggHB2 * cggHB;
    cggHB4 = cggHB3 * cggHB;

    dytHB2 = dytHB * dytHB;
    dytHB3 = dytHB2 * dytHB;
    dytHB4 = dytHB3 * dytHB;
    dytHB5 = dytHB4 * dytHB;
    dytHB6 = dytHB5 * dytHB;
    dytHB7 = dytHB6 * dytHB;
    dytHB8 = dytHB7 * dytHB;
    
    dKlambda2 = dKlambda * dKlambda;
    dKlambda3 = dKlambda2 * dKlambda;
    dKlambda4 = dKlambda3 * dKlambda;

    // The Chi2
    
    Chi2Tot = 2.0595082782796297e7 * cggHB2 - 3.6971136499764752e9 * cggHB3 + 1.7583900534677216e11 * cggHB4 
            - 630035.4483047676 * cggHB * dcZHB + 1.3588174266991532e8 * cggHB2 * dcZHB - 7.10364464231958e9 * cggHB3 * dcZHB 
            + 5311.651853836387 * dcZHB2 - 1.7067170379207395e6 * cggHB * dcZHB2 +  1.1851653627034137e8 * cggHB2 * dcZHB2 
            + 8180.119549200313 * dcZHB3 - 943018.2335425722 * cggHB * dcZHB3 + 3159.9135213745994 * dcZHB4 
            + 180518.97210352542 * cggHB * dKlambda - 2.8949546963646576e7 * cggHB2 * dKlambda - 5.501576225306801e8 * cggHB3 * dKlambda 
            + 1.5079027448500854e11 * cggHB4 * dKlambda - 2846.9365320948145 * dcZHB * dKlambda + 797208.485191074 * cggHB * dcZHB * dKlambda 
            - 4.978486710457227e6 * cggHB2 * dcZHB * dKlambda - 4.586348042437428e9 * cggHB3 * dcZHB * dKlambda - 6485.875373880575 * dcZHB2 * dKlambda 
            + 390177.86145601963 * cggHB * dcZHB2 * dKlambda + 5.056678567468029e7 * cggHB2 * dcZHB2 * dKlambda - 3291.6842405815532 * dcZHB3 * dKlambda 
            - 198301.99217208195 * cggHB * dcZHB3 * dKlambda + 399.29685823653153 * dKlambda2 - 95580.41780509672 * cggHB * dKlambda2 
            - 7.430874086734321e6 * cggHB2 * dKlambda2 + 7.720064658809748e8 * cggHB3 * dKlambda2 + 5.089872992160051e10 * cggHB4 * dKlambda2 
            + 1809.9095844013955 * dcZHB * dKlambda2 - 1150.4119995786175 * cggHB * dcZHB * dKlambda2 - 2.2786176268418655e7 * cggHB2 * dcZHB * dKlambda2 
            - 1.0351049455121036e9 * cggHB3 * dcZHB * dKlambda2 + 1362.5781363223641 * dcZHB2 * dKlambda2 + 170792.06609378837 * cggHB * dcZHB2 * dKlambda2 
            + 5.658917948194164e6 * cggHB2 * dcZHB2 * dKlambda2 - 178.77181321253659 * dKlambda3 - 11443.938844928987 * cggHB * dKlambda3 
            + 2.461878722072089e6 * cggHB2 * dKlambda3 + 2.821167791764089e8 * cggHB3 * dKlambda3 + 7.998289700049803e9 * cggHB4 * dKlambda3 
            - 267.7615464146533 * dcZHB * dKlambda3 - 52488.33374581051 * cggHB * dcZHB * dKlambda3 - 3.555711022595523e6 * cggHB2 * dcZHB * dKlambda3 
            - 8.149153208622633e7 * cggHB3 * dcZHB * dKlambda3 + 21.07398490236267 * dKlambda4 + 5735.3996792942135 * cggHB * dKlambda4 
            + 596986.3215027236 * cggHB2 * dKlambda4 + 2.773647081412465e7 * cggHB3 * dKlambda4 + 4.915460918180312e8 * cggHB4 * dKlambda4 
            + 740876.8879497008 * cggHB * dytHB - 1.938279550686329e8 * cggHB2 * dytHB + 1.1944585224312653e10 * cggHB3 * dytHB 
            - 12947.635844899749 * dcZHB * dytHB + 4.908519506685015e6 * cggHB * dcZHB * dytHB - 3.742271337006843e8 * cggHB2 * dcZHB * dytHB 
            - 33546.241370498166 * dcZHB2 * dytHB + 4.3134482870087875e6 * cggHB * dcZHB2 * dytHB - 18267.038917513022 * dcZHB3 * dytHB 
            + 3387.385955080094 * dKlambda * dytHB - 963072.1570381082 * cggHB * dKlambda * dytHB - 2.3453010760683898e7 * cggHB2 * dKlambda * dytHB 
            + 9.317798790237669e9 * cggHB3 * dKlambda * dytHB + 14461.190498065112 * dcZHB * dKlambda * dytHB - 276210.0620250288 * cggHB * dcZHB * dKlambda * dytHB 
            - 2.1850896154428744e8 * cggHB2 * dcZHB * dKlambda * dytHB + 7442.375770947524 * dcZHB2 * dKlambda * dytHB 
            + 1.6339998473341048e6 * cggHB * dcZHB2 * dKlambda * dytHB - 3291.6842405815532 * dcZHB3 * dKlambda * dytHB - 1559.6600507789517 * dKlambda2 * dytHB 
            - 212800.20942464058 * cggHB * dKlambda2 * dytHB + 3.499621075016396e7 * cggHB2 * dKlambda2 * dytHB + 2.9495867407085886e9 * cggHB3 * dKlambda2 * dytHB 
            - 132.54584108464164 * dcZHB * dKlambda2 * dytHB - 704650.5551856682 * cggHB * dcZHB * dKlambda2 * dytHB 
            - 4.6230021860231325e7 * cggHB2 * dcZHB * dKlambda2 * dytHB + 2725.1562726447282 * dcZHB2 * dKlambda2 * dytHB 
            + 170792.06609378837 * cggHB * dcZHB2 * dKlambda2 * dytHB - 174.87036642817392 * dKlambda3 * dytHB + 72002.66692264378 * cggHB * dKlambda3 * dytHB 
            + 1.2160354917437742e7 * cggHB2 * dKlambda3 * dytHB + 4.500393455278235e8 * cggHB3 * dKlambda3 * dytHB - 803.2846392439599 * dcZHB * dKlambda3 * dytHB 
            - 104976.66749162102 * cggHB * dcZHB * dKlambda3 * dytHB - 3.555711022595523e6 * cggHB2 * dcZHB * dKlambda3 * dytHB 
            + 84.29593960945068 * dKlambda4 * dytHB + 17206.19903788264 * cggHB * dKlambda4 * dytHB + 1.1939726430054472e6 * cggHB2 * dKlambda4 * dytHB 
            + 2.773647081412465e7 * cggHB3 * dKlambda4 * dytHB + 7985.615632692477 * dytHB2 - 4.312707242837639e6 * cggHB * dytHB2 
            + 4.446488644358661e8 * cggHB2 * dytHB2 - 5.669235052669609e9 * cggHB3 * dytHB2 + 59322.05816648064 * dcZHB * dytHB2 
            - 1.0048203483978426e7 * cggHB * dcZHB * dytHB2 + 2.009903412514487e8 * cggHB2 * dcZHB * dytHB2 + 64971.66315898899 * dcZHB2 * dytHB2 
            - 2.4669987769536236e6 * cggHB * dcZHB2 * dytHB2 + 11471.803789781865 * dcZHB3 * dytHB2 - 11811.249755773804 * dKlambda * dytHB2 
            + 431747.7364057698 * cggHB * dKlambda * dytHB2 + 2.2358583287946397e8 * cggHB2 * dKlambda * dytHB2 - 3.8910877145439386e9 * cggHB3 * dKlambda * dytHB2 
            - 16029.606555240167 * dcZHB * dKlambda * dytHB2 - 2.9253661324121524e6 * cggHB * dcZHB * dKlambda * dytHB2 
            + 8.987023921425158e7 * cggHB2 * dcZHB * dKlambda * dytHB2 + 4717.219498302798 * dcZHB2 * dKlambda * dytHB2 
            - 540895.9436706528 * cggHB * dcZHB2 * dKlambda * dytHB2 + 214.81067429237223 * dKlambda2 * dytHB2 + 567954.341114266 * cggHB * dKlambda2 * dytHB2 
            + 4.5123619667514816e7 * cggHB2 * dKlambda2 * dytHB2 - 9.277345617086976e8 * cggHB3 * dKlambda2 * dytHB2 
            - 3081.626211728115 * dcZHB * dKlambda2 * dytHB2 - 381097.4778098703 * cggHB * dcZHB * dKlambda2 * dytHB2 
            + 1.050966209735231e7 * cggHB2 * dcZHB * dKlambda2 * dytHB2 + 1362.5781363223641 * dcZHB2 * dKlambda2 * dytHB2 
            + 284.9520271687106 * dKlambda3 * dytHB2 + 127206.63260007375 * cggHB * dKlambda3 * dytHB2 + 6.267940600872645e6 * cggHB2 * dKlambda3 * dytHB2 
            - 7.655202990726441e7 * cggHB3 * dKlambda3 * dytHB2 - 803.2846392439599 * dcZHB * dKlambda3 * dytHB2 - 52488.33374581051 * cggHB * dcZHB * dKlambda3 * dytHB2 
            + 126.44390941417602 * dKlambda4 * dytHB2 + 17206.19903788264 * cggHB * dKlambda4 * dytHB2 + 596986.3215027236 * cggHB2 * dKlambda4 * dytHB2 
            - 37223.626257417236 * dytHB3 + 8.269994128894571e6 * cggHB * dytHB3 - 2.9221928856272686e8 * cggHB2 * dytHB3 - 105038.22976459829 * dcZHB * dytHB3 
            + 7.149383019204844e6 * cggHB * dcZHB * dytHB3 - 47474.492515326274 * dcZHB2 * dytHB3 + 11656.27418420629 * dKlambda * dytHB3 
            + 2.385352845620739e6 * cggHB * dKlambda * dytHB3 - 1.8438201632292444e8 * cggHB2 * dKlambda * dytHB3 - 8524.8765354653 * dcZHB * dKlambda * dytHB3 
            + 2.8867300035650665e6 * cggHB * dcZHB * dKlambda * dytHB3 - 9211.031646525304 * dcZHB2 * dKlambda * dytHB3 + 3263.1999469874036 * dKlambda2 * dytHB3 
            + 44138.45406924717 * cggHB * dKlambda2 * dytHB3 - 4.193837918690795e7 * cggHB2 * dKlambda2 * dytHB3 + 1474.023437403278 * dcZHB * dKlambda2 * dytHB3 
            + 322402.6653762193 * cggHB * dcZHB * dKlambda2 * dytHB3 + 116.36014794980927 * dKlambda3 * dytHB3 - 7370.4909474997985 * cggHB * dKlambda3 * dytHB3 
            - 3.4305355944930054e6 * cggHB2 * dKlambda3 * dytHB3 - 267.7615464146533 * dcZHB * dKlambda3 * dytHB3 + 84.29593960945068 * dKlambda4 * dytHB3 
            + 5735.3996792942135 * cggHB * dKlambda4 * dytHB3 + 66652.27308402126 * dytHB4 - 6.871040436399154e6 * cggHB * dytHB4 
            + 9.22099747455498e7 * cggHB2 * dytHB4 + 92021.78032189047 * dcZHB * dytHB4 - 2.257899878309953e6 * cggHB * dcZHB * dytHB4 
            + 16245.693309808961 * dcZHB2 * dytHB4 + 2838.4331580144003 * dKlambda * dytHB4 - 2.731422853592693e6 * cggHB * dKlambda * dytHB4 
            + 4.274439860749665e7 * cggHB2 * dKlambda * dytHB4 + 15892.926730807862 * dcZHB * dKlambda * dytHB4 - 515009.5486394962 * cggHB * dcZHB * dKlambda * dytHB4 
            - 1056.6073875703482 * dKlambda2 * dytHB4 - 482475.3464808796 * cggHB * dKlambda2 * dytHB4 + 5.170468004804585e6 * cggHB2 * dKlambda2 * dytHB4 
            + 2613.194223645355 * dcZHB * dKlambda2 * dytHB4 - 427.75818525652596 * dKlambda3 * dytHB4 - 51130.51778000078 * cggHB * dKlambda3 * dytHB4 
            + 21.07398490236267 * dKlambda4 * dytHB4 - 63203.969008703876 * dytHB5 + 3.151938475204292e6 * cggHB * dytHB5 - 42834.09620756765 * dcZHB * dytHB5 
            - 12524.979109927113 * dKlambda * dytHB5 + 1.3421161655790398e6 * cggHB * dKlambda * dytHB5 - 8919.930319126936 * dcZHB * dKlambda * dytHB5 
            - 849.49051561947 * dKlambda2 * dytHB5 + 158560.3321836832 * cggHB * dKlambda2 * dytHB5 - 263.0677528219873 * dKlambda3 * dytHB5 
            + 37913.4502786983 * dytHB6 - 712582.2268647491 * cggHB * dytHB6 + 10593.332328402174 * dcZHB * dytHB6 + 8514.598993531516 * dKlambda * dytHB6 
            - 169200.83566434312 * cggHB * dKlambda * dytHB6 + 1296.5492356304262 * dKlambda2 * dytHB6 - 13281.426292006341 * dytHB7 
            - 2976.898633587163 * dKlambda * dytHB7 + 2684.433665848417 * dytHB8;
    
    // To be used as Gaussian observable with mean=0, var=1 I must return the sqrt.
    return sqrt(Chi2Tot);
}

double NPSMEFTd6::AuxObs_NP15() const
{
    // To be used for some temporary observable
    return 0.0;
}

double NPSMEFTd6::AuxObs_NP16() const
{
    // To be used for some temporary observable
    return 0.0;
}

double NPSMEFTd6::AuxObs_NP17() const
{
    // To be used for some temporary observable
    return 0.0;
}

double NPSMEFTd6::AuxObs_NP18() const
{
    // To be used for some temporary observable
    return 0.0;
}

double NPSMEFTd6::AuxObs_NP19() const
{
    // To be used for some temporary observable
    return 0.0;
}

double NPSMEFTd6::AuxObs_NP20() const
{
    // To be used for some temporary observable
    return 0.0;
}

///////////////////////////////////////////////////////////////////////////////

double NPSMEFTd6::CLL_mu() const
{
    return (CLL_1122 + CLL_2211 + CiLL_1221 + CiLL_2112);
}

double NPSMEFTd6::CLL_tau() const
{
    return (CLL_1133 + CLL_3311 + CLL_1331 + CLL_3113);
}

double NPSMEFTd6::CLL_up() const
{
    return (CLQ1_1111-CLQ3_1111);
}

double NPSMEFTd6::CLL_down() const
{
    return (CLQ1_1111+CLQ3_1111);
}

double NPSMEFTd6::CLL_charm() const
{
    return (CLQ1_1122+CLQ1_2211-CLQ3_1122-CLQ3_2211);
}

double NPSMEFTd6::CLL_strange() const
{
    return (CLQ1_1122+CLQ1_2211+CLQ3_1122+CLQ3_2211);
}

double NPSMEFTd6::CLL_bottom() const
{
    return (CLQ1_1133+CLQ1_3311+CLQ3_1133+CLQ3_3311);
}

double NPSMEFTd6::CLR_mu() const
{
    return (CLe_1122+CLe_2211);
}

double NPSMEFTd6::CLR_tau() const
{
    return (CLe_1133+CLe_3311);
}

double NPSMEFTd6::CLR_up() const
{
    return (CLu_1111);
}

double NPSMEFTd6::CLR_down() const
{
    return (CLd_1111);
}

double NPSMEFTd6::CLR_charm() const
{
    return (CLu_1122+CLu_2211);
}

double NPSMEFTd6::CLR_strange() const
{
    return (CLd_1122+CLd_2211);
}

double NPSMEFTd6::CLR_bottom() const
{
    return (CLd_1133+CLd_3311);
}

double NPSMEFTd6::CRL_mu() const
{
    return (CLe_1122+CLe_2211);
}

double NPSMEFTd6::CRL_tau() const
{
    return (CLe_1133+CLe_3311);
}

double NPSMEFTd6::CRL_up() const
{
    return (CQe_1111);
}

double NPSMEFTd6::CRL_down() const
{
    return (CQe_1111);
}

double NPSMEFTd6::CRL_charm() const
{
    return (CQe_1122+CQe_2211);
}

double NPSMEFTd6::CRL_strange() const
{
    return (CQe_1122+CQe_2211);
}

double NPSMEFTd6::CRL_bottom() const
{
    return (CQe_1133+CQe_3311);
}

double NPSMEFTd6::CRR_mu() const
{
    return (Cee_1122+Cee_2211);
}

double NPSMEFTd6::CRR_tau() const
{
    return (Cee_1133+Cee_3311);
}


double NPSMEFTd6::CRR_up() const
{
    return (Ceu_1111);
}

double NPSMEFTd6::CRR_down() const
{
    return (Ced_1111);
}

double NPSMEFTd6::CRR_charm() const
{
    return (Ceu_1122+Ceu_2211);
}

double NPSMEFTd6::CRR_strange() const
{
    return (Ced_1122+Ced_2211);
}

double NPSMEFTd6::CRR_bottom() const
{
    return (Ced_1133+Ced_3311);
}
