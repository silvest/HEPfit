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
        = {"CG", "CW", "CHG", "CHW", "CHB", "CDHB", "CDHW", "CHWB", "CHD", "CHbox", "CH",
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
    "eHggint","eHggpar","eHWWint","eHWWpar","eHZZint","eHZZpar","eHZgaint","eHZgapar",
    "eHgagaint","eHgagapar","eHmumuint","eHmumupar","eHtautauint","eHtautaupar",
    "eHccint","eHccpar","eHbbint","eHbbpar",
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
        = {"CG", "CW", "CHG", "CHWHB_gaga", "CHWHB_gagaorth", "CDHB", "CDHW", "CHWB", "CHD", "CHbox", "CH",
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
    "eHggint","eHggpar","eHWWint","eHWWpar","eHZZint","eHZZpar","eHZgaint","eHZgapar",
    "eHgagaint","eHgagapar","eHmumuint","eHmumupar","eHtautauint","eHtautaupar",
    "eHccint","eHccpar","eHbbint","eHbbpar",
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
        = {"CG", "CW", "CHG", "CHW", "CHB", "CDHB", "CDHW", "CHWB", "CHD", "CHbox", "CH",
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
    "eHggint","eHggpar","eHWWint","eHWWpar","eHZZint","eHZZpar","eHZgaint","eHZgapar",
    "eHgagaint","eHgagapar","eHmumuint","eHmumupar","eHtautauint","eHtautaupar",
    "eHccint","eHccpar","eHbbint","eHbbpar",
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
        = {"CG", "CW", "CHG", "CHWHB_gaga", "CHWHB_gagaorth", "CDHB", "CDHW", "CHWB", "CHD", "CHbox", "CH",
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
    "eHggint","eHggpar","eHWWint","eHWWpar","eHZZint","eHZZpar","eHZgaint","eHZgapar",
    "eHgagaint","eHgagapar","eHmumuint","eHmumupar","eHtautauint","eHtautaupar",
    "eHccint","eHccpar","eHbbint","eHbbpar",
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
    FlagLoopHd6 = false;
    setModelLinearized();
    
    w_WW = gsl_integration_cquad_workspace_alloc(100);
    
    SMM.setObj((StandardModelMatching&) NPSMEFTd6M.getObj());

    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CG", boost::cref(CG)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CW", boost::cref(CW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHG", boost::cref(CHG)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHW", boost::cref(CHW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHB", boost::cref(CHB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHWHB_gaga", boost::cref(CHWHB_gaga)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHWHB_gagaorth", boost::cref(CHWHB_gagaorth)));        
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CDHB", boost::cref(CDHB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CDHW", boost::cref(CDHW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHWB", boost::cref(CHWB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHD", boost::cref(CHD)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHbox", boost::cref(CHbox)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CH", boost::cref(CH)));
    if (FlagLeptonUniversal) {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1", boost::cref(CHL1_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3", boost::cref(CHL3_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe", boost::cref(CHe_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_11r", boost::cref(CeH_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_22r", boost::cref(CeH_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_33r", boost::cref(CeH_33r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_11i", boost::cref(CeH_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_22i", boost::cref(CeH_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_33i", boost::cref(CeH_33i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLL", boost::cref(CLL_1221)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Cee", boost::cref(Cee_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLe", boost::cref(CLe_1111)));
    } else {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_11", boost::cref(CHL1_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_12r", boost::cref(CHL1_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_13r", boost::cref(CHL1_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_22", boost::cref(CHL1_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_23r", boost::cref(CHL1_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_33", boost::cref(CHL1_33)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_12i", boost::cref(CHL1_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_13i", boost::cref(CHL1_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_23i", boost::cref(CHL1_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_11", boost::cref(CHL3_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_12r", boost::cref(CHL3_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_13r", boost::cref(CHL3_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_22", boost::cref(CHL3_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_23r", boost::cref(CHL3_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_33", boost::cref(CHL3_33)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_12i", boost::cref(CHL3_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_13i", boost::cref(CHL3_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_23i", boost::cref(CHL3_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_11", boost::cref(CHe_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_12r", boost::cref(CHe_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_13r", boost::cref(CHe_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_22", boost::cref(CHe_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_23r", boost::cref(CHe_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_33", boost::cref(CHe_33)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_12i", boost::cref(CHe_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_13i", boost::cref(CHe_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_23i", boost::cref(CHe_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_11r", boost::cref(CeH_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_12r", boost::cref(CeH_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_13r", boost::cref(CeH_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_22r", boost::cref(CeH_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_23r", boost::cref(CeH_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_33r", boost::cref(CeH_33r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_11i", boost::cref(CeH_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_12i", boost::cref(CeH_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_13i", boost::cref(CeH_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_22i", boost::cref(CeH_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_23i", boost::cref(CeH_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_33i", boost::cref(CeH_33i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLL_1111", boost::cref(CLL_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLL_1221", boost::cref(CLL_1221)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLL_1122", boost::cref(CLL_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLL_1331", boost::cref(CLL_1331)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLL_1133", boost::cref(CLL_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Cee_1111", boost::cref(Cee_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Cee_1122", boost::cref(Cee_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Cee_1133", boost::cref(Cee_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLe_1111", boost::cref(CLe_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLe_1122", boost::cref(CLe_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLe_2211", boost::cref(CLe_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLe_1133", boost::cref(CLe_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLe_3311", boost::cref(CLe_3311)));       
    }
    if (FlagQuarkUniversal) {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1", boost::cref(CHQ1_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3", boost::cref(CHQ3_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu", boost::cref(CHu_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd", boost::cref(CHd_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_r", boost::cref(CHud_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_i", boost::cref(CHud_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_11r", boost::cref(CuH_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_22r", boost::cref(CuH_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_33r", boost::cref(CuH_33r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_11i", boost::cref(CuH_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_22i", boost::cref(CuH_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_33i", boost::cref(CuH_33i)));        
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_11r", boost::cref(CdH_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_22r", boost::cref(CdH_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_33r", boost::cref(CdH_33r)));        
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_11i", boost::cref(CdH_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_22i", boost::cref(CdH_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_33i", boost::cref(CdH_33i)));        
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_r", boost::cref(CuG_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_i", boost::cref(CuG_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_r", boost::cref(CuW_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_i", boost::cref(CuW_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_r", boost::cref(CuB_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_i", boost::cref(CuB_11i)));
    } else {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_11", boost::cref(CHQ1_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_12r", boost::cref(CHQ1_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_13r", boost::cref(CHQ1_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_22", boost::cref(CHQ1_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_23r", boost::cref(CHQ1_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_33", boost::cref(CHQ1_33)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_12i", boost::cref(CHQ1_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_13i", boost::cref(CHQ1_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_23i", boost::cref(CHQ1_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_11", boost::cref(CHQ3_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_12r", boost::cref(CHQ3_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_13r", boost::cref(CHQ3_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_22", boost::cref(CHQ3_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_23r", boost::cref(CHQ3_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_33", boost::cref(CHQ3_33)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_12i", boost::cref(CHQ3_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_13i", boost::cref(CHQ3_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_23i", boost::cref(CHQ3_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_11", boost::cref(CHu_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_12r", boost::cref(CHu_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_13r", boost::cref(CHu_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_22", boost::cref(CHu_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_23r", boost::cref(CHu_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_33", boost::cref(CHu_33)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_12i", boost::cref(CHu_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_13i", boost::cref(CHu_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_23i", boost::cref(CHu_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_11", boost::cref(CHd_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_12r", boost::cref(CHd_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_13r", boost::cref(CHd_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_22", boost::cref(CHd_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_23r", boost::cref(CHd_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_33", boost::cref(CHd_33)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_12i", boost::cref(CHd_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_13i", boost::cref(CHd_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_23i", boost::cref(CHd_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_11r", boost::cref(CHud_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_12r", boost::cref(CHud_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_13r", boost::cref(CHud_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_22r", boost::cref(CHud_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_23r", boost::cref(CHud_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_33r", boost::cref(CHud_33r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_11i", boost::cref(CHud_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_12i", boost::cref(CHud_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_13i", boost::cref(CHud_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_22i", boost::cref(CHud_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_23i", boost::cref(CHud_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_33i", boost::cref(CHud_33i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_11r", boost::cref(CuH_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_12r", boost::cref(CuH_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_13r", boost::cref(CuH_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_22r", boost::cref(CuH_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_23r", boost::cref(CuH_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_33r", boost::cref(CuH_33r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_11i", boost::cref(CuH_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_12i", boost::cref(CuH_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_13i", boost::cref(CuH_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_22i", boost::cref(CuH_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_23i", boost::cref(CuH_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_33i", boost::cref(CuH_33i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_11r", boost::cref(CdH_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_12r", boost::cref(CdH_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_13r", boost::cref(CdH_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_22r", boost::cref(CdH_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_23r", boost::cref(CdH_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_33r", boost::cref(CdH_33r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_11i", boost::cref(CdH_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_12i", boost::cref(CdH_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_13i", boost::cref(CdH_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_22i", boost::cref(CdH_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_23i", boost::cref(CdH_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_33i", boost::cref(CdH_33i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_11r", boost::cref(CuG_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_12r", boost::cref(CuG_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_13r", boost::cref(CuG_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_22r", boost::cref(CuG_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_23r", boost::cref(CuG_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_33r", boost::cref(CuG_33r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_11i", boost::cref(CuG_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_12i", boost::cref(CuG_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_13i", boost::cref(CuG_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_22i", boost::cref(CuG_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_23i", boost::cref(CuG_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_33i", boost::cref(CuG_33i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_11r", boost::cref(CuW_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_12r", boost::cref(CuW_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_13r", boost::cref(CuW_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_22r", boost::cref(CuW_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_23r", boost::cref(CuW_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_33r", boost::cref(CuW_33r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_11i", boost::cref(CuW_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_12i", boost::cref(CuW_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_13i", boost::cref(CuW_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_22i", boost::cref(CuW_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_23i", boost::cref(CuW_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_33i", boost::cref(CuW_33i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_11r", boost::cref(CuB_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_12r", boost::cref(CuB_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_13r", boost::cref(CuB_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_22r", boost::cref(CuB_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_23r", boost::cref(CuB_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_33r", boost::cref(CuB_33r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_11i", boost::cref(CuB_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_12i", boost::cref(CuB_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_13i", boost::cref(CuB_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_22i", boost::cref(CuB_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_23i", boost::cref(CuB_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_33i", boost::cref(CuB_33i)));
    }
    if(FlagLeptonUniversal && FlagQuarkUniversal){
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1", boost::cref(CLQ1_1111)));  
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3", boost::cref(CLQ3_1111))); 
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ceu", boost::cref(Ceu_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced", boost::cref(Ced_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLu", boost::cref(CLu_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd", boost::cref(CLd_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe", boost::cref(CQe_1111)));
    } else {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_1111", boost::cref(CLQ1_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_1122", boost::cref(CLQ1_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_2211", boost::cref(CLQ1_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_1221", boost::cref(CLQ1_1221)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_2112", boost::cref(CLQ1_2112)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_1133", boost::cref(CLQ1_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_3311", boost::cref(CLQ1_3311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_1331", boost::cref(CLQ1_1331)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_3113", boost::cref(CLQ1_3113)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_1123", boost::cref(CLQ1_1123)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_2223", boost::cref(CLQ1_2223)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_3323", boost::cref(CLQ1_3323)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_1132", boost::cref(CLQ1_1132)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_2232", boost::cref(CLQ1_2232)));       
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_3332", boost::cref(CLQ1_3332)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_1111", boost::cref(CLQ3_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_1122", boost::cref(CLQ3_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_2211", boost::cref(CLQ3_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_1221", boost::cref(CLQ3_1221)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_2112", boost::cref(CLQ3_2112)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_1133", boost::cref(CLQ3_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_3311", boost::cref(CLQ3_3311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_1331", boost::cref(CLQ3_1331)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_3113", boost::cref(CLQ3_3113)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_1123", boost::cref(CLQ3_1123)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_2223", boost::cref(CLQ3_2223)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_3323", boost::cref(CLQ3_3323)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_1132", boost::cref(CLQ3_1132)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_2232", boost::cref(CLQ3_2232)));       
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_3332", boost::cref(CLQ3_3332)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ceu_1111", boost::cref(Ceu_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ceu_1122", boost::cref(Ceu_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ceu_2211", boost::cref(Ceu_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ceu_1133", boost::cref(Ceu_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ceu_2233", boost::cref(Ceu_2233)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ceu_3311", boost::cref(Ceu_3311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_1111", boost::cref(Ced_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_1122", boost::cref(Ced_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_2211", boost::cref(Ced_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_1133", boost::cref(Ced_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_3311", boost::cref(Ced_3311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_1123", boost::cref(Ced_1123)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_2223", boost::cref(Ced_2223)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_3323", boost::cref(Ced_3323)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_1132", boost::cref(Ced_1132)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_2232", boost::cref(Ced_2232)));       
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_3332", boost::cref(Ced_3332)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLu_1111", boost::cref(CLu_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLu_1122", boost::cref(CLu_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLu_2211", boost::cref(CLu_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLu_1133", boost::cref(CLu_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLu_2233", boost::cref(CLu_2233)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLu_3311", boost::cref(CLu_3311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_1111", boost::cref(CLd_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_1122", boost::cref(CLd_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_2211", boost::cref(CLd_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_1133", boost::cref(CLd_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_3311", boost::cref(CLd_3311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_1123", boost::cref(CLd_1123)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_2223", boost::cref(CLd_2223)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_3323", boost::cref(CLd_3323)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_1132", boost::cref(CLd_1132)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_2232", boost::cref(CLd_2232)));       
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_3332", boost::cref(CLd_3332)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_1111", boost::cref(CQe_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_1122", boost::cref(CQe_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_2211", boost::cref(CQe_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_1133", boost::cref(CQe_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_3311", boost::cref(CQe_3311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_2311", boost::cref(CQe_2311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_2322", boost::cref(CQe_2322)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_2333", boost::cref(CQe_2333)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_3211", boost::cref(CQe_3211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_3222", boost::cref(CQe_3222)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_3233", boost::cref(CQe_3233)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLedQ_11", boost::cref(CLedQ_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLedQ_22", boost::cref(CLedQ_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CpLedQ_11", boost::cref(CpLedQ_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CpLedQ_22", boost::cref(CpLedQ_22)));
    }
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Lambda_NP", boost::cref(Lambda_NP)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BrHinv", boost::cref(BrHinv)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BrHexo", boost::cref(BrHexo)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("dg1Z", boost::cref(dg1Z)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("dKappaga", boost::cref(dKappaga)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("lambZ", boost::cref(lambZ)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eggFint", boost::cref(eggFint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eggFpar", boost::cref(eggFpar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettHint", boost::cref(ettHint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettHpar", boost::cref(ettHpar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBFint", boost::cref(eVBFint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBFpar", boost::cref(eVBFpar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWHint", boost::cref(eWHint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWHpar", boost::cref(eWHpar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZHint", boost::cref(eZHint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZHpar", boost::cref(eZHpar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eeeWBFint", boost::cref(eeeWBFint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eeeWBFpar", boost::cref(eeeWBFpar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eeeZHint", boost::cref(eeeZHint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eeeZHpar", boost::cref(eeeZHpar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eeettHint", boost::cref(eeettHint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eeettHpar", boost::cref(eeettHpar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHggint", boost::cref(eHggint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHggpar", boost::cref(eHggpar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHWWint", boost::cref(eHWWint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHWWpar", boost::cref(eHWWpar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHZZint", boost::cref(eHZZint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHZZpar", boost::cref(eHZZpar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHZgaint", boost::cref(eHZgaint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHZgapar", boost::cref(eHZgapar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHgagaint", boost::cref(eHgagaint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHgagapar", boost::cref(eHgagapar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHmumuint", boost::cref(eHmumuint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHmumupar", boost::cref(eHmumupar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHtautauint", boost::cref(eHtautauint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHtautaupar", boost::cref(eHtautaupar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHccint", boost::cref(eHccint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHccpar", boost::cref(eHccpar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHbbint", boost::cref(eHbbint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHbbpar", boost::cref(eHbbpar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_2_Hbox", boost::cref(eVBF_2_Hbox)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_2_HQ1_11", boost::cref(eVBF_2_HQ1_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_2_Hu_11", boost::cref(eVBF_2_Hu_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_2_Hd_11", boost::cref(eVBF_2_Hd_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_2_HQ3_11", boost::cref(eVBF_2_HQ3_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_2_HD", boost::cref(eVBF_2_HD)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_2_HB", boost::cref(eVBF_2_HB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_2_HW", boost::cref(eVBF_2_HW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_2_HWB", boost::cref(eVBF_2_HWB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_2_HG", boost::cref(eVBF_2_HG)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_2_DHB", boost::cref(eVBF_2_DHB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_2_DHW", boost::cref(eVBF_2_DHW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_2_DeltaGF", boost::cref(eVBF_2_DeltaGF)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_78_Hbox", boost::cref(eVBF_78_Hbox)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_78_HQ1_11", boost::cref(eVBF_78_HQ1_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_78_Hu_11", boost::cref(eVBF_78_Hu_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_78_Hd_11", boost::cref(eVBF_78_Hd_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_78_HQ3_11", boost::cref(eVBF_78_HQ3_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_78_HD", boost::cref(eVBF_78_HD)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_78_HB", boost::cref(eVBF_78_HB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_78_HW", boost::cref(eVBF_78_HW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_78_HWB", boost::cref(eVBF_78_HWB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_78_HG", boost::cref(eVBF_78_HG)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_78_DHB", boost::cref(eVBF_78_DHB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_78_DHW", boost::cref(eVBF_78_DHW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_78_DeltaGF", boost::cref(eVBF_78_DeltaGF)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_1314_Hbox", boost::cref(eVBF_1314_Hbox)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_1314_HQ1_11", boost::cref(eVBF_1314_HQ1_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_1314_Hu_11", boost::cref(eVBF_1314_Hu_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_1314_Hd_11", boost::cref(eVBF_1314_Hd_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_1314_HQ3_11", boost::cref(eVBF_1314_HQ3_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_1314_HD", boost::cref(eVBF_1314_HD)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_1314_HB", boost::cref(eVBF_1314_HB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_1314_HW", boost::cref(eVBF_1314_HW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_1314_HWB", boost::cref(eVBF_1314_HWB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_1314_HG", boost::cref(eVBF_1314_HG)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_1314_DHB", boost::cref(eVBF_1314_DHB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_1314_DHW", boost::cref(eVBF_1314_DHW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_1314_DeltaGF", boost::cref(eVBF_1314_DeltaGF)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_2_Hbox", boost::cref(eWH_2_Hbox)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_2_HQ3_11", boost::cref(eWH_2_HQ3_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_2_HD", boost::cref(eWH_2_HD)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_2_HW", boost::cref(eWH_2_HW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_2_HWB", boost::cref(eWH_2_HWB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_2_DHW", boost::cref(eWH_2_DHW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_2_DeltaGF", boost::cref(eWH_2_DeltaGF)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_78_Hbox", boost::cref(eWH_78_Hbox)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_78_HQ3_11", boost::cref(eWH_78_HQ3_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_78_HD", boost::cref(eWH_78_HD)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_78_HW", boost::cref(eWH_78_HW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_78_HWB", boost::cref(eWH_78_HWB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_78_DHW", boost::cref(eWH_78_DHW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_78_DeltaGF", boost::cref(eWH_78_DeltaGF)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_1314_Hbox", boost::cref(eWH_1314_Hbox)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_1314_HQ3_11", boost::cref(eWH_1314_HQ3_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_1314_HD", boost::cref(eWH_1314_HD)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_1314_HW", boost::cref(eWH_1314_HW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_1314_HWB", boost::cref(eWH_1314_HWB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_1314_DHW", boost::cref(eWH_1314_DHW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_1314_DeltaGF", boost::cref(eWH_1314_DeltaGF)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_2_Hbox", boost::cref(eZH_2_Hbox)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_2_HQ1_11", boost::cref(eZH_2_HQ1_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_2_Hu_11", boost::cref(eZH_2_Hu_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_2_Hd_11", boost::cref(eZH_2_Hd_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_2_HQ3_11", boost::cref(eZH_2_HQ3_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_2_HD", boost::cref(eZH_2_HD)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_2_HB", boost::cref(eZH_2_HB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_2_HW", boost::cref(eZH_2_HW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_2_HWB", boost::cref(eZH_2_HWB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_2_DHB", boost::cref(eZH_2_DHB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_2_DHW", boost::cref(eZH_2_DHW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_2_DeltaGF", boost::cref(eZH_2_DeltaGF)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_78_Hbox", boost::cref(eZH_78_Hbox)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_78_HQ1_11", boost::cref(eZH_78_HQ1_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_78_Hu_11", boost::cref(eZH_78_Hu_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_78_Hd_11", boost::cref(eZH_78_Hd_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_78_HQ3_11", boost::cref(eZH_78_HQ3_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_78_HD", boost::cref(eZH_78_HD)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_78_HB", boost::cref(eZH_78_HB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_78_HW", boost::cref(eZH_78_HW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_78_HWB", boost::cref(eZH_78_HWB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_78_DHB", boost::cref(eZH_78_DHB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_78_DHW", boost::cref(eZH_78_DHW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_78_DeltaGF", boost::cref(eZH_78_DeltaGF)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_1314_Hbox", boost::cref(eZH_1314_Hbox)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_1314_HQ1_11", boost::cref(eZH_1314_HQ1_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_1314_Hu_11", boost::cref(eZH_1314_Hu_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_1314_Hd_11", boost::cref(eZH_1314_Hd_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_1314_HQ3_11", boost::cref(eZH_1314_HQ3_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_1314_HD", boost::cref(eZH_1314_HD)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_1314_HB", boost::cref(eZH_1314_HB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_1314_HW", boost::cref(eZH_1314_HW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_1314_HWB", boost::cref(eZH_1314_HWB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_1314_DHB", boost::cref(eZH_1314_DHB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_1314_DHW", boost::cref(eZH_1314_DHW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_1314_DeltaGF", boost::cref(eZH_1314_DeltaGF)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH_2_HG", boost::cref(ettH_2_HG)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH_2_G", boost::cref(ettH_2_G)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH_2_uG_33r", boost::cref(ettH_2_uG_33r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH_2_DeltagHt", boost::cref(ettH_2_DeltagHt)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH_78_HG", boost::cref(ettH_78_HG)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH_78_G", boost::cref(ettH_78_G)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH_78_uG_33r", boost::cref(ettH_78_uG_33r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH_78_DeltagHt", boost::cref(ettH_78_DeltagHt)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH_1314_HG", boost::cref(ettH_1314_HG)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH_1314_G", boost::cref(ettH_1314_G)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH_1314_uG_33r", boost::cref(ettH_1314_uG_33r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH_1314_DeltagHt", boost::cref(ettH_1314_DeltagHt)));
    
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

    LambdaNP2 = Lambda_NP * Lambda_NP;
    v2_over_LambdaNP2 = v() * v() / LambdaNP2;
    aleMz = alphaMz();
    cW_tree = Mw_tree() / Mz;
    cW2_tree = cW_tree * cW_tree;
    sW2_tree = 1.0 - cW2_tree;
    sW_tree = sqrt(sW2_tree);

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
      
    if (FlagRotateCHWCHB) {
        CHW = sW2_tree * CHWHB_gaga - cW2_tree * CHWHB_gagaorth;
        CHB = cW2_tree * CHWHB_gaga + sW2_tree * CHWHB_gagaorth;
    } else {
        CHWHB_gaga = sW2_tree * CHW + cW2_tree * CHB;
        CHWHB_gagaorth = - cW2_tree * CHW + sW2_tree * CHB;
    }
    
    if (FlagFlavU3OfX) {

        CeH_11r = Yuke * CeH_11r;
        CeH_22r = Yukmu * CeH_22r;        
        CeH_33r = Yuktau * CeH_33r;
                
        CuH_11r = Yuku * CuH_11r;            
        CuH_22r = Yukc * CuH_22r;     
        CuH_33r = Yukt * CuH_33r;
        
        CdH_11r = Yukd * CdH_11r;
        CdH_22r = Yuks * CdH_22r;            
        CdH_33r = Yukb * CdH_33r;
        
        CuG_11r = Yuku * CuG_11r;  
        CuG_22r = Yukc * CuG_22r;
        CuG_33r = Yukt * CuG_33r;

        CuW_11r = Yuku * CuW_11r;
        CuW_22r = Yukc * CuW_22r;
        CuW_33r = Yukt * CuW_33r;

        CuB_11r = Yuku * CuB_11r;
        CuB_22r = Yukc * CuB_22r;                   
        CuB_33r = Yukt * CuB_33r;
    }
    
    if (!FlagLoopHd6) {
        cLHd6 = 0.;
    } else {
        cLHd6 = 1.;
    }

    delta_ZZ = (cW2_tree * CHW + sW2_tree * CHB + sW_tree * cW_tree * CHWB) * v2_over_LambdaNP2;
    delta_AA = (sW2_tree * CHW + cW2_tree * CHB - sW_tree * cW_tree * CHWB) * v2_over_LambdaNP2;
    delta_AZ = 2.0 * sW_tree * cW_tree * (CHW - CHB) * v2_over_LambdaNP2
            - (cW2_tree - sW2_tree) * CHWB * v2_over_LambdaNP2;
    delta_h = (-CHD / 4.0 + CHbox) * v2_over_LambdaNP2;
    
    NPSMEFTd6M.getObj().updateNPSMEFTd6Parameters();

    return (true);
}

void NPSMEFTd6::setParameter(const std::string name, const double& value)
{
    if (name.compare("CG") == 0)
        CG = value;
    else if (name.compare("CW") == 0)
        CW = value;
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
    else if (name.compare("CHWB") == 0)
        CHWB = value;
    else if (name.compare("CHD") == 0)
        CHD = value;
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
                          return false;
                    }
              }
        } else {
              for (int i = 0; i < NNPSMEFTd6Vars_LFU_QFU; i++) {
                    if (DPars.find(NPSMEFTd6Vars_LFU_QFU[i]) == DPars.end()) {
                          std::cout << "ERROR: Missing mandatory NPSMEFTd6_LFU_QFU parameter "
                          << NPSMEFTd6Vars_LFU_QFU[i] << std::endl;
                          return false;
                    }
              }
        }
          
        //} else if (FlagLeptonUniversal && !FlagQuarkUniversal) {
        //} else if (!FlagLeptonUniversal && FlagQuarkUniversal) {
    } else if (!FlagLeptonUniversal && !FlagQuarkUniversal) {
        if (FlagRotateCHWCHB) {
              for (int i = 0; i < NNPSMEFTd6Vars; i++) {
                    if (DPars.find(NPSMEFTd6VarsRot[i]) == DPars.end()) {
                          std::cout << "ERROR: Missing mandatory NPSMEFTd6 parameter"
                          << NPSMEFTd6VarsRot[i] << std::endl;
                          return false;
                    }
              }
        } else {
              for (int i = 0; i < NNPSMEFTd6Vars; i++) {
                    if (DPars.find(NPSMEFTd6Vars[i]) == DPars.end()) {
                          std::cout << "ERROR: Missing mandatory NPSMEFTd6 parameter"
                          << NPSMEFTd6Vars[i] << std::endl;
                          return false;
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
    } else if (name.compare("LoopHd6") == 0) {        
        FlagLoopHd6 = value;
        res = true;
    } else
        res = NPbase::setFlag(name, value);

    return (res);
}


////////////////////////////////////////////////////////////////////////

double NPSMEFTd6::CHF1_diag(const Particle F) const
{
    if (F.is("NEUTRINO_1") || F.is("ELECTRON"))
        return CHL1_11;
    else if (F.is("NEUTRINO_2") || F.is("MU"))
        return CHL1_22;
    else if (F.is("NEUTRINO_3") || F.is("TAU"))
        return CHL1_33;
    else if (F.is("UP") || F.is("DOWN"))
        return CHQ1_11;
    else if (F.is("CHARM") || F.is("STRANGE"))
        return CHQ1_22;
    else if (F.is("TOP") || F.is("BOTTOM"))
        return CHQ1_33;
    else
        throw std::runtime_error("NPSMEFTd6::CHF1_diag(): wrong argument");
}

double NPSMEFTd6::CHF3_diag(const Particle F) const
{
    if (F.is("NEUTRINO_1") || F.is("ELECTRON"))
        return CHL3_11;
    else if (F.is("NEUTRINO_2") || F.is("MU"))
        return CHL3_22;
    else if (F.is("NEUTRINO_3") || F.is("TAU"))
        return CHL3_33;
    else if (F.is("UP") || F.is("DOWN"))
        return CHQ3_11;
    else if (F.is("CHARM") || F.is("STRANGE"))
        return CHQ3_22;
    else if (F.is("TOP") || F.is("BOTTOM"))
        return CHQ3_33;
    else
        throw std::runtime_error("NPSMEFTd6::CHF3_diag(): wrong argument");
}

double NPSMEFTd6::CHf_diag(const Particle f) const
{
    if (f.is("NEUTRINO_1") || f.is("NEUTRINO_2") || f.is("NEUTRINO_3"))
        return 0.0;
    else if (f.is("ELECTRON"))
        return CHe_11;
    else if (f.is("MU"))
        return CHe_22;
    else if (f.is("TAU"))
        return CHe_33;
    else if (f.is("UP"))
        return CHu_11;
    else if (f.is("CHARM"))
        return CHu_22;
    else if (f.is("TOP"))
        return CHu_33;
    else if (f.is("DOWN"))
        return CHd_11;
    else if (f.is("STRANGE"))
        return CHd_22;
    else if (f.is("BOTTOM"))
        return CHd_33;
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
        return gslpp::complex(CeH_11r, CeH_11i, false);
    else if (f.is("MU"))
        return gslpp::complex(CeH_22r, CeH_22i, false);
    else if (f.is("TAU"))
        return gslpp::complex(CeH_33r, CeH_33i, false);
    else if (f.is("UP"))
        return gslpp::complex(CuH_11r, CuH_11i, false);
    else if (f.is("CHARM"))
        return gslpp::complex(CuH_22r, CuH_22i, false);
    else if (f.is("TOP"))
        return gslpp::complex(CuH_33r, CuH_33i, false);
    else if (f.is("DOWN"))
        return gslpp::complex(CdH_11r, CdH_11i, false);
    else if (f.is("STRANGE"))
        return gslpp::complex(CdH_22r, CdH_22i, false);
    else if (f.is("BOTTOM"))
        return gslpp::complex(CdH_33r, CdH_33i, false);
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
        return gslpp::complex(CuG_11r, CuG_11i, false);
    else if (f.is("CHARM"))
        return gslpp::complex(CuG_22r, CuG_22i, false);
    else if (f.is("TOP"))
        return gslpp::complex(CuG_33r, CuG_33i, false);
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
        return gslpp::complex(CuW_11r, CuW_11i, false);
    else if (f.is("CHARM"))
        return gslpp::complex(CuW_22r, CuW_22i, false);
    else if (f.is("TOP"))
        return gslpp::complex(CuW_33r, CuW_33i, false);
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
        return gslpp::complex(CuB_11r, CuB_11i, false);
    else if (f.is("CHARM"))
        return gslpp::complex(CuB_22r, CuB_22i, false);
    else if (f.is("TOP"))
        return gslpp::complex(CuB_33r, CuB_33i, false);
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
    return ((CHL3_11 + CHL3_22 - 0.5 * (CLL_1221 + CLL_2112)) * v2_over_LambdaNP2);
}

double NPSMEFTd6::obliqueS() const
{
    return (4.0 * sW_tree * cW_tree * CHWB / alphaMz() * v2_over_LambdaNP2);
}

double NPSMEFTd6::obliqueT() const
{
    return (-CHD / 2.0 / alphaMz() * v2_over_LambdaNP2);
}

double NPSMEFTd6::obliqueU() const
{
    return 0.0;
}

double NPSMEFTd6::deltaMz2() const
{
    return (2.0*(cW2_tree*CHW+sW2_tree*CHB+cW_tree*sW_tree*CHWB)*v2_over_LambdaNP2);
}

double NPSMEFTd6::Mw() const
{
    return (trueSM.Mw() - Mw_tree() / 4.0 / (cW2_tree - sW2_tree)
            *(4.0 * sW_tree * cW_tree * CHWB * v2_over_LambdaNP2
            + cW2_tree * CHD * v2_over_LambdaNP2
            + 2.0 * sW2_tree * DeltaGF()));
}

double NPSMEFTd6::deltaMw() const
{
    return (- 1.0 / 4.0 / (cW2_tree - sW2_tree)
            *(4.0 * sW_tree * cW_tree * CHWB * v2_over_LambdaNP2
            + cW2_tree * CHD * v2_over_LambdaNP2
            + 2.0 * sW2_tree * DeltaGF()));
}

double NPSMEFTd6::deltaMw2() const
{
    double dMW = deltaMw();
    
    return (dMW*dMW);
}

double NPSMEFTd6::deltaGamma_Wff(const Particle fi, const Particle fj) const
{
    double G0 = GF * pow(Mw(), 3.0) / 6.0 / sqrt(2.0) / M_PI;
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
            *(4.0 * sW_tree * cW_tree * CHWB * v2_over_LambdaNP2
            + cW2_tree * CHD * v2_over_LambdaNP2
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
    double G0 = GF * pow(Mw(), 3.0) / 6.0 / sqrt(2.0) / M_PI;
    double GammaW_tree = (3.0 + 2.0 * Nc) * G0;

    return (- 3.0 * GammaW_tree / 4.0 / (cW2_tree - sW2_tree)
            *(4.0 * sW_tree * cW_tree * CHWB * v2_over_LambdaNP2
            + cW2_tree * CHD * v2_over_LambdaNP2
            + 2.0 * (1.0 + cW2_tree) / 3.0 * DeltaGF())
            + 2.0 * G0 * (CHL3_11 + CHL3_22 + CHL3_33 + Nc*(CHQ3_11 + CHQ3_22)) * v2_over_LambdaNP2);          
//            + 2.0 * GammaW_tree / 3.0 * (CHL3_11 + CHQ3_11 + CHQ3_22) * v2_over_LambdaNP2);    
}

double NPSMEFTd6::GammaW() const
{
    return ( trueSM.GammaW() + deltaGamma_W() );
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

    NPindirect = -I3p / 4.0 * (CHD * v2_over_LambdaNP2 + 2.0 * DeltaGF())
                - Qp * sW2_tree / 4.0 / (cW2_tree - sW2_tree)
                *((4.0 * cW_tree / sW_tree * CHWB + CHD) * v2_over_LambdaNP2 + 2.0 * DeltaGF());

    double NPdirect = -0.5 * (CHF1 - 2.0 * I3p * CHF3) * v2_over_LambdaNP2;
    return (NPindirect + NPdirect);
}

double NPSMEFTd6::deltaGR_f(const Particle p) const
{
    double Qp = p.getCharge();
    double CHf = CHf_diag(p);
    double NPindirect;

    NPindirect = -Qp * sW2_tree / 4.0 / (cW2_tree - sW2_tree)
                *((4.0 * cW_tree / sW_tree * CHWB + CHD) * v2_over_LambdaNP2 + 2.0 * DeltaGF());

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
                * ((4.0 * sW_tree / cW_tree * CHWB + CHD) * v2_over_LambdaNP2 + 2.0 * DeltaGF());

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
    return (( 2.0 * CHW - sqrt( M_PI * aleMz ) * CDHW / sW_tree ) * v2_over_LambdaNP2 / v());
}

double NPSMEFTd6::deltaG2_hWW() const
{
    return ( - sqrt( M_PI * aleMz ) * ( CDHW / sW_tree ) * v2_over_LambdaNP2 / v());
}

double NPSMEFTd6::deltaG3_hWW() const
{
    double NPindirect;

    NPindirect = 2.0 * cW2_tree * Mz * Mz / v()
                * (delta_h - 1.0 / 2.0 / (cW2_tree - sW2_tree)
                * ((4.0 * sW_tree * cW_tree * CHWB + cW2_tree * CHD) * v2_over_LambdaNP2 + DeltaGF()));

    return NPindirect;
}

double NPSMEFTd6::deltaG1_hZZ() const
{
    return ( (delta_ZZ - 0.5 * sqrt( M_PI * aleMz ) * (CDHB / cW_tree + CDHW / sW_tree) * v2_over_LambdaNP2 )/ v());
}

double NPSMEFTd6::deltaG2_hZZ() const
{
    return ( - sqrt( M_PI * aleMz ) * ( CDHB / cW_tree + CDHW / sW_tree ) * v2_over_LambdaNP2 / v());
}

double NPSMEFTd6::deltaG3_hZZ() const
{
    double NPindirect = Mz * Mz / v() * (-0.5 * CHD * v2_over_LambdaNP2 + delta_h - 0.5 * DeltaGF());
    double NPdirect = Mz * Mz / v() * CHD * v2_over_LambdaNP2;
    return (NPindirect + NPdirect);
}

double NPSMEFTd6::deltaG1_hZA() const
{
    return ( (delta_AZ + 0.5 * sqrt( M_PI * aleMz ) * (CDHB / sW_tree - CDHW / cW_tree) * v2_over_LambdaNP2 )/ v());
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
    return ( sqrt( M_PI * aleMz ) * ( CDHB / sW_tree - CDHW / cW_tree ) * v2_over_LambdaNP2 / v());
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
    
    dg = -0.5 * DeltaGF() + 3.0 * delta_h - 2.0 * CH * v2_over_LambdaNP2 * v()*v()/mHl/mHl;

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
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        gslpp::complex tmp2 = tmpt +tmpb +tmpc + tmpHG;
        
        mu += tmp2.abs2();
        
//  Quadratic contribution from Higgs self-coupling
    mu = mu + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
        
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

        A1HH = 2.13;
        A2HH = 10.1;
        A3HH = 0.300;
        A4HH = 21.8;
        A5HH = 188;
        A6HH = -8.62;
        A7HH = -1.43;
        A8HH = 2.93;
        A9HH = 21.0;
        A10HH = 59.8;
        A11HH = -9.93;
        A12HH = -23.1;
        A13HH = 4.87;
        A14HH = 10.5;
        A15HH = 96.6;
        
    } else if (sqrt_s == 100.0) {

        A1HH = 1.95;
        A2HH = 11.2;
        A3HH = 0.229;
        A4HH = 16.0;
        A5HH = 386;
        A6HH = -8.32;
        A7HH = -1.18;
        A8HH = 2.55;
        A9HH = 16.9;
        A10HH = 52.4;
        A11HH = -7.49;
        A12HH = -17.3;
        A13HH = 3.55;
        A14HH = 8.46;
        A15HH = 87.5;
        
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::muggHH()");
    
    ct= 1.0 - 0.5 * DeltaGF() + delta_h - v() * CuH_33r * v2_over_LambdaNP2 / sqrt(2.0)/ mtpole;
    c2t = delta_h - 3.0 *v() * CuH_33r * v2_over_LambdaNP2 / 2.0 /sqrt(2.0)/ mtpole;
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
                +120936. * (1. + eVBF_2_Hbox ) * CHbox / LambdaNP2
                -9422.68 * (1. + eVBF_2_HQ1_11 ) * CHQ1_11 / LambdaNP2
                -10683.8 * (1. + eVBF_2_Hu_11 ) * CHu_11 / LambdaNP2
                +4055.59 * (1. + eVBF_2_Hd_11 ) * CHd_11 / LambdaNP2
                -229691. * (1. + eVBF_2_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -170093. * (1. + eVBF_2_HD ) * CHD / LambdaNP2
                +8971.22 * (1. + eVBF_2_HB ) * CHB / LambdaNP2
                -65827.6 * (1. + eVBF_2_HW ) * CHW / LambdaNP2
                -323514. * (1. + eVBF_2_HWB ) * CHWB / LambdaNP2
                +481332. * (1. + eVBF_2_HG ) * CHG / LambdaNP2
                +1255.16 * (1. + eVBF_2_DHB ) * CDHB / LambdaNP2
                -34956.7 * (1. + eVBF_2_DHW ) * CDHW / LambdaNP2
                -4.511 * (1. + eVBF_2_DeltaGF ) * DeltaGF()
                -3.481 * deltaMw()
                +45.398 * deltaMw2()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

        mu +=  0.0;

        }
        
    } else if (sqrt_s == 7.0) {
        
        C1 = 0.0065;

        mu += 
                +121582. * (1. + eVBF_78_Hbox ) * CHbox / LambdaNP2
                +13546.6 * (1. + eVBF_78_HQ1_11 ) * CHQ1_11 / LambdaNP2
                -27657.6 * (1. + eVBF_78_Hu_11 ) * CHu_11 / LambdaNP2
                +8892.12 * (1. + eVBF_78_Hd_11 ) * CHd_11 / LambdaNP2
                -411400. * (1. + eVBF_78_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -164286. * (1. + eVBF_78_HD ) * CHD / LambdaNP2
                -423.123 * (1. + eVBF_78_HB ) * CHB / LambdaNP2
                -89854. * (1. + eVBF_78_HW ) * CHW / LambdaNP2
                -312617. * (1. + eVBF_78_HWB ) * CHWB / LambdaNP2
                -82956.8 * (1. + eVBF_78_HG ) * CHG / LambdaNP2
                -279.08 * (1. + eVBF_78_DHB ) * CDHB / LambdaNP2
                -54861. * (1. + eVBF_78_DHW ) * CDHW / LambdaNP2
                -4.479 * (1. + eVBF_78_DeltaGF ) * DeltaGF()
                -3.22 * deltaMw()
                +36.828 * deltaMw2()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

        mu +=  0.0;

        }
        
    } else if (sqrt_s == 8.0) {
        
        C1 = 0.0065;

        mu += 
                +121042. * (1. + eVBF_78_Hbox ) * CHbox / LambdaNP2
                +12739.3 * (1. + eVBF_78_HQ1_11 ) * CHQ1_11 / LambdaNP2
                -28367.7 * (1. + eVBF_78_Hu_11 ) * CHu_11 / LambdaNP2
                +9134.21 * (1. + eVBF_78_Hd_11 ) * CHd_11 / LambdaNP2
                -423704. * (1. + eVBF_78_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -165182. * (1. + eVBF_78_HD ) * CHD / LambdaNP2
                -349.242 * (1. + eVBF_78_HB ) * CHB / LambdaNP2
                -87279.4 * (1. + eVBF_78_HW ) * CHW / LambdaNP2
                -313449. * (1. + eVBF_78_HWB ) * CHWB / LambdaNP2
                -69421.9 * (1. + eVBF_78_HG ) * CHG / LambdaNP2
                -373.338 * (1. + eVBF_78_DHB ) * CDHB / LambdaNP2
                -57028.1 * (1. + eVBF_78_DHW ) * CDHW / LambdaNP2
                -4.472 * (1. + eVBF_78_DeltaGF ) * DeltaGF()
                -3.138 * deltaMw()
                +7.472 * deltaMw2()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

        mu +=  0.0;

        }
    } else if (sqrt_s == 13.0) {
        
        C1 = 0.0064;

        mu += 
                +121798. * (1. + eVBF_1314_Hbox ) * CHbox / LambdaNP2
                +10339.7 * (1. + eVBF_1314_HQ1_11 ) * CHQ1_11 / LambdaNP2
                -30827.2 * (1. + eVBF_1314_Hu_11 ) * CHu_11 / LambdaNP2
                +10564.3 * (1. + eVBF_1314_Hd_11 ) * CHd_11 / LambdaNP2
                -466270. * (1. + eVBF_1314_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -164119. * (1. + eVBF_1314_HD ) * CHD / LambdaNP2
                -61.471 * (1. + eVBF_1314_HB ) * CHB / LambdaNP2
                -82985.3 * (1. + eVBF_1314_HW ) * CHW / LambdaNP2
                -313815. * (1. + eVBF_1314_HWB ) * CHWB / LambdaNP2
                -36554. * (1. + eVBF_1314_HG ) * CHG / LambdaNP2
                -725.694 * (1. + eVBF_1314_DHB ) * CDHB / LambdaNP2
                -65253.4 * (1. + eVBF_1314_DHW ) * CDHW / LambdaNP2
                -4.474 * (1. + eVBF_1314_DeltaGF ) * DeltaGF()
                -3.109 * deltaMw()
                -13.397 * deltaMw2()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 14.0) {
        
        C1 = 0.0064;

        mu += 
                +120948. * (1. + eVBF_1314_Hbox ) * CHbox / LambdaNP2
                +9896.36 * (1. + eVBF_1314_HQ1_11 ) * CHQ1_11 / LambdaNP2
                -31371. * (1. + eVBF_1314_Hu_11 ) * CHu_11 / LambdaNP2
                +10716.4 * (1. + eVBF_1314_Hd_11 ) * CHd_11 / LambdaNP2
                -473497. * (1. + eVBF_1314_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -164672. * (1. + eVBF_1314_HD ) * CHD / LambdaNP2
                -60.253 * (1. + eVBF_1314_HB ) * CHB / LambdaNP2
                -83504.9 * (1. + eVBF_1314_HW ) * CHW / LambdaNP2
                -314059. * (1. + eVBF_1314_HWB ) * CHWB / LambdaNP2
                -33627.6 * (1. + eVBF_1314_HG ) * CHG / LambdaNP2
                -775.959 * (1. + eVBF_1314_DHB ) * CDHB / LambdaNP2
                -66336.3 * (1. + eVBF_1314_DHW ) * CDHW / LambdaNP2
                -4.474 * (1. + eVBF_1314_DeltaGF ) * DeltaGF()
                -3.193 * deltaMw()
                -27.218 * deltaMw2()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }

    } else if (sqrt_s == 100.0) {
        
        C1 = 0.0; // N.A.

        mu += 
                +121714. * CHbox / LambdaNP2
                -2261.73 * CHQ1_11 / LambdaNP2
                -42045.4 * CHu_11 / LambdaNP2
                +17539.2 * CHd_11 / LambdaNP2
                -674206. * CHQ3_11 / LambdaNP2
                -163344. * CHD / LambdaNP2
                +71.488 * CHB / LambdaNP2
                -90808.2 * CHW / LambdaNP2
                -312544. * CHWB / LambdaNP2
                -8165.65 * CHG / LambdaNP2
                -2615.48 * CDHB / LambdaNP2
                -96539.6 * CDHW / LambdaNP2
                -4.452 * DeltaGF()
                -2.949 * deltaMw()
                +14.174 * deltaMw2()
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
    
    if (FlagQuadraticTerms) {
//  Quadratic contribution from Higgs self-coupling
        mu = mu + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();        
    }

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
                +119630. * CHbox / LambdaNP2
                -501300. * CHQ3_11 / LambdaNP2
                -200890. * CHD / LambdaNP2
                +11852.5 * CHB / LambdaNP2
                -131586. * CHW / LambdaNP2
                -361991. * CHWB / LambdaNP2
                -18894.5 * CDHB / LambdaNP2
                -69025.4 * CDHW / LambdaNP2
                +23773.1 * CW / LambdaNP2
                -4.629 * DeltaGF()
                -5.637 * deltaMw()
                +20.854 * deltaMw2()
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
    
    if (FlagQuadraticTerms) {
//  Quadratic contribution from Higgs self-coupling
        mu = mu + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();        
    }

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
                +121120. * CHbox / LambdaNP2
                -138682. * CHL3_11 / LambdaNP2
                -203727. * CHD / LambdaNP2
                -24699.7 * CHW / LambdaNP2
                -379830. * CHWB / LambdaNP2
                -18173.7 * CDHW / LambdaNP2
                -4.716 * DeltaGF()
                -5.665 * deltaMw()
                +20.817 * deltaMw2()
                ;
        
//        if (FlagQuadraticTerms) {
//            //Add contributions that are quadratic in the effective coefficients
//
//        }
          
    } else if (sqrt_s == 0.250) {
        
        C1 = 0.0064;

        mu += 
                +121142. * CHbox / LambdaNP2
                -147357. * CHL3_11 / LambdaNP2
                -203726. * CHD / LambdaNP2
                -26559.2 * CHW / LambdaNP2
                -379797. * CHWB / LambdaNP2
                -19265.3 * CDHW / LambdaNP2
                -4.717 * DeltaGF()
                -5.593 * deltaMw()
                +21.406 * deltaMw2()
                ;
        
//        if (FlagQuadraticTerms) {
//            //Add contributions that are quadratic in the effective coefficients
//
//        }
        
    } else if (sqrt_s == 0.350) {
        
        C1 = 0.0062;

        mu += 
                +121107. * CHbox / LambdaNP2
                -219582. * CHL3_11 / LambdaNP2
                -203717. * CHD / LambdaNP2
                -39722.3 * CHW / LambdaNP2
                -379795. * CHWB / LambdaNP2
                -28864.2 * CDHW / LambdaNP2
                -4.714 * DeltaGF()
                -5.13 * deltaMw()
                +14.296 * deltaMw2()
                ;
        
//        if (FlagQuadraticTerms) {
//            //Add contributions that are quadratic in the effective coefficients
//
//        }
//      
        
    } else if (sqrt_s == 0.365) {
        
        C1 = 0.0062; // Use the same as 350 GeV

        mu += 
                +121071. * CHbox / LambdaNP2
                -228452. * CHL3_11 / LambdaNP2
                -203725. * CHD / LambdaNP2
                -40966.9 * CHW / LambdaNP2
                -379798. * CHWB / LambdaNP2
                -30110.4 * CDHW / LambdaNP2
                -4.714 * DeltaGF()
                -5.08 * deltaMw()
                +13.196 * deltaMw2()
                ;
        
    } else if (sqrt_s == 0.380) {
        
        C1 = 0.0062; // Use the same as 350 GeV

        mu += 
                +121001. * CHbox / LambdaNP2
                -237126. * CHL3_11 / LambdaNP2
                -203726. * CHD / LambdaNP2
                -42070.9 * CHW / LambdaNP2
                -379788. * CHWB / LambdaNP2
                -31352.7 * CDHW / LambdaNP2
                -4.714 * DeltaGF()
                -5.044 * deltaMw()
                +12.295 * deltaMw2()
                ;

    } else if (sqrt_s == 0.500) {
        
        C1 = 0.0061;
        
        mu += 
                +121063. * CHbox / LambdaNP2
                -295115. * CHL3_11 / LambdaNP2
                -203679. * CHD / LambdaNP2
                -47539.5 * CHW / LambdaNP2
                -379773. * CHWB / LambdaNP2
                -39825.1 * CDHW / LambdaNP2
                -4.715 * DeltaGF()
                -4.817 * deltaMw()
                +12.532 * deltaMw2()
                ;
        
    } else if (sqrt_s == 1.0) {
        
        C1 = 0.0059;

        mu += 
                +120960. * CHbox / LambdaNP2
                -442647. * CHL3_11 / LambdaNP2
                -203748. * CHD / LambdaNP2
                -49375.4 * CHW / LambdaNP2
                -379685. * CHWB / LambdaNP2
                -63503.9 * CDHW / LambdaNP2
                -4.712 * DeltaGF()
                -4.481 * deltaMw()
                +5.949 * deltaMw2()
                ;
  
    } else if (sqrt_s == 1.4) {
        
        C1 = 0.0058;

        mu += 
                +121118. * CHbox / LambdaNP2
                -515189. * CHL3_11 / LambdaNP2
                -203684. * CHD / LambdaNP2
                -46619.5 * CHW / LambdaNP2
                -379667. * CHWB / LambdaNP2
                -75747.8 * CDHW / LambdaNP2
                -4.714 * DeltaGF()
                -4.391 * deltaMw()
                +8.567 * deltaMw2()
                ;
        
    } else if (sqrt_s == 1.5) {
        
        C1 = 0.0058;// Use the same as 1400 GeV

        mu += 
                +121200. * CHbox / LambdaNP2
                -530152. * CHL3_11 / LambdaNP2
                -203649. * CHD / LambdaNP2
                -45921.3 * CHW / LambdaNP2
                -379591. * CHWB / LambdaNP2
                -78241.3 * CDHW / LambdaNP2
                -4.715 * DeltaGF()
                -4.38 * deltaMw()
                +12.793 * deltaMw2()
                ;

        
    } else if (sqrt_s == 3.0) {
        
        C1 = 0.0057;
        
        mu += 
                +121321. * CHbox / LambdaNP2
                -684382. * CHL3_11 / LambdaNP2
                -203585. * CHD / LambdaNP2
                -38239. * CHW / LambdaNP2
                -379518. * CHWB / LambdaNP2
                -104465. * CDHW / LambdaNP2
                -4.714 * DeltaGF()
                -4.258 * deltaMw()
                +10.482 * deltaMw2()
                ;
        

    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::mueeWBF()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += eeeWBFint + eeeWBFpar;
    
//  Linear contribution from Higgs self-coupling
    mu = mu + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
    
    if (FlagQuadraticTerms) {
//  Quadratic contribution from Higgs self-coupling
        mu = mu + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();        
    }

    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}


double NPSMEFTd6::mueeWBFPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{

//  Pure WBF, hence only initiated by LH fermions. No difference between polarizations at the linear level.
//  Expand like other functions when quadratic terms are included
    
    return mueeWBF(sqrt_s);
}

double NPSMEFTd6::mueeZBF(const double sqrt_s) const
{
    double mu = 1.0;
    
    double C1 = 0.0;
    
    if (sqrt_s == 0.240) {
        
    C1 = 0.0070;

        mu += 
                +121661. * CHbox / LambdaNP2
                +489617. * CHL1_11 / LambdaNP2
                -357163. * CHe_11 / LambdaNP2
                +489617. * CHL3_11 / LambdaNP2
                -39217.8 * CHD / LambdaNP2
                +1525468. * CHB / LambdaNP2
                +378019. * CHW / LambdaNP2
                +215983. * CHWB / LambdaNP2
                -6554.11 * CDHB / LambdaNP2
                +1175.47 * CDHW / LambdaNP2
                -3.161 * DeltaGF()
                ;
        
    } else if (sqrt_s == 0.250) {
        
    C1 = 0.0070;

        mu += 
                +122144. * CHbox / LambdaNP2
                +444406. * CHL1_11 / LambdaNP2
                -315727. * CHe_11 / LambdaNP2
                +444406. * CHL3_11 / LambdaNP2
                -41440.8 * CHD / LambdaNP2
                +1186855. * CHB / LambdaNP2
                +301913. * CHW / LambdaNP2
                +98540.5 * CHWB / LambdaNP2
                -5766.35 * CDHB / LambdaNP2
                +294.724 * CDHW / LambdaNP2
                -3.279 * DeltaGF()
                ;
        
    } else if (sqrt_s == 0.350) {
        
    C1 = 0.0069;

        mu += 
                +121556. * CHbox / LambdaNP2
                +46354.9 * CHL1_11 / LambdaNP2
                -251.929 * CHe_11 / LambdaNP2
                +46354.9 * CHL3_11 / LambdaNP2
                -43426.2 * CHD / LambdaNP2
                +450512. * CHB / LambdaNP2
                +166493. * CHW / LambdaNP2
                -198898. * CHWB / LambdaNP2
                -4408.76 * CDHB / LambdaNP2
                -17005.2 * CDHW / LambdaNP2
                -3.427 * DeltaGF()
                ;        
        
    } else if (sqrt_s == 0.365) {
        
    C1 = 0.0069; // use same as 350 GeV

        mu += 
                +121067. * CHbox / LambdaNP2
                +9887.64 * CHL1_11 / LambdaNP2
                +27809. * CHe_11 / LambdaNP2
                +9887.64 * CHL3_11 / LambdaNP2
                -43174.2 * CHD / LambdaNP2
                +417865. * CHB / LambdaNP2
                +154270. * CHW / LambdaNP2
                -201517. * CHWB / LambdaNP2
                -4943.82 * CDHB / LambdaNP2
                -19213.5 * CDHW / LambdaNP2
                -3.423 * DeltaGF()
                ;
        
    } else if (sqrt_s == 0.380) {
        
    C1 = 0.0069; // use same as 350 GeV

        mu += 
                +121214. * CHbox / LambdaNP2
                -22289.7 * CHL1_11 / LambdaNP2
                +52903.2 * CHe_11 / LambdaNP2
                -22289.7 * CHL3_11 / LambdaNP2
                -43137.3 * CHD / LambdaNP2
                +388336. * CHB / LambdaNP2
                +140923. * CHW / LambdaNP2
                -202884. * CHWB / LambdaNP2
                -5363.69 * CDHB / LambdaNP2
                -21404.2 * CDHW / LambdaNP2
                -3.418 * DeltaGF()
                ;
        
    } else if (sqrt_s == 0.500) {
        
    C1 = 0.0067;

        mu += 
                +121453. * CHbox / LambdaNP2
                -185326. * CHL1_11 / LambdaNP2
                +178925. * CHe_11 / LambdaNP2
                -185326. * CHL3_11 / LambdaNP2
                -42051.6 * CHD / LambdaNP2
                +236945. * CHB / LambdaNP2
                +67833.5 * CHW / LambdaNP2
                -178623. * CHWB / LambdaNP2
                -8004.61 * CDHB / LambdaNP2
                -33567.3 * CDHW / LambdaNP2
                -3.416 * DeltaGF()
                ;
        
    } else if (sqrt_s == 1.0) {
        
    C1 = 0.0065;

        mu += 
                +121062. * CHbox / LambdaNP2
                -409543. * CHL1_11 / LambdaNP2
                +356730. * CHe_11 / LambdaNP2
                -409543. * CHL3_11 / LambdaNP2
                -42133.9 * CHD / LambdaNP2
                +69851. * CHB / LambdaNP2
                -14416.8 * CHW / LambdaNP2
                -113198. * CHWB / LambdaNP2
                -18688.4 * CDHB / LambdaNP2
                -61696. * CDHW / LambdaNP2
                -3.405 * DeltaGF()
                ;
        
    } else if (sqrt_s == 1.4) {
        
    C1 = 0.0065;

        mu += 
                +120749. * CHbox / LambdaNP2
                -493617. * CHL1_11 / LambdaNP2
                +426669. * CHe_11 / LambdaNP2
                -493617. * CHL3_11 / LambdaNP2
                -42486.9 * CHD / LambdaNP2
                +34633.1 * CHB / LambdaNP2
                -27609.6 * CHW / LambdaNP2
                -97014.2 * CHWB / LambdaNP2
                -23942.2 * CDHB / LambdaNP2
                -74940.3 * CDHW / LambdaNP2
                -3.405 * DeltaGF()
                ;
        
    } else if (sqrt_s == 1.5) {
        
        C1 = 0.0065;// Use the same as 1400 GeV

        mu += 
                +120587. * CHbox / LambdaNP2
                -510290. * CHL1_11 / LambdaNP2
                +440504. * CHe_11 / LambdaNP2
                -510290. * CHL3_11 / LambdaNP2
                -42529.6 * CHD / LambdaNP2
                +30448.1 * CHB / LambdaNP2
                -30741.2 * CHW / LambdaNP2
                -95903.3 * CHWB / LambdaNP2
                -25074.9 * CDHB / LambdaNP2
                -77634.5 * CDHW / LambdaNP2
                -3.401 * DeltaGF()
                ;
        
    } else if (sqrt_s == 3.0) {
        
    C1 = 0.0063;
        
        mu += 
                +120474. * CHbox / LambdaNP2
                -677185. * CHL1_11 / LambdaNP2
                +582037. * CHe_11 / LambdaNP2
                -677185. * CHL3_11 / LambdaNP2
                -42541.3 * CHD / LambdaNP2
                +6810.6 * CHB / LambdaNP2
                -32994.5 * CHW / LambdaNP2
                -78012.3 * CHWB / LambdaNP2
                -36250. * CDHB / LambdaNP2
                -105734. * CDHW / LambdaNP2
                -3.405 * DeltaGF()
                ;

    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZBF()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    //(Assume similar to WBF.)
    mu += eeeWBFint + eeeWBFpar;
    
//  Linear contribution from Higgs self-coupling
    mu = mu + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
    
    if (FlagQuadraticTerms) {
//  Quadratic contribution from Higgs self-coupling
        mu = mu + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();        
    }

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
                +121531. * CHbox / LambdaNP2
                +58943.5 * CHL1_11 / LambdaNP2
                -939512. * CHe_11 / LambdaNP2
                +58943.5 * CHL3_11 / LambdaNP2
                +77442.6 * CHD / LambdaNP2
                +2082256. * CHB / LambdaNP2
                +108043. * CHW / LambdaNP2
                +1362693. * CHWB / LambdaNP2
                +40385. * CDHB / LambdaNP2
                -21886. * CDHW / LambdaNP2
                +0.563 * DeltaGF()
                ;            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +122065. * CHbox / LambdaNP2
                +905327. * CHL1_11 / LambdaNP2
                -55689. * CHe_11 / LambdaNP2
                +905327. * CHL3_11 / LambdaNP2
                -124548. * CHD / LambdaNP2
                +905057. * CHB / LambdaNP2
                +540185. * CHW / LambdaNP2
                -329708. * CHWB / LambdaNP2
                -37296.9 * CDHB / LambdaNP2
                +20497.1 * CDHW / LambdaNP2
                -5.854 * DeltaGF()
                ; 
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121947. * CHbox / LambdaNP2
                +88774.4 * CHL1_11 / LambdaNP2
                -753269. * CHe_11 / LambdaNP2
                +88774.4 * CHL3_11 / LambdaNP2
                +54593.2 * CHD / LambdaNP2
                +2101955. * CHB / LambdaNP2
                +182237. * CHW / LambdaNP2
                +972861. * CHWB / LambdaNP2
                +29346.2 * CDHB / LambdaNP2
                -18562.1 * CDHW / LambdaNP2
                -0.206 * DeltaGF()
                ; 
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +122265. * CHbox / LambdaNP2
                +785643. * CHL1_11 / LambdaNP2
                -66907.6 * CHe_11 / LambdaNP2
                +785643. * CHL3_11 / LambdaNP2
                -107673. * CHD / LambdaNP2
                +1115316. * CHB / LambdaNP2
                +521873. * CHW / LambdaNP2
                -331727. * CHWB / LambdaNP2
                -32442.4 * CDHB / LambdaNP2
                +15348.7 * CDHW / LambdaNP2
                -5.334 * DeltaGF()
                ; 
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZBFPol()");
        }
        
    } else if (sqrt_s == 0.250) {
        
        C1 = 0.0070;
        
        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121054. * CHbox / LambdaNP2
                +51113. * CHL1_11 / LambdaNP2
                -851357. * CHe_11 / LambdaNP2
                +51113. * CHL3_11 / LambdaNP2
                +76762.9 * CHD / LambdaNP2
                +1629614. * CHB / LambdaNP2
                +72741.6 * CHW / LambdaNP2
                +1130834. * CHWB / LambdaNP2
                +34381.7 * CDHB / LambdaNP2
                -19876.5 * CDHW / LambdaNP2
                +0.563 * DeltaGF()
                ;            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121471. * CHbox / LambdaNP2
                +824294. * CHL1_11 / LambdaNP2
                -45066.5 * CHe_11 / LambdaNP2
                +824294. * CHL3_11 / LambdaNP2
                -128864. * CHD / LambdaNP2
                +644513. * CHB / LambdaNP2
                +425051. * CHW / LambdaNP2
                -383720. * CHWB / LambdaNP2
                -32434.3 * CDHB / LambdaNP2
                +15329.4 * CDHW / LambdaNP2
                -6.022 * DeltaGF()
                ; 
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121494. * CHbox / LambdaNP2
                +77372.1 * CHL1_11 / LambdaNP2
                -676199. * CHe_11 / LambdaNP2
                +77372.1 * CHL3_11 / LambdaNP2
                +53294.7 * CHD / LambdaNP2
                +1668830. * CHB / LambdaNP2
                +145010. * CHW / LambdaNP2
                +772902. * CHWB / LambdaNP2
                +23910.6 * CDHB / LambdaNP2
                -16890.6 * CDHW / LambdaNP2
                -0.226 * DeltaGF()
                ; 
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121947. * CHbox / LambdaNP2
                +713174. * CHL1_11 / LambdaNP2
                -53393.3 * CHe_11 / LambdaNP2
                +713174. * CHL3_11 / LambdaNP2
                -111120. * CHD / LambdaNP2
                +843388. * CHB / LambdaNP2
                +417838. * CHW / LambdaNP2
                -386753. * CHWB / LambdaNP2
                -27915.7 * CDHB / LambdaNP2
                +11946.5 * CDHW / LambdaNP2
                -5.496 * DeltaGF()
                ; 
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZBFPol()");
        }
        
    } else if (sqrt_s == 0.350) {
        
        C1 = 0.0069;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121674. * CHbox / LambdaNP2
                -47420.2 * CHL1_11 / LambdaNP2
                -172088. * CHe_11 / LambdaNP2
                -47420.2 * CHL3_11 / LambdaNP2
                +59728. * CHD / LambdaNP2
                +544205. * CHB / LambdaNP2
                +83604.4 * CHW / LambdaNP2
                +435393. * CHWB / LambdaNP2
                -24800.4 * CDHB / LambdaNP2
                -4583.09 * CDHW / LambdaNP2
                -0.05 * DeltaGF()
                ;            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121541. * CHbox / LambdaNP2
                +197618. * CHL1_11 / LambdaNP2
                +42238.9 * CHe_11 / LambdaNP2
                +197618. * CHL3_11 / LambdaNP2
                -124376. * CHD / LambdaNP2
                +181154. * CHB / LambdaNP2
                +195329. * CHW / LambdaNP2
                -505800. * CHWB / LambdaNP2
                +13082.6 * CDHB / LambdaNP2
                -26607.4 * CDHW / LambdaNP2
                -6.096 * DeltaGF()
                ; 
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121760. * CHbox / LambdaNP2
                -62853. * CHL1_11 / LambdaNP2
                -83019.6 * CHe_11 / LambdaNP2
                -62853. * CHL3_11 / LambdaNP2
                +34395.4 * CHD / LambdaNP2
                +623389. * CHB / LambdaNP2
                +123932. * CHW / LambdaNP2
                +181789. * CHWB / LambdaNP2
                -20420. * CDHB / LambdaNP2
                -7820.42 * CDHW / LambdaNP2
                -0.875 * DeltaGF()
                ; 
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121557. * CHbox / LambdaNP2
                +131443. * CHL1_11 / LambdaNP2
                +63326.7 * CHe_11 / LambdaNP2
                +131443. * CHL3_11 / LambdaNP2
                -103038. * CHD / LambdaNP2
                +323596. * CHB / LambdaNP2
                +201676. * CHW / LambdaNP2
                -491019. * CHWB / LambdaNP2
                +7992.43 * CDHB / LambdaNP2
                -24283.6 * CDHW / LambdaNP2
                -5.391 * DeltaGF()
                ; 
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZBFPol()");
        } 
    
    } else if (sqrt_s == 0.365) {
        
        C1 = 0.0069; // Use same as 350 GeV

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121458. * CHbox / LambdaNP2
                -58695.1 * CHL1_11 / LambdaNP2
                -109686. * CHe_11 / LambdaNP2
                -58695.1 * CHL3_11 / LambdaNP2
                +58496.7 * CHD / LambdaNP2
                +489137. * CHB / LambdaNP2
                +80751.3 * CHW / LambdaNP2
                +410304. * CHWB / LambdaNP2
                -30918.3 * CDHB / LambdaNP2
                -3571.31 * CDHW / LambdaNP2
                -0.085 * DeltaGF()
                ;            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121152. * CHbox / LambdaNP2
                +136019. * CHL1_11 / LambdaNP2
                +50762. * CHe_11 / LambdaNP2
                +136019. * CHL3_11 / LambdaNP2
                -123859. * CHD / LambdaNP2
                +165799. * CHB / LambdaNP2
                +176652. * CHW / LambdaNP2
                -504889. * CHWB / LambdaNP2
                +16920.7 * CDHB / LambdaNP2
                -31414.1 * CDHW / LambdaNP2
                -6.076 * DeltaGF()
                ; 
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121193. * CHbox / LambdaNP2
                -76905.7 * CHL1_11 / LambdaNP2
                -32264.3 * CHe_11 / LambdaNP2
                -76905.7 * CHL3_11 / LambdaNP2
                +33650.3 * CHD / LambdaNP2
                +573505. * CHB / LambdaNP2
                +117937. * CHW / LambdaNP2
                +166382. * CHWB / LambdaNP2
                -25012.1 * CDHB / LambdaNP2
                -7703.47 * CDHW / LambdaNP2
                -0.911 * DeltaGF()
                ; 
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121177. * CHbox / LambdaNP2
                +77981.5 * CHL1_11 / LambdaNP2
                +74274.1 * CHe_11 / LambdaNP2
                +77981.5 * CHL3_11 / LambdaNP2
                -102068. * CHD / LambdaNP2
                +305730. * CHB / LambdaNP2
                +183682. * CHW / LambdaNP2
                -487770. * CHWB / LambdaNP2
                +10624.8 * CDHB / LambdaNP2
                -28092.3 * CDHW / LambdaNP2
                -5.366 * DeltaGF()
                ; 
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZBFPol()");
        } 
    
    } else if (sqrt_s == 0.380) {
        
        C1 = 0.0069; // Use same as 350 GeV

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121392. * CHbox / LambdaNP2
                -68799.8 * CHL1_11 / LambdaNP2
                -54383.2 * CHe_11 / LambdaNP2
                -68799.8 * CHL3_11 / LambdaNP2
                +57427.7 * CHD / LambdaNP2
                +439155. * CHB / LambdaNP2
                +76978.2 * CHW / LambdaNP2
                +392293. * CHWB / LambdaNP2
                -36175.9 * CDHB / LambdaNP2
                -3193.74 * CDHW / LambdaNP2
                -0.11 * DeltaGF()
                ;            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121306. * CHbox / LambdaNP2
                +80159.7 * CHL1_11 / LambdaNP2
                +58002.2 * CHe_11 / LambdaNP2
                +80159.7 * CHL3_11 / LambdaNP2
                -123524. * CHD / LambdaNP2
                +151617. * CHB / LambdaNP2
                +154342. * CHW / LambdaNP2
                -500961. * CHWB / LambdaNP2
                +20509.9 * CDHB / LambdaNP2
                -35718.1 * CDHW / LambdaNP2
                -6.064 * DeltaGF()
                ; 
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121171. * CHbox / LambdaNP2
                -89494.3 * CHL1_11 / LambdaNP2
                +11882.3 * CHe_11 / LambdaNP2
                -89494.3 * CHL3_11 / LambdaNP2
                +32430.1 * CHD / LambdaNP2
                +524620. * CHB / LambdaNP2
                +111520. * CHW / LambdaNP2
                +156122. * CHWB / LambdaNP2
                -29271.1 * CDHB / LambdaNP2
                -8056.8 * CDHW / LambdaNP2
                -0.928 * DeltaGF()
                ; 
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121286. * CHbox / LambdaNP2
                +30046.7 * CHL1_11 / LambdaNP2
                +84014. * CHe_11 / LambdaNP2
                +30046.7 * CHL3_11 / LambdaNP2
                -101539. * CHD / LambdaNP2
                +286981. * CHB / LambdaNP2
                +164662. * CHW / LambdaNP2
                -480410. * CHWB / LambdaNP2
                +13149.6 * CDHB / LambdaNP2
                -31886.7 * CDHW / LambdaNP2
                -5.346 * DeltaGF()
                ; 
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZBFPol()");
        } 
    
    } else if (sqrt_s == 0.500) {
        
        C1 = 0.0067;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121372. * CHbox / LambdaNP2
                -121062. * CHL1_11 / LambdaNP2
                +224754. * CHe_11 / LambdaNP2
                -121062. * CHL3_11 / LambdaNP2
                +55161.7 * CHD / LambdaNP2
                +201238. * CHB / LambdaNP2
                +52456.6 * CHW / LambdaNP2
                +335517. * CHWB / LambdaNP2
                -63733.4 * CDHB / LambdaNP2
                -2379.21 * CDHW / LambdaNP2
                -0.207 * DeltaGF()
                ;            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121399. * CHbox / LambdaNP2
                -200849. * CHL1_11 / LambdaNP2
                +96427.7 * CHe_11 / LambdaNP2
                -200849. * CHL3_11 / LambdaNP2
                -121178. * CHD / LambdaNP2
                +83220.9 * CHB / LambdaNP2
                +42832.2 * CHW / LambdaNP2
                -464173. * CHWB / LambdaNP2
                +37654.2 * CDHB / LambdaNP2
                -59029.6 * CDHW / LambdaNP2
                -6.025 * DeltaGF()
                ; 
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121435. * CHbox / LambdaNP2
                -154953. * CHL1_11 / LambdaNP2
                +235326. * CHe_11 / LambdaNP2
                -154953. * CHL3_11 / LambdaNP2
                +30472. * CHD / LambdaNP2
                +298145. * CHB / LambdaNP2
                +75047.6 * CHW / LambdaNP2
                +137304. * CHWB / LambdaNP2
                -49636.1 * CDHB / LambdaNP2
                -10277.1 * CDHW / LambdaNP2
                -1.027 * DeltaGF()
                ; 
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121468. * CHbox / LambdaNP2
                -208577. * CHL1_11 / LambdaNP2
                +134790. * CHe_11 / LambdaNP2
                -208577. * CHL3_11 / LambdaNP2
                -98708.1 * CHD / LambdaNP2
                +190310. * CHB / LambdaNP2
                +62321.4 * CHW / LambdaNP2
                -429412. * CHWB / LambdaNP2
                +24628.2 * CDHB / LambdaNP2
                -51722.9 * CDHW / LambdaNP2
                -5.287 * DeltaGF()
                ; 
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZBFPol()");
        } 
    
    } else if (sqrt_s == 1.0) {
        
        C1 = 0.0065;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121044. * CHbox / LambdaNP2
                -206156. * CHL1_11 / LambdaNP2
                +586357. * CHe_11 / LambdaNP2
                -206156. * CHL3_11 / LambdaNP2
                +54157.3 * CHD / LambdaNP2
                -30839.6 * CHB / LambdaNP2
                +18110.3 * CHW / LambdaNP2
                +345253. * CHWB / LambdaNP2
                -108488. * CDHB / LambdaNP2
                -12324.2 * CDHW / LambdaNP2
                -0.229 * DeltaGF()
                ;            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121085. * CHbox / LambdaNP2
                -565700. * CHL1_11 / LambdaNP2
                +157498. * CHe_11 / LambdaNP2
                -565700. * CHL3_11 / LambdaNP2
                -120795. * CHD / LambdaNP2
                +7953.6 * CHB / LambdaNP2
                -79908.9 * CHW / LambdaNP2
                -402278. * CHWB / LambdaNP2
                +54805.3 * CDHB / LambdaNP2
                -101988. * CDHW / LambdaNP2
                -6.001 * DeltaGF()
                ; 
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +120996. * CHbox / LambdaNP2
                -263143. * CHL1_11 / LambdaNP2
                +533190. * CHe_11 / LambdaNP2
                -263143. * CHL3_11 / LambdaNP2
                +29434.5 * CHD / LambdaNP2
                +63176.5 * CHB / LambdaNP2
                +26728.5 * CHW / LambdaNP2
                +184228. * CHWB / LambdaNP2
                -85487.1 * CDHB / LambdaNP2
                -24906.1 * CDHW / LambdaNP2
                -1.044 * DeltaGF()
                ; 
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121114. * CHbox / LambdaNP2
                -524119. * CHL1_11 / LambdaNP2
                +218758. * CHe_11 / LambdaNP2
                -524119. * CHL3_11 / LambdaNP2
                -98164. * CHD / LambdaNP2
                +74694.7 * CHB / LambdaNP2
                -49060.4 * CHW / LambdaNP2
                -348619. * CHWB / LambdaNP2
                +33861.6 * CDHB / LambdaNP2
                -90369.8 * CDHW / LambdaNP2
                -5.256 * DeltaGF()
                ; 
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZBFPol()");
        } 
    
    } else if (sqrt_s == 1.4) {
        
        C1 = 0.0065;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +120762. * CHbox / LambdaNP2
                -242720. * CHL1_11 / LambdaNP2
                +714345. * CHe_11 / LambdaNP2
                -242720. * CHL3_11 / LambdaNP2
                +53823.3 * CHD / LambdaNP2
                -64876.7 * CHB / LambdaNP2
                +9362.37 * CHW / LambdaNP2
                +355440. * CHWB / LambdaNP2
                -127361. * CDHB / LambdaNP2
                -18147.3 * CDHW / LambdaNP2
                -0.228 * DeltaGF()
                ;            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +120818. * CHbox / LambdaNP2
                -692905. * CHL1_11 / LambdaNP2
                +184416. * CHe_11 / LambdaNP2
                -692905. * CHL3_11 / LambdaNP2
                -121143. * CHD / LambdaNP2
                -4989.81 * CHB / LambdaNP2
                -93241.6 * CHW / LambdaNP2
                -392394. * CHWB / LambdaNP2
                +60556.9 * CDHB / LambdaNP2
                -121409. * CDHW / LambdaNP2
                -6.003 * DeltaGF()
                ; 
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +120773. * CHbox / LambdaNP2
                -309806. * CHL1_11 / LambdaNP2
                +643900. * CHe_11 / LambdaNP2
                -309806. * CHL3_11 / LambdaNP2
                +29091.1 * CHD / LambdaNP2
                +22438.3 * CHB / LambdaNP2
                +16021.7 * CHW / LambdaNP2
                +202496. * CHWB / LambdaNP2
                -100775. * CDHB / LambdaNP2
                -32830.8 * CDHW / LambdaNP2
                -1.043 * DeltaGF()
                ; 
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +120795. * CHbox / LambdaNP2
                -637584. * CHL1_11 / LambdaNP2
                +256188. * CHe_11 / LambdaNP2
                -637584. * CHL3_11 / LambdaNP2
                -98543.3 * CHD / LambdaNP2
                +49040.2 * CHB / LambdaNP2
                -63051.7 * CHW / LambdaNP2
                -332850. * CHWB / LambdaNP2
                +36510.1 * CDHB / LambdaNP2
                -108018. * CDHW / LambdaNP2
                -5.256 * DeltaGF()
                ; 
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZBFPol()");
        } 
        
    } else if (sqrt_s == 1.5) {
        
        C1 = 0.0065;// Use the same as 1400 GeV

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +120570. * CHbox / LambdaNP2
                -250340. * CHL1_11 / LambdaNP2
                +739684. * CHe_11 / LambdaNP2
                -250340. * CHL3_11 / LambdaNP2
                +53685.8 * CHD / LambdaNP2
                -71192.9 * CHB / LambdaNP2
                +9743.41 * CHW / LambdaNP2
                +357556. * CHWB / LambdaNP2
                -131206. * CDHB / LambdaNP2
                -19448. * CDHW / LambdaNP2
                -0.224 * DeltaGF()
                ;            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +120602. * CHbox / LambdaNP2
                -718001. * CHL1_11 / LambdaNP2
                +189852. * CHe_11 / LambdaNP2
                -718001. * CHL3_11 / LambdaNP2
                -121214. * CHD / LambdaNP2
                -6057.91 * CHB / LambdaNP2
                -95148.1 * CHW / LambdaNP2
                -390958. * CHWB / LambdaNP2
                +61690.7 * CDHB / LambdaNP2
                -125382. * CDHW / LambdaNP2
                -5.997 * DeltaGF()
                ; 
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +120563. * CHbox / LambdaNP2
                -319378. * CHL1_11 / LambdaNP2
                +665765. * CHe_11 / LambdaNP2
                -319378. * CHL3_11 / LambdaNP2
                +29010.7 * CHD / LambdaNP2
                +14190.4 * CHB / LambdaNP2
                +16080. * CHW / LambdaNP2
                +205187. * CHWB / LambdaNP2
                -103927. * CDHB / LambdaNP2
                -34420.2 * CDHW / LambdaNP2
                -1.04 * DeltaGF()
                ; 
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +120607. * CHbox / LambdaNP2
                -659879. * CHL1_11 / LambdaNP2
                +263841. * CHe_11 / LambdaNP2
                -659879. * CHL3_11 / LambdaNP2
                -98617.3 * CHD / LambdaNP2
                +46418.4 * CHB / LambdaNP2
                -64166.6 * CHW / LambdaNP2
                -330855. * CHWB / LambdaNP2
                +36774.5 * CDHB / LambdaNP2
                -111573. * CDHW / LambdaNP2
                -5.253 * DeltaGF()
                ; 
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZBFPol()");
        }
    
    } else if (sqrt_s == 3.0) {
        
        C1 = 0.0063;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +120539. * CHbox / LambdaNP2
                -327096. * CHL1_11 / LambdaNP2
                +988310. * CHe_11 / LambdaNP2
                -327096. * CHL3_11 / LambdaNP2
                +53758.1 * CHD / LambdaNP2
                -79161. * CHB / LambdaNP2
                +3856.87 * CHW / LambdaNP2
                +369878. * CHWB / LambdaNP2
                -170059. * CDHB / LambdaNP2
                -32235.8 * CDHW / LambdaNP2
                -0.226 * DeltaGF()
                ;            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +120565. * CHbox / LambdaNP2
                -961658. * CHL1_11 / LambdaNP2
                +247947. * CHe_11 / LambdaNP2
                -961658. * CHL3_11 / LambdaNP2
                -121230. * CHD / LambdaNP2
                -10752.9 * CHB / LambdaNP2
                -92123.7 * CHW / LambdaNP2
                -391807. * CHWB / LambdaNP2
                +73242.2 * CDHB / LambdaNP2
                -165690. * CDHW / LambdaNP2
                -6.002 * DeltaGF()
                ; 
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +120534. * CHbox / LambdaNP2
                -417962. * CHL1_11 / LambdaNP2
                +884851. * CHe_11 / LambdaNP2
                -417962. * CHL3_11 / LambdaNP2
                +29065.5 * CHD / LambdaNP2
                -10885.4 * CHB / LambdaNP2
                +8249.25 * CHW / LambdaNP2
                +228820. * CHWB / LambdaNP2
                -135851. * CDHB / LambdaNP2
                -51177.2 * CDHW / LambdaNP2
                -1.04 * DeltaGF()
                ; 
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +120480. * CHbox / LambdaNP2
                -880604. * CHL1_11 / LambdaNP2
                +344657. * CHe_11 / LambdaNP2
                -880604. * CHL3_11 / LambdaNP2
                -98656.8 * CHD / LambdaNP2
                +28681.4 * CHB / LambdaNP2
                -66216.6 * CHW / LambdaNP2
                -320715. * CHWB / LambdaNP2
                +41721.6 * CDHB / LambdaNP2
                -148698. * CDHW / LambdaNP2
                -5.256 * DeltaGF()
                ; 
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
    
    if (FlagQuadraticTerms) {
//  Quadratic contribution from Higgs self-coupling
        mu = mu + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();        
    }

    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::muepWBF(const double sqrt_s) const
{
    double mu = 1.0;
    
    if (sqrt_s == 1.3) {
        
        mu += 
                +121790. * CHbox / LambdaNP2
                -161604. * CHL3_11 / LambdaNP2
                -161282. * CHQ3_11 / LambdaNP2
                -203141. * CHD / LambdaNP2
                -88171.6 * CHW / LambdaNP2
                -377218. * CHWB / LambdaNP2
                -37738.9 * CDHW / LambdaNP2
                -4.676 * DeltaGF()
                -4.916 * deltaMw()
                +3.627 * deltaMw2()
                ;
        
//        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

//        }
        
    } else if (sqrt_s == 1.8) {
        
        mu += 
                +121867. * CHbox / LambdaNP2
                -182643. * CHL3_11 / LambdaNP2
                -181961. * CHQ3_11 / LambdaNP2
                -202400. * CHD / LambdaNP2
                -78295.8 * CHW / LambdaNP2
                -377193. * CHWB / LambdaNP2
                -45757.3 * CDHW / LambdaNP2
                -4.672 * DeltaGF()
                -4.637 * deltaMw()
                +38.639 * deltaMw2()
                ;
        
//        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

//        }
        
    } else if (sqrt_s == 3.5) {

        mu += 
                +121250. * CHbox / LambdaNP2
                -216885. * CHL3_11 / LambdaNP2
                -218544. * CHQ3_11 / LambdaNP2
                -202390. * CHD / LambdaNP2
                -64783.2 * CHW / LambdaNP2
                -377727. * CHWB / LambdaNP2
                -60431.2 * CDHW / LambdaNP2
                -4.688 * DeltaGF()
                -4.573 * deltaMw()
                +14.015 * deltaMw2()
                ;
        
//        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

//        }
          
    } else if (sqrt_s == 5.0) {

        mu += 
                +119662. * CHbox / LambdaNP2
                -237868. * CHL3_11 / LambdaNP2
                -236470. * CHQ3_11 / LambdaNP2
                -203294. * CHD / LambdaNP2
                -60911. * CHW / LambdaNP2
                -378045. * CHWB / LambdaNP2
                -67483.7 * CDHW / LambdaNP2
                -4.667 * DeltaGF()
                -4.437 * deltaMw()
                -26.106 * deltaMw2()
                ;
        
//        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

//        }
        
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::muepWBF()");

    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::muepZBF(const double sqrt_s) const
{
    double mu = 1.0;

    if (sqrt_s == 1.3) {
        
        mu += 
                +121280. * CHbox / LambdaNP2
                -152367. * CHL1_11 / LambdaNP2
                +32200. * CHQ1_11 / LambdaNP2
                +124934. * CHe_11 / LambdaNP2
                -42209.5 * CHu_11 / LambdaNP2
                +12445.7 * CHd_11 / LambdaNP2
                -152367. * CHL3_11 / LambdaNP2
                -165343. * CHQ3_11 / LambdaNP2
                -173922. * CHD / LambdaNP2
                -34636.2 * CHB / LambdaNP2
                -121438. * CHW / LambdaNP2
                -74939.1 * CHWB / LambdaNP2
                -5454.93 * CDHB / LambdaNP2
                -39349.6 * CDHW / LambdaNP2
                -3.719 * DeltaGF()
                ;
        
//        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

//        }
        
    } else if (sqrt_s == 1.8) {
        
        mu += 
                +120218. * CHbox / LambdaNP2
                -173566. * CHL1_11 / LambdaNP2
                +26307.1 * CHQ1_11 / LambdaNP2
                +142600. * CHe_11 / LambdaNP2
                -47449. * CHu_11 / LambdaNP2
                +14356.2 * CHd_11 / LambdaNP2
                -173566. * CHL3_11 / LambdaNP2
                -188606. * CHQ3_11 / LambdaNP2
                -174301. * CHD / LambdaNP2
                -19800. * CHB / LambdaNP2
                -103254. * CHW / LambdaNP2
                -89049.2 * CHWB / LambdaNP2
                -8304.85 * CDHB / LambdaNP2
                -48942.9 * CDHW / LambdaNP2
                -3.714 * DeltaGF()
                ;
        
//        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

//        }
        
    } else if (sqrt_s == 3.5) {

        mu +=              
                +123119. * CHbox / LambdaNP2
                -206981. * CHL1_11 / LambdaNP2
                +18620.9 * CHQ1_11 / LambdaNP2
                +177706. * CHe_11 / LambdaNP2
                -53822. * CHu_11 / LambdaNP2
                +20491.5 * CHd_11 / LambdaNP2
                -206981. * CHL3_11 / LambdaNP2
                -227549. * CHQ3_11 / LambdaNP2
                -172298. * CHD / LambdaNP2
                -6887.17 * CHB / LambdaNP2
                -79245. * CHW / LambdaNP2
                -103223. * CHWB / LambdaNP2
                -9863.11 * CDHB / LambdaNP2
                -61304.3 * CDHW / LambdaNP2
                -3.721 * DeltaGF()
                ;
        
//        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

//        }
          
    } else if (sqrt_s == 5.0) {

        mu += 
                +121709. * CHbox / LambdaNP2
                -225267. * CHL1_11 / LambdaNP2
                +13471.8 * CHQ1_11 / LambdaNP2
                +193542. * CHe_11 / LambdaNP2
                -57640.9 * CHu_11 / LambdaNP2
                +22573. * CHd_11 / LambdaNP2
                -225267. * CHL3_11 / LambdaNP2
                -247738. * CHQ3_11 / LambdaNP2
                -172768. * CHD / LambdaNP2
                -4524.89 * CHB / LambdaNP2
                -71935.4 * CHW / LambdaNP2
                -104998. * CHWB / LambdaNP2
                -11877.8 * CDHB / LambdaNP2
                -69467.3 * CDHW / LambdaNP2
                -3.71 * DeltaGF()
                ;
        
//        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

//        }
        
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::muepZBF()");

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
                +121173. * (1. + eWH_2_Hbox ) * CHbox / LambdaNP2
                +1566788. * (1. + eWH_2_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -160914. * (1. + eWH_2_HD ) * CHD / LambdaNP2
                +860916. * (1. + eWH_2_HW ) * CHW / LambdaNP2
                -286409. * (1. + eWH_2_HWB ) * CHWB / LambdaNP2
                +134641. * (1. + eWH_2_DHW ) * CDHW / LambdaNP2
                -3.31 * (1. + eWH_2_DeltaGF ) * DeltaGF()
                -2.199 * deltaMw()
                +0.804 * deltaMw2()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
        
    } else if (sqrt_s == 7.0) {
        
        C1 = 0.0106;

        mu += 
                +121015. * (1. + eWH_78_Hbox ) * CHbox / LambdaNP2
                +1792020. * (1. + eWH_78_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -159689. * (1. + eWH_78_HD ) * CHD / LambdaNP2
                +881065. * (1. + eWH_78_HW ) * CHW / LambdaNP2
                -283895. * (1. + eWH_78_HWB ) * CHWB / LambdaNP2
                +168173. * (1. + eWH_78_DHW ) * CDHW / LambdaNP2
                -3.273 * (1. + eWH_78_DeltaGF ) * DeltaGF()
                -2.143 * deltaMw()
                -2.312 * deltaMw2()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
        
    } else if (sqrt_s == 8.0) {
        
        C1 = 0.0105;

        mu += 
                +121226. * (1. + eWH_78_Hbox ) * CHbox / LambdaNP2
                +1830192. * (1. + eWH_78_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -159543. * (1. + eWH_78_HD ) * CHD / LambdaNP2
                +884671. * (1. + eWH_78_HW ) * CHW / LambdaNP2
                -283662. * (1. + eWH_78_HWB ) * CHWB / LambdaNP2
                +174061. * (1. + eWH_78_DHW ) * CDHW / LambdaNP2
                -3.278 * (1. + eWH_78_DeltaGF ) * DeltaGF()
                -2.147 * deltaMw()
                +6.049 * deltaMw2()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }

    } else if (sqrt_s == 13.0) {
        
        C1 = 0.0103;

        mu += 
                +120439. * (1. + eWH_1314_Hbox ) * CHbox / LambdaNP2
                +1953200. * (1. + eWH_1314_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -159847. * (1. + eWH_1314_HD ) * CHD / LambdaNP2
                +892264. * (1. + eWH_1314_HW ) * CHW / LambdaNP2
                -283830. * (1. + eWH_1314_HWB ) * CHWB / LambdaNP2
                +192168. * (1. + eWH_1314_DHW ) * CDHW / LambdaNP2
                -3.269 * (1. + eWH_1314_DeltaGF ) * DeltaGF()
                -2.101 * deltaMw()
                -15.282 * deltaMw2()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
      
    } else if (sqrt_s == 14.0) {
        
        C1 = 0.0103;

        mu += 
                +120284. * (1. + eWH_1314_Hbox ) * CHbox / LambdaNP2
                +1971011. * (1. + eWH_1314_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -159830. * (1. + eWH_1314_HD ) * CHD / LambdaNP2
                +893216. * (1. + eWH_1314_HW ) * CHW / LambdaNP2
                -283818. * (1. + eWH_1314_HWB ) * CHWB / LambdaNP2
                +194877. * (1. + eWH_1314_DHW ) * CDHW / LambdaNP2
                -3.272 * (1. + eWH_1314_DeltaGF ) * DeltaGF()
                -2.103 * deltaMw()
                -15.576 * deltaMw2()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
          
    } else if (sqrt_s == 100.0) {
        
        C1 = 0.0; // N.A. 

        mu += 
                +121319. * CHbox / LambdaNP2
                +2294991. * CHQ3_11 / LambdaNP2
                -159242. * CHD / LambdaNP2
                +908130. * CHW / LambdaNP2
                -282574. * CHWB / LambdaNP2
                +245406. * CDHW / LambdaNP2
                -3.259 * DeltaGF()
                -2.047 * deltaMw()
                +0. * deltaMw2()
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
    
    if (FlagQuadraticTerms) {
//  Quadratic contribution from Higgs self-coupling
        mu = mu + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();        
    }

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
                +121197. * (1. + eZH_2_Hbox ) * CHbox / LambdaNP2
                -810445. * (1. + eZH_2_HQ1_11 ) * CHQ1_11 / LambdaNP2
                +529340. * (1. + eZH_2_Hu_11 ) * CHu_11 / LambdaNP2
                -69410.3 * (1. + eZH_2_Hd_11 ) * CHd_11 / LambdaNP2
                +1567161. * (1. + eZH_2_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -16992.5 * (1. + eZH_2_HD ) * CHD / LambdaNP2
                +79314.5 * (1. + eZH_2_HB ) * CHB / LambdaNP2
                +711710. * (1. + eZH_2_HW ) * CHW / LambdaNP2
                +189054. * (1. + eZH_2_HWB ) * CHWB / LambdaNP2
                +9774.73 * (1. + eZH_2_DHB ) * CDHB / LambdaNP2
                +130777. * (1. + eZH_2_DHW ) * CDHW / LambdaNP2
                -2.535 * (1. + eZH_2_DeltaGF ) * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
        
    } else if (sqrt_s == 7.0) {
        
        C1 = 0.0123;

        mu += 
                +121069. * (1. + eZH_78_Hbox ) * CHbox / LambdaNP2
                -182215. * (1. + eZH_78_HQ1_11 ) * CHQ1_11 / LambdaNP2
                +421780. * (1. + eZH_78_Hu_11 ) * CHu_11 / LambdaNP2
                -139169. * (1. + eZH_78_Hd_11 ) * CHd_11 / LambdaNP2
                +1712111. * (1. + eZH_78_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -15395.4 * (1. + eZH_78_HD ) * CHD / LambdaNP2
                +87094.9 * (1. + eZH_78_HB ) * CHB / LambdaNP2
                +717388. * (1. + eZH_78_HW ) * CHW / LambdaNP2
                +203105. * (1. + eZH_78_HWB ) * CHWB / LambdaNP2
                +17532.4 * (1. + eZH_78_DHB ) * CDHB / LambdaNP2
                +152950. * (1. + eZH_78_DHW ) * CDHW / LambdaNP2
                -2.502 * (1. + eZH_78_DeltaGF ) * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
        
    } else if (sqrt_s == 8.0) {
        
        C1 = 0.0122;

        mu += 
                +121334. * (1. + eZH_78_Hbox ) * CHbox / LambdaNP2
                -176804. * (1. + eZH_78_HQ1_11 ) * CHQ1_11 / LambdaNP2
                +428587. * (1. + eZH_78_Hu_11 ) * CHu_11 / LambdaNP2
                -142508. * (1. + eZH_78_Hd_11 ) * CHd_11 / LambdaNP2
                +1747367. * (1. + eZH_78_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -15002.7 * (1. + eZH_78_HD ) * CHD / LambdaNP2
                +87781.5 * (1. + eZH_78_HB ) * CHB / LambdaNP2
                +721405. * (1. + eZH_78_HW ) * CHW / LambdaNP2
                +204540. * (1. + eZH_78_HWB ) * CHWB / LambdaNP2
                +18868.6 * (1. + eZH_78_DHB ) * CDHB / LambdaNP2
                +158432. * (1. + eZH_78_DHW ) * CDHW / LambdaNP2
                -2.507 * (1. + eZH_78_DeltaGF ) * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }

    } else if (sqrt_s == 13.0) {
        
        C1 = 0.0119;

        mu += 
                +121374. * (1. + eZH_1314_Hbox ) * CHbox / LambdaNP2
                -152273. * (1. + eZH_1314_HQ1_11 ) * CHQ1_11 / LambdaNP2
                +448168. * (1. + eZH_1314_Hu_11 ) * CHu_11 / LambdaNP2
                -155999. * (1. + eZH_1314_Hd_11 ) * CHd_11 / LambdaNP2
                +1862364. * (1. + eZH_1314_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -15185. * (1. + eZH_1314_HD ) * CHD / LambdaNP2
                +88937.9 * (1. + eZH_1314_HB ) * CHB / LambdaNP2
                +728207. * (1. + eZH_1314_HW ) * CHW / LambdaNP2
                +207857. * (1. + eZH_1314_HWB ) * CHWB / LambdaNP2
                +21647.4 * (1. + eZH_1314_DHB ) * CDHB / LambdaNP2
                +175015. * (1. + eZH_1314_DHW ) * CDHW / LambdaNP2
                -2.506 * (1. + eZH_1314_DeltaGF ) * DeltaGF()
                ;

        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
        
    } else if (sqrt_s == 14.0) {
        
        C1 = 0.0118;

        mu += 
                +121437. * (1. + eZH_1314_Hbox ) * CHbox / LambdaNP2
                -147580. * (1. + eZH_1314_HQ1_11 ) * CHQ1_11 / LambdaNP2
                +450628. * (1. + eZH_1314_Hu_11 ) * CHu_11 / LambdaNP2
                -157625. * (1. + eZH_1314_Hd_11 ) * CHd_11 / LambdaNP2
                +1878132. * (1. + eZH_1314_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -15299.4 * (1. + eZH_1314_HD ) * CHD / LambdaNP2
                +88761.8 * (1. + eZH_1314_HB ) * CHB / LambdaNP2
                +729243. * (1. + eZH_1314_HW ) * CHW / LambdaNP2
                +207707. * (1. + eZH_1314_HWB ) * CHWB / LambdaNP2
                +21958.9 * (1. + eZH_1314_DHB ) * CDHB / LambdaNP2
                +177300. * (1. + eZH_1314_DHW ) * CDHW / LambdaNP2
                -2.507 * (1. + eZH_1314_DeltaGF ) * DeltaGF()
                ;

        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;

        }
          
    } else if (sqrt_s == 100.0) {
        
        C1 = 0.0; // N.A.
 
        mu += 
                +121269. * CHbox / LambdaNP2
                +90.68 * CHQ1_11 / LambdaNP2
                +484275. * CHu_11 / LambdaNP2
                -197878. * CHd_11 / LambdaNP2
                +2175601. * CHQ3_11 / LambdaNP2
                -14992.4 * CHD / LambdaNP2
                +91707.3 * CHB / LambdaNP2
                +741805. * CHW / LambdaNP2
                +215319. * CHWB / LambdaNP2
                +31435.6 * CDHB / LambdaNP2
                +223843. * CDHW / LambdaNP2
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
    
    if (FlagQuadraticTerms) {
//  Quadratic contribution from Higgs self-coupling
        mu = mu + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();        
    }

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
                +121263. * CHbox / LambdaNP2
                +898682. * CHL1_11 / LambdaNP2
                -767820. * CHe_11 / LambdaNP2
                +898682. * CHL3_11 / LambdaNP2
                -6046.36 * CHD / LambdaNP2
                +122439. * CHB / LambdaNP2
                +540057. * CHW / LambdaNP2
                +231063. * CHWB / LambdaNP2
                +17593.2 * CDHB / LambdaNP2
                +53409.5 * CDHW / LambdaNP2
                -2.2 * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
            }

    } else if (sqrt_s == 0.250) {
        
        C1 = 0.015;

        mu += 
                +121263. * CHbox / LambdaNP2
                +975101. * CHL1_11 / LambdaNP2
                -833750. * CHe_11 / LambdaNP2
                +975101. * CHL3_11 / LambdaNP2
                -6046.36 * CHD / LambdaNP2
                +128443. * CHB / LambdaNP2
                +568273. * CHW / LambdaNP2
                +244206. * CHWB / LambdaNP2
                +19818.6 * CDHB / LambdaNP2
                +60127.6 * CDHW / LambdaNP2
                -2.2 * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 0.350) {
        
        C1 = 0.0057;

        mu += 
                +121283. * CHbox / LambdaNP2
                +1911340. * CHL1_11 / LambdaNP2
                -1640958. * CHe_11 / LambdaNP2
                +1911340. * CHL3_11 / LambdaNP2
                -6009.52 * CHD / LambdaNP2
                +173183. * CHB / LambdaNP2
                +785843. * CHW / LambdaNP2
                +344494. * CHWB / LambdaNP2
                +59158.7 * CDHB / LambdaNP2
                +167954. * CDHW / LambdaNP2
                -2.201 * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 0.365) {
        
        C1 = 0.0057; // Use same as 350 GeV

        mu += 
                +121243. * CHbox / LambdaNP2
                +2078482. * CHL1_11 / LambdaNP2
                -1785085. * CHe_11 / LambdaNP2
                +2078482. * CHL3_11 / LambdaNP2
                -6010.65 * CHD / LambdaNP2
                +178173. * CHB / LambdaNP2
                +809806. * CHW / LambdaNP2
                +355487. * CHWB / LambdaNP2
                +67662.7 * CDHB / LambdaNP2
                +190194. * CDHW / LambdaNP2
                -2.201 * DeltaGF()
                ;
        
    } else if (sqrt_s == 0.380) {
        
        C1 = 0.0057; // Use same as 350 GeV

        mu += 
                +121281. * CHbox / LambdaNP2
                +2253013. * CHL1_11 / LambdaNP2
                -1934557. * CHe_11 / LambdaNP2
                +2253013. * CHL3_11 / LambdaNP2
                -6026.37 * CHD / LambdaNP2
                +182674. * CHB / LambdaNP2
                +832109. * CHW / LambdaNP2
                +365819. * CHWB / LambdaNP2
                +76742. * CDHB / LambdaNP2
                +214030. * CDHW / LambdaNP2
                -2.202 * DeltaGF()
                ;        
        
    } else if (sqrt_s == 0.500) {
        
        C1 = 0.00099;

        mu += 
                +121264. * CHbox / LambdaNP2
                +3900384. * CHL1_11 / LambdaNP2
                -3350136. * CHe_11 / LambdaNP2
                +3900384. * CHL3_11 / LambdaNP2
                -6019.22 * CHD / LambdaNP2
                +209229. * CHB / LambdaNP2
                +959942. * CHW / LambdaNP2
                +425112. * CHWB / LambdaNP2
                +169841. * CDHB / LambdaNP2
                +455437. * CDHW / LambdaNP2
                -2.202 * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 1.0) {
        
        C1 = -0.0012;

        mu += 
                +121274. * CHbox / LambdaNP2
                +15601820. * CHL1_11 / LambdaNP2
                -13395670. * CHe_11 / LambdaNP2
                +15601820. * CHL3_11 / LambdaNP2
                -6040.16 * CHD / LambdaNP2
                +243960. * CHB / LambdaNP2
                +1128805. * CHW / LambdaNP2
                +503138. * CHWB / LambdaNP2
                +899357. * CDHB / LambdaNP2
                +2321619. * CDHW / LambdaNP2
                -2.202 * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 1.4) {
        
        C1 = -0.0011;

        mu += 
                +121283. * CHbox / LambdaNP2
                +30579278. * CHL1_11 / LambdaNP2
                -26253064. * CHe_11 / LambdaNP2
                +30579278. * CHL3_11 / LambdaNP2
                -6010.77 * CHD / LambdaNP2
                +250804. * CHB / LambdaNP2
                +1161208. * CHW / LambdaNP2
                +518040. * CHWB / LambdaNP2
                +1848758. * CDHB / LambdaNP2
                +4747422. * CDHW / LambdaNP2
                -2.203 * DeltaGF()
                ;
        
    } else if (sqrt_s == 1.5) {
        
        C1 = -0.0011;// Use the same as 1400 GeV

        mu += 
                +121262. * CHbox / LambdaNP2
                +35102329. * CHL1_11 / LambdaNP2
                -30135878. * CHe_11 / LambdaNP2
                +35102329. * CHL3_11 / LambdaNP2
                -6034.22 * CHD / LambdaNP2
                +251576. * CHB / LambdaNP2
                +1165634. * CHW / LambdaNP2
                +519954. * CHWB / LambdaNP2
                +2132554. * CDHB / LambdaNP2
                +5481906. * CDHW / LambdaNP2
                -2.203 * DeltaGF()
                ;
        
    } else if (sqrt_s == 3.0) {
        
        C1 = -0.00054;
        
        mu += 
                +121279. * CHbox / LambdaNP2
                +140413697. * CHL1_11 / LambdaNP2
                -120540988. * CHe_11 / LambdaNP2
                +140413697. * CHL3_11 / LambdaNP2
                -6012.61 * CHD / LambdaNP2
                +257222. * CHB / LambdaNP2
                +1188444. * CHW / LambdaNP2
                +530503. * CHWB / LambdaNP2
                +8839419. * CDHB / LambdaNP2
                +22583370. * CDHW / LambdaNP2
                -2.202 * DeltaGF()
                ;
        
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZH()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += eeeZHint + eeeZHpar;
    
//  Linear contribution from Higgs self-coupling
    mu = mu + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
    
    if (FlagQuadraticTerms) {
//  Quadratic contribution from Higgs self-coupling
        mu = mu + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();        
    }
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::mueeZHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    double mu = 1.0;
    
    double C1 = 0.0;

    if (sqrt_s == 0.240) {
        
        C1 = 0.017;
        
        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121260. * CHbox / LambdaNP2
                +117191. * CHL1_11 / LambdaNP2
                -1681596. * CHe_11 / LambdaNP2
                +117191. * CHL3_11 / LambdaNP2
                +74555.1 * CHD / LambdaNP2
                +528105. * CHB / LambdaNP2
                +134403. * CHW / LambdaNP2
                +872560. * CHWB / LambdaNP2
                +137571. * CDHB / LambdaNP2
                -12321.5 * CDHW / LambdaNP2
                +0.459 * DeltaGF()
                ;            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121254. * CHbox / LambdaNP2
                +1495015. * CHL1_11 / LambdaNP2
                -76567.2 * CHe_11 / LambdaNP2
                +1495015. * CHL3_11 / LambdaNP2
                -67582.1 * CHD / LambdaNP2
                -187104. * CHB / LambdaNP2
                +849552. * CHW / LambdaNP2
                -258537. * CHWB / LambdaNP2
                -73970.1 * CDHB / LambdaNP2
                +103582. * CDHW / LambdaNP2
                -4.23 * DeltaGF()
                ; 
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121256. * CHbox / LambdaNP2
                +204529. * CHL1_11 / LambdaNP2
                -1578998. * CHe_11 / LambdaNP2
                +204529. * CHL3_11 / LambdaNP2
                +65548.7 * CHD / LambdaNP2
                +482729. * CHB / LambdaNP2
                +179733. * CHW / LambdaNP2
                +800870. * CHWB / LambdaNP2
                +124170. * CDHB / LambdaNP2
                -5016.48 * CDHW / LambdaNP2
                +0.162 * DeltaGF()
                ; 
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121264. * CHbox / LambdaNP2
                +1442776. * CHL1_11 / LambdaNP2
                -137405. * CHe_11 / LambdaNP2
                +1442776. * CHL3_11 / LambdaNP2
                -62167.6 * CHD / LambdaNP2
                -159988. * CHB / LambdaNP2
                +822448. * CHW / LambdaNP2
                -215639. * CHWB / LambdaNP2
                -65950.1 * CDHB / LambdaNP2
                +99206.1 * CDHW / LambdaNP2
                -4.052 * DeltaGF()
                ; 
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZHPol()");
        }
        
    } else if (sqrt_s == 0.250) {
        
        C1 = 0.015;
        
        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121264. * CHbox / LambdaNP2
                +127210. * CHL1_11 / LambdaNP2
                -1824910. * CHe_11 / LambdaNP2
                +127210. * CHL3_11 / LambdaNP2
                +74597.1 * CHD / LambdaNP2
                +560319. * CHB / LambdaNP2
                +136129. * CHW / LambdaNP2
                +902676. * CHWB / LambdaNP2
                +154358. * CDHB / LambdaNP2
                -13612.9 * CDHW / LambdaNP2
                +0.459 * DeltaGF()
                ;            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121257. * CHbox / LambdaNP2
                +1622228. * CHL1_11 / LambdaNP2
                -83107. * CHe_11 / LambdaNP2
                +1622228. * CHL3_11 / LambdaNP2
                -67554.3 * CHD / LambdaNP2
                -201409. * CHB / LambdaNP2
                +898116. * CHW / LambdaNP2
                -258306. * CHWB / LambdaNP2
                -82898. * CDHB / LambdaNP2
                +116421. * CDHW / LambdaNP2
                -4.23 * DeltaGF()
                ; 
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121309. * CHbox / LambdaNP2
                +221930. * CHL1_11 / LambdaNP2
                -1714047. * CHe_11 / LambdaNP2
                +221930. * CHL3_11 / LambdaNP2
                +65599.6 * CHD / LambdaNP2
                +512136. * CHB / LambdaNP2
                +184424. * CHW / LambdaNP2
                +829145. * CHWB / LambdaNP2
                +139369. * CDHB / LambdaNP2
                -5351.17 * CDHW / LambdaNP2
                +0.162 * DeltaGF()
                ; 
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121269. * CHbox / LambdaNP2
                +1565559. * CHL1_11 / LambdaNP2
                -148908. * CHe_11 / LambdaNP2
                +1565559. * CHL3_11 / LambdaNP2
                -62170. * CHD / LambdaNP2
                -172540. * CHB / LambdaNP2
                +869218. * CHW / LambdaNP2
                -214299. * CHWB / LambdaNP2
                -73929.8 * CDHB / LambdaNP2
                +111494. * CDHW / LambdaNP2
                -4.053 * DeltaGF()
                ; 
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZHPol()");
        }
        
    } else if (sqrt_s == 0.350) {
        
        C1 = 0.0057;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121274. * CHbox / LambdaNP2
                +249309. * CHL1_11 / LambdaNP2
                -3576996. * CHe_11 / LambdaNP2
                +249309. * CHL3_11 / LambdaNP2
                +74596.5 * CHD / LambdaNP2
                +812491. * CHB / LambdaNP2
                +146212. * CHW / LambdaNP2
                +1135161. * CHWB / LambdaNP2
                +395085. * CDHB / LambdaNP2
                -16140.8 * CDHW / LambdaNP2
                +0.458 * DeltaGF()
                ;            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121289. * CHbox / LambdaNP2
                +3179548. * CHL1_11 / LambdaNP2
                -163347. * CHe_11 / LambdaNP2
                +3179548. * CHL3_11 / LambdaNP2
                -67524.8 * CHD / LambdaNP2
                -314653. * CHB / LambdaNP2
                +1273817. * CHW / LambdaNP2
                -258947. * CHWB / LambdaNP2
                -197137. * CDHB / LambdaNP2
                +308384. * CDHW / LambdaNP2
                -4.231 * DeltaGF()
                ; 
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121304. * CHbox / LambdaNP2
                +434952. * CHL1_11 / LambdaNP2
                -3360980. * CHe_11 / LambdaNP2
                +434952. * CHL3_11 / LambdaNP2
                +65624.7 * CHD / LambdaNP2
                +741142. * CHB / LambdaNP2
                +217654. * CHW / LambdaNP2
                +1046799. * CHWB / LambdaNP2
                +357606. * CDHB / LambdaNP2
                +4440.1 * CDHW / LambdaNP2
                +0.161 * DeltaGF()
                ; 
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121259. * CHbox / LambdaNP2
                +3068356. * CHL1_11 / LambdaNP2
                -292427. * CHe_11 / LambdaNP2
                +3068356. * CHL3_11 / LambdaNP2
                -62160.7 * CHD / LambdaNP2
                -271962. * CHB / LambdaNP2
                +1231171. * CHW / LambdaNP2
                -206112. * CHWB / LambdaNP2
                -174718. * CDHB / LambdaNP2
                +296046. * CDHW / LambdaNP2
                -4.053 * DeltaGF()
                ; 
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZHPol()");
        } 
    
    } else if (sqrt_s == 0.365) {
        
        C1 = 0.0057; // Use same as 350 GeV

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121270. * CHbox / LambdaNP2
                +271098. * CHL1_11 / LambdaNP2
                -3890169. * CHe_11 / LambdaNP2
                +271098. * CHL3_11 / LambdaNP2
                +74554. * CHD / LambdaNP2
                +840573. * CHB / LambdaNP2
                +147108. * CHW / LambdaNP2
                +1160947. * CHWB / LambdaNP2
                +442125. * CDHB / LambdaNP2
                -15038.8 * CDHW / LambdaNP2
                +0.459 * DeltaGF()
                ;            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121238. * CHbox / LambdaNP2
                +3457848. * CHL1_11 / LambdaNP2
                -177584. * CHe_11 / LambdaNP2
                +3457848. * CHL3_11 / LambdaNP2
                -67578.3 * CHD / LambdaNP2
                -327391. * CHB / LambdaNP2
                +1315671. * CHW / LambdaNP2
                -259142. * CHWB / LambdaNP2
                -218241. * CDHB / LambdaNP2
                +346804. * CDHW / LambdaNP2
                -4.231 * DeltaGF()
                ; 
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121251. * CHbox / LambdaNP2
                +472985. * CHL1_11 / LambdaNP2
                -3655203. * CHe_11 / LambdaNP2
                +472985. * CHL3_11 / LambdaNP2
                +65559.4 * CHD / LambdaNP2
                +766585. * CHB / LambdaNP2
                +221202. * CHW / LambdaNP2
                +1070933. * CHWB / LambdaNP2
                +400293. * CDHB / LambdaNP2
                +7914.02 * CDHW / LambdaNP2
                +0.161 * DeltaGF()
                ; 
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121238. * CHbox / LambdaNP2
                +3336984. * CHL1_11 / LambdaNP2
                -317944. * CHe_11 / LambdaNP2
                +3336984. * CHL3_11 / LambdaNP2
                -62188.9 * CHD / LambdaNP2
                -283174. * CHB / LambdaNP2
                +1271272. * CHW / LambdaNP2
                -205330. * CHWB / LambdaNP2
                -193153. * CDHB / LambdaNP2
                +333078. * CDHW / LambdaNP2
                -4.053 * DeltaGF()
                ; 
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZHPol()");
        } 
    
    } else if (sqrt_s == 0.380) {
        
        C1 = 0.0057; // Use same as 350 GeV

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121228. * CHbox / LambdaNP2
                +293860. * CHL1_11 / LambdaNP2
                -4216491. * CHe_11 / LambdaNP2
                +293860. * CHL3_11 / LambdaNP2
                +74561.4 * CHD / LambdaNP2
                +866754. * CHB / LambdaNP2
                +147982. * CHW / LambdaNP2
                +1184912. * CHWB / LambdaNP2
                +492018. * CDHB / LambdaNP2
                -13596.5 * CDHW / LambdaNP2
                +0.459 * DeltaGF()
                ;            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121226. * CHbox / LambdaNP2
                +3747707. * CHL1_11 / LambdaNP2
                -192650. * CHe_11 / LambdaNP2
                +3747707. * CHL3_11 / LambdaNP2
                -67608.3 * CHD / LambdaNP2
                -339193. * CHB / LambdaNP2
                +1354040. * CHW / LambdaNP2
                -259321. * CHWB / LambdaNP2
                -240311. * CDHB / LambdaNP2
                +387710. * CDHW / LambdaNP2
                -4.23 * DeltaGF()
                ; 
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121325. * CHbox / LambdaNP2
                +512707. * CHL1_11 / LambdaNP2
                -3961665. * CHe_11 / LambdaNP2
                +512707. * CHL3_11 / LambdaNP2
                +65601.7 * CHD / LambdaNP2
                +790306. * CHB / LambdaNP2
                +224394. * CHW / LambdaNP2
                +1093297. * CHWB / LambdaNP2
                +445530. * CDHB / LambdaNP2
                +11860.4 * CDHW / LambdaNP2
                +0.161 * DeltaGF()
                ; 
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121273. * CHbox / LambdaNP2
                +3617032. * CHL1_11 / LambdaNP2
                -344629. * CHe_11 / LambdaNP2
                +3617032. * CHL3_11 / LambdaNP2
                -62148.3 * CHD / LambdaNP2
                -293491. * CHB / LambdaNP2
                +1308558. * CHW / LambdaNP2
                -204594. * CHWB / LambdaNP2
                -212514. * CDHB / LambdaNP2
                +372554. * CDHW / LambdaNP2
                -4.053 * DeltaGF()
                ; 
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZHPol()");
        } 
    
    } else if (sqrt_s == 0.500) {
        
        C1 = 0.00099;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121268. * CHbox / LambdaNP2
                +508715. * CHL1_11 / LambdaNP2
                -7299333. * CHe_11 / LambdaNP2
                +508715. * CHL3_11 / LambdaNP2
                +74603.6 * CHD / LambdaNP2
                +1018069. * CHB / LambdaNP2
                +151257. * CHW / LambdaNP2
                +1323862. * CHWB / LambdaNP2
                +985604. * CDHB / LambdaNP2
                +8362.16 * CDHW / LambdaNP2
                +0.458 * DeltaGF()
                ;            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121273. * CHbox / LambdaNP2
                +6488707. * CHL1_11 / LambdaNP2
                -332950. * CHe_11 / LambdaNP2
                +6488707. * CHL3_11 / LambdaNP2
                -67530.9 * CHD / LambdaNP2
                -408101. * CHB / LambdaNP2
                +1576859. * CHW / LambdaNP2
                -260777. * CHWB / LambdaNP2
                -452746. * CDHB / LambdaNP2
                +796569. * CDHW / LambdaNP2
                -4.231 * DeltaGF()
                ; 
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121280. * CHbox / LambdaNP2
                +887632. * CHL1_11 / LambdaNP2
                -6858533. * CHe_11 / LambdaNP2
                +887632. * CHL3_11 / LambdaNP2
                +65606.6 * CHD / LambdaNP2
                +927745. * CHB / LambdaNP2
                +241619. * CHW / LambdaNP2
                +1223535. * CHWB / LambdaNP2
                +894441. * CDHB / LambdaNP2
                +58317. * CDHW / LambdaNP2
                +0.161 * DeltaGF()
                ; 
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121268. * CHbox / LambdaNP2
                +6262095. * CHL1_11 / LambdaNP2
                -597046. * CHe_11 / LambdaNP2
                +6262095. * CHL3_11 / LambdaNP2
                -62148.8 * CHD / LambdaNP2
                -353914. * CHB / LambdaNP2
                +1522841. * CHW / LambdaNP2
                -200684. * CHWB / LambdaNP2
                -398214. * CDHB / LambdaNP2
                +766821. * CDHW / LambdaNP2
                -4.054 * DeltaGF()
                ; 
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZHPol()");
        } 
    
    } else if (sqrt_s == 1.0) {
        
        C1 = -0.0012;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121236. * CHbox / LambdaNP2
                +2034785. * CHL1_11 / LambdaNP2
                -29195703. * CHe_11 / LambdaNP2
                +2034785. * CHL3_11 / LambdaNP2
                +74612.7 * CHD / LambdaNP2
                +1218284. * CHB / LambdaNP2
                +154779. * CHW / LambdaNP2
                +1507673. * CHWB / LambdaNP2
                +4701988. * CDHB / LambdaNP2
                +239404. * CDHW / LambdaNP2
                +0.458 * DeltaGF()
                ;            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121298. * CHbox / LambdaNP2
                +25954994. * CHL1_11 / LambdaNP2
                -1333713. * CHe_11 / LambdaNP2
                +25954994. * CHL3_11 / LambdaNP2
                -67536.7 * CHD / LambdaNP2
                -499699. * CHB / LambdaNP2
                +1872177. * CHW / LambdaNP2
                -263454. * CHWB / LambdaNP2
                -1999387. * CDHB / LambdaNP2
                +3910434. * CDHW / LambdaNP2
                -4.233 * DeltaGF()
                ; 
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121307. * CHbox / LambdaNP2
                +3550656. * CHL1_11 / LambdaNP2
                -27432206. * CHe_11 / LambdaNP2
                +3550656. * CHL3_11 / LambdaNP2
                +65607.4 * CHD / LambdaNP2
                +1109435. * CHB / LambdaNP2
                +263679. * CHW / LambdaNP2
                +1395519. * CHWB / LambdaNP2
                +4277336. * CDHB / LambdaNP2
                +472106. * CDHW / LambdaNP2
                +0.159 * DeltaGF()
                ; 
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121327. * CHbox / LambdaNP2
                +25048839. * CHL1_11 / LambdaNP2
                -2390358. * CHe_11 / LambdaNP2
                +25048839. * CHL3_11 / LambdaNP2
                -62132.7 * CHD / LambdaNP2
                -434824. * CHB / LambdaNP2
                +1807095. * CHW / LambdaNP2
                -196264. * CHWB / LambdaNP2
                -1746222. * CDHB / LambdaNP2
                +3771341. * CDHW / LambdaNP2
                -4.056 * DeltaGF()
                ; 
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZHPol()");
        } 
    
    } else if (sqrt_s == 1.4) {
        
        C1 = -0.0011;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121277. * CHbox / LambdaNP2
                +3988231. * CHL1_11 / LambdaNP2
                -57226150. * CHe_11 / LambdaNP2
                +3988231. * CHL3_11 / LambdaNP2
                +74608.5 * CHD / LambdaNP2
                +1256970. * CHB / LambdaNP2
                +155358. * CHW / LambdaNP2
                +1542655. * CHWB / LambdaNP2
                +9506894. * CDHB / LambdaNP2
                +553431. * CDHW / LambdaNP2
                +0.457 * DeltaGF()
                ;            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121314. * CHbox / LambdaNP2
                +50871646. * CHL1_11 / LambdaNP2
                -2614134. * CHe_11 / LambdaNP2
                +50871646. * CHL3_11 / LambdaNP2
                -67535.5 * CHD / LambdaNP2
                -516385. * CHB / LambdaNP2
                +1928805. * CHW / LambdaNP2
                -264072. * CHWB / LambdaNP2
                -3989947. * CDHB / LambdaNP2
                +7948308. * CDHW / LambdaNP2
                -4.233 * DeltaGF()
                ; 
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121250. * CHbox / LambdaNP2
                +6958750. * CHL1_11 / LambdaNP2
                -53762500. * CHe_11 / LambdaNP2
                +6958750. * CHL3_11 / LambdaNP2
                +65589.3 * CHD / LambdaNP2
                +1144464. * CHB / LambdaNP2
                +267732. * CHW / LambdaNP2
                +1428214. * CHWB / LambdaNP2
                +8650536. * CDHB / LambdaNP2
                +1021964. * CDHW / LambdaNP2
                +0.16 * DeltaGF()
                ; 
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121278. * CHbox / LambdaNP2
                +49094486. * CHL1_11 / LambdaNP2
                -4685522. * CHe_11 / LambdaNP2
                +49094486. * CHL3_11 / LambdaNP2
                -62150.9 * CHD / LambdaNP2
                -450090. * CHB / LambdaNP2
                +1861602. * CHW / LambdaNP2
                -195621. * CHWB / LambdaNP2
                -3478338. * CDHB / LambdaNP2
                +7668095. * CDHW / LambdaNP2
                -4.055 * DeltaGF()
                ; 
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZHPol()");
        } 
        
    } else if (sqrt_s == 1.5) {
        
        C1 = -0.0011;// Use the same as 1400 GeV

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121268. * CHbox / LambdaNP2
                +4578315. * CHL1_11 / LambdaNP2
                -65691823. * CHe_11 / LambdaNP2
                +4578315. * CHL3_11 / LambdaNP2
                +74595.2 * CHD / LambdaNP2
                +1262261. * CHB / LambdaNP2
                +155435. * CHW / LambdaNP2
                +1547379. * CHWB / LambdaNP2
                +10961322. * CDHB / LambdaNP2
                +649157. * CDHW / LambdaNP2
                +0.457 * DeltaGF()
                ;            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121277. * CHbox / LambdaNP2
                +58398883. * CHL1_11 / LambdaNP2
                -3000385. * CHe_11 / LambdaNP2
                +58398883. * CHL3_11 / LambdaNP2
                -67535.8 * CHD / LambdaNP2
                -518798. * CHB / LambdaNP2
                +1936613. * CHW / LambdaNP2
                -264171. * CHWB / LambdaNP2
                -4590136. * CDHB / LambdaNP2
                +9169803. * CDHW / LambdaNP2
                -4.233 * DeltaGF()
                ; 
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121289. * CHbox / LambdaNP2
                +7988570. * CHL1_11 / LambdaNP2
                -61718691. * CHe_11 / LambdaNP2
                +7988570. * CHL3_11 / LambdaNP2
                +65599. * CHD / LambdaNP2
                +1149083. * CHB / LambdaNP2
                +268317. * CHW / LambdaNP2
                +1432777. * CHWB / LambdaNP2
                +9972576. * CDHB / LambdaNP2
                +1188554. * CDHW / LambdaNP2
                +0.16 * DeltaGF()
                ; 
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121259. * CHbox / LambdaNP2
                +56356946. * CHL1_11 / LambdaNP2
                -5378233. * CHe_11 / LambdaNP2
                +56356946. * CHL3_11 / LambdaNP2
                -62168.7 * CHD / LambdaNP2
                -452149. * CHB / LambdaNP2
                +1869136. * CHW / LambdaNP2
                -195562. * CHWB / LambdaNP2
                -4000306. * CDHB / LambdaNP2
                +8846432. * CDHW / LambdaNP2
                -4.055 * DeltaGF()
                ; 
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZHPol()");
        }
    
    } else if (sqrt_s == 3.0) {
        
        C1 = -0.00054;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121320. * CHbox / LambdaNP2
                +18314161. * CHL1_11 / LambdaNP2
                -262773345. * CHe_11 / LambdaNP2
                +18314161. * CHL3_11 / LambdaNP2
                +74663.6 * CHD / LambdaNP2
                +1289569. * CHB / LambdaNP2
                +155612. * CHW / LambdaNP2
                +1572580. * CHWB / LambdaNP2
                +44806408. * CDHB / LambdaNP2
                +2877519. * CDHW / LambdaNP2
                +0.456 * DeltaGF()
                ;            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121305. * CHbox / LambdaNP2
                +233598342. * CHL1_11 / LambdaNP2
                -12002450. * CHe_11 / LambdaNP2
                +233598342. * CHL3_11 / LambdaNP2
                -67507.7 * CHD / LambdaNP2
                -531387. * CHB / LambdaNP2
                +1976750. * CHW / LambdaNP2
                -264661. * CHWB / LambdaNP2
                -18587969. * CDHB / LambdaNP2
                +37618569. * CDHW / LambdaNP2
                -4.233 * DeltaGF()
                ; 
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121225. * CHbox / LambdaNP2
                +31953446. * CHL1_11 / LambdaNP2
                -246870182. * CHe_11 / LambdaNP2
                +31953446. * CHL3_11 / LambdaNP2
                +65576.5 * CHD / LambdaNP2
                +1173703. * CHB / LambdaNP2
                +270983. * CHW / LambdaNP2
                +1456032. * CHWB / LambdaNP2
                +40783748. * CDHB / LambdaNP2
                +5077924. * CDHW / LambdaNP2
                +0.16 * DeltaGF()
                ; 
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121248. * CHbox / LambdaNP2
                +225427310. * CHL1_11 / LambdaNP2
                -21505526. * CHe_11 / LambdaNP2
                +225427310. * CHL3_11 / LambdaNP2
                -62193.4 * CHD / LambdaNP2
                -463403. * CHB / LambdaNP2
                +1907593. * CHW / LambdaNP2
                -195017. * CHWB / LambdaNP2
                -16188019. * CDHB / LambdaNP2
                +36299719. * CDHW / LambdaNP2
                -4.054 * DeltaGF()
                ; 
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZHPol()");
        } 
    
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZHPol()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += eeeZHint + eeeZHpar;
    
//  Linear contribution from Higgs self-coupling
    mu = mu + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
    
    if (FlagQuadraticTerms) {
//  Quadratic contribution from Higgs self-coupling
        mu = mu + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();        
    }

    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
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
                +566792. * (1. + ettH_2_uG_33r ) * CuG_33r / LambdaNP2
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
                +811678. * (1. + ettH_78_uG_33r ) * CuG_33r / LambdaNP2
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
                +825379. * (1. + ettH_78_uG_33r ) * CuG_33r / LambdaNP2
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
                +860470. * (1. + ettH_1314_uG_33r ) * CuG_33r / LambdaNP2
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
                +863670. * (1. + ettH_1314_uG_33r ) * CuG_33r / LambdaNP2
                -2.839 * (1. + ettH_1314_DeltagHt ) * deltaG_hff(quarks[TOP]).real()
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
                +880378. * CuG_33r / LambdaNP2
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
    
    if (FlagQuadraticTerms) {
//  Quadratic contribution from Higgs self-coupling
        mu = mu + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();        
    }

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
                +121901. * CHbox / LambdaNP2
                +84038.2 * CHL1_11 / LambdaNP2
                +41671.2 * CHe_11 / LambdaNP2
                -31418.2 * CHu_11 / LambdaNP2
                +84038.2 * CHL3_11 / LambdaNP2
                -121791. * CuH_33r / LambdaNP2
                -59467.6 * CHD / LambdaNP2
                +138929. * CHB / LambdaNP2
                +130909. * CHW / LambdaNP2
                -253030. * CHWB / LambdaNP2
                -1757.66 * CDHB / LambdaNP2
                +1501.34 * CDHW / LambdaNP2
                +1386027. * CuW_33r / LambdaNP2
                +1698012. * CuB_33r / LambdaNP2
                -1.965 * DeltaGF()
                -1.187 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 1.0) {
        
        C1 = 0.017;

        mu += 
                +122013. * CHbox / LambdaNP2
                +889282. * CHL1_11 / LambdaNP2
                -543424. * CHe_11 / LambdaNP2
                -8240.83 * CHu_11 / LambdaNP2
                +889282. * CHL3_11 / LambdaNP2
                -116099. * CuH_33r / LambdaNP2
                -60351.9 * CHD / LambdaNP2
                +352804. * CHB / LambdaNP2
                +361918. * CHW / LambdaNP2
                -397547. * CHWB / LambdaNP2
                +37326.1 * CDHB / LambdaNP2
                +113772. * CDHW / LambdaNP2
                +2758980. * CuW_33r / LambdaNP2
                +3462941. * CuB_33r / LambdaNP2
                -2.08 * DeltaGF()
                -2.575 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  0.0;
        }
        
    } else if (sqrt_s == 1.4) {
        
        C1 = 0.0094;
       
        mu += 
                +122081. * CHbox / LambdaNP2
                +2544832. * CHL1_11 / LambdaNP2
                -1901938. * CHe_11 / LambdaNP2
                +3241.73 * CHu_11 / LambdaNP2
                +2544832. * CHL3_11 / LambdaNP2
                -112208. * CuH_33r / LambdaNP2
                -60340.4 * CHD / LambdaNP2
                +464967. * CHB / LambdaNP2
                +487659. * CHW / LambdaNP2
                -471053. * CHWB / LambdaNP2
                +134900. * CDHB / LambdaNP2
                +371767. * CDHW / LambdaNP2
                +3804096. * CuW_33r / LambdaNP2
                +4800265. * CuB_33r / LambdaNP2
                -2.139 * DeltaGF()
                -3.203 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ;
        
    } else if (sqrt_s == 1.5) {
        
        C1 = 0.0094;// Use the same as 1400 GeV

        mu += 
                +122173. * CHbox / LambdaNP2
                +3117293. * CHL1_11 / LambdaNP2
                -2378233. * CHe_11 / LambdaNP2
                +5531.15 * CHu_11 / LambdaNP2
                +3117293. * CHL3_11 / LambdaNP2
                -111274. * CuH_33r / LambdaNP2
                -60192. * CHD / LambdaNP2
                +487962. * CHB / LambdaNP2
                +513503. * CHW / LambdaNP2
                -485782. * CHWB / LambdaNP2
                +170734. * CDHB / LambdaNP2
                +462665. * CDHW / LambdaNP2
                +4068326. * CuW_33r / LambdaNP2
                +5138930. * CuB_33r / LambdaNP2
                -2.149 * DeltaGF()
                -3.325 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ;
        
    } else if (sqrt_s == 3.0) {
        
        C1 = 0.0037;
        
        mu += 
                +121915. * CHbox / LambdaNP2
                +19529668. * CHL1_11 / LambdaNP2
                -16356276. * CHe_11 / LambdaNP2
                +23142.9 * CHu_11 / LambdaNP2
                +19529668. * CHL3_11 / LambdaNP2
                -104011. * CuH_33r / LambdaNP2
                -58710.4 * CHD / LambdaNP2
                +697868. * CHB / LambdaNP2
                +751003. * CHW / LambdaNP2
                -625171. * CHWB / LambdaNP2
                +1204441. * CDHB / LambdaNP2
                +3111413. * CDHW / LambdaNP2
                +8604912. * CuW_33r / LambdaNP2
                +10946841. * CuB_33r / LambdaNP2
                -2.224 * DeltaGF()
                -4.279 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ;
        
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::mueettH()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += eeettHint + eeettHpar;
    
//  Linear contribution from Higgs self-coupling
    mu = mu + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
    
    if (FlagQuadraticTerms) {
//  Quadratic contribution from Higgs self-coupling
        mu = mu + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();        
    }

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
                +121861. * CHbox / LambdaNP2
                +14207.9 * CHL1_11 / LambdaNP2
                +124191. * CHe_11 / LambdaNP2
                +112591. * CHu_11 / LambdaNP2
                +14207.9 * CHL3_11 / LambdaNP2
                -123399. * CuH_33r / LambdaNP2
                -12437.7 * CHD / LambdaNP2
                +249779. * CHB / LambdaNP2
                +18912.8 * CHW / LambdaNP2
                -109936. * CHWB / LambdaNP2
                -5170.73 * CDHB / LambdaNP2
                +3167.65 * CDHW / LambdaNP2
                +174267. * CuW_33r / LambdaNP2
                +3032981. * CuB_33r / LambdaNP2
                -0.388 * DeltaGF()
                +3.51 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ;            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121809. * CHbox / LambdaNP2
                +116253. * CHL1_11 / LambdaNP2
                +3415.4 * CHe_11 / LambdaNP2
                -98311.8 * CHu_11 / LambdaNP2
                +116253. * CHL3_11 / LambdaNP2
                -121117. * CuH_33r / LambdaNP2
                -81321.2 * CHD / LambdaNP2
                +87352.2 * CHB / LambdaNP2
                +182702. * CHW / LambdaNP2
                -319427. * CHWB / LambdaNP2
                -21.616 * CDHB / LambdaNP2
                +799.81 * CDHW / LambdaNP2
                +1948272. * CuW_33r / LambdaNP2
                +1078489. * CuB_33r / LambdaNP2
                -2.697 * DeltaGF()
                -3.37 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ;
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121837. * CHbox / LambdaNP2
                +24323.6 * CHL1_11 / LambdaNP2
                +111998. * CHe_11 / LambdaNP2
                +91391.1 * CHu_11 / LambdaNP2
                +24323.6 * CHL3_11 / LambdaNP2
                -123203. * CuH_33r / LambdaNP2
                -19404.2 * CHD / LambdaNP2
                +233452. * CHB / LambdaNP2
                +35310.2 * CHW / LambdaNP2
                -131019. * CHWB / LambdaNP2
                -4810.06 * CDHB / LambdaNP2
                +2842.31 * CDHW / LambdaNP2
                +351790. * CuW_33r / LambdaNP2
                +2837005. * CuB_33r / LambdaNP2
                -0.617 * DeltaGF()
                +2.818 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ; 
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121814. * CHbox / LambdaNP2
                +113858. * CHL1_11 / LambdaNP2
                +6221.44 * CHe_11 / LambdaNP2
                -93321.6 * CHu_11 / LambdaNP2
                +113858. * CHL3_11 / LambdaNP2
                -121180. * CuH_33r / LambdaNP2
                -79695. * CHD / LambdaNP2
                +91201.9 * CHB / LambdaNP2
                +178853. * CHW / LambdaNP2
                -314513. * CHWB / LambdaNP2
                -137.642 * CDHB / LambdaNP2
                +853.383 * CDHW / LambdaNP2
                +1906734. * CuW_33r / LambdaNP2
                +1124181. * CuB_33r / LambdaNP2
                -2.642 * DeltaGF()
                -3.21 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ; 
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueettHPol()");
        }
          
    } else if (sqrt_s == 1.0) {
        
        C1 = 0.017;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +122269. * CHbox / LambdaNP2
                +148925. * CHL1_11 / LambdaNP2
                -1516295. * CHe_11 / LambdaNP2
                +181376. * CHu_11 / LambdaNP2
                +148925. * CHL3_11 / LambdaNP2
                -115721. * CuH_33r / LambdaNP2
                -9966.97 * CHD / LambdaNP2
                +648027. * CHB / LambdaNP2
                +58990.6 * CHW / LambdaNP2
                -166947. * CHWB / LambdaNP2
                +258446. * CDHB / LambdaNP2
                +27641. * CDHW / LambdaNP2
                +416063. * CuW_33r / LambdaNP2
                +5771745. * CuB_33r / LambdaNP2
                -0.426 * DeltaGF()
                +3.026 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ;            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +122212. * CHbox / LambdaNP2
                +1266376. * CHL1_11 / LambdaNP2
                -47326.8 * CHe_11 / LambdaNP2
                -104685. * CHu_11 / LambdaNP2
                +1266376. * CHL3_11 / LambdaNP2
                -116193. * CuH_33r / LambdaNP2
                -85861. * CHD / LambdaNP2
                +202732. * CHB / LambdaNP2
                +516612. * CHW / LambdaNP2
                -514723. * CHWB / LambdaNP2
                -75504.5 * CDHB / LambdaNP2
                +158356. * CDHW / LambdaNP2
                +3954267. * CuW_33r / LambdaNP2
                +2288387. * CuB_33r / LambdaNP2
                -2.929 * DeltaGF()
                -5.432 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ; 
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +122564. * CHbox / LambdaNP2
                +252265. * CHL1_11 / LambdaNP2
                -1381101. * CHe_11 / LambdaNP2
                +155161. * CHu_11 / LambdaNP2
                +252265. * CHL3_11 / LambdaNP2
                -115358. * CuH_33r / LambdaNP2
                -16813.1 * CHD / LambdaNP2
                +607466. * CHB / LambdaNP2
                +101359. * CHW / LambdaNP2
                -198737. * CHWB / LambdaNP2
                +227834. * CDHB / LambdaNP2
                +39939.6 * CDHW / LambdaNP2
                +742520. * CuW_33r / LambdaNP2
                +5453267. * CuB_33r / LambdaNP2
                -0.659 * DeltaGF()
                +2.273 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ; 
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +122380. * CHbox / LambdaNP2
                +1238124. * CHL1_11 / LambdaNP2
                -84811.2 * CHe_11 / LambdaNP2
                -97259.2 * CHu_11 / LambdaNP2
                +1238124. * CHL3_11 / LambdaNP2
                -116044. * CuH_33r / LambdaNP2
                -83798.9 * CHD / LambdaNP2
                +214128. * CHB / LambdaNP2
                +505118. * CHW / LambdaNP2
                -505830. * CHWB / LambdaNP2
                -66814.1 * CDHB / LambdaNP2
                +155075. * CDHW / LambdaNP2
                +3863710. * CuW_33r / LambdaNP2
                +2378351. * CuB_33r / LambdaNP2
                -2.867 * DeltaGF()
                -5.212 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ; 
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueettHPol()");
        } 
        
    } else if (sqrt_s == 1.4) {
        
        C1 = 0.0094;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121945. * CHbox / LambdaNP2
                +416437. * CHL1_11 / LambdaNP2
                -5198451. * CHe_11 / LambdaNP2
                +211446. * CHu_11 / LambdaNP2
                +416437. * CHL3_11 / LambdaNP2
                -110413. * CuH_33r / LambdaNP2
                -8089.5 * CHD / LambdaNP2
                +852065. * CHB / LambdaNP2
                +78915.7 * CHW / LambdaNP2
                -191411. * CHWB / LambdaNP2
                +881670. * CDHB / LambdaNP2
                +72289.2 * CDHW / LambdaNP2
                +588296. * CuW_33r / LambdaNP2
                +7812392. * CuB_33r / LambdaNP2
                -0.441 * DeltaGF()
                +2.819 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ;            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +122124. * CHbox / LambdaNP2
                +3668482. * CHL1_11 / LambdaNP2
                -164738. * CHe_11 / LambdaNP2
                -106285. * CHu_11 / LambdaNP2
                +3668482. * CHL3_11 / LambdaNP2
                -112775. * CuH_33r / LambdaNP2
                -87497.2 * CHD / LambdaNP2
                +261266. * CHB / LambdaNP2
                +703789. * CHW / LambdaNP2
                -618584. * CHWB / LambdaNP2
                -257636. * CDHB / LambdaNP2
                +530202. * CDHW / LambdaNP2
                +5501929. * CuW_33r / LambdaNP2
                +3213842. * CuB_33r / LambdaNP2
                -3.038 * DeltaGF()
                -6.378 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ; 
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121843. * CHbox / LambdaNP2
                +706068. * CHL1_11 / LambdaNP2
                -4748505. * CHe_11 / LambdaNP2
                +182964. * CHu_11 / LambdaNP2
                +706068. * CHL3_11 / LambdaNP2
                -110672. * CuH_33r / LambdaNP2
                -15249.5 * CHD / LambdaNP2
                +798771. * CHB / LambdaNP2
                +134415. * CHW / LambdaNP2
                -229663. * CHWB / LambdaNP2
                +779863. * CDHB / LambdaNP2
                +112951. * CDHW / LambdaNP2
                +1026697. * CuW_33r / LambdaNP2
                +7402171. * CuB_33r / LambdaNP2
                -0.673 * DeltaGF()
                +1.996 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ; 
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +122069. * CHbox / LambdaNP2
                +3581543. * CHL1_11 / LambdaNP2
                -298692. * CHe_11 / LambdaNP2
                -97874.3 * CHu_11 / LambdaNP2
                +3581543. * CHL3_11 / LambdaNP2
                -112737. * CuH_33r / LambdaNP2
                -85431.2 * CHD / LambdaNP2
                +276629. * CHB / LambdaNP2
                +687136. * CHW / LambdaNP2
                -607155. * CHWB / LambdaNP2
                -227375. * CDHB / LambdaNP2
                +517945. * CDHW / LambdaNP2
                +5370183. * CuW_33r / LambdaNP2
                +3335906. * CuB_33r / LambdaNP2
                -2.969 * DeltaGF()
                -6.138 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ; 
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueettHPol()");
        } 
        
    } else if (sqrt_s == 1.5) {
        
        C1 = 0.0094;// Use the same as 1400 GeV

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +121854. * CHbox / LambdaNP2
                +507190. * CHL1_11 / LambdaNP2
                -6475118. * CHe_11 / LambdaNP2
                +216935. * CHu_11 / LambdaNP2
                +507190. * CHL3_11 / LambdaNP2
                -109820. * CuH_33r / LambdaNP2
                -7568.59 * CHD / LambdaNP2
                +893094. * CHB / LambdaNP2
                +82781.5 * CHW / LambdaNP2
                -196556. * CHWB / LambdaNP2
                +1099527. * CDHB / LambdaNP2
                +87228. * CDHW / LambdaNP2
                +630747. * CuW_33r / LambdaNP2
                +8328477. * CuB_33r / LambdaNP2
                -0.442 * DeltaGF()
                +2.756 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ;            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +121994. * CHbox / LambdaNP2
                +4501280. * CHL1_11 / LambdaNP2
                -206085. * CHe_11 / LambdaNP2
                -106381. * CHu_11 / LambdaNP2
                +4501280. * CHL3_11 / LambdaNP2
                -112104. * CuH_33r / LambdaNP2
                -87805.6 * CHD / LambdaNP2
                +273106. * CHB / LambdaNP2
                +741955. * CHW / LambdaNP2
                -639545. * CHWB / LambdaNP2
                -322155. * CDHB / LambdaNP2
                +661931. * CDHW / LambdaNP2
                +5892414. * CuW_33r / LambdaNP2
                +3448015. * CuB_33r / LambdaNP2
                -3.057 * DeltaGF()
                -6.552 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ; 
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +121793. * CHbox / LambdaNP2
                +861242. * CHL1_11 / LambdaNP2
                -5919951. * CHe_11 / LambdaNP2
                +188249. * CHu_11 / LambdaNP2
                +861242. * CHL3_11 / LambdaNP2
                -109696. * CuH_33r / LambdaNP2
                -14806.7 * CHD / LambdaNP2
                +837632. * CHB / LambdaNP2
                +141142. * CHW / LambdaNP2
                -235907. * CHWB / LambdaNP2
                +973107. * CDHB / LambdaNP2
                +138331. * CDHW / LambdaNP2
                +1097452. * CuW_33r / LambdaNP2
                +7895510. * CuB_33r / LambdaNP2
                -0.673 * DeltaGF()
                +1.935 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ; 
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +122029. * CHbox / LambdaNP2
                +4394189. * CHL1_11 / LambdaNP2
                -373205. * CHe_11 / LambdaNP2
                -97750.6 * CHu_11 / LambdaNP2
                +4394189. * CHL3_11 / LambdaNP2
                -112024. * CuH_33r / LambdaNP2
                -85643.3 * CHD / LambdaNP2
                +289620. * CHB / LambdaNP2
                +724463. * CHW / LambdaNP2
                -627885. * CHWB / LambdaNP2
                -284076. * CDHB / LambdaNP2
                +646658. * CDHW / LambdaNP2
                +5753330. * CuW_33r / LambdaNP2
                +3578793. * CuB_33r / LambdaNP2
                -2.989 * DeltaGF()
                -6.311 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ; 
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueettHPol()");
        }
             
    } else if (sqrt_s == 3.0) {
        
        C1 = 0.0037;

        if (Pol_em == 80. && Pol_ep == -30.){
            mu += 
                +122442. * CHbox / LambdaNP2
                +3092340. * CHL1_11 / LambdaNP2
                -43264264. * CHe_11 / LambdaNP2
                +259622. * CHu_11 / LambdaNP2
                +3092340. * CHL3_11 / LambdaNP2
                -100510. * CuH_33r / LambdaNP2
                -3230.01 * CHD / LambdaNP2
                +1267548. * CHB / LambdaNP2
                +118886. * CHW / LambdaNP2
                -247164. * CHWB / LambdaNP2
                +7397753. * CDHB / LambdaNP2
                +510206. * CDHW / LambdaNP2
                +1343630. * CuW_33r / LambdaNP2
                +17234081. * CuB_33r / LambdaNP2
                -0.459 * DeltaGF()
                +2.453 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ;            
        } else if (Pol_em == -80. && Pol_ep == 30.){
            mu += 
                +122230. * CHbox / LambdaNP2
                +28686134. * CHL1_11 / LambdaNP2
                -1435177. * CHe_11 / LambdaNP2
                -108195. * CHu_11 / LambdaNP2
                +28686134. * CHL3_11 / LambdaNP2
                -105858. * CuH_33r / LambdaNP2
                -89803.1 * CHD / LambdaNP2
                +381886. * CHB / LambdaNP2
                +1102843. * CHW / LambdaNP2
                -834821. * CHWB / LambdaNP2
                -2237555. * CDHB / LambdaNP2
                +4557030. * CDHW / LambdaNP2
                +12639913. * CuW_33r / LambdaNP2
                +7455995. * CuB_33r / LambdaNP2
                -3.212 * DeltaGF()
                -8.009 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ; 
        } else if (Pol_em == 80. && Pol_ep == 0.){
            mu += 
                +122688. * CHbox / LambdaNP2
                +5271741. * CHL1_11 / LambdaNP2
                -39707692. * CHe_11 / LambdaNP2
                +228729. * CHu_11 / LambdaNP2
                +5271741. * CHL3_11 / LambdaNP2
                -100891. * CuH_33r / LambdaNP2
                -10526.3 * CHD / LambdaNP2
                +1192421. * CHB / LambdaNP2
                +202915. * CHW / LambdaNP2
                -296939. * CHWB / LambdaNP2
                +6582510. * CDHB / LambdaNP2
                +853895. * CDHW / LambdaNP2
                +2303644. * CuW_33r / LambdaNP2
                +16407287. * CuB_33r / LambdaNP2
                -0.693 * DeltaGF()
                +1.565 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ; 
        } else if (Pol_em == -80. && Pol_ep == 0.){
            mu += 
                +121781. * CHbox / LambdaNP2
                +27966374. * CHL1_11 / LambdaNP2
                -2597153. * CHe_11 / LambdaNP2
                -98089.4 * CHu_11 / LambdaNP2
                +27966374. * CHL3_11 / LambdaNP2
                -105885. * CuH_33r / LambdaNP2
                -87600.3 * CHD / LambdaNP2
                +406305. * CHB / LambdaNP2
                +1075086. * CHW / LambdaNP2
                -818808. * CHWB / LambdaNP2
                -1967062. * CDHB / LambdaNP2
                +4442109. * CDHW / LambdaNP2
                +12322125. * CuW_33r / LambdaNP2
                +7728315. * CuB_33r / LambdaNP2
                -3.134 * DeltaGF()
                -7.724 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ; 
        } else {
            throw std::runtime_error("Bad argument in NPSMEFTd6::mueettHPol()");
        } 
    
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::mueettHPol()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += eeettHint + eeettHpar;
    
//  Linear contribution from Higgs self-coupling
    mu = mu + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
    
    if (FlagQuadraticTerms) {
//  Quadratic contribution from Higgs self-coupling
        mu = mu + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();        
    }

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
    
    Br += deltaGammaHggRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        Br += - deltaGammaHggRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHggRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPSMEFTd6::BrHWWRatio() const
{
    
    return BrHWW4fRatio();

}

double NPSMEFTd6::BrHWlvRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHWlvRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        Br += - deltaGammaHWlvRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHWlvRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;
}
    
double NPSMEFTd6::BrHWW2l2vRatio() const
{
    double deltaBRratio;
    
    deltaBRratio = deltaGamma_Wff(leptons[NEUTRINO_1], leptons[ELECTRON]) 
            + deltaGamma_Wff(leptons[NEUTRINO_2], leptons[MU]);
    
    deltaBRratio = deltaBRratio / 
            ( trueSM.GammaW(leptons[NEUTRINO_1], leptons[ELECTRON]) 
            + trueSM.GammaW(leptons[NEUTRINO_2], leptons[MU]) );
    
    deltaBRratio = deltaBRratio - deltaGamma_W() / trueSM.GammaW();
    
    return ( BrHWlvRatio() + deltaBRratio );
}

double NPSMEFTd6::BrHWjjRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHWjjRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        Br += - deltaGammaHWjjRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHWjjRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;
}

double NPSMEFTd6::BrHWW4jRatio() const
{
    double deltaBRratio;
    
    deltaBRratio = deltaGamma_Wff(quarks[UP], quarks[DOWN]) 
            + deltaGamma_Wff(quarks[CHARM], quarks[STRANGE]);
    
    deltaBRratio = deltaBRratio / 
            ( trueSM.GammaW(quarks[UP], quarks[DOWN]) 
            + trueSM.GammaW(quarks[CHARM], quarks[STRANGE]) );
    
    deltaBRratio = deltaBRratio - deltaGamma_W() / trueSM.GammaW();
    
    return ( BrHWjjRatio() + deltaBRratio );
}

double NPSMEFTd6::BrHWffRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHWffRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        Br += - deltaGammaHWffRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHWffRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;
}


double NPSMEFTd6::BrHWW4fRatio() const
{    
    return BrHWffRatio();
}

double NPSMEFTd6::BrHZZRatio() const
{    
    return BrHZZ4fRatio();
}

double NPSMEFTd6::BrHZllRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHZllRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        Br += - deltaGammaHZllRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHZllRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;
}

double NPSMEFTd6::BrHZZ4lRatio() const
{
    double deltaBRratio;
    
    deltaBRratio = deltaGamma_Zf(leptons[ELECTRON]) 
            + deltaGamma_Zf(leptons[MU]);
    
    deltaBRratio = deltaBRratio / 
            ( trueSM.GammaZ(leptons[ELECTRON]) + trueSM.GammaZ(leptons[MU]) );
    
    deltaBRratio = deltaBRratio - deltaGamma_Z() / trueSM.Gamma_Z();
    
    return ( BrHZllRatio() + deltaBRratio );
}

double NPSMEFTd6::BrHZvvRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHZvvRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        Br += - deltaGammaHZvvRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHZvvRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;
}

double NPSMEFTd6::BrHZZ4vRatio() const
{
    double deltaBRratio;
    
    deltaBRratio = deltaGamma_Zf(leptons[NEUTRINO_1]) 
            + deltaGamma_Zf(leptons[NEUTRINO_2])
            + deltaGamma_Zf(leptons[NEUTRINO_3]);
    
    deltaBRratio = deltaBRratio / 
            ( trueSM.GammaZ(leptons[NEUTRINO_1]) 
            + trueSM.GammaZ(leptons[NEUTRINO_2]) 
            + trueSM.GammaZ(leptons[NEUTRINO_3]) );
    
    deltaBRratio = deltaBRratio - deltaGamma_Z() / trueSM.Gamma_Z();
    
    return ( BrHZvvRatio() + deltaBRratio );
}

double NPSMEFTd6::BrHZuuRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHZuuRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        Br += - deltaGammaHZuuRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHZuuRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
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
    
    Br += deltaGammaHZddRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        Br += - deltaGammaHZddRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHZddRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
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
    
    Br += deltaGammaHZffRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        Br += - deltaGammaHZffRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHZffRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;
}

double NPSMEFTd6::BrHZZ4fRatio() const
{    
    return BrHZffRatio();
}

double NPSMEFTd6::BrHZgaRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHZgaRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        Br += - deltaGammaHZgaRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHZgaRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
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

double NPSMEFTd6::BrHgagaRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHgagaRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        Br += - deltaGammaHgagaRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHgagaRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}
      
double NPSMEFTd6::BrHmumuRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHmumuRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        Br += - deltaGammaHmumuRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHmumuRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPSMEFTd6::BrHtautauRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHtautauRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        Br += - deltaGammaHtautauRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHtautauRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPSMEFTd6::BrHccRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHccRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        Br += - deltaGammaHccRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHccRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPSMEFTd6::BrHbbRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHbbRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        Br += - deltaGammaHbbRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHbbRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPSMEFTd6::computeGammaTotalRatio() const
{
    double width = 1.0;

    width += deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaTotalRatio2();
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
    
    double C1 = 0.0066;

    width += deltaGammaHggRatio1();
    
//  Linear contribution from Higgs self-coupling
    width = width + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHggRatio2();  
            
//  Quadratic contribution from Higgs self-coupling
            width = width + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio(); 
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHggRatio1() const
{
    double dwidth = 0.0;
    
    dwidth = ( +37523820. * CHG / LambdaNP2
            + cLHd6 * (
            +121249. * CHbox / LambdaNP2
            -121249. * CHL3_11 / LambdaNP2
            +173342. * CuH_22r / LambdaNP2
            -129305. * CuH_33r / LambdaNP2
            +248515. * CdH_33r / LambdaNP2
            -30312.3 * CHD / LambdaNP2
            -60624.6 * DeltaGF() / v() / v() ) 
            );

    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHggint + eHggpar;
    
    return dwidth;
}

double NPSMEFTd6::deltaGammaHggRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    return ( 0.0 );
    
}

double NPSMEFTd6::GammaHWWRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;
    
    double C1 = 0.0073;

    width += deltaGammaHWWRatio1();
    
//  Linear contribution from Higgs self-coupling
    width = width + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHWWRatio2();
            
//  Quadratic contribution from Higgs self-coupling
            width = width + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio(); 
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHWWRatio1() const
{    
    double dwidth = 0.0;

    dwidth = deltaGammaHWffRatio1();
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHWWint + eHWWpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHWWRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    return deltaGammaHWffRatio2();
    
}

double NPSMEFTd6::GammaHWlvRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;
    
    double C1 = 0.0073;

    width += deltaGammaHWlvRatio1();
    
//  Linear contribution from Higgs self-coupling
    width = width + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHWlvRatio2();
            
//  Quadratic contribution from Higgs self-coupling
            width = width + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio(); 
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHWlvRatio1() const
{    
    double dwidth = 0.0;

    dwidth = ( +121875. * CHbox / LambdaNP2
                +18351.9 * (1.0/2.0) * ( CHL3_11 + CHL3_22 ) / LambdaNP2
                -159873. * CHD / LambdaNP2
                -91288.7 * CHW / LambdaNP2
                -284689. * CHWB / LambdaNP2
                +37703.7 * CDHW / LambdaNP2
                -3.292 * DeltaGF()
                -15.14 * deltaMw()
                +115.287 * deltaMw2()  );
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHWWint + eHWWpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHWlvRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    return ( 0.0 );
    
}

double NPSMEFTd6::GammaHWjjRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;
    
    double C1 = 0.0073;

    width += deltaGammaHWjjRatio1();
    
//  Linear contribution from Higgs self-coupling
    width = width + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHWjjRatio2();
            
//  Quadratic contribution from Higgs self-coupling
            width = width + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio(); 
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHWjjRatio1() const
{    
    double dwidth = 0.0;

    dwidth = ( +121611. * CHbox / LambdaNP2
                +17701.4 * (1.0/2.0) * ( CHQ3_11 + CHQ3_22 ) / LambdaNP2
                -159273. * CHD / LambdaNP2
                -91021.6 * CHW / LambdaNP2
                -282574. * CHWB / LambdaNP2
                +37917.5 * CDHW / LambdaNP2
                -3.259 * DeltaGF()
                -15.198 * deltaMw()
                +121.807 * deltaMw2()  );
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHWWint + eHWWpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHWjjRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    return ( 0.0 );
    
}

double NPSMEFTd6::GammaHWffRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;
    
    double C1 = 0.0073;

    width += deltaGammaHWffRatio1();
    
//  Linear contribution from Higgs self-coupling
    width = width + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHWffRatio2();
            
//  Quadratic contribution from Higgs self-coupling
            width = width + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio(); 
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHWffRatio1() const
{    
    double dwidth = 0.0;

    dwidth = ( +121288. * CHbox / LambdaNP2
                +5395.21 * (1.0/3.0) * ( CHL3_11 + CHL3_22 + CHL3_33 ) / LambdaNP2
                +11680.9 * (1.0/2.0) * ( CHQ3_11 + CHQ3_22 ) / LambdaNP2
                -159787. * CHD / LambdaNP2
                -91509.1 * CHW / LambdaNP2
                -283092. * CHWB / LambdaNP2
                +37845.1 * CDHW / LambdaNP2
                -3.259 * DeltaGF()
                -15.196 * deltaMw()
                +106.857 * deltaMw2()  );
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHWWint + eHWWpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHWffRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    return ( 0.0 );
    
}

double NPSMEFTd6::GammaHZZRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;
    
    double C1 = 0.0083;

    width += deltaGammaHZZRatio1();
    
//  Linear contribution from Higgs self-coupling
    width = width + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHZZRatio2();   
            
//  Quadratic contribution from Higgs self-coupling
            width = width + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHZZRatio1() const
{
    double dwidth = 0.0;

    dwidth = deltaGammaHZffRatio1();
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHZZint + eHZZpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHZZRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    return deltaGammaHZffRatio2();
    
}

double NPSMEFTd6::GammaHZllRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;
    
    double C1 = 0.0083;

    width += deltaGammaHZllRatio1();
    
//  Linear contribution from Higgs self-coupling
    width = width + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHZllRatio2();   
            
//  Quadratic contribution from Higgs self-coupling
            width = width + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHZllRatio1() const
{
    double dwidth = 0.0;

    dwidth = ( +121715. * CHbox / LambdaNP2
                +8726.9 * (1.0/2.0) * ( CHL1_11 + CHL1_22 ) / LambdaNP2
                -7315.2 * (1.0/2.0) * ( CHe_11 + CHe_22 ) / LambdaNP2
                +8726.9 * (1.0/2.0) * ( CHL3_11 + CHL3_22 ) / LambdaNP2
                -5544.15 * CHD / LambdaNP2
                -13560.9 * CHB / LambdaNP2
                -45585.2 * CHW / LambdaNP2
                -53507.9 * CHWB / LambdaNP2
                +16829.2 * CDHB / LambdaNP2
                +30766.6 * CDHW / LambdaNP2
                -2.204 * DeltaGF()  );
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHZZint + eHZZpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHZllRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    return ( 0.0 );
    
}

double NPSMEFTd6::GammaHZvvRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;
    
    double C1 = 0.0083;

    width += deltaGammaHZvvRatio1();
    
//  Linear contribution from Higgs self-coupling
    width = width + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHZvvRatio2();   
            
//  Quadratic contribution from Higgs self-coupling
            width = width + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHZvvRatio1() const
{
    double dwidth = 0.0;

    dwidth = ( +121530. * CHbox / LambdaNP2
                -7943.34 * (1.0/3.0) * ( CHL1_11 + CHL1_22 + CHL1_33 ) / LambdaNP2
                +7943.34 * (1.0/3.0) * ( CHL3_11 + CHL3_22 + CHL3_33 ) / LambdaNP2
                -229.41 * CHD / LambdaNP2
                -13535.2 * CHB / LambdaNP2
                -45480.6 * CHW / LambdaNP2
                -24891. * CHWB / LambdaNP2
                +16833. * CDHB / LambdaNP2
                +30597.6 * CDHW / LambdaNP2
                -2. * DeltaGF()  );
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHZZint + eHZZpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHZvvRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    return ( 0.0 );
    
}

double NPSMEFTd6::GammaHZuuRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;
    
    double C1 = 0.0083;

    width += deltaGammaHZuuRatio1();
    
//  Linear contribution from Higgs self-coupling
    width = width + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHZuuRatio2();   
            
//  Quadratic contribution from Higgs self-coupling
            width = width + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHZuuRatio1() const
{
    double dwidth = 0.0;

    dwidth = ( +121512. * CHbox / LambdaNP2
                -9648.28 * (1.0/2.0) * ( CHQ1_11 + CHQ1_22 ) / LambdaNP2
                +4218.6 * (1.0/2.0) * ( CHu_11 + CHu_22 ) / LambdaNP2
                +9648.28 * (1.0/2.0) * ( CHQ3_11 + CHQ3_22 ) / LambdaNP2
                -17762.5 * CHD / LambdaNP2
                -13473.3 * CHB / LambdaNP2
                -45667.9 * CHW / LambdaNP2
                -110057. * CHWB / LambdaNP2
                +16854.2 * CDHB / LambdaNP2
                +30781.7 * CDHW / LambdaNP2
                -2.6 * DeltaGF()  );
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHZZint + eHZZpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHZuuRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    return ( 0.0 );
    
}

double NPSMEFTd6::GammaHZddRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;
    
    double C1 = 0.0083;

    width += deltaGammaHZddRatio1();
    
//  Linear contribution from Higgs self-coupling
    width = width + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHZddRatio2();   
            
//  Quadratic contribution from Higgs self-coupling
            width = width + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHZddRatio1() const
{
    double dwidth = 0.0;

    dwidth = ( +121756. * CHbox / LambdaNP2
                +9252.73 * (1.0/3.0) * ( CHQ1_11 + CHQ1_22 + CHQ1_33 ) / LambdaNP2
                -1471.03 * (1.0/3.0) * ( CHd_11 + CHd_22 + CHd_33 ) / LambdaNP2
                +9252.73 * (1.0/3.0) * ( CHQ3_11 + CHQ3_22 + CHQ3_33 ) / LambdaNP2
                -12714.3 * CHD / LambdaNP2
                -13589.3 * CHB / LambdaNP2
                -45689.4 * CHW / LambdaNP2
                -85582.3 * CHWB / LambdaNP2
                +17007.1 * CDHB / LambdaNP2
                +30733.1 * CDHW / LambdaNP2
                -2.427 * DeltaGF()  );
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHZZint + eHZZpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHZddRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    return ( 0.0 );
    
}

double NPSMEFTd6::GammaHZffRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;
    
    double C1 = 0.0083;

    width += deltaGammaHZffRatio1();
    
//  Linear contribution from Higgs self-coupling
    width = width + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHZffRatio2();   
            
//  Quadratic contribution from Higgs self-coupling
            width = width + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHZffRatio1() const
{
    double dwidth = 0.0;

    dwidth = ( +121551. * CHbox / LambdaNP2
                -824.482 * (1.0/3.0) * ( CHL1_11 + CHL1_22 + CHL1_33 ) / LambdaNP2
                +1840.54 * (1.0/12.0) * ( 5.0 * CHQ1_11 + 5.0 * CHQ1_22 + 2.0 * CHQ1_33 - CHQ3_11 - CHQ3_22 + 2.0 * CHQ3_33 ) / LambdaNP2
                -795.383 * (1.0/3.0) * ( CHe_11 + CHe_22 + CHe_33 ) / LambdaNP2
                +1069.4 * (1.0/2.0) * ( CHu_11 + CHu_22 ) / LambdaNP2
                -579.563 * (1.0/3.0) * ( CHd_11 + CHd_22 + CHd_33 ) / LambdaNP2
                +3164.56 * (1.0/3.0) * ( CHL3_11 + CHL3_22 + CHL3_33 ) / LambdaNP2
                +6413.99 * (-1.0/12.0) * ( CHQ1_11 + CHQ1_22 - 2.0 * CHQ1_33 - 5.0 * CHQ3_11 - 5.0 * CHQ3_22 - 2.0 * CHQ3_33) / LambdaNP2
                -10839.5 * CHD / LambdaNP2
                -14222.3 * CHB / LambdaNP2
                -45455.6 * CHW / LambdaNP2
                -75343.1 * CHWB / LambdaNP2
                +16804.9 * CDHB / LambdaNP2
                +30421. * CDHW / LambdaNP2
                -2.356 * DeltaGF()  );
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHZZint + eHZZpar;
    
    return dwidth;
    
}

double NPSMEFTd6::deltaGammaHZffRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    return ( 0.0 );
    
}

double NPSMEFTd6::GammaHZgaRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;
    
    double C1 = 0.0;

    width += deltaGammaHZgaRatio1();
    
//  Linear contribution from Higgs self-coupling
    width = width + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHZgaRatio2();      
            
//  Quadratic contribution from Higgs self-coupling
            width = width + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHZgaRatio1() const
{
    double dwidth = 0.0;

//  It does not include modifications of Zff vertices or MW
    dwidth = ( +14894896. * CHB / LambdaNP2
            -14894896. * CHW / LambdaNP2
            +9508500. * CHWB / LambdaNP2
            -2869606. * CDHB / LambdaNP2
            +1572591. * CDHW / LambdaNP2
            + cLHd6 * (
            +120004. * CHbox / LambdaNP2
            -374.02 * CeH_33r / LambdaNP2
            -2953.25 * CuH_22r / LambdaNP2
            +6645.39 * CuH_33r / LambdaNP2
            -6121.81 * CdH_33r / LambdaNP2
            -120353. * CHD / LambdaNP2
            -198058. * CHWB / LambdaNP2
            -114271. * DeltaGF() / v() / v() )
             );
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHZgaint + eHZgapar;
    
    return dwidth;
}

double NPSMEFTd6::deltaGammaHZgaRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    return ( 0.0 );
    
}

double NPSMEFTd6::GammaHgagaRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;
    
    double C1 = 0.0049;

    width += deltaGammaHgagaRatio1();
    
//  Linear contribution from Higgs self-coupling
    width = width + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHgagaRatio2();  
            
//  Quadratic contribution from Higgs self-coupling
            width = width + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio(); 
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHgagaRatio1() const
{
    double dwidth = 0.0;

//  It does not include modifications of MW
    dwidth = ( -48316360. * CHB / LambdaNP2
            -14510470. * CHW / LambdaNP2
            +26478163. * CHWB / LambdaNP2 
            + cLHd6 * (
            +119767. * CHbox / LambdaNP2
            -42566.9 * CeH_33r / LambdaNP2
            -48869.9 * CuH_22r / LambdaNP2
            +32119.8 * CuH_33r / LambdaNP2
            -18429. * CdH_33r / LambdaNP2
            -137455. * CHD / LambdaNP2
            -235675. * CHWB / LambdaNP2
            -124461. * DeltaGF() / v() / v() )
            );
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHgagaint + eHgagapar;
    
    return dwidth;
}

double NPSMEFTd6::deltaGammaHgagaRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    return ( 0.0 );
    
}
      
double NPSMEFTd6::GammaHmumuRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;
    
    double C1 = 0.0;

    width += deltaGammaHmumuRatio1();
    
//  Linear contribution from Higgs self-coupling
    width = width + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHmumuRatio2();
            
//  Quadratic contribution from Higgs self-coupling
            width = width + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio(); 
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHmumuRatio1() const
{
    double dwidth = 0.0;
    
    dwidth = ( +121249. * CHbox / LambdaNP2
            -199794752. * CeH_22r / LambdaNP2
            -30312.3 * CHD / LambdaNP2
            -60624.6 * DeltaGF() / v() / v() ); 
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHmumuint + eHmumupar;
    
    return dwidth;
}

double NPSMEFTd6::deltaGammaHmumuRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    return ( 0.0 );
    
}

double NPSMEFTd6::GammaHtautauRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;
    
    double C1 = 0.0;

    width += deltaGammaHtautauRatio1();
    
//  Linear contribution from Higgs self-coupling
    width = width + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHtautauRatio2();    
            
//  Quadratic contribution from Higgs self-coupling
            width = width + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHtautauRatio1() const
{
    double dwidth = 0.0;

    dwidth = ( +121249. * CHbox / LambdaNP2
            -11880769. * CeH_33r / LambdaNP2
            -30312.3 * CHD / LambdaNP2
            -60624.6 * DeltaGF() / v() / v() );
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHtautauint + eHtautaupar;
    
    return dwidth;
}

double NPSMEFTd6::deltaGammaHtautauRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    return ( 0.0 );
    
}

double NPSMEFTd6::GammaHccRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;
    
    double C1 = 0.0;

    width += deltaGammaHccRatio1();
    
//  Linear contribution from Higgs self-coupling
    width = width + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHccRatio2();   
            
//  Quadratic contribution from Higgs self-coupling
            width = width + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio();
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHccRatio1() const
{
    double dwidth = 0.0;
    
    if (FlagLoopHd6) {
        
        dwidth = ( +121249. * CHbox / LambdaNP2
            -16422087. * CuH_22r / LambdaNP2
            -993.221 * CuH_33r / LambdaNP2
            -30312.3 * CHD / LambdaNP2
            -60624.6 * DeltaGF() / v() / v() );
        
    } else {
        
        dwidth = ( +121249. * CHbox / LambdaNP2
            -16556853. * CuH_22r / LambdaNP2
            -30312.3 * CHD / LambdaNP2
            -60624.6 * DeltaGF() / v() / v() );  
    }
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHccint + eHccpar;
    
    return dwidth;
}

double NPSMEFTd6::deltaGammaHccRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    return ( 0.0 );
    
}

double NPSMEFTd6::GammaHbbRatio() const
{
      // SM (1). Intrinsic + parametric theory relative errors (free pars) included in deltaGammaHXXRatio1
    double width = 1.0;
    
    double C1 = 0.0;
    
    width += deltaGammaHbbRatio1();
    
//  Linear contribution from Higgs self-coupling
    width = width + cLHd6*(C1 + 2.0*dZH)*deltaG_hhhRatio();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHbbRatio2();    
            
//  Quadratic contribution from Higgs self-coupling
            width = width + cLHd6*dZH*deltaG_hhhRatio()*deltaG_hhhRatio(); 
        }
    
    return width;
}

double NPSMEFTd6::deltaGammaHbbRatio1() const
{
    double dwidth = 0.0;
    
    
    if (FlagLoopHd6) {
        
        dwidth = ( +121249. * CHbox / LambdaNP2
            -558.741 * CuH_33r / LambdaNP2
            -5027111. * CdH_33r / LambdaNP2
            -30312.3 * CHD / LambdaNP2
            -60624.6 * DeltaGF() / v() / v() );
        
    } else {
        
        dwidth = ( +121249. * CHbox / LambdaNP2
            -5050236. * CdH_33r / LambdaNP2
            -30312.3 * CHD / LambdaNP2
            -60624.6 * DeltaGF() / v() / v() );  
    }
    
    // SM (1) + intrinsic + parametric theory relative errors (free pars)    
    dwidth += eHbbint + eHbbpar;
    
    return dwidth;
}

double NPSMEFTd6::deltaGammaHbbRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    return ( 0.0 );        
    
}

double NPSMEFTd6::Br_H_exo() const
{    
    return BrHexo;
}

double NPSMEFTd6::Br_H_inv() const
{    
    return BrHinv;
}


double NPSMEFTd6::BrHvisRatio() const
{    
    double Br = 1.0;
    double dvis1,dvis2,delta2SM;

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
    
    Br += dvis1 - deltaGammaTotalRatio1();
    
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
        Br += - dvis1 * deltaGammaTotalRatio1()
                + dvis2 - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;
}


///////////////////////SPECIAL OBSERVABLES/////////////////////////

double NPSMEFTd6::muttHZbbboost(const double sqrt_s) const
{
    /* Ratios of BR with the SM*/
    double BrHbbrat = BrHbbRatio();
    double BrZbbSM = (trueSM.GammaZ(quarks[BOTTOM]))/trueSM.Gamma_Z();
    double BrZbbrat = BR_Zf(quarks[BOTTOM])/BrZbbSM;

    gslpp::complex dKappa_t = deltaG_hff(quarks[TOP]) / (-mtpole / v());    
    double dkt = dKappa_t.real();
    
    double dgV = deltaGV_f(quarks[TOP]);
    double dgA = deltaGA_f(quarks[TOP]);
    double gLSM = quarks[TOP].getIsospin() 
    - (quarks[TOP].getCharge())*sW2_tree;
    double gRSM = - (quarks[TOP].getCharge())*sW2_tree;
    
    double dgL = 0.5*(dgV + dgA)/gLSM;
    double dgR = 0.5*(dgV - dgA)/gRSM;

    /* VERY CRUDE APPROX. */
    double dsigmarat;
    
    dsigmarat = 1.0 + 
            2.0 * dkt -
            2.0 * (gLSM*gLSM*dgL + gRSM*gRSM*dgR)/(gLSM*gLSM + gRSM*gRSM);
    
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

///////////////////////////////////////////////////////////////////////////////

double NPSMEFTd6::deltag1ZNP() const
{
      double NPdirect, NPindirect;
      
      /*    From own calculations. Agrees with with LHCHXWG-INT-2015-001 for common interactions */
      NPdirect = sW_tree / sqrt( 4.0 * M_PI * aleMz );
      NPdirect = - NPdirect * (Mz * Mz / v () / v() ) * CDHW * v2_over_LambdaNP2;
      
      NPindirect = - 1.0 / (cW2_tree-sW2_tree);
      
      NPindirect = NPindirect * (sW_tree * CHWB / cW_tree 
              + 0.25 * CHD ) * v2_over_LambdaNP2
              + 0.5 * NPindirect * DeltaGF() ;
      
      return NPdirect + NPindirect + dg1Z ;
}
      
double NPSMEFTd6::deltaKgammaNP() const
{
      double NPdirect;

      /*    Translate from LHCHXWG-INT-2015-001: Checked with own calculations  */
      NPdirect = sqrt( 4.0 * M_PI * aleMz ) / 4.0 / sW2_tree;
      
      NPdirect = NPdirect * ( (4.0 * sW_tree * cW_tree / sqrt( 4.0 * M_PI * aleMz ) ) * CHWB 
              - sW_tree * CDHW 
              - cW_tree * CDHB ) * v2_over_LambdaNP2;
      
      return NPdirect + dKappaga ;
}
      
double NPSMEFTd6::lambdaZNP() const
{
      double NPdirect;

      /*    Translate from LHCHXWG-INT-2015-001: Checked with own calculations  */
      NPdirect = - (3.0 / 2.0) * (sqrt( 4.0 * M_PI * aleMz ) / sW_tree) * CW * v2_over_LambdaNP2;

      return NPdirect + lambZ ;
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
    double g1Z,g1ga,kZ,kga,lambZ,lambga,g4Z,g4ga,g5Z,g5ga,ktZ,ktga,lambtZ,lambtga;

//  TGC present in the SM     
    g1Z=1.0 + deltag1ZNP();
    g1ga=1.0;
    kZ=1.0 + deltag1ZNP() - (sW2_tree/cW2_tree) * deltaKgammaNP();
    kga=1.0 + deltaKgammaNP();
//  TGC not present in the SM
    lambZ=lambdaZNP(); //Check normalization
    lambga=lambZ;
    g4Z=0.0;
    g4ga=0.0;
    g5Z=0.0;
    g5ga=0.0;
    ktZ=0.0;
    ktga=0.0;
    lambtZ=0.0;
    lambtga=0.0;
    
    double f3Z, f3ga;
    
    f3Z = g1Z + kZ + lambZ;    
    f3ga = g1ga + kga + lambga;
    
 // Kinematic factors
    double beta, gamma, gamma2;
    
    beta = sqrt(1.0 - 4.0 * mw * mw / s);
    gamma = sqrt_sGeV/(2.0 * mw);
    gamma2= gamma*gamma;
    
//  J=1 Subamplitudes: Z
    gslpp::complex AZpp, AZmm, AZp0, AZm0, AZ0p, AZ0m, AZ00;
    
    AZpp = gslpp::complex(g1Z + 2.0* gamma2* lambZ, (ktZ + lambtZ - 2.0*lambtZ)/beta , false);
    AZmm = gslpp::complex(g1Z + 2.0* gamma2* lambZ, -(ktZ + lambtZ - 2.0*lambtZ)/beta , false);
    AZp0 = gslpp::complex(f3Z + beta * g5Z , -g4Z + (ktZ-lambtZ)/beta , false);
    AZp0 = gamma * AZp0;
    AZm0 = gslpp::complex(f3Z - beta * g5Z , -g4Z - (ktZ-lambtZ)/beta , false);
    AZm0 = gamma * AZm0;
    AZ0p = gslpp::complex(f3Z - beta * g5Z , g4Z + (ktZ-lambtZ)/beta , false);
    AZ0p = gamma * AZ0p;
    AZ0m = gslpp::complex(f3Z + beta * g5Z , g4Z - (ktZ-lambtZ)/beta , false);
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
    
    Agapp = gslpp::complex(g1ga + 2.0* gamma2* lambga, (ktga + lambtga - 2.0*lambtga)/beta , false);
    Agamm = gslpp::complex(g1ga + 2.0* gamma2* lambga, -(ktga + lambtga - 2.0*lambtga)/beta , false);
    Agap0 = gslpp::complex(f3ga + beta * g5ga , -g4ga + (ktga-lambtga)/beta , false);
    Agap0 = gamma * Agap0;
    Agam0 = gslpp::complex(f3ga - beta * g5ga , -g4ga - (ktga-lambtga)/beta , false);
    Agam0 = gamma * Agam0;
    Aga0p = gslpp::complex(f3ga - beta * g5ga , g4ga + (ktga-lambtga)/beta , false);
    Aga0p = gamma * Aga0p;
    Aga0m = gslpp::complex(f3ga + beta * g5ga , g4ga - (ktga-lambtga)/beta , false);
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

///////////////////////////////////////////////////////////////////////////////

double NPSMEFTd6::CLL_mu() const
{
    return (CLL_1122 + CLL_2211 + CLL_1221 + CLL_2112);
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
