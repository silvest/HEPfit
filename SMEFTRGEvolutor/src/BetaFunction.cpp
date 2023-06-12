

int RGESolver::func(double logmu, const double y[], double f[], void* params) {
    int i, j, a, b, l, k, d, w, count, v;
    int p, r, s, t;
    double loop_factor = 1.
            / (16. * 3.1415926535 * 3.1415926535);
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

    double yeR[NG][NG], yeI[NG][NG], //yukawas
            ydR[NG][NG], ydI[NG][NG],
            yuR[NG][NG], yuI[NG][NG];
    double yedagR[NG][NG], yedagI[NG][NG],
            yddagR[NG][NG], yddagI[NG][NG],
            yudagR[NG][NG], yudagI[NG][NG];

    for (i = 0; i < NG; i ++) {
        for (j = 0; j < NG; j ++) {
            a = 0;
            yuR[i][j] = y[c];
            yudagR[j][i] = y[c];
            a ++;
            yuI[i][j] = y[c + DF];
            yudagI[j][i] = - y[c + DF];
            a ++;
            ydR[i][j] = y[c + 2 * DF];
            yddagR[j][i] = y[c + 2 * DF];
            a ++;
            ydI[i][j] = y[c + 3 * DF];
            yddagI[j][i] = - y[c + 3 * DF];
            a ++;
            yeR[i][j] = y[c + 4 * DF];
            yedagR[j][i] = y[c + 4 * DF];
            a ++;
            yeI[i][j] = y[c + 5 * DF];
            yedagI[j][i] = - y[c + 5 * DF];
            a ++;
            c ++;
        }
    }
    c += (2. * Nyukawa - 1.) * DF;



    //SMEFT class 1 
    double CG = y[c];
    double CGT = y[c + 1];
    double CW = y[c + 2];
    double CWT = y[c + 3];
    c += N1;

    //SMEFT class 2-3 
    double CH = y[c];
    double CHBOX = y[c + 1];
    double CHD = y[c + 2];
    c += N23;
    //SMEFT class 4 
    double CHG = y[c];
    double CHB = y[c + 1];
    double CHW = y[c + 2];
    double CHWB = y[c + 3];

    double CHGT = y[c + 4];
    double CHBT = y[c + 5];
    double CHWT = y[c + 6];
    double CHWBT = y[c + 7];
    c += N4;
    //SMEFT class 5
    double CeHR[NG * NG], CeHI[NG * NG], CuHR[NG * NG],
            CuHI[NG * NG], CdHR[NG * NG], CdHI[NG * NG];

    for (i = 0; i < NG; i ++) {
        for (j = 0; j < NG; j ++) {
            a = 0;
            WC1_set(CuHR, i, j, y[c + a * DF]);
            a ++;
            WC1_set(CuHI, i, j, y[c + a * DF]);
            a ++;
            WC1_set(CdHR, i, j, y[c + a * DF]);
            a ++;
            WC1_set(CdHI, i, j, y[c + a * DF]);
            a ++;
            WC1_set(CeHR, i, j, y[c + a * DF]);
            a ++;
            WC1_set(CeHI, i, j, y[c + a * DF]);
            a ++;
            c ++;
        }
    }
    c += (N5 * 2 - 1) * DF;


    //SMEFT class 6 
    double CeWR[NG * NG], CeWI[NG * NG];
    double CeBR[NG * NG], CeBI[NG * NG];
    double CuGR[NG * NG], CuGI[NG * NG];
    double CuWR[NG * NG], CuWI[NG * NG], CuBR[NG * NG], CuBI[NG * NG];
    double CdGR[NG * NG], CdGI[NG * NG], CdWR[NG * NG], CdWI[NG * NG],
            CdBR[NG * NG], CdBI[NG * NG];
    for (i = 0; i < NG; i ++) {
        for (j = 0; j < NG; j ++) {
            a = 0;
            WC1_set(CeWR, i, j, y[c + a * DF]);
            a ++;
            WC1_set(CeWI, i, j, y[c + a * DF]);
            a ++;
            WC1_set(CeBR, i, j, y[c + a * DF]);
            a ++;
            WC1_set(CeBI, i, j, y[c + a * DF]);
            a ++;
            WC1_set(CuGR, i, j, y[c + a * DF]);
            a ++;
            WC1_set(CuGI, i, j, y[c + a * DF]);
            a ++;
            WC1_set(CuWR, i, j, y[c + a * DF]);
            a ++;
            WC1_set(CuWI, i, j, y[c + a * DF]);
            a ++;
            WC1_set(CuBR, i, j, y[c + a * DF]);
            a ++;
            WC1_set(CuBI, i, j, y[c + a * DF]);
            a ++;
            WC1_set(CdGR, i, j, y[c + a * DF]);
            a ++;
            WC1_set(CdGI, i, j, y[c + a * DF]);
            a ++;
            WC1_set(CdWR, i, j, y[c + a * DF]);
            a ++;
            WC1_set(CdWI, i, j, y[c + a * DF]);
            a ++;
            WC1_set(CdBR, i, j, y[c + a * DF]);
            a ++;
            WC1_set(CdBI, i, j, y[c + a * DF]);
            a ++;
            c ++;
        }
    }
    c += (N6 * 2 - 1) * DF;

    //SMEFT class 7
    double CHl1R[NG * NG], CHl3R[NG * NG], CHl1I[NG * NG], CHl3I[NG * NG];
    double CHeR[NG * NG], CHeI[NG * NG];
    double CHq1R[NG * NG], CHq1I[NG * NG], CHq3R[NG * NG], CHq3I[NG * NG];
    double CHuR[NG * NG], CHuI[NG * NG], CHdR[NG * NG], CHdI[NG * NG];
    double CHudR[NG * NG], CHudI[NG * NG];

    {
        for (i = 0; i < DWC2R; i ++) {
            WC2R_set(CHl1R, WC2R_indices[i][0], WC2R_indices[i][1], y[c]);
            c ++;
        }
        for (i = 0; i < DWC2I; i ++) {
            WC2I_set(CHl1I, WC2I_indices[i][0], WC2I_indices[i][1], y[c]);
            c ++;
        }
        for (i = 0; i < DWC2R; i ++) {
            WC2R_set(CHl3R, WC2R_indices[i][0], WC2R_indices[i][1], y[c]);
            c ++;
        }
        for (i = 0; i < DWC2I; i ++) {
            WC2I_set(CHl3I, WC2I_indices[i][0], WC2I_indices[i][1], y[c]);
            c ++;
        }
        for (i = 0; i < DWC2R; i ++) {
            WC2R_set(CHeR, WC2R_indices[i][0], WC2R_indices[i][1], y[c]);
            c ++;
        }
        for (i = 0; i < DWC2I; i ++) {
            WC2I_set(CHeI, WC2I_indices[i][0], WC2I_indices[i][1], y[c]);
            c ++;
        }

        for (i = 0; i < DWC2R; i ++) {
            WC2R_set(CHq1R, WC2R_indices[i][0], WC2R_indices[i][1], y[c]);
            c ++;
        }
        for (i = 0; i < DWC2I; i ++) {
            WC2I_set(CHq1I, WC2I_indices[i][0], WC2I_indices[i][1], y[c]);
            c ++;
        }
        for (i = 0; i < DWC2R; i ++) {
            WC2R_set(CHq3R, WC2R_indices[i][0], WC2R_indices[i][1], y[c]);
            c ++;
        }
        for (i = 0; i < DWC2I; i ++) {
            WC2I_set(CHq3I, WC2I_indices[i][0], WC2I_indices[i][1], y[c]);
            c ++;
        }

        for (i = 0; i < DWC2R; i ++) {
            WC2R_set(CHuR, WC2R_indices[i][0], WC2R_indices[i][1], y[c]);
            c ++;
        }
        for (i = 0; i < DWC2I; i ++) {
            WC2I_set(CHuI, WC2I_indices[i][0], WC2I_indices[i][1], y[c]);
            c ++;
        }

        for (i = 0; i < DWC2R; i ++) {
            WC2R_set(CHdR, WC2R_indices[i][0], WC2R_indices[i][1], y[c]);
            c ++;
        }
        for (i = 0; i < DWC2I; i ++) {
            WC2I_set(CHdI, WC2I_indices[i][0], WC2I_indices[i][1], y[c]);
            c ++;
        }

        for (i = 0; i < NG; i ++) {
            for (j = 0; j < NG; j ++) {
                WC1_set(CHudR, i, j, y[c]);
                c ++;
            }
        }
        for (i = 0; i < NG; i ++) {
            for (j = 0; j < NG; j ++) {
                WC1_set(CHudI, i, j, y[c]);
                c ++;
            }
        }
    }

    //SMEFT class 8_LLLL
    double CllR[NG * NG * NG * NG], CllI[NG * NG * NG * NG], Cqq1R[NG * NG * NG * NG], Cqq1I[NG * NG * NG * NG],
            Cqq3R[NG * NG * NG * NG], Cqq3I[NG * NG * NG * NG];

    double Clq1R[NG * NG * NG * NG];
    double Clq1I[NG * NG * NG * NG];
    double Clq3R[NG * NG * NG * NG];
    double Clq3I[NG * NG * NG * NG];

    for (a = 0; a < DWC6R; a ++) {
        WC6R_set(CllR, WC6R_indices[a][0], WC6R_indices[a][1],
                WC6R_indices[a][2], WC6R_indices[a][3], y[c]);
        c ++;
    }
    for (a = 0; a < DWC6I; a ++) {
        WC6I_set(CllI, WC6I_indices[a][0], WC6I_indices[a][1],
                WC6I_indices[a][2], WC6I_indices[a][3], y[c]);
        c ++;
    }
    for (a = 0; a < DWC6R; a ++) {
        WC6R_set(Cqq1R, WC6R_indices[a][0], WC6R_indices[a][1],
                WC6R_indices[a][2], WC6R_indices[a][3], y[c]);
        c ++;
    }
    for (a = 0; a < DWC6I; a ++) {
        WC6I_set(Cqq1I, WC6I_indices[a][0], WC6I_indices[a][1],
                WC6I_indices[a][2], WC6I_indices[a][3], y[c]);
        c ++;
    }
    for (a = 0; a < DWC6R; a ++) {
        WC6R_set(Cqq3R, WC6R_indices[a][0], WC6R_indices[a][1],
                WC6R_indices[a][2], WC6R_indices[a][3], y[c]);
        c ++;
    }
    for (a = 0; a < DWC6I; a ++) {
        WC6I_set(Cqq3I, WC6I_indices[a][0], WC6I_indices[a][1],
                WC6I_indices[a][2], WC6I_indices[a][3], y[c]);
        c ++;
    }
    for (a = 0; a < DWC7R; a ++) {
        WC7R_set(Clq1R, WC7R_indices[a][0], WC7R_indices[a][1],
                WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
        c ++;
    }
    for (a = 0; a < DWC7I; a ++) {
        WC7I_set(Clq1I, WC7I_indices[a][0], WC7I_indices[a][1],
                WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
        c ++;
    }
    for (a = 0; a < DWC7R; a ++) {
        WC7R_set(Clq3R, WC7R_indices[a][0], WC7R_indices[a][1],
                WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
        c ++;
    }
    for (a = 0; a < DWC7I; a ++) {
        WC7I_set(Clq3I, WC7I_indices[a][0], WC7I_indices[a][1],
                WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
        c ++;
    }

    //SMEFT class 8_RRRR
    double CuuR[NG * NG * NG * NG];
    double CuuI[NG * NG * NG * NG];
    double CeeR[NG * NG * NG * NG];
    double CeeI[NG * NG * NG * NG];
    double CddR[NG * NG * NG * NG];
    double CddI[NG * NG * NG * NG];
    double CeuR[NG * NG * NG * NG];
    double CeuI[NG * NG * NG * NG];
    double CedR[NG * NG * NG * NG];
    double CedI[NG * NG * NG * NG];
    double Cud1R[NG * NG * NG * NG];
    double Cud1I[NG * NG * NG * NG];
    double Cud8R[NG * NG * NG * NG];
    double Cud8I[NG * NG * NG * NG];



    for (a = 0; a < DWC8R; a ++) {
        WC8R_set(CeeR, WC8R_indices[a][0], WC8R_indices[a][1],
                WC8R_indices[a][2], WC8R_indices[a][3], y[c]);
        c ++;
    }
    for (a = 0; a < DWC8I; a ++) {
        WC8I_set(CeeI, WC8I_indices[a][0], WC8I_indices[a][1],
                WC8I_indices[a][2], WC8I_indices[a][3], y[c]);
        c ++;
    }

    for (a = 0; a < DWC6R; a ++) {
        WC6R_set(CuuR, WC6R_indices[a][0], WC6R_indices[a][1],
                WC6R_indices[a][2], WC6R_indices[a][3], y[c]);
        c ++;
    }
    for (a = 0; a < DWC6I; a ++) {
        WC6I_set(CuuI, WC6I_indices[a][0], WC6I_indices[a][1],
                WC6I_indices[a][2], WC6I_indices[a][3], y[c]);
        c ++;
    }

    for (a = 0; a < DWC6R; a ++) {
        WC6R_set(CddR, WC6R_indices[a][0], WC6R_indices[a][1],
                WC6R_indices[a][2], WC6R_indices[a][3], y[c]);
        c ++;
    }
    for (a = 0; a < DWC6I; a ++) {
        WC6I_set(CddI, WC6I_indices[a][0], WC6I_indices[a][1],
                WC6I_indices[a][2], WC6I_indices[a][3], y[c]);
        c ++;
    }



    for (a = 0; a < DWC7R; a ++) {
        WC7R_set(CeuR, WC7R_indices[a][0], WC7R_indices[a][1],
                WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
        c ++;

    }
    for (a = 0; a < DWC7I; a ++) {
        WC7I_set(CeuI, WC7I_indices[a][0], WC7I_indices[a][1],
                WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
        c ++;
    }


    for (a = 0; a < DWC7R; a ++) {
        WC7R_set(CedR, WC7R_indices[a][0], WC7R_indices[a][1],
                WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
        c ++;
    }
    for (a = 0; a < DWC7I; a ++) {
        WC7I_set(CedI, WC7I_indices[a][0], WC7I_indices[a][1],
                WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
        c ++;
    }


    for (a = 0; a < DWC7R; a ++) {
        WC7R_set(Cud1R, WC7R_indices[a][0], WC7R_indices[a][1],
                WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
        c ++;
    }
    for (a = 0; a < DWC7I; a ++) {
        WC7I_set(Cud1I, WC7I_indices[a][0], WC7I_indices[a][1],
                WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
        c ++;
    }


    for (a = 0; a < DWC7R; a ++) {
        WC7R_set(Cud8R, WC7R_indices[a][0], WC7R_indices[a][1],
                WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
        c ++;
    }
    for (a = 0; a < DWC7I; a ++) {
        WC7I_set(Cud8I, WC7I_indices[a][0], WC7I_indices[a][1],
                WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
        c ++;
    }

    //SMEFT class 8_LLRR
    double CleR[NG * NG * NG * NG];
    double CleI[NG * NG * NG * NG];
    double CluR[NG * NG * NG * NG];
    double CluI[NG * NG * NG * NG];
    double CldR[NG * NG * NG * NG];
    double CldI[NG * NG * NG * NG];
    double CqeR[NG * NG * NG * NG];
    double CqeI[NG * NG * NG * NG];
    double Cqu1R[NG * NG * NG * NG];
    double Cqu1I[NG * NG * NG * NG];
    double Cqu8R[NG * NG * NG * NG];
    double Cqu8I[NG * NG * NG * NG];
    double Cqd1R[NG * NG * NG * NG];
    double Cqd1I[NG * NG * NG * NG];
    double Cqd8R[NG * NG * NG * NG];
    double Cqd8I[NG * NG * NG * NG];



    {
        for (a = 0; a < DWC7R; a ++) {
            WC7R_set(CleR, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
            c ++;
        }
        for (a = 0; a < DWC7I; a ++) {
            WC7I_set(CleI, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
            c ++;
        }
        for (a = 0; a < DWC7R; a ++) {
            WC7R_set(CluR, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
            c ++;
        }
        for (a = 0; a < DWC7I; a ++) {
            WC7I_set(CluI, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
            c ++;
        }

        for (a = 0; a < DWC7R; a ++) {
            WC7R_set(CldR, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
            c ++;
        }
        for (a = 0; a < DWC7I; a ++) {
            WC7I_set(CldI, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
            c ++;
        }
        for (a = 0; a < DWC7R; a ++) {
            WC7R_set(CqeR, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
            c ++;
        }
        for (a = 0; a < DWC7I; a ++) {
            WC7I_set(CqeI, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
            c ++;
        }

        for (a = 0; a < DWC7R; a ++) {
            WC7R_set(Cqu1R, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
            c ++;
        }
        for (a = 0; a < DWC7I; a ++) {
            WC7I_set(Cqu1I, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
            c ++;
        }
        for (a = 0; a < DWC7R; a ++) {
            WC7R_set(Cqu8R, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
            c ++;
        }
        for (a = 0; a < DWC7I; a ++) {
            WC7I_set(Cqu8I, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
            c ++;
        }

        for (a = 0; a < DWC7R; a ++) {
            WC7R_set(Cqd1R, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
            c ++;
        }
        for (a = 0; a < DWC7I; a ++) {
            WC7I_set(Cqd1I, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
            c ++;
        }
        for (a = 0; a < DWC7R; a ++) {
            WC7R_set(Cqd8R, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
            c ++;
        }
        for (a = 0; a < DWC7I; a ++) {
            WC7I_set(Cqd8I, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
            c ++;
        }
    }


    //SMEFT class 8_LRRL
    double CledqR[NG * NG * NG * NG], CledqI[NG * NG * NG * NG];

    for (i = 0; i < NG; i ++) {
        for (j = 0; j < NG; j ++) {
            for (k = 0; k < NG; k ++) {
                for (l = 0; l < NG; l ++) {
                    CledqR[NG * NG * NG * i + NG * NG * j + NG * k + l] = y[c];
                    CledqI[NG * NG * NG * i + NG * NG * j + NG * k + l] = y[c + NG * NG * NG * NG];
                    c ++;
                }
            }
        }
    }
    c += NG * NG * NG*NG;

    //SMEFT class 8_LRLR

    double Cquqd1R[NG * NG * NG * NG], Cquqd1I[NG * NG * NG * NG],
            Cquqd8R[NG * NG * NG * NG], Cquqd8I[NG * NG * NG * NG],
            Clequ1R[NG * NG * NG * NG], Clequ1I[NG * NG * NG * NG],
            Clequ3R[NG * NG * NG * NG], Clequ3I[NG * NG * NG * NG];

    for (i = 0; i < NG; i ++) {
        for (j = 0; j < NG; j ++) {
            for (k = 0; k < NG; k ++) {
                for (l = 0; l < NG; l ++) {
                    a = 0;
                    WC5_set(Cquqd1R, i, j, k, l, y[c + a * NG * NG * NG * NG]);
                    a ++;
                    WC5_set(Cquqd1I, i, j, k, l, y[c + a * NG * NG * NG * NG]);
                    a ++;
                    WC5_set(Cquqd8R, i, j, k, l, y[c + a * NG * NG * NG * NG]);
                    a ++;
                    WC5_set(Cquqd8I, i, j, k, l, y[c + a * NG * NG * NG * NG]);
                    a ++;
                    WC5_set(Clequ1R, i, j, k, l, y[c + a * NG * NG * NG * NG]);
                    a ++;
                    WC5_set(Clequ1I, i, j, k, l, y[c + a * NG * NG * NG * NG]);
                    a ++;
                    WC5_set(Clequ3R, i, j, k, l, y[c + a * NG * NG * NG * NG]);
                    a ++;
                    WC5_set(Clequ3I, i, j, k, l, y[c + a * NG * NG * NG * NG]);
                    a ++;
                    c ++;
                }
            }
        }
    }

    c += NG * NG * NG * NG * (2 * N8_LRLR - 1);

    //Auxiliary quantities: 

    //Powers and products of gauge couplings
    double g12 = g1*g1;
    double g22 = g2*g2;
    double g32 = g3*g3;
    double g1g3 = g1*g3;
    double g1g2 = g1*g2;
    double g2g3 = g2*g3;

    double gammaH = 0.; //Higgs wavefunction normalization  
    double H = 0.; //scalar yukawa-trace dependent appearing in lambda RGE
    double gqR[NG][NG] = {
        {0.}
    }; //wavefunction ren. const. of q 
    double gqI[NG][NG] = {
        {0.}
    };
    double glR[NG][NG] = {
        {0.}
    }; //wavefunction ren. const. of l 
    double glI[NG][NG] = {
        {0.}
    };
    double guR[NG][NG] = {
        {0.}
    }; //wavefunction ren. const. of u
    double guI[NG][NG] = {
        {0.}
    };
    double gdR[NG][NG] = {
        {0.}
    }; //wavefunction ren. const. of d 
    double gdI[NG][NG] = {
        {0.}
    };
    double geR[NG][NG] = {
        {0.}
    }; //wavefunction ren. const. of e 
    double geI[NG][NG] = {
        {0.}
    };
    double yudyuR[NG][NG] = {
        {0.}
    }; //yu^dag yu 
    double yudyuI[NG][NG] = {
        {0.}
    };
    double yddydR[NG][NG] = {
        {0.}
    }; //yd^dag yd 
    double yddydI[NG][NG] = {
        {0.}
    };
    double yedyeR[NG][NG] = {
        {0.}
    }; //ye^dag ye 
    double yedyeI[NG][NG] = {
        {0.}
    };
    double ydyudR[NG][NG] = {
        {0.}
    }; //yd yu^dag
    double ydyudI[NG][NG] = {
        {0.}
    };
    double yuyddR[NG][NG] = {
        {0.}
    }; //yu yd^dag
    double yuyddI[NG][NG] = {
        {0.}
    };

    //3 Yukawa matrices product

    double yuyudyuR[NG][NG] = {
        {0.}
    }; //yu yu^dag yu
    double yuyudyuI[NG][NG] = {
        {0.}
    };
    double yudyuyudR[NG][NG] = {
        {0.}
    }; //yu^dag yu yu^dag
    double yudyuyudI[NG][NG] = {
        {0.}
    };
    double ydyddydR[NG][NG] = {
        {0.}
    }; //yd yd^dag yd
    double ydyddydI[NG][NG] = {
        {0.}
    };
    double yddydyddR[NG][NG] = {
        {0.}
    }; //yd^dag yd yd^dag
    double yddydyddI[NG][NG] = {
        {0.}
    };
    double yeyedyeR[NG][NG] = {
        {0.}
    }; //ye ye^dag ye
    double yeyedyeI[NG][NG] = {
        {0.}
    };
    double yedyeyedR[NG][NG] = {
        {0.}
    }; //ye^dag ye ye^dag
    double yedyeyedI[NG][NG] = {
        {0.}
    };

    double yddydyudR[NG][NG] = {
        {0.}
    }; // yd^dag yd yu^dag
    double yddydyudI[NG][NG] = {
        {0.}
    };
    double yudyuyddR[NG][NG] = {
        {0.}
    }; // yu^dag yu yd^dag
    double yudyuyddI[NG][NG] = {
        {0.}
    };


    //Real numbers defined in (A.3) https://arxiv.org/abs/1310.4838 
    double eta1 = 0.;
    double eta2 = 0.;
    double eta3 = 0.;
    double eta4 = 0.;
    double eta5 = 0.;

    //Complex matrices defined in (A.4) https://arxiv.org/abs/1310.4838 
    double xieR[NG][NG] = {
        {0.}
    };
    double xieI[NG][NG] = {
        {0.}
    };
    double xiuR[NG][NG] = {
        {0.}
    };
    double xiuI[NG][NG] = {
        {0.}
    };
    double xidR[NG][NG] = {
        {0.}
    };
    double xidI[NG][NG] = {
        {0.}
    };


    //Products of Yukawas 
    for (i = 0; i < NG; i ++) {
        for (j = 0; j < NG; j ++) {
            for (a = 0; a < NG; a ++) {
                guR[i][j] += yuR[i][a] * yudagR[a][j] - yuI[i][a] * yudagI[a][j];
                guI[i][j] += yuI[i][a] * yudagR[a][j] + yuR[i][a] * yudagI[a][j];
                yudyuR[i][j] += yudagR[i][a] * yuR[a][j] - yudagI[i][a] * yuI[a][j];
                yudyuI[i][j] += yudagI[i][a] * yuR[a][j] + yudagR[i][a] * yuI[a][j];

                gdR[i][j] += ydR[i][a] * yddagR[a][j] - ydI[i][a] * yddagI[a][j];
                gdI[i][j] += ydI[i][a] * yddagR[a][j] + ydR[i][a] * yddagI[a][j];
                yddydR[i][j] += yddagR[i][a] * ydR[a][j] - yddagI[i][a] * ydI[a][j];
                yddydI[i][j] += yddagI[i][a] * ydR[a][j] + yddagR[i][a] * ydI[a][j];

                geR[i][j] += yeR[i][a] * yedagR[a][j] - yeI[i][a] * yedagI[a][j];
                geI[i][j] += yeI[i][a] * yedagR[a][j] + yeR[i][a] * yedagI[a][j];
                yedyeR[i][j] += yedagR[i][a] * yeR[a][j] - yedagI[i][a] * yeI[a][j];
                yedyeI[i][j] += yedagI[i][a] * yeR[a][j] + yedagR[i][a] * yeI[a][j];

                ydyudR[i][j] += ydR[i][a] * yudagR[a][j] - ydI[i][a] * yudagI[a][j];
                ydyudI[i][j] += ydI[i][a] * yudagR[a][j] + ydR[i][a] * yudagI[a][j];
            }
        }
    }


    for (i = 0; i < NG; i ++) {
        for (j = 0; j < NG; j ++) {
            gqR[i][j] = 0.5 * (yudyuR[i][j] + yddydR[i][j]);
            gqI[i][j] = 0.5 * (yudyuI[i][j] + yddydI[i][j]);
            glR[i][j] = 0.5 * yedyeR[i][j];
            glI[i][j] = 0.5 * yedyeI[i][j];

            for (a = 0; a < NG; a ++) {
                yuyudyuR[i][j] += guR[i][a] * yuR[a][j] - guI[i][a] * yuI[a][j];
                yuyudyuI[i][j] += guI[i][a] * yuR[a][j] + guR[i][a] * yuI[a][j];
                ydyddydR[i][j] += gdR[i][a] * ydR[a][j] - gdI[i][a] * ydI[a][j];
                ydyddydI[i][j] += gdI[i][a] * ydR[a][j] + gdR[i][a] * ydI[a][j];
                yeyedyeR[i][j] += geR[i][a] * yeR[a][j] - geI[i][a] * yeI[a][j];
                yeyedyeI[i][j] += geI[i][a] * yeR[a][j] + geR[i][a] * yeI[a][j];

                yddydyudR[i][j] += yddydR[i][a] * yudagR[a][j] - yddydI[i][a] * yudagI[a][j];
                yddydyudI[i][j] += yddydI[i][a] * yudagR[a][j] + yddydR[i][a] * yudagI[a][j];
                yudyuyddR[i][j] += yudyuR[i][a] * yddagR[a][j] - yudyuI[i][a] * yddagI[a][j];
                yudyuyddI[i][j] += yudyuI[i][a] * yddagR[a][j] + yudyuR[i][a] * yddagI[a][j];

            }
        }
    }

    for (i = 0; i < NG; i ++) {
        for (j = 0; j < NG; j ++) {
            yudyuyudR[i][j] = yuyudyuR[j][i];
            yudyuyudI[i][j] = - yuyudyuI[j][i];
            yedyeyedR[i][j] = yeyedyeR[j][i];
            yedyeyedI[i][j] = - yeyedyeI[j][i];
            yddydyddR[i][j] = ydyddydR[j][i];
            yddydyddI[i][j] = - ydyddydI[j][i];
            yuyddR[i][j] = ydyudR[j][i];
            yuyddI[i][j] = - ydyudI[j][i];

        }
    }

    for (r = 0; r < NG; r ++) {
        for (s = 0; s < NG; s ++) {
            eta1 += + 0.5 * NC * (
                    + WC1(CdHR, r, s) * ydR[s][r] - WC1(CdHI, r, s) * ydI[s][r]
                    + WC1(CdHR, r, s) * yddagR[r][s] + WC1(CdHI, r, s) * yddagI[r][s]
                    )
                    + 0.5 * NC * (
                    + WC1(CuHR, r, s) * yuR[s][r] - WC1(CuHI, r, s) * yuI[s][r]
                    + WC1(CuHR, r, s) * yudagR[r][s] + WC1(CuHI, r, s) * yudagI[r][s]
                    )
                    + 0.5 * (
                    + WC1(CeHR, r, s) * yeR[s][r] - WC1(CeHI, r, s) * yeI[s][r]
                    + WC1(CeHR, r, s) * yedagR[r][s] + WC1(CeHI, r, s) * yedagI[r][s]
                    );

            eta2 += - 2. * NC * (
                    + WC2R(CHq3R, r, s)*(yudyuR[s][r] + yddydR[s][r])
                    - WC2I(CHq3I, r, s)*(yudyuI[s][r] + yddydI[s][r])
                    )
                    + NC * (
                    + WC1(CHudR, r, s) * ydyudR[s][r] - WC1(CHudI, r, s) * ydyudI[s][r]
                    + WC1(CHudR, s, r) * yuyddR[s][r] + WC1(CHudI, s, r) * yuyddI[s][r]
                    )
                    - 2. * (
                    + WC2R(CHl3R, r, s)*(yedyeR[s][r])
                    - WC2I(CHl3I, r, s)*(yedyeI[s][r])
                    );

            eta3 += + NC * (
                    + WC2R(CHq1R, r, s)*(yddydR[s][r] - yudyuR[s][r])
                    - WC2I(CHq1I, r, s)*(yddydI[s][r] - yudyuI[s][r])
                    )
                    + 3. * NC * (
                    + WC2R(CHq3R, r, s)*(yudyuR[s][r] + yddydR[s][r])
                    - WC2I(CHq3I, r, s)*(yudyuI[s][r] + yddydI[s][r])
                    )
                    + NC * (
                    + WC2R(CHuR, r, s) * guR[s][r] - WC2I(CHuI, r, s) * guI[s][r])
                    - NC * (
                    + WC2R(CHdR, r, s) * gdR[s][r] - WC2I(CHdI, r, s) * gdI[s][r])
                    - NC * (
                    + WC1(CHudR, r, s) * ydyudR[s][r] - WC1(CHudI, r, s) * ydyudI[s][r]
                    + WC1(CHudR, s, r) * yuyddR[s][r] + WC1(CHudI, s, r) * yuyddI[s][r]
                    )
                    + (
                    + (3. * WC2R(CHl3R, r, s) + WC2R(CHl1R, r, s)) * yedyeR[s][r]
                    -(3. * WC2I(CHl3I, r, s) + WC2I(CHl1I, r, s)) * yedyeI[s][r]
                    )
                    -(WC2R(CHeR, r, s) * geR[s][r] - WC2I(CHeI, r, s) * geI[s][r]);

            eta4 += 4. * NC * (
                    + WC2R(CHq1R, r, s)*(yddydR[s][r] - yudyuR[s][r])
                    - WC2I(CHq1I, r, s)*(yddydI[s][r] - yudyuI[s][r])
                    )
                    + 4. * NC * (
                    + WC2R(CHuR, r, s) * guR[s][r] - WC2I(CHuI, r, s) * guI[s][r])
                    - 4. * NC * (
                    + WC2R(CHdR, r, s) * gdR[s][r] - WC2I(CHdI, r, s) * gdI[s][r])
                    + 2. * NC * (
                    + WC1(CHudR, r, s) * ydyudR[s][r] - WC1(CHudI, r, s) * ydyudI[s][r]
                    + WC1(CHudR, s, r) * yuyddR[s][r] + WC1(CHudI, s, r) * yuyddI[s][r]
                    )
                    + 4. * (
                    + WC2R(CHl1R, r, s) * yedyeR[s][r] - WC2I(CHl1R, r, s) * yedyeI[s][r]
                    )
                    - 4. *
                    (WC2R(CHeR, r, s) * geR[s][r] - WC2I(CHeI, r, s) * geI[s][r]);


            eta5 +=
                    + 0.5 * NC * (
                    + (+ WC1(CdHR, r, s) * ydI[s][r] + WC1(CdHI, r, s) * ydR[s][r])
                    -(+ WC1(CdHR, r, s) * yddagI[r][s] - WC1(CdHI, r, s) * yddagR[r][s])
                    )
                    - 0.5 * NC * (
                    + (+ WC1(CuHR, r, s) * yuI[s][r] + WC1(CuHI, r, s) * yuR[s][r])
                    -(+ WC1(CuHR, r, s) * yudagI[r][s] - WC1(CuHI, r, s) * yudagR[r][s])
                    )
                    + 0.5 * (
                    + (+ WC1(CeHR, r, s) * yeI[s][r] + WC1(CeHI, r, s) * yeR[s][r])
                    -(+ WC1(CeHR, r, s) * yedagI[r][s] - WC1(CeHI, r, s) * yedagR[r][s])
                    );
        }
    }


    //gammaH
    for (i = 0; i < NG; i ++) {
        gammaH += yedyeR[i][i] + NC * (yudyuR[i][i] + yddydR[i][i]);
        for (j = 0; j < NG; j ++) {
            H += yedyeR[i][j] * yedyeR[j][i] - yedyeI[i][j] * yedyeI[j][i]
                    + NC * (yudyuR[i][j] * yudyuR[j][i] - yudyuI[i][j] * yudyuI[j][i]
                    + yddydR[i][j] * yddydR[j][i] - yddydI[i][j] * yddydI[j][i]);
        }
    }

    //xiB defined in (C.1) https://arxiv.org/abs/1312.2014
    double xiB = + FOUR_THIRDS * Yh * (CHBOX + CHD)
            + EIGHT_THIRDS * (
            + 2. * Yl * (WC2R(CHl1R, 0, 0) + WC2R(CHl1R, 1, 1) + WC2R(CHl1R, 2, 2))
            + 2. * Yq * NC * (WC2R(CHq1R, 0, 0) + WC2R(CHq1R, 1, 1) + WC2R(CHq1R, 2, 2))
            + Ye * (WC2R(CHeR, 0, 0) + WC2R(CHeR, 1, 1) + WC2R(CHeR, 2, 2))
            + Yu * NC * (WC2R(CHuR, 0, 0) + WC2R(CHuR, 1, 1) + WC2R(CHuR, 2, 2))
            + Yd * NC * (WC2R(CHdR, 0, 0) + WC2R(CHdR, 1, 1) + WC2R(CHdR, 2, 2))
            );

    //xiu,xid,xid
    for (p = 0; p < NG; p ++) {
        for (t = 0; t < NG; t ++) {
            //2 summed indices 
            for (r = 0; r < NG; r ++) {
                for (s = 0; s < NG; s ++) {
                    xieR[p][t] += 2. * (WC7R(CleR, p, r, s, t) * yeR[s][r] + WC7I(CleI, p, r, s, t) * yeI[s][r])
                            - NC * (WC5(CledqR, p, t, s, r) * ydR[s][r] + WC5(CledqI, p, t, s, r) * ydI[s][r])
                            + NC * (WC5(Clequ1R, p, t, s, r) * yuR[r][s] - WC5(Clequ1I, p, t, s, r) * yuI[r][s]);
                    xieI[p][t] += 2. * (- WC7R(CleR, p, r, s, t) * yeI[s][r] + WC7I(CleI, p, r, s, t) * yeR[s][r])
                            - NC * (- WC5(CledqR, p, t, s, r) * ydI[s][r] + WC5(CledqI, p, t, s, r) * ydR[s][r])
                            + NC * (WC5(Clequ1R, p, t, s, r) * yuI[r][s] + WC5(Clequ1I, p, t, s, r) * yuR[r][s]);

                    xidR[p][t] += 2. * ((WC7R(Cqd1R, p, r, s, t) + cF3 * WC7R(Cqd8R, p, r, s, t)) * ydR[s][r]
                            +(WC7I(Cqd1I, p, r, s, t) + cF3 * WC7I(Cqd8I, p, r, s, t)) * ydI[s][r])
                            -(
                            (NC * WC5(Cquqd1R, s, r, p, t) + 0.5 * WC5(Cquqd1R, p, r, s, t)
                            + 0.5 * cF3 * WC5(Cquqd8R, p, r, s, t)) * yuR[r][s]

                            -(NC * WC5(Cquqd1I, s, r, p, t) + 0.5 * WC5(Cquqd1I, p, r, s, t)
                            + 0.5 * cF3 * WC5(Cquqd8I, p, r, s, t)) * yuI[r][s]
                            )
                            -(WC5(CledqR, s, r, t, p) * yeR[r][s] - WC5(CledqI, s, r, t, p) * yeI[r][s]);

                    xidI[p][t] += 2. * (- (WC7R(Cqd1R, p, r, s, t) + cF3 * WC7R(Cqd8R, p, r, s, t)) * ydI[s][r]
                            +(WC7I(Cqd1I, p, r, s, t) + cF3 * WC7I(Cqd8I, p, r, s, t)) * ydR[s][r])
                            -(
                            (NC * WC5(Cquqd1R, s, r, p, t) + 0.5 * WC5(Cquqd1R, p, r, s, t)
                            + 0.5 * cF3 * WC5(Cquqd8R, p, r, s, t)) * yuI[r][s]
                            +(NC * WC5(Cquqd1I, s, r, p, t) + 0.5 * WC5(Cquqd1I, p, r, s, t)
                            + 0.5 * cF3 * WC5(Cquqd8I, p, r, s, t)) * yuR[r][s]
                            )
                            -(- WC5(CledqR, s, r, t, p) * yeI[r][s] - WC5(CledqI, s, r, t, p) * yeR[r][s]);

                    xiuR[p][t] += 2. * (
                            (WC7R(Cqu1R, p, r, s, t) + cF3 * WC7R(Cqu8R, p, r, s, t)) * yuR[s][r]
                            +(WC7I(Cqu1I, p, r, s, t) + cF3 * WC7I(Cqu8I, p, r, s, t)) * yuI[s][r]
                            )
                            -(
                            (NC * WC5(Cquqd1R, p, t, s, r) + 0.5 * WC5(Cquqd1R, s, t, p, r)
                            + 0.5 * cF3 * WC5(Cquqd8R, s, t, p, r)) * ydR[r][s]
                            -(NC * WC5(Cquqd1I, p, t, s, r) + 0.5 * WC5(Cquqd1I, s, t, p, r)
                            + 0.5 * cF3 * WC5(Cquqd8I, s, t, p, r)) * ydI[r][s]
                            )
                            +(WC5(Clequ1R, s, r, p, t) * yeR[r][s] - WC5(Clequ1I, s, r, p, t) * yeI[r][s]);
                    xiuI[p][t] += 2. * (- (WC7R(Cqu1R, p, r, s, t) + cF3 * WC7R(Cqu8R, p, r, s, t)) * yuI[s][r]
                            +(WC7I(Cqu1I, p, r, s, t) + cF3 * WC7I(Cqu8I, p, r, s, t)) * yuR[s][r])
                            -((NC * WC5(Cquqd1R, p, t, s, r) + 0.5 * WC5(Cquqd1R, s, t, p, r)
                            + 0.5 * cF3 * WC5(Cquqd8R, s, t, p, r)) * ydI[r][s]
                            +(NC * WC5(Cquqd1I, p, t, s, r) + 0.5 * WC5(Cquqd1I, s, t, p, r)
                            + 0.5 * cF3 * WC5(Cquqd8I, s, t, p, r)) * ydR[r][s])
                            +(WC5(Clequ1R, s, r, p, t) * yeI[r][s] + WC5(Clequ1I, s, r, p, t) * yeR[r][s]);
                }
            }
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
        f[c] = (- b02 * g22 //SM
                - 4. * mh2 * CHW //SMEFT
                ) * g2 * loop_factor;
        c ++;
        //g1
        f[c] = (- b01 * g12 //SM
                - 4. * mh2 * CHB //SMEFT
                ) * g1 * loop_factor;
        c ++;
        //g3
        f[c] = (- b03 * g32 //SM
                - 4. * mh2 * CHG //SMEFT
                ) * g3 * loop_factor;
        c ++;

        //lambda
        f[c] = lambda * (24. * lambda - 3. * g12 - 9. * g22 + 4. * gammaH)
                + 0.375 * g12 * g12 + 0.75 * g12 * g22
                + 1.125 * g22 * g22 - 2. * H //SM
                + mh2 * (12. * CH + (- 32. * lambda + TEN_THIRDS * g22) * CHBOX
                + (12. * lambda - 1.5 * g22 + 6. * g12 * Yh2) * CHD + 2. * (eta1 + eta2)
                + 12. * g22 * cF2 * CHW + 12. * g12 * Yh2 * CHB + 6. * g1g2 * Yh * CHWB
                + FOUR_THIRDS * g22 * (
                (WC2R(CHl3R, 0, 0) + WC2R(CHl3R, 1, 1) + WC2R(CHl3R, 2, 2)) +
                NC * (WC2R(CHq3R, 0, 0) + WC2R(CHq3R, 1, 1) + WC2R(CHq3R, 2, 2)))
                ) //SMEFT
                ;
        f[c] *= loop_factor;
        c ++;
        //mh2
        f[c] = mh2 * (12. * lambda + 2. * gammaH - 1.5 * g12 - 4.5 * g22 //SM  
                + mh2 * (2. * CHD - 4. * CHBOX) //SMEFT
                ) * loop_factor;
        c ++;
    }

    //---------RGE YUKAWA --------
    //SM beta functions can be found in
    // https://arxiv.org/abs/hep-ph/0207271
    for (i = 0; i < NG; i ++) {
        for (j = 0; j < NG; j ++) {
            //Entries without matrix products:      
            //yuR  
            f[c ] = (gammaH - (17. / 12.) * g12 - 2.25 * g22 - 8. * g32) * yuR[i][j]
                    +(3. * WC1(CuHR, j, i) - CHBOX * yuR[i][j]
                    + 0.5 * CHD * yuR[i][j]) * mh2
                    ;
            //yuI
            f[c + DF] = (gammaH - (17. / 12.) * g12 - 2.25 * g22 - 8. * g32) * yuI[i][j]
                    +(- 3. * WC1(CuHI, j, i) - CHBOX * yuI[i][j]
                    + 0.5 * CHD * yuI[i][j]) * mh2
                    ;
            //ydR     
            f[c + 2 * DF] = (gammaH - (5. / 12.) * g12 - 2.25 * g22
                    - 8. * g32) * ydR[i][j]
                    +(3. * WC1(CdHR, j, i) - CHBOX * ydR[i][j]
                    + 0.5 * CHD * ydR[i][j]) * mh2
                    ;
            //ydI  
            f[c + 3 * DF] = (gammaH - (5. / 12.) * g12 - 2.25 * g22
                    - 8. * g32) * ydI[i][j]
                    +(- 3. * WC1(CdHI, j, i) - CHBOX * ydI[i][j]
                    + 0.5 * CHD * ydI[i][j]) * mh2
                    ;
            //yeR
            f[c + 4 * DF] = (gammaH - 2.25 * g22 - 3.75 * g12) * yeR[i][j]
                    +(3. * WC1(CeHR, j, i) - CHBOX * yeR[i][j]
                    + 0.5 * CHD * yeR[i][j]) * mh2
                    ;
            //yeI
            f[c + 5 * DF] = (gammaH - 2.25 * g22 - 3.75 * g12) * yeI[i][j]
                    +(- 3. * WC1(CeHI, j, i) - CHBOX * yeI[i][j]
                    + 0.5 * CHD * yeI[i][j]) * mh2
                    ;
            //Entries with 1 matrix product (1 summed index)
            for (b = 0; b < NG; b ++) {
                //yuR
                f[c ] +=
                        + 1.5 * (yuR[i][b]*(yudyuR[b][j] - yddydR[b][j])
                        - yuI[i][b]*(yudyuI[b][j] - yddydI[b][j]))
                        +(- (yuR[i][b]*(WC2R(CHq1R, b, j) - 3. * WC2R(CHq3R, b, j))
                        - yuI[i][b]*(WC2I(CHq1I, b, j) - 3. * WC2I(CHq3I, b, j)))
                        +(WC2R(CHuR, i, b) * yuR[b][j] - WC2I(CHuI, i, b) * yuI[b][j])
                        - (WC1(CHudR, i, b) * ydR[b][j] - WC1(CHudI, i, b) * ydI[b][j])
                        ) * mh2
                        ;
                //yuI
                f[c + DF] += 1.5 * (yuI[i][b]*(yudyuR[b][j] - yddydR[b][j])
                        + yuR[i][b]*(yudyuI[b][j] - yddydI[b][j]))
                        +(- (yuI[i][b]*(WC2R(CHq1R, b, j) - 3. * WC2R(CHq3R, b, j))
                        + yuR[i][b]*(WC2I(CHq1I, b, j) - 3. * WC2I(CHq3I, b, j)))
                        +(WC2R(CHuR, i, b) * yuI[b][j] + WC2I(CHuI, i, b) * yuR[b][j])
                        - (WC1(CHudR, i, b) * ydI[b][j] + WC1(CHudI, i, b) * ydR[b][j])
                        ) * mh2
                        ;
                //ydR       
                f[c + 2 * DF] += 1.5 * (ydR[i][b]*(yddydR[b][j] - yudyuR[b][j])
                        - ydI[i][b]*(yddydI[b][j] - yudyuI[b][j]))
                        +(ydR[i][b]*(WC2R(CHq1R, b, j) + 3. * WC2R(CHq3R, b, j))
                        - ydI[i][b]*(WC2I(CHq1I, b, j) + 3. * WC2I(CHq3I, b, j))
                        -(WC2R(CHdR, i, b) * ydR[b][j] - WC2I(CHdI, i, b) * ydI[b][j])
                        -(yuR[b][j] * WC1(CHudR, b, i) + yuI[b][j] * WC1(CHudI, b, i))
                        ) * mh2
                        ;
                //ydI 
                f[c + 3 * DF] += 1.5 * (
                        ydI[i][b]*(yddydR[b][j] - yudyuR[b][j])
                        + ydR[i][b]*(yddydI[b][j] - yudyuI[b][j]))
                        +(ydI[i][b]*(WC2R(CHq1R, b, j) + 3. * WC2R(CHq3R, b, j))
                        + ydR[i][b]*(WC2I(CHq1I, b, j) + 3. * WC2I(CHq3I, b, j))
                        -(WC2R(CHdR, i, b) * ydI[b][j] + WC2I(CHdI, i, b) * ydR[b][j])
                        -(yuI[b][j] * WC1(CHudR, b, i) - yuR[b][j] * WC1(CHudI, b, i))
                        ) * mh2
                        ;
                //yeR
                f[c + 4 * DF] += 1.5 * (yeR[i][b]*(yedyeR[b][j])
                        - yeI[i][b]*(yedyeI[b][j]))
                        +(yeR[i][b]*(WC2R(CHl1R, b, j) + 3. * WC2R(CHl3R, b, j))
                        - yeI[i][b]*(WC2I(CHl1I, b, j) + 3. * WC2I(CHl3I, b, j))
                        -(WC2R(CHeR, i, b) * yeR[b][j] - WC2I(CHeI, i, b) * yeI[b][j])
                        ) * mh2
                        ;
                //yeI 
                f[c + 5 * DF] += 1.5 * (yeI[i][b]*(yedyeR[b][j])
                        + yeR[i][b]*(yedyeI[b][j]))
                        +(yeI[i][b]*(WC2R(CHl1R, b, j) + 3. * WC2R(CHl3R, b, j))
                        + yeR[i][b]*(WC2I(CHl1I, b, j) + 3. * WC2I(CHl3I, b, j))
                        -(WC2R(CHeR, i, b) * yeI[b][j] + WC2I(CHeI, i, b) * yeR[b][j])
                        ) * mh2
                        ;

                for (a = 0; a < NG; a ++) {
                    //yuR
                    f[c] += (- 2. * (
                            (WC7R(Cqu1R, a, j, i, b) + cF3 * WC7R(Cqu8R, a, j, i, b)) * yuR[b][a]
                            -(WC7I(Cqu1I, a, j, i, b) + cF3 * WC7I(Cqu8I, a, j, i, b)) * yuI[b][a])
                            -(WC5(Clequ1R, a, b, j, i) * yeR[b][a] - WC5(Clequ1I, a, b, j, i) * yeI[b][a])
                            + NC * (WC5(Cquqd1R, j, i, a, b) * ydR[b][a] - WC5(Cquqd1I, j, i, a, b) * ydI[b][a])
                            + 0.5 * (
                            (WC5(Cquqd1R, a, i, j, b) + cF3 * WC5(Cquqd8R, a, i, j, b)) * ydR[b][a]
                            -(WC5(Cquqd1I, a, i, j, b) + cF3 * WC5(Cquqd8I, a, i, j, b)) * ydI[b][a]
                            )
                            ) * mh2;
                    //yuI 
                    f[c + DF] += (- 2. * (
                            (WC7R(Cqu1R, a, j, i, b) + cF3 * WC7R(Cqu8R, a, j, i, b)) * yuI[b][a]
                            +(WC7I(Cqu1I, a, j, i, b) + cF3 * WC7I(Cqu8I, a, j, i, b)) * yuR[b][a])
                            -(- WC5(Clequ1R, a, b, j, i) * yeI[b][a] - WC5(Clequ1I, a, b, j, i) * yeR[b][a])
                            + NC * (- WC5(Cquqd1R, j, i, a, b) * ydI[b][a] - WC5(Cquqd1I, j, i, a, b) * ydR[b][a])
                            + 0.5 * (
                            - (WC5(Cquqd1R, a, i, j, b) + cF3 * WC5(Cquqd8R, a, i, j, b)) * ydI[b][a]
                            -(WC5(Cquqd1I, a, i, j, b) + cF3 * WC5(Cquqd8I, a, i, j, b)) * ydR[b][a])
                            ) * mh2;
                    //ydR
                    f[c + 2 * DF] += (- 2. * (
                            (WC7R(Cqd1R, a, j, i, b) + cF3 * WC7R(Cqd8R, a, j, i, b)) * ydR[b][a]
                            -(WC7I(Cqd1I, a, j, i, b) + cF3 * WC7I(Cqd8I, a, j, i, b)) * ydI[b][a])
                            +(WC5(CledqR, a, b, i, j) * yeR[b][a] - WC5(CledqI, a, b, i, j) * yeI[b][a])
                            + NC * (WC5(Cquqd1R, a, b, j, i) * yuR[b][a] - WC5(Cquqd1I, a, b, j, i) * yuI[b][a])
                            + 0.5 * ((WC5(Cquqd1R, j, a, b, i) + cF3 * WC5(Cquqd8R, j, a, b, i)) * yuR[a][b]
                            -(WC5(Cquqd1I, j, a, b, i) + cF3 * WC5(Cquqd8I, j, a, b, i)) * yuI[a][b])
                            ) * mh2;
                    //ydI
                    f[c + 3 * DF] += (- 2. * (
                            (WC7R(Cqd1R, a, j, i, b) + cF3 * WC7R(Cqd8R, a, j, i, b)) * ydI[b][a]
                            +(WC7I(Cqd1I, a, j, i, b) + cF3 * WC7I(Cqd8I, a, j, i, b)) * ydR[b][a])
                            +(WC5(CledqR, a, b, i, j) * yeI[b][a] + WC5(CledqI, a, b, i, j) * yeR[b][a])
                            + NC * (- WC5(Cquqd1R, a, b, j, i) * yuI[b][a] - WC5(Cquqd1I, a, b, j, i) * yuR[b][a])
                            + 0.5 * (
                            - (WC5(Cquqd1R, j, a, b, i) + cF3 * WC5(Cquqd8R, j, a, b, i)) * yuI[a][b]
                            -(WC5(Cquqd1I, j, a, b, i) + cF3 * WC5(Cquqd8I, j, a, b, i)) * yuR[a][b])
                            ) * mh2;

                    //yeR
                    f[c + 4 * DF] += (
                            - 2. * (WC7R(CleR, a, j, i, b) * yeR[b][a] - WC7I(CleI, a, j, i, b) * yeI[b][a])
                            + NC * (WC5(CledqR, j, i, a, b) * ydR[a][b] + WC5(CledqI, j, i, a, b) * ydI[a][b])
                            - NC * (WC5(Clequ1R, j, i, a, b) * yuR[b][a] - WC5(Clequ1I, j, i, a, b) * yuI[b][a])
                            ) * mh2;
                    //yeI
                    f[c + 5 * DF] += (
                            - 2. * (WC7R(CleR, a, j, i, b) * yeI[b][a] + WC7I(CleI, a, j, i, b) * yeR[b][a])
                            + NC * (WC5(CledqR, j, i, a, b) * ydI[a][b] - WC5(CledqI, j, i, a, b) * ydR[a][b])
                            - NC * (- WC5(Clequ1R, j, i, a, b) * yuI[b][a] - WC5(Clequ1I, j, i, a, b) * yuR[b][a])
                            ) * mh2;
                }
            }


            for (a = 0; a < (2 * Nyukawa); a ++) {
                f[c + a * DF] *= loop_factor;
            }
            c ++;
        }
    }
    c += (2 * Nyukawa - 1) * DF;


    //SMEFT RGE are in 
    //RGE 1 : https://arxiv.org/abs/1308.2627
    //RGE 2 : https://arxiv.org/abs/1310.4838
    //RGE 3 : https://arxiv.org/abs/1312.2014

    //----------------------RGE SMEFT class 1-3--------------------------
    {
        f[c] = (12. * cA3 - 3. * b03) * g32 * CG * loop_factor; //CG
        f[c + 1] = (12. * cA3 - 3. * b03) * g32 * CGT * loop_factor; //CGT
        f[c + 2] = (12. * cA2 - 3. * b02) * g22 * CW * loop_factor; //CW
        f[c + 3] = (12. * cA2 - 3. * b02) * g22 * CWT * loop_factor; //CWT
        c += N1;

        //CH
        f[c] = lambda * (108. * CH + 8. * eta1 + 8. * eta2) //RGEI 
                - lambda * lambda * (160. * CHBOX - 48. * CHD)
                + 6. * gammaH*CH;
        for (i = 0; i < NG; i ++) {
            for (j = 0; j < NG; j ++) {
                f[c] += - 4. * NC * (
                        2. * (WC1(CuHR, i, j) * yuyudyuR[j][i] - WC1(CuHI, i, j) * yuyudyuI[j][i])
                        + 2. * (WC1(CdHR, i, j) * ydyddydR[j][i] - WC1(CdHI, i, j) * ydyddydI[j][i]))
                        - 4. * (2. * (WC1(CeHR, i, j) * yeyedyeR[j][i] - WC1(CeHI, i, j) * yeyedyeI[j][i]));
            }
        }

        //RGE 3 
        f[c] += (- 13.5 * g22 - 4.5 * g12) * CH +
                lambda * ((40. / 3.) * g22 * CHBOX + (- 6. * g22 + 24. * g12 * Yh2) * CHD)
                - 0.75 * (4. * Yh2 * g12 + g22) * (4. * Yh2 * g12 + g22) * CHD
                + 12. * lambda * (3. * g22 * CHW + 4. * g12 * Yh2 * CHB + 2 * g1g2 * Yh * CHWB)
                -(12. * g12 * g22 * Yh2 + 9. * g22 * g22) * CHW
                - g12 * Yh2 * (48. * g12 * Yh2 + 12. * g22) * CHB
                - g1g2 * Yh * (24. * g12 * Yh2 + 6. * g22) * CHWB
                + (16. / 3.) * lambda * g22 * (
                WC2R(CHl3R, 0, 0) + WC2R(CHl3R, 1, 1) + WC2R(CHl3R, 2, 2)
                + NC * (WC2R(CHq3R, 0, 0) + WC2R(CHq3R, 1, 1) + WC2R(CHq3R, 2, 2)));

        //CHBOX
        f[c + 1] = 24. * lambda * CHBOX - (4. * g22 + (16. / 3.) * g12 * Yh2) * CHBOX
                + (20. / 3.) * g12 * Yh2 * CHD
                + 2. * g22 * (WC2R(CHl3R, 0, 0) + WC2R(CHl3R, 1, 1) + WC2R(CHl3R, 2, 2)
                + NC * (WC2R(CHq3R, 0, 0) + WC2R(CHq3R, 1, 1) + WC2R(CHq3R, 2, 2)))
                + FOUR_THIRDS * g12 * Yh * (NC * (
                Yu * (WC2R(CHuR, 0, 0) + WC2R(CHuR, 1, 1) + WC2R(CHuR, 2, 2))
                + Yd * (WC2R(CHdR, 0, 0) + WC2R(CHdR, 1, 1) + WC2R(CHdR, 2, 2))
                + 2. * Yq * (WC2R(CHq1R, 0, 0) + WC2R(CHq1R, 1, 1) + WC2R(CHq1R, 2, 2)))
                + Ye * (WC2R(CHeR, 0, 0) + WC2R(CHeR, 1, 1) + WC2R(CHeR, 2, 2))
                + 2. * Yl * (WC2R(CHl1R, 0, 0) + WC2R(CHl1R, 1, 1) + WC2R(CHl1R, 2, 2)))
                - 2. * eta3 + 4. * gammaH*CHBOX;

        //CHD
        f[c + 2] = 12. * lambda * CHD + (80. / 3.) * g12 * Yh2 * CHBOX
                + (4.5 * g22 - TEN_THIRDS * g12 * Yh2) * CHD
                + (16. / 3.) * g12 * Yh * (NC * (
                + Yu * (WC2R(CHuR, 0, 0) + WC2R(CHuR, 1, 1) + WC2R(CHuR, 2, 2))
                + Yd * (WC2R(CHdR, 0, 0) + WC2R(CHdR, 1, 1) + WC2R(CHdR, 2, 2))
                + 2. * Yq * (WC2R(CHq1R, 0, 0) + WC2R(CHq1R, 1, 1) + WC2R(CHq1R, 2, 2)))
                + Ye * (WC2R(CHeR, 0, 0) + WC2R(CHeR, 1, 1) + WC2R(CHeR, 2, 2))
                + 2. * Yl * (WC2R(CHl1R, 0, 0) + WC2R(CHl1R, 1, 1) + WC2R(CHl1R, 2, 2)))
                - 2. * eta4 + 4. * gammaH*CHD;

        f[c] *= loop_factor;
        f[c + 1] *= loop_factor;
        f[c + 2] *= loop_factor;

        c += N23;
    }
    //----------------------RGE SMEFT class 4--------------------------

    {
        //Entries without matrix products

        //CHG and CHGT
        f[c] = (12. * lambda + 2. * gammaH
                - 6. * Yh2 * g12 - 4.5 * g22
                - 2. * b03 * g32) * CHG;
        f[c + 4] = (12. * lambda + 2. * gammaH
                - 6. * Yh2 * g12 - 4.5 * g22
                - 2. * b03 * g32) * CHGT;

        //CHW and CHWT
        f[c + 2] = (12. * lambda + 2. * gammaH
                - 6. * Yh2 * g12 -
                g22 * (2.5 + 2. * b02)) * CHW
                - 15. * g22 * g2 * CW + 2. * g1g2 * Yh*CHWB;
        f[c + 6] = (12. * lambda + 2. * gammaH
                - 6. * Yh2 * g12 -
                g22 * (2.5 + 2. * b02)) * CHWT
                - 15. * g22 * g2 * CWT + 2. * g1g2 * Yh*CHWBT;

        //CHB and CHBT
        f[c + 1] = (12. * lambda + 2. * gammaH
                + (2. * Yh2 - 2. * b01) * g12
                - 4.5 * g22) * CHB
                + 6. * g1g2 * Yh * CHWB;
        f[c + 5] = (12. * lambda + 2. * gammaH
                + (2. * Yh2 - 2. * b01) * g12
                - 4.5 * g22) * CHBT
                + 6. * g1g2 * Yh * CHWBT;

        //CHWB and CHWBT
        f[c + 3] = (4. * lambda + 2. * gammaH
                - g12 * (2. * Yh2 + b01)
                + g22 * (4.5 - b02)) * CHWB
                + 6. * g1 * g22 * Yh * CW
                + 4. * g1g2 * Yh * (CHB + CHW);

        f[c + 7] = (4. * lambda + 2. * gammaH
                - g12 * (2 * Yh2 + b01)
                + g22 * (4.5 - b02)) * CHWBT
                + 6. * g1 * g22 * Yh * CWT
                + 4. * g1g2 * Yh * (CHBT + CHWT);


        //Entries with traces (2 summed indices)
        for (a = 0; a < NG; a ++) {
            for (b = 0; b < NG; b ++) {
                //CHG and CHGT
                f[c] += - 2. * g3 * (2. * (WC1(CuGR, a, b) * yuR[b][a] - WC1(CuGI, a, b) * yuI[b][a]))
                        - 2. * g3 * (2. * (WC1(CdGR, a, b) * ydR[b][a] - WC1(CdGI, a, b) * ydI[b][a]));
                f[c + 4] += - 2. * g3 * (2. * (WC1(CuGI, a, b) * yuR[b][a] + WC1(CuGR, a, b) * yuI[b][a]))
                        - 2. * g3 * (2. * (WC1(CdGI, a, b) * ydR[b][a] + WC1(CdGR, a, b) * ydI[b][a]));

                //CHW and CHWT
                f[c + 2 ] += - g2 * NC * (2. * (WC1(CuWR, a, b) * yuR[b][a] - WC1(CuWI, a, b) * yuI[b][a]))
                        - g2 * NC * (2. * (WC1(CdWR, a, b) * ydR[b][a] - WC1(CdWI, a, b) * ydI[b][a]))
                        - g2 * (2. * (WC1(CeWR, a, b) * yeR[b][a] - WC1(CeWI, a, b) * yeI[b][a]));
                f[c + 6 ] += - g2 * NC * (2. * (WC1(CuWI, a, b) * yuR[b][a] + WC1(CuWR, a, b) * yuI[b][a]))
                        - g2 * NC * (2. * (WC1(CdWI, a, b) * ydR[b][a] + WC1(CdWR, a, b) * ydI[b][a]))
                        - g2 * (2. * (WC1(CeWI, a, b) * yeR[b][a] + WC1(CeWR, a, b) * yeI[b][a]));

                //CHB and CHBT
                f[c + 1 ] += - 2. * g1 * NC * (Yu + Yq) * (2. * (WC1(CuBR, a, b) * yuR[b][a] - WC1(CuBI, a, b) * yuI[b][a]))
                        - 2. * g1 * (Yl + Ye) * (2. * (WC1(CeBR, a, b) * yeR[b][a] - WC1(CeBI, a, b) * yeI[b][a]))
                        - 2. * g1 * NC * (Yd + Yq) * (2. * (WC1(CdBR, a, b) * ydR[b][a] - WC1(CdBI, a, b) * ydI[b][a]));
                f[c + 5 ] +=
                        - 2. * g1 * NC * (Yu + Yq) *(2. * (WC1(CuBI, a, b) * yuR[b][a] + WC1(CuBR, a, b) * yuI[b][a]))
                        - 2. * g1 * (Yl + Ye) *(2. * (WC1(CeBI, a, b) * yeR[b][a] + WC1(CeBR, a, b) * yeI[b][a]))
                        - 2. * g1 * NC * (Yd + Yq) *(2. * (WC1(CdBI, a, b) * ydR[b][a] + WC1(CdBR, a, b) * ydI[b][a]));

                //CHWB and CHWBT
                f[c + 3 ] += g2 * NC * (2. * (WC1(CuBR, a, b) * yuR[b][a] - WC1(CuBI, a, b) * yuI[b][a]))
                        - g2 * NC * (2. * (WC1(CdBR, a, b) * ydR[b][a] - WC1(CdBI, a, b) * ydI[b][a]))
                        - g2 * (2. * (WC1(CeBR, a, b) * yeR[b][a] - WC1(CeBI, a, b) * yeI[b][a]))
                        + 2. * g1 * (Yq + Yu) * NC * (2. * (WC1(CuWR, a, b) * yuR[b][a] - WC1(CuWI, a, b) * yuI[b][a]))
                        - 2. * g1 * (Yq + Yd) * NC * (2. * (WC1(CdWR, a, b) * ydR[b][a] - WC1(CdWI, a, b) * ydI[b][a]))
                        - 2. * g1 * (Ye + Yl) * (2. * (WC1(CeWR, a, b) * yeR[b][a] - WC1(CeWI, a, b) * yeI[b][a]));

                f[c + 7 ] += g2 * NC * (2. * (WC1(CuBR, a, b) * yuI[b][a] + WC1(CuBI, a, b) * yuR[b][a]))
                        - g2 * NC * (2. * (WC1(CdBR, a, b) * ydI[b][a] + WC1(CdBI, a, b) * ydR[b][a]))
                        - g2 * (2. * (WC1(CeBR, a, b) * yeI[b][a] + WC1(CeBI, a, b) * yeR[b][a]))
                        + 2. * g1 * (Yq + Yu) * NC * (2. * (WC1(CuWR, a, b) * yuI[b][a] + WC1(CuWI, a, b) * yuR[b][a]))
                        - 2. * g1 * (Yq + Yd) * NC * (2. * (WC1(CdWR, a, b) * ydI[b][a] + WC1(CdWI, a, b) * ydR[b][a]))
                        - 2. * g1 * (Ye + Yl) * (2. * (WC1(CeWR, a, b) * yeI[b][a] + WC1(CeWI, a, b) * yeR[b][a]));
            }
        }


        for (i = 0; i < N4; i ++) {
            f[c + i] *= loop_factor;
        }
        c += N4;
    }

    //----------------------RGE SMEFT class 5--------------------------
    for (r = 0; r < NG; r ++) {
        for (s = 0; s < NG; s ++) {
            count = 0;
            //CuH : Re and Im
            f[c + count * DF] = lambda * (24. * WC1(CuHR, r, s)
                    - 4. * yudagR[r][s] * CHBOX + 2. * yudagR[r][s] * CHD) //RGE 1
                    + 2. * ((eta1 + eta2) * yudagR[r][s] + eta5 * yudagI[r][s])
                    + yudyuyudR[r][s]*(CHD - 6. * CHBOX) + 3. * gammaH * WC1(CuHR, r, s)
                    //RGE 2
                    + yudagR[r][s]*(TEN_THIRDS * g22 * CHBOX
                    + (- 1.5 * g22 + 6. * g12 * Yh2) * CHD)
                    -(3. * (3. * Yq2 + 3. * Yu2 - 4. * YuYq) * g12
                    + 6.75 * g22 + 6. * cF3 * g32) * WC1(CuHR, r, s)
                    + 3. * (yudagR[r][s]*(8. * g32 * cF3 * CHG + 3. * g22 * CHW
                    + 4. * (Yh2 + 2. * YuYq) * g12 * CHB - 2. * Yq * g1g2 * CHWB)
                    - yudagI[r][s]*(8. * g32 * cF3 * CHGT + 3. * g22 * CHWT
                    + 4. * (Yh2 + 2. * YuYq) * g12 * CHBT - 2. * Yq * g1g2 * CHWBT))
                    - 6. * (4. * g12 * g1 * Yh2 * Yu + 4. * g12 * g1 * Yh2 * Yq
                    - g22 * g1 * Yh) * WC1(CuBR, r, s)
                    + 3. * (4. * g12 * g2 * YhYu + 4. * g12 * g2 * YhYq - 3. * g22 * g2) * WC1(CuWR, r, s)
                    ; //RGE 3
            count ++;
            f[c + count * DF] = lambda * (24. * WC1(CuHI, r, s)
                    - 4. * yudagI[r][s] * CHBOX + 2. * yudagI[r][s] * CHD) //RGE 1
                    + 2. * ((eta1 + eta2) * yudagI[r][s] - eta5 * yudagR[r][s])
                    + yudyuyudI[r][s]*(CHD - 6. * CHBOX) + 3. * gammaH * WC1(CuHI, r, s)
                    //RGE 2
                    + yudagI[r][s]*(TEN_THIRDS * g22 * CHBOX
                    + (- 1.5 * g22 + 6. * g12 * Yh2) * CHD)
                    -(3. * (3. * Yq2 + 3. * Yu2 - 4. * YuYq) * g12
                    + 6.75 * g22 + 6. * cF3 * g32) * WC1(CuHI, r, s)
                    + 3. * (yudagI[r][s]*(8. * g32 * cF3 * CHG + 3. * g22 * CHW
                    + 4. * (Yh2 + 2. * YuYq) * g12 * CHB - 2. * Yq * g1g2 * CHWB)
                    + yudagR[r][s]*(8. * g32 * cF3 * CHGT + 3. * g22 * CHWT
                    + 4. * (Yh2 + 2. * YuYq) * g12 * CHBT - 2. * Yq * g1g2 * CHWBT))
                    - 6. * (4. * g12 * g1 * Yh2 * Yu + 4. * g12 * g1 * Yh2 * Yq
                    - g22 * g1 * Yh) * WC1(CuBI, r, s)
                    + 3. * (4. * g12 * g2 * YhYu + 4. * g12 * g2 * YhYq
                    - 3. * g22 * g2) * WC1(CuWI, r, s)
                    ; //RGE 3
            count ++;
            //CdH : Re and Im
            f[c + count * DF] = lambda * (24. * WC1(CdHR, r, s)
                    - 4. * yddagR[r][s] * CHBOX + 2. * yddagR[r][s] * CHD)
                    //RGE 1
                    + 2. * ((eta1 + eta2) * yddagR[r][s] - eta5 * yddagI[r][s])
                    + yddydyddR[r][s]*(CHD - 6. * CHBOX) + 3. * gammaH * WC1(CdHR, r, s)
                    //RGE 2
                    + yddagR[r][s]*(TEN_THIRDS * g22 * CHBOX + (- 1.5 * g22 + 6. * Yh2 * g12) * CHD)
                    -(3. * (3. * Yq2 + 3. * Yd2 - 4. * YdYq) * g12
                    + 6.75 * g22 + 6. * cF3 * g32) * WC1(CdHR, r, s)
                    + 3. * (
                    yddagR[r][s]*(8. * cF3 * g32 * CHG + 3. * g22 * CHW
                    + 4. * (Yh2 + 2. * YdYq) * g12 * CHB + 2. * Yq * g1g2 * CHWB)
                    - yddagI[r][s]*(8. * cF3 * g32 * CHGT + 3. * g22 * CHWT
                    + 4. * (Yh2 + 2. * YdYq) * g12 * CHBT + 2. * Yq * g1g2 * CHWBT)
                    )
                    - 6. * (4. * g12 * g1 * Yh2 * Yd + 4. * g12 * g1 * Yh2 * Yq
                    + g22 * g1 * Yh) * WC1(CdBR, r, s)
                    - 3. * (4. * g12 * g2 * YhYd + 4. * g12 * g2 * YhYq
                    + 3. * g22 * g2) * WC1(CdWR, r, s)
                    ; //RGE 3
            count ++;

            f[c + count * DF] = lambda * (24. * WC1(CdHI, r, s)
                    - 4. * yddagI[r][s] * CHBOX + 2. * yddagI[r][s] * CHD)
                    //RGE 1
                    + 2. * ((eta1 + eta2) * yddagI[r][s] + eta5 * yddagR[r][s])
                    + yddydyddI[r][s]*(CHD - 6. * CHBOX) + 3. * gammaH * WC1(CdHI, r, s)
                    //RGE 2
                    + yddagI[r][s]*(TEN_THIRDS * g22 * CHBOX + (- 1.5 * g22 + 6. * Yh2 * g12) * CHD)
                    -(3. * (3. * Yq2 + 3. * Yd2 - 4. * YdYq) * g12
                    + 6.75 * g22 + 6. * cF3 * g32) * WC1(CdHI, r, s)
                    + 3. * (
                    yddagI[r][s]*(8. * cF3 * g32 * CHG + 3. * g22 * CHW
                    + 4. * (Yh2 + 2. * YdYq) * g12 * CHB + 2. * Yq * g1g2 * CHWB)
                    + yddagR[r][s]*(8. * cF3 * g32 * CHGT + 3. * g22 * CHWT
                    + 4. * (Yh2 + 2. * YdYq) * g12 * CHBT + 2. * Yq * g1g2 * CHWBT)
                    )
                    - 6. * (4. * g12 * g1 * Yh2 * Yd + 4. * g12 * g1 * Yh2 * Yq
                    + g22 * g1 * Yh) * WC1(CdBI, r, s)
                    - 3. * (4. * g12 * g2 * YhYd + 4. * g12 * g2 * YhYq
                    + 3. * g22 * g2) * WC1(CdWI, r, s)
                    ; //RGE 3
            count ++;
            //CeH : Re and Im
            f[c + count * DF] = lambda * (24. * WC1(CeHR, r, s)
                    - 4. * yedagR[r][s] * CHBOX + 2. * yedagR[r][s] * CHD) //RGE 1
                    + 2. * ((eta1 + eta2) * yedagR[r][s] - eta5 * yedagI[r][s])
                    + yedyeyedR[r][s]*(CHD - 6. * CHBOX) + 3. * gammaH * WC1(CeHR, r, s)
                    //RGE 2
                    + yedagR[r][s]*(TEN_THIRDS * g22 * CHBOX + (- 1.5 * g22 + 6. * g12 * Yh2) * CHD)
                    -(3. * (3. * Yl2 + 3. * Ye2 - 4. * YeYl) * g12
                    + 6.75 * g22) * WC1(CeHR, r, s)
                    + 3. * (yedagR[r][s]*(3. * g22 * CHW
                    + 4. * (Yh2 + 2. * YeYl) * g12 * CHB
                    + 2. * g1g2 * Yl * CHWB)
                    - yedagI[r][s]*(3. * g22 * CHWT
                    + 4. * (Yh2 + 2. * YeYl) * g12 * CHBT
                    + 2. * g1g2 * Yl * CHWBT))
                    - 6. * (4. * g12 * g1 * Yh2 * Ye + 4. * g12 * g1 * Yh2 * Yl
                    + g22 * g1 * Yh) * WC1(CeBR, r, s)
                    - 3. * (4. * g12 * g2 * YhYe + 4. * g12 * g2 * YhYl
                    + 3. * g22 * g2) * WC1(CeWR, r, s);
            //RGE 3
            count ++;
            f[c + count * DF] = lambda * (24. * WC1(CeHI, r, s)
                    - 4. * yedagI[r][s] * CHBOX + 2. * yedagI[r][s] * CHD) //RGE 1
                    + 2. * ((eta1 + eta2) * yedagI[r][s] + eta5 * yedagR[r][s])
                    + yedyeyedI[r][s]*(CHD - 6. * CHBOX) + 3. * gammaH * WC1(CeHI, r, s)
                    //RGE 2
                    + yedagI[r][s]*(TEN_THIRDS * g22 * CHBOX + (- 1.5 * g22 + 6. * g12 * Yh2) * CHD)
                    -(3. * (3. * Yl2 + 3. * Ye2 - 4. * YeYl) * g12
                    + 6.75 * g22) * WC1(CeHI, r, s)
                    + 3. * (yedagI[r][s]*(3. * g22 * CHW
                    + 4. * (Yh2 + 2. * YeYl) * g12 * CHB
                    + 2. * g1g2 * Yl * CHWB)
                    + yedagR[r][s]*(3. * g22 * CHWT
                    + 4. * (Yh2 + 2. * YeYl) * g12 * CHBT
                    + 2. * g1g2 * Yl * CHWBT))
                    - 6. * (4. * g12 * g1 * Yh2 * Ye + 4. * g12 * g1 * Yh2 * Yl
                    + g22 * g1 * Yh) * WC1(CeBI, r, s)
                    - 3. * (4. * g12 * g2 * YhYe + 4. * g12 * g2 * YhYl
                    + 3. * g22 * g2) * WC1(CeWI, r, s)
                    ;
            //RGE 3
            count ++;
            //Entries with 1 summed index
            for (t = 0; t < NG; t ++) {
                count = 0;
                //CuH : Re and Im
                f[c + count * DF] += lambda * (- 4. * (WC2R(CHq1R, r, t) * yudagR[t][s]
                        - WC2I(CHq1I, r, t) * yudagI[t][s])
                        + 12. * (WC2R(CHq3R, r, t) * yudagR[t][s] - WC2I(CHq3I, r, t) * yudagI[t][s])
                        + 4. * (yudagR[r][t] * WC2R(CHuR, t, s) - yudagI[r][t] * WC2I(CHuI, t, s))
                        - 4. * (yddagR[r][t] * WC1(CHudR, s, t) + yddagI[r][t] * WC1(CHudI, s, t))) //RGE1
                        - 2. * (WC2R(CHq1R, r, t) * yudyuyudR[t][s] - WC2I(CHq1I, r, t) * yudyuyudI[t][s])
                        + 6. * (WC2R(CHq3R, r, t) * yddydyudR[t][s] - WC2I(CHq3I, r, t) * yddydyudI[t][s])
                        + 2. * (yudyuyudR[r][t] * WC2R(CHuR, t, s) - yudyuyudI[r][t] * WC2I(CHuI, t, s))
                        - 2. * (yddydyddR[r][t] * WC1(CHudR, s, t) + yddydyddI[r][t] * WC1(CHudI, s, t))
                        + 4. * (WC1(CuHR, r, t) * guR[t][s] - WC1(CuHI, r, t) * guI[t][s])
                        + 5. * (yudyuR[r][t] * WC1(CuHR, t, s) - yudyuI[r][t] * WC1(CuHI, t, s))
                        -(WC1(CdHR, r, t) * ydyudR[t][s] - WC1(CdHI, r, t) * ydyudI[t][s])
                        - 2. * (yddydR[r][t] * WC1(CuHR, t, s) - yddydI[r][t] * WC1(CuHI, t, s))
                        + (gqR[r][t] * WC1(CuHR, t, s) - gqI[r][t] * WC1(CuHI, t, s))
                        +(WC1(CuHR, r, t) * guR[t][s] - WC1(CuHI, r, t) * guI[t][s])
                        //RGE 2
                        - 12. * g2 * (yddydR[r][t] * WC1(CuWR, t, s) - yddydI[r][t] * WC1(CuWI, t, s))
                        - 6. * g2 * (WC1(CdWR, r, t) * ydyudR[t][s] - WC1(CdWI, r, t) * ydyudI[t][s])
                        - 3. * ((4. * g3 * cF3 * WC1(CuGR, r, t) + g2 * WC1(CuWR, r, t)
                        +(3. * Yu + Yd) * g1 * WC1(CuBR, r, t)) * guR[t][s]
                        -(4. * g3 * cF3 * WC1(CuGI, r, t) + g2 * WC1(CuWI, r, t)
                        +(3. * Yu + Yd) * WC1(CuBI, r, t)) * guI[t][s])
                        - 3. * (yudyuR[r][t]*(4. * cF3 * g3 * WC1(CuGR, t, s)
                        - g2 * WC1(CuWR, t, s) + 2. * (Yq + Yu) * g1 * WC1(CuBR, t, s))
                        - yudyuI[r][t]*(4. * cF3 * g3 * WC1(CuGI, t, s)
                        - g2 * WC1(CuWI, t, s) + 2. * (Yq + Yu) * g1 * WC1(CuBI, t, s)))
                        -(3. * g22 - 12. * g12 * YhYq)*(
                        yudagR[r][t] * WC2R(CHuR, t, s) - yudagI[r][t] * WC2I(CHuI, t, s))
                        + 3. * g22 * (yddagR[r][t] * WC1(CHudR, s, t)
                        + yddagI[r][t] * WC1(CHudI, s, t))
                        + 12. * g12 * YhYu * (WC2R(CHq1R, r, t) * yudagR[t][s]
                        - WC2I(CHq1I, r, t) * yudagI[t][s])
                        - 12. * g12 * YhYu * (WC2R(CHq3R, r, t) * yudagR[t][s]
                        - WC2I(CHq3I, r, t) * yudagI[t][s])
                        + FOUR_THIRDS * g22 * yudagR[r][s]*(
                        WC2R(CHl3R, t, t) + NC * WC2R(CHq3R, t, t))
                        ; //RGE 3



                count ++;
                f[c + count * DF] += lambda * (- 4. * (WC2R(CHq1R, r, t) * yudagI[t][s]
                        + WC2I(CHq1I, r, t) * yudagR[t][s])
                        + 12. * (WC2R(CHq3R, r, t) * yudagI[t][s] + WC2I(CHq3I, r, t) * yudagR[t][s])
                        + 4. * (yudagI[r][t] * WC2R(CHuR, t, s) + yudagR[r][t] * WC2I(CHuI, t, s))
                        - 4. * (yddagI[r][t] * WC1(CHudR, s, t) - yddagR[r][t] * WC1(CHudI, s, t))) //RGE1
                        - 2. * (WC2R(CHq1R, r, t) * yudyuyudI[t][s] + WC2I(CHq1I, r, t) * yudyuyudR[t][s])
                        + 6. * (WC2R(CHq3R, r, t) * yddydyudI[t][s]+(WC2I(CHq3I, r, t) * yddydyudR[t][s]))
                        + 2. * (yudyuyudI[r][t] * WC2R(CHuR, t, s) + yudyuyudR[r][t] * WC2I(CHuI, t, s))
                        - 2. * (yddydyddI[r][t] * WC1(CHudR, s, t) - yddydyddR[r][t] * WC1(CHudI, s, t))
                        + 4. * (WC1(CuHR, r, t) * guI[t][s] + WC1(CuHI, r, t) * guR[t][s])
                        + 5. * (yudyuI[r][t] * WC1(CuHR, t, s) + yudyuR[r][t] * WC1(CuHI, t, s))
                        -(WC1(CdHR, r, t) * ydyudI[t][s] + WC1(CdHI, r, t) * ydyudR[t][s])
                        - 2. * (yddydI[r][t] * WC1(CuHR, t, s) + yddydR[r][t] * WC1(CuHI, t, s))
                        + (gqI[r][t] * WC1(CuHR, t, s) + gqR[r][t] * WC1(CuHI, t, s))
                        +(WC1(CuHR, r, t) * guI[t][s] + WC1(CuHI, r, t) * guR[t][s])
                        //RGE 2
                        - 12. * g2 * (yddydI[r][t] * WC1(CuWR, t, s) + yddydR[r][t] * WC1(CuWI, t, s))
                        - 6. * g2 * (WC1(CdWR, r, t) * ydyudI[t][s] + WC1(CdWI, r, t) * ydyudR[t][s])
                        - 3. * ((4. * g3 * cF3 * WC1(CuGR, r, t) + g2 * WC1(CuWR, r, t)
                        +(3. * Yu + Yd) * WC1(CuBR, r, t)) * guI[t][s]
                        +(4. * g3 * cF3 * WC1(CuGI, r, t) + g2 * WC1(CuWI, r, t)
                        +(3. * Yu + Yd) * g1 * WC1(CuBI, r, t)) * guR[t][s])
                        - 3. * (yudyuI[r][t]*(4. * cF3 * g3 * WC1(CuGR, t, s)
                        - g2 * WC1(CuWR, t, s) + 2. * (Yq + Yu) * g1 * WC1(CuBR, t, s))
                        + yudyuR[r][t]*(4. * cF3 * g3 * WC1(CuGI, t, s)
                        - g2 * WC1(CuWI, t, s) + 2. * (Yq + Yu) * g1 * WC1(CuBI, t, s)))
                        -(3. * g22 - 12. * g12 * YhYq)*(
                        yudagI[r][t] * WC2R(CHuR, t, s) + yudagR[r][t] * WC2I(CHuI, t, s))
                        + 3. * g22 * (yddagI[r][t] * WC1(CHudR, s, t)
                        - yddagR[r][t] * WC1(CHudI, s, t))
                        + 12. * g12 * YhYu * (WC2R(CHq1R, r, t) * yudagI[t][s]
                        + WC2I(CHq1I, r, t) * yudagR[t][s])
                        - 12. * g12 * YhYu * (WC2R(CHq3R, r, t) * yudagI[t][s]
                        + WC2I(CHq3I, r, t) * yudagR[t][s])
                        + FOUR_THIRDS * g22 * yudagI[r][s]*(
                        WC2R(CHl3R, t, t) + NC * WC2R(CHq3R, t, t))
                        ; //RGE 3
                count ++;

                //CdH : Re and Im
                f[c + count * DF] += lambda * (
                        + 4. * (WC2R(CHq1R, r, t) * yddagR[t][s] - WC2I(CHq1I, r, t) * yddagI[t][s])
                        + 12. * (WC2R(CHq3R, r, t) * yddagR[t][s] - WC2I(CHq3I, r, t) * yddagI[t][s])
                        - 4. * (yddagR[r][t] * WC2R(CHdR, t, s) - yddagI[r][t] * WC2I(CHdI, t, s))
                        - 4. * (yudagR[r][t] * WC1(CHudR, t, s) - yudagI[r][t] * WC1(CHudI, t, s)))
                        //RGE 1
                        + 2. * (WC2R(CHq1R, r, t) * yddydyddR[t][s] - WC2I(CHq1I, r, t) * yddydyddI[t][s])
                        + 6. * (WC2R(CHq3R, r, t) * yudyuyddR[t][s] - WC2I(CHq3I, r, t) * yudyuyddI[t][s])
                        - 2. * (yddydyddR[r][t] * WC2R(CHdR, t, s) - yddydyddI[r][t] * WC2I(CHdI, t, s))
                        - 2. * (yudyuyudR[r][t] * WC1(CHudR, t, s) - yudyuyudI[r][t] * WC1(CHudI, t, s))
                        + 4. * (WC1(CdHR, r, t) * gdR[t][s] - WC1(CdHI, r, t) * gdI[t][s])
                        + 5. * (yddydR[r][t] * WC1(CdHR, t, s) - yddydI[r][t] * WC1(CdHI, t, s))
                        -(WC1(CuHR, r, t) * yuyddR[t][s] - WC1(CuHI, r, t) * yuyddI[t][s])
                        - 2. * (yudyuR[r][t] * WC1(CdHR, t, s) - yudyuI[r][t] * WC1(CdHI, t, s))
                        +(gqR[r][t] * WC1(CdHR, t, s) - gqI[r][t] * WC1(CdHI, t, s))
                        +(WC1(CdHR, r, t) * gdR[t][s] - WC1(CdHI, r, t) * gdI[t][s])
                        //RGE 2
                        - 12. * g2 * (yudyuR[r][t] * WC1(CdWR, t, s) - yudyuI[r][t] * WC1(CdWI, t, s))
                        - 6. * g2 * (WC1(CuWR, r, t) * yuyddR[t][s] - WC1(CuWI, r, t) * yuyddI[t][s])
                        - 3. * (
                        (4. * cF3 * g3 * WC1(CdGR, r, t)
                        + g2 * WC1(CdWR, r, t)+(3. * Yd + Yu) * g1 * WC1(CdBR, r, t)) * gdR[t][s]
                        -(4. * cF3 * g3 * WC1(CdGI, r, t)
                        + g2 * WC1(CdWI, r, t)+(3. * Yd + Yu) * g1 * WC1(CdBI, r, t)) * gdI[t][s]
                        )
                        - 3. * (
                        yddydR[r][t]*(4. * cF3 * g3 * WC1(CdGR, t, s)
                        - g2 * WC1(CdWR, t, s) + 2. * (Yq + Yd) * g1 * WC1(CdBR, t, s))
                        - yddydI[r][t]*(4. * cF3 * g3 * WC1(CdGI, t, s)
                        - g2 * WC1(CdWI, t, s) + 2. * (Yq + Yd) * g1 * WC1(CdBI, t, s))
                        )
                        +(3. * g22 + 12. * g12 * YhYq)*
                        (yddagR[r][t] * WC2R(CHdR, t, s) - yddagI[r][t] * WC2I(CHdI, t, s))
                        + 3. * g22 * (yudagR[r][t] * WC1(CHudR, t, s)
                        - yudagI[r][t] * WC1(CHudI, t, s))
                        + 12. * g12 * YhYd * (WC2R(CHq1R, r, t) * yddagR[t][s]
                        - WC2I(CHq1I, r, t) * yddagI[t][s])
                        + 12. * g12 * YhYd * (WC2R(CHq3R, r, t) * yddagR[t][s]
                        - WC2I(CHq3I, r, t) * yddagI[t][s])
                        + FOUR_THIRDS * g22 * yddagR[r][s]*
                        (WC2R(CHl3R, t, t) + NC * WC2R(CHq3R, t, t))
                        ; //RGE 3


                count ++;
                f[c + count * DF] += lambda * (
                        + 4. * (WC2R(CHq1R, r, t) * yddagI[t][s] + WC2I(CHq1I, r, t) * yddagR[t][s])
                        + 12. * (WC2R(CHq3R, r, t) * yddagI[t][s] + WC2I(CHq3I, r, t) * yddagR[t][s])
                        - 4. * (yddagI[r][t] * WC2R(CHdR, t, s) + yddagR[r][t] * WC2I(CHdI, t, s))
                        - 4. * (yudagI[r][t] * WC1(CHudR, t, s) + yudagR[r][t] * WC1(CHudI, t, s)))
                        //RGE 1
                        + 2. * (WC2R(CHq1R, r, t) * yddydyddI[t][s] + WC2I(CHq1I, r, t) * yddydyddR[t][s])
                        + 6. * (WC2R(CHq3R, r, t) * yudyuyddI[t][s] + WC2I(CHq3I, r, t) * yudyuyddR[t][s])
                        - 2. * (yddydyddI[r][t] * WC2R(CHdR, t, s) + yddydyddR[r][t] * WC2I(CHdI, t, s))
                        - 2. * (yudyuyudI[r][t] * WC1(CHudR, t, s) + yudyuyudR[r][t] * WC1(CHudI, t, s))
                        + 4. * (WC1(CdHR, r, t) * gdI[t][s] + WC1(CdHI, r, t) * gdR[t][s])
                        + 5. * (yddydI[r][t] * WC1(CdHR, t, s) + yddydR[r][t] * WC1(CdHI, t, s))
                        -(WC1(CuHR, r, t) * yuyddI[t][s] + WC1(CuHI, r, t) * yuyddR[t][s])
                        - 2. * (yudyuI[r][t] * WC1(CdHR, t, s) + yudyuR[r][t] * WC1(CdHI, t, s))
                        +(gqI[r][t] * WC1(CdHR, t, s) + gqR[r][t] * WC1(CdHI, t, s))
                        +(WC1(CdHR, r, t) * gdI[t][s] + WC1(CdHI, r, t) * gdR[t][s])
                        //RGE 2
                        - 12. * g2 * (yudyuI[r][t] * WC1(CdWR, t, s) + yudyuR[r][t] * WC1(CdWI, t, s))
                        - 6. * g2 * (WC1(CuWR, r, t) * yuyddI[t][s] + WC1(CuWI, r, t) * yuyddR[t][s])
                        - 3. * (
                        (4. * cF3 * g3 * WC1(CdGR, r, t)
                        + g2 * WC1(CdWR, r, t)+(3. * Yd + Yu) * g1 * WC1(CdBR, r, t)) * gdI[t][s]
                        +(4. * cF3 * g3 * WC1(CdGI, r, t)
                        + g2 * WC1(CdWI, r, t)+(3. * Yd + Yu) * g1 * WC1(CdBI, r, t)) * gdR[t][s]
                        )
                        - 3. * (
                        yddydI[r][t]*(4. * cF3 * g3 * WC1(CdGR, t, s)
                        - g2 * WC1(CdWR, t, s) + 2. * (Yq + Yd) * g1 * WC1(CdBR, t, s))
                        + yddydR[r][t]*(4. * cF3 * g3 * WC1(CdGI, t, s)
                        - g2 * WC1(CdWI, t, s) + 2. * (Yq + Yd) * g1 * WC1(CdBI, t, s))
                        )
                        +(3. * g22 + 12. * g12 * YhYq)*
                        (yddagI[r][t] * WC2R(CHdR, t, s) + yddagR[r][t] * WC2I(CHdI, t, s))
                        + 3. * g22 * (yudagI[r][t] * WC1(CHudR, t, s)
                        + yudagR[r][t] * WC1(CHudI, t, s))
                        + 12. * g12 * YhYd * (WC2R(CHq1R, r, t) * yddagI[t][s]
                        + WC2I(CHq1I, r, t) * yddagR[t][s])
                        + 12. * g12 * YhYd * (WC2R(CHq3R, r, t) * yddagI[t][s]
                        + WC2I(CHq3I, r, t) * yddagR[t][s])
                        + FOUR_THIRDS * g22 * yddagI[r][s]*
                        (WC2R(CHl3R, t, t) + NC * WC2R(CHq3R, t, t))
                        ; //RGE 3

                count ++;
                //CeH : Re and Im 
                f[c + count * DF] += lambda * (4. * (WC2R(CHl1R, r, t) * yedagR[t][s] - WC2I(CHl1I, r, t) * yedagI[t][s])
                        + 12. * (WC2R(CHl3R, r, t) * yedagR[t][s] - WC2I(CHl3I, r, t) * yedagI[t][s])
                        - 4. * (yedagR[r][t] * WC2R(CHeR, t, s) - yedagI[r][t] * WC2I(CHeI, t, s))) //RGE1
                        + 2. * (WC2R(CHl1R, r, t) * yedyeyedR[t][s] - WC2I(CHl1I, r, t) * yedyeyedI[t][s])
                        - 2. * (yedyeyedR[r][t] * WC2R(CHeR, t, s) - yedyeyedI[r][t] * WC2I(CHeI, t, s))
                        + 5. * (yedyeR[r][t] * WC1(CeHR, t, s) - yedyeI[r][t] * WC1(CeHI, t, s))
                        + 4. * (WC1(CeHR, r, t) * geR[t][s] - WC1(CeHI, r, t) * geI[t][s])
                        +(glR[r][t] * WC1(CeHR, t, s) - glI[r][t] * WC1(CeHI, t, s))
                        +(WC1(CeHR, r, t) * geR[t][s] - WC1(CeHI, r, t) * geI[t][s])
                        //RGE 2
                        - 3. * ((3. * g1 * Ye * WC1(CeBR, r, t) + g2 * WC1(CeWR, r, t)) * geR[t][s]
                        -(3. * g1 * Ye * WC1(CeBI, r, t) + g2 * WC1(CeWI, r, t)) * geI[t][s])
                        - 3. * (yedyeR[r][t]*(2. * g1 * (Yl + Ye) * WC1(CeBR, t, s) - g2 * WC1(CeWR, t, s))
                        - yedyeI[r][t]*(2. * g1 * (Yl + Ye) * WC1(CeBI, t, s) - g2 * WC1(CeWI, t, s)))
                        +(3. * g22 + 12. * g12 * YhYl)*
                        (yedagR[r][t] * WC2R(CHeR, t, s) - yedagI[r][t] * WC2I(CHeI, t, s))
                        + 12. * g12 * YhYe * (WC2R(CHl1R, r, t) * yedagR[t][s]
                        - WC2I(CHl1I, r, t) * yedagI[t][s])
                        + 12. * g12 * YhYe * (WC2R(CHl3R, r, t) * yedagR[t][s]
                        - WC2I(CHl3I, r, t) * yedagI[t][s])
                        + FOUR_THIRDS * g22 * yedagR[r][s]*
                        (WC2R(CHl3R, t, t) + NC * WC2R(CHq3R, t, t))
                        ; //RGE 3

                count ++;
                f[c + count * DF] += lambda * (4. * (WC2R(CHl1R, r, t) * yedagI[t][s] + WC2I(CHl1I, r, t) * yedagR[t][s])
                        + 12. * (WC2R(CHl3R, r, t) * yedagI[t][s] + WC2I(CHl3I, r, t) * yedagR[t][s])
                        - 4. * (yedagI[r][t] * WC2R(CHeR, t, s) + yedagR[r][t] * WC2I(CHeI, t, s))) //RGE1
                        + 2. * (WC2R(CHl1R, r, t) * yedyeyedI[t][s] + WC2I(CHl1I, r, t) * yedyeyedR[t][s])
                        - 2. * (yedyeyedI[r][t] * WC2R(CHeR, t, s) + yedyeyedR[r][t] * WC2I(CHeI, t, s))
                        + 5. * (yedyeI[r][t] * WC1(CeHR, t, s) + yedyeR[r][t] * WC1(CeHI, t, s))
                        + 4. * (WC1(CeHR, r, t) * geI[t][s] + WC1(CeHI, r, t) * geR[t][s])
                        +(glI[r][t] * WC1(CeHR, t, s) + glR[r][t] * WC1(CeHI, t, s))
                        +(WC1(CeHR, r, t) * geI[t][s] + WC1(CeHI, r, t) * geR[t][s])
                        //RGE 2
                        - 3. * ((3. * g1 * Ye * WC1(CeBR, r, t) + g2 * WC1(CeWR, r, t)) * geI[t][s]
                        +(3. * g1 * Ye * WC1(CeBI, r, t) + g2 * WC1(CeWI, r, t)) * geR[t][s])
                        - 3. * (yedyeI[r][t]*(2. * g1 * (Yl + Ye) * WC1(CeBR, t, s) - g2 * WC1(CeWR, t, s))
                        + yedyeR[r][t]*(2. * g1 * (Yl + Ye) * WC1(CeBI, t, s) - g2 * WC1(CeWI, t, s)))
                        +(3. * g22 + 12. * g12 * YhYl)*
                        (yedagI[r][t] * WC2R(CHeR, t, s) + yedagR[r][t] * WC2I(CHeI, t, s))
                        + 12. * g12 * YhYe * (WC2R(CHl1R, r, t) * yedagI[t][s]
                        + WC2I(CHl1I, r, t) * yedagR[t][s])
                        + 12. * g12 * YhYe * (WC2R(CHl3R, r, t) * yedagI[t][s]
                        + WC2I(CHl3I, r, t) * yedagR[t][s])
                        + FOUR_THIRDS * g22 * yedagI[r][s]*
                        (WC2R(CHl3R, t, t) + NC * WC2R(CHq3R, t, t))
                        ; //RGE 3



                count ++;

                //entries with 2 summed indices 
                for (p = 0; p < NG; p ++) {
                    count = 0;
                    //CuH : Re and Im 
                    f[c + count * DF] += lambda * (- 8. * (WC7R(Cqu1R, r, p, t, s) * yudagR[p][t] - WC7I(Cqu1I, r, p, t, s) * yudagI[p][t])
                            - 8. * cF3 * (WC7R(Cqu8R, r, p, t, s) * yudagR[p][t] - WC7I(Cqu8I, r, p, t, s) * yudagI[p][t])
                            - 4. * (WC5(Clequ1R, p, t, r, s) * yeR[t][p] - WC5(Clequ1I, p, t, r, s) * yeI[t][p])
                            + 4. * NC * (WC5(Cquqd1R, r, s, p, t) * ydR[t][p] - WC5(Cquqd1I, r, s, p, t) * ydI[t][p])
                            + 2. * (WC5(Cquqd1R, p, s, r, t) * ydR[t][p] - WC5(Cquqd1I, p, s, r, t) * ydI[t][p])
                            + 2. * cF3 * (WC5(Cquqd8R, p, s, r, t) * ydR[t][p] - WC5(Cquqd8I, p, s, r, t) * ydI[t][p])
                            ) //RGE 1
                            + 8. * ((WC7R(Cqu1R, r, p, t, s) + cF3 * WC7R(Cqu8R, r, p, t, s)) * yudyuyudR[p][t]
                            -(WC7I(Cqu1I, r, p, t, s) + cF3 * WC7I(Cqu8I, r, p, t, s)) * yudyuyudI[p][t])
                            - 2. * ((2. * NC * WC5(Cquqd1R, r, s, t, p) + WC5(Cquqd1R, t, s, r, p)
                            + cF3 * WC5(Cquqd8R, t, s, r, p)) * ydyddydR[p][t]
                            -(2. * NC * WC5(Cquqd1I, r, s, t, p) + WC5(Cquqd1I, t, s, r, p)
                            + cF3 * WC5(Cquqd8I, t, s, r, p)) * ydyddydI[p][t])
                            + 4. * (WC5(Clequ1R, t, p, r, s) * yeyedyeR[p][t] - WC5(Clequ1I, t, p, r, s) * yeyedyeI[p][t])
                            - 2. * ((yddagR[r][t] * WC1(CdHR, p, t) + yddagI[r][t] * WC1(CdHI, p, t)) * yudagR[p][s]
                            -(yddagI[r][t] * WC1(CdHR, p, t) - yddagR[r][t] * WC1(CdHI, p, t)) * yudagI[p][s]);
                    ; //RGE 2
                    count ++;
                    f[c + count * DF] += lambda * (- 8. * (WC7R(Cqu1R, r, p, t, s) * yudagI[p][t] + WC7I(Cqu1I, r, p, t, s) * yudagR[p][t])
                            - 8. * cF3 * (WC7R(Cqu8R, r, p, t, s) * yudagI[p][t] + WC7I(Cqu8I, r, p, t, s) * yudagR[p][t])
                            - 4. * (WC5(Clequ1R, p, t, r, s) * yeI[t][p] + WC5(Clequ1I, p, t, r, s) * yeR[t][p])
                            + 4. * NC * (WC5(Cquqd1R, r, s, p, t) * ydI[t][p] + WC5(Cquqd1I, r, s, p, t) * ydR[t][p])
                            + 2. * (WC5(Cquqd1R, p, s, r, t) * ydI[t][p] + WC5(Cquqd1I, p, s, r, t) * ydR[t][p])
                            + 2. * cF3 * (WC5(Cquqd8R, p, s, r, t) * ydI[t][p] + WC5(Cquqd8I, p, s, r, t) * ydR[t][p])
                            ) //RGE 1
                            + 8. * ((WC7R(Cqu1R, r, p, t, s) + cF3 * WC7R(Cqu8R, r, p, t, s)) * yudyuyudI[p][t]
                            +(WC7I(Cqu1I, r, p, t, s) + cF3 * WC7I(Cqu8I, r, p, t, s)) * yudyuyudR[p][t])
                            - 2. * ((2. * NC * WC5(Cquqd1R, r, s, t, p) + WC5(Cquqd1R, t, s, r, p)
                            + cF3 * WC5(Cquqd8R, t, s, r, p)) * ydyddydI[p][t]
                            +(2. * NC * WC5(Cquqd1I, r, s, t, p) + WC5(Cquqd1I, t, s, r, p)
                            + cF3 * WC5(Cquqd8I, t, s, r, p)) * ydyddydR[p][t])
                            + 4. * (WC5(Clequ1R, t, p, r, s) * yeyedyeI[p][t] + WC5(Clequ1I, t, p, r, s) * yeyedyeR[p][t])
                            - 2. * ((yddagR[r][t] * WC1(CdHR, p, t) + yddagI[r][t] * WC1(CdHI, p, t)) * yudagI[p][s]
                            +(yddagI[r][t] * WC1(CdHR, p, t) - yddagR[r][t] * WC1(CdHI, p, t)) * yudagR[p][s]);
                    ; //RGE 2
                    count ++;

                    //CdH : Re and Im
                    f[c + count * DF] += lambda * (
                            - 8. * (WC7R(Cqd1R, r, p, t, s) * yddagR[p][t] - WC7I(Cqd1I, r, p, t, s) * yddagI[p][t])
                            - 8. * cF3 * (WC7R(Cqd8R, r, p, t, s) * yddagR[p][t] - WC7I(Cqd8I, r, p, t, s) * yddagI[p][t])
                            + 4. * (WC5(CledqR, p, t, s, r) * yedagR[p][t] + WC5(CledqI, p, t, s, r) * yedagI[p][t])
                            + 4. * NC * (WC5(Cquqd1R, p, t, r, s) * yuR[t][p] - WC5(Cquqd1I, p, t, r, s) * yuI[t][p])
                            + 2. * (WC5(Cquqd1R, r, t, p, s) * yuR[t][p] - WC5(Cquqd1I, r, t, p, s) * yuI[t][p])
                            + 2. * cF3 * (WC5(Cquqd8R, r, t, p, s) * yuR[t][p] - WC5(Cquqd8I, r, t, p, s) * yuI[t][p])
                            ) //RGE 1
                            + 8. * ((WC7R(Cqd1R, r, p, t, s) + cF3 * WC7R(Cqd8R, r, p, t, s)) * yddydyddR[p][t]
                            -(WC7I(Cqd1I, r, p, t, s) + cF3 * WC7I(Cqd8I, r, p, t, s)) * yddydyddI[p][t])
                            - 4. * (WC5(CledqR, p, t, s, r) * yedyeyedR[p][t]
                            + WC5(CledqI, p, t, s, r) * yedyeyedI[p][t])
                            - 2. * (
                            (2. * NC * WC5(Cquqd1R, t, p, r, s) + WC5(Cquqd1R, r, p, t, s)
                            + cF3 * WC5(Cquqd8R, r, p, t, s)) * yuyudyuR[p][t]
                            -(2. * NC * WC5(Cquqd1I, t, p, r, s) + WC5(Cquqd1I, r, p, t, s)
                            + cF3 * WC5(Cquqd8I, r, p, t, s)) * yuyudyuI[p][t]
                            )
                            - 2. * ((yudagR[r][t] * WC1(CuHR, p, t) + yudagI[r][t] * WC1(CuHI, p, t)) * yddagR[p][s]
                            -(yudagI[r][t] * WC1(CuHR, p, t) - yudagR[r][t] * WC1(CuHI, p, t)) * yddagI[p][s])
                            ; //RGE 2
                    count ++;

                    f[c + count * DF] += lambda * (
                            - 8. * (WC7R(Cqd1R, r, p, t, s) * yddagI[p][t] + WC7I(Cqd1I, r, p, t, s) * yddagR[p][t])
                            - 8. * cF3 * (WC7R(Cqd8R, r, p, t, s) * yddagI[p][t] + WC7I(Cqd8I, r, p, t, s) * yddagR[p][t])
                            + 4. * (WC5(CledqR, p, t, s, r) * yedagI[p][t] - WC5(CledqI, p, t, s, r) * yedagR[p][t])
                            + 4. * NC * (WC5(Cquqd1R, p, t, r, s) * yuI[t][p] + WC5(Cquqd1I, p, t, r, s) * yuR[t][p])
                            + 2. * (WC5(Cquqd1R, r, t, p, s) * yuI[t][p] + WC5(Cquqd1I, r, t, p, s) * yuR[t][p])
                            + 2. * cF3 * (WC5(Cquqd8R, r, t, p, s) * yuI[t][p] + WC5(Cquqd8I, r, t, p, s) * yuR[t][p])
                            ) //RGE 1
                            + 8. * ((WC7R(Cqd1R, r, p, t, s) + cF3 * WC7R(Cqd8R, r, p, t, s)) * yddydyddI[p][t]
                            +(WC7I(Cqd1I, r, p, t, s) + cF3 * WC7I(Cqd8I, r, p, t, s)) * yddydyddR[p][t])
                            - 4. * (WC5(CledqR, p, t, s, r) * yedyeyedI[p][t]
                            - WC5(CledqI, p, t, s, r) * yedyeyedR[p][t])
                            - 2. * (
                            (2. * NC * WC5(Cquqd1R, t, p, r, s) + WC5(Cquqd1R, r, p, t, s)
                            + cF3 * WC5(Cquqd8R, r, p, t, s)) * yuyudyuI[p][t]
                            +(2. * NC * WC5(Cquqd1I, t, p, r, s) + WC5(Cquqd1I, r, p, t, s)
                            + cF3 * WC5(Cquqd8I, r, p, t, s)) * yuyudyuR[p][t]
                            )
                            - 2. * ((yudagR[r][t] * WC1(CuHR, p, t) + yudagI[r][t] * WC1(CuHI, p, t)) * yddagI[p][s]
                            +(yudagI[r][t] * WC1(CuHR, p, t) - yudagR[r][t] * WC1(CuHI, p, t)) * yddagR[p][s])
                            ; //RGE 2
                    count ++;




                    //CeH : Re and Im
                    f[c + count * DF] += lambda * (- 8. * (WC7R(CleR, r, p, t, s) * yedagR[p][t] -
                            WC7I(CleI, r, p, t, s) * yedagI[p][t])
                            + 4. * NC * (WC5(CledqR, r, s, p, t) * yddagR[t][p] - WC5(CledqI, r, s, p, t) * yddagI[t][p])
                            - 4. * NC * (WC5(Clequ1R, r, s, p, t) * yuR[t][p] - WC5(Clequ1I, r, s, p, t) * yuI[t][p])
                            ) //RGE 1
                            + 8. * (WC7R(CleR, r, p, t, s) * yedyeyedR[p][t] - WC7I(CleI, r, p, t, s) * yedyeyedI[p][t])
                            - 4. * NC * (WC5(CledqR, r, s, p, t) * yddydyddR[t][p] - WC5(CledqI, r, s, p, t) * yddydyddI[t][p])
                            + 4. * NC * (WC5(Clequ1R, r, s, t, p) * yuyudyuR[p][t] - WC5(Clequ1I, r, s, t, p) * yuyudyuI[p][t])
                            ; //RGE 2

                    count ++;
                    f[c + count * DF] += lambda * (- 8. * (WC7R(CleR, r, p, t, s) * yedagI[p][t] +
                            WC7I(CleI, r, p, t, s) * yedagR[p][t])
                            + 4. * NC * (WC5(CledqR, r, s, p, t) * yddagI[t][p] + WC5(CledqI, r, s, p, t) * yddagR[t][p])
                            - 4. * NC * (WC5(Clequ1R, r, s, p, t) * yuI[t][p] + WC5(Clequ1I, r, s, p, t) * yuR[t][p])
                            ) //RGE 1
                            + 8. * (WC7R(CleR, r, p, t, s) * yedyeyedI[p][t] + WC7I(CleI, r, p, t, s) * yedyeyedR[p][t])
                            - 4. * NC * (WC5(CledqR, r, s, p, t) * yddydyddI[t][p] + WC5(CledqI, r, s, p, t) * yddydyddR[t][p])
                            + 4. * NC * (WC5(Clequ1R, r, s, t, p) * yuyudyuI[p][t] + WC5(Clequ1I, r, s, t, p) * yuyudyuR[p][t])
                            ; //RGE 2
                    count ++;

                }


            }


            for (l = 0; l < 2 * N5; l ++) {
                f[c + l * DF] *= loop_factor;
            }
            c ++;
        }
    }
    c += (N5 * 2 - 1) * DF;



    //----------------------RGE SMEFT class 6--------------------------
    for (i = 0; i < NG; i ++) {
        for (j = 0; j < NG; j ++) {
            //Entries without  matrix products 
            l = 0;
            //CeW : Re and Im;
            f[c + l * DF] = ((3. * cF2 - b02) * g22 +
                    (- 3. * Ye2 + 8. * YeYl - 3. * Yl2) * g12) * WC1(CeWR, i, j)
                    + g1g2 * (3. * Yl - Ye) * WC1(CeBR, i, j)
                    - g2 * (yeR[j][i] * CHW + yeI[j][i] * CHWT)
                    - g1 * (Yl + Ye)*(yeR[j][i] * CHWB + yeI[j][i] * CHWBT)
                    + gammaH * WC1(CeWR, i, j);
            l ++;
            f[c + l * DF] = ((3. * cF2 - b02) * g22 +
                    (- 3. * Ye2 + 8. * YeYl - 3. * Yl2) * g12) * WC1(CeWI, i, j)
                    + g1g2 * (3. * Yl - Ye) * WC1(CeBI, i, j)
                    - g2 * (yeR[j][i] * CHWT - yeI[j][i] * CHW)
                    - g1 * (Yl + Ye)*(yeR[j][i] * CHWBT - yeI[j][i] * CHWB)
                    + gammaH * WC1(CeWI, i, j);
            l ++;
            //CeB : Re and Im;
            f[c + l * DF] = (- 3. * cF2 * g22 +
                    (3. * Ye2 + 4. * YeYl + 3. * Yl2 - b01) * g12) * WC1(CeBR, i, j)
                    + 4. * cF2 * g1g2 * (3. * Yl - Ye) * WC1(CeWR, i, j)
                    - 2. * g1 * (Yl + Ye)*(yeR[j][i] * CHB + yeI[j][i] * CHBT)
                    - 1.5 * g2 * (yeR[j][i] * CHWB + yeI[j][i] * CHWBT)
                    + gammaH * WC1(CeBR, i, j);
            l ++;
            f[c + l * DF] = (- 3. * cF2 * g22 +
                    (3. * Ye2 + 4. * YeYl + 3. * Yl2 - b01) * g12) * WC1(CeBI, i, j)
                    + 4. * cF2 * g1g2 * (3. * Yl - Ye) * WC1(CeWI, i, j)
                    - 2. * g1 * (Yl + Ye)*(yeR[j][i] * CHBT - yeI[j][i] * CHB)
                    - 1.5 * g2 * (yeR[j][i] * CHWBT - yeI[j][i] * CHWB)
                    + gammaH * WC1(CeBI, i, j);
            l ++;
            //CuG: Re and Im
            f[c + l * DF] = ((10. * cF3 - 4. * cA3 - b03) * g32 - 3. * cF2 * g22
                    + (- 3. * Yu2 + 8. * YuYq - 3. * Yq2) * g12) * WC1(CuGR, i, j)
                    + 8. * cF2 * g2g3 * WC1(CuWR, i, j) + 4. * g1g3 * (Yu + Yq) * WC1(CuBR, i, j)
                    - 4. * g3 * (yuR[j][i] * CHG + yuI[j][i] * CHGT)
                    + 3. * g32 * cA3 * (yuR[j][i] * CG + yuI[j][i] * CGT)
                    + gammaH * WC1(CuGR, i, j);
            l ++;
            f[c + l * DF ] = ((10. * cF3 - 4. * cA3 - b03) * g32 - 3 * cF2 * g22
                    + (- 3. * Yu2 + 8. * YuYq - 3. * Yq2) * g12) * WC1(CuGI, i, j)
                    + 8. * cF2 * g2g3 * WC1(CuWI, i, j) + 4. * g1g3 * (Yu + Yq) * WC1(CuBI, i, j)
                    - 4. * g3 * (yuR[j][i] * CHGT - yuI[j][i] * CHG)
                    + 3. * g32 * cA3 * (yuR[j][i] * CGT - yuI[j][i] * CG)
                    + gammaH * WC1(CuGI, i, j);
            l ++;
            //CuW: Re and Im
            f[c + l * DF ] = (2. * cF3 * g32 + (3. * cF2 - b02) * g22
                    + (- 3. * Yu2 + 8. * YuYq - 3. * Yq2) * g12) * WC1(CuWR, i, j)
                    + 2. * cF3 * g2g3 * WC1(CuGR, i, j) + g1g2 * (3. * Yq - Yu) * WC1(CuBR, i, j)
                    - g2 * (yuR[j][i] * CHW + yuI[j][i] * CHWT)
                    + g1 * (Yq + Yu) * (yuR[j][i] * CHWB + yuI[j][i] * CHWBT)
                    + gammaH * WC1(CuWR, i, j);
            l ++;
            f[c + l * DF ] = (2. * cF3 * g32 + (3. * cF2 - b02) * g22
                    + (- 3. * Yu2 + 8. * YuYq - 3. * Yq2) * g12) * WC1(CuWI, i, j)
                    + 2. * cF3 * g2g3 * WC1(CuGI, i, j) + g1g2 * (3. * Yq - Yu) * WC1(CuBI, i, j)
                    - g2 * (yuR[j][i] * CHWT - yuI[j][i] * CHW)
                    + g1 * (Yq + Yu) * (yuR[j][i] * CHWBT - yuI[j][i] * CHWB)
                    + gammaH * WC1(CuWI, i, j);
            l ++;
            //CuB: Re and Im
            f[c + l * DF ] = (2. * cF3 * g32 - 3. * cF2 * g22 +
                    (3. * Yu2 + 4. * YuYq + 3. * Yq2 - b01) * g12) * WC1(CuBR, i, j)
                    + 4. * cF2 * g1g2 * (3. * Yq - Yu) * WC1(CuWR, i, j)
                    + 4. * cF3 * g1g3 * (Yq + Yu) * WC1(CuGR, i, j)
                    -(yuR[j][i]*(2. * g1 * (Yq + Yu) * CHB - 1.5 * g2 * CHWB)
                    + yuI[j][i]*(2. * g1 * (Yq + Yu) * CHBT - 1.5 * g2 * CHWBT))
                    + gammaH * WC1(CuBR, i, j);
            l ++;
            f[c + l * DF ] = (2. * cF3 * g32 - 3. * cF2 * g22 +
                    (3. * Yu2 + 4. * YuYq + 3. * Yq2 - b01) * g12) * WC1(CuBI, i, j)
                    + 4. * cF2 * g1g2 * (3. * Yq - Yu) * WC1(CuWI, i, j)
                    + 4. * cF3 * g1g3 * (Yq + Yu) * WC1(CuGI, i, j)
                    - (yuR[j][i]*(2. * g1 * (Yq + Yu) * CHBT - 1.5 * g2 * CHWBT)
                    - yuI[j][i]*(2. * g1 * (Yq + Yu) * CHB - 1.5 * g2 * CHWB))
                    + gammaH * WC1(CuBI, i, j);
            l ++;


            //CdG : Re and Im;
            f[c + l * DF] = ((10. * cF3 - 4. * cA3 - b03) * g32 - 3. * cF2 * g22
                    + (- 3. * Yd2 + 8. * YdYq - 3. * Yq2) * g12) * WC1(CdGR, i, j)
                    + 8. * cF2 * g2g3 * WC1(CdWR, i, j) + 4. * g1g3 * (Yd + Yq) * WC1(CdBR, i, j)
                    - 4. * g3 * (ydR[j][i] * CHG + ydI[j][i] * CHGT)
                    + 3. * g32 * cA3 * (ydR[j][i] * CG + ydI[j][i] * CGT)
                    + gammaH * WC1(CdGR, i, j);
            l ++;
            f[c + l * DF] = ((10. * cF3 - 4. * cA3 - b03) * g32 - 3. * cF2 * g22
                    + (- 3. * Yd2 + 8. * YdYq - 3. * Yq2) * g12) * WC1(CdGI, i, j)
                    + 8. * cF2 * g2g3 * WC1(CdWI, i, j) + 4. * g1g3 * (Yd + Yq) * WC1(CdBI, i, j)
                    - 4. * g3 * (ydR[j][i] * CHGT - ydI[j][i] * CHG)
                    + 3. * g32 * cA3 * (ydR[j][i] * CGT - ydI[j][i] * CG)
                    + gammaH * WC1(CdGI, i, j);

            l ++;
            //CdW : Re and Im;
            f[c + l * DF] = (2. * cF3 * g32 + (3. * cF2 - b02) * g22
                    + (- 3. * Yd2 + 8. * YdYq - 3. * Yq2) * g12) * WC1(CdWR, i, j)
                    + 2. * cF3 * g2g3 * WC1(CdGR, i, j) + g1g2 * (3. * Yq - Yd) * WC1(CdBR, i, j)
                    - g2 * (ydR[j][i] * CHW + ydI[j][i] * CHWT)
                    - g1 * (Yq + Yd)*(ydR[j][i] * CHWB + ydI[j][i] * CHWBT)
                    + gammaH * WC1(CdWR, i, j);
            l ++;
            f[c + l * DF] = (2. * cF3 * g32 + (3. * cF2 - b02) * g22
                    + (- 3. * Yd2 + 8. * YdYq - 3. * Yq2) * g12) * WC1(CdWI, i, j)
                    + 2. * cF3 * g2g3 * WC1(CdGI, i, j) + g1g2 * (3. * Yq - Yd) * WC1(CdBI, i, j)
                    - g2 * (ydR[j][i] * CHWT - ydI[j][i] * CHW)
                    - g1 * (Yq + Yd)*(ydR[j][i] * CHWBT - ydI[j][i] * CHWB)
                    + gammaH * WC1(CdWI, i, j);

            l ++;
            //CdB : Re and Im;
            f[c + l * DF] = (2. * cF3 * g32 - 3. * cF2 * g22 +
                    (3. * Yd2 + 4. * YdYq + 3. * Yq2 - b01) * g12) * WC1(CdBR, i, j)
                    + 4. * cF2 * g1g2 * (3. * Yq - Yd) * WC1(CdWR, i, j)
                    + 4. * cF3 * g1g3 * (Yq + Yd) * WC1(CdGR, i, j)
                    - 2. * g1 * (Yq + Yd)*(ydR[j][i] * CHB + ydI[j][i] * CHBT)
                    - 1.5 * g2 * (ydR[j][i] * CHWB + ydI[j][i] * CHWBT)
                    + gammaH * WC1(CdBR, i, j);
            l ++;
            f[c + l * DF] = (2. * cF3 * g32 - 3. * cF2 * g22 +
                    (3. * Yd2 + 4. * YdYq + 3. * Yq2 - b01) * g12) * WC1(CdBI, i, j)
                    + 4. * cF2 * g1g2 * (3. * Yq - Yd) * WC1(CdWI, i, j)
                    + 4. * cF3 * g1g3 * (Yq + Yd) * WC1(CdGI, i, j)
                    - 2. * g1 * (Yq + Yd)*(ydR[j][i] * CHBT - ydI[j][i] * CHB)
                    - 1.5 * g2 * (ydR[j][i] * CHWBT - ydI[j][i] * CHWB)
                    + gammaH * WC1(CdBI, i, j);

            l ++;
            //Entries with 1 matrix product (1 summed index)

            for (a = 0; a < NG; a ++) {
                l = 0;
                //CeW : Re and Im;
                f[c + l * DF] += 2. * (WC1(CeWR, i, a) * geR[a][j] - WC1(CeWI, i, a) * geI[a][j])
                        +(glR[i][a] * WC1(CeWR, a, j) - glI[i][a] * WC1(CeWI, a, j));
                l ++;
                f[c + l * DF] += 2. * (WC1(CeWR, i, a) * geI[a][j] + WC1(CeWI, i, a) * geR[a][j])
                        +(glR[i][a] * WC1(CeWI, a, j) + glI[i][a] * WC1(CeWR, a, j));
                l ++;

                //CeB : Re and Im;
                f[c + l * DF] += 2. * (WC1(CeBR, i, a) * geR[a][j] - WC1(CeBI, i, a) * geI[a][j])
                        + 2. * (yedyeR[i][a] * WC1(CeBR, a, j) - yedyeI[i][a] * WC1(CeBI, a, j))
                        +((glR[i][a] * WC1(CeBR, a, j) - glI[i][a] * WC1(CeBI, a, j)));
                l ++;
                f[c + l * DF] += 2. * (WC1(CeBR, i, a) * geI[a][j] + WC1(CeBI, i, a) * geR[a][j])
                        + 2. * (yedyeR[i][a] * WC1(CeBI, a, j) + yedyeI[i][a] * WC1(CeBR, a, j))
                        +((glR[i][a] * WC1(CeBI, a, j) + glI[i][a] * WC1(CeBR, a, j)));
                l ++;

                //CuG: Re and Im
                f[c + l * DF] += 2. * ((yudyuR[i][a] - yddydR[i][a]) * WC1(CuGR, a, j)- (yudyuI[i][a] - yddydI[i][a]) * WC1(CuGI, a, j))
                        + 2. * (WC1(CuGR, i, a) * guR[a][j] - WC1(CuGI, i, a) * guI[a][j])
                        + (gqR[i][a] * WC1(CuGR, a, j) - gqI[i][a] * WC1(CuGI, a, j))
                        -(WC1(CdGR, i, a) * ydyudR[a][j] - WC1(CdGI, i, a) * ydyudI[a][j]);
                l ++;
                f[c + l * DF] += 2. * ((yudyuR[i][a] - yddydR[i][a]) * WC1(CuGI, a, j)+ (yudyuI[i][a] - yddydI[i][a]) * WC1(CuGR, a, j))
                        + 2. * (WC1(CuGR, i, a) * guI[a][j] + WC1(CuGI, i, a) * guR[a][j])
                        + (gqR[i][a] * WC1(CuGI, a, j) + gqI[i][a] * WC1(CuGR, a, j))
                        -(WC1(CdGI, i, a) * ydyudR[a][j] + WC1(CdGR, i, a) * ydyudI[a][j]);
                l ++;

                //CuW: Re and Im
                f[c + l * DF] += 2. * (yddydR[i][a] * WC1(CuWR, a, j) - yddydI[i][a] * WC1(CuWI, a, j))
                        + 2. * (WC1(CuWR, i, a) * guR[a][j] - WC1(CuWI, i, a) * guI[a][j])
                        + (gqR[i][a] * WC1(CuWR, a, j) - gqI[i][a] * WC1(CuWI, a, j))
                        -(WC1(CdWR, i, a) * ydyudR[a][j] - WC1(CdWI, i, a) * ydyudI[a][j]);
                l ++;

                f[c + l * DF] += 2. * (yddydI[i][a] * WC1(CuWR, a, j) + yddydR[i][a] * WC1(CuWI, a, j))
                        + 2. * (WC1(CuWR, i, a) * guI[a][j] + WC1(CuWI, i, a) * guR[a][j])
                        + (gqI[i][a] * WC1(CuWR, a, j) + gqR[i][a] * WC1(CuWI, a, j))
                        -(WC1(CdWR, i, a) * ydyudI[a][j] + WC1(CdWI, i, a) * ydyudR[a][j]);
                l ++;

                //CuB: Re and Im
                f[c + l * DF] += 2. * ((yudyuR[i][a] - yddydR[i][a]) * WC1(CuBR, a, j)- (yudyuI[i][a] - yddydI[i][a]) * WC1(CuBI, a, j))
                        + 2. * (WC1(CuBR, i, a) * guR[a][j] - WC1(CuBI, i, a) * guI[a][j])
                        + (gqR[i][a] * WC1(CuBR, a, j) - gqI[i][a] * WC1(CuBI, a, j))
                        -(WC1(CdBR, i, a) * ydyudR[a][j] - WC1(CdBI, i, a) * ydyudI[a][j]);
                l ++;


                f[c + l * DF] += 2. * ((yudyuR[i][a] - yddydR[i][a]) * WC1(CuBI, a, j) + (yudyuI[i][a] - yddydI[i][a]) * WC1(CuBR, a, j))
                        + 2. * (WC1(CuBR, i, a) * guI[a][j] + WC1(CuBI, i, a) * guR[a][j])
                        + (gqR[i][a] * WC1(CuBI, a, j) + gqI[i][a] * WC1(CuBR, a, j))
                        -(WC1(CdBR, i, a) * ydyudI[a][j] + WC1(CdBI, i, a) * ydyudR[a][j]);
                l ++;

                //CdG: Re and Im 
                f[c + l * DF] += - 2. * ((yudyuR[i][a] - yddydR[i][a]) * WC1(CdGR, a, j)-(yudyuI[i][a] - yddydI[i][a]) * WC1(CdGI, a, j))
                        -(WC1(CuGR, i, a) * yuyddR[a][j] - WC1(CuGI, i, a) * yuyddI[a][j]) +
                        2. * (WC1(CdGR, i, a) * gdR[a][j] - WC1(CdGI, i, a) * gdI[a][j])
                        +(gqR[i][a] * WC1(CdGR, a, j) - gqI[i][a] * WC1(CdGI, a, j));
                l ++;
                f[c + l * DF] += - 2. * ((yudyuR[i][a] - yddydR[i][a]) * WC1(CdGI, a, j)+(yudyuI[i][a] - yddydI[i][a]) * WC1(CdGR, a, j))
                        -(WC1(CuGR, i, a) * yuyddI[a][j] + WC1(CuGI, i, a) * yuyddR[a][j]) +
                        2. * (WC1(CdGR, i, a) * gdI[a][j] + WC1(CdGI, i, a) * gdR[a][j])
                        +(gqR[i][a] * WC1(CdGI, a, j) + gqI[i][a] * WC1(CdGR, a, j));
                l ++;

                //CdW: Re and Im 
                f[c + l * DF] += 2. * (yudyuR[i][a] * WC1(CdWR, a, j) - yudyuI[i][a] * WC1(CdWI, a, j))
                        -(WC1(CuWR, i, a) * yuyddR[a][j] - WC1(CuWI, i, a) * yuyddI[a][j])
                        + 2. * (WC1(CdWR, i, a) * gdR[a][j] - WC1(CdWI, i, a) * gdI[a][j])
                        +(gqR[i][a] * WC1(CdWR, a, j) - gqI[i][a] * WC1(CdWI, a, j));

                l ++;
                f[c + l * DF] += 2. * (yudyuR[i][a] * WC1(CdWI, a, j) + yudyuI[i][a] * WC1(CdWR, a, j))
                        -(WC1(CuWR, i, a) * yuyddI[a][j] + WC1(CuWI, i, a) * yuyddR[a][j])
                        + 2. * (WC1(CdWR, i, a) * gdI[a][j] + WC1(CdWI, i, a) * gdR[a][j])
                        +(gqR[i][a] * WC1(CdWI, a, j) + gqI[i][a] * WC1(CdWR, a, j));
                l ++;

                // CdB: Re and Im 
                f[c + l * DF] += - 2. * ((yudyuR[i][a] - yddydR[i][a]) * WC1(CdBR, a, j)-(yudyuI[i][a] - yddydI[i][a]) * WC1(CdBI, a, j))
                        -(WC1(CuBR, i, a) * yuyddR[a][j] - WC1(CuBI, i, a) * yuyddI[a][j])
                        + 2. * (WC1(CdBR, i, a) * gdR[a][j] - WC1(CdBI, i, a) * gdI[a][j])
                        +(gqR[i][a] * WC1(CdBR, a, j) - gqI[i][a] * WC1(CdBI, a, j));
                l ++;
                f[c + l * DF] += - 2. * ((yudyuI[i][a] - yddydI[i][a]) * WC1(CdBR, a, j)+(yudyuR[i][a] - yddydR[i][a]) * WC1(CdBI, a, j))
                        -(WC1(CuBR, i, a) * yuyddI[a][j] + WC1(CuBI, i, a) * yuyddR[a][j])
                        + 2. * (WC1(CdBR, i, a) * gdI[a][j] + WC1(CdBI, i, a) * gdR[a][j])
                        +(gqI[i][a] * WC1(CdBR, a, j) + gqR[i][a] * WC1(CdBI, a, j));
                l ++;

                //Entries with 2 matrix products (2 summed indices) 
                for (b = 0; b < NG; b ++) {
                    l = 0;
                    //CeW : Re and Im;
                    f[c + l * DF] += - 2. * g2 * NC * (WC5(Clequ3R, i, j, a, b) * yuR[b][a]
                            - WC5(Clequ3I, i, j, a, b) * yuI[b][a]);
                    l ++;
                    f[c + l * DF] += - 2. * g2 * NC * (WC5(Clequ3R, i, j, a, b) * yuI[b][a]
                            + WC5(Clequ3I, i, j, a, b) * yuR[b][a]);
                    l ++;
                    //CeB : Re and Im;
                    f[c + l * DF] += 4. * g1 * NC * (Yu + Yq)*
                            (WC5(Clequ3R, i, j, a, b) * yuR[b][a]
                            - WC5(Clequ3I, i, j, a, b) * yuI[b][a]);
                    l ++;
                    f[c + l * DF] += 4. * g1 * NC * (Yu + Yq)*
                            (WC5(Clequ3R, i, j, a, b) * yuI[b][a]
                            + WC5(Clequ3I, i, j, a, b) * yuR[b][a]);
                    l ++;
                    //CuG: Re and Im
                    f[c + l * DF] += - g3 * (
                            (WC5(Cquqd1R, a, j, i, b) -(0.5 / NC) * WC5(Cquqd8R, a, j, i, b)) * ydR[b][a]
                            - (WC5(Cquqd1I, a, j, i, b) -(0.5 / NC) * WC5(Cquqd8I, a, j, i, b)) * ydI[b][a]);
                    l ++;
                    f[c + l * DF] += - g3 * (
                            (WC5(Cquqd1R, a, j, i, b) -(0.5 / NC) * WC5(Cquqd8R, a, j, i, b)) * ydI[b][a]
                            + (WC5(Cquqd1I, a, j, i, b) -(0.5 / NC) * WC5(Cquqd8I, a, j, i, b)) * ydR[b][a]);
                    l ++;
                    //CuW: Re and Im
                    f[c + l * DF] += - 2. * g2 *
                            (WC5(Clequ3R, a, b, i, j) * yeR[b][a] - WC5(Clequ3I, a, b, i, j) * yeI[b][a])
                            + 0.25 * g2 * (
                            (WC5(Cquqd1R, a, j, i, b) + cF3 * WC5(Cquqd8R, a, j, i, b)) * ydR[b][a]
                            - (WC5(Cquqd1I, a, j, i, b) + cF3 * WC5(Cquqd8I, a, j, i, b)) * ydI[b][a]
                            );
                    l ++;
                    f[c + l * DF] += - 2. * g2 *
                            (WC5(Clequ3R, a, b, i, j) * yeI[b][a] + WC5(Clequ3I, a, b, i, j) * yeR[b][a])
                            + 0.25 * g2 * (
                            (WC5(Cquqd1R, a, j, i, b) + cF3 * WC5(Cquqd8R, a, j, i, b)) * ydI[b][a]
                            + (WC5(Cquqd1I, a, j, i, b) + cF3 * WC5(Cquqd8I, a, j, i, b)) * ydR[b][a]
                            );
                    l ++;
                    //CuB: Re and Im
                    f[c + l * DF] += 4. * g1 * (Ye + Yl)*(
                            WC5(Clequ3R, a, b, i, j) * yeR[b][a] - WC5(Clequ3I, a, b, i, j) * yeI[b][a])
                            - 0.5 * g1 * (Yd + Yq)*(
                            (WC5(Cquqd1R, a, j, i, b) + cF3 * WC5(Cquqd8R, a, j, i, b)) * ydR[b][a]
                            - (WC5(Cquqd1I, a, j, i, b) + cF3 * WC5(Cquqd8I, a, j, i, b)) * ydI[b][a]
                            );
                    l ++;
                    f[c + l * DF] += 4. * g1 * (Ye + Yl)*(
                            WC5(Clequ3R, a, b, i, j) * yeI[b][a] + WC5(Clequ3I, a, b, i, j) * yeR[b][a])
                            - 0.5 * g1 * (Yd + Yq)*(
                            (WC5(Cquqd1R, a, j, i, b) + cF3 * WC5(Cquqd8R, a, j, i, b)) * ydI[b][a]
                            + (WC5(Cquqd1I, a, j, i, b) + cF3 * WC5(Cquqd8I, a, j, i, b)) * ydR[b][a]
                            );
                    l ++;
                    //CdG: Re and Im 
                    f[c + l * DF] += - g3 * (
                            (WC5(Cquqd1R, i, b, a, j)-(0.5 / NC) * WC5(Cquqd8R, i, b, a, j)) * yuR[b][a]
                            -(WC5(Cquqd1I, i, b, a, j)-(0.5 / NC) * WC5(Cquqd8I, i, b, a, j)) * yuI[b][a]);
                    l ++;
                    f[c + l * DF] += - g3 * (
                            (WC5(Cquqd1R, i, b, a, j)-(0.5 / NC) * WC5(Cquqd8R, i, b, a, j)) * yuI[b][a]
                            +(WC5(Cquqd1I, i, b, a, j)-(0.5 / NC) * WC5(Cquqd8I, i, b, a, j)) * yuR[b][a]);
                    l ++;
                    //CdW: Re and Im 
                    f[c + l * DF] += 0.25 * g2 * (
                            (WC5(Cquqd1R, i, b, a, j) + cF3 * WC5(Cquqd8R, i, b, a, j)) * yuR[b][a]
                            - (WC5(Cquqd1I, i, b, a, j) + cF3 * WC5(Cquqd8I, i, b, a, j)) * yuI[b][a]
                            );
                    l ++;
                    f[c + l * DF] += 0.25 * g2 * (
                            (WC5(Cquqd1R, i, b, a, j) + cF3 * WC5(Cquqd8R, i, b, a, j)) * yuI[b][a]
                            + (WC5(Cquqd1I, i, b, a, j) + cF3 * WC5(Cquqd8I, i, b, a, j)) * yuR[b][a]
                            );
                    l ++;
                    //CdB: Re and Im 
                    f[c + l * DF] += - 0.5 * g1 * (Yu + Yq)*(
                            (WC5(Cquqd1R, i, b, a, j) + cF3 * WC5(Cquqd8R, i, b, a, j)) * yuR[b][a]
                            -(WC5(Cquqd1I, i, b, a, j) + cF3 * WC5(Cquqd8I, i, b, a, j)) * yuI[b][a]
                            );
                    l ++;
                    f[c + l * DF] += - 0.5 * g1 * (Yu + Yq)*(
                            (WC5(Cquqd1R, i, b, a, j) + cF3 * WC5(Cquqd8R, i, b, a, j)) * yuI[b][a]
                            +(WC5(Cquqd1I, i, b, a, j) + cF3 * WC5(Cquqd8I, i, b, a, j)) * yuR[b][a]
                            );
                    l ++;
                }
            }

            for (l = 0; l < 2 * N6; l ++) {
                f[c + l * DF] *= loop_factor;
            }

            c ++;
        }
    }
    c += (N6 * 2 - 1) * DF;


    //----------------------RGE SMEFT class 7--------------------------
    {

        //CHl1R
        for (i = 0; i < DWC2R; i ++) {
            p = WC2R_indices[i][0];
            r = WC2R_indices[i][1];
            //RGE 2
            f[c] = 0.;
            f[c] = - 0.5 * yedyeR[p][r] *(CHBOX + CHD) + 2. * gammaH * WC2R(CHl1R, p, r);
            for (t = 0; t < NG; t ++) {
                f[c] += 1.5 * (yedyeR[p][t] *(WC2R(CHl1R, t, r) + 3. * WC2R(CHl3R, t, r))
                        - yedyeI[p][t] *(WC2I(CHl1I, t, r) + 3. * WC2I(CHl3I, t, r)))
                        + 1.5 * ((WC2R(CHl1R, p, t) + 3. * WC2R(CHl3R, p, t)) * yedyeR[t][r]
                        -(WC2I(CHl1I, p, t) + 3. * WC2I(CHl3I, p, t)) * yedyeI[t][r])
                        +(glR[p][t] * WC2R(CHl1R, t, r) - glI[p][t] * WC2I(CHl1I, t, r))
                        +(WC2R(CHl1R, p, t) * glR[t][r] - WC2I(CHl1I, p, t) * glI[t][r]);

                for (s = 0; s < NG; s ++) {
                    f[c] += - ((yedagR[p][s] * WC2R(CHeR, s, t) - yedagI[p][s] * WC2I(CHeI, s, t)) * yeR[t][r]
                            -(yedagI[p][s] * WC2R(CHeR, s, t) + yedagR[p][s] * WC2I(CHeI, s, t)) * yeI[t][r])
                            + 2. * (WC7R(CleR, p, r, s, t) * geR[t][s] - WC7I(CleI, p, r, s, t) * geI[t][s])
                            +((- 2. * WC6R(CllR, p, r, s, t) - 2. * WC6R(CllR, s, t, p, r)
                            - WC6R(CllR, p, t, s, r) - WC6R(CllR, s, r, p, t)) * yedyeR[t][s]
                            -(- 2. * WC6I(CllI, p, r, s, t) - 2. * WC6I(CllI, s, t, p, r)
                            - WC6I(CllI, p, t, s, r) - WC6R(CllI, s, r, p, t)) * yedyeI[t][s])
                            - 2. * NC * (WC7R(Clq1R, p, r, s, t)*(yddydR[t][s] - yudyuR[t][s])
                            - WC7I(Clq1I, p, r, s, t)*(yddydI[t][s] - yudyuI[t][s]))
                            - 2. * NC * (WC7R(CluR, p, r, s, t) * guR[t][s] - WC7I(CluI, p, r, s, t) * guI[t][s])
                            + 2. * NC * (WC7R(CldR, p, r, s, t) * gdR[t][s] - WC7I(CldI, p, r, s, t) * gdI[t][s]);
                }
            }
            //Necessary since [RGE2] and [RGE3] use different indices.
            r = WC2R_indices[i][0];
            s = WC2R_indices[i][1];
            //RGE 3
            f[c] += 0.5 * xiB * g12 * delta[r][s] * Yl
                    + FOUR_THIRDS * g12 * Yh2 * WC2R(CHl1R, r, s);
            for (w = 0; w < NG; w ++) {
                f[c] += FOUR_THIRDS * g12 * Yh * NC * Yd * WC7R(CldR, r, s, w, w)
                        + FOUR_THIRDS * g12 * YhYe * WC7R(CleR, r, s, w, w)
                        + FOUR_THIRDS * g12 * YhYl * (
                        2. * WC6R(CllR, r, s, w, w) + WC6R(CllR, r, w, w, s)
                        + WC6R(CllR, w, s, r, w) + 2. * WC6R(CllR, w, w, r, s))
                        + EIGHT_THIRDS * g12 * NC * YhYq * WC7R(Clq1R, r, s, w, w)
                        + FOUR_THIRDS * g12 * NC * YhYu * WC7R(CluR, r, s, w, w);
            }

            f[c] *= loop_factor;
            c ++;
        }
        //CHl1I
        for (i = 0; i < DWC2I; i ++) {
            p = WC2I_indices[i][0];
            r = WC2I_indices[i][1];
            //RGE 2
            f[c] = - 0.5 * yedyeI[p][r] *(CHBOX + CHD) + 2. * gammaH * WC2I(CHl1I, p, r);
            for (t = 0; t < NG; t ++) {
                f[c] += 1.5 * (yedyeI[p][t] *(WC2R(CHl1R, t, r) + 3. * WC2R(CHl3R, t, r))
                        + yedyeR[p][t] *(WC2I(CHl1I, t, r) + 3. * WC2I(CHl3I, t, r)))
                        + 1.5 * ((WC2R(CHl1R, p, t) + 3. * WC2R(CHl3R, p, t)) * yedyeI[t][r]
                        +(WC2I(CHl1I, p, t) + 3. * WC2I(CHl3I, p, t)) * yedyeR[t][r])
                        +(glI[p][t] * WC2R(CHl1R, t, r) + glR[p][t] * WC2I(CHl1I, t, r))
                        +(WC2R(CHl1R, p, t) * glI[t][r] + WC2I(CHl1I, p, t) * glR[t][r]);

                for (s = 0; s < NG; s ++) {
                    f[c] += - ((yedagR[p][s] * WC2R(CHeR, s, t) - yedagI[p][s] * WC2I(CHeI, s, t)) * yeI[t][r]
                            +(yedagI[p][s] * WC2R(CHeR, s, t) + yedagR[p][s] * WC2I(CHeI, s, t)) * yeR[t][r])
                            + 2. * (WC7R(CleR, p, r, s, t) * geI[t][s] + WC7I(CleI, p, r, s, t) * geR[t][s])
                            +((- 2. * WC6R(CllR, p, r, s, t) - 2. * WC6R(CllR, s, t, p, r)
                            - WC6R(CllR, p, t, s, r) - WC6R(CllR, s, r, p, t)) * yedyeI[t][s]
                            +(- 2. * WC6I(CllI, p, r, s, t) - 2. * WC6I(CllI, s, t, p, r)
                            - WC6I(CllI, p, t, s, r) - WC6R(CllI, s, r, p, t)) * yedyeR[t][s])
                            - 2. * NC * (WC7R(Clq1R, p, r, s, t)*(yddydI[t][s] - yudyuI[t][s])
                            + WC7I(Clq1I, p, r, s, t)*(yddydR[t][s] - yudyuR[t][s]))
                            - 2. * NC * (WC7R(CluR, p, r, s, t) * guI[t][s] + WC7I(CluI, p, r, s, t) * guR[t][s])
                            + 2. * NC * (WC7R(CldR, p, r, s, t) * gdI[t][s] + WC7I(CldI, p, r, s, t) * gdR[t][s]);
                }
            }

            //Necessary since RGE 2 and RGE 3 use different indices.
            r = WC2I_indices[i][0];
            s = WC2I_indices[i][1];
            //RGE 3
            f[c] += FOUR_THIRDS * g12 * Yh2 * WC2I(CHl1I, r, s);
            for (w = 0; w < NG; w ++) {
                f[c] += FOUR_THIRDS * g12 * Yh * NC * Yd * WC7I(CldI, r, s, w, w)
                        + FOUR_THIRDS * g12 * YhYe * WC7I(CleI, r, s, w, w)
                        + FOUR_THIRDS * g12 * YhYl * (
                        2. * WC6I(CllI, r, s, w, w) + WC6I(CllI, r, w, w, s)
                        + WC6I(CllI, w, s, r, w) + 2. * WC6I(CllI, w, w, r, s))
                        + EIGHT_THIRDS * g12 * NC * YhYq * WC7I(Clq1I, r, s, w, w)
                        + FOUR_THIRDS * g12 * NC * YhYu * WC7I(CluI, r, s, w, w);
            }


            f[c] *= loop_factor;
            c ++;
        }



        //CHl3R
        for (i = 0; i < DWC2R; i ++) {
            p = WC2R_indices[i][0];
            r = WC2R_indices[i][1];
            //RGE 2
            f[c] = - 0.5 * yedyeR[p][r] * CHBOX + 2. * gammaH * WC2R(CHl3R, p, r);
            for (t = 0; t < NG; t ++) {
                f[c] += 0.5 * (yedyeR[p][t]*(3. * WC2R(CHl1R, t, r) + WC2R(CHl3R, t, r))
                        - yedyeI[p][t]*(3. * WC2I(CHl1I, t, r) + WC2I(CHl3I, t, r)))
                        + 0.5 * ((3. * WC2R(CHl1R, p, t) + WC2R(CHl3R, p, t)) * yedyeR[t][r]
                        -(3. * WC2I(CHl1I, p, t) + WC2I(CHl3I, p, t)) * yedyeI[t][r])
                        +(glR[p][t] * WC2R(CHl3R, t, r) - glI[p][t] * WC2I(CHl3I, t, r))
                        +(WC2R(CHl3R, p, t) * glR[t][r] - WC2I(CHl3I, p, t) * glI[t][r]);
                for (s = 0; s < NG; s ++) {
                    f[c] += - ((WC6R(CllR, p, t, s, r) + WC6R(CllR, s, r, p, t)) * yedyeR[t][s]
                            -(WC6I(CllI, p, t, s, r) + WC6I(CllI, s, r, p, t)) * yedyeI[t][s])
                            - 2. * NC * (WC7R(Clq3R, p, r, s, t)*(yddydR[t][s] + yudyuR[t][s])
                            - WC7I(Clq3I, p, r, s, t)*(yddydI[t][s] + yudyuI[t][s]));
                }
            }
            //Necessary since RGE 2 and RGE 3 use different indices.
            r = WC2R_indices[i][0];
            s = WC2R_indices[i][1];
            //RGE 3
            f[c] += ONE_SIXTH * g22 * CHBOX * delta[r][s] + ONE_THIRD * g22 * WC2R(CHl3R, r, s)
                    - 6. * g22 * WC2R(CHl3R, r, s);
            for (w = 0; w < NG; w ++) {
                f[c] += TWO_THIRDS * g22 * delta[r][s]*
                        (WC2R(CHl3R, w, w) + NC * WC2R(CHq3R, w, w))
                        + ONE_THIRD * g22 * (WC6R(CllR, r, w, w, s) + WC6R(CllR, w, s, r, w))
                        + TWO_THIRDS * g22 * NC * WC7R(Clq3R, r, s, w, w);
            }
            f[c] *= loop_factor;
            c ++;
        }


        //CHl3I
        for (i = 0; i < DWC2I; i ++) {
            p = WC2I_indices[i][0];
            r = WC2I_indices[i][1];
            //RGE 2
            f[c] = - 0.5 * yedyeI[p][r] * CHBOX + 2. * gammaH * WC2I(CHl3I, p, r);
            for (t = 0; t < NG; t ++) {
                f[c] += 0.5 * (yedyeI[p][t]*(3. * WC2R(CHl1R, t, r) + WC2R(CHl3R, t, r))
                        + yedyeR[p][t]*(3. * WC2I(CHl1I, t, r) + WC2I(CHl3I, t, r)))
                        + 0.5 * ((3. * WC2R(CHl1R, p, t) + WC2R(CHl3R, p, t)) * yedyeI[t][r]
                        +(3. * WC2I(CHl1I, p, t) + WC2I(CHl3I, p, t)) * yedyeR[t][r])
                        +(glI[p][t] * WC2R(CHl3R, t, r) + glR[p][t] * WC2I(CHl3I, t, r))
                        +(WC2R(CHl3R, p, t) * glI[t][r] + WC2I(CHl3I, p, t) * glR[t][r]);
                for (s = 0; s < NG; s ++) {
                    f[c] += - ((WC6R(CllR, p, t, s, r) + WC6R(CllR, s, r, p, t)) * yedyeI[t][s]
                            +(WC6I(CllI, p, t, s, r) + WC6I(CllI, s, r, p, t)) * yedyeR[t][s])
                            - 2. * NC * (WC7R(Clq3R, p, r, s, t)*(yddydI[t][s] + yudyuI[t][s])
                            + WC7I(Clq3I, p, r, s, t)*(yddydR[t][s] + yudyuR[t][s]));
                }
            }
            //Necessary since RGE 2 and RGE 3 use different indices.
            r = WC2I_indices[i][0];
            s = WC2I_indices[i][1];
            //RGE 3  
            f[c] += ONE_THIRD * g22 * WC2I(CHl3I, r, s)
                    - 6. * g22 * WC2I(CHl3I, r, s);
            for (w = 0; w < NG; w ++) {
                f[c] += ONE_THIRD * g22 * (WC6I(CllI, r, w, w, s) + WC6I(CllI, w, s, r, w))
                        + TWO_THIRDS * g22 * NC * WC7I(Clq3I, r, s, w, w);
            }
            f[c] *= loop_factor;
            c ++;
        }


        //CHeR
        for (i = 0; i < DWC2R; i ++) {
            p = WC2R_indices[i][0];
            r = WC2R_indices[i][1];
            //RGE 2
            f[c] = geR[p][r]*(CHBOX + CHD) + 2. * gammaH * WC2R(CHeR, p, r);
            for (t = 0; t < NG; t ++) {
                f[c] += 4. * (geR[p][t] * WC2R(CHeR, t, r) - geI[p][t] * WC2I(CHeI, t, r))
                        + 4. * (WC2R(CHeR, p, t) * geR[t][r] - WC2I(CHeI, p, t) * geI[t][r]);

                for (s = 0; s < NG; s ++) {
                    f[c] += - 2. * ((yeR[p][s] * WC2R(CHl1R, s, t) - yeI[p][s] * WC2I(CHl1I, s, t)) * yedagR[t][r]
                            - (yeI[p][s] * WC2R(CHl1R, s, t) + yeR[p][s] * WC2I(CHl1I, s, t)) * yedagI[t][r])
                            - 2. * (WC7R(CleR, s, t, p, r) * yedyeR[t][s] - WC7I(CleI, s, t, p, r) * yedyeI[t][s])
                            + 2. * ((WC8R(CeeR, p, r, s, t) + WC8R(CeeR, s, t, p, r)
                            + WC8R(CeeR, p, t, s, r) + WC8R(CeeR, s, r, p, t)) * geR[t][s]
                            -(WC8I(CeeI, p, r, s, t) + WC8I(CeeI, s, t, p, r)
                            + WC8I(CeeI, p, t, s, r) + WC8I(CeeI, s, r, p, t)) * geI[t][s])
                            - 2. * NC * (WC7R(CeuR, p, r, s, t) * guR[t][s] - WC7I(CeuI, p, r, s, t) * guI[t][s])
                            + 2. * NC * (WC7R(CedR, p, r, s, t) * gdR[t][s] - WC7I(CedI, p, r, s, t) * gdI[t][s])
                            - 2. * NC * (WC7R(CqeR, s, t, p, r)*(yddydR[t][s] - yudyuR[t][s])
                            - WC7I(CqeI, s, t, p, r)*(yddydI[t][s] - yudyuI[t][s]));
                }
            }
            //Necessary since RGE 2 and RGE 3 use different indices.
            r = WC2R_indices[i][0];
            s = WC2R_indices[i][1];
            //RGE 3
            f[c] += 0.5 * xiB * g12 * delta[r][s] * Ye
                    + FOUR_THIRDS * g12 * Yh2 * WC2R(CHeR, r, s);
            for (w = 0; w < NG; w ++) {
                f[c] += FOUR_THIRDS * NC * g12 * YhYd * WC7R(CedR, r, s, w, w)
                        + FOUR_THIRDS * g12 * YhYe *
                        (WC8R(CeeR, r, s, w, w) + WC8R(CeeR, r, w, w, s)
                        + WC8R(CeeR, w, s, r, w) + WC8R(CeeR, w, w, r, s))
                        + FOUR_THIRDS * g12 * NC * YhYu * WC7R(CeuR, r, s, w, w)
                        + EIGHT_THIRDS * g12 * YhYl * WC7R(CleR, w, w, r, s)
                        + EIGHT_THIRDS * g12 * NC * YhYq * WC7R(CqeR, w, w, r, s);
            }
            f[c] *= loop_factor;
            c ++;
        }


        //CHeI
        for (i = 0; i < DWC2I; i ++) {
            p = WC2I_indices[i][0];
            r = WC2I_indices[i][1];
            //RGE 2
            f[c] = geI[p][r]*(CHBOX + CHD) + 2. * gammaH * WC2I(CHeI, p, r);
            for (t = 0; t < NG; t ++) {
                f[c] += 4. * (geI[p][t] * WC2R(CHeR, t, r) + geR[p][t] * WC2I(CHeI, t, r))
                        + 4. * (WC2R(CHeR, p, t) * geI[t][r] + WC2I(CHeI, p, t) * geR[t][r]);

                for (s = 0; s < NG; s ++) {
                    f[c] += - 2. * ((yeR[p][s] * WC2R(CHl1R, s, t) - yeI[p][s] * WC2I(CHl1I, s, t)) * yedagI[t][r]
                            +(yeI[p][s] * WC2R(CHl1R, s, t) + yeR[p][s] * WC2I(CHl1I, s, t)) * yedagR[t][r])
                            - 2. * (WC7R(CleR, s, t, p, r) * yedyeI[t][s] + WC7I(CleI, s, t, p, r) * yedyeR[t][s])
                            + 2. * ((WC8R(CeeR, p, r, s, t) + WC8R(CeeR, s, t, p, r)
                            + WC8R(CeeR, p, t, s, r) + WC8R(CeeR, s, r, p, t)) * geI[t][s]
                            +(WC8I(CeeI, p, r, s, t) + WC8I(CeeI, s, t, p, r)
                            + WC8I(CeeI, p, t, s, r) + WC8I(CeeI, s, r, p, t)) * geR[t][s])
                            - 2. * NC * (WC7R(CeuR, p, r, s, t) * guI[t][s] + WC7I(CeuI, p, r, s, t) * guR[t][s])
                            + 2. * NC * (WC7R(CedR, p, r, s, t) * gdI[t][s] + WC7I(CedI, p, r, s, t) * gdR[t][s])
                            - 2. * NC * (WC7R(CqeR, s, t, p, r)*(yddydI[t][s] - yudyuI[t][s])
                            + WC7I(CqeI, s, t, p, r)*(yddydR[t][s] - yudyuR[t][s]));
                }
            }
            //Necessary since RGE 2 and RGE 3 use different indices.
            r = WC2I_indices[i][0];
            s = WC2I_indices[i][1];
            //RGE 3
            f[c] += FOUR_THIRDS * g12 * Yh2 * WC2I(CHeI, r, s);
            for (w = 0; w < NG; w ++) {
                f[c] += FOUR_THIRDS * NC * g12 * YhYd * WC7I(CedI, r, s, w, w)
                        + FOUR_THIRDS * g12 * YhYe *
                        (WC8I(CeeI, r, s, w, w) + WC8I(CeeI, r, w, w, s)
                        + WC8I(CeeI, w, s, r, w) + WC8I(CeeI, w, w, r, s))
                        + FOUR_THIRDS * g12 * NC * YhYu * WC7I(CeuI, r, s, w, w)
                        + EIGHT_THIRDS * g12 * YhYl * WC7I(CleI, w, w, r, s)
                        + EIGHT_THIRDS * g12 * NC * YhYq * WC7I(CqeI, w, w, r, s);
            }
            f[c] *= loop_factor;
            c ++;
        }

        //CHq1R
        for (i = 0; i < DWC2R; i ++) {
            p = WC2R_indices[i][0];
            r = WC2R_indices[i][1];
            //RGE 2
            f[c] = 0.5 * (yudyuR[p][r] - yddydR[p][r])*(CHBOX + CHD) + 2. * gammaH * WC2R(CHq1R, p, r);
            for (t = 0; t < NG; t ++) {
                f[c] += 1.5 * ((yddydR[p][t] + yudyuR[p][t]) * WC2R(CHq1R, t, r)
                        -(yddydI[p][t] + yudyuI[p][t]) * WC2I(CHq1I, t, r))
                        + 1.5 * (WC2R(CHq1R, p, t)*(yddydR[t][r] + yudyuR[t][r])
                        - WC2I(CHq1I, p, t)*(yddydI[t][r] + yudyuI[t][r]))
                        + 4.5 * ((yddydR[p][t] - yudyuR[p][t]) * WC2R(CHq3R, t, r)
                        -(yddydI[p][t] - yudyuI[p][t]) * WC2I(CHq3I, t, r))
                        + 4.5 * (WC2R(CHq3R, p, t)*(yddydR[t][r] - yudyuR[t][r])
                        - WC2I(CHq3I, p, t)*(yddydI[t][r] - yudyuI[t][r]))
                        +(gqR[p][t] * WC2R(CHq1R, t, r) - gqI[p][t] * WC2I(CHq1I, t, r))
                        +(WC2R(CHq1R, p, t) * gqR[t][r] - WC2I(CHq1I, p, t) * gqI[t][r])
                        ;
                for (s = 0; s < NG; s ++) {
                    f[c] += - ((yudagR[p][s] * WC2R(CHuR, s, t) - yudagI[p][s] * WC2I(CHuI, s, t)) * yuR[t][r]
                            -(yudagI[p][s] * WC2R(CHuR, s, t) + yudagR[p][s] * WC2I(CHuI, s, t)) * yuI[t][r])
                            -((yddagR[p][s] * WC2R(CHdR, s, t) - yddagI[p][s] * WC2I(CHdI, s, t)) * ydR[t][r]
                            -(yddagI[p][s] * WC2R(CHdR, s, t) + yddagR[p][s] * WC2I(CHdI, s, t)) * ydI[t][r])
                            + 2. * (WC7R(CqeR, p, r, s, t) * geR[t][s] - WC7I(CqeI, p, r, s, t) * geI[t][s])
                            - 2. * (WC7R(Clq1R, s, t, p, r) * yedyeR[t][s] - WC7I(Clq1I, s, t, p, r) * yedyeI[t][s])
                            -(
                            (2. * NC * WC6R(Cqq1R, p, r, s, t) + 2. * NC * WC6R(Cqq1R, s, t, p, r)
                            + WC6R(Cqq1R, p, t, s, r) + WC6R(Cqq1R, s, r, p, t)
                            + 3. * WC6R(Cqq3R, p, t, s, r) + 3. * WC6R(Cqq3R, s, r, p, t))*(yddydR[t][s] - yudyuR[t][s])
                            -(2. * NC * WC6I(Cqq1I, p, r, s, t) + 2. * NC * WC6I(Cqq1I, s, t, p, r)
                            + WC6I(Cqq1I, p, t, s, r) + WC6I(Cqq1I, s, r, p, t)
                            + 3. * WC6I(Cqq3I, p, t, s, r) + 3. * WC6I(Cqq3I, s, r, p, t))*(yddydI[t][s] - yudyuI[t][s])
                            )
                            - 2. * NC * (WC7R(Cqu1R, p, r, s, t) * guR[t][s] - WC7I(Cqu1I, p, r, s, t) * guI[t][s])
                            + 2. * NC * (WC7R(Cqd1R, p, r, s, t) * gdR[t][s] - WC7I(Cqd1I, p, r, s, t) * gdI[t][s]);

                }
            }
            //Necessary since RGE 2 and RGE 3 use different indices.
            r = WC2R_indices[i][0];
            s = WC2R_indices[i][1];
            f[c] += 0.5 * xiB * g12 * delta[r][s] * Yq
                    + FOUR_THIRDS * g12 * Yh2 * WC2R(CHq1R, r, s);
            for (w = 0; w < NG; w ++) {
                f[c] += EIGHT_THIRDS * g12 * YhYl * WC7R(Clq1R, w, w, r, s)
                        + FOUR_THIRDS * g12 * NC * YhYd * WC7R(Cqd1R, r, s, w, w)
                        + FOUR_THIRDS * g12 * YhYe * WC7R(CqeR, r, s, w, w)
                        + FOUR_THIRDS * g12 * YhYq * (
                        2. * NC * WC6R(Cqq1R, r, s, w, w) + WC6R(Cqq1R, r, w, w, s)
                        + WC6R(Cqq1R, w, s, r, w) + 2. * NC * WC6R(Cqq1R, w, w, r, s))
                        + 4. * g12 * YhYq *
                        (WC6R(Cqq3R, r, w, w, s) + WC6R(Cqq3R, w, s, r, w))
                        + FOUR_THIRDS * g12 * NC * YhYu * WC7R(Cqu1R, r, s, w, w);
            }

            f[c] *= loop_factor;
            c ++;
        }



        //CHq1I
        for (i = 0; i < DWC2I; i ++) {
            p = WC2I_indices[i][0];
            r = WC2I_indices[i][1];
            //RGE 2
            f[c] = 0.5 * (yudyuI[p][r] - yddydI[p][r])*(CHBOX + CHD) + 2. * gammaH * WC2I(CHq1I, p, r);
            for (t = 0; t < NG; t ++) {
                f[c] += 1.5 * ((yddydI[p][t] + yudyuI[p][t]) * WC2R(CHq1R, t, r)
                        +(yddydR[p][t] + yudyuR[p][t]) * WC2I(CHq1I, t, r))
                        + 1.5 * (WC2R(CHq1R, p, t)*(yddydI[t][r] + yudyuI[t][r])
                        + WC2I(CHq1I, p, t)*(yddydR[t][r] + yudyuR[t][r]))
                        + 4.5 * ((yddydI[p][t] - yudyuI[p][t]) * WC2R(CHq3R, t, r)
                        +(yddydR[p][t] - yudyuR[p][t]) * WC2I(CHq3I, t, r))
                        + 4.5 * (WC2R(CHq3R, p, t)*(yddydI[t][r] - yudyuI[t][r])
                        + WC2I(CHq3I, p, t)*(yddydR[t][r] - yudyuR[t][r]))
                        +(gqI[p][t] * WC2R(CHq1R, t, r) + gqR[p][t] * WC2I(CHq1I, t, r))
                        +(WC2R(CHq1R, p, t) * gqI[t][r] + WC2I(CHq1I, p, t) * gqR[t][r])
                        ;
                for (s = 0; s < NG; s ++) {
                    f[c] += - ((yudagR[p][s] * WC2R(CHuR, s, t) - yudagI[p][s] * WC2I(CHuI, s, t)) * yuI[t][r]
                            +(yudagI[p][s] * WC2R(CHuR, s, t) + yudagR[p][s] * WC2I(CHuI, s, t)) * yuR[t][r])
                            -((yddagR[p][s] * WC2R(CHdR, s, t) - yddagI[p][s] * WC2I(CHdI, s, t)) * ydI[t][r]
                            +(yddagI[p][s] * WC2R(CHdR, s, t) + yddagR[p][s] * WC2I(CHdI, s, t)) * ydR[t][r])
                            + 2. * (WC7R(CqeR, p, r, s, t) * geI[t][s] + WC7I(CqeI, p, r, s, t) * geR[t][s])
                            - 2. * (WC7R(Clq1R, s, t, p, r) * yedyeI[t][s] + WC7I(Clq1I, s, t, p, r) * yedyeR[t][s])
                            -(
                            (2. * NC * WC6R(Cqq1R, p, r, s, t) + 2. * NC * WC6R(Cqq1R, s, t, p, r)
                            + WC6R(Cqq1R, p, t, s, r) + WC6R(Cqq1R, s, r, p, t)
                            + 3. * WC6R(Cqq3R, p, t, s, r) + 3. * WC6R(Cqq3R, s, r, p, t))*(yddydI[t][s] - yudyuI[t][s])
                            +(2. * NC * WC6I(Cqq1I, p, r, s, t) + 2. * NC * WC6I(Cqq1I, s, t, p, r)
                            + WC6I(Cqq1I, p, t, s, r) + WC6I(Cqq1I, s, r, p, t)
                            + 3. * WC6I(Cqq3I, p, t, s, r) + 3. * WC6I(Cqq3I, s, r, p, t))*(yddydR[t][s] - yudyuR[t][s])
                            )
                            - 2. * NC * (WC7R(Cqu1R, p, r, s, t) * guI[t][s] + WC7I(Cqu1I, p, r, s, t) * guR[t][s])
                            + 2. * NC * (WC7R(Cqd1R, p, r, s, t) * gdI[t][s] + WC7I(Cqd1I, p, r, s, t) * gdR[t][s]);

                }
            }
            //Necessary since RGE 2 and RGE 3 use different indices.
            r = WC2I_indices[i][0];
            s = WC2I_indices[i][1];
            f[c] += FOUR_THIRDS * g12 * Yh2 * WC2I(CHq1I, r, s);
            for (w = 0; w < NG; w ++) {
                f[c] += EIGHT_THIRDS * g12 * YhYl * WC7I(Clq1I, w, w, r, s)
                        + FOUR_THIRDS * g12 * NC * YhYd * WC7I(Cqd1I, r, s, w, w)
                        + FOUR_THIRDS * g12 * YhYe * WC7I(CqeI, r, s, w, w)
                        + FOUR_THIRDS * g12 * YhYq * (
                        2. * NC * WC6I(Cqq1I, r, s, w, w) + WC6I(Cqq1I, r, w, w, s)
                        + WC6I(Cqq1I, w, s, r, w) + 2. * NC * WC6I(Cqq1I, w, w, r, s))
                        + 4. * g12 * YhYq *
                        (WC6I(Cqq3I, r, w, w, s) + WC6I(Cqq3I, w, s, r, w))
                        + FOUR_THIRDS * g12 * NC * YhYu * WC7I(Cqu1I, r, s, w, w);
            }
            f[c] *= loop_factor;
            c ++;
        }


        //CHq3R 
        for (i = 0; i < DWC2R; i ++) {
            p = WC2R_indices[i][0];
            r = WC2R_indices[i][1];
            //RGE 2
            f[c] = - 0.5 * (yudyuR[p][r] + yddydR[p][r]) * CHBOX + 2. * gammaH * WC2R(CHq3R, p, r);
            for (t = 0; t < NG; t ++) {
                f[c] += 1.5 * ((yddydR[p][t] - yudyuR[p][t]) * WC2R(CHq1R, t, r)
                        -(yddydI[p][t] - yudyuI[p][t]) * WC2I(CHq1I, t, r))
                        + 1.5 * (WC2R(CHq1R, p, t)*(yddydR[t][r] - yudyuR[t][r])
                        - WC2I(CHq1I, p, t)*(yddydI[t][r] - yudyuI[t][r]))

                        + 0.5 * ((yddydR[p][t] + yudyuR[p][t]) * WC2R(CHq3R, t, r)
                        -(yddydI[p][t] + yudyuI[p][t]) * WC2I(CHq3I, t, r))
                        + 0.5 * (WC2R(CHq3R, p, t)*(yddydR[t][r] + yudyuR[t][r])
                        - WC2I(CHq3I, p, t)*(yddydI[t][r] + yudyuI[t][r]))

                        +(gqR[p][t] * WC2R(CHq3R, t, r) - gqI[p][t] * WC2I(CHq3I, t, r))
                        +(WC2R(CHq3R, p, t) * gqR[t][r] - WC2I(CHq3I, p, t) * gqI[t][r]);


                for (s = 0; s < NG; s ++) {
                    f[c] += - (
                            (2. * NC * WC6R(Cqq3R, p, r, s, t) + 2. * NC * WC6R(Cqq3R, s, t, p, r)
                            + WC6R(Cqq1R, p, t, s, r) + WC6R(Cqq1R, s, r, p, t)
                            - WC6R(Cqq3R, p, t, s, r) - WC6R(Cqq3R, s, r, p, t))*(yddydR[t][s] + yudyuR[t][s])
                            -(2. * NC * WC6I(Cqq3I, p, r, s, t) + 2. * NC * WC6I(Cqq3I, s, t, p, r)
                            + WC6I(Cqq1I, p, t, s, r) + WC6I(Cqq1I, s, r, p, t)
                            - WC6I(Cqq3I, p, t, s, r) - WC6I(Cqq3I, s, r, p, t))* (yddydI[t][s] + yudyuI[t][s])
                            )
                            - 2. * (WC7R(Clq3R, s, t, p, r) * yedyeR[t][s] - WC7I(Clq3I, s, t, p, r) * yedyeI[t][s]);
                }
            }
            //Necessary since RGE 2 and RGE 3 use different indices.
            r = WC2R_indices[i][0];
            s = WC2R_indices[i][1];
            //RGE 3
            f[c] += ONE_SIXTH * g22 * CHBOX * delta[r][s]
                    + ONE_THIRD * g22 * WC2R(CHq3R, r, s) - 6. * g22 * WC2R(CHq3R, r, s);
            for (w = 0; w < NG; w ++) {
                f[c] += TWO_THIRDS * g22 * delta[r][s]*
                        (WC2R(CHl3R, w, w) + NC * WC2R(CHq3R, w, w))
                        + TWO_THIRDS * g22 * WC7R(Clq3R, w, w, r, s)
                        + ONE_THIRD * g22 * (WC6R(Cqq1R, r, w, w, s) + WC6R(Cqq1R, w, s, r, w))
                        + ONE_THIRD * g22 *
                        (2. * NC * WC6R(Cqq3R, r, s, w, w) - WC6R(Cqq3R, r, w, w, s)
                        - WC6R(Cqq3R, w, s, r, w) + 2. * NC * WC6R(Cqq3R, w, w, r, s));
            }
            f[c] *= loop_factor;
            c ++;
        }

        //CHq3I
        for (i = 0; i < DWC2I; i ++) {
            p = WC2I_indices[i][0];
            r = WC2I_indices[i][1];
            //RGE 2
            f[c] = - 0.5 * (yudyuI[p][r] + yddydI[p][r]) * CHBOX + 2. * gammaH * WC2I(CHq3I, p, r);
            for (t = 0; t < NG; t ++) {
                f[c] += 1.5 * ((yddydI[p][t] - yudyuI[p][t]) * WC2R(CHq1R, t, r)
                        +(yddydR[p][t] - yudyuR[p][t]) * WC2I(CHq1I, t, r))
                        + 1.5 * (WC2R(CHq1R, p, t)*(yddydI[t][r] - yudyuI[t][r])
                        + WC2I(CHq1I, p, t)*(yddydR[t][r] - yudyuR[t][r]))

                        + 0.5 * ((yddydI[p][t] + yudyuI[p][t]) * WC2R(CHq3R, t, r)
                        +(yddydR[p][t] + yudyuR[p][t]) * WC2I(CHq3I, t, r))
                        + 0.5 * (WC2R(CHq3R, p, t)*(yddydI[t][r] + yudyuI[t][r])
                        + WC2I(CHq3I, p, t)*(yddydR[t][r] + yudyuR[t][r]))

                        +(gqI[p][t] * WC2R(CHq3R, t, r) + gqR[p][t] * WC2I(CHq3I, t, r))
                        +(WC2R(CHq3R, p, t) * gqI[t][r] + WC2I(CHq3I, p, t) * gqR[t][r]);


                for (s = 0; s < NG; s ++) {
                    f[c] += - (
                            (2. * NC * WC6R(Cqq3R, p, r, s, t) + 2. * NC * WC6R(Cqq3R, s, t, p, r)
                            + WC6R(Cqq1R, p, t, s, r) + WC6R(Cqq1R, s, r, p, t)
                            - WC6R(Cqq3R, p, t, s, r) - WC6R(Cqq3R, s, r, p, t))*(yddydI[t][s] + yudyuI[t][s])
                            +(2. * NC * WC6I(Cqq3I, p, r, s, t) + 2. * NC * WC6I(Cqq3I, s, t, p, r)
                            + WC6I(Cqq1I, p, t, s, r) + WC6I(Cqq1I, s, r, p, t)
                            - WC6I(Cqq3I, p, t, s, r) - WC6I(Cqq3I, s, r, p, t))* (yddydR[t][s] + yudyuR[t][s])
                            )
                            - 2. * (WC7R(Clq3R, s, t, p, r) * yedyeI[t][s] + WC7I(Clq3I, s, t, p, r) * yedyeR[t][s]);
                }
            }
            //Necessary since RGE 2 and RGE 3 use different indices.
            r = WC2I_indices[i][0];
            s = WC2I_indices[i][1];
            //RGE 3
            f[c] += ONE_THIRD * g22 * WC2I(CHq3I, r, s) - 6. * g22 * WC2I(CHq3I, r, s);
            for (w = 0; w < NG; w ++) {
                f[c] += TWO_THIRDS * g22 * WC7I(Clq3I, w, w, r, s)
                        + ONE_THIRD * g22 *
                        (WC6I(Cqq1I, r, w, w, s) + WC6I(Cqq1I, w, s, r, w))
                        + ONE_THIRD * g22 *
                        (2. * NC * WC6I(Cqq3I, r, s, w, w) - WC6I(Cqq3I, r, w, w, s)
                        - WC6I(Cqq3I, w, s, r, w) + 2. * NC * WC6I(Cqq3I, w, w, r, s));
            }
            f[c] *= loop_factor;
            c ++;
        }


        //CHuR
        for (i = 0; i < DWC2R; i ++) {
            p = WC2R_indices[i][0];
            r = WC2R_indices[i][1];
            //RGE 2
            f[c] = - guR[p][r]*(CHBOX + CHD) + 2. * gammaH * WC2R(CHuR, p, r);
            for (t = 0; t < NG; t ++) {
                f[c] += 4. * (guR[p][t] * WC2R(CHuR, t, r) - guI[p][t] * WC2I(CHuI, t, r))
                        + 4. * (WC2R(CHuR, p, t) * guR[t][r] - WC2I(CHuI, p, t) * guI[t][r])
                        +(yuyddR[p][t] * WC1(CHudR, r, t) + yuyddI[p][t] * WC1(CHudI, r, t))
                        +(WC1(CHudR, p, t) * ydyudR[t][r] - WC1(CHudI, p, t) * ydyudI[t][r])
                        ;
                for (s = 0; s < NG; s ++) {
                    f[c] += - 2. * ((yuR[p][s] * WC2R(CHq1R, s, t) - yuI[p][s] * WC2I(CHq1I, s, t)) * yudagR[t][r]
                            -(yuI[p][s] * WC2R(CHq1R, s, t) + yuR[p][s] * WC2I(CHq1I, s, t)) * yudagI[t][r])
                            - 2. * (
                            (NC * WC6R(CuuR, p, r, s, t) + NC * WC6R(CuuR, s, t, p, r)
                            + WC6R(CuuR, p, t, s, r) + WC6R(CuuR, s, r, p, t)) * guR[t][s]
                            -(NC * WC6I(CuuI, p, r, s, t) + NC * WC6I(CuuI, s, t, p, r)
                            + WC6I(CuuI, p, t, s, r) + WC6I(CuuI, s, r, p, t)) * guI[t][s])
                            + 2. * (WC7R(CeuR, s, t, p, r) * geR[t][s] - WC7I(CeuI, s, t, p, r) * geI[t][s])
                            - 2. * (WC7R(CluR, s, t, p, r) * yedyeR[t][s] - WC7I(CluI, s, t, p, r) * yedyeI[t][s])
                            + 2. * NC * (WC7R(Cud1R, p, r, s, t) * gdR[t][s] - WC7I(Cud1I, p, r, s, t) * gdI[t][s])
                            - 2. * NC * (WC7R(Cqu1R, s, t, p, r)*(yddydR[t][s] - yudyuR[t][s])
                            - WC7I(Cqu1I, s, t, p, r)*(yddydI[t][s] - yudyuI[t][s]));
                }
            }

            //Necessary since RGE 2 and RGE 3 use different indices.
            r = WC2R_indices[i][0];
            s = WC2R_indices[i][1];
            //RGE 3
            f[c] += 0.5 * xiB * g12 * delta[r][s] * Yu
                    + FOUR_THIRDS * g12 * Yh2 * WC2R(CHuR, r, s);
            for (w = 0; w < NG; w ++) {
                f[c] += FOUR_THIRDS * g12 * YhYe * WC7R(CeuR, w, w, r, s)
                        + EIGHT_THIRDS * g12 * YhYl * WC7R(CluR, w, w, r, s)
                        + EIGHT_THIRDS * g12 * NC * YhYq * WC7R(Cqu1R, w, w, r, s)
                        + FOUR_THIRDS * g12 * NC * YhYd * WC7R(Cud1R, r, s, w, w)
                        + FOUR_THIRDS * g12 * YhYu * (
                        NC * WC6R(CuuR, r, s, w, w) + WC6R(CuuR, r, w, w, s)
                        + WC6R(CuuR, w, s, r, w) + NC * WC6R(CuuR, w, w, r, s));

            }
            f[c] *= loop_factor;
            c ++;
        }
        //CHuI
        for (i = 0; i < DWC2I; i ++) {
            p = WC2I_indices[i][0];
            r = WC2I_indices[i][1];
            //RGE 2
            f[c] = - guI[p][r]*(CHBOX + CHD) + 2. * gammaH * WC2I(CHuI, p, r);
            for (t = 0; t < NG; t ++) {
                f[c] += 4. * (guI[p][t] * WC2R(CHuR, t, r) + guR[p][t] * WC2I(CHuI, t, r))
                        + 4. * (WC2R(CHuR, p, t) * guI[t][r] + WC2I(CHuI, p, t) * guR[t][r])
                        +(yuyddI[p][t] * WC1(CHudR, r, t) - yuyddR[p][t] * WC1(CHudI, r, t))
                        +(WC1(CHudR, p, t) * ydyudI[t][r] + WC1(CHudI, p, t) * ydyudR[t][r])
                        ;
                for (s = 0; s < NG; s ++) {
                    f[c] += - 2. * ((yuR[p][s] * WC2R(CHq1R, s, t) - yuI[p][s] * WC2I(CHq1I, s, t)) * yudagI[t][r]
                            +(yuI[p][s] * WC2R(CHq1R, s, t) + yuR[p][s] * WC2I(CHq1I, s, t)) * yudagR[t][r])
                            - 2. * (
                            (NC * WC6R(CuuR, p, r, s, t) + NC * WC6R(CuuR, s, t, p, r)
                            + WC6R(CuuR, p, t, s, r) + WC6R(CuuR, s, r, p, t)) * guI[t][s]
                            +(NC * WC6I(CuuI, p, r, s, t) + NC * WC6I(CuuI, s, t, p, r)
                            + WC6I(CuuI, p, t, s, r) + WC6I(CuuI, s, r, p, t)) * guR[t][s])
                            + 2. * (WC7R(CeuR, s, t, p, r) * geI[t][s] + WC7I(CeuI, s, t, p, r) * geR[t][s])
                            - 2. * (WC7R(CluR, s, t, p, r) * yedyeI[t][s] + WC7I(CluI, s, t, p, r) * yedyeR[t][s])
                            + 2. * NC * (WC7R(Cud1R, p, r, s, t) * gdI[t][s] + WC7I(Cud1I, p, r, s, t) * gdR[t][s])
                            - 2. * NC * (WC7R(Cqu1R, s, t, p, r)*(yddydI[t][s] - yudyuI[t][s])
                            + WC7I(Cqu1I, s, t, p, r)*(yddydR[t][s] - yudyuR[t][s]));
                }
            }
            //Necessary since RGE 2 and RGE 3 use different indices.
            r = WC2I_indices[i][0];
            s = WC2I_indices[i][1];
            //RGE 3
            f[c] += FOUR_THIRDS * g12 * Yh2 * WC2I(CHuI, r, s);
            for (w = 0; w < NG; w ++) {
                f[c] += FOUR_THIRDS * g12 * YhYe * WC7I(CeuI, w, w, r, s)
                        + EIGHT_THIRDS * g12 * YhYl * WC7I(CluI, w, w, r, s)
                        + EIGHT_THIRDS * g12 * NC * YhYq * WC7I(Cqu1I, w, w, r, s)
                        + FOUR_THIRDS * g12 * NC * YhYd * WC7I(Cud1I, r, s, w, w)
                        + FOUR_THIRDS * g12 * YhYu * (
                        NC * WC6I(CuuI, r, s, w, w) + WC6I(CuuI, r, w, w, s)
                        + WC6I(CuuI, w, s, r, w) + NC * WC6I(CuuI, w, w, r, s));

            }
            f[c] *= loop_factor;
            c ++;
        }




        //CHdR
        for (i = 0; i < DWC2R; i ++) {
            p = WC2R_indices[i][0];
            r = WC2R_indices[i][1];
            //RGE 2
            f[c] = gdR[p][r]*(CHBOX + CHD) + 2. * gammaH * WC2R(CHdR, p, r);
            for (t = 0; t < NG; t ++) {
                f[c] += 4. * (gdR[p][t] * WC2R(CHdR, t, r) - gdI[p][t] * WC2I(CHdI, t, r))
                        + 4. * (WC2R(CHdR, p, t) * gdR[t][r] - WC2I(CHdI, p, t) * gdI[t][r])
                        -(ydyudR[p][t] * WC1(CHudR, t, r) - ydyudI[p][t] * WC1(CHudI, t, r))
                        -(WC1(CHudR, t, p) * yuyddR[t][r] + WC1(CHudI, t, p) * yuyddI[t][r]);
                for (s = 0; s < NG; s ++) {
                    f[c] += - 2. * ((ydR[p][s] * WC2R(CHq1R, s, t) - ydI[p][s] * WC2I(CHq1I, s, t)) * yddagR[t][r]
                            -(ydR[p][s] * WC2I(CHq1I, s, t) + ydI[p][s] * WC2R(CHq1R, s, t)) * yddagI[t][r])
                            + 2. * (
                            (NC * WC6R(CddR, p, r, s, t) + NC * WC6R(CddR, s, t, p, r)
                            + WC6R(CddR, p, t, s, r) + WC6R(CddR, s, r, p, t)) * gdR[t][s]
                            -(NC * WC6I(CddI, p, r, s, t) + NC * WC6I(CddI, s, t, p, r)
                            + WC6I(CddI, p, t, s, r) + WC6I(CddI, s, r, p, t)) * gdI[t][s]
                            )
                            + 2. * (WC7R(CedR, s, t, p, r) * geR[t][s] - WC7I(CedI, s, t, p, r) * geI[t][s])
                            - 2. * (WC7R(CldR, s, t, p, r) * yedyeR[t][s] - WC7I(CldI, s, t, p, r) * yedyeI[t][s])
                            - 2. * NC * (WC7R(Cud1R, s, t, p, r) * guR[t][s] - WC7I(Cud1I, s, t, p, r) * guI[t][s])
                            - 2. * NC * (WC7R(Cqd1R, s, t, p, r)*(yddydR[t][s] - yudyuR[t][s])
                            - WC7I(Cqd1I, s, t, p, r)*(yddydI[t][s] - yudyuI[t][s]))
                            ;
                }
            }
            //Necessary since RGE 2 and RGE 3 use different indices.
            r = WC2R_indices[i][0];
            s = WC2R_indices[i][1];
            //RGE 3
            f[c] += 0.5 * xiB * g12 * delta[r][s] * Yd
                    + FOUR_THIRDS * g12 * Yh2 * WC2R(CHdR, r, s);
            for (w = 0; w < NG; w ++) {
                f[c] += FOUR_THIRDS * g12 * YhYd * (
                        NC * WC6R(CddR, r, s, w, w) + WC6R(CddR, r, w, w, s)
                        + WC6R(CddR, w, s, r, w) + NC * WC6R(CddR, w, w, r, s))
                        + FOUR_THIRDS * g12 * YhYe * WC7R(CedR, w, w, r, s)
                        + EIGHT_THIRDS * g12 * YhYl * WC7R(CldR, w, w, r, s)
                        + EIGHT_THIRDS * g12 * NC * YhYq * WC7R(Cqd1R, w, w, r, s)
                        + FOUR_THIRDS * g12 * NC * YhYu * WC7R(Cud1R, w, w, r, s);
            }
            f[c] *= loop_factor;
            c ++;
        }

        //CHdI  
        for (i = 0; i < DWC2I; i ++) {
            p = WC2I_indices[i][0];
            r = WC2I_indices[i][1];
            //RGE 2
            f[c] = gdI[p][r]*(CHBOX + CHD) + 2. * gammaH * WC2I(CHdI, p, r);
            for (t = 0; t < NG; t ++) {
                f[c] += 4. * (gdI[p][t] * WC2R(CHdR, t, r) + gdR[p][t] * WC2I(CHdI, t, r))
                        + 4. * (WC2R(CHdR, p, t) * gdI[t][r] + WC2I(CHdI, p, t) * gdR[t][r])
                        -(ydyudI[p][t] * WC1(CHudR, t, r) + ydyudR[p][t] * WC1(CHudI, t, r))
                        -(WC1(CHudR, t, p) * yuyddI[t][r] - WC1(CHudI, t, p) * yuyddR[t][r]);
                for (s = 0; s < NG; s ++) {
                    f[c] += - 2. * ((ydR[p][s] * WC2R(CHq1R, s, t) - ydI[p][s] * WC2I(CHq1I, s, t)) * yddagI[t][r]
                            +(ydR[p][s] * WC2I(CHq1I, s, t) + ydI[p][s] * WC2R(CHq1R, s, t)) * yddagR[t][r])
                            + 2. * (
                            (NC * WC6R(CddR, p, r, s, t) + NC * WC6R(CddR, s, t, p, r)
                            + WC6R(CddR, p, t, s, r) + WC6R(CddR, s, r, p, t)) * gdI[t][s]
                            +(NC * WC6I(CddI, p, r, s, t) + NC * WC6I(CddI, s, t, p, r)
                            + WC6I(CddI, p, t, s, r) + WC6I(CddI, s, r, p, t)) * gdR[t][s]
                            )
                            + 2. * (WC7R(CedR, s, t, p, r) * geI[t][s] + WC7I(CedI, s, t, p, r) * geR[t][s])
                            - 2. * (WC7R(CldR, s, t, p, r) * yedyeI[t][s] + WC7I(CldI, s, t, p, r) * yedyeR[t][s])
                            - 2. * NC * (WC7R(Cud1R, s, t, p, r) * guI[t][s] + WC7I(Cud1I, s, t, p, r) * guR[t][s])
                            - 2. * NC * (WC7R(Cqd1R, s, t, p, r)*(yddydI[t][s] - yudyuI[t][s])
                            + WC7I(Cqd1I, s, t, p, r)*(yddydR[t][s] - yudyuR[t][s]))
                            ;
                }
            }
            //Necessary since RGE 2 and RGE 3 use different indices.
            r = WC2I_indices[i][0];
            s = WC2I_indices[i][1];
            //RGE 3
            f[c] += FOUR_THIRDS * g12 * Yh2 * WC2I(CHdI, r, s);
            for (w = 0; w < NG; w ++) {
                f[c] += FOUR_THIRDS * g12 * YhYd * (
                        NC * WC6I(CddI, r, s, w, w) + WC6I(CddI, r, w, w, s)
                        + WC6I(CddI, w, s, r, w) + NC * WC6I(CddI, w, w, r, s))
                        + FOUR_THIRDS * g12 * YhYe * WC7I(CedI, w, w, r, s)
                        + EIGHT_THIRDS * g12 * YhYl * WC7I(CldI, w, w, r, s)
                        + EIGHT_THIRDS * g12 * NC * YhYq * WC7I(Cqd1I, w, w, r, s)
                        + FOUR_THIRDS * g12 * NC * YhYu * WC7I(Cud1I, w, w, r, s);
            }
            f[c] *= loop_factor;
            c ++;
        }

        for (p = 0; p < NG; p ++) {
            for (r = 0; r < NG; r ++) {
                //CHudR
                f[c] = yuyddR[p][r]*(2. * CHBOX - CHD) + 2. * gammaH * WC1(CHudR, p, r);
                for (t = 0; t < NG; t ++) {
                    f[c] += - 2. * (yuyddR[p][t] * WC2R(CHdR, t, r) - yuyddI[p][t] * WC2I(CHdI, t, r))
                            + 2. * (WC2R(CHuR, p, t) * yuyddR[r][t] - WC2I(CHuI, p, t) * yuyddI[r][t])
                            + 3. * (guR[p][t] * WC1(CHudR, t, r) - guI[p][t] * WC1(CHudI, t, r))
                            + 3. * (WC1(CHudR, p, t) * gdR[t][r] - WC1(CHudI, p, t) * gdI[t][r])
                            ;
                    for (s = 0; s < NG; s ++) {
                        f[c] += 4. * (
                                (WC7R(Cud1R, p, t, s, r) + cF3 * WC7R(Cud8R, p, t, s, r)) * yuyddR[t][s]
                                -(WC7I(Cud1I, p, t, s, r) + cF3 * WC7I(Cud8I, p, t, s, r)) * yuyddI[t][s]
                                );
                    }
                }
                //RGE 3 
                f[c] += - 3. * g12 * (Yu - Yd)*(Yu - Yd) * WC1(CHudR, p, r);

                f[c] *= loop_factor;

                //CHudI
                f[c + DF] = yuyddI[p][r]*(2. * CHBOX - CHD) + 2. * gammaH * WC1(CHudI, p, r);
                for (t = 0; t < NG; t ++) {
                    f[c + DF] += - 2. * (yuyddI[p][t] * WC2R(CHdR, t, r) + yuyddR[p][t] * WC2I(CHdI, t, r))
                            + 2. * (WC2R(CHuR, p, t) * yuyddI[r][t] + WC2I(CHuI, p, t) * yuyddR[r][t])
                            + 3. * (guI[p][t] * WC1(CHudR, t, r) + guR[p][t] * WC1(CHudI, t, r))
                            + 3. * (WC1(CHudR, p, t) * gdI[t][r] + WC1(CHudI, p, t) * gdR[t][r])
                            ;
                    for (s = 0; s < NG; s ++) {
                        f[c + DF] += 4. * (
                                (WC7R(Cud1R, p, t, s, r) + cF3 * WC7R(Cud8R, p, t, s, r)) * yuyddI[t][s]
                                +(WC7I(Cud1I, p, t, s, r) + cF3 * WC7I(Cud8I, p, t, s, r)) * yuyddR[t][s]
                                );
                    }
                }
                //RGE 3 
                f[c + DF] += - 3. * g12 * (Yu - Yd)*(Yu - Yd) * WC1(CHudI, p, r);


                f[c + DF] *= loop_factor;
                c ++;
            }
        }
        c += DF;

    }
    //----------------------RGE SMEFT class 8_LLLL----------------------
    //CllR
    for (d = 0; d < DWC6R; d ++) {
        p = WC6R_indices[d][0];
        r = WC6R_indices[d][1];
        s = WC6R_indices[d][2];
        t = WC6R_indices[d][3];

        f[c] = - 0.5 * (yedyeR[p][r] * WC2R(CHl1R, s, t) - yedyeI[p][r] * WC2I(CHl1I, s, t))
                - 0.5 * (yedyeR[s][t] * WC2R(CHl1R, p, r) - yedyeI[s][t] * WC2I(CHl1I, p, r))
                + 0.5 * (yedyeR[p][r] * WC2R(CHl3R, s, t) - yedyeI[p][r] * WC2I(CHl3I, s, t))
                + 0.5 * (yedyeR[s][t] * WC2R(CHl3R, p, r) - yedyeI[s][t] * WC2I(CHl3I, p, r))
                -(yedyeR[s][r] * WC2R(CHl3R, p, t) - yedyeI[s][r] * WC2I(CHl3I, p, t))
                -(yedyeR[p][t] * WC2R(CHl3R, s, r) - yedyeI[p][t] * WC2I(CHl3I, s, r))//RGE 2
                + TWO_THIRDS * g12 * YhYl * (WC2R(CHl1R, s, t) * delta[p][r] + WC2R(CHl1R, p, r) * delta[s][t])
                - ONE_SIXTH * g22 * (WC2R(CHl3R, s, t) * delta[p][r] + WC2R(CHl3R, p, r) * delta[s][t])
                + ONE_THIRD * g22 * (WC2R(CHl3R, s, r) * delta[p][t] + WC2R(CHl3R, p, t) * delta[r][s])
                + 6. * g22 * WC6R(CllR, p, t, s, r) - 3. * (g22 - 4. * Yl2 * g12) * WC6R(CllR, p, r, s, t)
                //RGE 3
                ;
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] += (glR[p][v] * WC6R(CllR, v, r, s, t) - glI[p][v] * WC6I(CllI, v, r, s, t))
                    +(glR[s][v] * WC6R(CllR, p, r, v, t) - glI[s][v] * WC6I(CllI, p, r, v, t))
                    +(WC6R(CllR, p, v, s, t) * glR[v][r] - WC6I(CllI, p, v, s, t) * glI[v][r])
                    +(WC6R(CllR, p, r, s, v) * glR[v][t] - WC6I(CllI, p, r, s, v) * glI[v][t])
                    //RGE 2
                    + FOUR_THIRDS * g12 * Yl2 * (
                    WC6R(CllR, p, r, w, w) * delta[s][t] + WC6R(CllR, s, t, w, w) * delta[p][r]
                    + WC6R(CllR, w, w, s, t) * delta[p][r] + WC6R(CllR, w, w, p, r) * delta[s][t])
                    + TWO_THIRDS * g12 * Yl2 * (
                    WC6R(CllR, p, w, w, r) * delta[s][t] + WC6R(CllR, s, w, w, t) * delta[p][r]
                    + WC6R(CllR, w, r, p, w) * delta[s][t] + WC6R(CllR, w, t, s, w) * delta[p][r])
                    - ONE_SIXTH * g22 * (
                    WC6R(CllR, p, w, w, r) * delta[s][t] + WC6R(CllR, s, w, w, t) * delta[p][r]
                    + WC6R(CllR, p, w, w, r) * delta[s][t] + WC6R(CllR, w, t, s, w) * delta[p][r])
                    + ONE_THIRD * g22 * (
                    WC6R(CllR, s, w, w, r) * delta[p][t] + WC6R(CllR, p, w, w, t) * delta[r][s]
                    + WC6R(CllR, w, r, s, w) * delta[p][t] + WC6R(CllR, w, t, p, w) * delta[r][s])
                    + FOUR_THIRDS * g12 * NC * YlYq * (
                    WC7R(Clq1R, p, r, w, w) * delta[s][t] + WC7R(Clq1R, s, t, w, w) * delta[p][r])
                    + ONE_THIRD * g22 * NC *
                    (- WC7R(Clq3R, p, r, w, w) * delta[s][t] - WC7R(Clq3R, s, t, w, w) * delta[p][r]
                    + 2. * WC7R(Clq3R, s, r, w, w) * delta[p][t] + 2. * WC7R(Clq3R, p, t, w, w) * delta[r][s])
                    + TWO_THIRDS * g12 * NC * YuYl * (
                    WC7R(CluR, p, r, w, w) * delta[s][t] + WC7R(CluR, s, t, w, w) * delta[p][r])
                    + TWO_THIRDS * g12 * NC * YdYl * (
                    WC7R(CldR, p, r, w, w) * delta[s][t] + WC7R(CldR, s, t, w, w) * delta[p][r])
                    + TWO_THIRDS * g12 * YeYl * (
                    WC7R(CleR, p, r, w, w) * delta[s][t] + WC7R(CleR, s, t, w, w) * delta[p][r]);
            //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] += - 0.5 * (
                        + (yedagR[s][v] * yeR[w][t] - yedagI[s][v] * yeI[w][t]) * WC7R(CleR, p, r, v, w)
                        -(yedagR[s][v] * yeI[w][t] + yedagI[s][v] * yeR[w][t]) * WC7I(CleI, p, r, v, w))
                        - 0.5 * ((yedagR[p][v] * yeR[w][r] - yedagI[p][v] * yeI[w][r]) * WC7R(CleR, s, t, v, w)
                        -(yedagR[p][v] * yeI[w][r] + yedagI[p][v] * yeR[w][r]) * WC7I(CleI, s, t, v, w))
                        //RGE 2
                        ;
            }
        }

        f[c] *= loop_factor;
        c ++;
    }
    //CllI
    for (d = 0; d < DWC6I; d ++) {
        p = WC6I_indices[d][0];
        r = WC6I_indices[d][1];
        s = WC6I_indices[d][2];
        t = WC6I_indices[d][3];

        f[c] = - 0.5 * (yedyeI[p][r] * WC2R(CHl1R, s, t) + yedyeR[p][r] * WC2I(CHl1I, s, t))
                - 0.5 * (yedyeI[s][t] * WC2R(CHl1R, p, r) + yedyeR[s][t] * WC2I(CHl1I, p, r))
                + 0.5 * (yedyeI[p][r] * WC2R(CHl3R, s, t) + yedyeR[p][r] * WC2I(CHl3I, s, t))
                + 0.5 * (yedyeI[s][t] * WC2R(CHl3R, p, r) + yedyeR[s][t] * WC2I(CHl3I, p, r))
                -(yedyeI[s][r] * WC2R(CHl3R, p, t) + yedyeR[s][r] * WC2I(CHl3I, p, t))
                -(yedyeI[p][t] * WC2R(CHl3R, s, r) + yedyeR[p][t] * WC2I(CHl3I, s, r))//RGE 2
                + TWO_THIRDS * g12 * YhYl * (WC2I(CHl1I, s, t) * delta[p][r] + WC2I(CHl1I, p, r) * delta[s][t])
                - ONE_SIXTH * g22 * (WC2I(CHl3I, s, t) * delta[p][r] + WC2I(CHl3I, p, r) * delta[s][t])
                + ONE_THIRD * g22 * (WC2I(CHl3I, s, r) * delta[p][t] + WC2I(CHl3I, p, t) * delta[r][s])
                + 6. * g22 * WC6I(CllI, p, t, s, r) - 3. * (g22 - 4. * Yl2 * g12) * WC6I(CllI, p, r, s, t)
                //RGE 3
                ;
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] += (glI[p][v] * WC6R(CllR, v, r, s, t) + glR[p][v] * WC6I(CllI, v, r, s, t))
                    +(glI[s][v] * WC6R(CllR, p, r, v, t) + glR[s][v] * WC6I(CllI, p, r, v, t))
                    +(WC6R(CllR, p, v, s, t) * glI[v][r] + WC6I(CllI, p, v, s, t) * glR[v][r])
                    +(WC6R(CllR, p, r, s, v) * glI[v][t] + WC6I(CllI, p, r, s, v) * glR[v][t])
                    //RGE 2
                    + FOUR_THIRDS * g12 * Yl2 * (
                    WC6I(CllI, p, r, w, w) * delta[s][t] + WC6I(CllI, s, t, w, w) * delta[p][r]
                    + WC6I(CllI, w, w, s, t) * delta[p][r] + WC6I(CllI, w, w, p, r) * delta[s][t])
                    + TWO_THIRDS * g12 * Yl2 * (
                    WC6I(CllI, p, w, w, r) * delta[s][t] + WC6I(CllI, s, w, w, t) * delta[p][r]
                    + WC6I(CllI, w, r, p, w) * delta[s][t] + WC6I(CllI, w, t, s, w) * delta[p][r])
                    - ONE_SIXTH * g22 * (
                    WC6I(CllI, p, w, w, r) * delta[s][t] + WC6I(CllI, s, w, w, t) * delta[p][r]
                    + WC6I(CllI, p, w, w, r) * delta[s][t] + WC6I(CllI, w, t, s, w) * delta[p][r])
                    + ONE_THIRD * g22 * (
                    WC6I(CllI, s, w, w, r) * delta[p][t] + WC6I(CllI, p, w, w, t) * delta[r][s]
                    + WC6I(CllI, w, r, s, w) * delta[p][t] + WC6I(CllI, w, t, p, w) * delta[r][s])
                    + FOUR_THIRDS * g12 * NC * YlYq * (
                    WC7I(Clq1I, p, r, w, w) * delta[s][t] + WC7I(Clq1I, s, t, w, w) * delta[p][r])
                    + ONE_THIRD * g22 * NC *
                    (- WC7I(Clq3I, p, r, w, w) * delta[s][t] - WC7I(Clq3I, s, t, w, w) * delta[p][r]
                    + 2. * WC7I(Clq3I, s, r, w, w) * delta[p][t] + 2. * WC7I(Clq3I, p, t, w, w) * delta[r][s])
                    + TWO_THIRDS * g12 * NC * YuYl * (
                    WC7I(CluI, p, r, w, w) * delta[s][t] + WC7I(CluI, s, t, w, w) * delta[p][r])
                    + TWO_THIRDS * g12 * NC * YdYl * (
                    WC7I(CldI, p, r, w, w) * delta[s][t] + WC7I(CldI, s, t, w, w) * delta[p][r])
                    + TWO_THIRDS * g12 * YeYl * (
                    WC7I(CleI, p, r, w, w) * delta[s][t] + WC7I(CleI, s, t, w, w) * delta[p][r]);
            //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] += - 0.5 * ((yedagR[s][v] * yeR[w][t] - yedagI[s][v] * yeI[w][t]) * WC7I(CleI, p, r, v, w)
                        +(yedagR[s][v] * yeI[w][t] + yedagI[s][v] * yeR[w][t]) * WC7R(CleR, p, r, v, w))
                        - 0.5 * ((yedagR[p][v] * yeR[w][r] - yedagI[p][v] * yeI[w][r]) * WC7I(CleI, s, t, v, w)
                        +(yedagR[p][v] * yeI[w][r] + yedagI[p][v] * yeR[w][r]) * WC7R(CleR, s, t, v, w))
                        //RGE 2
                        ;

            }
        }

        f[c] *= loop_factor;
        c ++;
    }
    //Cqq1R
    for (d = 0; d < DWC6R; d ++) {
        p = WC6R_indices[d][0];
        r = WC6R_indices[d][1];
        s = WC6R_indices[d][2];
        t = WC6R_indices[d][3];
        f[c] = 0.5 * ((yudyuR[p][r] - yddydR[p][r]) * WC2R(CHq1R, s, t)
                -(yudyuI[p][r] - yddydI[p][r]) * WC2I(CHq1I, s, t))
                + 0.5 * ((yudyuR[s][t] - yddydR[s][t]) * WC2R(CHq1R, p, r)
                -(yudyuI[s][t] - yddydI[s][t]) * WC2I(CHq1I, p, r))//RGE 2
                + TWO_THIRDS * g12 * YhYq * (WC2R(CHq1R, s, t) * delta[p][r] + WC2R(CHq1R, p, r) * delta[s][t])
                + 3. * g32 * WC6R(Cqq1R, p, t, s, r) + 9. * g32 * WC6R(Cqq3R, p, t, s, r)
                + 9. * g22 * WC6R(Cqq3R, p, r, s, t)
                -(6. / NC)*(g32 - 2. * NC * Yq2 * g12) * WC6R(Cqq1R, p, r, s, t)//RGE 3
                ;
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] += (gqR[p][v] * WC6R(Cqq1R, v, r, s, t) - gqI[p][v] * WC6I(Cqq1I, v, r, s, t))
                    +(gqR[s][v] * WC6R(Cqq1R, p, r, v, t) - gqI[s][v] * WC6I(Cqq1I, p, r, v, t))
                    +(WC6R(Cqq1R, p, v, s, t) * gqR[v][r] - WC6I(Cqq1I, p, v, s, t) * gqI[v][r])
                    +(WC6R(Cqq1R, p, r, s, v) * gqR[v][t] - WC6I(Cqq1I, p, r, s, v) * gqI[v][t])
                    //RGE 2
                    + FOUR_THIRDS * g12 * YlYq * (WC7R(Clq1R, w, w, s, t) * delta[p][r]
                    + WC7R(Clq1R, w, w, p, r) * delta[s][t])
                    + FOUR_THIRDS * g12 * Yq2 * NC * (
                    WC6R(Cqq1R, p, r, w, w) * delta[s][t] + WC6R(Cqq1R, s, t, w, w) * delta[p][r]
                    + WC6R(Cqq1R, w, w, s, t) * delta[p][r] + WC6R(Cqq1R, w, w, p, r) * delta[s][t])
                    + TWO_THIRDS * g12 * Yq2 * (
                    WC6R(Cqq1R, p, w, w, r) * delta[s][t] + WC6R(Cqq1R, s, w, w, t) * delta[p][r]
                    + WC6R(Cqq1R, w, r, p, w) * delta[s][t] + WC6R(Cqq1R, w, t, s, w) * delta[p][r])
                    + ONE_SIXTH * g32 * (
                    WC6R(Cqq1R, s, w, w, r) * delta[p][t] + WC6R(Cqq1R, p, w, w, t) * delta[r][s]
                    + WC6R(Cqq1R, w, r, s, w) * delta[p][t] + WC6R(Cqq1R, w, t, p, w) * delta[r][s])
                    - g32 * (ONE_THIRD / NC)*(
                    WC6R(Cqq1R, p, w, w, r) * delta[s][t] + WC6R(Cqq1R, s, w, w, t) * delta[p][r]
                    + WC6R(Cqq1R, w, r, p, w) * delta[s][t] + WC6R(Cqq1R, w, t, s, w) * delta[p][r])
                    + 2. * g12 * Yq2 * (
                    WC6R(Cqq3R, p, w, w, r) * delta[s][t] + WC6R(Cqq3R, s, w, w, t) * delta[p][r]
                    + WC6R(Cqq3R, w, r, p, w) * delta[s][t] + WC6R(Cqq3R, w, t, s, w) * delta[p][r])
                    + 0.5 * g32 * (
                    WC6R(Cqq3R, s, w, w, r) * delta[p][t] + WC6R(Cqq3R, p, w, w, t) * delta[r][s]
                    + WC6R(Cqq3R, w, r, s, w) * delta[p][t] + WC6R(Cqq3R, w, t, p, w) * delta[r][s])
                    -(1. / NC) * g32 * (
                    WC6R(Cqq3R, p, w, w, r) * delta[s][t] + WC6R(Cqq3R, s, w, w, t) * delta[p][r]
                    + WC6R(Cqq3R, w, r, p, w) * delta[s][t] + WC6R(Cqq3R, w, t, s, w) * delta[p][r])
                    + TWO_THIRDS * g12 * NC * YuYq *
                    (WC7R(Cqu1R, p, r, w, w) * delta[s][t] + WC7R(Cqu1R, s, t, w, w) * delta[p][r])
                    + TWO_THIRDS * g12 * NC * YdYq *
                    (WC7R(Cqd1R, p, r, w, w) * delta[s][t] + WC7R(Cqd1R, s, t, w, w) * delta[p][r])
                    +(g32 / 12.)*(
                    WC7R(Cqu8R, s, r, w, w) * delta[p][t] + WC7R(Cqu8R, p, t, w, w) * delta[r][s]
                    -(2. / NC) * WC7R(Cqu8R, p, r, w, w) * delta[s][t]-(2. / NC) * WC7R(Cqu8R, s, t, w, w) * delta[p][r])
                    +(g32 / 12.)*(
                    WC7R(Cqd8R, s, r, w, w) * delta[p][t] + WC7R(Cqd8R, p, t, w, w) * delta[r][s]
                    -(2. / NC) * WC7R(Cqd8R, p, r, w, w) * delta[s][t]-(2. / NC) * WC7R(Cqd8R, s, t, w, w) * delta[p][r])
                    + TWO_THIRDS * g12 * YeYq *
                    (WC7R(CqeR, p, r, w, w) * delta[s][t] + WC7R(CqeR, s, t, w, w) * delta[p][r])
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] += + (0.25 / NC)*(
                        (yudagR[p][v] * yuR[w][r] - yudagI[p][v] * yuI[w][r]) * WC7R(Cqu8R, s, t, v, w)
                        -(yudagR[p][v] * yuI[w][r] + yudagI[p][v] * yuR[w][r]) * WC7I(Cqu8I, s, t, v, w)
                        +(yudagR[s][v] * yuR[w][t] - yudagI[s][v] * yuI[w][t]) * WC7R(Cqu8R, p, r, v, w)
                        -(yudagR[s][v] * yuI[w][t] + yudagI[s][v] * yuR[w][t]) * WC7I(Cqu8I, p, r, v, w))
                        +(0.25 / NC)*(
                        (yddagR[p][v] * ydR[w][r] - yddagI[p][v] * ydI[w][r]) * WC7R(Cqd8R, s, t, v, w)
                        -(yddagR[p][v] * ydI[w][r] + yddagI[p][v] * ydR[w][r]) * WC7I(Cqd8I, s, t, v, w)
                        +(yddagR[s][v] * ydR[w][t] - yddagI[s][v] * ydI[w][t]) * WC7R(Cqd8R, p, r, v, w)
                        -(yddagR[s][v] * ydI[w][t] + yddagI[s][v] * ydR[w][t]) * WC7I(Cqd8I, p, r, v, w))

                        - 0.125 * (
                        (yudagR[p][v] * yuR[w][t] - yudagI[p][v] * yuI[w][t]) * WC7R(Cqu8R, s, r, v, w)
                        -(yudagR[p][v] * yuI[w][t] + yudagI[p][v] * yuR[w][t]) * WC7I(Cqu8I, s, r, v, w)
                        +(yudagR[s][v] * yuR[w][r] - yudagI[s][v] * yuI[w][r]) * WC7R(Cqu8R, p, t, v, w)
                        -(yudagR[s][v] * yuI[w][r] + yudagI[s][v] * yuR[w][r]) * WC7I(Cqu8I, p, t, v, w))
                        - 0.125 * (
                        (yddagR[p][v] * ydR[w][t] - yddagI[p][v] * ydI[w][t]) * WC7R(Cqd8R, s, r, v, w)
                        -(yddagR[p][v] * ydI[w][t] + yddagI[p][v] * ydR[w][t]) * WC7I(Cqd8I, s, r, v, w)
                        +(yddagR[s][v] * ydR[w][r] - yddagI[s][v] * ydI[w][r]) * WC7R(Cqd8R, p, t, v, w)
                        -(yddagR[s][v] * ydI[w][r] + yddagI[s][v] * ydR[w][r]) * WC7I(Cqd8I, p, t, v, w))

                        +(0.0625 / NC) *(
                        (ydR[w][t] * yuR[v][r] - ydI[w][t] * yuI[v][r]) * WC5(Cquqd8R, p, v, s, w)
                        -(ydI[w][t] * yuR[v][r] + ydR[w][t] * yuI[v][r]) * WC5(Cquqd8I, p, v, s, w)
                        +(ydR[w][r] * yuR[v][t] - ydI[w][r] * yuI[v][t]) * WC5(Cquqd8R, s, v, p, w)
                        -(ydI[w][r] * yuR[v][t] + ydR[w][r] * yuI[v][t]) * WC5(Cquqd8I, s, v, p, w))
                        +(0.0625 / NC) *(
                        (yddagR[s][w] * yudagR[p][v] - yddagI[s][w] * yudagI[p][v]) * WC5(Cquqd8R, r, v, t, w)
                        +(yddagI[s][w] * yudagR[p][v] + yddagR[s][w] * yudagI[p][v]) * WC5(Cquqd8I, r, v, t, w)
                        +(yddagR[p][w] * yudagR[s][v] - yddagI[p][w] * yudagI[s][v]) * WC5(Cquqd8R, t, v, r, w)
                        +(yddagI[p][w] * yudagR[s][v] + yddagR[p][w] * yudagI[s][v]) * WC5(Cquqd8I, t, v, r, w))

                        + 0.0625 * (
                        (ydR[w][t] * yuR[v][r] - ydI[w][t] * yuI[v][r]) * WC5(Cquqd8R, s, v, p, w)
                        -(ydI[w][t] * yuR[v][r] + ydR[w][t] * yuI[v][r]) * WC5(Cquqd8I, s, v, p, w)
                        +(ydR[w][r] * yuR[v][t] - ydI[w][r] * yuI[v][t]) * WC5(Cquqd8R, p, v, s, w)
                        -(ydI[w][r] * yuR[v][t] + ydR[w][r] * yuI[v][t]) * WC5(Cquqd8I, p, v, s, w))
                        + 0.0625 * (
                        (yddagR[s][w] * yudagR[p][v] - yddagI[s][w] * yudagI[p][v]) * WC5(Cquqd8R, t, v, r, w)
                        +(yddagI[s][w] * yudagR[p][v] + yddagR[s][w] * yudagI[p][v]) * WC5(Cquqd8I, t, v, r, w)
                        +(yddagR[p][w] * yudagR[s][v] - yddagI[p][w] * yudagI[s][v]) * WC5(Cquqd8R, r, v, t, w)
                        +(yddagI[p][w] * yudagR[s][v] + yddagR[p][w] * yudagI[s][v]) * WC5(Cquqd8I, r, v, t, w))

                        - 0.5 * ((yudagR[p][v] * yuR[w][r] - yudagI[p][v] * yuI[w][r]) * WC7R(Cqu1R, s, t, v, w)
                        -(yudagI[p][v] * yuR[w][r] + yudagR[p][v] * yuI[w][r]) * WC7I(Cqu1I, s, t, v, w))
                        - 0.5 * ((yddagR[p][v] * ydR[w][r] - yddagI[p][v] * ydI[w][r]) * WC7R(Cqd1R, s, t, v, w)
                        -(yddagI[s][v] * ydR[w][r] + yddagR[s][v] * ydI[w][r]) * WC7I(Cqd1I, s, t, v, w))
                        - 0.5 * ((yudagR[s][v] * yuR[w][t] - yudagI[s][v] * yuI[w][t]) * WC7R(Cqu1R, p, r, v, w)
                        -(yudagI[s][v] * yuR[w][t] + yudagR[s][v] * yuI[w][t]) * WC7I(Cqu1I, p, r, v, w))
                        - 0.5 * ((yddagR[s][v] * ydR[w][t] - yddagI[s][v] * ydI[w][t]) * WC7R(Cqd1R, p, r, v, w)
                        -(yddagI[s][v] * ydR[w][t] + yddagR[s][v] * ydI[w][t]) * WC7I(Cqd1I, p, r, v, w))

                        - 0.125 * ((ydR[w][t] * yuR[v][r] - ydI[w][t] * yuI[v][r]) * WC5(Cquqd1R, p, v, s, w)
                        -(ydI[w][t] * yuR[v][r] + ydR[w][t] * yuI[v][r]) * WC5(Cquqd1I, p, v, s, w))
                        - 0.125 * ((yddagR[s][w] * yudagR[p][v] - yddagI[s][w] * yudagI[p][v]) * WC5(Cquqd1R, r, v, t, w)
                        +(yddagI[s][w] * yudagR[p][v] + yddagR[s][w] * yudagI[p][v]) * WC5(Cquqd1I, r, v, t, w))
                        - 0.125 * ((ydR[w][r] * yuR[v][t] - ydI[w][r] * yuI[v][t]) * WC5(Cquqd1R, s, v, p, w)
                        -(ydI[w][r] * yuR[v][t] + ydR[w][r] * yuI[v][t]) * WC5(Cquqd1I, s, v, p, w))
                        - 0.125 * ((yddagR[p][w] * yudagR[s][v] - yddagI[p][w] * yudagI[s][v]) * WC5(Cquqd1R, t, v, r, w)
                        +(yddagI[p][w] * yudagR[s][v] + yddagR[p][w] * yudagI[s][v]) * WC5(Cquqd1I, t, v, r, w))
                        //RGE 2
                        ;
            }
        }

        f[c] *= loop_factor;
        c ++;
    }


    //Cqq1I
    for (d = 0; d < DWC6I; d ++) {
        p = WC6I_indices[d][0];
        r = WC6I_indices[d][1];
        s = WC6I_indices[d][2];
        t = WC6I_indices[d][3];
        f[c] = 0.5 * ((yudyuI[p][r] - yddydI[p][r]) * WC2R(CHq1R, s, t)
                +(yudyuR[p][r] - yddydR[p][r]) * WC2I(CHq1I, s, t))
                + 0.5 * ((yudyuI[s][t] - yddydI[s][t]) * WC2R(CHq1R, p, r)
                +(yudyuR[s][t] - yddydR[s][t]) * WC2I(CHq1I, p, r))//RGE 2
                + TWO_THIRDS * g12 * YhYq * (WC2I(CHq1I, s, t) * delta[p][r] + WC2I(CHq1I, p, r) * delta[s][t])
                + 3. * g32 * WC6I(Cqq1I, p, t, s, r) + 9. * g32 * WC6I(Cqq3I, p, t, s, r)
                + 9. * g22 * WC6I(Cqq3I, p, r, s, t)
                -(6. / NC)*(g32 - 2. * NC * Yq2 * g12) * WC6I(Cqq1I, p, r, s, t)//RGE 3
                ;
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] += (gqI[p][v] * WC6R(Cqq1R, v, r, s, t) + gqR[p][v] * WC6I(Cqq1I, v, r, s, t))
                    +(gqI[s][v] * WC6R(Cqq1R, p, r, v, t) + gqR[s][v] * WC6I(Cqq1I, p, r, v, t))
                    +(WC6R(Cqq1R, p, v, s, t) * gqI[v][r] + WC6I(Cqq1I, p, v, s, t) * gqR[v][r])
                    +(WC6R(Cqq1R, p, r, s, v) * gqI[v][t] + WC6I(Cqq1I, p, r, s, v) * gqR[v][t]) //RGE 2
                    //RGE 2
                    + FOUR_THIRDS * g12 * YlYq * (WC7I(Clq1I, w, w, s, t) * delta[p][r]
                    + WC7I(Clq1I, w, w, p, r) * delta[s][t])
                    + FOUR_THIRDS * g12 * Yq2 * NC * (
                    WC6I(Cqq1I, p, r, w, w) * delta[s][t] + WC6I(Cqq1I, s, t, w, w) * delta[p][r]
                    + WC6I(Cqq1I, w, w, s, t) * delta[p][r] + WC6I(Cqq1I, w, w, p, r) * delta[s][t])
                    + TWO_THIRDS * g12 * Yq2 * (
                    WC6I(Cqq1I, p, w, w, r) * delta[s][t] + WC6I(Cqq1I, s, w, w, t) * delta[p][r]
                    + WC6I(Cqq1I, w, r, p, w) * delta[s][t] + WC6I(Cqq1I, w, t, s, w) * delta[p][r])
                    + ONE_SIXTH * g32 * (
                    WC6I(Cqq1I, s, w, w, r) * delta[p][t] + WC6I(Cqq1I, p, w, w, t) * delta[r][s]
                    + WC6I(Cqq1I, w, r, s, w) * delta[p][t] + WC6I(Cqq1I, w, t, p, w) * delta[r][s])
                    - g32 * (ONE_THIRD / NC)*(
                    WC6I(Cqq1I, p, w, w, r) * delta[s][t] + WC6I(Cqq1I, s, w, w, t) * delta[p][r]
                    + WC6I(Cqq1I, w, r, p, w) * delta[s][t] + WC6I(Cqq1I, w, t, s, w) * delta[p][r])
                    + 2. * g12 * Yq2 * (
                    WC6I(Cqq3I, p, w, w, r) * delta[s][t] + WC6I(Cqq3I, s, w, w, t) * delta[p][r]
                    + WC6I(Cqq3I, w, r, p, w) * delta[s][t] + WC6I(Cqq3I, w, t, s, w) * delta[p][r])
                    + 0.5 * g32 * (
                    WC6I(Cqq3I, s, w, w, r) * delta[p][t] + WC6I(Cqq3I, p, w, w, t) * delta[r][s]
                    + WC6I(Cqq3I, w, r, s, w) * delta[p][t] + WC6I(Cqq3I, w, t, p, w) * delta[r][s])
                    -(1. / NC) * g32 * (
                    WC6I(Cqq3I, p, w, w, r) * delta[s][t] + WC6I(Cqq3I, s, w, w, t) * delta[p][r]
                    + WC6I(Cqq3I, w, r, p, w) * delta[s][t] + WC6I(Cqq3I, w, t, s, w) * delta[p][r])
                    + TWO_THIRDS * g12 * NC * YuYq *
                    (WC7I(Cqu1I, p, r, w, w) * delta[s][t] + WC7I(Cqu1I, s, t, w, w) * delta[p][r])
                    + TWO_THIRDS * g12 * NC * YdYq *
                    (WC7I(Cqd1I, p, r, w, w) * delta[s][t] + WC7I(Cqd1I, s, t, w, w) * delta[p][r])
                    +(g32 / 12.)*(
                    WC7I(Cqu8I, s, r, w, w) * delta[p][t] + WC7I(Cqu8I, p, t, w, w) * delta[r][s]
                    -(2. / NC) * WC7I(Cqu8I, p, r, w, w) * delta[s][t]-(2. / NC) * WC7I(Cqu8I, s, t, w, w) * delta[p][r])
                    +(g32 / 12.)*(
                    WC7I(Cqd8I, s, r, w, w) * delta[p][t] + WC7I(Cqd8I, p, t, w, w) * delta[r][s]
                    -(2. / NC) * WC7I(Cqd8I, p, r, w, w) * delta[s][t]-(2. / NC) * WC7I(Cqd8I, s, t, w, w) * delta[p][r])
                    + TWO_THIRDS * g12 * YeYq *
                    (WC7I(CqeI, p, r, w, w) * delta[s][t] + WC7I(CqeI, s, t, w, w) * delta[p][r])
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] += + (0.25 / NC)*(
                        (yudagR[p][v] * yuR[w][r] - yudagI[p][v] * yuI[w][r]) * WC7I(Cqu8I, s, t, v, w)
                        +(yudagR[p][v] * yuI[w][r] + yudagI[p][v] * yuR[w][r]) * WC7R(Cqu8R, s, t, v, w)
                        +(yudagR[s][v] * yuR[w][t] - yudagI[s][v] * yuI[w][t]) * WC7I(Cqu8I, p, r, v, w)
                        +(yudagR[s][v] * yuI[w][t] + yudagI[s][v] * yuR[w][t]) * WC7R(Cqu8R, p, r, v, w))
                        +(0.25 / NC)*(
                        (yddagR[p][v] * ydR[w][r] - yddagI[p][v] * ydI[w][r]) * WC7I(Cqd8I, s, t, v, w)
                        +(yddagR[p][v] * ydI[w][r] + yddagI[p][v] * ydR[w][r]) * WC7R(Cqd8R, s, t, v, w)
                        +(yddagR[s][v] * ydR[w][t] - yddagI[s][v] * ydI[w][t]) * WC7I(Cqd8I, p, r, v, w)
                        +(yddagR[s][v] * ydI[w][t] + yddagI[s][v] * ydR[w][t]) * WC7R(Cqd8R, p, r, v, w))

                        - 0.125 * (
                        (yudagR[p][v] * yuR[w][t] - yudagI[p][v] * yuI[w][t]) * WC7I(Cqu8I, s, r, v, w)
                        +(yudagR[p][v] * yuI[w][t] + yudagI[p][v] * yuR[w][t]) * WC7R(Cqu8R, s, r, v, w)
                        +(yudagR[s][v] * yuR[w][r] - yudagI[s][v] * yuI[w][r]) * WC7I(Cqu8I, p, t, v, w)
                        +(yudagR[s][v] * yuI[w][r] + yudagI[s][v] * yuR[w][r]) * WC7R(Cqu8R, p, t, v, w))
                        - 0.125 * (
                        (yddagR[p][v] * ydR[w][t] - yddagI[p][v] * ydI[w][t]) * WC7I(Cqd8I, s, r, v, w)
                        +(yddagR[p][v] * ydI[w][t] + yddagI[p][v] * ydR[w][t]) * WC7R(Cqd8R, s, r, v, w)
                        +(yddagR[s][v] * ydR[w][r] - yddagI[s][v] * ydI[w][r]) * WC7I(Cqd8I, p, t, v, w)
                        +(yddagR[s][v] * ydI[w][r] + yddagI[s][v] * ydR[w][r]) * WC7R(Cqd8R, p, t, v, w))

                        +(0.0625 / NC) *(
                        (ydR[w][t] * yuR[v][r] - ydI[w][t] * yuI[v][r]) * WC5(Cquqd8I, p, v, s, w)
                        +(ydI[w][t] * yuR[v][r] + ydR[w][t] * yuI[v][r]) * WC5(Cquqd8R, p, v, s, w)
                        +(ydR[w][r] * yuR[v][t] - ydI[w][r] * yuI[v][t]) * WC5(Cquqd8I, s, v, p, w)
                        +(ydI[w][r] * yuR[v][t] + ydR[w][r] * yuI[v][t]) * WC5(Cquqd8R, s, v, p, w))
                        +(0.0625 / NC) *(
                        - (yddagR[s][w] * yudagR[p][v] - yddagI[s][w] * yudagI[p][v]) * WC5(Cquqd8I, r, v, t, w)
                        +(yddagI[s][w] * yudagR[p][v] + yddagR[s][w] * yudagI[p][v]) * WC5(Cquqd8R, r, v, t, w)
                        -(yddagR[p][w] * yudagR[s][v] - yddagI[p][w] * yudagI[s][v]) * WC5(Cquqd8I, t, v, r, w)
                        +(yddagI[p][w] * yudagR[s][v] + yddagR[p][w] * yudagI[s][v]) * WC5(Cquqd8R, t, v, r, w))

                        + 0.0625 * (
                        (ydR[w][t] * yuR[v][r] - ydI[w][t] * yuI[v][r]) * WC5(Cquqd8I, s, v, p, w)
                        +(ydI[w][t] * yuR[v][r] + ydR[w][t] * yuI[v][r]) * WC5(Cquqd8R, s, v, p, w)
                        +(ydR[w][r] * yuR[v][t] - ydI[w][r] * yuI[v][t]) * WC5(Cquqd8I, p, v, s, w)
                        +(ydI[w][r] * yuR[v][t] + ydR[w][r] * yuI[v][t]) * WC5(Cquqd8R, p, v, s, w))
                        + 0.0625 * (
                        - (yddagR[s][w] * yudagR[p][v] - yddagI[s][w] * yudagI[p][v]) * WC5(Cquqd8I, t, v, r, w)
                        +(yddagI[s][w] * yudagR[p][v] + yddagR[s][w] * yudagI[p][v]) * WC5(Cquqd8R, t, v, r, w)
                        -(yddagR[p][w] * yudagR[s][v] - yddagI[p][w] * yudagI[s][v]) * WC5(Cquqd8I, r, v, t, w)
                        +(yddagI[p][w] * yudagR[s][v] + yddagR[p][w] * yudagI[s][v]) * WC5(Cquqd8R, r, v, t, w))

                        - 0.5 * ((yudagR[p][v] * yuR[w][r] - yudagI[p][v] * yuI[w][r]) * WC7I(Cqu1I, s, t, v, w)
                        +(yudagI[p][v] * yuR[w][r] + yudagR[p][v] * yuI[w][r]) * WC7R(Cqu1R, s, t, v, w))
                        - 0.5 * ((yddagR[p][v] * ydR[w][r] - yddagI[p][v] * ydI[w][r]) * WC7I(Cqd1I, s, t, v, w)
                        +(yddagI[s][v] * ydR[w][r] + yddagR[s][v] * ydI[w][r]) * WC7R(Cqd1R, s, t, v, w))
                        - 0.5 * ((yudagR[s][v] * yuR[w][t] - yudagI[s][v] * yuI[w][t]) * WC7I(Cqu1I, p, r, v, w)
                        +(yudagI[s][v] * yuR[w][t] + yudagR[s][v] * yuI[w][t]) * WC7R(Cqu1R, p, r, v, w))
                        - 0.5 * ((yddagR[s][v] * ydR[w][t] - yddagI[s][v] * ydI[w][t]) * WC7I(Cqd1I, p, r, v, w)
                        +(yddagI[s][v] * ydR[w][t] + yddagR[s][v] * ydI[w][t]) * WC7R(Cqd1R, p, r, v, w))

                        - 0.125 * ((ydR[w][t] * yuR[v][r] - ydI[w][t] * yuI[v][r]) * WC5(Cquqd1I, p, v, s, w)
                        +(ydI[w][t] * yuR[v][r] + ydR[w][t] * yuI[v][r]) * WC5(Cquqd1R, p, v, s, w))
                        - 0.125 * (- (yddagR[s][w] * yudagR[p][v] - yddagI[s][w] * yudagI[p][v]) * WC5(Cquqd1I, r, v, t, w)
                        +(yddagI[s][w] * yudagR[p][v] + yddagR[s][w] * yudagI[p][v]) * WC5(Cquqd1R, r, v, t, w))
                        - 0.125 * ((ydR[w][r] * yuR[v][t] - ydI[w][r] * yuI[v][t]) * WC5(Cquqd1I, s, v, p, w)
                        +(ydI[w][r] * yuR[v][t] + ydR[w][r] * yuI[v][t]) * WC5(Cquqd1R, s, v, p, w))
                        - 0.125 * (- (yddagR[p][w] * yudagR[s][v] - yddagI[p][w] * yudagI[s][v]) * WC5(Cquqd1I, t, v, r, w)
                        +(yddagI[p][w] * yudagR[s][v] + yddagR[p][w] * yudagI[s][v]) * WC5(Cquqd1R, t, v, r, w))
                        ;

            }
        }
        f[c] *= loop_factor;
        c ++;
    }
    //Cqq3R
    for (d = 0; d < DWC6R; d ++) {
        p = WC6R_indices[d][0];
        r = WC6R_indices[d][1];
        s = WC6R_indices[d][2];
        t = WC6R_indices[d][3];
        f[c] = - 0.5 * ((yudyuR[p][r] + yddydR[p][r]) * WC2R(CHq3R, s, t)
                - (yudyuI[p][r] + yddydI[p][r]) * WC2I(CHq3I, s, t))
                - 0.5 * ((yudyuR[s][t] + yddydR[s][t]) * WC2R(CHq3R, p, r)
                - (yudyuI[s][t] + yddydI[s][t]) * WC2I(CHq3I, p, r)) //RGE 2
                + ONE_SIXTH * g22 * (WC2R(CHq3R, s, t) * delta[p][r] + WC2R(CHq3R, p, r) * delta[s][t])
                - 3. * g32 * WC6R(Cqq3R, p, t, s, r)-(6. / NC) * g32 * WC6R(Cqq3R, p, r, s, t)
                - 6. * g22 * WC6R(Cqq3R, p, r, s, t) + 12. * Yq2 * g12 * WC6R(Cqq3R, p, r, s, t)
                + 3. * g32 * WC6R(Cqq1R, p, t, s, r) + 3. * g22 * WC6R(Cqq1R, p, r, s, t) //RGE3
                ;


        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index

            f[c] += (gqR[p][v] * WC6R(Cqq3R, v, r, s, t) - gqI[p][v] * WC6I(Cqq3I, v, r, s, t))
                    + (gqR[s][v] * WC6R(Cqq3R, p, r, v, t) - gqI[s][v] * WC6I(Cqq3I, p, r, v, t))
                    +(WC6R(Cqq3R, p, v, s, t) * gqR[v][r] - WC6I(Cqq3I, p, v, s, t) * gqI[v][r])
                    +(WC6R(Cqq3R, p, r, s, v) * gqR[v][t] - WC6I(Cqq3I, p, r, s, v) * gqI[v][t])
                    //RGE 2
                    + ONE_THIRD * g22 * (
                    WC7R(Clq3R, w, w, s, t) * delta[p][r] + WC7R(Clq3R, w, w, p, r) * delta[s][t])
                    + ONE_SIXTH * g22 * (
                    + WC6R(Cqq1R, s, w, w, t) * delta[p][r] + WC6R(Cqq1R, w, t, s, w) * delta[p][r]
                    + WC6R(Cqq1R, p, w, w, r) * delta[s][t] + WC6R(Cqq1R, w, r, p, w) * delta[s][t])
                    + ONE_SIXTH * g32 * (
                    WC6R(Cqq1R, p, w, w, t) * delta[r][s] + WC6R(Cqq1R, s, w, w, r) * delta[p][t]
                    + WC6R(Cqq1R, w, t, p, w) * delta[r][s] + WC6R(Cqq1R, w, r, s, w) * delta[p][t])
                    + ONE_THIRD * g22 * NC * (
                    WC6R(Cqq3R, p, r, w, w) * delta[s][t] + WC6R(Cqq3R, s, t, w, w) * delta[p][r]
                    + WC6R(Cqq3R, w, w, s, t) * delta[p][r] + WC6R(Cqq3R, w, w, p, r) * delta[s][t])
                    - ONE_SIXTH * g22 * (
                    WC6R(Cqq3R, p, w, w, r) * delta[s][t] + WC6R(Cqq3R, s, w, w, t) * delta[p][r]
                    + WC6R(Cqq3R, w, r, p, w) * delta[s][t] + WC6R(Cqq3R, w, t, s, w) * delta[p][r])
                    + 0.5 * g32 * (
                    WC6R(Cqq3R, p, w, w, t) * delta[r][s] + WC6R(Cqq3R, s, w, w, r) * delta[p][t]
                    + WC6R(Cqq3R, w, t, p, w) * delta[r][s] + WC6R(Cqq3R, w, r, s, w) * delta[p][t])
                    + (1. / 12.) * g32 * (
                    + WC7R(Cqu8R, p, t, w, w) * delta[r][s] + WC7R(Cqu8R, s, r, w, w) * delta[p][t]
                    + WC7R(Cqd8R, p, t, w, w) * delta[r][s] + WC7R(Cqd8R, s, r, w, w) * delta[p][t])
                    //RGE 3
                    ;
            for (w = 0; w < NG; w ++) {
                f[c] += - 0.125 * (
                        + (yudagR[p][v] * yuR[w][t] - yudagI[p][v] * yuI[w][t]) * WC7R(Cqu8R, s, r, v, w)
                        -(yudagI[p][v] * yuR[w][t] + yudagR[p][v] * yuI[w][t]) * WC7I(Cqu8I, s, r, v, w)
                        +(yudagR[s][v] * yuR[w][r] - yudagI[s][v] * yuI[w][r]) * WC7R(Cqu8R, p, t, v, w)
                        -(yudagI[s][v] * yuR[w][r] + yudagR[s][v] * yuI[w][r]) * WC7I(Cqu8I, p, t, v, w)
                        )
                        - 0.125 * (
                        + (yddagR[p][v] * ydR[w][t] - yddagI[p][v] * ydI[w][t]) * WC7R(Cqd8R, s, r, v, w)
                        -(yddagI[p][v] * ydR[w][t] + yddagR[p][v] * ydI[w][t]) * WC7I(Cqd8I, s, r, v, w)
                        +(yddagR[s][v] * ydR[w][r] - yddagI[s][v] * ydI[w][r]) * WC7R(Cqd8R, p, t, v, w)
                        -(yddagI[s][v] * ydR[w][r] + yddagR[s][v] * ydI[w][r]) * WC7I(Cqd8I, p, t, v, w)
                        )

                        -(0.0625 / NC)*(
                        + (ydR[w][t] * yuR[v][r] - ydI[w][t] * yuI[v][r]) * WC5(Cquqd8R, p, v, s, w)
                        -(ydI[w][t] * yuR[v][r] + ydR[w][t] * yuI[v][r]) * WC5(Cquqd8I, p, v, s, w)
                        +(ydR[w][r] * yuR[v][t] - ydI[w][r] * yuI[v][t]) * WC5(Cquqd8R, s, v, p, w)
                        -(ydI[w][r] * yuR[v][t] + ydR[w][r] * yuI[v][t]) * WC5(Cquqd8I, s, v, p, w)
                        )
                        -(0.0625 / NC)*(
                        + (yddagR[s][w] * yudagR[p][v] - yddagI[s][w] * yudagI[p][v]) * WC5(Cquqd8R, r, v, t, w)
                        +(yddagI[s][w] * yudagR[p][v] + yddagR[s][w] * yudagI[p][v]) * WC5(Cquqd8I, r, v, t, w)
                        +(yddagR[p][w] * yudagR[s][v] - yddagI[p][w] * yudagI[s][v]) * WC5(Cquqd8R, t, v, r, w)
                        +(yddagI[p][w] * yudagR[s][v] + yddagR[p][w] * yudagI[s][v]) * WC5(Cquqd8I, t, v, r, w)
                        )

                        - 0.0625 * (
                        + (ydR[w][t] * yuR[v][r] - ydI[w][t] * yuI[v][r]) * WC5(Cquqd8R, s, v, p, w)
                        -(ydI[w][t] * yuR[v][r] + ydR[w][t] * yuI[v][r]) * WC5(Cquqd8I, s, v, p, w)
                        +(ydR[w][r] * yuR[v][t] - ydI[w][r] * yuI[v][t]) * WC5(Cquqd8R, p, v, s, w)
                        -(ydI[w][r] * yuR[v][t] + ydR[w][r] * yuI[v][t]) * WC5(Cquqd8I, p, v, s, w)
                        )
                        - 0.0625 * (
                        + (yddagR[s][w] * yudagR[p][v] - yddagI[s][w] * yudagI[p][v]) * WC5(Cquqd8R, t, v, r, w)
                        +(yddagI[s][w] * yudagR[p][v] + yddagR[s][w] * yudagI[p][v]) * WC5(Cquqd8I, t, v, r, w)
                        +(yddagR[p][w] * yudagR[s][v] - yddagI[p][w] * yudagI[s][v]) * WC5(Cquqd8R, r, v, t, w)
                        +(yddagI[p][w] * yudagR[s][v] + yddagR[p][w] * yudagI[s][v]) * WC5(Cquqd8I, r, v, t, w)
                        )

                        + 0.125 * (
                        + (ydR[w][t] * yuR[v][r] - ydI[w][t] * yuI[v][r]) * WC5(Cquqd1R, p, v, s, w)
                        -(ydI[w][t] * yuR[v][r] + ydR[w][t] * yuI[v][r]) * WC5(Cquqd1I, p, v, s, w)
                        )
                        + 0.125 * (
                        + (yddagR[s][w] * yudagR[p][v] - yddagI[s][w] * yudagI[p][v]) * WC5(Cquqd1R, r, v, t, w)
                        +(yddagI[s][w] * yudagR[p][v] + yddagR[s][w] * yudagI[p][v]) * WC5(Cquqd1I, r, v, t, w)
                        )
                        + 0.125 * (
                        + (ydR[w][r] * yuR[v][t] - ydI[w][r] * yuI[v][t]) * WC5(Cquqd1R, s, v, p, w)
                        -(ydI[w][r] * yuR[v][t] + ydR[w][r] * yuI[v][t]) * WC5(Cquqd1I, s, v, p, w)
                        )
                        + 0.125 * (
                        + (yddagR[p][w] * yudagR[s][v] - yddagI[p][w] * yudagI[s][v]) * WC5(Cquqd1R, t, v, r, w)
                        +(yddagI[p][w] * yudagR[s][v] + yddagR[p][w] * yudagI[s][v]) * WC5(Cquqd1I, t, v, r, w)
                        )
                        //RGE 2
                        ;
            }
        }
        f[c] *= loop_factor;
        c ++;


    }
    //Cqq3I
    for (d = 0; d < DWC6I; d ++) {
        p = WC6I_indices[d][0];
        r = WC6I_indices[d][1];
        s = WC6I_indices[d][2];
        t = WC6I_indices[d][3];
        f[c] = - 0.5 * ((yudyuI[p][r] + yddydI[p][r]) * WC2R(CHq3R, s, t)
                + (yudyuR[p][r] + yddydR[p][r]) * WC2I(CHq3I, s, t))
                - 0.5 * ((yudyuI[s][t] + yddydI[s][t]) * WC2R(CHq3R, p, r)
                + (yudyuR[s][t] + yddydR[s][t]) * WC2I(CHq3I, p, r)) //RGE 2
                + ONE_SIXTH * g22 * (WC2I(CHq3I, s, t) * delta[p][r] + WC2I(CHq3I, p, r) * delta[s][t])
                - 3. * g32 * WC6I(Cqq3I, p, t, s, r)-(6. / NC) * g32 * WC6I(Cqq3I, p, r, s, t)
                - 6. * g22 * WC6I(Cqq3I, p, r, s, t) + 12. * Yq2 * g12 * WC6I(Cqq3I, p, r, s, t)
                + 3. * g32 * WC6I(Cqq1I, p, t, s, r) + 3. * g22 * WC6I(Cqq1I, p, r, s, t) //RGE3
                ;
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] += (gqI[p][v] * WC6R(Cqq3R, v, r, s, t) + gqR[p][v] * WC6I(Cqq3I, v, r, s, t))
                    + (gqI[s][v] * WC6R(Cqq3R, p, r, v, t) + gqR[s][v] * WC6I(Cqq3I, p, r, v, t))
                    +(WC6R(Cqq3R, p, v, s, t) * gqI[v][r] + WC6I(Cqq3I, p, v, s, t) * gqR[v][r])
                    +(WC6R(Cqq3R, p, r, s, v) * gqI[v][t] + WC6I(Cqq3I, p, r, s, v) * gqR[v][t])
                    //RGE 2
                    + ONE_THIRD * g22 * (
                    WC7I(Clq3I, w, w, s, t) * delta[p][r] + WC7I(Clq3I, w, w, p, r) * delta[s][t])
                    + ONE_SIXTH * g22 * (
                    + WC6I(Cqq1I, s, w, w, t) * delta[p][r] + WC6I(Cqq1I, w, t, s, w) * delta[p][r]
                    + WC6I(Cqq1I, p, w, w, r) * delta[s][t] + WC6I(Cqq1I, w, r, p, w) * delta[s][t])
                    + ONE_SIXTH * g32 * (
                    WC6I(Cqq1I, p, w, w, t) * delta[r][s] + WC6I(Cqq1I, s, w, w, r) * delta[p][t]
                    + WC6I(Cqq1I, w, t, p, w) * delta[r][s] + WC6I(Cqq1I, w, r, s, w) * delta[p][t])
                    + ONE_THIRD * g22 * NC * (
                    WC6I(Cqq3I, p, r, w, w) * delta[s][t] + WC6I(Cqq3I, s, t, w, w) * delta[p][r]
                    + WC6I(Cqq3I, w, w, s, t) * delta[p][r] + WC6I(Cqq3I, w, w, p, r) * delta[s][t])
                    - ONE_SIXTH * g22 * (
                    WC6I(Cqq3I, p, w, w, r) * delta[s][t] + WC6I(Cqq3I, s, w, w, t) * delta[p][r]
                    + WC6I(Cqq3I, w, r, p, w) * delta[s][t] + WC6I(Cqq3I, w, t, s, w) * delta[p][r])
                    + 0.5 * g32 * (
                    WC6I(Cqq3I, p, w, w, t) * delta[r][s] + WC6I(Cqq3I, s, w, w, r) * delta[p][t]
                    + WC6I(Cqq3I, w, t, p, w) * delta[r][s] + WC6I(Cqq3I, w, r, s, w) * delta[p][t])
                    + (1. / 12.) * g32 * (
                    + WC7I(Cqu8I, p, t, w, w) * delta[r][s] + WC7I(Cqu8I, s, r, w, w) * delta[p][t]
                    + WC7I(Cqd8I, p, t, w, w) * delta[r][s] + WC7I(Cqd8I, s, r, w, w) * delta[p][t])
                    //RGE 3
                    ;

            for (w = 0; w < NG; w ++) {
                f[c] += - 0.125 * (
                        + (yudagR[p][v] * yuR[w][t] - yudagI[p][v] * yuI[w][t]) * WC7I(Cqu8I, s, r, v, w)
                        +(yudagI[p][v] * yuR[w][t] + yudagR[p][v] * yuI[w][t]) * WC7R(Cqu8R, s, r, v, w)
                        +(yudagR[s][v] * yuR[w][r] - yudagI[s][v] * yuI[w][r]) * WC7I(Cqu8I, p, t, v, w)
                        +(yudagI[s][v] * yuR[w][r] + yudagR[s][v] * yuI[w][r]) * WC7R(Cqu8R, p, t, v, w)
                        )
                        - 0.125 * (
                        + (yddagR[p][v] * ydR[w][t] - yddagI[p][v] * ydI[w][t]) * WC7I(Cqd8I, s, r, v, w)
                        +(yddagI[p][v] * ydR[w][t] + yddagR[p][v] * ydI[w][t]) * WC7R(Cqd8R, s, r, v, w)
                        +(yddagR[s][v] * ydR[w][r] - yddagI[s][v] * ydI[w][r]) * WC7I(Cqd8I, p, t, v, w)
                        + (yddagI[s][v] * ydR[w][r] + yddagR[s][v] * ydI[w][r]) * WC7R(Cqd8R, p, t, v, w)
                        )

                        -(0.0625 / NC)*(
                        + (ydR[w][t] * yuR[v][r] - ydI[w][t] * yuI[v][r]) * WC5(Cquqd8I, p, v, s, w)
                        +(ydI[w][t] * yuR[v][r] + ydR[w][t] * yuI[v][r]) * WC5(Cquqd8R, p, v, s, w)
                        +(ydR[w][r] * yuR[v][t] - ydI[w][r] * yuI[v][t]) * WC5(Cquqd8I, s, v, p, w)
                        +(ydI[w][r] * yuR[v][t] + ydR[w][r] * yuI[v][t]) * WC5(Cquqd8R, s, v, p, w)
                        )
                        -(0.0625 / NC)*(
                        - (yddagR[s][w] * yudagR[p][v] - yddagI[s][w] * yudagI[p][v]) * WC5(Cquqd8I, r, v, t, w)
                        +(yddagI[s][w] * yudagR[p][v] + yddagR[s][w] * yudagI[p][v]) * WC5(Cquqd8R, r, v, t, w)
                        -(yddagR[p][w] * yudagR[s][v] - yddagI[p][w] * yudagI[s][v]) * WC5(Cquqd8I, t, v, r, w)
                        +(yddagI[p][w] * yudagR[s][v] + yddagR[p][w] * yudagI[s][v]) * WC5(Cquqd8R, t, v, r, w)
                        )

                        - 0.0625 * (
                        + (ydR[w][t] * yuR[v][r] - ydI[w][t] * yuI[v][r]) * WC5(Cquqd8I, s, v, p, w)
                        +(ydI[w][t] * yuR[v][r] + ydR[w][t] * yuI[v][r]) * WC5(Cquqd8R, s, v, p, w)
                        +(ydR[w][r] * yuR[v][t] - ydI[w][r] * yuI[v][t]) * WC5(Cquqd8I, p, v, s, w)
                        +(ydI[w][r] * yuR[v][t] + ydR[w][r] * yuI[v][t]) * WC5(Cquqd8R, p, v, s, w)
                        )
                        - 0.0625 * (
                        - (yddagR[s][w] * yudagR[p][v] - yddagI[s][w] * yudagI[p][v]) * WC5(Cquqd8I, t, v, r, w)
                        +(yddagI[s][w] * yudagR[p][v] + yddagR[s][w] * yudagI[p][v]) * WC5(Cquqd8R, t, v, r, w)
                        -(yddagR[p][w] * yudagR[s][v] - yddagI[p][w] * yudagI[s][v]) * WC5(Cquqd8I, r, v, t, w)
                        +(yddagI[p][w] * yudagR[s][v] + yddagR[p][w] * yudagI[s][v]) * WC5(Cquqd8R, r, v, t, w)
                        )

                        + 0.125 * (
                        + (ydR[w][t] * yuR[v][r] - ydI[w][t] * yuI[v][r]) * WC5(Cquqd1I, p, v, s, w)
                        +(ydI[w][t] * yuR[v][r] + ydR[w][t] * yuI[v][r]) * WC5(Cquqd1R, p, v, s, w)
                        )
                        + 0.125 * (
                        - (yddagR[s][w] * yudagR[p][v] - yddagI[s][w] * yudagI[p][v]) * WC5(Cquqd1I, r, v, t, w)
                        +(yddagI[s][w] * yudagR[p][v] + yddagR[s][w] * yudagI[p][v]) * WC5(Cquqd1R, r, v, t, w)
                        )
                        + 0.125 * (
                        + (ydR[w][r] * yuR[v][t] - ydI[w][r] * yuI[v][t]) * WC5(Cquqd1I, s, v, p, w)
                        +(ydI[w][r] * yuR[v][t] + ydR[w][r] * yuI[v][t]) * WC5(Cquqd1R, s, v, p, w)
                        )
                        + 0.125 * (
                        - (yddagR[p][w] * yudagR[s][v] - yddagI[p][w] * yudagI[s][v]) * WC5(Cquqd1I, t, v, r, w)
                        +(yddagI[p][w] * yudagR[s][v] + yddagR[p][w] * yudagI[s][v]) * WC5(Cquqd1R, t, v, r, w)
                        )
                        //RGE 2
                        ;
            }
        }
        f[c] *= loop_factor;
        c ++;
    }

    //Clq1R
    for (d = 0; d < DWC7R; d ++) {
        p = WC7R_indices[d][0];
        r = WC7R_indices[d][1];
        s = WC7R_indices[d][2];
        t = WC7R_indices[d][3];
        f[c] =
                - (yedyeR[p][r] * WC2R(CHq1R, s, t) - yedyeI[p][r] * WC2I(CHq1I, s, t))
                +(
                + (yudyuR[s][t] - yddydR[s][t]) * WC2R(CHl1R, p, r)
                -(yudyuI[s][t] - yddydI[s][t]) * WC2I(CHl1I, p, r))
                //RGE 2
                + FOUR_THIRDS * g12 * YhYl * WC2R(CHq1R, s, t) * delta[p][r]
                + FOUR_THIRDS * g12 * YhYq * WC2R(CHl1R, p, r) * delta[s][t]
                + 12. * YlYq * g12 * WC7R(Clq1R, p, r, s, t)
                + 9. * g22 * WC7R(Clq3R, p, r, s, t)
                ; //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] +=
                    + (glR[p][v] * WC7R(Clq1R, v, r, s, t) - glI[p][v] * WC7I(Clq1I, v, r, s, t))
                    +(gqR[s][v] * WC7R(Clq1R, p, r, v, t) - gqI[s][v] * WC7I(Clq1I, p, r, v, t))
                    +(WC7R(Clq1R, p, v, s, t) * glR[v][r] - WC7I(Clq1I, p, v, s, t) * glI[v][r])
                    +(WC7R(Clq1R, p, r, s, v) * gqR[v][t] - WC7I(Clq1I, p, r, s, v) * gqI[v][t])
                    //RGE 2
                    + EIGHT_THIRDS * g12 * YlYq * delta[s][t]*(
                    + WC6R(CllR, p, r, w, w) + WC6R(CllR, w, w, p, r)
                    )
                    + FOUR_THIRDS * g12 * YlYq * delta[s][t]*(
                    + WC6R(CllR, p, w, w, r) + WC6R(CllR, w, r, p, w)
                    )
                    + EIGHT_THIRDS * g12 * NC * Yq2 * WC7R(Clq1R, p, r, w, w) * delta[s][t]
                    + EIGHT_THIRDS * g12 * Yl2 * WC7R(Clq1R, w, w, s, t) * delta[p][r]
                    + EIGHT_THIRDS * g12 * NC * YlYq * delta[p][r]*(
                    + WC6R(Cqq1R, s, t, w, w) + WC6R(Cqq1R, w, w, s, t)
                    )
                    + FOUR_THIRDS * g12 * YlYq * delta[p][r]*(
                    + WC6R(Cqq1R, s, w, w, t) + WC6R(Cqq1R, w, t, s, w)
                    )
                    + 4. * g12 * YlYq * delta[p][r]* (
                    + WC6R(Cqq3R, s, w, w, t) + WC6R(Cqq3R, w, t, s, w)
                    )
                    + FOUR_THIRDS * g12 * NC * YuYl * WC7R(Cqu1R, s, t, w, w) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YdYl * WC7R(Cqd1R, s, t, w, w) * delta[p][r]
                    + FOUR_THIRDS * g12 * YeYl * WC7R(CqeR, s, t, w, w) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YuYq * WC7R(CluR, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * NC * YdYq * WC7R(CldR, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * YeYq * WC7R(CleR, p, r, w, w) * delta[s][t]

                    //RGE 3
                    ;
            for (w = 0; w < NG; w ++) {
                f[c] += 0.25 * (
                        + (yuR[w][t] * yeR[v][r] - yuI[w][t] * yeI[v][r]) * WC5(Clequ1R, p, v, s, w)
                        -(yuI[w][t] * yeR[v][r] + yuR[w][t] * yeI[v][r]) * WC5(Clequ1I, p, v, s, w)
                        )
                        + 0.25 * (
                        + (yudagR[s][w] * yedagR[p][v] - yudagI[s][w] * yedagI[p][v]) * WC5(Clequ1R, r, v, t, w)
                        +(yudagI[s][w] * yedagR[p][v] + yudagR[s][w] * yedagI[p][v]) * WC5(Clequ1I, r, v, t, w)
                        )
                        -(
                        + (yudagR[s][v] * yuR[w][t] - yudagI[s][v] * yuI[w][t]) * WC7R(CluR, p, r, v, w)
                        -(yudagI[s][v] * yuR[w][t] + yudagR[s][v] * yuI[w][t]) * WC7I(CluI, p, r, v, w)
                        )
                        -(
                        + (yddagR[s][v] * ydR[w][t] - yddagI[s][v] * ydI[w][t]) * WC7R(CldR, p, r, v, w)
                        -(yddagI[s][v] * ydR[w][t] + yddagR[s][v] * ydI[w][t]) * WC7I(CldI, p, r, v, w)
                        )
                        -(
                        + (yeR[p][v] * yeR[w][r] - yeI[p][v] * yeI[w][r]) * WC7R(CqeR, s, t, v, w)
                        -(yeI[p][v] * yeR[w][r] + yeR[p][v] * yeI[w][r]) * WC7I(CqeI, s, t, v, w)
                        )
                        + 0.25 * (
                        + (yddagR[s][w] * yeR[v][r] - yddagI[s][w] * yeI[v][r]) * WC5(CledqR, p, v, w, t)
                        -(yddagI[s][w] * yeR[v][r] + yddagR[s][w] * yeI[v][r]) * WC5(CledqI, p, v, w, t)
                        +(yedagR[p][v] * ydR[w][t] - yedagI[p][v] * ydI[w][t]) * WC5(CledqR, r, v, w, s)
                        +(yedagI[p][v] * ydR[w][t] + yedagR[p][v] * ydI[w][t]) * WC5(CledqI, r, v, w, s)
                        )
                        - 3. * (
                        + (yeR[v][r] * yuR[w][t] - yeI[v][r] * yuI[w][t]) * WC5(Clequ3R, p, v, s, w)
                        -(yeI[v][r] * yuR[w][t] + yeR[v][r] * yuI[w][t]) * WC5(Clequ3I, p, v, s, w)
                        +(yedagR[p][v] * yudagR[s][w] - yedagI[p][v] * yudagI[s][w]) * WC5(Clequ3R, r, v, t, w)
                        +(yedagI[p][v] * yudagR[s][w] + yedagR[p][v] * yudagI[s][w]) * WC5(Clequ3I, r, v, t, w)
                        )
                        ; //RGE 2
            }
        }
        f[c] *= loop_factor;
        c ++;
    }


    //Clq1I
    for (d = 0; d < DWC7I; d ++) {
        p = WC7I_indices[d][0];
        r = WC7I_indices[d][1];
        s = WC7I_indices[d][2];
        t = WC7I_indices[d][3];
        f[c] =
                - (yedyeI[p][r] * WC2R(CHq1R, s, t) + yedyeR[p][r] * WC2I(CHq1I, s, t))
                +(
                + (yudyuI[s][t] - yddydI[s][t]) * WC2R(CHl1R, p, r)
                +(yudyuR[s][t] - yddydR[s][t]) * WC2I(CHl1I, p, r))
                //RGE 2
                + FOUR_THIRDS * g12 * YhYl * WC2I(CHq1I, s, t) * delta[p][r]
                + FOUR_THIRDS * g12 * YhYq * WC2I(CHl1I, p, r) * delta[s][t]
                + 12. * YlYq * g12 * WC7I(Clq1I, p, r, s, t)
                + 9. * g22 * WC7I(Clq3I, p, r, s, t)
                ; //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] +=
                    + (glI[p][v] * WC7R(Clq1R, v, r, s, t) + glR[p][v] * WC7I(Clq1I, v, r, s, t))
                    +(gqI[s][v] * WC7R(Clq1R, p, r, v, t) + gqR[s][v] * WC7I(Clq1I, p, r, v, t))
                    +(WC7R(Clq1R, p, v, s, t) * glI[v][r] + WC7I(Clq1I, p, v, s, t) * glR[v][r])
                    +(WC7R(Clq1R, p, r, s, v) * gqI[v][t] + WC7I(Clq1I, p, r, s, v) * gqR[v][t])
                    //RGE 2
                    + EIGHT_THIRDS * g12 * YlYq * delta[s][t]*(
                    + WC6I(CllI, p, r, w, w) + WC6I(CllI, w, w, p, r)
                    )
                    + FOUR_THIRDS * g12 * YlYq * delta[s][t]*(
                    + WC6I(CllI, p, w, w, r) + WC6I(CllI, w, r, p, w)
                    )
                    + EIGHT_THIRDS * g12 * NC * Yq2 * WC7I(Clq1I, p, r, w, w) * delta[s][t]
                    + EIGHT_THIRDS * g12 * Yl2 * WC7I(Clq1I, w, w, s, t) * delta[p][r]
                    + EIGHT_THIRDS * g12 * NC * YlYq * delta[p][r]*(
                    + WC6I(Cqq1I, s, t, w, w) + WC6I(Cqq1I, w, w, s, t)
                    )
                    + FOUR_THIRDS * g12 * YlYq * delta[p][r]*(
                    + WC6I(Cqq1I, s, w, w, t) + WC6I(Cqq1I, w, t, s, w)
                    )
                    + 4. * g12 * YlYq * delta[p][r]* (
                    + WC6I(Cqq3I, s, w, w, t) + WC6I(Cqq3I, w, t, s, w)
                    )
                    + FOUR_THIRDS * g12 * NC * YuYl * WC7I(Cqu1I, s, t, w, w) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YdYl * WC7I(Cqd1I, s, t, w, w) * delta[p][r]
                    + FOUR_THIRDS * g12 * YeYl * WC7I(CqeI, s, t, w, w) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YuYq * WC7I(CluI, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * NC * YdYq * WC7I(CldI, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * YeYq * WC7I(CleI, p, r, w, w) * delta[s][t]

                    //RGE 3
                    ;
            for (w = 0; w < NG; w ++) {
                f[c] += 0.25 * (
                        + (yuR[w][t] * yeR[v][r] - yuI[w][t] * yeI[v][r]) * WC5(Clequ1I, p, v, s, w)
                        +(yuI[w][t] * yeR[v][r] + yuR[w][t] * yeI[v][r]) * WC5(Clequ1R, p, v, s, w)
                        )
                        + 0.25 * (
                        - (yudagR[s][w] * yedagR[p][v] - yudagI[s][w] * yedagI[p][v]) * WC5(Clequ1I, r, v, t, w)
                        +(yudagI[s][w] * yedagR[p][v] + yudagR[s][w] * yedagI[p][v]) * WC5(Clequ1R, r, v, t, w)
                        )
                        -(
                        + (yudagR[s][v] * yuR[w][t] - yudagI[s][v] * yuI[w][t]) * WC7I(CluI, p, r, v, w)
                        +(yudagI[s][v] * yuR[w][t] + yudagR[s][v] * yuI[w][t]) * WC7R(CluR, p, r, v, w)
                        )
                        -(
                        + (yddagR[s][v] * ydR[w][t] - yddagI[s][v] * ydI[w][t]) * WC7I(CldI, p, r, v, w)
                        +(yddagI[s][v] * ydR[w][t] + yddagR[s][v] * ydI[w][t]) * WC7R(CldR, p, r, v, w)
                        )
                        -(
                        + (yeR[p][v] * yeR[w][r] - yeI[p][v] * yeI[w][r]) * WC7I(CqeI, s, t, v, w)
                        +(yeI[p][v] * yeR[w][r] + yeR[p][v] * yeI[w][r]) * WC7R(CqeR, s, t, v, w)
                        )
                        + 0.25 * (
                        + (yddagR[s][w] * yeR[v][r] - yddagI[s][w] * yeI[v][r]) * WC5(CledqI, p, v, w, t)
                        +(yddagI[s][w] * yeR[v][r] + yddagR[s][w] * yeI[v][r]) * WC5(CledqR, p, v, w, t)
                        -(yedagR[p][v] * ydR[w][t] - yedagI[p][v] * ydI[w][t]) * WC5(CledqI, r, v, w, s)
                        +(yedagI[p][v] * ydR[w][t] + yedagR[p][v] * ydI[w][t]) * WC5(CledqR, r, v, w, s)
                        )
                        - 3. * (
                        + (yeR[v][r] * yuR[w][t] - yeI[v][r] * yuI[w][t]) * WC5(Clequ3I, p, v, s, w)
                        +(yeI[v][r] * yuR[w][t] + yeR[v][r] * yuI[w][t]) * WC5(Clequ3R, p, v, s, w)
                        -(yedagR[p][v] * yudagR[s][w] - yedagI[p][v] * yudagI[s][w]) * WC5(Clequ3I, r, v, t, w)
                        +(yedagI[p][v] * yudagR[s][w] + yedagR[p][v] * yudagI[s][w]) * WC5(Clequ3R, r, v, t, w)
                        )
                        ; //RGE 2
            }
        }
        f[c] *= loop_factor;
        c ++;
    }



    //Clq3R
    for (d = 0; d < DWC7R; d ++) {
        p = WC7R_indices[d][0];
        r = WC7R_indices[d][1];
        s = WC7R_indices[d][2];
        t = WC7R_indices[d][3];
        f[c] =
                - (yedyeR[p][r] * WC2R(CHq3R, s, t) - yedyeI[p][r] * WC2I(CHq3I, s, t))
                -(
                + (yudyuR[s][t] + yddydR[s][t]) * WC2R(CHl3R, p, r)
                -(yudyuI[s][t] + yddydI[s][t]) * WC2I(CHl3I, p, r)
                )
                //RGE 2
                + ONE_THIRD * g22 * WC2R(CHq3R, s, t) * delta[p][r]
                + ONE_THIRD * g22 * WC2R(CHl3R, p, r) * delta[s][t]
                + 3. * g22 * WC7R(Clq1R, p, r, s, t)
                - 6. * (g22 - 2. * YlYq * g12) * WC7R(Clq3R, p, r, s, t)
                ; //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] +=
                    + (glR[p][v] * WC7R(Clq3R, v, r, s, t) - glI[p][v] * WC7I(Clq3I, v, r, s, t))
                    +(gqR[s][v] * WC7R(Clq3R, p, r, v, t) - gqI[s][v] * WC7I(Clq3I, p, r, v, t))
                    +(WC7R(Clq3R, p, v, s, t) * glR[v][r] - WC7I(Clq3I, p, v, s, t) * glI[v][r])
                    +(WC7R(Clq3R, p, r, s, v) * gqR[v][t] - WC7I(Clq3I, p, r, s, v) * gqI[v][t])
                    //RGE 2
                    + TWO_THIRDS * g22 * NC * WC7R(Clq3R, p, r, w, w) * delta[s][t]
                    + TWO_THIRDS * g22 * WC7R(Clq3R, w, w, s, t) * delta[p][r]
                    + ONE_THIRD * g22 * delta[p][r] * (
                    + WC6R(Cqq1R, s, w, w, t) + WC6R(Cqq1R, w, t, s, w)
                    )
                    + TWO_THIRDS * g22 * NC * delta[p][r]*(
                    + WC6R(Cqq3R, s, t, w, w) + WC6R(Cqq3R, w, w, s, t)
                    )
                    - ONE_THIRD * g22 * delta[p][r]*(
                    + WC6R(Cqq3R, s, w, w, t) + WC6R(Cqq3R, w, t, s, w)
                    )
                    + ONE_THIRD * g22 * delta[s][t]*(
                    + WC6R(CllR, p, w, w, r) + WC6R(CllR, w, r, p, w)
                    )
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] +=
                        - 0.25 * (
                        + (yuR[w][t] * yeR[v][r] - yuI[w][t] * yeI[v][r]) * WC5(Clequ1R, p, v, s, w)
                        -(yuI[w][t] * yeR[v][r] + yuR[w][t] * yeI[v][r]) * WC5(Clequ1I, p, v, s, w)
                        )
                        - 0.25 * (
                        + (yudagR[s][w] * yedagR[p][v] - yudagI[s][w] * yedagI[p][v]) * WC5(Clequ1R, r, v, t, w)
                        +(yudagI[s][w] * yedagR[p][v] + yudagR[s][w] * yedagI[p][v]) * WC5(Clequ1I, r, v, t, w)
                        )
                        + 0.25 * (
                        + (yddagR[s][w] * yeR[v][r] - yddagI[s][w] * yeI[v][r]) * WC5(CledqR, p, v, w, t)
                        -(yddagI[s][w] * yeR[v][r] + yddagR[s][w] * yeI[v][r]) * WC5(CledqI, p, v, w, t)
                        +(yedagR[p][v] * ydR[w][t] - yedagI[p][v] * ydI[w][t]) * WC5(CledqR, r, v, w, s)
                        +(yedagI[p][v] * ydR[w][t] + yedagR[p][v] * ydI[w][t]) * WC5(CledqI, r, v, w, s)
                        )
                        + 3. * (
                        + (yeR[v][r] * yuR[w][t] - yeI[v][r] * yuI[w][t]) * WC5(Clequ3R, p, v, s, w)
                        -(yeI[v][r] * yuR[w][t] + yeR[v][r] * yuI[w][t]) * WC5(Clequ3I, p, v, s, w)
                        +(yedagR[p][v] * yudagR[s][w] - yedagI[p][v] * yudagI[s][w]) * WC5(Clequ3R, r, v, t, w)
                        +(yedagI[p][v] * yudagR[s][w] + yedagR[p][v] * yudagI[s][w]) * WC5(Clequ3I, r, v, t, w)
                        )
                        ; //RGE 2
            }
        }
        f[c] *= loop_factor;
        c ++;
    }

    //Clq3I
    for (d = 0; d < DWC7I; d ++) {
        p = WC7I_indices[d][0];
        r = WC7I_indices[d][1];
        s = WC7I_indices[d][2];
        t = WC7I_indices[d][3];
        f[c] =
                - (yedyeI[p][r] * WC2R(CHq3R, s, t) + yedyeR[p][r] * WC2I(CHq3I, s, t))
                -(
                + (yudyuI[s][t] + yddydI[s][t]) * WC2R(CHl3R, p, r)
                +(yudyuR[s][t] + yddydR[s][t]) * WC2I(CHl3I, p, r)
                )
                //RGE 2
                + ONE_THIRD * g22 * WC2I(CHq3I, s, t) * delta[p][r]
                + ONE_THIRD * g22 * WC2I(CHl3I, p, r) * delta[s][t]
                + 3. * g22 * WC7I(Clq1I, p, r, s, t)
                - 6. * (g22 - 2. * YlYq * g12) * WC7I(Clq3I, p, r, s, t)
                ; //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] +=
                    + (glI[p][v] * WC7R(Clq3R, v, r, s, t) + glR[p][v] * WC7I(Clq3I, v, r, s, t))
                    +(gqI[s][v] * WC7R(Clq3R, p, r, v, t) + gqR[s][v] * WC7I(Clq3I, p, r, v, t))
                    +(WC7R(Clq3R, p, v, s, t) * glI[v][r] + WC7I(Clq3I, p, v, s, t) * glR[v][r])
                    +(WC7R(Clq3R, p, r, s, v) * gqI[v][t] + WC7I(Clq3I, p, r, s, v) * gqR[v][t])
                    //RGE 2
                    + TWO_THIRDS * g22 * NC * WC7I(Clq3I, p, r, w, w) * delta[s][t]
                    + TWO_THIRDS * g22 * WC7I(Clq3I, w, w, s, t) * delta[p][r]
                    + ONE_THIRD * g22 * delta[p][r] * (
                    + WC6I(Cqq1I, s, w, w, t) + WC6I(Cqq1I, w, t, s, w)
                    )
                    + TWO_THIRDS * g22 * NC * delta[p][r]*(
                    + WC6I(Cqq3I, s, t, w, w) + WC6I(Cqq3I, w, w, s, t)
                    )
                    - ONE_THIRD * g22 * delta[p][r]*(
                    + WC6I(Cqq3I, s, w, w, t) + WC6I(Cqq3I, w, t, s, w)
                    )
                    + ONE_THIRD * g22 * delta[s][t]*(
                    + WC6I(CllI, p, w, w, r) + WC6I(CllI, w, r, p, w)
                    )
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] +=
                        - 0.25 * (
                        + (yuR[w][t] * yeR[v][r] - yuI[w][t] * yeI[v][r]) * WC5(Clequ1I, p, v, s, w)
                        +(yuI[w][t] * yeR[v][r] + yuR[w][t] * yeI[v][r]) * WC5(Clequ1R, p, v, s, w)
                        )
                        - 0.25 * (
                        - (yudagR[s][w] * yedagR[p][v] - yudagI[s][w] * yedagI[p][v]) * WC5(Clequ1I, r, v, t, w)
                        +(yudagI[s][w] * yedagR[p][v] + yudagR[s][w] * yedagI[p][v]) * WC5(Clequ1R, r, v, t, w)
                        )
                        + 0.25 * (
                        + (yddagR[s][w] * yeR[v][r] - yddagI[s][w] * yeI[v][r]) * WC5(CledqI, p, v, w, t)
                        +(yddagI[s][w] * yeR[v][r] + yddagR[s][w] * yeI[v][r]) * WC5(CledqR, p, v, w, t)
                        -(yedagR[p][v] * ydR[w][t] - yedagI[p][v] * ydI[w][t]) * WC5(CledqI, r, v, w, s)
                        +(yedagI[p][v] * ydR[w][t] + yedagR[p][v] * ydI[w][t]) * WC5(CledqR, r, v, w, s)
                        )
                        + 3. * (
                        + (yeR[v][r] * yuR[w][t] - yeI[v][r] * yuI[w][t]) * WC5(Clequ3I, p, v, s, w)
                        +(yeI[v][r] * yuR[w][t] + yeR[v][r] * yuI[w][t]) * WC5(Clequ3R, p, v, s, w)
                        -(yedagR[p][v] * yudagR[s][w] - yedagI[p][v] * yudagI[s][w]) * WC5(Clequ3I, r, v, t, w)
                        +(yedagI[p][v] * yudagR[s][w] + yedagR[p][v] * yudagI[s][w]) * WC5(Clequ3R, r, v, t, w)
                        )
                        ; //RGE 2
            }
        }
        f[c] *= loop_factor;
        c ++;
    }




    //----------------------RGE SMEFT class 8_RRRR----------------------
    //CeeR
    for (d = 0; d < DWC8R; d ++) {
        p = WC8R_indices[d][0];
        r = WC8R_indices[d][1];
        s = WC8R_indices[d][2];
        t = WC8R_indices[d][3];
        f[c] = (geR[p][r] * WC2R(CHeR, s, t) - geI[p][r] * WC2I(CHeI, s, t))
                +(geR[s][t] * WC2R(CHeR, p, r) - geI[s][t] * WC2I(CHeI, p, r)) //RGE 2
                + TWO_THIRDS * g12 * YhYe * WC2R(CHeR, s, t) * delta[p][r]
                + TWO_THIRDS * g12 * YhYe * WC2R(CHeR, p, r) * delta[s][t]
                + 12. * Ye2 * g12 * WC8R(CeeR, p, r, s, t); //RGE 3

        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index

            f[c] += + (geR[p][v] * WC8R(CeeR, v, r, s, t) - geI[p][v] * WC8I(CeeI, v, r, s, t))
                    +(geR[s][v] * WC8R(CeeR, p, r, v, t) - geI[s][v] * WC8I(CeeI, p, r, v, t))
                    +(WC8R(CeeR, p, v, s, t) * geR[v][r] - WC8I(CeeI, p, v, s, t) * geI[v][r])
                    +(WC8R(CeeR, p, r, s, v) * geR[v][t] - WC8I(CeeI, p, r, s, v) * geI[v][t])
                    //RGE 2
                    + FOUR_THIRDS * g12 * YeYl * (
                    WC7R(CleR, w, w, p, r) * delta[s][t] + WC7R(CleR, w, w, s, t) * delta[p][r])
                    + FOUR_THIRDS * g12 * NC * YeYq * (
                    WC7R(CqeR, w, w, p, r) * delta[s][t] + WC7R(CqeR, w, w, s, t) * delta[p][r])
                    + TWO_THIRDS * g12 * NC * YuYe * (
                    WC7R(CeuR, p, r, w, w) * delta[s][t] + WC7R(CeuR, s, t, w, w) * delta[p][r])
                    + TWO_THIRDS * g12 * NC * YdYe * (
                    WC7R(CedR, p, r, w, w) * delta[s][t] + WC7R(CedR, s, t, w, w) * delta[p][r])

                    + TWO_THIRDS * g12 * Ye2 * (
                    + WC8R(CeeR, p, r, w, w) * delta[s][t] + WC8R(CeeR, s, t, w, w) * delta[p][r]
                    + WC8R(CeeR, w, w, p, r) * delta[s][t] + WC8R(CeeR, w, w, s, t) * delta[p][r]
                    + WC8R(CeeR, p, w, w, r) * delta[s][t] + WC8R(CeeR, s, w, w, t) * delta[p][r]
                    + WC8R(CeeR, w, t, s, w) * delta[p][r] + WC8R(CeeR, w, r, p, w) * delta[s][t])
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] += - (
                        (yedagR[w][r] * yeR[p][v] - yedagI[w][r] * yeI[p][v]) * WC7R(CleR, v, w, s, t)
                        -(yedagI[w][r] * yeR[p][v] + yedagR[w][r] * yeI[p][v]) * WC7I(CleI, v, w, s, t)
                        )
                        -(
                        (yedagR[w][t] * yeR[s][v] - yedagI[w][t] * yeI[s][v]) * WC7R(CleR, v, w, p, r)
                        -(yedagI[w][t] * yeR[s][v] + yedagR[w][t] * yeI[s][v]) * WC7I(CleI, v, w, p, r)
                        )
                        //RGE 2 
                        ;
            }
        }
        f[c] *= loop_factor;
        c ++;
    }
    //CeeI
    for (d = 0; d < DWC8I; d ++) {
        p = WC8I_indices[d][0];
        r = WC8I_indices[d][1];
        s = WC8I_indices[d][2];
        t = WC8I_indices[d][3];

        f[c] = (geI[p][r] * WC2R(CHeR, s, t) + geR[p][r] * WC2I(CHeI, s, t))
                +(geI[s][t] * WC2R(CHeR, p, r) + geR[s][t] * WC2I(CHeI, p, r)) //RGE 2
                + TWO_THIRDS * g12 * YhYe * WC2I(CHeI, s, t) * delta[p][r]
                + TWO_THIRDS * g12 * YhYe * WC2I(CHeI, p, r) * delta[s][t]
                + 12. * Ye2 * g12 * WC8I(CeeI, p, r, s, t); //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index

            f[c] += + (geI[p][v] * WC8R(CeeR, v, r, s, t) + geR[p][v] * WC8I(CeeI, v, r, s, t))
                    +(geI[s][v] * WC8R(CeeR, p, r, v, t) + geR[s][v] * WC8I(CeeI, p, r, v, t))
                    +(WC8R(CeeR, p, v, s, t) * geI[v][r] + WC8I(CeeI, p, v, s, t) * geR[v][r])
                    +(WC8R(CeeR, p, r, s, v) * geI[v][t] + WC8I(CeeI, p, r, s, v) * geR[v][t])
                    //RGE 2
                    + FOUR_THIRDS * g12 * YeYl * (
                    WC7I(CleI, w, w, p, r) * delta[s][t] + WC7I(CleI, w, w, s, t) * delta[p][r])
                    + FOUR_THIRDS * g12 * NC * YeYq * (
                    WC7I(CqeI, w, w, p, r) * delta[s][t] + WC7I(CqeI, w, w, s, t) * delta[p][r])
                    + TWO_THIRDS * g12 * NC * YuYe * (
                    WC7I(CeuI, p, r, w, w) * delta[s][t] + WC7I(CeuI, s, t, w, w) * delta[p][r])
                    + TWO_THIRDS * g12 * NC * YdYe * (
                    WC7I(CedI, p, r, w, w) * delta[s][t] + WC7I(CedI, s, t, w, w) * delta[p][r])

                    + TWO_THIRDS * g12 * Ye2 * (
                    + WC8I(CeeI, p, r, w, w) * delta[s][t] + WC8I(CeeI, s, t, w, w) * delta[p][r]
                    + WC8I(CeeI, w, w, p, r) * delta[s][t] + WC8I(CeeI, w, w, s, t) * delta[p][r]
                    + WC8I(CeeI, p, w, w, r) * delta[s][t] + WC8I(CeeI, s, w, w, t) * delta[p][r]
                    + WC8I(CeeI, w, t, s, w) * delta[p][r] + WC8I(CeeI, w, r, p, w) * delta[s][t])
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] += - (
                        (yedagR[w][r] * yeR[p][v] - yedagI[w][r] * yeI[p][v]) * WC7I(CleI, v, w, s, t)
                        +(yedagI[w][r] * yeR[p][v] + yedagR[w][r] * yeI[p][v]) * WC7R(CleR, v, w, s, t)
                        )
                        -(
                        (yedagR[w][t] * yeR[s][v] - yedagI[w][t] * yeI[s][v]) * WC7I(CleI, v, w, p, r)
                        +(yedagI[w][t] * yeR[s][v] + yedagR[w][t] * yeI[s][v]) * WC7R(CleR, v, w, p, r)
                        )
                        //RGE 2 
                        ;
            }

        }
        f[c] *= loop_factor;
        c ++;
    }
    //CuuR
    for (d = 0; d < DWC6R; d ++) {
        p = WC6R_indices[d][0];
        r = WC6R_indices[d][1];
        s = WC6R_indices[d][2];
        t = WC6R_indices[d][3];
        f[c] = - (guR[p][r] * WC2R(CHuR, s, t) - guI[p][r] * WC2I(CHuI, s, t))
                -(guR[s][t] * WC2R(CHuR, p, r) - guI[s][t] * WC2I(CHuI, p, r))
                //RGE 2
                + TWO_THIRDS * g12 * YhYu *
                (WC2R(CHuR, s, t) * delta[p][r] + WC2R(CHuR, p, r) * delta[s][t])
                + 6. * g32 * WC6R(CuuR, p, t, s, r) - 6. * (g32 / NC - 2. * Yu2 * g12) * WC6R(CuuR, p, r, s, t)
                ; //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] += + (guR[p][v] * WC6R(CuuR, v, r, s, t) - guI[p][v] * WC6I(CuuI, v, r, s, t))
                    +(guR[s][v] * WC6R(CuuR, p, r, v, t) - guI[s][v] * WC6I(CuuI, p, r, v, t))
                    +(WC6R(CuuR, p, v, s, t) * guR[v][r] - WC6I(CuuI, p, v, s, t) * guI[v][r])
                    +(WC6R(CuuR, p, r, s, v) * guR[v][t] - WC6I(CuuI, p, r, s, v) * guI[v][t])
                    //RGE 2
                    + TWO_THIRDS * g12 * YuYe * (
                    WC7R(CeuR, w, w, s, t) * delta[p][r] + WC7R(CeuR, w, w, p, r) * delta[s][t])
                    + FOUR_THIRDS * g12 * YuYl * (
                    WC7R(CluR, w, w, p, r) * delta[s][t] + WC7R(CluR, w, w, s, t) * delta[p][r])
                    + FOUR_THIRDS * g12 * NC * YuYq * (
                    WC7R(Cqu1R, w, w, s, t) * delta[p][r] + WC7R(Cqu1R, w, w, p, r) * delta[s][t])
                    + ONE_THIRD * g32 * (
                    WC7R(Cqu8R, w, w, p, t) * delta[r][s] + WC7R(Cqu8R, w, w, s, r) * delta[p][t])
                    -(ONE_THIRD / NC) * g32 * (
                    WC7R(Cqu8R, w, w, s, t) * delta[p][r] + WC7R(Cqu8R, w, w, p, r) * delta[s][t])
                    + TWO_THIRDS * g12 * NC * Yu2 * (
                    WC6R(CuuR, p, r, w, w) * delta[s][t] + WC6R(CuuR, s, t, w, w) * delta[p][r]
                    + WC6R(CuuR, w, w, p, r) * delta[s][t] + WC6R(CuuR, w, w, s, t) * delta[p][r])
                    + TWO_THIRDS * g12 * Yu2 * (
                    WC6R(CuuR, p, w, w, r) * delta[s][t] + WC6R(CuuR, s, w, w, t) * delta[p][r]
                    + WC6R(CuuR, w, r, p, w) * delta[s][t] + WC6R(CuuR, w, t, s, w) * delta[p][r])
                    + ONE_THIRD * g32 * (
                    WC6R(CuuR, p, w, w, t) * delta[r][s] + WC6R(CuuR, s, w, w, r) * delta[p][t]
                    + WC6R(CuuR, w, t, p, w) * delta[r][s] + WC6R(CuuR, w, r, s, w) * delta[p][t])
                    -(ONE_THIRD / NC) * g32 * (
                    WC6R(CuuR, p, w, w, r) * delta[s][t] + WC6R(CuuR, s, w, w, t) * delta[p][r]
                    + WC6R(CuuR, w, r, p, w) * delta[s][t] + WC6R(CuuR, w, t, s, w) * delta[p][r])
                    + TWO_THIRDS * g12 * NC * YuYd * (
                    WC7R(Cud1R, p, r, w, w) * delta[s][t] + WC7R(Cud1R, s, t, w, w) * delta[p][r])
                    + ONE_SIXTH * g32 * (
                    WC7R(Cud8R, p, t, w, w) * delta[r][s] + WC7R(Cud8R, s, r, w, w) * delta[p][t])
                    -(ONE_SIXTH / NC) * g32 *
                    (WC7R(Cud8R, p, r, w, w) * delta[s][t] + WC7R(Cud8R, s, t, w, w) * delta[p][r])
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] += - (
                        + (yudagR[w][r] * yuR[p][v] - yudagI[w][r] * yuI[p][v]) * WC7R(Cqu1R, v, w, s, t)
                        -(yudagI[w][r] * yuR[p][v] + yudagR[w][r] * yuI[p][v]) * WC7I(Cqu1I, v, w, s, t))
                        -(
                        + (yudagR[w][t] * yuR[s][v] - yudagI[w][t] * yuI[s][v]) * WC7R(Cqu1R, v, w, p, r)
                        -(yudagI[w][t] * yuR[s][v] + yudagR[w][t] * yuI[s][v]) * WC7I(Cqu1I, v, w, p, r))

                        +(0.5 / NC)*(
                        + (yuR[p][v] * yudagR[w][r] - yuI[p][v] * yudagI[w][r]) * WC7R(Cqu8R, v, w, s, t)
                        -(yuI[p][v] * yudagR[w][r] + yuR[p][v] * yudagI[w][r]) * WC7I(Cqu8I, v, w, s, t))
                        +(0.5 / NC)*(
                        + (yuR[s][v] * yudagR[w][t] - yuI[s][v] * yudagI[w][t]) * WC7R(Cqu8R, v, w, p, r)
                        -(yuI[s][v] * yudagR[w][t] + yuR[s][v] * yudagI[w][t]) * WC7I(Cqu8I, v, w, p, r))

                        - 0.5 * (
                        + (yudagR[w][r] * yuR[s][v] - yudagI[w][r] * yuI[s][v]) * WC7R(Cqu8R, v, w, p, t)
                        -(yudagI[w][r] * yuR[s][v] + yudagR[w][r] * yuI[s][v]) * WC7I(Cqu8I, v, w, p, t))
                        - 0.5 * (
                        + (yudagR[w][t] * yuR[p][v] - yudagI[w][t] * yuI[p][v]) * WC7R(Cqu8R, v, w, s, r)
                        -(yudagI[w][t] * yuR[p][v] + yudagR[w][t] * yuI[p][v]) * WC7I(Cqu8I, v, w, s, r))

                        //RGE 2
                        ;
            }
        }

        f[c] *= loop_factor;
        c ++;
    }
    //CuuI
    for (d = 0; d < DWC6I; d ++) {
        p = WC6I_indices[d][0];
        r = WC6I_indices[d][1];
        s = WC6I_indices[d][2];
        t = WC6I_indices[d][3];
        f[c] = - (guI[p][r] * WC2R(CHuR, s, t) + guR[p][r] * WC2I(CHuI, s, t))
                -(guI[s][t] * WC2R(CHuR, p, r) + guR[s][t] * WC2I(CHuI, p, r))
                //RGE 2
                + TWO_THIRDS * g12 * YhYu *
                (WC2I(CHuI, s, t) * delta[p][r] + WC2I(CHuI, p, r) * delta[s][t])
                + 6. * g32 * WC6I(CuuI, p, t, s, r) - 6. * (g32 / NC - 2. * Yu2 * g12) * WC6I(CuuI, p, r, s, t)
                ; //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] += + (guI[p][v] * WC6R(CuuR, v, r, s, t) + guR[p][v] * WC6I(CuuI, v, r, s, t))
                    +(guI[s][v] * WC6R(CuuR, p, r, v, t) + guR[s][v] * WC6I(CuuI, p, r, v, t))
                    +(WC6R(CuuR, p, v, s, t) * guI[v][r] + WC6I(CuuI, p, v, s, t) * guR[v][r])
                    +(WC6R(CuuR, p, r, s, v) * guI[v][t] + WC6I(CuuI, p, r, s, v) * guR[v][t])
                    //RGE 2
                    + TWO_THIRDS * g12 * YuYe * (
                    WC7I(CeuI, w, w, s, t) * delta[p][r] + WC7I(CeuI, w, w, p, r) * delta[s][t])
                    + FOUR_THIRDS * g12 * YuYl * (
                    WC7I(CluI, w, w, p, r) * delta[s][t] + WC7I(CluI, w, w, s, t) * delta[p][r])
                    + FOUR_THIRDS * g12 * NC * YuYq * (
                    WC7I(Cqu1I, w, w, s, t) * delta[p][r] + WC7I(Cqu1I, w, w, p, r) * delta[s][t])
                    + ONE_THIRD * g32 * (
                    WC7I(Cqu8I, w, w, p, t) * delta[r][s] + WC7I(Cqu8I, w, w, s, r) * delta[p][t])
                    -(ONE_THIRD / NC) * g32 * (
                    WC7I(Cqu8I, w, w, s, t) * delta[p][r] + WC7I(Cqu8I, w, w, p, r) * delta[s][t])
                    + TWO_THIRDS * g12 * NC * Yu2 * (
                    WC6I(CuuI, p, r, w, w) * delta[s][t] + WC6I(CuuI, s, t, w, w) * delta[p][r]
                    + WC6I(CuuI, w, w, p, r) * delta[s][t] + WC6I(CuuI, w, w, s, t) * delta[p][r])
                    + TWO_THIRDS * g12 * Yu2 * (
                    WC6I(CuuI, p, w, w, r) * delta[s][t] + WC6I(CuuI, s, w, w, t) * delta[p][r]
                    + WC6I(CuuI, w, r, p, w) * delta[s][t] + WC6I(CuuI, w, t, s, w) * delta[p][r])
                    + ONE_THIRD * g32 * (
                    WC6I(CuuI, p, w, w, t) * delta[r][s] + WC6I(CuuI, s, w, w, r) * delta[p][t]
                    + WC6I(CuuI, w, t, p, w) * delta[r][s] + WC6I(CuuI, w, r, s, w) * delta[p][t])
                    -(ONE_THIRD / NC) * g32 * (
                    WC6I(CuuI, p, w, w, r) * delta[s][t] + WC6I(CuuI, s, w, w, t) * delta[p][r]
                    + WC6I(CuuI, w, r, p, w) * delta[s][t] + WC6I(CuuI, w, t, s, w) * delta[p][r])
                    + TWO_THIRDS * g12 * NC * YuYd * (
                    WC7I(Cud1I, p, r, w, w) * delta[s][t] + WC7I(Cud1I, s, t, w, w) * delta[p][r])
                    + ONE_SIXTH * g32 * (
                    WC7I(Cud8I, p, t, w, w) * delta[r][s] + WC7I(Cud8I, s, r, w, w) * delta[p][t])
                    -(ONE_SIXTH / NC) * g32 *
                    (WC7I(Cud8I, p, r, w, w) * delta[s][t] + WC7I(Cud8I, s, t, w, w) * delta[p][r])
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] += - (
                        + (yudagR[w][r] * yuR[p][v] - yudagI[w][r] * yuI[p][v]) * WC7I(Cqu1I, v, w, s, t)
                        +(yudagI[w][r] * yuR[p][v] + yudagR[w][r] * yuI[p][v]) * WC7R(Cqu1R, v, w, s, t))
                        -(
                        + (yudagR[w][t] * yuR[s][v] - yudagI[w][t] * yuI[s][v]) * WC7I(Cqu1I, v, w, p, r)
                        +(yudagI[w][t] * yuR[s][v] + yudagR[w][t] * yuI[s][v]) * WC7R(Cqu1R, v, w, p, r))

                        +(0.5 / NC)*(
                        + (yuR[p][v] * yudagR[w][r] - yuI[p][v] * yudagI[w][r]) * WC7I(Cqu8I, v, w, s, t)
                        +(yuI[p][v] * yudagR[w][r] + yuR[p][v] * yudagI[w][r]) * WC7R(Cqu8R, v, w, s, t))
                        +(0.5 / NC)*(
                        + (yuR[s][v] * yudagR[w][t] - yuI[s][v] * yudagI[w][t]) * WC7I(Cqu8I, v, w, p, r)
                        +(yuI[s][v] * yudagR[w][t] + yuR[s][v] * yudagI[w][t]) * WC7R(Cqu8R, v, w, p, r))

                        - 0.5 * (
                        + (yudagR[w][r] * yuR[s][v] - yudagI[w][r] * yuI[s][v]) * WC7I(Cqu8I, v, w, p, t)
                        +(yudagI[w][r] * yuR[s][v] + yudagR[w][r] * yuI[s][v]) * WC7R(Cqu8R, v, w, p, t))
                        - 0.5 * (
                        + (yudagR[w][t] * yuR[p][v] - yudagI[w][t] * yuI[p][v]) * WC7I(Cqu8I, v, w, s, r)
                        +(yudagI[w][t] * yuR[p][v] + yudagR[w][t] * yuI[p][v]) * WC7R(Cqu8R, v, w, s, r))

                        //RGE 2
                        ;
            }
        }
        f[c] *= loop_factor;
        c ++;
    }
    //CddR
    for (d = 0; d < DWC6R; d ++) {
        p = WC6R_indices[d][0];
        r = WC6R_indices[d][1];
        s = WC6R_indices[d][2];
        t = WC6R_indices[d][3];
        f[c] = + (gdR[p][r] * WC2R(CHdR, s, t) - gdI[p][r] * WC2I(CHdI, s, t))
                +(gdR[s][t] * WC2R(CHdR, p, r) - gdI[s][t] * WC2I(CHdI, p, r))
                //RGE 2
                + TWO_THIRDS * g12 * YhYd *
                (WC2R(CHdR, s, t) * delta[p][r] + WC2R(CHdR, p, r) * delta[s][t])
                + 6. * g32 * WC6R(CddR, p, t, s, r) - 6. * (g32 / NC - 2. * Yd2 * g12) * WC6R(CddR, p, r, s, t)
                ; //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] += + (gdR[p][v] * WC6R(CddR, v, r, s, t) - gdI[p][v] * WC6I(CddI, v, r, s, t))
                    +(gdR[s][v] * WC6R(CddR, p, r, v, t) - gdI[s][v] * WC6I(CddI, p, r, v, t))
                    +(WC6R(CddR, p, v, s, t) * gdR[v][r] - WC6I(CddI, p, v, s, t) * gdI[v][r])
                    +(WC6R(CddR, p, r, s, v) * gdR[v][t] - WC6I(CddI, p, r, s, v) * gdI[v][t])
                    //RGE 2
                    + TWO_THIRDS * g12 * NC * Yd2 * (
                    WC6R(CddR, p, r, w, w) * delta[s][t] + WC6R(CddR, s, t, w, w) * delta[p][r]
                    + WC6R(CddR, w, w, p, r) * delta[s][t] + WC6R(CddR, w, w, s, t) * delta[p][r])
                    + TWO_THIRDS * g12 * Yd2 * (
                    WC6R(CddR, p, w, w, r) * delta[s][t] + WC6R(CddR, s, w, w, t) * delta[p][r]
                    + WC6R(CddR, w, t, s, w) * delta[p][r] + WC6R(CddR, w, r, p, w) * delta[s][t])
                    + ONE_THIRD * g32 * (
                    WC6R(CddR, p, w, w, t) * delta[r][s] + WC6R(CddR, s, w, w, r) * delta[p][t]
                    + WC6R(CddR, w, t, p, w) * delta[r][s] + WC6R(CddR, w, r, s, w) * delta[p][t])
                    -(ONE_THIRD / NC) * g32 * (
                    WC6R(CddR, p, w, w, r) * delta[s][t] + WC6R(CddR, s, w, w, t) * delta[p][r]
                    + WC6R(CddR, w, t, s, w) * delta[p][r] + WC6R(CddR, w, r, p, w) * delta[s][t])

                    + FOUR_THIRDS * g12 * YdYl * (
                    WC7R(CldR, w, w, p, r) * delta[s][t] + WC7R(CldR, w, w, s, t) * delta[p][r])
                    + FOUR_THIRDS * g12 * NC * YdYq * (
                    WC7R(Cqd1R, w, w, p, r) * delta[s][t] + WC7R(Cqd1R, w, w, s, t) * delta[p][r])
                    + ONE_THIRD * g32 * (
                    WC7R(Cqd8R, w, w, s, r) * delta[p][t] + WC7R(Cqd8R, w, w, p, t) * delta[r][s])
                    -(ONE_THIRD / NC) * g32 * (
                    WC7R(Cqd8R, w, w, p, r) * delta[s][t] + WC7R(Cqd8R, w, w, s, t) * delta[p][r])
                    + TWO_THIRDS * g12 * YdYe * (
                    WC7R(CedR, w, w, p, r) * delta[s][t] + WC7R(CedR, w, w, s, t) * delta[p][r])
                    + TWO_THIRDS * g12 * NC * YuYd * (
                    WC7R(Cud1R, w, w, p, r) * delta[s][t] + WC7R(Cud1R, w, w, s, t) * delta[p][r])
                    + ONE_SIXTH * g32 * (
                    WC7R(Cud8R, w, w, p, t) * delta[r][s] + WC7R(Cud8R, w, w, s, r) * delta[p][t])
                    -(ONE_SIXTH / NC) * g32 * (
                    WC7R(Cud8R, w, w, p, r) * delta[s][t] + WC7R(Cud8R, w, w, s, t) * delta[p][r])
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] += - (
                        + (yddagR[w][r] * ydR[p][v] - yddagI[w][r] * ydI[p][v]) * WC7R(Cqd1R, v, w, s, t)
                        -(yddagI[w][r] * ydR[p][v] + yddagR[w][r] * ydI[p][v]) * WC7I(Cqd1I, v, w, s, t))
                        -(
                        + (yddagR[w][t] * ydR[s][v] - yddagI[w][t] * ydI[s][v]) * WC7R(Cqd1R, v, w, p, r)
                        -(yddagI[w][t] * ydR[s][v] + yddagR[w][t] * ydI[s][v]) * WC7I(Cqd1I, v, w, p, r))

                        +(0.5 / NC)*(
                        + (ydR[p][v] * yddagR[w][r] - ydI[p][v] * yddagI[w][r]) * WC7R(Cqd8R, v, w, s, t)
                        -(ydI[p][v] * yddagR[w][r] + ydR[p][v] * yddagI[w][r]) * WC7I(Cqd8I, v, w, s, t))
                        +(0.5 / NC)*(
                        + (ydR[s][v] * yddagR[w][t] - ydI[s][v] * yddagI[w][t]) * WC7R(Cqd8R, v, w, p, r)
                        -(ydI[s][v] * yddagR[w][t] + ydR[s][v] * yddagI[w][t]) * WC7I(Cqd8I, v, w, p, r))

                        - 0.5 * (
                        + (yddagR[w][r] * ydR[s][v] - yddagI[w][r] * ydI[s][v]) * WC7R(Cqd8R, v, w, p, t)
                        -(yddagI[w][r] * ydR[s][v] + yddagR[w][r] * ydI[s][v]) * WC7I(Cqd8I, v, w, p, t))
                        - 0.5 * (
                        + (yddagR[w][t] * ydR[p][v] - yddagI[w][t] * ydI[p][v]) * WC7R(Cqd8R, v, w, s, r)
                        -(yddagI[w][t] * ydR[p][v] + yddagR[w][t] * ydI[p][v]) * WC7I(Cqd8I, v, w, s, r))

                        //RGE 2
                        ;
            }
        }

        f[c] *= loop_factor;
        c ++;
    }
    //CddI
    for (d = 0; d < DWC6I; d ++) {
        p = WC6I_indices[d][0];
        r = WC6I_indices[d][1];
        s = WC6I_indices[d][2];
        t = WC6I_indices[d][3];
        f[c] = + (gdI[p][r] * WC2R(CHdR, s, t) + gdR[p][r] * WC2I(CHdI, s, t))
                +(gdI[s][t] * WC2R(CHdR, p, r) + gdR[s][t] * WC2I(CHdI, p, r))
                //RGE 2
                + TWO_THIRDS * g12 * YhYd *
                (WC2I(CHdI, s, t) * delta[p][r] + WC2I(CHdI, p, r) * delta[s][t])
                + 6. * g32 * WC6I(CddI, p, t, s, r) - 6. * (g32 / NC - 2. * Yd2 * g12) * WC6I(CddI, p, r, s, t)
                ; //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] += + (gdI[p][v] * WC6R(CddR, v, r, s, t) + gdR[p][v] * WC6I(CddI, v, r, s, t))
                    +(gdI[s][v] * WC6R(CddR, p, r, v, t) + gdR[s][v] * WC6I(CddI, p, r, v, t))
                    +(WC6R(CddR, p, v, s, t) * gdI[v][r] + WC6I(CddI, p, v, s, t) * gdR[v][r])
                    +(WC6R(CddR, p, r, s, v) * gdI[v][t] + WC6I(CddI, p, r, s, v) * gdR[v][t])
                    //RGE 2
                    + TWO_THIRDS * g12 * NC * Yd2 * (
                    WC6I(CddI, p, r, w, w) * delta[s][t] + WC6I(CddI, s, t, w, w) * delta[p][r]
                    + WC6I(CddI, w, w, p, r) * delta[s][t] + WC6I(CddI, w, w, s, t) * delta[p][r])
                    + TWO_THIRDS * g12 * Yd2 * (
                    WC6I(CddI, p, w, w, r) * delta[s][t] + WC6I(CddI, s, w, w, t) * delta[p][r]
                    + WC6I(CddI, w, t, s, w) * delta[p][r] + WC6I(CddI, w, r, p, w) * delta[s][t])
                    + ONE_THIRD * g32 * (
                    WC6I(CddI, p, w, w, t) * delta[r][s] + WC6I(CddI, s, w, w, r) * delta[p][t]
                    + WC6I(CddI, w, t, p, w) * delta[r][s] + WC6I(CddI, w, r, s, w) * delta[p][t])
                    -(ONE_THIRD / NC) * g32 * (
                    WC6I(CddI, p, w, w, r) * delta[s][t] + WC6I(CddI, s, w, w, t) * delta[p][r]
                    + WC6I(CddI, w, t, s, w) * delta[p][r] + WC6I(CddI, w, r, p, w) * delta[s][t])

                    + FOUR_THIRDS * g12 * YdYl * (
                    WC7I(CldI, w, w, p, r) * delta[s][t] + WC7I(CldI, w, w, s, t) * delta[p][r])
                    + FOUR_THIRDS * g12 * NC * YdYq * (
                    WC7I(Cqd1I, w, w, p, r) * delta[s][t] + WC7I(Cqd1I, w, w, s, t) * delta[p][r])
                    + ONE_THIRD * g32 * (
                    WC7I(Cqd8I, w, w, s, r) * delta[p][t] + WC7I(Cqd8I, w, w, p, t) * delta[r][s])
                    -(ONE_THIRD / NC) * g32 * (
                    WC7I(Cqd8I, w, w, p, r) * delta[s][t] + WC7I(Cqd8I, w, w, s, t) * delta[p][r])
                    + TWO_THIRDS * g12 * YdYe * (
                    WC7I(CedI, w, w, p, r) * delta[s][t] + WC7I(CedI, w, w, s, t) * delta[p][r])
                    + TWO_THIRDS * g12 * NC * YuYd * (
                    WC7I(Cud1I, w, w, p, r) * delta[s][t] + WC7I(Cud1I, w, w, s, t) * delta[p][r])
                    + ONE_SIXTH * g32 * (
                    WC7I(Cud8I, w, w, p, t) * delta[r][s] + WC7I(Cud8I, w, w, s, r) * delta[p][t])
                    -(ONE_SIXTH / NC) * g32 * (
                    WC7I(Cud8I, w, w, p, r) * delta[s][t] + WC7I(Cud8I, w, w, s, t) * delta[p][r])
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] += - (
                        + (yddagR[w][r] * ydR[p][v] - yddagI[w][r] * ydI[p][v]) * WC7I(Cqd1I, v, w, s, t)
                        +(yddagI[w][r] * ydR[p][v] + yddagR[w][r] * ydI[p][v]) * WC7R(Cqd1R, v, w, s, t))
                        -(
                        + (yddagR[w][t] * ydR[s][v] - yddagI[w][t] * ydI[s][v]) * WC7I(Cqd1I, v, w, p, r)
                        +(yddagI[w][t] * ydR[s][v] + yddagR[w][t] * ydI[s][v]) * WC7R(Cqd1R, v, w, p, r))

                        +(0.5 / NC)*(
                        + (ydR[p][v] * yddagR[w][r] - ydI[p][v] * yddagI[w][r]) * WC7I(Cqd8I, v, w, s, t)
                        +(ydI[p][v] * yddagR[w][r] + ydR[p][v] * yddagI[w][r]) * WC7R(Cqd8R, v, w, s, t))
                        +(0.5 / NC)*(
                        + (ydR[s][v] * yddagR[w][t] - ydI[s][v] * yddagI[w][t]) * WC7I(Cqd8I, v, w, p, r)
                        +(ydI[s][v] * yddagR[w][t] + ydR[s][v] * yddagI[w][t]) * WC7R(Cqd8R, v, w, p, r))

                        - 0.5 * (
                        + (yddagR[w][r] * ydR[s][v] - yddagI[w][r] * ydI[s][v]) * WC7I(Cqd8I, v, w, p, t)
                        +(yddagI[w][r] * ydR[s][v] + yddagR[w][r] * ydI[s][v]) * WC7R(Cqd8R, v, w, p, t))
                        - 0.5 * (
                        + (yddagR[w][t] * ydR[p][v] - yddagI[w][t] * ydI[p][v]) * WC7I(Cqd8I, v, w, s, r)
                        +(yddagI[w][t] * ydR[p][v] + yddagR[w][t] * ydI[p][v]) * WC7R(Cqd8R, v, w, s, r))

                        //RGE 2
                        ;
            }
        }
        f[c] *= loop_factor;
        c ++;
    }
    //CeuR
    for (d = 0; d < DWC7R; d ++) {
        p = WC7R_indices[d][0];
        r = WC7R_indices[d][1];
        s = WC7R_indices[d][2];
        t = WC7R_indices[d][3];
        f[c] = 2. * (geR[p][r] * WC2R(CHuR, s, t) - geI[p][r] * WC2I(CHuI, s, t))
                - 2. * (guR[s][t] * WC2R(CHeR, p, r) - guI[s][t] * WC2I(CHeI, p, r))
                //RGE 2
                + FOUR_THIRDS * g12 * YhYe * WC2R(CHuR, s, t) * delta[p][r]
                + FOUR_THIRDS * g12 * YhYu * WC2R(CHeR, p, r) * delta[s][t]
                + 12. * g12 * YuYe * WC7R(CeuR, p, r, s, t); //RGE 3

        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] += (geR[p][v] * WC7R(CeuR, v, r, s, t) - geI[p][v] * WC7I(CeuI, v, r, s, t))
                    +(guR[s][v] * WC7R(CeuR, p, r, v, t) - guI[s][v] * WC7I(CeuI, p, r, v, t))
                    +(WC7R(CeuR, p, v, s, t) * geR[v][r] - WC7I(CeuI, p, v, s, t) * geI[v][r])
                    +(WC7R(CeuR, p, r, s, v) * guR[v][t] - WC7I(CeuI, p, r, s, v) * guI[v][t])
                    //RGE 2
                    + EIGHT_THIRDS * g12 * YuYl * WC7R(CleR, w, w, p, r) * delta[s][t]
                    + EIGHT_THIRDS * g12 * YeYl * WC7R(CluR, w, w, s, t) * delta[p][r]
                    + EIGHT_THIRDS * g12 * NC * YuYq * WC7R(CqeR, w, w, p, r) * delta[s][t]
                    + EIGHT_THIRDS * g12 * NC * YeYq * WC7R(Cqu1R, w, w, s, t) * delta[p][r]

                    + FOUR_THIRDS * g12 * NC * YuYe * (
                    + WC6R(CuuR, s, t, w, w) * delta[p][r] + WC6R(CuuR, w, w, s, t) * delta[p][r]
                    )
                    + FOUR_THIRDS * g12 * YuYe * (
                    + WC6R(CuuR, s, w, w, t) * delta[p][r] + WC6R(CuuR, w, t, s, w) * delta[p][r]
                    )
                    + FOUR_THIRDS * g12 * YuYe * (
                    + WC8R(CeeR, p, r, w, w) * delta[s][t] + WC8R(CeeR, w, w, p, r) * delta[s][t]
                    + WC8R(CeeR, p, w, w, r) * delta[s][t] + WC8R(CeeR, w, r, p, w) * delta[s][t]
                    )
                    + FOUR_THIRDS * g12 * NC * Yu2 * WC7R(CeuR, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * Ye2 * WC7R(CeuR, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YdYe * WC7R(Cud1R, s, t, w, w) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YuYd * WC7R(CedR, p, r, w, w) * delta[s][t]
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] += (
                        + (yeR[p][v] * yuR[s][w] - yeI[p][v] * yuI[s][w]) * WC5(Clequ1R, v, r, w, t)
                        - (yeI[p][v] * yuR[s][w] + yeR[p][v] * yuI[s][w]) * WC5(Clequ1I, v, r, w, t)
                        +(yedagR[v][r] * yudagR[w][t] - yedagI[v][r] * yudagI[w][t]) * WC5(Clequ1R, v, p, w, s)
                        +(yedagI[v][r] * yudagR[w][t] + yedagR[v][r] * yudagI[w][t]) * WC5(Clequ1I, v, p, w, s)
                        )
                        - 2. * (
                        + (yeR[p][v] * yedagR[w][r] - yeI[p][v] * yedagI[w][r]) * WC7R(CluR, v, w, s, t)
                        -(yeI[p][v] * yedagR[w][r] + yeR[p][v] * yedagI[w][r]) * WC7I(CluI, v, w, s, t)
                        )
                        - 12. * (
                        + (yeR[p][v] * yuR[s][w] - yeI[p][v] * yuI[s][w]) * WC5(Clequ3R, v, r, w, t)
                        - (yeI[p][v] * yuR[s][w] + yeR[p][v] * yuI[s][w]) * WC5(Clequ3I, v, r, w, t)
                        +(yedagR[v][r] * yudagR[w][t] - yedagI[v][r] * yudagI[w][t]) * WC5(Clequ3R, v, p, w, s)
                        +(yedagI[v][r] * yudagR[w][t] + yedagR[v][r] * yudagI[w][t]) * WC5(Clequ3I, v, p, w, s)
                        )
                        - 2. * (
                        + (yuR[s][v] * yudagR[w][t] - yuI[s][v] * yudagI[w][t]) * WC7R(CqeR, v, w, p, r)
                        - -(yuI[s][v] * yudagR[w][t] + yuR[s][v] * yudagI[w][t]) * WC7I(CqeI, v, w, p, r)
                        )
                        ; //RGE 2
            }
        }
        f[c] *= loop_factor;
        c ++;
    }
    //CeuI
    for (d = 0; d < DWC7I; d ++) {
        p = WC7I_indices[d][0];
        r = WC7I_indices[d][1];
        s = WC7I_indices[d][2];
        t = WC7I_indices[d][3];
        f[c] = 2. * (geI[p][r] * WC2R(CHuR, s, t) + geR[p][r] * WC2I(CHuI, s, t))
                - 2. * (guI[s][t] * WC2R(CHeR, p, r) + guR[s][t] * WC2I(CHeI, p, r))
                //RGE 2
                + FOUR_THIRDS * g12 * YhYe * WC2I(CHuI, s, t) * delta[p][r]
                + FOUR_THIRDS * g12 * YhYu * WC2I(CHeI, p, r) * delta[s][t]
                + 12. * g12 * YuYe * WC7I(CeuI, p, r, s, t); //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] += (geI[p][v] * WC7R(CeuR, v, r, s, t) + geR[p][v] * WC7I(CeuI, v, r, s, t))
                    +(guI[s][v] * WC7R(CeuR, p, r, v, t) + guR[s][v] * WC7I(CeuI, p, r, v, t))
                    +(WC7R(CeuR, p, v, s, t) * geI[v][r] + WC7I(CeuI, p, v, s, t) * geR[v][r])
                    +(WC7R(CeuR, p, r, s, v) * guI[v][t] + WC7I(CeuI, p, r, s, v) * guR[v][t])
                    //RGE 2
                    + EIGHT_THIRDS * g12 * YuYl * WC7I(CleI, w, w, p, r) * delta[s][t]
                    + EIGHT_THIRDS * g12 * YeYl * WC7I(CluI, w, w, s, t) * delta[p][r]
                    + EIGHT_THIRDS * g12 * NC * YuYq * WC7I(CqeI, w, w, p, r) * delta[s][t]
                    + EIGHT_THIRDS * g12 * NC * YeYq * WC7I(Cqu1I, w, w, s, t) * delta[p][r]

                    + FOUR_THIRDS * g12 * NC * YuYe * (
                    + WC6I(CuuI, s, t, w, w) * delta[p][r] + WC6I(CuuI, w, w, s, t) * delta[p][r]
                    )
                    + FOUR_THIRDS * g12 * YuYe * (
                    + WC6I(CuuI, s, w, w, t) * delta[p][r] + WC6I(CuuI, w, t, s, w) * delta[p][r]
                    )
                    + FOUR_THIRDS * g12 * YuYe * (
                    + WC8I(CeeI, p, r, w, w) * delta[s][t] + WC8I(CeeI, w, w, p, r) * delta[s][t]
                    + WC8I(CeeI, p, w, w, r) * delta[s][t] + WC8I(CeeI, w, r, p, w) * delta[s][t]
                    )
                    + FOUR_THIRDS * g12 * NC * Yu2 * WC7I(CeuI, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * Ye2 * WC7I(CeuI, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YdYe * WC7I(Cud1I, s, t, w, w) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YuYd * WC7I(CedI, p, r, w, w) * delta[s][t]
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] += (
                        + (yeR[p][v] * yuR[s][w] - yeI[p][v] * yuI[s][w]) * WC5(Clequ1I, v, r, w, t)
                        +(yeI[p][v] * yuR[s][w] + yeR[p][v] * yuI[s][w]) * WC5(Clequ1R, v, r, w, t)
                        -(yedagR[v][r] * yudagR[w][t] - yedagI[v][r] * yudagI[w][t]) * WC5(Clequ1I, v, p, w, s)
                        +(yedagI[v][r] * yudagR[w][t] + yedagR[v][r] * yudagI[w][t]) * WC5(Clequ1R, v, p, w, s)
                        )
                        - 2. * (
                        + (yeR[p][v] * yedagR[w][r] - yeI[p][v] * yedagI[w][r]) * WC7I(CluI, v, w, s, t)
                        +(yeI[p][v] * yedagR[w][r] + yeR[p][v] * yedagI[w][r]) * WC7R(CluR, v, w, s, t)
                        )
                        - 12. * (
                        + (yeR[p][v] * yuR[s][w] - yeI[p][v] * yuI[s][w]) * WC5(Clequ3I, v, r, w, t)
                        + (yeI[p][v] * yuR[s][w] + yeR[p][v] * yuI[s][w]) * WC5(Clequ3R, v, r, w, t)
                        - (yedagR[v][r] * yudagR[w][t] - yedagI[v][r] * yudagI[w][t]) * WC5(Clequ3I, v, p, w, s)
                        + (yedagI[v][r] * yudagR[w][t] + yedagR[v][r] * yudagI[w][t]) * WC5(Clequ3R, v, p, w, s)
                        )
                        - 2. * (
                        + (yuR[s][v] * yudagR[w][t] - yuI[s][v] * yudagI[w][t]) * WC7I(CqeI, v, w, p, r)
                        +(yuI[s][v] * yudagR[w][t] + yuR[s][v] * yudagI[w][t]) * WC7R(CqeR, v, w, p, r)
                        )
                        ; //RGE 2
            }
        }
        f[c] *= loop_factor;
        c ++;
    }
    //CedR
    for (d = 0; d < DWC7R; d ++) {
        p = WC7R_indices[d][0];
        r = WC7R_indices[d][1];
        s = WC7R_indices[d][2];
        t = WC7R_indices[d][3];

        f[c] = 2. * (geR[p][r] * WC2R(CHdR, s, t) - geI[p][r] * WC2I(CHdI, s, t))
                + 2. * (gdR[s][t] * WC2R(CHeR, p, r) - gdI[s][t] * WC2I(CHeI, p, r))
                //RGE 2
                + FOUR_THIRDS * g12 * YhYe * WC2R(CHdR, s, t) * delta[p][r]
                + FOUR_THIRDS * g12 * YhYd * WC2R(CHeR, p, r) * delta[s][t]
                + 12. * YdYe * g12 * WC7R(CedR, p, r, s, t)
                ; //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] += (geR[p][v] * WC7R(CedR, v, r, s, t) - geI[p][v] * WC7I(CedI, v, r, s, t))
                    +(gdR[s][v] * WC7R(CedR, p, r, v, t) - gdI[s][v] * WC7I(CedI, p, r, v, t))
                    +(WC7R(CedR, p, v, s, t) * geR[v][r] - WC7I(CedI, p, v, s, t) * geI[v][r])
                    +(WC7R(CedR, p, r, s, v) * gdR[v][t] - WC7I(CedI, p, r, s, v) * gdI[v][t])
                    //RGE 2
                    + FOUR_THIRDS * g12 * NC * YdYe * (
                    + WC6R(CddR, s, t, w, w) * delta[p][r] + WC6R(CddR, w, w, s, t) * delta[p][r]
                    )
                    + FOUR_THIRDS * g12 * YdYe * (
                    + WC6R(CddR, s, w, w, t) * delta[p][r] + WC6R(CddR, w, t, s, w) * delta[p][r]
                    )
                    + FOUR_THIRDS * g12 * YdYe * (
                    + WC8R(CeeR, p, r, w, w) * delta[s][t] + WC8R(CeeR, w, w, p, r) * delta[s][t]
                    + WC8R(CeeR, p, w, w, r) * delta[s][t] + WC8R(CeeR, w, r, p, w) * delta[s][t]
                    )
                    + EIGHT_THIRDS * g12 * YeYl * WC7R(CldR, w, w, s, t) * delta[p][r]
                    + EIGHT_THIRDS * g12 * YdYl * WC7R(CleR, w, w, p, r) * delta[s][t]
                    + EIGHT_THIRDS * g12 * NC * YeYq * WC7R(Cqd1R, w, w, s, t) * delta[p][r]
                    + EIGHT_THIRDS * g12 * NC * YdYq * WC7R(CqeR, w, w, p, r) * delta[s][t]

                    + FOUR_THIRDS * g12 * NC * YuYd * WC7R(CeuR, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * NC * YuYe * WC7R(Cud1R, w, w, s, t) * delta[p][r]

                    + FOUR_THIRDS * g12 * NC * Yd2 * WC7R(CedR, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * Ye2 * WC7R(CedR, w, w, s, t) * delta[p][r]
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] += - 2. * (
                        + (yeR[p][v] * yedagR[w][r] - yeI[p][v] * yedagI[w][r]) * WC7R(CldR, v, w, s, t)
                        -(yeI[p][v] * yedagR[w][r] + yeR[p][v] * yedagI[w][r]) * WC7I(CldI, v, w, s, t)
                        )
                        - 2. * (
                        + (ydR[s][v] * yddagR[w][t] - ydI[s][v] * yddagI[w][t]) * WC7R(CqeR, v, w, p, r)
                        -(ydI[s][v] * yddagR[w][t] + ydR[s][v] * yddagI[w][t]) * WC7I(CqeI, v, w, p, r)
                        )
                        +(
                        + (yeR[p][v] * yddagR[w][t] - yeI[p][v] * yddagI[w][t]) * WC5(CledqR, v, r, s, w)
                        -(yeI[p][v] * yddagR[w][t] + yeR[p][v] * yddagI[w][t]) * WC5(CledqI, v, r, s, w)
                        )
                        +(
                        + (yedagR[v][r] * ydR[s][w] - yedagI[v][r] * ydI[s][w]) * WC5(CledqR, v, p, t, w)
                        +(yedagI[v][r] * ydR[s][w] + yedagR[v][r] * ydI[s][w]) * WC5(CledqI, v, p, t, w)
                        )
                        //RGE 2
                        ;
            }
        }
        f[c] *= loop_factor;
        c ++;
    }
    //CedI
    for (d = 0; d < DWC7I; d ++) {
        p = WC7I_indices[d][0];
        r = WC7I_indices[d][1];
        s = WC7I_indices[d][2];
        t = WC7I_indices[d][3];

        f[c] = 2. * (geI[p][r] * WC2R(CHdR, s, t) + geR[p][r] * WC2I(CHdI, s, t))
                + 2. * (gdI[s][t] * WC2R(CHeR, p, r) + gdR[s][t] * WC2I(CHeI, p, r))
                //RGE 2
                + FOUR_THIRDS * g12 * YhYe * WC2I(CHdI, s, t) * delta[p][r]
                + FOUR_THIRDS * g12 * YhYd * WC2I(CHeI, p, r) * delta[s][t]
                + 12. * YdYe * g12 * WC7I(CedI, p, r, s, t)
                ; //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index

            f[c] += (geI[p][v] * WC7R(CedR, v, r, s, t) + geR[p][v] * WC7I(CedI, v, r, s, t))
                    +(gdI[s][v] * WC7R(CedR, p, r, v, t) + gdR[s][v] * WC7I(CedI, p, r, v, t))
                    +(WC7R(CedR, p, v, s, t) * geI[v][r] + WC7I(CedI, p, v, s, t) * geR[v][r])
                    +(WC7R(CedR, p, r, s, v) * gdI[v][t] + WC7I(CedI, p, r, s, v) * gdR[v][t])
                    //RGE 2
                    + FOUR_THIRDS * g12 * NC * YdYe * (
                    + WC6I(CddI, s, t, w, w) * delta[p][r] + WC6I(CddI, w, w, s, t) * delta[p][r]
                    )
                    + FOUR_THIRDS * g12 * YdYe * (
                    + WC6I(CddI, s, w, w, t) * delta[p][r] + WC6I(CddI, w, t, s, w) * delta[p][r]
                    )
                    + FOUR_THIRDS * g12 * YdYe * (
                    + WC8I(CeeI, p, r, w, w) * delta[s][t] + WC8I(CeeI, w, w, p, r) * delta[s][t]
                    + WC8I(CeeI, p, w, w, r) * delta[s][t] + WC8I(CeeI, w, r, p, w) * delta[s][t]
                    )
                    + EIGHT_THIRDS * g12 * YeYl * WC7I(CldI, w, w, s, t) * delta[p][r]
                    + EIGHT_THIRDS * g12 * YdYl * WC7I(CleI, w, w, p, r) * delta[s][t]
                    + EIGHT_THIRDS * g12 * NC * YeYq * WC7I(Cqd1I, w, w, s, t) * delta[p][r]
                    + EIGHT_THIRDS * g12 * NC * YdYq * WC7I(CqeI, w, w, p, r) * delta[s][t]

                    + FOUR_THIRDS * g12 * NC * YuYd * WC7I(CeuI, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * NC * YuYe * WC7I(Cud1I, w, w, s, t) * delta[p][r]

                    + FOUR_THIRDS * g12 * NC * Yd2 * WC7I(CedI, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * Ye2 * WC7I(CedI, w, w, s, t) * delta[p][r]
                    ; //RGE 3
            ;
            for (w = 0; w < NG; w ++) {
                f[c] += - 2. * (
                        + (yeR[p][v] * yedagR[w][r] - yeI[p][v] * yedagI[w][r]) * WC7I(CldI, v, w, s, t)
                        +(yeI[p][v] * yedagR[w][r] + yeR[p][v] * yedagI[w][r]) * WC7R(CldR, v, w, s, t)
                        )
                        - 2. * (
                        + (ydR[s][v] * yddagR[w][t] - ydI[s][v] * yddagI[w][t]) * WC7I(CqeI, v, w, p, r)
                        +(ydI[s][v] * yddagR[w][t] + ydR[s][v] * yddagI[w][t]) * WC7R(CqeR, v, w, p, r)
                        )
                        +(
                        + (yeR[p][v] * yddagR[w][t] - yeI[p][v] * yddagI[w][t]) * WC5(CledqI, v, r, s, w)
                        +(yeI[p][v] * yddagR[w][t] + yeR[p][v] * yddagI[w][t]) * WC5(CledqR, v, r, s, w)
                        )
                        +(
                        - (yedagR[v][r] * ydR[s][w] - yedagI[v][r] * ydI[s][w]) * WC5(CledqI, v, p, t, w)
                        +(yedagI[v][r] * ydR[s][w] + yedagR[v][r] * ydI[s][w]) * WC5(CledqR, v, p, t, w)
                        )
                        //RGE 2
                        ;
            }
        }
        f[c] *= loop_factor;
        c ++;
    }
    //Cud1R
    for (d = 0; d < DWC7R; d ++) {
        p = WC7R_indices[d][0];
        r = WC7R_indices[d][1];
        s = WC7R_indices[d][2];
        t = WC7R_indices[d][3];
        f[c] =
                - 2. * (guR[p][r] * WC2R(CHdR, s, t) - guI[p][r] * WC2I(CHdI, s, t))
                + 2. * (gdR[s][t] * WC2R(CHuR, p, r) - gdI[s][t] * WC2I(CHuI, p, r))
                +(2. / NC)*
                (ydyudR[s][r] * WC1(CHudR, p, t) - ydyudI[s][r] * WC1(CHudI, p, t))
                +(2. / NC)*
                (yuyddR[p][t] * WC1(CHudR, r, s) + yuyddI[p][t] * WC1(CHudI, r, s))
                //RGE 2
                + FOUR_THIRDS * g12 * YhYu * WC2R(CHdR, s, t) * delta[p][r]
                + FOUR_THIRDS * g12 * YhYd * WC2R(CHuR, p, r) * delta[s][t]
                + 12. * YuYd * g12 * WC7R(Cud1R, p, r, s, t)
                + 3. * ((NC2 - 1) / NC2) * g32 * WC7R(Cud8R, p, r, s, t)
                //RGE 3
                ;
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] +=
                    + (guR[p][v] * WC7R(Cud1R, v, r, s, t) - guI[p][v] * WC7I(Cud1I, v, r, s, t))
                    +(gdR[s][v] * WC7R(Cud1R, p, r, v, t) - gdI[s][v] * WC7I(Cud1I, p, r, v, t))
                    +(WC7R(Cud1R, p, v, s, t) * guR[v][r] - WC7I(Cud1I, p, v, s, t) * guI[v][r])
                    +(WC7R(Cud1R, p, r, s, v) * gdR[v][t] - WC7I(Cud1I, p, r, s, v) * gdI[v][t])
                    //RGE 2
                    + FOUR_THIRDS * g12 * NC * YuYd * delta[s][t]*(
                    WC6R(CuuR, p, r, w, w) + WC6R(CuuR, w, w, p, r))
                    + FOUR_THIRDS * g12 * YuYd * delta[s][t]*(
                    WC6R(CuuR, p, w, w, r) + WC6R(CuuR, w, r, p, w))
                    + FOUR_THIRDS * g12 * NC * YuYd * delta[p][r]*(
                    WC6R(CddR, s, t, w, w) + WC6R(CddR, w, w, s, t))
                    + FOUR_THIRDS * g12 * YuYd * delta[p][r]*(
                    WC6R(CddR, s, w, w, t) + WC6R(CddR, w, t, s, w))
                    + EIGHT_THIRDS * g12 * NC * YdYq * WC7R(Cqu1R, w, w, p, r) * delta[s][t]
                    + EIGHT_THIRDS * g12 * NC * YuYq * WC7R(Cqd1R, w, w, s, t) * delta[p][r]
                    + EIGHT_THIRDS * g12 * YdYl * WC7R(CluR, w, w, p, r) * delta[s][t]
                    + EIGHT_THIRDS * g12 * YuYl * WC7R(CldR, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * Yd2 * WC7R(Cud1R, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * NC * Yu2 * WC7R(Cud1R, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * YdYe * WC7R(CeuR, w, w, p, r) * delta[s][t]
                    + FOUR_THIRDS * g12 * YuYe * WC7R(CedR, w, w, s, t) * delta[p][r]
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] += (1. / NC)*(
                        + (ydR[s][v] * yuR[p][w] - ydI[s][v] * yuI[p][w]) * WC5(Cquqd1R, v, r, w, t)
                        -(ydI[s][v] * yuR[p][w] + ydR[s][v] * yuI[p][w]) * WC5(Cquqd1I, v, r, w, t)
                        +(yddagR[v][t] * yudagR[w][r] - yddagI[v][t] * yudagI[w][r]) * WC5(Cquqd1R, v, p, w, s)
                        +(yddagI[v][t] * yudagR[w][r] + yddagR[v][t] * yudagI[w][r]) * WC5(Cquqd1I, v, p, w, s)
                        )
                        -(
                        + (ydR[s][w] * yuR[p][v] - ydI[s][w] * yuI[p][v]) * WC5(Cquqd1R, v, r, w, t)
                        -(ydI[s][w] * yuR[p][v] + ydR[s][w] * yuI[p][v]) * WC5(Cquqd1I, v, r, w, t)
                        )
                        -(
                        + (yddagR[w][t] * yudagR[v][r] - yddagI[w][t] * yudagI[v][r]) * WC5(Cquqd1R, v, p, w, s)
                        +(yddagI[w][t] * yudagR[v][r] + yddagR[w][t] * yudagI[v][r]) * WC5(Cquqd1I, v, p, w, s)
                        )
                        + 0.5 * ((NC2 - 1.) / NC2)*(
                        + (ydR[s][v] * yuR[p][w] - ydI[s][v] * yuI[p][w]) * WC5(Cquqd8R, v, r, w, t)
                        -(ydI[s][v] * yuR[p][w] + ydR[s][v] * yuI[p][w]) * WC5(Cquqd8I, v, r, w, t)
                        +(yddagR[v][t] * yudagR[w][r] - yddagI[v][t] * yudagI[w][r]) * WC5(Cquqd8R, v, p, w, s)
                        +(yddagI[v][t] * yudagR[w][r] + yddagR[v][t] * yudagI[w][r]) * WC5(Cquqd8I, v, p, w, s)
                        )
                        - 2. * (
                        + (yuR[p][v] * yudagR[w][r] - yuI[p][v] * yudagI[w][r]) * WC7R(Cqd1R, v, w, s, t)
                        -(yuI[p][v] * yudagR[w][r] + yuR[p][v] * yudagI[w][r]) * WC7I(Cqd1I, v, w, s, t)
                        )
                        - 2. * (
                        + (ydR[s][v] * yddagR[w][t] - ydI[s][v] * yddagI[w][t]) * WC7R(Cqu1R, v, w, p, r)
                        -(ydI[s][v] * yddagR[w][t] + ydR[s][v] * yddagI[w][t]) * WC7I(Cqu1I, v, w, p, r)
                        )
                        //RGE 2 
                        ;
            }
        }
        f[c] *= loop_factor;
        c ++;
    }
    //Cud1I
    for (d = 0; d < DWC7I; d ++) {
        p = WC7I_indices[d][0];
        r = WC7I_indices[d][1];
        s = WC7I_indices[d][2];
        t = WC7I_indices[d][3];
        f[c] =
                - 2. * (guI[p][r] * WC2R(CHdR, s, t) + guR[p][r] * WC2I(CHdI, s, t))
                + 2. * (gdI[s][t] * WC2R(CHuR, p, r) + gdR[s][t] * WC2I(CHuI, p, r))
                +(2. / NC)*
                (ydyudI[s][r] * WC1(CHudR, p, t) + ydyudR[s][r] * WC1(CHudI, p, t))
                +(2. / NC)*
                (yuyddI[p][t] * WC1(CHudR, r, s) - yuyddR[p][t] * WC1(CHudI, r, s))
                //RGE 2
                + FOUR_THIRDS * g12 * YhYu * WC2I(CHdI, s, t) * delta[p][r]
                + FOUR_THIRDS * g12 * YhYd * WC2I(CHuI, p, r) * delta[s][t]
                + 12. * YuYd * g12 * WC7I(Cud1I, p, r, s, t)
                + 3. * ((NC2 - 1) / NC2) * g32 * WC7I(Cud8I, p, r, s, t)
                //RGE 3
                ;
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] +=
                    + (guI[p][v] * WC7R(Cud1R, v, r, s, t) + guR[p][v] * WC7I(Cud1I, v, r, s, t))
                    +(gdI[s][v] * WC7R(Cud1R, p, r, v, t) + gdR[s][v] * WC7I(Cud1I, p, r, v, t))
                    +(WC7R(Cud1R, p, v, s, t) * guI[v][r] + WC7I(Cud1I, p, v, s, t) * guR[v][r])
                    +(WC7R(Cud1R, p, r, s, v) * gdI[v][t] + WC7I(Cud1I, p, r, s, v) * gdR[v][t])
                    //RGE 2
                    + FOUR_THIRDS * g12 * NC * YuYd * delta[s][t]*(
                    WC6I(CuuI, p, r, w, w) + WC6I(CuuI, w, w, p, r))
                    + FOUR_THIRDS * g12 * YuYd * delta[s][t]*(
                    WC6I(CuuI, p, w, w, r) + WC6I(CuuI, w, r, p, w))
                    + FOUR_THIRDS * g12 * NC * YuYd * delta[p][r]*(
                    WC6I(CddI, s, t, w, w) + WC6I(CddI, w, w, s, t))
                    + FOUR_THIRDS * g12 * YuYd * delta[p][r]*(
                    WC6I(CddI, s, w, w, t) + WC6I(CddI, w, t, s, w))
                    + EIGHT_THIRDS * g12 * NC * YdYq * WC7I(Cqu1I, w, w, p, r) * delta[s][t]
                    + EIGHT_THIRDS * g12 * NC * YuYq * WC7I(Cqd1I, w, w, s, t) * delta[p][r]
                    + EIGHT_THIRDS * g12 * YdYl * WC7I(CluI, w, w, p, r) * delta[s][t]
                    + EIGHT_THIRDS * g12 * YuYl * WC7I(CldI, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * Yd2 * WC7I(Cud1I, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * NC * Yu2 * WC7I(Cud1I, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * YdYe * WC7I(CeuI, w, w, p, r) * delta[s][t]
                    + FOUR_THIRDS * g12 * YuYe * WC7I(CedI, w, w, s, t) * delta[p][r]
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] +=
                        + (1. / NC)*(
                        + (ydR[s][v] * yuR[p][w] - ydI[s][v] * yuI[p][w]) * WC5(Cquqd1I, v, r, w, t)
                        +(ydI[s][v] * yuR[p][w] + ydR[s][v] * yuI[p][w]) * WC5(Cquqd1R, v, r, w, t)
                        -(yddagR[v][t] * yudagR[w][r] - yddagI[v][t] * yudagI[w][r]) * WC5(Cquqd1I, v, p, w, s)
                        +(yddagI[v][t] * yudagR[w][r] + yddagR[v][t] * yudagI[w][r]) * WC5(Cquqd1R, v, p, w, s)
                        )
                        -(
                        + (ydR[s][w] * yuR[p][v] - ydI[s][w] * yuI[p][v]) * WC5(Cquqd1I, v, r, w, t)
                        +(ydI[s][w] * yuR[p][v] + ydR[s][w] * yuI[p][v]) * WC5(Cquqd1R, v, r, w, t)
                        )
                        -(
                        - (yddagR[w][t] * yudagR[v][r] - yddagI[w][t] * yudagI[v][r]) * WC5(Cquqd1I, v, p, w, s)
                        +(yddagI[w][t] * yudagR[v][r] + yddagR[w][t] * yudagI[v][r]) * WC5(Cquqd1R, v, p, w, s)
                        )
                        + 0.5 * ((NC2 - 1.) / NC2)*(
                        + (ydR[s][v] * yuR[p][w] - ydI[s][v] * yuI[p][w]) * WC5(Cquqd8I, v, r, w, t)
                        +(ydI[s][v] * yuR[p][w] + ydR[s][v] * yuI[p][w]) * WC5(Cquqd8R, v, r, w, t)
                        -(yddagR[v][t] * yudagR[w][r] - yddagI[v][t] * yudagI[w][r]) * WC5(Cquqd8I, v, p, w, s)
                        +(yddagI[v][t] * yudagR[w][r] + yddagR[v][t] * yudagI[w][r]) * WC5(Cquqd8R, v, p, w, s)
                        )
                        - 2. * (
                        + (yuR[p][v] * yudagR[w][r] - yuI[p][v] * yudagI[w][r]) * WC7I(Cqd1I, v, w, s, t)
                        +(yuI[p][v] * yudagR[w][r] + yuR[p][v] * yudagI[w][r]) * WC7R(Cqd1R, v, w, s, t)
                        )
                        - 2. * (
                        + (ydR[s][v] * yddagR[w][t] - ydI[s][v] * yddagI[w][t]) * WC7I(Cqu1I, v, w, p, r)
                        +(ydI[s][v] * yddagR[w][t] + ydR[s][v] * yddagI[w][t]) * WC7R(Cqu1R, v, w, p, r)
                        )
                        //RGE 2 
                        ;
            }
        }
        f[c] *= loop_factor;
        c ++;
    }
    //Cud8R
    for (d = 0; d < DWC7R; d ++) {
        p = WC7R_indices[d][0];
        r = WC7R_indices[d][1];
        s = WC7R_indices[d][2];
        t = WC7R_indices[d][3];
        f[c] =
                + 4. * (ydyudR[s][r] * WC1(CHudR, p, t) - ydyudI[s][r] * WC1(CHudI, p, t))
                + 4. * (yuyddR[p][t] * WC1(CHudR, r, s) + yuyddI[p][t] * WC1(CHudI, r, s))
                //RGE 2
                + 12. * (YuYd * g12 - g32 / NC) * WC7R(Cud8R, p, r, s, t)
                + 12. * g32 * WC7R(Cud1R, p, r, s, t)
                ; //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] +=
                    + (guR[p][v] * WC7R(Cud8R, v, r, s, t) - guI[p][v] * WC7I(Cud8I, v, r, s, t))
                    +(gdR[s][v] * WC7R(Cud8R, p, r, v, t) - gdI[s][v] * WC7I(Cud8I, p, r, v, t))
                    +(WC7R(Cud8R, p, v, s, t) * guR[v][r] - WC7I(Cud8I, p, v, s, t) * guI[v][r])
                    +(WC7R(Cud8R, p, r, s, v) * gdR[v][t] - WC7I(Cud8I, p, r, s, v) * gdI[v][t])
                    //RGE 2
                    + FOUR_THIRDS * g32 * delta[s][t]*(
                    + WC6R(CuuR, p, w, w, r) + WC6R(CuuR, w, r, p, w)
                    )
                    + FOUR_THIRDS * g32 * delta[p][r]*(
                    + WC6R(CddR, s, w, w, t) + WC6R(CddR, w, t, s, w)
                    )
                    + FOUR_THIRDS * g32 * WC7R(Cqu8R, w, w, p, r) * delta[s][t]
                    + FOUR_THIRDS * g32 * WC7R(Cqd8R, w, w, s, t) * delta[p][r]
                    + TWO_THIRDS * g32 * (
                    + WC7R(Cud8R, p, r, w, w) * delta[s][t] + WC7R(Cud8R, w, w, s, t) * delta[p][r]
                    )
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] += + 2. * (
                        + (ydR[s][v] * yuR[p][w] - ydI[s][v] * yuI[p][w]) * WC5(Cquqd1R, v, r, w, t)
                        -(ydI[s][v] * yuR[p][w] + ydR[s][v] * yuI[p][w]) * WC5(Cquqd1I, v, r, w, t)
                        +(yddagR[v][t] * yudagR[w][r] - yddagI[v][t] * yudagI[w][r]) * WC5(Cquqd1R, v, p, w, s)
                        +(yddagI[v][t] * yudagR[w][r] + yddagR[v][t] * yudagI[w][r]) * WC5(Cquqd1I, v, p, w, s)
                        )
                        - 2. * (
                        + (yuR[p][v] * yudagR[w][r] - yuI[p][v] * yudagI[w][r]) * WC7R(Cqd8R, v, w, s, t)
                        -(yuI[p][v] * yudagR[w][r] + yuR[p][v] * yudagI[w][r]) * WC7I(Cqd8I, v, w, s, t)
                        )
                        - 2. * (
                        + (ydR[s][v] * yddagR[w][t] - ydI[s][v] * yddagI[w][t]) * WC7R(Cqu8R, v, w, p, r)
                        -(ydI[s][v] * yddagR[w][t] + ydR[s][v] * yddagI[w][t]) * WC7I(Cqu8I, v, w, p, r)
                        )
                        -(1. / NC)*(
                        + (ydR[s][v] * yuR[p][w] - ydI[s][v] * yuI[p][w]) * WC5(Cquqd8R, v, r, w, t)
                        -(ydI[s][v] * yuR[p][w] + ydR[s][v] * yuI[p][w]) * WC5(Cquqd8I, v, r, w, t)
                        +(yddagR[v][t] * yudagR[w][r] - yddagI[v][t] * yudagI[w][r]) * WC5(Cquqd8R, v, p, w, s)
                        +(yddagI[v][t] * yudagR[w][r] + yddagR[v][t] * yudagI[w][r]) * WC5(Cquqd8I, v, p, w, s)
                        )
                        -(
                        + (ydR[s][w] * yuR[p][v] - ydI[s][w] * yuI[p][v]) * WC5(Cquqd8R, v, r, w, t)
                        - (ydI[s][w] * yuR[p][v] + ydR[s][w] * yuI[p][v]) * WC5(Cquqd8I, v, r, w, t)
                        +(yddagR[w][t] * yudagR[v][r] - yddagI[w][t] * yudagI[v][r]) * WC5(Cquqd8R, v, p, w, s)
                        +(yddagI[w][t] * yudagR[v][r] + yddagR[w][t] * yudagI[v][r]) * WC5(Cquqd8I, v, p, w, s)
                        )
                        ; //RGE 2
            }
        }
        f[c] *= loop_factor;
        c ++;
    }
    //Cud8I
    for (d = 0; d < DWC7I; d ++) {
        p = WC7I_indices[d][0];
        r = WC7I_indices[d][1];
        s = WC7I_indices[d][2];
        t = WC7I_indices[d][3];
        f[c] =
                + 4. * (ydyudI[s][r] * WC1(CHudR, p, t) + ydyudR[s][r] * WC1(CHudI, p, t))
                + 4. * (yuyddI[p][t] * WC1(CHudR, r, s) - yuyddR[p][t] * WC1(CHudI, r, s))
                //RGE 2
                + 12. * (YuYd * g12 - g32 / NC) * WC7I(Cud8I, p, r, s, t)
                + 12. * g32 * WC7I(Cud1I, p, r, s, t)
                ; //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] +=
                    + (guI[p][v] * WC7R(Cud8R, v, r, s, t) + guR[p][v] * WC7I(Cud8I, v, r, s, t))
                    +(gdI[s][v] * WC7R(Cud8R, p, r, v, t) + gdR[s][v] * WC7I(Cud8I, p, r, v, t))
                    +(WC7R(Cud8R, p, v, s, t) * guI[v][r] + WC7I(Cud8I, p, v, s, t) * guR[v][r])
                    +(WC7R(Cud8R, p, r, s, v) * gdI[v][t] + WC7I(Cud8I, p, r, s, v) * gdR[v][t])
                    //RGE 2
                    + FOUR_THIRDS * g32 * delta[s][t]*(
                    + WC6I(CuuI, p, w, w, r) + WC6I(CuuI, w, r, p, w)
                    )
                    + FOUR_THIRDS * g32 * delta[p][r]*(
                    + WC6I(CddI, s, w, w, t) + WC6I(CddI, w, t, s, w)
                    )
                    + FOUR_THIRDS * g32 * WC7I(Cqu8I, w, w, p, r) * delta[s][t]
                    + FOUR_THIRDS * g32 * WC7I(Cqd8I, w, w, s, t) * delta[p][r]
                    + TWO_THIRDS * g32 * (
                    + WC7I(Cud8I, p, r, w, w) * delta[s][t] + WC7I(Cud8I, w, w, s, t) * delta[p][r]
                    )
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] += + 2. * (
                        + (ydR[s][v] * yuR[p][w] - ydI[s][v] * yuI[p][w]) * WC5(Cquqd1I, v, r, w, t)
                        +(ydI[s][v] * yuR[p][w] + ydR[s][v] * yuI[p][w]) * WC5(Cquqd1R, v, r, w, t)
                        -(yddagR[v][t] * yudagR[w][r] - yddagI[v][t] * yudagI[w][r]) * WC5(Cquqd1I, v, p, w, s)
                        +(yddagI[v][t] * yudagR[w][r] + yddagR[v][t] * yudagI[w][r]) * WC5(Cquqd1R, v, p, w, s)
                        )
                        - 2. * (
                        + (yuR[p][v] * yudagR[w][r] - yuI[p][v] * yudagI[w][r]) * WC7I(Cqd8I, v, w, s, t)
                        +(yuI[p][v] * yudagR[w][r] + yuR[p][v] * yudagI[w][r]) * WC7R(Cqd8R, v, w, s, t)
                        )
                        - 2. * (
                        + (ydR[s][v] * yddagR[w][t] - ydI[s][v] * yddagI[w][t]) * WC7I(Cqu8I, v, w, p, r)
                        +(ydI[s][v] * yddagR[w][t] + ydR[s][v] * yddagI[w][t]) * WC7R(Cqu8R, v, w, p, r)
                        )
                        -(1. / NC)*(
                        + (ydR[s][v] * yuR[p][w] - ydI[s][v] * yuI[p][w]) * WC5(Cquqd8I, v, r, w, t)
                        +(ydI[s][v] * yuR[p][w] + ydR[s][v] * yuI[p][w]) * WC5(Cquqd8R, v, r, w, t)
                        -(yddagR[v][t] * yudagR[w][r] - yddagI[v][t] * yudagI[w][r]) * WC5(Cquqd8I, v, p, w, s)
                        +(yddagI[v][t] * yudagR[w][r] + yddagR[v][t] * yudagI[w][r]) * WC5(Cquqd8R, v, p, w, s)
                        )
                        -(
                        + (ydR[s][w] * yuR[p][v] - ydI[s][w] * yuI[p][v]) * WC5(Cquqd8I, v, r, w, t)
                        + (ydI[s][w] * yuR[p][v] + ydR[s][w] * yuI[p][v]) * WC5(Cquqd8R, v, r, w, t)
                        -(yddagR[w][t] * yudagR[v][r] - yddagI[w][t] * yudagI[v][r]) * WC5(Cquqd8I, v, p, w, s)
                        +(yddagI[w][t] * yudagR[v][r] + yddagR[w][t] * yudagI[v][r]) * WC5(Cquqd8R, v, p, w, s)
                        )
                        ; //RGE 2

            }
        }
        f[c] *= loop_factor;
        c ++;
    }



    //----------------------RGE SMEFT class 8_LLRR----------------------

    //CleR
    for (d = 0; d < DWC7R; d ++) {
        p = WC7R_indices[d][0];
        r = WC7R_indices[d][1];
        s = WC7R_indices[d][2];
        t = WC7R_indices[d][3];
        f[c] =
                + (yeR[s][r] * xieR[p][t] - yeI[s][r] * xieI[p][t])
                +(yedagR[p][t] * xieR[r][s] + yedagI[p][t] * xieI[r][s])
                -(yedyeR[p][r] * WC2R(CHeR, s, t) - yedyeI[p][r] * WC2I(CHeI, s, t))
                + 2. * (geR[s][t] * WC2R(CHl1R, p, r) - geI[s][t] * WC2I(CHl1I, p, r))
                //RGE 2
                + FOUR_THIRDS * g12 * YhYl * WC2R(CHeR, s, t) * delta[p][r]
                + FOUR_THIRDS * g12 * YhYe * WC2R(CHl1R, p, r) * delta[s][t]
                - 12. * YeYl * g12 * WC7R(CleR, p, r, s, t)
                ; //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] +=
                    + (glR[p][v] * WC7R(CleR, v, r, s, t) - glI[p][v] * WC7I(CleI, v, r, s, t))
                    +(geR[s][v] * WC7R(CleR, p, r, v, t) - geI[s][v] * WC7I(CleI, p, r, v, t))
                    +(WC7R(CleR, p, v, s, t) * glR[v][r] - WC7I(CleI, p, v, s, t) * glI[v][r])
                    +(WC7R(CleR, p, r, s, v) * geR[v][t] - WC7I(CleI, p, r, s, v) * geI[v][t])
                    //RGE 2
                    + EIGHT_THIRDS * g12 * YeYl * delta[s][t]*(
                    + WC6R(CllR, p, r, w, w) + WC6R(CllR, w, w, p, r)
                    )
                    + FOUR_THIRDS * g12 * YeYl * delta[s][t]*(
                    + WC6R(CllR, p, w, w, r) + WC6R(CllR, w, r, p, w)
                    )
                    + EIGHT_THIRDS * g12 * NC * YeYq * WC7R(Clq1R, p, r, w, w) * delta[s][t]
                    + EIGHT_THIRDS * g12 * NC * YlYq * WC7R(CqeR, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * Ye2 * WC7R(CleR, p, r, w, w) * delta[s][t]
                    + EIGHT_THIRDS * g12 * Yl2 * WC7R(CleR, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YuYe * WC7R(CluR, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * NC * YdYe * WC7R(CldR, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * NC * YuYl * WC7R(CeuR, s, t, w, w) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YdYl * WC7R(CedR, s, t, w, w) * delta[p][r]
                    + FOUR_THIRDS * g12 * YeYl * delta[p][r]* (
                    + WC8R(CeeR, s, t, w, w) + WC8R(CeeR, s, w, w, t)
                    + WC8R(CeeR, w, t, s, w) + WC8R(CeeR, w, w, s, t)
                    )
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] +=
                        - (
                        + (yedagR[p][v] * yeR[w][r] - yedagI[p][v] * yeI[w][r]) * WC8R(CeeR, v, t, s, w)
                        -(yedagI[p][v] * yeR[w][r] + yedagR[p][v] * yeI[w][r]) * WC8I(CeeI, v, t, s, w)
                        )
                        -(
                        + (yedagR[p][w] * yeR[v][r] - yedagI[p][w] * yeI[v][r]) * WC8R(CeeR, w, t, s, v)
                        -(yedagI[p][w] * yeR[v][r] + yedagR[p][w] * yeI[v][r]) * WC8I(CeeI, w, t, s, v)
                        )
                        - 2. * (
                        + (yedagR[p][v] * yeR[w][r] - yedagI[p][v] * yeI[w][r]) * WC8R(CeeR, v, w, s, t)
                        -(yedagI[p][v] * yeR[w][r] + yedagR[p][v] * yeI[w][r]) * WC8I(CeeI, v, w, s, t)
                        )
                        +(
                        + (yedagR[p][w] * yeR[s][v] - yedagI[p][w] * yeI[s][v]) * WC7R(CleR, v, r, w, t)
                        -(yedagI[p][w] * yeR[s][v] + yedagR[p][w] * yeI[s][v]) * WC7I(CleI, v, r, w, t)
                        )
                        -(
                        + (yedagR[w][t] * yeR[s][v] - yedagI[w][t] * yeI[s][v]) * WC6R(CllR, p, w, v, r)
                        -(yedagI[w][t] * yeR[s][v] + yedagR[w][t] * yeI[s][v]) * WC6I(CllI, p, w, v, r)
                        )
                        -(
                        + (yedagR[v][t] * yeR[s][w] - yedagI[v][t] * yeI[s][w]) * WC6R(CllR, p, v, w, r)
                        -(yedagI[v][t] * yeR[s][w] + yedagR[v][t] * yeI[s][w]) * WC6I(CllI, p, v, w, r)
                        )
                        - 4. * (
                        + (yedagR[w][t] * yeR[s][v] - yedagI[w][t] * yeI[s][v]) * WC6R(CllR, p, r, v, w)
                        -(yedagI[w][t] * yeR[s][v] + yedagR[w][t] * yeI[s][v]) * WC6I(CllI, p, r, v, w)
                        )
                        +(
                        + (yedagR[v][t] * yeR[w][r] - yedagI[v][t] * yeI[w][r]) * WC7R(CleR, p, v, s, w)
                        -(yedagI[v][t] * yeR[w][r] + yedagR[v][t] * yeI[w][r]) * WC7I(CleI, p, v, s, w)
                        )
                        //RGE 2
                        ;
            }
        }
        f[c] *= loop_factor;
        c ++;
    }
    //CleI
    for (d = 0; d < DWC7I; d ++) {
        p = WC7I_indices[d][0];
        r = WC7I_indices[d][1];
        s = WC7I_indices[d][2];
        t = WC7I_indices[d][3];
        f[c] =
                + (yeI[s][r] * xieR[p][t] + yeR[s][r] * xieI[p][t])
                +(yedagI[p][t] * xieR[r][s] - yedagR[p][t] * xieI[r][s])
                -(yedyeI[p][r] * WC2R(CHeR, s, t) + yedyeR[p][r] * WC2I(CHeI, s, t))
                + 2. * (geI[s][t] * WC2R(CHl1R, p, r) + geR[s][t] * WC2I(CHl1I, p, r))
                //RGE 2
                + FOUR_THIRDS * g12 * YhYl * WC2I(CHeI, s, t) * delta[p][r]
                + FOUR_THIRDS * g12 * YhYe * WC2I(CHl1I, p, r) * delta[s][t]
                - 12. * YeYl * g12 * WC7I(CleI, p, r, s, t)
                ; //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] +=
                    + (glI[p][v] * WC7R(CleR, v, r, s, t) + glR[p][v] * WC7I(CleI, v, r, s, t))
                    +(geI[s][v] * WC7R(CleR, p, r, v, t) + geR[s][v] * WC7I(CleI, p, r, v, t))
                    +(WC7R(CleR, p, v, s, t) * glI[v][r] + WC7I(CleI, p, v, s, t) * glR[v][r])
                    +(WC7R(CleR, p, r, s, v) * geI[v][t] + WC7I(CleI, p, r, s, v) * geR[v][t])
                    //RGE 2
                    + EIGHT_THIRDS * g12 * YeYl * delta[s][t]*(
                    + WC6I(CllI, p, r, w, w) + WC6I(CllI, w, w, p, r)
                    )
                    + FOUR_THIRDS * g12 * YeYl * delta[s][t]*(
                    + WC6I(CllI, p, w, w, r) + WC6I(CllI, w, r, p, w)
                    )
                    + EIGHT_THIRDS * g12 * NC * YeYq * WC7I(Clq1I, p, r, w, w) * delta[s][t]
                    + EIGHT_THIRDS * g12 * NC * YlYq * WC7I(CqeI, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * Ye2 * WC7I(CleI, p, r, w, w) * delta[s][t]
                    + EIGHT_THIRDS * g12 * Yl2 * WC7I(CleI, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YuYe * WC7I(CluI, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * NC * YdYe * WC7I(CldI, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * NC * YuYl * WC7I(CeuI, s, t, w, w) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YdYl * WC7I(CedI, s, t, w, w) * delta[p][r]
                    + FOUR_THIRDS * g12 * YeYl * delta[p][r]* (
                    + WC8I(CeeI, s, t, w, w) + WC8I(CeeI, s, w, w, t)
                    + WC8I(CeeI, w, t, s, w) + WC8I(CeeI, w, w, s, t)
                    )
                    ; //RGE 3

            for (w = 0; w < NG; w ++) {
                f[c] += - (
                        + (yedagR[p][v] * yeR[w][r] - yedagI[p][v] * yeI[w][r]) * WC8I(CeeI, v, t, s, w)
                        +(yedagI[p][v] * yeR[w][r] + yedagR[p][v] * yeI[w][r]) * WC8R(CeeR, v, t, s, w)
                        )
                        -(
                        + (yedagR[p][w] * yeR[v][r] - yedagI[p][w] * yeI[v][r]) * WC8I(CeeI, w, t, s, v)
                        +(yedagI[p][w] * yeR[v][r] + yedagR[p][w] * yeI[v][r]) * WC8R(CeeR, w, t, s, v)
                        )
                        - 2. * (
                        + (yedagR[p][v] * yeR[w][r] - yedagI[p][v] * yeI[w][r]) * WC8I(CeeI, v, w, s, t)
                        +(yedagI[p][v] * yeR[w][r] + yedagR[p][v] * yeI[w][r]) * WC8R(CeeR, v, w, s, t)
                        )
                        +(
                        + (yedagR[p][w] * yeR[s][v] - yedagI[p][w] * yeI[s][v]) * WC7I(CleI, v, r, w, t)
                        +(yedagI[p][w] * yeR[s][v] + yedagR[p][w] * yeI[s][v]) * WC7R(CleR, v, r, w, t)
                        )
                        -(
                        + (yedagR[w][t] * yeR[s][v] - yedagI[w][t] * yeI[s][v]) * WC6I(CllI, p, w, v, r)
                        +(yedagI[w][t] * yeR[s][v] + yedagR[w][t] * yeI[s][v]) * WC6R(CllR, p, w, v, r)
                        )
                        -(
                        + (yedagR[v][t] * yeR[s][w] - yedagI[v][t] * yeI[s][w]) * WC6I(CllI, p, v, w, r)
                        +(yedagI[v][t] * yeR[s][w] + yedagR[v][t] * yeI[s][w]) * WC6R(CllR, p, v, w, r)
                        )
                        - 4. * (
                        + (yedagR[w][t] * yeR[s][v] - yedagI[w][t] * yeI[s][v]) * WC6I(CllI, p, r, v, w)
                        +(yedagI[w][t] * yeR[s][v] + yedagR[w][t] * yeI[s][v]) * WC6R(CllR, p, r, v, w)
                        )
                        +(
                        + (yedagR[v][t] * yeR[w][r] - yedagI[v][t] * yeI[w][r]) * WC7I(CleI, p, v, s, w)
                        +(yedagI[v][t] * yeR[w][r] + yedagR[v][t] * yeI[w][r]) * WC7R(CleR, p, v, s, w)
                        )
                        //RGE 2
                        ;
            }
        }
        f[c] *= loop_factor;
        c ++;
    }




    //CluR
    for (d = 0; d < DWC7R; d ++) {
        p = WC7R_indices[d][0];
        r = WC7R_indices[d][1];
        s = WC7R_indices[d][2];
        t = WC7R_indices[d][3];
        f[c] =
                - (yedyeR[p][r] * WC2R(CHuR, s, t) - yedyeI[p][r] * WC2I(CHuI, s, t))
                - 2. * (guR[s][t] * WC2R(CHl1R, p, r) - guI[s][t] * WC2I(CHl1I, p, r))
                //RGE 2
                + FOUR_THIRDS * g12 * YhYl * WC2R(CHuR, s, t) * delta[p][r]
                + FOUR_THIRDS * g12 * YhYu * WC2R(CHl1R, p, r) * delta[s][t]
                - 12. * g12 * YuYl * WC7R(CluR, p, r, s, t)
                ; //RGE 3

        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] +=
                    + (glR[p][v] * WC7R(CluR, v, r, s, t) - glI[p][v] * WC7I(CluI, v, r, s, t))
                    +(guR[s][v] * WC7R(CluR, p, r, v, t) - guI[s][v] * WC7I(CluI, p, r, v, t))
                    +(WC7R(CluR, p, v, s, t) * glR[v][r] - WC7I(CluI, p, v, s, t) * glI[v][r])
                    +(WC7R(CluR, p, r, s, v) * guR[v][t] - WC7I(CluI, p, r, s, v) * guI[v][t])
                    //RGE 2
                    + EIGHT_THIRDS * g12 * YuYl * delta[s][t]*(
                    + WC6R(CllR, p, r, w, w) + WC6R(CllR, w, w, p, r)
                    )
                    + FOUR_THIRDS * g12 * YuYl * delta[s][t]*(
                    + WC6R(CllR, p, w, w, r) + WC6R(CllR, w, r, p, w)
                    )
                    + EIGHT_THIRDS * g12 * NC * YuYq * WC7R(Clq1R, p, r, w, w) * delta[s][t]
                    + EIGHT_THIRDS * g12 * NC * YlYq * WC7R(Cqu1R, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * Yu2 * WC7R(CluR, p, r, w, w) * delta[s][t]
                    + EIGHT_THIRDS * g12 * Yl2 * WC7R(CluR, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YuYd * WC7R(CldR, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * YuYe * WC7R(CleR, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * NC * YdYl * WC7R(Cud1R, s, t, w, w) * delta[p][r]
                    + FOUR_THIRDS * g12 * YeYl * WC7R(CeuR, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YuYl * delta[p][r]* (
                    + WC6R(CuuR, s, t, w, w) + WC6R(CuuR, w, w, s, t)
                    )
                    + FOUR_THIRDS * g12 * YuYl * delta[p][r]* (
                    + WC6R(CuuR, s, w, w, t) + WC6R(CuuR, w, t, s, w)
                    )
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] +=
                        - 0.5 * (
                        + (yeR[v][r] * yuR[s][w] - yeI[v][r] * yuI[s][w]) * WC5(Clequ1R, p, v, w, t)
                        -(yeI[v][r] * yuR[s][w] + yeR[v][r] * yuI[s][w]) * WC5(Clequ1I, p, v, w, t)
                        +(yedagR[p][v] * yudagR[w][t] - yedagI[p][v] * yudagI[w][t]) * WC5(Clequ1R, r, v, w, s)
                        +(yedagI[p][v] * yudagR[w][t] + yedagR[p][v] * yudagI[w][t]) * WC5(Clequ1I, r, v, w, s)
                        )
                        - 2. * (
                        + (yuR[s][v] * yudagR[w][t] - yuI[s][v] * yudagI[w][t]) * WC7R(Clq1R, p, r, v, w)
                        -(yuI[s][v] * yudagR[w][t] + yuR[s][v] * yudagI[w][t]) * WC7I(Clq1I, p, r, v, w)
                        )
                        - 6. * (
                        + (yeR[v][r] * yuR[s][w] - yeI[v][r] * yuI[s][w]) * WC5(Clequ3R, p, v, w, t)
                        -(yeI[v][r] * yuR[s][w] + yeR[v][r] * yuI[s][w]) * WC5(Clequ3I, p, v, w, t)
                        +(yedagR[p][v] * yudagR[w][t] - yedagI[p][v] * yudagI[w][t]) * WC5(Clequ3R, r, v, w, s)
                        +(yedagI[p][v] * yudagR[w][t] + yedagR[p][v] * yudagI[w][t]) * WC5(Clequ3I, r, v, w, s)
                        )
                        -(
                        + (yeR[w][r] * yedagR[p][v] - yeI[w][r] * yedagI[p][v]) * WC7R(CeuR, v, w, s, t)
                        -(yeI[w][r] * yedagR[p][v] + yeR[w][r] * yedagI[p][v]) * WC7I(CeuI, v, w, s, t)
                        )
                        ; //RGE 2
            }
        }
        f[c] *= loop_factor;
        c ++;
    }
    //CluI
    for (d = 0; d < DWC7I; d ++) {
        p = WC7I_indices[d][0];
        r = WC7I_indices[d][1];
        s = WC7I_indices[d][2];
        t = WC7I_indices[d][3];
        f[c] =
                - (yedyeI[p][r] * WC2R(CHuR, s, t) + yedyeR[p][r] * WC2I(CHuI, s, t))
                - 2. * (guI[s][t] * WC2R(CHl1R, p, r) + guR[s][t] * WC2I(CHl1I, p, r))
                //RGE 2
                + FOUR_THIRDS * g12 * YhYl * WC2I(CHuI, s, t) * delta[p][r]
                + FOUR_THIRDS * g12 * YhYu * WC2I(CHl1I, p, r) * delta[s][t]
                - 12. * g12 * YuYl * WC7I(CluI, p, r, s, t)
                ; //RGE 3

        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] +=
                    + (glI[p][v] * WC7R(CluR, v, r, s, t) + glR[p][v] * WC7I(CluI, v, r, s, t))
                    +(guI[s][v] * WC7R(CluR, p, r, v, t) + guR[s][v] * WC7I(CluI, p, r, v, t))
                    +(WC7R(CluR, p, v, s, t) * glI[v][r] + WC7I(CluI, p, v, s, t) * glR[v][r])
                    +(WC7R(CluR, p, r, s, v) * guI[v][t] + WC7I(CluI, p, r, s, v) * guR[v][t])
                    //RGE 2
                    + EIGHT_THIRDS * g12 * YuYl * delta[s][t]*(
                    + WC6I(CllI, p, r, w, w) + WC6I(CllI, w, w, p, r)
                    )
                    + FOUR_THIRDS * g12 * YuYl * delta[s][t]*(
                    + WC6I(CllI, p, w, w, r) + WC6I(CllI, w, r, p, w)
                    )
                    + EIGHT_THIRDS * g12 * NC * YuYq * WC7I(Clq1I, p, r, w, w) * delta[s][t]
                    + EIGHT_THIRDS * g12 * NC * YlYq * WC7I(Cqu1I, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * Yu2 * WC7I(CluI, p, r, w, w) * delta[s][t]
                    + EIGHT_THIRDS * g12 * Yl2 * WC7I(CluI, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YuYd * WC7I(CldI, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * YuYe * WC7I(CleI, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * NC * YdYl * WC7I(Cud1I, s, t, w, w) * delta[p][r]
                    + FOUR_THIRDS * g12 * YeYl * WC7I(CeuI, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YuYl * delta[p][r]* (
                    + WC6I(CuuI, s, t, w, w) + WC6I(CuuI, w, w, s, t)
                    )
                    + FOUR_THIRDS * g12 * YuYl * delta[p][r]* (
                    + WC6I(CuuI, s, w, w, t) + WC6I(CuuI, w, t, s, w)
                    )
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] +=
                        - 0.5 * (
                        + (yeR[v][r] * yuR[s][w] - yeI[v][r] * yuI[s][w]) * WC5(Clequ1I, p, v, w, t)
                        +(yeI[v][r] * yuR[s][w] + yeR[v][r] * yuI[s][w]) * WC5(Clequ1R, p, v, w, t)
                        -(yedagR[p][v] * yudagR[w][t] - yedagI[p][v] * yudagI[w][t]) * WC5(Clequ1I, r, v, w, s)
                        +(yedagI[p][v] * yudagR[w][t] + yedagR[p][v] * yudagI[w][t]) * WC5(Clequ1R, r, v, w, s)
                        )
                        - 2. * (
                        + (yuR[s][v] * yudagR[w][t] - yuI[s][v] * yudagI[w][t]) * WC7I(Clq1I, p, r, v, w)
                        +(yuI[s][v] * yudagR[w][t] + yuR[s][v] * yudagI[w][t]) * WC7R(Clq1R, p, r, v, w)
                        )
                        - 6. * (
                        + (yeR[v][r] * yuR[s][w] - yeI[v][r] * yuI[s][w]) * WC5(Clequ3I, p, v, w, t)
                        +(yeI[v][r] * yuR[s][w] + yeR[v][r] * yuI[s][w]) * WC5(Clequ3R, p, v, w, t)
                        -(yedagR[p][v] * yudagR[w][t] - yedagI[p][v] * yudagI[w][t]) * WC5(Clequ3I, r, v, w, s)
                        +(yedagI[p][v] * yudagR[w][t] + yedagR[p][v] * yudagI[w][t]) * WC5(Clequ3R, r, v, w, s)
                        )
                        -(
                        + (yeR[w][r] * yedagR[p][v] - yeI[w][r] * yedagI[p][v]) * WC7I(CeuI, v, w, s, t)
                        +(yeI[w][r] * yedagR[p][v] + yeR[w][r] * yedagI[p][v]) * WC7R(CeuR, v, w, s, t)
                        )
                        ; //RGE 2
            }
        }
        f[c] *= loop_factor;
        c ++;
    }

    //CldR
    for (d = 0; d < DWC7R; d ++) {
        p = WC7R_indices[d][0];
        r = WC7R_indices[d][1];
        s = WC7R_indices[d][2];
        t = WC7R_indices[d][3];
        f[c] =
                - (yedyeR[p][r] * WC2R(CHdR, s, t) - yedyeI[p][r] * WC2I(CHdI, s, t))
                + 2. * (gdR[s][t] * WC2R(CHl1R, p, r) - gdI[s][t] * WC2I(CHl1I, p, r))
                //RGE 2
                + FOUR_THIRDS * g12 * YhYl * WC2R(CHdR, s, t) * delta[p][r]
                + FOUR_THIRDS * g12 * YhYd * WC2R(CHl1R, p, r) * delta[s][t]
                - 12. * g12 * YdYl * WC7R(CldR, p, r, s, t)
                ; //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] +=
                    + (glR[p][v] * WC7R(CldR, v, r, s, t) - glI[p][v] * WC7I(CldI, v, r, s, t))
                    +(gdR[s][v] * WC7R(CldR, p, r, v, t) - gdI[s][v] * WC7I(CldI, p, r, v, t))
                    +(WC7R(CldR, p, v, s, t) * glR[v][r] - WC7I(CldI, p, v, s, t) * glI[v][r])
                    +(WC7R(CldR, p, r, s, v) * gdR[v][t] - WC7I(CldI, p, r, s, v) * gdI[v][t])
                    //RGE 2
                    + EIGHT_THIRDS * g12 * YdYl * delta[s][t]*(
                    + WC6R(CllR, p, r, w, w) + WC6R(CllR, w, w, p, r)
                    )
                    + FOUR_THIRDS * g12 * YdYl * delta[s][t]*(
                    + WC6R(CllR, p, w, w, r) + WC6R(CllR, w, r, p, w)
                    )
                    + EIGHT_THIRDS * g12 * NC * YdYq * WC7R(Clq1R, p, r, w, w) * delta[s][t]
                    + EIGHT_THIRDS * g12 * NC * YlYq * WC7R(Cqd1R, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * Yd2 * WC7R(CldR, p, r, w, w) * delta[s][t]
                    + EIGHT_THIRDS * g12 * Yl2 * WC7R(CldR, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YuYd * WC7R(CluR, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * YdYe * WC7R(CleR, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * NC * YuYl * WC7R(Cud1R, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * YeYl * WC7R(CedR, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YdYl * delta[p][r]*(
                    + WC6R(CddR, s, t, w, w) + WC6R(CddR, w, w, s, t)
                    )
                    + FOUR_THIRDS * g12 * YdYl * delta[p][r]*(
                    + WC6R(CddR, s, w, w, t) + WC6R(CddR, w, t, s, w)
                    )
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] +=
                        - 0.5 * (
                        + (yeR[v][r] * yddagR[w][t] - yeI[v][r] * yddagI[w][t]) * WC5(CledqR, p, v, s, w)
                        -(yeI[v][r] * yddagR[w][t] + yeR[v][r] * yddagI[w][t]) * WC5(CledqI, p, v, s, w)
                        +(yedagR[p][v] * ydR[s][w] - yedagI[p][v] * ydI[s][w]) * WC5(CledqR, r, v, t, w)
                        +(yedagI[p][v] * ydR[s][w] + yedagR[p][v] * ydI[s][w]) * WC5(CledqI, r, v, t, w)
                        )
                        - 2. * (
                        + (ydR[s][v] * yddagR[w][t] - ydI[s][v] * yddagI[w][t]) * WC7R(Clq1R, p, r, v, w)
                        -(ydI[s][v] * yddagR[w][t] + ydR[s][v] * yddagI[w][t]) * WC7I(Clq1I, p, r, v, w)
                        )
                        - (
                        + (yeR[w][r] * yedagR[p][v] - yeI[w][r] * yedagI[p][v]) * WC7R(CedR, v, w, s, t)
                        -(yeI[w][r] * yedagR[p][v] + yeR[w][r] * yedagI[p][v]) * WC7I(CedI, v, w, s, t)
                        )

                        ; //RGE 2
            }
        }
        f[c] *= loop_factor;
        c ++;
    }
    //CldI
    for (d = 0; d < DWC7I; d ++) {
        p = WC7I_indices[d][0];
        r = WC7I_indices[d][1];
        s = WC7I_indices[d][2];
        t = WC7I_indices[d][3];
        f[c] =
                - (yedyeI[p][r] * WC2R(CHdR, s, t) + yedyeR[p][r] * WC2I(CHdI, s, t))
                + 2. * (gdI[s][t] * WC2R(CHl1R, p, r) + gdR[s][t] * WC2I(CHl1I, p, r))
                //RGE 2
                + FOUR_THIRDS * g12 * YhYl * WC2I(CHdI, s, t) * delta[p][r]
                + FOUR_THIRDS * g12 * YhYd * WC2I(CHl1I, p, r) * delta[s][t]
                - 12. * g12 * YdYl * WC7I(CldI, p, r, s, t)
                ; //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] +=
                    + (glI[p][v] * WC7R(CldR, v, r, s, t) + glR[p][v] * WC7I(CldI, v, r, s, t))
                    +(gdI[s][v] * WC7R(CldR, p, r, v, t) + gdR[s][v] * WC7I(CldI, p, r, v, t))
                    +(WC7R(CldR, p, v, s, t) * glI[v][r] + WC7I(CldI, p, v, s, t) * glR[v][r])
                    +(WC7R(CldR, p, r, s, v) * gdI[v][t] + WC7I(CldI, p, r, s, v) * gdR[v][t])
                    //RGE 2
                    + EIGHT_THIRDS * g12 * YdYl * delta[s][t]*(
                    + WC6I(CllI, p, r, w, w) + WC6I(CllI, w, w, p, r)
                    )
                    + FOUR_THIRDS * g12 * YdYl * delta[s][t]*(
                    + WC6I(CllI, p, w, w, r) + WC6I(CllI, w, r, p, w)
                    )
                    + EIGHT_THIRDS * g12 * NC * YdYq * WC7I(Clq1I, p, r, w, w) * delta[s][t]
                    + EIGHT_THIRDS * g12 * NC * YlYq * WC7I(Cqd1I, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * Yd2 * WC7I(CldI, p, r, w, w) * delta[s][t]
                    + EIGHT_THIRDS * g12 * Yl2 * WC7I(CldI, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YuYd * WC7I(CluI, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * YdYe * WC7I(CleI, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * NC * YuYl * WC7I(Cud1I, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * YeYl * WC7I(CedI, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YdYl * delta[p][r]*(
                    + WC6I(CddI, s, t, w, w) + WC6I(CddI, w, w, s, t)
                    )
                    + FOUR_THIRDS * g12 * YdYl * delta[p][r]*(
                    + WC6I(CddI, s, w, w, t) + WC6I(CddI, w, t, s, w)
                    )
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] +=
                        - 0.5 * (
                        + (yeR[v][r] * yddagR[w][t] - yeI[v][r] * yddagI[w][t]) * WC5(CledqI, p, v, s, w)
                        +(yeI[v][r] * yddagR[w][t] + yeR[v][r] * yddagI[w][t]) * WC5(CledqR, p, v, s, w)
                        -(yedagR[p][v] * ydR[s][w] - yedagI[p][v] * ydI[s][w]) * WC5(CledqI, r, v, t, w)
                        +(yedagI[p][v] * ydR[s][w] + yedagR[p][v] * ydI[s][w]) * WC5(CledqR, r, v, t, w)
                        )
                        - 2. * (
                        + (ydR[s][v] * yddagR[w][t] - ydI[s][v] * yddagI[w][t]) * WC7I(Clq1I, p, r, v, w)
                        +(ydI[s][v] * yddagR[w][t] + ydR[s][v] * yddagI[w][t]) * WC7R(Clq1R, p, r, v, w)
                        )
                        - (
                        + (yeR[w][r] * yedagR[p][v] - yeI[w][r] * yedagI[p][v]) * WC7I(CedI, v, w, s, t)
                        +(yeI[w][r] * yedagR[p][v] + yeR[w][r] * yedagI[p][v]) * WC7R(CedR, v, w, s, t)
                        )

                        ; //RGE 2
            }
        }
        f[c] *= loop_factor;
        c ++;
    }


    //CqeR
    for (d = 0; d < DWC7R; d ++) {
        p = WC7R_indices[d][0];
        r = WC7R_indices[d][1];
        s = WC7R_indices[d][2];
        t = WC7R_indices[d][3];
        f[c] = (
                + (yudyuR[p][r] - yddydR[p][r]) * WC2R(CHeR, s, t)
                -(yudyuI[p][r] - yddydI[p][r]) * WC2I(CHeI, s, t)
                )
                + 2. * (geR[s][t] * WC2R(CHq1R, p, r) - geI[s][t] * WC2I(CHq1I, p, r))
                //RGE 2
                + FOUR_THIRDS * g12 * YhYq * WC2R(CHeR, s, t) * delta[p][r]
                + FOUR_THIRDS * g12 * YhYe * WC2R(CHq1R, p, r) * delta[s][t]
                - 12. * YeYq * g12 * WC7R(CqeR, p, r, s, t)
                ; //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] +=
                    + (gqR[p][v] * WC7R(CqeR, v, r, s, t) - gqI[p][v] * WC7I(CqeI, v, r, s, t))
                    + (geR[s][v] * WC7R(CqeR, p, r, v, t) - geI[s][v] * WC7I(CqeI, p, r, v, t))
                    +(WC7R(CqeR, p, v, s, t) * gqR[v][r] - WC7I(CqeI, p, v, s, t) * gqI[v][r])
                    +(WC7R(CqeR, p, r, s, v) * geR[v][t] - WC7I(CqeI, p, r, s, v) * geI[v][t])
                    //RGE 2
                    + EIGHT_THIRDS * g12 * NC * YeYq * delta[s][t]*(
                    + WC6R(Cqq1R, p, r, w, w) + WC6R(Cqq1R, w, w, p, r)
                    )
                    + FOUR_THIRDS * g12 * YeYq * delta[s][t]*(
                    + WC6R(Cqq1R, p, w, w, r) + WC6R(Cqq1R, w, r, p, w)
                    )
                    + 4. * g12 * YeYq * delta[s][t]*(
                    + WC6R(Cqq3R, p, w, w, r) + WC6R(Cqq3R, w, r, p, w)
                    )
                    + EIGHT_THIRDS * g12 * YeYl * WC7R(Clq1R, w, w, p, r) * delta[s][t]
                    + EIGHT_THIRDS * g12 * YlYq * WC7R(CleR, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * Ye2 * WC7R(CqeR, p, r, w, w) * delta[s][t]
                    + EIGHT_THIRDS * g12 * NC * Yq2 * WC7R(CqeR, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YuYe * WC7R(Cqu1R, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * NC * YdYe * WC7R(Cqd1R, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * NC * YuYq * WC7R(CeuR, s, t, w, w) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YdYq * WC7R(CedR, s, t, w, w) * delta[p][r]
                    + FOUR_THIRDS * g12 * YeYq * delta[p][r]*(
                    + WC8R(CeeR, s, t, w, w) + WC8R(CeeR, w, w, s, t)
                    + WC8R(CeeR, s, w, w, t) + WC8R(CeeR, w, t, s, w)
                    )
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] += - 0.5 * (
                        + (yddagR[p][w] * yeR[s][v] - yddagI[p][w] * yeI[s][v]) * WC5(CledqR, v, t, w, r)
                        -(yddagI[p][w] * yeR[s][v] + yddagR[p][w] * yeI[s][v]) * WC5(CledqI, v, t, w, r)
                        +(yedagR[v][t] * ydR[w][r] - yedagI[v][t] * ydI[w][r]) * WC5(CledqR, v, s, w, p)
                        +(yedagI[v][t] * ydR[w][r] + yedagR[v][t] * ydI[w][r]) * WC5(CledqI, v, s, w, p)
                        )
                        - 2. * (
                        + (yeR[s][v] * yedagR[w][t] - yeI[s][v] * yedagI[w][t]) * WC7R(Clq1R, v, w, p, r)
                        -(yeI[s][v] * yedagR[w][t] + yeR[s][v] * yedagI[w][t]) * WC7I(Clq1I, v, w, p, r)
                        )
                        - 0.5 * (
                        + (yuR[w][r] * yeR[s][v] - yuI[w][r] * yeI[s][v]) * WC5(Clequ1R, v, t, p, w)
                        -(yuI[w][r] * yeR[s][v] + yuR[w][r] * yeI[s][v]) * WC5(Clequ1I, v, t, p, w)
                        +(yedagR[v][t] * yudagR[p][w] - yedagI[v][t] * yudagI[p][w]) * WC5(Clequ1R, v, s, r, w)
                        +(yedagI[v][t] * yudagR[p][w] + yedagR[v][t] * yudagI[p][w]) * WC5(Clequ1I, v, s, r, w)
                        )
                        -(
                        + (ydR[w][r] * yddagR[p][v] - ydI[w][r] * yddagI[p][v]) * WC7R(CedR, s, t, v, w)
                        -(ydI[w][r] * yddagR[p][v] + ydR[w][r] * yddagI[p][v]) * WC7I(CedI, s, t, v, w)
                        )
                        - 6. * (
                        + (yuR[w][r] * yeR[s][v] - yuI[w][r] * yeI[s][v]) * WC5(Clequ3R, v, t, p, w)
                        -(yuI[w][r] * yeR[s][v] + yuR[w][r] * yeI[s][v]) * WC5(Clequ3I, v, t, p, w)
                        +(yedagR[v][t] * yudagR[p][w] - yedagI[v][t] * yudagI[p][w]) * WC5(Clequ3R, v, s, r, w)
                        +(yedagI[v][t] * yudagR[p][w] + yedagR[v][t] * yudagI[p][w]) * WC5(Clequ3I, v, s, r, w)
                        )
                        -(
                        + (yuR[w][r] * yudagR[p][v] - yuI[w][r] * yudagI[p][v]) * WC7R(CeuR, s, t, v, w)
                        -(yuI[w][r] * yudagR[p][v] + yuR[w][r] * yudagI[p][v]) * WC7I(CeuI, s, t, v, w)
                        )
                        ; //RGE 2
            }
        }


        f[c] *= loop_factor;
        c ++;
    }

    //CqeI
    for (d = 0; d < DWC7I; d ++) {
        p = WC7I_indices[d][0];
        r = WC7I_indices[d][1];
        s = WC7I_indices[d][2];
        t = WC7I_indices[d][3];
        f[c] = (
                + (yudyuI[p][r] - yddydI[p][r]) * WC2R(CHeR, s, t)
                +(yudyuR[p][r] - yddydR[p][r]) * WC2I(CHeI, s, t)
                )
                + 2. * (geI[s][t] * WC2R(CHq1R, p, r) + geR[s][t] * WC2I(CHq1I, p, r))
                //RGE 2
                + FOUR_THIRDS * g12 * YhYq * WC2I(CHeI, s, t) * delta[p][r]
                + FOUR_THIRDS * g12 * YhYe * WC2I(CHq1I, p, r) * delta[s][t]
                - 12. * YeYq * g12 * WC7I(CqeI, p, r, s, t)
                ; //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] +=
                    + (gqI[p][v] * WC7R(CqeR, v, r, s, t) + gqR[p][v] * WC7I(CqeI, v, r, s, t))
                    + (geI[s][v] * WC7R(CqeR, p, r, v, t) + geR[s][v] * WC7I(CqeI, p, r, v, t))
                    +(WC7R(CqeR, p, v, s, t) * gqI[v][r] + WC7I(CqeI, p, v, s, t) * gqR[v][r])
                    +(WC7R(CqeR, p, r, s, v) * geI[v][t] + WC7I(CqeI, p, r, s, v) * geR[v][t])
                    //RGE 2
                    + EIGHT_THIRDS * g12 * NC * YeYq * delta[s][t]*(
                    + WC6I(Cqq1I, p, r, w, w) + WC6I(Cqq1I, w, w, p, r)
                    )
                    + FOUR_THIRDS * g12 * YeYq * delta[s][t]*(
                    + WC6I(Cqq1I, p, w, w, r) + WC6I(Cqq1I, w, r, p, w)
                    )
                    + 4. * g12 * YeYq * delta[s][t]*(
                    + WC6I(Cqq3I, p, w, w, r) + WC6I(Cqq3I, w, r, p, w)
                    )
                    + EIGHT_THIRDS * g12 * YeYl * WC7I(Clq1I, w, w, p, r) * delta[s][t]
                    + EIGHT_THIRDS * g12 * YlYq * WC7I(CleI, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * Ye2 * WC7I(CqeI, p, r, w, w) * delta[s][t]
                    + EIGHT_THIRDS * g12 * NC * Yq2 * WC7I(CqeI, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YuYe * WC7I(Cqu1I, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * NC * YdYe * WC7I(Cqd1I, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * NC * YuYq * WC7I(CeuI, s, t, w, w) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YdYq * WC7I(CedI, s, t, w, w) * delta[p][r]
                    + FOUR_THIRDS * g12 * YeYq * delta[p][r]*(
                    + WC8I(CeeI, s, t, w, w) + WC8I(CeeI, w, w, s, t)
                    + WC8I(CeeI, s, w, w, t) + WC8I(CeeI, w, t, s, w)
                    )
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] += - 0.5 * (
                        + (yddagR[p][w] * yeR[s][v] - yddagI[p][w] * yeI[s][v]) * WC5(CledqI, v, t, w, r)
                        +(yddagI[p][w] * yeR[s][v] + yddagR[p][w] * yeI[s][v]) * WC5(CledqR, v, t, w, r)
                        -(yedagR[v][t] * ydR[w][r] - yedagI[v][t] * ydI[w][r]) * WC5(CledqI, v, s, w, p)
                        +(yedagI[v][t] * ydR[w][r] + yedagR[v][t] * ydI[w][r]) * WC5(CledqR, v, s, w, p)
                        )
                        - 2. * (
                        + (yeR[s][v] * yedagR[w][t] - yeI[s][v] * yedagI[w][t]) * WC7I(Clq1I, v, w, p, r)
                        +(yeI[s][v] * yedagR[w][t] + yeR[s][v] * yedagI[w][t]) * WC7R(Clq1R, v, w, p, r)
                        )
                        - 0.5 * (
                        + (yuR[w][r] * yeR[s][v] - yuI[w][r] * yeI[s][v]) * WC5(Clequ1I, v, t, p, w)
                        +(yuI[w][r] * yeR[s][v] + yuR[w][r] * yeI[s][v]) * WC5(Clequ1R, v, t, p, w)
                        -(yedagR[v][t] * yudagR[p][w] - yedagI[v][t] * yudagI[p][w]) * WC5(Clequ1I, v, s, r, w)
                        +(yedagI[v][t] * yudagR[p][w] + yedagR[v][t] * yudagI[p][w]) * WC5(Clequ1R, v, s, r, w)
                        )
                        -(
                        + (ydR[w][r] * yddagR[p][v] - ydI[w][r] * yddagI[p][v]) * WC7I(CedI, s, t, v, w)
                        +(ydI[w][r] * yddagR[p][v] + ydR[w][r] * yddagI[p][v]) * WC7R(CedR, s, t, v, w)
                        )
                        - 6. * (
                        + (yuR[w][r] * yeR[s][v] - yuI[w][r] * yeI[s][v]) * WC5(Clequ3I, v, t, p, w)
                        +(yuI[w][r] * yeR[s][v] + yuR[w][r] * yeI[s][v]) * WC5(Clequ3R, v, t, p, w)
                        -(yedagR[v][t] * yudagR[p][w] - yedagI[v][t] * yudagI[p][w]) * WC5(Clequ3I, v, s, r, w)
                        +(yedagI[v][t] * yudagR[p][w] + yedagR[v][t] * yudagI[p][w]) * WC5(Clequ3R, v, s, r, w)
                        )
                        -(
                        + (yuR[w][r] * yudagR[p][v] - yuI[w][r] * yudagI[p][v]) * WC7I(CeuI, s, t, v, w)
                        +(yuI[w][r] * yudagR[p][v] + yuR[w][r] * yudagI[p][v]) * WC7R(CeuR, s, t, v, w)
                        )
                        ; //RGE 2
            }
        }
        f[c] *= loop_factor;
        c ++;
    }

    //Cqu1R
    for (d = 0; d < DWC7R; d ++) {
        p = WC7R_indices[d][0];
        r = WC7R_indices[d][1];
        s = WC7R_indices[d][2];
        t = WC7R_indices[d][3];
        f[c] =
                + (1. / NC)*(yuR[s][r] * xiuR[p][t] - yuI[s][r] * xiuI[p][t])
                +(1. / NC)*(yudagR[p][t] * xiuR[r][s] + yudagI[p][t] * xiuI[r][s])
                +(
                + (yudyuR[p][r] - yddydR[p][r]) * WC2R(CHuR, s, t)
                -(yudyuI[p][r] - yddydI[p][r]) * WC2I(CHuI, s, t)
                )
                - 2. * (guR[s][t] * WC2R(CHq1R, p, r) - guI[s][t] * WC2I(CHq1I, p, r))
                //RGE 2
                + FOUR_THIRDS * g12 * YhYq * WC2R(CHuR, s, t) * delta[p][r]
                + FOUR_THIRDS * g12 * YhYu * WC2R(CHq1R, p, r) * delta[s][t]
                - 12. * YuYq * g12 * WC7R(Cqu1R, p, r, s, t)
                - 3. * ((NC2 - 1.) / NC2) * g32 * WC7R(Cqu8R, p, r, s, t)
                ; //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] +=
                    + (gqR[p][v] * WC7R(Cqu1R, v, r, s, t) - gqI[p][v] * WC7I(Cqu1I, v, r, s, t))
                    +(guR[s][v] * WC7R(Cqu1R, p, r, v, t) - guI[s][v] * WC7I(Cqu1I, p, r, v, t))
                    +(WC7R(Cqu1R, p, v, s, t) * gqR[v][r] - WC7I(Cqu1I, p, v, s, t) * gqI[v][r])
                    +(WC7R(Cqu1R, p, r, s, v) * guR[v][t] - WC7I(Cqu1I, p, r, s, v) * guI[v][t])
                    //RGE 2
                    + EIGHT_THIRDS * g12 * NC * YuYq * delta[s][t]*
                    (+ WC6R(Cqq1R, p, r, w, w) + WC6R(Cqq1R, w, w, p, r)
                    )
                    + FOUR_THIRDS * g12 * YuYq * delta[s][t]*
                    (+ WC6R(Cqq1R, p, w, w, r) + WC6R(Cqq1R, w, r, p, w)
                    )
                    + 4. * g12 * YuYq * delta[s][t]* (
                    + WC6R(Cqq3R, p, w, w, r) + WC6R(Cqq3R, w, r, p, w)
                    )
                    + EIGHT_THIRDS * g12 * YuYl * WC7R(Clq1R, w, w, p, r) * delta[s][t]
                    + FOUR_THIRDS * g12 * YuYe * WC7R(CqeR, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * NC * YuYd * WC7R(Cqd1R, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * NC * Yu2 * WC7R(Cqu1R, p, r, w, w) * delta[s][t]

                    + EIGHT_THIRDS * g12 * NC * Yq2 * WC7R(Cqu1R, w, w, s, t) * delta[p][r]
                    + EIGHT_THIRDS * g12 * YlYq * WC7R(CluR, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * YeYq * WC7R(CeuR, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YdYq * WC7R(Cud1R, s, t, w, w) * delta[p][r]

                    + FOUR_THIRDS * g12 * NC * YuYq * delta[p][r]*(
                    + WC6R(CuuR, s, t, w, w) + WC6R(CuuR, w, w, s, t)
                    )
                    + FOUR_THIRDS * g12 * YuYq * delta[p][r]*(
                    + WC6R(CuuR, s, w, w, t) + WC6R(CuuR, w, t, s, w)
                    )
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] += + (1. / NC)*(
                        + (yudagR[p][w] * yuR[s][v] - yudagI[p][w] * yuI[s][v]) * WC7R(Cqu1R, v, r, w, t)
                        -(yudagI[p][w] * yuR[s][v] + yudagR[p][w] * yuI[s][v]) * WC7I(Cqu1I, v, r, w, t)
                        +(yudagR[v][t] * yuR[w][r] - yudagI[v][t] * yuI[w][r]) * WC7R(Cqu1R, p, v, s, w)
                        -(yudagI[v][t] * yuR[w][r] + yudagR[v][t] * yuI[w][r]) * WC7I(Cqu1I, p, v, s, w)
                        +(ydR[w][r] * yuR[s][v] - ydI[w][r] * yuI[s][v]) * WC5(Cquqd1R, p, t, v, w)
                        -(ydI[w][r] * yuR[s][v] + ydR[w][r] * yuI[s][v]) * WC5(Cquqd1I, p, t, v, w)
                        +(yddagR[p][w] * yudagR[v][t] - yddagI[p][w] * yudagI[v][t]) * WC5(Cquqd1R, r, s, v, w)
                        +(yddagI[p][w] * yudagR[v][t] + yddagR[p][w] * yudagI[v][t]) * WC5(Cquqd1I, r, s, v, w)
                        )
                        -(0.5 / NC2)*(
                        + (yudagR[p][w] * yuR[s][v] - yudagI[p][w] * yuI[s][v]) * WC7R(Cqu8R, v, r, w, t)
                        -(yudagI[p][w] * yuR[s][v] + yudagR[p][w] * yuI[s][v]) * WC7I(Cqu8I, v, r, w, t)
                        +(yudagR[v][t] * yuR[w][r] - yudagI[v][t] * yuI[w][r]) * WC7R(Cqu8R, p, v, s, w)
                        -(yudagI[v][t] * yuR[w][r] + yudagR[v][t] * yuI[w][r]) * WC7I(Cqu8I, p, v, s, w)
                        +(ydR[w][r] * yuR[s][v] - ydI[w][r] * yuI[s][v]) * WC5(Cquqd8R, p, t, v, w)
                        -(ydI[w][r] * yuR[s][v] + ydR[w][r] * yuI[s][v]) * WC5(Cquqd8I, p, t, v, w)
                        +(yddagR[p][w] * yudagR[v][t] - yddagI[p][w] * yudagI[v][t]) * WC5(Cquqd8R, r, s, v, w)
                        +(yddagI[p][w] * yudagR[v][t] + yddagR[p][w] * yudagI[v][t]) * WC5(Cquqd8I, r, s, v, w)
                        )
                        -(2. / NC)*(
                        + (yudagR[v][t] * yuR[s][w] - yudagI[v][t] * yuI[s][w]) * WC6R(Cqq1R, p, v, w, r)
                        -(yudagI[v][t] * yuR[s][w] + yudagR[v][t] * yuI[s][w]) * WC6I(Cqq1I, p, v, w, r)
                        +(yudagR[p][v] * yuR[w][r] - yudagI[p][v] * yuI[w][r]) * WC6R(CuuR, v, t, s, w)
                        -(yudagI[p][v] * yuR[w][r] + yudagR[p][v] * yuI[w][r]) * WC6I(CuuI, v, t, s, w)
                        )
                        -(6. / NC)*(
                        + (yudagR[v][t] * yuR[s][w] - yudagI[v][t] * yuI[s][w]) * WC6R(Cqq3R, p, v, w, r)
                        -(yudagI[v][t] * yuR[s][w] + yudagR[v][t] * yuI[s][w]) * WC6I(Cqq3I, p, v, w, r)
                        )
                        + 0.5 * (
                        + (yudagR[p][w] * yuR[s][v] - yudagI[p][w] * yuI[s][v]) * WC7R(Cqu8R, v, r, w, t)
                        -(yudagI[p][w] * yuR[s][v] + yudagR[p][w] * yuI[s][v]) * WC7I(Cqu8I, v, r, w, t)
                        +(yudagR[v][t] * yuR[w][r] - yudagI[v][t] * yuI[w][r]) * WC7R(Cqu8R, p, v, s, w)
                        -(yudagI[v][t] * yuR[w][r] + yudagR[v][t] * yuI[w][r]) * WC7I(Cqu8I, p, v, s, w)
                        )
                        + 0.5 * (
                        + (yuR[s][v] * ydR[w][r] - yuI[s][v] * ydI[w][r]) * WC5(Cquqd1R, v, t, p, w)
                        -(yuI[s][v] * ydR[w][r] + yuR[s][v] * ydI[w][r]) * WC5(Cquqd1I, v, t, p, w)
                        +(yddagR[p][w] * yudagR[v][t] - yddagI[p][w] * yudagI[v][t]) * WC5(Cquqd1R, v, s, r, w)
                        +(yddagI[p][w] * yudagR[v][t] + yddagR[p][w] * yudagI[v][t]) * WC5(Cquqd1I, v, s, r, w)
                        +(yuR[s][v] * ydR[w][r] - yuI[s][v] * ydI[w][r]) * WC5(Cquqd8R, p, t, v, w)
                        -(yuI[s][v] * ydR[w][r] + yuR[s][v] * ydI[w][r]) * WC5(Cquqd8I, p, t, v, w)
                        +(yddagR[p][w] * yudagR[v][t] - yddagI[p][w] * yudagI[v][t]) * WC5(Cquqd8R, r, s, v, w)
                        +(yddagI[p][w] * yudagR[v][t] + yddagR[p][w] * yudagI[v][t]) * WC5(Cquqd8I, r, s, v, w)
                        )

                        - 4. * (
                        + (yudagR[w][t] * yuR[s][v] - yudagI[w][t] * yuI[s][v]) * WC6R(Cqq1R, p, r, v, w)
                        -(yudagI[w][t] * yuR[s][v] + yudagR[w][t] * yuI[s][v]) * WC6I(Cqq1I, p, r, v, w)
                        )
                        - 2. * (
                        + (yudagR[p][v] * yuR[w][r] - yudagI[p][v] * yuI[w][r]) * WC6R(CuuR, v, w, s, t)
                        -(yudagI[p][v] * yuR[w][r] + yudagR[p][v] * yuI[w][r]) * WC6I(CuuI, v, w, s, t)
                        )
                        -(
                        + (yddagR[p][v] * ydR[w][r] - yddagI[p][v] * ydI[w][r]) * WC7R(Cud1R, s, t, v, w)
                        -(yddagI[p][v] * ydR[w][r] + yddagR[p][v] * ydI[w][r]) * WC7I(Cud1I, s, t, v, w)
                        )



                        ; //RGE 2
            }

        }
        f[c] *= loop_factor;
        c ++;
    }

    //Cqu1I
    for (d = 0; d < DWC7I; d ++) {
        p = WC7I_indices[d][0];
        r = WC7I_indices[d][1];
        s = WC7I_indices[d][2];
        t = WC7I_indices[d][3];
        f[c] =
                + (1. / NC)*(yuI[s][r] * xiuR[p][t] + yuR[s][r] * xiuI[p][t])
                +(1. / NC)*(yudagI[p][t] * xiuR[r][s] - yudagR[p][t] * xiuI[r][s])
                +(
                + (yudyuI[p][r] - yddydI[p][r]) * WC2R(CHuR, s, t)
                +(yudyuR[p][r] - yddydR[p][r]) * WC2I(CHuI, s, t)
                )
                - 2. * (guI[s][t] * WC2R(CHq1R, p, r) + guR[s][t] * WC2I(CHq1I, p, r))
                //RGE 2
                + FOUR_THIRDS * g12 * YhYq * WC2I(CHuI, s, t) * delta[p][r]
                + FOUR_THIRDS * g12 * YhYu * WC2I(CHq1I, p, r) * delta[s][t]
                - 12. * YuYq * g12 * WC7I(Cqu1I, p, r, s, t)
                - 3. * ((NC2 - 1.) / NC2) * g32 * WC7I(Cqu8I, p, r, s, t)
                ; //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] +=
                    + (gqI[p][v] * WC7R(Cqu1R, v, r, s, t) + gqR[p][v] * WC7I(Cqu1I, v, r, s, t))
                    +(guI[s][v] * WC7R(Cqu1R, p, r, v, t) + guR[s][v] * WC7I(Cqu1I, p, r, v, t))
                    +(WC7R(Cqu1R, p, v, s, t) * gqI[v][r] + WC7I(Cqu1I, p, v, s, t) * gqR[v][r])
                    +(WC7R(Cqu1R, p, r, s, v) * guI[v][t] + WC7I(Cqu1I, p, r, s, v) * guR[v][t])
                    //RGE 2
                    + EIGHT_THIRDS * g12 * NC * YuYq * delta[s][t]*
                    (+ WC6I(Cqq1I, p, r, w, w) + WC6I(Cqq1I, w, w, p, r)
                    )
                    + FOUR_THIRDS * g12 * YuYq * delta[s][t]*
                    (+ WC6I(Cqq1I, p, w, w, r) + WC6I(Cqq1I, w, r, p, w)
                    )
                    + 4. * g12 * YuYq * delta[s][t]* (
                    + WC6I(Cqq3I, p, w, w, r) + WC6I(Cqq3I, w, r, p, w)
                    )
                    + EIGHT_THIRDS * g12 * YuYl * WC7I(Clq1I, w, w, p, r) * delta[s][t]
                    + FOUR_THIRDS * g12 * YuYe * WC7I(CqeI, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * NC * YuYd * WC7I(Cqd1I, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * NC * Yu2 * WC7I(Cqu1I, p, r, w, w) * delta[s][t]

                    + EIGHT_THIRDS * g12 * NC * Yq2 * WC7I(Cqu1I, w, w, s, t) * delta[p][r]
                    + EIGHT_THIRDS * g12 * YlYq * WC7I(CluI, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * YeYq * WC7I(CeuI, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YdYq * WC7I(Cud1I, s, t, w, w) * delta[p][r]

                    + FOUR_THIRDS * g12 * NC * YuYq * delta[p][r]*(
                    + WC6I(CuuI, s, t, w, w) + WC6I(CuuI, w, w, s, t)
                    )
                    + FOUR_THIRDS * g12 * YuYq * delta[p][r]*(
                    + WC6I(CuuI, s, w, w, t) + WC6I(CuuI, w, t, s, w)
                    )
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] += + (1. / NC)*(
                        + (yudagR[p][w] * yuR[s][v] - yudagI[p][w] * yuI[s][v]) * WC7I(Cqu1I, v, r, w, t)
                        +(yudagI[p][w] * yuR[s][v] + yudagR[p][w] * yuI[s][v]) * WC7R(Cqu1R, v, r, w, t)
                        +(yudagR[v][t] * yuR[w][r] - yudagI[v][t] * yuI[w][r]) * WC7I(Cqu1I, p, v, s, w)
                        +(yudagI[v][t] * yuR[w][r] + yudagR[v][t] * yuI[w][r]) * WC7R(Cqu1R, p, v, s, w)
                        +(ydR[w][r] * yuR[s][v] - ydI[w][r] * yuI[s][v]) * WC5(Cquqd1I, p, t, v, w)
                        +(ydI[w][r] * yuR[s][v] + ydR[w][r] * yuI[s][v]) * WC5(Cquqd1R, p, t, v, w)
                        -(yddagR[p][w] * yudagR[v][t] - yddagI[p][w] * yudagI[v][t]) * WC5(Cquqd1I, r, s, v, w)
                        +(yddagI[p][w] * yudagR[v][t] + yddagR[p][w] * yudagI[v][t]) * WC5(Cquqd1R, r, s, v, w)
                        )
                        -(0.5 / NC2)*(
                        + (yudagR[p][w] * yuR[s][v] - yudagI[p][w] * yuI[s][v]) * WC7I(Cqu8I, v, r, w, t)
                        +(yudagI[p][w] * yuR[s][v] + yudagR[p][w] * yuI[s][v]) * WC7R(Cqu8R, v, r, w, t)
                        +(yudagR[v][t] * yuR[w][r] - yudagI[v][t] * yuI[w][r]) * WC7I(Cqu8I, p, v, s, w)
                        +(yudagI[v][t] * yuR[w][r] + yudagR[v][t] * yuI[w][r]) * WC7R(Cqu8R, p, v, s, w)
                        +(ydR[w][r] * yuR[s][v] - ydI[w][r] * yuI[s][v]) * WC5(Cquqd8I, p, t, v, w)
                        +(ydI[w][r] * yuR[s][v] + ydR[w][r] * yuI[s][v]) * WC5(Cquqd8R, p, t, v, w)
                        -(yddagR[p][w] * yudagR[v][t] - yddagI[p][w] * yudagI[v][t]) * WC5(Cquqd8I, r, s, v, w)
                        +(yddagI[p][w] * yudagR[v][t] + yddagR[p][w] * yudagI[v][t]) * WC5(Cquqd8R, r, s, v, w)
                        )
                        -(2. / NC)*(
                        + (yudagR[v][t] * yuR[s][w] - yudagI[v][t] * yuI[s][w]) * WC6I(Cqq1I, p, v, w, r)
                        +(yudagI[v][t] * yuR[s][w] + yudagR[v][t] * yuI[s][w]) * WC6R(Cqq1R, p, v, w, r)
                        +(yudagR[p][v] * yuR[w][r] - yudagI[p][v] * yuI[w][r]) * WC6I(CuuI, v, t, s, w)
                        +(yudagI[p][v] * yuR[w][r] + yudagR[p][v] * yuI[w][r]) * WC6R(CuuR, v, t, s, w)
                        )
                        -(6. / NC)*(
                        + (yudagR[v][t] * yuR[s][w] - yudagI[v][t] * yuI[s][w]) * WC6I(Cqq3I, p, v, w, r)
                        +(yudagI[v][t] * yuR[s][w] + yudagR[v][t] * yuI[s][w]) * WC6R(Cqq3R, p, v, w, r)
                        )
                        + 0.5 * (
                        + (yudagR[p][w] * yuR[s][v] - yudagI[p][w] * yuI[s][v]) * WC7I(Cqu8I, v, r, w, t)
                        +(yudagI[p][w] * yuR[s][v] + yudagR[p][w] * yuI[s][v]) * WC7R(Cqu8R, v, r, w, t)
                        +(yudagR[v][t] * yuR[w][r] - yudagI[v][t] * yuI[w][r]) * WC7I(Cqu8I, p, v, s, w)
                        +(yudagI[v][t] * yuR[w][r] + yudagR[v][t] * yuI[w][r]) * WC7R(Cqu8R, p, v, s, w)
                        )
                        + 0.5 * (
                        + (yuR[s][v] * ydR[w][r] - yuI[s][v] * ydI[w][r]) * WC5(Cquqd1I, v, t, p, w)
                        +(yuI[s][v] * ydR[w][r] + yuR[s][v] * ydI[w][r]) * WC5(Cquqd1R, v, t, p, w)
                        -(yddagR[p][w] * yudagR[v][t] - yddagI[p][w] * yudagI[v][t]) * WC5(Cquqd1I, v, s, r, w)
                        +(yddagI[p][w] * yudagR[v][t] + yddagR[p][w] * yudagI[v][t]) * WC5(Cquqd1R, v, s, r, w)
                        +(yuR[s][v] * ydR[w][r] - yuI[s][v] * ydI[w][r]) * WC5(Cquqd8I, p, t, v, w)
                        +(yuI[s][v] * ydR[w][r] + yuR[s][v] * ydI[w][r]) * WC5(Cquqd8R, p, t, v, w)
                        -(yddagR[p][w] * yudagR[v][t] - yddagI[p][w] * yudagI[v][t]) * WC5(Cquqd8I, r, s, v, w)
                        +(yddagI[p][w] * yudagR[v][t] + yddagR[p][w] * yudagI[v][t]) * WC5(Cquqd8R, r, s, v, w)
                        )

                        - 4. * (
                        + (yudagR[w][t] * yuR[s][v] - yudagI[w][t] * yuI[s][v]) * WC6I(Cqq1I, p, r, v, w)
                        +(yudagI[w][t] * yuR[s][v] + yudagR[w][t] * yuI[s][v]) * WC6R(Cqq1R, p, r, v, w)
                        )
                        - 2. * (
                        + (yudagR[p][v] * yuR[w][r] - yudagI[p][v] * yuI[w][r]) * WC6I(CuuI, v, w, s, t)
                        +(yudagI[p][v] * yuR[w][r] + yudagR[p][v] * yuI[w][r]) * WC6R(CuuR, v, w, s, t)
                        )
                        -(
                        + (yddagR[p][v] * ydR[w][r] - yddagI[p][v] * ydI[w][r]) * WC7I(Cud1I, s, t, v, w)
                        +(yddagI[p][v] * ydR[w][r] + yddagR[p][v] * ydI[w][r]) * WC7R(Cud1R, s, t, v, w)
                        )
                        ; //RGE 2
            }
        }
        f[c] *= loop_factor;
        c ++;
    }

    //Cqu8R
    for (d = 0; d < DWC7R; d ++) {
        p = WC7R_indices[d][0];
        r = WC7R_indices[d][1];
        s = WC7R_indices[d][2];
        t = WC7R_indices[d][3];
        f[c] =
                + 2. * (yuR[s][r] * xiuR[p][t] - yuI[s][r] * xiuI[p][t])
                + 2. * (yudagR[p][t] * xiuR[r][s] + yudagI[p][t] * xiuI[r][s])
                //RGE 2
                -(12. * YuYq * g12 + 6. * (NC - (2. / NC)) * g32) * WC7R(Cqu8R, p, r, s, t)
                - 12. * g32 * WC7R(Cqu1R, p, r, s, t)
                ; //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] +=
                    + (gqR[p][v] * WC7R(Cqu8R, v, r, s, t) - gqI[p][v] * WC7I(Cqu8I, v, r, s, t))
                    +(guR[s][v] * WC7R(Cqu8R, p, r, v, t) - guI[s][v] * WC7I(Cqu8I, p, r, v, t))
                    +(WC7R(Cqu8R, p, v, s, t) * gqR[v][r] - WC7I(Cqu8I, p, v, s, t) * gqI[v][r])
                    +(WC7R(Cqu8R, p, r, s, v) * guR[v][t] - WC7I(Cqu8I, p, r, s, v) * guI[v][t])
                    //RGE 2
                    + FOUR_THIRDS * g32 * delta[s][t]*(
                    + WC6R(Cqq1R, p, w, w, r) + WC6R(Cqq1R, w, r, p, w)
                    )
                    + 4. * g32 * delta[s][t]*(
                    + WC6R(Cqq3R, p, w, w, r) + WC6R(Cqq3R, w, r, p, w)
                    )
                    + TWO_THIRDS * g32 * WC7R(Cqu8R, p, r, w, w) * delta[s][t]
                    + TWO_THIRDS * g32 * WC7R(Cqd8R, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g32 * WC7R(Cqu8R, w, w, s, t) * delta[p][r]
                    + TWO_THIRDS * g32 * WC7R(Cud8R, s, t, w, w) * delta[p][r]
                    + FOUR_THIRDS * g32 * delta[p][r]*(
                    + WC6R(CuuR, s, w, w, t) + WC6R(CuuR, w, t, s, w)
                    )
                    ; // RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] +=
                        - (1. / NC)*(
                        + (yudagR[p][w] * yuR[s][v] - yudagI[p][w] * yuI[s][v]) * WC7R(Cqu8R, v, r, w, t)
                        -(yudagI[p][w] * yuR[s][v] + yudagR[p][w] * yuI[s][v]) * WC7I(Cqu8I, v, r, w, t)
                        +(yudagR[v][t] * yuR[w][r] - yudagI[v][t] * yuI[w][r]) * WC7R(Cqu8R, p, v, s, w)
                        -(yudagI[v][t] * yuR[w][r] + yudagR[v][t] * yuI[w][r]) * WC7I(Cqu8I, p, v, s, w)
                        +(ydR[w][r] * yuR[s][v] - ydI[w][r] * yuI[s][v]) * WC5(Cquqd8R, p, t, v, w)
                        -(ydI[w][r] * yuR[s][v] + ydR[w][r] * yuI[s][v]) * WC5(Cquqd8I, p, t, v, w)
                        +(yddagR[p][w] * yudagR[v][t] - yddagI[p][w] * yudagI[v][t]) * WC5(Cquqd8R, r, s, v, w)
                        +(yddagI[p][w] * yudagR[v][t] + yddagR[p][w] * yudagI[v][t]) * WC5(Cquqd8I, r, s, v, w)
                        )

                        + 2. * (
                        (
                        + (ydR[w][r] * yuR[s][v] - ydI[w][r] * yuI[s][v]) * WC5(Cquqd1R, p, t, v, w)
                        -(ydI[w][r] * yuR[s][v] + ydR[w][r] * yuI[s][v]) * WC5(Cquqd1I, p, t, v, w)
                        +(yddagR[p][w] * yudagR[v][t] - yddagI[p][w] * yudagI[v][t]) * WC5(Cquqd1R, r, s, v, w)
                        +(yddagI[p][w] * yudagR[v][t] + yddagR[p][w] * yudagI[v][t]) * WC5(Cquqd1I, r, s, v, w)
                        )
                        + 0.25 * (
                        + (ydR[w][r] * yuR[s][v] - ydI[w][r] * yuI[s][v]) * WC5(Cquqd8R, v, t, p, w)
                        -(ydI[w][r] * yuR[s][v] + ydR[w][r] * yuI[s][v]) * WC5(Cquqd8I, v, t, p, w)
                        +(yddagR[p][w] * yudagR[v][t] - yddagI[p][w] * yudagI[v][t]) * WC5(Cquqd8R, v, s, r, w)
                        +(yddagI[p][w] * yudagR[v][t] + yddagR[p][w] * yudagI[v][t]) * WC5(Cquqd8I, v, s, r, w)
                        )
                        )

                        - 2. * (
                        + 2. * (
                        + (yudagR[v][t] * yuR[s][w] - yudagI[v][t] * yuI[s][w]) * WC6R(Cqq1R, p, v, w, r)
                        -(yudagI[v][t] * yuR[s][w] + yudagR[v][t] * yuI[s][w]) * WC6I(Cqq1I, p, v, w, r)
                        )
                        -(
                        + (yudagR[p][w] * yuR[s][v] - yudagI[p][w] * yuI[s][v]) * WC7R(Cqu1R, v, r, w, t)
                        -(yudagI[p][w] * yuR[s][v] + yudagR[p][w] * yuI[s][v]) * WC7I(Cqu1I, v, r, w, t)
                        )
                        -(
                        + (yudagR[v][t] * yuR[w][r] - yudagI[v][t] * yuI[w][r]) * WC7R(Cqu1R, p, v, s, w)
                        -(yudagI[v][t] * yuR[w][r] + yudagR[v][t] * yuI[w][r]) * WC7I(Cqu1I, p, v, s, w)
                        )
                        )

                        - 4. * (
                        + (yudagR[p][v] * yuR[w][r] - yudagI[p][v] * yuI[w][r]) * WC6R(CuuR, v, t, s, w)
                        -(yudagI[p][v] * yuR[w][r] + yudagR[p][v] * yuI[w][r]) * WC6I(CuuI, v, t, s, w)
                        )
                        - 12. * (
                        + (yudagR[v][t] * yuR[s][w] - yudagI[v][t] * yuI[s][w]) * WC6R(Cqq3R, p, v, w, r)
                        -(yudagI[v][t] * yuR[s][w] + yudagR[v][t] * yuI[s][w]) * WC6I(Cqq3I, p, v, w, r)
                        )
                        -(
                        + (yddagR[p][v] * ydR[w][r] - yddagI[p][v] * ydI[w][r]) * WC7R(Cud8R, s, t, v, w)
                        -(yddagI[p][v] * ydR[w][r] + yddagR[p][v] * ydI[w][r]) * WC7I(Cud8I, s, t, v, w)
                        )

                        ; //RGE 2
            }
        }
        f[c] *= loop_factor;
        c ++;
    }

    //Cqu8I
    for (d = 0; d < DWC7I; d ++) {
        p = WC7I_indices[d][0];
        r = WC7I_indices[d][1];
        s = WC7I_indices[d][2];
        t = WC7I_indices[d][3];
        f[c] =
                + 2. * (yuI[s][r] * xiuR[p][t] + yuR[s][r] * xiuI[p][t])
                + 2. * (yudagI[p][t] * xiuR[r][s] - yudagR[p][t] * xiuI[r][s])
                //RGE 2
                -(12. * YuYq * g12 + 6. * (NC - (2. / NC)) * g32) * WC7I(Cqu8I, p, r, s, t)
                - 12. * g32 * WC7I(Cqu1I, p, r, s, t)
                ; //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] +=
                    + (gqI[p][v] * WC7R(Cqu8R, v, r, s, t) + gqR[p][v] * WC7I(Cqu8I, v, r, s, t))
                    +(guI[s][v] * WC7R(Cqu8R, p, r, v, t) + guR[s][v] * WC7I(Cqu8I, p, r, v, t))
                    +(WC7R(Cqu8R, p, v, s, t) * gqI[v][r] + WC7I(Cqu8I, p, v, s, t) * gqR[v][r])
                    +(WC7R(Cqu8R, p, r, s, v) * guI[v][t] + WC7I(Cqu8I, p, r, s, v) * guR[v][t])
                    //RGE 2
                    + FOUR_THIRDS * g32 * delta[s][t]*(
                    + WC6I(Cqq1I, p, w, w, r) + WC6I(Cqq1I, w, r, p, w)
                    )
                    + 4. * g32 * delta[s][t]*(
                    + WC6I(Cqq3I, p, w, w, r) + WC6I(Cqq3I, w, r, p, w)
                    )
                    + TWO_THIRDS * g32 * WC7I(Cqu8I, p, r, w, w) * delta[s][t]
                    + TWO_THIRDS * g32 * WC7I(Cqd8I, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g32 * WC7I(Cqu8I, w, w, s, t) * delta[p][r]
                    + TWO_THIRDS * g32 * WC7I(Cud8I, s, t, w, w) * delta[p][r]
                    + FOUR_THIRDS * g32 * delta[p][r]*(
                    + WC6I(CuuI, s, w, w, t) + WC6I(CuuI, w, t, s, w)
                    )
                    ; // RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] += - (1. / NC)*(
                        + (yudagR[p][w] * yuR[s][v] - yudagI[p][w] * yuI[s][v]) * WC7I(Cqu8I, v, r, w, t)
                        +(yudagI[p][w] * yuR[s][v] + yudagR[p][w] * yuI[s][v]) * WC7R(Cqu8R, v, r, w, t)
                        +(yudagR[v][t] * yuR[w][r] - yudagI[v][t] * yuI[w][r]) * WC7I(Cqu8I, p, v, s, w)
                        +(yudagI[v][t] * yuR[w][r] + yudagR[v][t] * yuI[w][r]) * WC7R(Cqu8R, p, v, s, w)
                        +(ydR[w][r] * yuR[s][v] - ydI[w][r] * yuI[s][v]) * WC5(Cquqd8I, p, t, v, w)
                        +(ydI[w][r] * yuR[s][v] + ydR[w][r] * yuI[s][v]) * WC5(Cquqd8R, p, t, v, w)
                        -(yddagR[p][w] * yudagR[v][t] - yddagI[p][w] * yudagI[v][t]) * WC5(Cquqd8I, r, s, v, w)
                        +(yddagI[p][w] * yudagR[v][t] + yddagR[p][w] * yudagI[v][t]) * WC5(Cquqd8R, r, s, v, w)
                        )

                        + 2. * (
                        (
                        + (ydR[w][r] * yuR[s][v] - ydI[w][r] * yuI[s][v]) * WC5(Cquqd1I, p, t, v, w)
                        +(ydI[w][r] * yuR[s][v] + ydR[w][r] * yuI[s][v]) * WC5(Cquqd1R, p, t, v, w)
                        -(yddagR[p][w] * yudagR[v][t] - yddagI[p][w] * yudagI[v][t]) * WC5(Cquqd1I, r, s, v, w)
                        +(yddagI[p][w] * yudagR[v][t] + yddagR[p][w] * yudagI[v][t]) * WC5(Cquqd1R, r, s, v, w)
                        )
                        + 0.25 * (
                        + (ydR[w][r] * yuR[s][v] - ydI[w][r] * yuI[s][v]) * WC5(Cquqd8I, v, t, p, w)
                        +(ydI[w][r] * yuR[s][v] + ydR[w][r] * yuI[s][v]) * WC5(Cquqd8R, v, t, p, w)
                        -(yddagR[p][w] * yudagR[v][t] - yddagI[p][w] * yudagI[v][t]) * WC5(Cquqd8I, v, s, r, w)
                        +(yddagI[p][w] * yudagR[v][t] + yddagR[p][w] * yudagI[v][t]) * WC5(Cquqd8R, v, s, r, w)
                        )
                        )

                        - 2. * (
                        + 2. * (
                        + (yudagR[v][t] * yuR[s][w] - yudagI[v][t] * yuI[s][w]) * WC6I(Cqq1I, p, v, w, r)
                        +(yudagI[v][t] * yuR[s][w] + yudagR[v][t] * yuI[s][w]) * WC6R(Cqq1R, p, v, w, r)
                        )
                        -(
                        + (yudagR[p][w] * yuR[s][v] - yudagI[p][w] * yuI[s][v]) * WC7I(Cqu1I, v, r, w, t)
                        +(yudagI[p][w] * yuR[s][v] + yudagR[p][w] * yuI[s][v]) * WC7R(Cqu1R, v, r, w, t)
                        )
                        -(
                        + (yudagR[v][t] * yuR[w][r] - yudagI[v][t] * yuI[w][r]) * WC7I(Cqu1I, p, v, s, w)
                        +(yudagI[v][t] * yuR[w][r] + yudagR[v][t] * yuI[w][r]) * WC7R(Cqu1R, p, v, s, w)
                        )
                        )

                        - 4. * (
                        + (yudagR[p][v] * yuR[w][r] - yudagI[p][v] * yuI[w][r]) * WC6I(CuuI, v, t, s, w)
                        +(yudagI[p][v] * yuR[w][r] + yudagR[p][v] * yuI[w][r]) * WC6R(CuuR, v, t, s, w)
                        )
                        - 12. * (
                        + (yudagR[v][t] * yuR[s][w] - yudagI[v][t] * yuI[s][w]) * WC6I(Cqq3I, p, v, w, r)
                        +(yudagI[v][t] * yuR[s][w] + yudagR[v][t] * yuI[s][w]) * WC6R(Cqq3R, p, v, w, r)
                        )
                        -(
                        + (yddagR[p][v] * ydR[w][r] - yddagI[p][v] * ydI[w][r]) * WC7I(Cud8I, s, t, v, w)
                        +(yddagI[p][v] * ydR[w][r] + yddagR[p][v] * ydI[w][r]) * WC7R(Cud8R, s, t, v, w)
                        )

                        ; //RGE 2
            }
        }
        f[c] *= loop_factor;
        c ++;
    }


    //Cqd1R
    for (d = 0; d < DWC7R; d ++) {
        p = WC7R_indices[d][0];
        r = WC7R_indices[d][1];
        s = WC7R_indices[d][2];
        t = WC7R_indices[d][3];
        f[c] =
                + (1. / NC)*(ydR[s][r] * xidR[p][t] - ydI[s][r] * xidI[p][t])
                +(1. / NC)*(yddagR[p][t] * xidR[r][s] + yddagI[p][t] * xidI[r][s])
                +(
                + (yudyuR[p][r] - yddydR[p][r]) * WC2R(CHdR, s, t)
                -(yudyuI[p][r] - yddydI[p][r]) * WC2I(CHdI, s, t)
                )
                + 2. * (gdR[s][t] * WC2R(CHq1R, p, r) - gdI[s][t] * WC2I(CHq1I, p, r))
                //RGE 2
                + FOUR_THIRDS * g12 * YhYq * WC2R(CHdR, s, t) * delta[p][r]
                + FOUR_THIRDS * g12 * YhYd * WC2R(CHq1R, p, r) * delta[s][t]
                - 12. * YdYq * g12 * WC7R(Cqd1R, p, r, s, t)
                - 3. * ((NC2 - 1.) / NC2) * g32 * WC7R(Cqd8R, p, r, s, t)
                ; //RGE 3


        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] +=
                    + (gqR[p][v] * WC7R(Cqd1R, v, r, s, t) - gqI[p][v] * WC7I(Cqd1I, v, r, s, t))
                    +(gdR[s][v] * WC7R(Cqd1R, p, r, v, t) - gdI[s][v] * WC7I(Cqd1I, p, r, v, t))
                    +(WC7R(Cqd1R, p, v, s, t) * gqR[v][r] - WC7I(Cqd1I, p, v, s, t) * gqI[v][r])
                    +(WC7R(Cqd1R, p, r, s, v) * gdR[v][t] - WC7I(Cqd1I, p, r, s, v) * gdI[v][t])
                    //RGE 2
                    + EIGHT_THIRDS * g12 * NC * YdYq * delta[s][t]*(
                    + WC6R(Cqq1R, p, r, w, w) + WC6R(Cqq1R, w, w, p, r)
                    )
                    + FOUR_THIRDS * g12 * YdYq * delta[s][t]*(
                    + WC6R(Cqq1R, p, w, w, r) + WC6R(Cqq1R, w, r, p, w)
                    )
                    + 4. * g12 * YdYq * delta[s][t]*(
                    + WC6R(Cqq3R, p, w, w, r) + WC6R(Cqq3R, w, r, p, w)
                    )
                    + EIGHT_THIRDS * g12 * YdYl * WC7R(Clq1R, w, w, p, r) * delta[s][t]
                    + FOUR_THIRDS * g12 * YdYe * WC7R(CqeR, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * NC * YuYd * WC7R(Cqu1R, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * NC * Yd2 * WC7R(Cqd1R, p, r, w, w) * delta[s][t]
                    + EIGHT_THIRDS * g12 * NC * Yq2 * WC7R(Cqd1R, w, w, s, t) * delta[p][r]
                    + EIGHT_THIRDS * g12 * YlYq * WC7R(CldR, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * YeYq * WC7R(CedR, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YuYq * WC7R(Cud1R, w, w, s, t) * delta[p][r]

                    + FOUR_THIRDS * g12 * NC * YdYq * delta[p][r]*(
                    + WC6R(CddR, s, t, w, w) + WC6R(CddR, w, w, s, t)
                    )
                    + FOUR_THIRDS * g12 * YdYq * delta[p][r]*(
                    + WC6R(CddR, s, w, w, t) + WC6R(CddR, w, t, s, w)
                    )
                    ; //RGE 3


            for (w = 0; w < NG; w ++) {
                f[c] += + (1. / NC)*(
                        + (yddagR[p][w] * ydR[s][v] - yddagI[p][w] * ydI[s][v]) * WC7R(Cqd1R, v, r, w, t)
                        -(yddagI[p][w] * ydR[s][v] + yddagR[p][w] * ydI[s][v]) * WC7I(Cqd1I, v, r, w, t)
                        +(yddagR[v][t] * ydR[w][r] - yddagI[v][t] * ydI[w][r]) * WC7R(Cqd1R, p, v, s, w)
                        -(yddagI[v][t] * ydR[w][r] + yddagR[v][t] * ydI[w][r]) * WC7I(Cqd1I, p, v, s, w)
                        +(yuR[w][r] * ydR[s][v] - yuI[w][r] * ydI[s][v]) * WC5(Cquqd1R, v, w, p, t)
                        -(yuI[w][r] * ydR[s][v] + yuR[w][r] * ydI[s][v]) * WC5(Cquqd1I, v, w, p, t)
                        +(yudagR[p][w] * yddagR[v][t] - yudagI[p][w] * yddagI[v][t]) * WC5(Cquqd1R, v, w, r, s)
                        +(yudagI[p][w] * yddagR[v][t] + yudagR[p][w] * yddagI[v][t]) * WC5(Cquqd1I, v, w, r, s)
                        )

                        -(0.5 / NC2)*(
                        + (yddagR[p][w] * ydR[s][v] - yddagI[p][w] * ydI[s][v]) * WC7R(Cqd8R, v, r, w, t)
                        -(yddagI[p][w] * ydR[s][v] + yddagR[p][w] * ydI[s][v]) * WC7I(Cqd8I, v, r, w, t)
                        +(yddagR[v][t] * ydR[w][r] - yddagI[v][t] * ydI[w][r]) * WC7R(Cqd8R, p, v, s, w)
                        -(yddagI[v][t] * ydR[w][r] + yddagR[v][t] * ydI[w][r]) * WC7I(Cqd8I, p, v, s, w)
                        +(yuR[w][r] * ydR[s][v] - yuI[w][r] * ydI[s][v]) * WC5(Cquqd8R, v, w, p, t)
                        -(yuI[w][r] * ydR[s][v] + yuR[w][r] * ydI[s][v]) * WC5(Cquqd8I, v, w, p, t)
                        +(yudagR[p][w] * yddagR[v][t] - yudagI[p][w] * yddagI[v][t]) * WC5(Cquqd8R, v, w, r, s)
                        +(yudagI[p][w] * yddagR[v][t] + yudagR[p][w] * yddagI[v][t]) * WC5(Cquqd8I, v, w, r, s)
                        )

                        -(2. / NC)*(
                        + (yddagR[v][t] * ydR[s][w] - yddagI[v][t] * ydI[s][w]) * WC6R(Cqq1R, p, v, w, r)
                        -(yddagI[v][t] * ydR[s][w] + yddagR[v][t] * ydI[s][w]) * WC6I(Cqq1I, p, v, w, r)
                        +(yddagR[p][v] * ydR[w][r] - yddagI[p][v] * ydI[w][r]) * WC6R(CddR, v, t, s, w)
                        -(yddagI[p][v] * ydR[w][r] + yddagR[p][v] * ydI[w][r]) * WC6I(CddI, v, t, s, w)
                        )

                        -(6. / NC)*(
                        + (yddagR[v][t] * ydR[s][w] - yddagI[v][t] * ydI[s][w]) * WC6R(Cqq3R, p, v, w, r)
                        -(yddagI[v][t] * ydR[s][w] + yddagR[v][t] * ydI[s][w]) * WC6I(Cqq3I, p, v, w, r)
                        )

                        + 0.5 * (
                        + (yddagR[p][w] * ydR[s][v] - yddagI[p][w] * ydI[s][v]) * WC7R(Cqd8R, v, r, w, t)
                        -(yddagI[p][w] * ydR[s][v] + yddagR[p][w] * ydI[s][v]) * WC7I(Cqd8I, v, r, w, t)
                        +(yddagR[v][t] * ydR[w][r] - yddagI[v][t] * ydI[w][r]) * WC7R(Cqd8R, p, v, s, w)
                        -(yddagI[v][t] * ydR[w][r] + yddagR[v][t] * ydI[w][r]) * WC7I(Cqd8I, p, v, s, w)
                        )

                        + 0.5 * (
                        + (ydR[s][w] * yuR[v][r] - ydI[s][w] * yuI[v][r]) * WC5(Cquqd1R, p, v, w, t)
                        -(ydI[s][w] * yuR[v][r] + ydR[s][w] * yuI[v][r]) * WC5(Cquqd1I, p, v, w, t)
                        +(yudagR[p][v] * yddagR[w][t] - yudagI[p][v] * yddagI[w][t]) * WC5(Cquqd1R, r, v, w, s)
                        +(yudagI[p][v] * yddagR[w][t] + yudagR[p][v] * yddagI[w][t]) * WC5(Cquqd1I, r, v, w, s)
                        +(ydR[s][v] * yuR[w][r] - ydI[s][v] * yuI[w][r]) * WC5(Cquqd8R, v, w, p, t)
                        -(ydI[s][v] * yuR[w][r] + ydR[s][v] * yuI[w][r]) * WC5(Cquqd8I, v, w, p, t)
                        +(yudagR[p][w] * yddagR[v][t] - yudagI[p][w] * yddagI[v][t]) * WC5(Cquqd8R, v, w, r, s)
                        +(yudagI[p][w] * yddagR[v][t] + yudagR[p][w] * yddagI[v][t]) * WC5(Cquqd8I, v, w, r, s)
                        )

                        - 4. * (
                        + (yddagR[w][t] * ydR[s][v] - yddagI[w][t] * ydI[s][v]) * WC6R(Cqq1R, p, r, v, w)
                        -(yddagI[w][t] * ydR[s][v] + yddagR[w][t] * ydI[s][v]) * WC6I(Cqq1I, p, r, v, w)
                        )
                        - 2. * (
                        + (yddagR[p][v] * ydR[w][r] - yddagI[p][v] * ydI[w][r]) * WC6R(CddR, v, w, s, t)
                        -(yddagI[p][v] * ydR[w][r] + yddagR[p][v] * ydI[w][r]) * WC6I(CddI, v, w, s, t)
                        )

                        -(
                        + (yudagR[p][v] * yuR[w][r] - yudagI[p][v] * yuI[w][r]) * WC7R(Cud1R, v, w, s, t)
                        -(yudagI[p][v] * yuR[w][r] + yudagR[p][v] * yuI[w][r]) * WC7I(Cud1I, v, w, s, t)
                        )
                        ; //RGE 2
            }
        }
        f[c] *= loop_factor;
        c ++;
    }


    //Cqd1I
    for (d = 0; d < DWC7I; d ++) {
        p = WC7I_indices[d][0];
        r = WC7I_indices[d][1];
        s = WC7I_indices[d][2];
        t = WC7I_indices[d][3];
        f[c] =
                + (1. / NC)*(ydI[s][r] * xidR[p][t] + ydR[s][r] * xidI[p][t])
                +(1. / NC)*(yddagI[p][t] * xidR[r][s] - yddagR[p][t] * xidI[r][s])
                +(
                + (yudyuI[p][r] - yddydI[p][r]) * WC2R(CHdR, s, t)
                +(yudyuR[p][r] - yddydR[p][r]) * WC2I(CHdI, s, t)
                )
                + 2. * (gdI[s][t] * WC2R(CHq1R, p, r) + gdR[s][t] * WC2I(CHq1I, p, r))
                //RGE 2
                + FOUR_THIRDS * g12 * YhYq * WC2I(CHdI, s, t) * delta[p][r]
                + FOUR_THIRDS * g12 * YhYd * WC2I(CHq1I, p, r) * delta[s][t]
                - 12. * YdYq * g12 * WC7I(Cqd1I, p, r, s, t)
                - 3. * ((NC2 - 1.) / NC2) * g32 * WC7I(Cqd8I, p, r, s, t)
                ; //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] +=
                    + (gqI[p][v] * WC7R(Cqd1R, v, r, s, t) + gqR[p][v] * WC7I(Cqd1I, v, r, s, t))
                    +(gdI[s][v] * WC7R(Cqd1R, p, r, v, t) + gdR[s][v] * WC7I(Cqd1I, p, r, v, t))
                    +(WC7R(Cqd1R, p, v, s, t) * gqI[v][r] + WC7I(Cqd1I, p, v, s, t) * gqR[v][r])
                    +(WC7R(Cqd1R, p, r, s, v) * gdI[v][t] + WC7I(Cqd1I, p, r, s, v) * gdR[v][t])
                    //RGE 2
                    + EIGHT_THIRDS * g12 * NC * YdYq * delta[s][t]*(
                    + WC6I(Cqq1I, p, r, w, w) + WC6I(Cqq1I, w, w, p, r)
                    )
                    + FOUR_THIRDS * g12 * YdYq * delta[s][t]*(
                    + WC6I(Cqq1I, p, w, w, r) + WC6I(Cqq1I, w, r, p, w)
                    )
                    + 4. * g12 * YdYq * delta[s][t]*(
                    + WC6I(Cqq3I, p, w, w, r) + WC6I(Cqq3I, w, r, p, w)
                    )
                    + EIGHT_THIRDS * g12 * YdYl * WC7I(Clq1I, w, w, p, r) * delta[s][t]
                    + FOUR_THIRDS * g12 * YdYe * WC7I(CqeI, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * NC * YuYd * WC7I(Cqu1I, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g12 * NC * Yd2 * WC7I(Cqd1I, p, r, w, w) * delta[s][t]
                    + EIGHT_THIRDS * g12 * NC * Yq2 * WC7I(Cqd1I, w, w, s, t) * delta[p][r]
                    + EIGHT_THIRDS * g12 * YlYq * WC7I(CldI, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * YeYq * WC7I(CedI, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g12 * NC * YuYq * WC7I(Cud1I, w, w, s, t) * delta[p][r]

                    + FOUR_THIRDS * g12 * NC * YdYq * delta[p][r]*(
                    + WC6I(CddI, s, t, w, w) + WC6I(CddI, w, w, s, t)
                    )
                    + FOUR_THIRDS * g12 * YdYq * delta[p][r]*(
                    + WC6I(CddI, s, w, w, t) + WC6I(CddI, w, t, s, w)
                    )
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] += + (1. / NC)*(
                        + (yddagR[p][w] * ydR[s][v] - yddagI[p][w] * ydI[s][v]) * WC7I(Cqd1I, v, r, w, t)
                        +(yddagI[p][w] * ydR[s][v] + yddagR[p][w] * ydI[s][v]) * WC7R(Cqd1R, v, r, w, t)
                        +(yddagR[v][t] * ydR[w][r] - yddagI[v][t] * ydI[w][r]) * WC7I(Cqd1I, p, v, s, w)
                        +(yddagI[v][t] * ydR[w][r] + yddagR[v][t] * ydI[w][r]) * WC7R(Cqd1R, p, v, s, w)
                        +(yuR[w][r] * ydR[s][v] - yuI[w][r] * ydI[s][v]) * WC5(Cquqd1I, v, w, p, t)
                        +(yuI[w][r] * ydR[s][v] + yuR[w][r] * ydI[s][v]) * WC5(Cquqd1R, v, w, p, t)
                        -(yudagR[p][w] * yddagR[v][t] - yudagI[p][w] * yddagI[v][t]) * WC5(Cquqd1I, v, w, r, s)
                        +(yudagI[p][w] * yddagR[v][t] + yudagR[p][w] * yddagI[v][t]) * WC5(Cquqd1R, v, w, r, s)
                        )

                        -(0.5 / NC2)*(
                        + (yddagR[p][w] * ydR[s][v] - yddagI[p][w] * ydI[s][v]) * WC7I(Cqd8I, v, r, w, t)
                        +(yddagI[p][w] * ydR[s][v] + yddagR[p][w] * ydI[s][v]) * WC7R(Cqd8R, v, r, w, t)
                        +(yddagR[v][t] * ydR[w][r] - yddagI[v][t] * ydI[w][r]) * WC7I(Cqd8I, p, v, s, w)
                        +(yddagI[v][t] * ydR[w][r] + yddagR[v][t] * ydI[w][r]) * WC7R(Cqd8R, p, v, s, w)
                        +(yuR[w][r] * ydR[s][v] - yuI[w][r] * ydI[s][v]) * WC5(Cquqd8I, v, w, p, t)
                        +(yuI[w][r] * ydR[s][v] + yuR[w][r] * ydI[s][v]) * WC5(Cquqd8R, v, w, p, t)
                        -(yudagR[p][w] * yddagR[v][t] - yudagI[p][w] * yddagI[v][t]) * WC5(Cquqd8I, v, w, r, s)
                        +(yudagI[p][w] * yddagR[v][t] + yudagR[p][w] * yddagI[v][t]) * WC5(Cquqd8R, v, w, r, s)
                        )

                        -(2. / NC)*(
                        + (yddagR[v][t] * ydR[s][w] - yddagI[v][t] * ydI[s][w]) * WC6I(Cqq1I, p, v, w, r)
                        +(yddagI[v][t] * ydR[s][w] + yddagR[v][t] * ydI[s][w]) * WC6R(Cqq1R, p, v, w, r)
                        +(yddagR[p][v] * ydR[w][r] - yddagI[p][v] * ydI[w][r]) * WC6I(CddI, v, t, s, w)
                        +(yddagI[p][v] * ydR[w][r] + yddagR[p][v] * ydI[w][r]) * WC6R(CddR, v, t, s, w)
                        )

                        -(6. / NC)*(
                        + (yddagR[v][t] * ydR[s][w] - yddagI[v][t] * ydI[s][w]) * WC6I(Cqq3I, p, v, w, r)
                        +(yddagI[v][t] * ydR[s][w] + yddagR[v][t] * ydI[s][w]) * WC6R(Cqq3R, p, v, w, r)
                        )

                        + 0.5 * (
                        + (yddagR[p][w] * ydR[s][v] - yddagI[p][w] * ydI[s][v]) * WC7I(Cqd8I, v, r, w, t)
                        +(yddagI[p][w] * ydR[s][v] + yddagR[p][w] * ydI[s][v]) * WC7R(Cqd8R, v, r, w, t)
                        +(yddagR[v][t] * ydR[w][r] - yddagI[v][t] * ydI[w][r]) * WC7I(Cqd8I, p, v, s, w)
                        +(yddagI[v][t] * ydR[w][r] + yddagR[v][t] * ydI[w][r]) * WC7R(Cqd8R, p, v, s, w)
                        )

                        + 0.5 * (
                        + (ydR[s][w] * yuR[v][r] - ydI[s][w] * yuI[v][r]) * WC5(Cquqd1I, p, v, w, t)
                        +(ydI[s][w] * yuR[v][r] + ydR[s][w] * yuI[v][r]) * WC5(Cquqd1R, p, v, w, t)
                        -(yudagR[p][v] * yddagR[w][t] - yudagI[p][v] * yddagI[w][t]) * WC5(Cquqd1I, r, v, w, s)
                        +(yudagI[p][v] * yddagR[w][t] + yudagR[p][v] * yddagI[w][t]) * WC5(Cquqd1R, r, v, w, s)
                        +(ydR[s][v] * yuR[w][r] - ydI[s][v] * yuI[w][r]) * WC5(Cquqd8I, v, w, p, t)
                        +(ydI[s][v] * yuR[w][r] + ydR[s][v] * yuI[w][r]) * WC5(Cquqd8R, v, w, p, t)
                        -(yudagR[p][w] * yddagR[v][t] - yudagI[p][w] * yddagI[v][t]) * WC5(Cquqd8I, v, w, r, s)
                        +(yudagI[p][w] * yddagR[v][t] + yudagR[p][w] * yddagI[v][t]) * WC5(Cquqd8R, v, w, r, s)
                        )

                        - 4. * (
                        + (yddagR[w][t] * ydR[s][v] - yddagI[w][t] * ydI[s][v]) * WC6I(Cqq1I, p, r, v, w)
                        +(yddagI[w][t] * ydR[s][v] + yddagR[w][t] * ydI[s][v]) * WC6R(Cqq1R, p, r, v, w)
                        )
                        - 2. * (
                        + (yddagR[p][v] * ydR[w][r] - yddagI[p][v] * ydI[w][r]) * WC6I(CddI, v, w, s, t)
                        +(yddagI[p][v] * ydR[w][r] + yddagR[p][v] * ydI[w][r]) * WC6R(CddR, v, w, s, t)
                        )

                        -(
                        + (yudagR[p][v] * yuR[w][r] - yudagI[p][v] * yuI[w][r]) * WC7I(Cud1I, v, w, s, t)
                        +(yudagI[p][v] * yuR[w][r] + yudagR[p][v] * yuI[w][r]) * WC7R(Cud1R, v, w, s, t)
                        )
                        ; //RGE 2
            }
        }
        f[c] *= loop_factor;
        c ++;
    }




    //Cqd8R
    for (d = 0; d < DWC7R; d ++) {
        p = WC7R_indices[d][0];
        r = WC7R_indices[d][1];
        s = WC7R_indices[d][2];
        t = WC7R_indices[d][3];
        f[c] =
                + 2. * (ydR[s][r] * xidR[p][t] - ydI[s][r] * xidI[p][t])
                + 2. * (yddagR[p][t] * xidR[r][s] + yddagI[p][t] * xidI[r][s])
                //RGE 2
                -(12. * YdYq * g12 + 6. * (NC - (2. / NC)) * g32) * WC7R(Cqd8R, p, r, s, t)
                - 12. * g32 * WC7R(Cqd1R, p, r, s, t)
                ; //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] +=
                    + (gqR[p][v] * WC7R(Cqd8R, v, r, s, t) - gqI[p][v] * WC7I(Cqd8I, v, r, s, t))
                    +(gdR[s][v] * WC7R(Cqd8R, p, r, v, t) - gdI[s][v] * WC7I(Cqd8I, p, r, v, t))
                    +(WC7R(Cqd8R, p, v, s, t) * gqR[v][r] - WC7I(Cqd8I, p, v, s, t) * gqI[v][r])
                    +(WC7R(Cqd8R, p, r, s, v) * gdR[v][t] - WC7I(Cqd8I, p, r, s, v) * gdI[v][t])
                    //RGE 2
                    + FOUR_THIRDS * g32 * delta[s][t]*(
                    + WC6R(Cqq1R, p, w, w, r) + WC6R(Cqq1R, w, r, p, w)
                    )
                    + 4. * g32 * delta[s][t]*(
                    + WC6R(Cqq3R, p, w, w, r) + WC6R(Cqq3R, w, r, p, w)
                    )
                    + TWO_THIRDS * g32 * WC7R(Cqu8R, p, r, w, w) * delta[s][t]
                    + TWO_THIRDS * g32 * WC7R(Cqd8R, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g32 * WC7R(Cqd8R, w, w, s, t) * delta[p][r]
                    + TWO_THIRDS * g32 * WC7R(Cud8R, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g32 * delta[p][r]*(
                    + WC6R(CddR, s, w, w, t) + WC6R(CddR, w, t, s, w)
                    )
                    ; //RGE 3
            for (w = 0; w < NG; w ++) {
                f[c] += - (1. / NC)*(
                        + (yddagR[p][w] * ydR[s][v] - yddagI[p][w] * ydI[s][v]) * WC7R(Cqd8R, v, r, w, t)
                        -(yddagI[p][w] * ydR[s][v] + yddagR[p][w] * ydI[s][v]) * WC7I(Cqd8I, v, r, w, t)
                        +(yddagR[v][t] * ydR[w][r] - yddagI[v][t] * ydI[w][r]) * WC7R(Cqd8R, p, v, s, w)
                        -(yddagI[v][t] * ydR[w][r] + yddagR[v][t] * ydI[w][r]) * WC7I(Cqd8I, p, v, s, w)
                        +(yuR[w][r] * ydR[s][v] - yuI[w][r] * ydI[s][v]) * WC5(Cquqd8R, v, w, p, t)
                        -(yuI[w][r] * ydR[s][v] + yuR[w][r] * ydI[s][v]) * WC5(Cquqd8I, v, w, p, t)
                        +(yudagR[p][w] * yddagR[v][t] - yudagI[p][w] * yddagI[v][t]) * WC5(Cquqd8R, v, w, r, s)
                        +(yudagI[p][w] * yddagR[v][t] + yudagR[p][w] * yddagI[v][t]) * WC5(Cquqd8I, v, w, r, s)
                        )
                        + 2. * (
                        + (
                        + (yuR[w][r] * ydR[s][v] - yuI[w][r] * ydI[s][v]) * WC5(Cquqd1R, v, w, p, t)
                        -(yuI[w][r] * ydR[s][v] + yuR[w][r] * ydI[s][v]) * WC5(Cquqd1I, v, w, p, t)
                        +(yudagR[p][w] * yddagR[v][t] - yudagI[p][w] * yddagI[v][t]) * WC5(Cquqd1R, v, w, r, s)
                        +(yudagI[p][w] * yddagR[v][t] + yudagR[p][w] * yddagI[v][t]) * WC5(Cquqd1I, v, w, r, s)
                        )
                        + 0.25 * (
                        + (yuR[v][r] * ydR[s][w] - yuI[v][r] * ydI[s][w]) * WC5(Cquqd8R, p, v, w, t)
                        -(yuI[v][r] * ydR[s][w] + yuR[v][r] * ydI[s][w]) * WC5(Cquqd8I, p, v, w, t)
                        +(yudagR[p][v] * yddagR[w][t] - yudagI[p][v] * yddagI[w][t]) * WC5(Cquqd8R, r, v, w, s)
                        +(yudagI[p][v] * yddagR[w][t] + yudagR[p][v] * yddagI[w][t]) * WC5(Cquqd8I, r, v, w, s)
                        )
                        )

                        - 2. * (
                        + 2. * (
                        + (yddagR[v][t] * ydR[s][w] - yddagI[v][t] * ydI[s][w]) * WC6R(Cqq1R, p, v, w, r)
                        -(yddagI[v][t] * ydR[s][w] + yddagR[v][t] * ydI[s][w]) * WC6I(Cqq1I, p, v, w, r)
                        )
                        -(
                        + (yddagR[p][w] * ydR[s][v] - yddagI[p][w] * ydI[s][v]) * WC7R(Cqd1R, v, r, w, t)
                        -(yddagI[p][w] * ydR[s][v] + yddagR[p][w] * ydI[s][v]) * WC7I(Cqd1I, v, r, w, t)
                        )
                        -(
                        + (yddagR[v][t] * ydR[w][r] - yddagI[v][t] * ydI[w][r]) * WC7R(Cqd1R, p, v, s, w)
                        -(yddagI[v][t] * ydR[w][r] + yddagR[v][t] * ydI[w][r]) * WC7I(Cqd1I, p, v, s, w)
                        )
                        )
                        - 4. * (
                        + (yddagR[p][v] * ydR[w][r] - yddagI[p][v] * ydI[w][r]) * WC6R(CddR, v, t, s, w)
                        -(yddagI[p][v] * ydR[w][r] + yddagR[p][v] * ydI[w][r]) * WC6I(CddI, v, t, s, w)
                        )
                        - 12. * (
                        + (yddagR[v][t] * ydR[s][w] - yddagI[v][t] * ydI[s][w]) * WC6R(Cqq3R, p, v, w, r)
                        -(yddagI[v][t] * ydR[s][w] + yddagR[v][t] * ydI[s][w]) * WC6I(Cqq3I, p, v, w, r)
                        )
                        -(
                        + (yudagR[p][v] * yuR[w][r] - yudagI[p][v] * yuI[w][r]) * WC7R(Cud8R, v, w, s, t)
                        -(yudagI[p][v] * yuR[w][r] + yudagR[p][v] * yuI[w][r]) * WC7I(Cud8I, v, w, s, t)
                        )
                        ; //RGE 2
            }
        }
        f[c] *= loop_factor;
        c ++;
    }
    //Cqd8I
    for (d = 0; d < DWC7I; d ++) {
        p = WC7I_indices[d][0];
        r = WC7I_indices[d][1];
        s = WC7I_indices[d][2];
        t = WC7I_indices[d][3];
        f[c] =
                + 2. * (ydI[s][r] * xidR[p][t] + ydR[s][r] * xidI[p][t])
                + 2. * (yddagI[p][t] * xidR[r][s] - yddagR[p][t] * xidI[r][s])
                //RGE 2
                -(12. * YdYq * g12 + 6. * (NC - (2. / NC)) * g32) * WC7I(Cqd8I, p, r, s, t)
                - 12. * g32 * WC7I(Cqd1I, p, r, s, t)
                ; //RGE 3
        for (v = 0; v < NG; v ++) {
            w = v; //Necessary since RGE 3 uses w as summed index
            f[c] +=
                    + (gqI[p][v] * WC7R(Cqd8R, v, r, s, t) + gqR[p][v] * WC7I(Cqd8I, v, r, s, t))
                    +(gdI[s][v] * WC7R(Cqd8R, p, r, v, t) + gdR[s][v] * WC7I(Cqd8I, p, r, v, t))
                    +(WC7R(Cqd8R, p, v, s, t) * gqI[v][r] + WC7I(Cqd8I, p, v, s, t) * gqR[v][r])
                    +(WC7R(Cqd8R, p, r, s, v) * gdI[v][t] + WC7I(Cqd8I, p, r, s, v) * gdR[v][t])
                    //RGE 2
                    + FOUR_THIRDS * g32 * delta[s][t]*(
                    + WC6I(Cqq1I, p, w, w, r) + WC6I(Cqq1I, w, r, p, w)
                    )
                    + 4. * g32 * delta[s][t]*(
                    + WC6I(Cqq3I, p, w, w, r) + WC6I(Cqq3I, w, r, p, w)
                    )
                    + TWO_THIRDS * g32 * WC7I(Cqu8I, p, r, w, w) * delta[s][t]
                    + TWO_THIRDS * g32 * WC7I(Cqd8I, p, r, w, w) * delta[s][t]
                    + FOUR_THIRDS * g32 * WC7I(Cqd8I, w, w, s, t) * delta[p][r]
                    + TWO_THIRDS * g32 * WC7I(Cud8I, w, w, s, t) * delta[p][r]
                    + FOUR_THIRDS * g32 * delta[p][r]*(
                    + WC6I(CddI, s, w, w, t) + WC6I(CddI, w, t, s, w)
                    )
                    ;
            for (w = 0; w < NG; w ++) {
                f[c] += - (1. / NC)*(
                        + (yddagR[p][w] * ydR[s][v] - yddagI[p][w] * ydI[s][v]) * WC7I(Cqd8I, v, r, w, t)
                        +(yddagI[p][w] * ydR[s][v] + yddagR[p][w] * ydI[s][v]) * WC7R(Cqd8R, v, r, w, t)
                        +(yddagR[v][t] * ydR[w][r] - yddagI[v][t] * ydI[w][r]) * WC7I(Cqd8I, p, v, s, w)
                        +(yddagI[v][t] * ydR[w][r] + yddagR[v][t] * ydI[w][r]) * WC7R(Cqd8R, p, v, s, w)
                        +(yuR[w][r] * ydR[s][v] - yuI[w][r] * ydI[s][v]) * WC5(Cquqd8I, v, w, p, t)
                        +(yuI[w][r] * ydR[s][v] + yuR[w][r] * ydI[s][v]) * WC5(Cquqd8R, v, w, p, t)
                        -(yudagR[p][w] * yddagR[v][t] - yudagI[p][w] * yddagI[v][t]) * WC5(Cquqd8I, v, w, r, s)
                        +(yudagI[p][w] * yddagR[v][t] + yudagR[p][w] * yddagI[v][t]) * WC5(Cquqd8R, v, w, r, s)
                        )
                        + 2. * (
                        + (
                        + (yuR[w][r] * ydR[s][v] - yuI[w][r] * ydI[s][v]) * WC5(Cquqd1I, v, w, p, t)
                        +(yuI[w][r] * ydR[s][v] + yuR[w][r] * ydI[s][v]) * WC5(Cquqd1R, v, w, p, t)
                        -(yudagR[p][w] * yddagR[v][t] - yudagI[p][w] * yddagI[v][t]) * WC5(Cquqd1I, v, w, r, s)
                        +(yudagI[p][w] * yddagR[v][t] + yudagR[p][w] * yddagI[v][t]) * WC5(Cquqd1R, v, w, r, s)
                        )
                        + 0.25 * (
                        + (yuR[v][r] * ydR[s][w] - yuI[v][r] * ydI[s][w]) * WC5(Cquqd8I, p, v, w, t)
                        +(yuI[v][r] * ydR[s][w] + yuR[v][r] * ydI[s][w]) * WC5(Cquqd8R, p, v, w, t)
                        -(yudagR[p][v] * yddagR[w][t] - yudagI[p][v] * yddagI[w][t]) * WC5(Cquqd8I, r, v, w, s)
                        +(yudagI[p][v] * yddagR[w][t] + yudagR[p][v] * yddagI[w][t]) * WC5(Cquqd8R, r, v, w, s)
                        )
                        )

                        - 2. * (
                        + 2. * (
                        + (yddagR[v][t] * ydR[s][w] - yddagI[v][t] * ydI[s][w]) * WC6I(Cqq1I, p, v, w, r)
                        +(yddagI[v][t] * ydR[s][w] + yddagR[v][t] * ydI[s][w]) * WC6R(Cqq1R, p, v, w, r)
                        )
                        -(
                        + (yddagR[p][w] * ydR[s][v] - yddagI[p][w] * ydI[s][v]) * WC7I(Cqd1I, v, r, w, t)
                        +(yddagI[p][w] * ydR[s][v] + yddagR[p][w] * ydI[s][v]) * WC7R(Cqd1R, v, r, w, t)
                        )
                        -(
                        + (yddagR[v][t] * ydR[w][r] - yddagI[v][t] * ydI[w][r]) * WC7I(Cqd1I, p, v, s, w)
                        +(yddagI[v][t] * ydR[w][r] + yddagR[v][t] * ydI[w][r]) * WC7R(Cqd1R, p, v, s, w)
                        )
                        )

                        - 4. * (
                        + (yddagR[p][v] * ydR[w][r] - yddagI[p][v] * ydI[w][r]) * WC6I(CddI, v, t, s, w)
                        +(yddagI[p][v] * ydR[w][r] + yddagR[p][v] * ydI[w][r]) * WC6R(CddR, v, t, s, w)
                        )
                        - 12. * (
                        + (yddagR[v][t] * ydR[s][w] - yddagI[v][t] * ydI[s][w]) * WC6I(Cqq3I, p, v, w, r)
                        +(yddagI[v][t] * ydR[s][w] + yddagR[v][t] * ydI[s][w]) * WC6R(Cqq3R, p, v, w, r)
                        )
                        -(
                        + (yudagR[p][v] * yuR[w][r] - yudagI[p][v] * yuI[w][r]) * WC7I(Cud8I, v, w, s, t)
                        +(yudagI[p][v] * yuR[w][r] + yudagR[p][v] * yuI[w][r]) * WC7R(Cud8R, v, w, s, t)
                        )
                        ; //RGE 2
            }
        }
        f[c] *= loop_factor;
        c ++;
    }


    //----------------------RGE SMEFT class 8_LRRL----------------------
    for (p = 0; p < NG; p ++) {
        for (r = 0; r < NG; r ++) {
            for (s = 0; s < NG; s ++) {
                for (t = 0; t < NG; t ++) {
                    //CledqR
                    f[c] = - 2. * (
                            + (ydR[s][t] * xieR[p][r] - ydI[s][t] * xieI[p][r])
                            +(yedagR[p][r] * xidR[t][s] + yedagI[p][r] * xidI[t][s])
                            ) //RGE 2
                            -(
                            6. * (Yd * (Yq - Ye) + Ye * (Ye + Yq)) * g12
                            + 3. * (NC - (1. / NC)) * g32
                            ) * WC5(CledqR, p, r, s, t); //RGE 3
                    //CledqI
                    f[c + NG * NG * NG * NG] = - 2. * (
                            + (ydI[s][t] * xieR[p][r] + ydR[s][t] * xieI[p][r])
                            +(yedagI[p][r] * xidR[t][s] - yedagR[p][r] * xidI[t][s])
                            ) //RGE 2
                            -(
                            6. * (Yd * (Yq - Ye) + Ye * (Ye + Yq)) * g12
                            + 3. * (NC - (1. / NC)) * g32
                            ) * WC5(CledqI, p, r, s, t); //RGE 3

                    //Entries with 1 summed index
                    for (v = 0; v < NG; v ++) {
                        //CledqR
                        f[c] +=
                                + (glR[p][v] * WC5(CledqR, v, r, s, t) - glI[p][v] * WC5(CledqI, v, r, s, t))
                                +(gdR[s][v] * WC5(CledqR, p, r, v, t) - gdI[s][v] * WC5(CledqI, p, r, v, t))
                                +(WC5(CledqR, p, v, s, t) * geR[v][r] - WC5(CledqI, p, v, s, t) * geI[v][r])
                                +(WC5(CledqR, p, r, s, v) * gqR[v][t] - WC5(CledqI, p, r, s, v) * gqI[v][t])
                                ; //RGE 2
                        //CledqI
                        f[c + NG * NG * NG * NG] +=
                                + (glI[p][v] * WC5(CledqR, v, r, s, t) + glR[p][v] * WC5(CledqI, v, r, s, t))
                                +(gdI[s][v] * WC5(CledqR, p, r, v, t) + gdR[s][v] * WC5(CledqI, p, r, v, t))
                                +(WC5(CledqR, p, v, s, t) * geI[v][r] + WC5(CledqI, p, v, s, t) * geR[v][r])
                                +(WC5(CledqR, p, r, s, v) * gqI[v][t] + WC5(CledqI, p, r, s, v) * gqR[v][t])
                                ; //RGE 2
                        for (w = 0; w < NG; w ++) {
                            //CledqR
                            f[c] += 2. * (
                                    + (yedagR[p][v] * ydR[w][t] - yedagI[p][v] * ydI[w][t]) * WC7R(CedR, v, r, s, w)
                                    -(yedagI[p][v] * ydR[w][t] + yedagR[p][v] * ydI[w][t]) * WC7I(CedI, v, r, s, w))
                                    - 2. * (
                                    + (yedagR[v][r] * ydR[w][t] - yedagI[v][r] * ydI[w][t]) * WC7R(CldR, p, v, s, w)
                                    -(yedagI[v][r] * ydR[w][t] + yedagR[v][r] * ydI[w][t]) * WC7I(CldI, p, v, s, w))
                                    + 2. * (
                                    + (yedagR[v][r] * ydR[s][w] - yedagI[v][r] * ydI[s][w]) * WC7R(Clq1R, p, v, w, t)
                                    -(yedagI[v][r] * ydR[s][w] + yedagR[v][r] * ydI[s][w]) * WC7I(Clq1I, p, v, w, t))
                                    + 6. * (
                                    + (yedagR[v][r] * ydR[s][w] - yedagI[v][r] * ydI[s][w]) * WC7R(Clq3R, p, v, w, t)
                                    -(yedagI[v][r] * ydR[s][w] + yedagR[v][r] * ydI[s][w]) * WC7I(Clq3I, p, v, w, t))
                                    - 2. * (
                                    + (yedagR[p][w] * ydR[s][v] - yedagI[p][w] * ydI[s][v]) * WC7R(CqeR, v, t, w, r)
                                    -(yedagI[p][w] * ydR[s][v] + yedagR[p][w] * ydI[s][v]) * WC7I(CqeI, v, t, w, r))
                                    + 2. * (
                                    + (ydR[s][v] * yuR[w][t] - ydI[s][v] * yuI[w][t]) * WC5(Clequ1R, p, r, v, w)
                                    -(ydI[s][v] * yuR[w][t] + ydR[s][v] * yuI[w][t]) * WC5(Clequ1I, p, r, v, w))
                                    ; //RGE 2
                            //CledqI
                            f[c + NG * NG * NG * NG] += 2. * (
                                    + (yedagR[p][v] * ydR[w][t] - yedagI[p][v] * ydI[w][t]) * WC7I(CedI, v, r, s, w)
                                    +(yedagI[p][v] * ydR[w][t] + yedagR[p][v] * ydI[w][t]) * WC7R(CedR, v, r, s, w))
                                    - 2. * (
                                    + (yedagR[v][r] * ydR[w][t] - yedagI[v][r] * ydI[w][t]) * WC7I(CldI, p, v, s, w)
                                    +(yedagI[v][r] * ydR[w][t] + yedagR[v][r] * ydI[w][t]) * WC7R(CldR, p, v, s, w))
                                    + 2. * (
                                    + (yedagR[v][r] * ydR[s][w] - yedagI[v][r] * ydI[s][w]) * WC7I(Clq1I, p, v, w, t)
                                    +(yedagI[v][r] * ydR[s][w] + yedagR[v][r] * ydI[s][w]) * WC7R(Clq1R, p, v, w, t))
                                    + 6. * (
                                    + (yedagR[v][r] * ydR[s][w] - yedagI[v][r] * ydI[s][w]) * WC7I(Clq3I, p, v, w, t)
                                    +(yedagI[v][r] * ydR[s][w] + yedagR[v][r] * ydI[s][w]) * WC7R(Clq3R, p, v, w, t))
                                    - 2. * (
                                    + (yedagR[p][w] * ydR[s][v] - yedagI[p][w] * ydI[s][v]) * WC7I(CqeI, v, t, w, r)
                                    +(yedagI[p][w] * ydR[s][v] + yedagR[p][w] * ydI[s][v]) * WC7R(CqeR, v, t, w, r))
                                    + 2. * (
                                    + (ydR[s][v] * yuR[w][t] - ydI[s][v] * yuI[w][t]) * WC5(Clequ1I, p, r, v, w)
                                    +(ydI[s][v] * yuR[w][t] + ydR[s][v] * yuI[w][t]) * WC5(Clequ1R, p, r, v, w))
                                    ; //RGE 2
                        }
                    }

                    f[c] *= loop_factor;
                    f[c + NG * NG * NG * NG] *= loop_factor;
                    c ++;
                }
            }
        }
    }
    c += NG * NG * NG * NG;


    //----------------------RGE SMEFT class 8_LRLR----------------------
    for (p = 0; p < NG; p ++) {
        for (r = 0; r < NG; r ++) {
            for (s = 0; s < NG; s ++) {
                for (t = 0; t < NG; t ++) {
                    count = 0;
                    //Cquqd1R
                    f[c + count * NG * NG * NG * NG] =
                            - 2. * (yudagR[p][r] * xidR[s][t] - yudagI[p][r] * xidI[s][t])
                            - 2. * (yddagR[s][t] * xiuR[p][r] - yddagI[s][t] * xiuI[p][r])
                            //RGE 2
                            + 4. * g1 * (Yq + Yu)*
                            (WC1(CdBR, s, t) * yudagR[p][r] - WC1(CdBI, s, t) * yudagI[p][r])
                            - 6. * g2 *
                            (WC1(CdWR, s, t) * yudagR[p][r] - WC1(CdWI, s, t) * yudagI[p][r])
                            -(8. / NC) * g1 * (Yq + Yu)*
                            (WC1(CdBR, p, t) * yudagR[s][r] - WC1(CdBI, p, t) * yudagI[s][r])
                            +(12. / NC) * g2 *
                            (WC1(CdWR, p, t) * yudagR[s][r] - WC1(CdWI, p, t) * yudagI[s][r])
                            - 8. * ((NC2 - 1.) / NC2) * g3 *
                            (WC1(CdGR, p, t) * yudagR[s][r] - WC1(CdGI, p, t) * yudagI[s][r])

                            + 4. * g1 * (Yq + Yd)*
                            (WC1(CuBR, p, r) * yddagR[s][t] - WC1(CuBI, p, r) * yddagI[s][t])
                            - 6. * g2 *
                            (WC1(CuWR, p, r) * yddagR[s][t] - WC1(CuWI, p, r) * yddagI[s][t])
                            -(8. / NC) * g1 * (Yq + Yd)*
                            (WC1(CuBR, s, r) * yddagR[p][t] - WC1(CuBI, s, r) * yddagI[p][t])
                            +(12. / NC) * g2 *
                            (WC1(CuWR, s, r) * yddagR[p][t] - WC1(CuWI, s, r) * yddagI[p][t])
                            - 8. * ((NC2 - 1.) / NC2) * g3 *
                            (WC1(CuGR, s, r) * yddagR[p][t] - WC1(CuGI, s, r) * yddagI[p][t])

                            - 0.5 * (
                            (3. * Yd2 + 2. * YuYd + 3. * Yu2) * g12 + 3. * g22
                            + 12. * (NC - (1. / NC)) * g32) * WC5(Cquqd1R, p, r, s, t)

                            -(1. / NC)*(
                            + (3. * Yd2 + 10 * YuYd + 3. * Yu2) * g12 - 3. * g22
                            + 8. * (NC - (1. / NC)) * g32
                            ) * WC5(Cquqd1R, s, r, p, t)

                            - 0.5 * (1. - (1. / NC2))*(
                            (3. * Yd2 + 10. * YuYd + 3. * Yu2) * g12
                            - 3. * g22 + 4. * (NC - (2. / NC)) * g32) * WC5(Cquqd8R, s, r, p, t)
                            + 2. * (1. - (1. / NC2)) * g32 * WC5(Cquqd8R, p, r, s, t)
                            //RGE 3
                            ;
                    count ++;
                    //Cquqd1I
                    f[c + count * NG * NG * NG * NG] =
                            - 2. * (yudagI[p][r] * xidR[s][t] + yudagR[p][r] * xidI[s][t])
                            - 2. * (yddagI[s][t] * xiuR[p][r] + yddagR[s][t] * xiuI[p][r])
                            //RGE 2
                            + 4. * g1 * (Yq + Yu)*
                            (WC1(CdBR, s, t) * yudagI[p][r] + WC1(CdBI, s, t) * yudagR[p][r])
                            - 6. * g2 *
                            (WC1(CdWR, s, t) * yudagI[p][r] + WC1(CdWI, s, t) * yudagR[p][r])
                            -(8. / NC) * g1 * (Yq + Yu)*
                            (WC1(CdBR, p, t) * yudagI[s][r] + WC1(CdBI, p, t) * yudagR[s][r])
                            +(12. / NC) * g2 *
                            (WC1(CdWR, p, t) * yudagI[s][r] + WC1(CdWI, p, t) * yudagR[s][r])
                            - 8. * ((NC2 - 1.) / NC2) * g3 *
                            (WC1(CdGR, p, t) * yudagI[s][r] + WC1(CdGI, p, t) * yudagR[s][r])

                            + 4. * g1 * (Yq + Yd)*
                            (WC1(CuBR, p, r) * yddagI[s][t] + WC1(CuBI, p, r) * yddagR[s][t])
                            - 6. * g2 *
                            (WC1(CuWR, p, r) * yddagI[s][t] + WC1(CuWI, p, r) * yddagR[s][t])
                            -(8. / NC) * g1 * (Yq + Yd)*
                            (WC1(CuBR, s, r) * yddagI[p][t] + WC1(CuBI, s, r) * yddagR[p][t])
                            +(12. / NC) * g2 *
                            (WC1(CuWR, s, r) * yddagI[p][t] + WC1(CuWI, s, r) * yddagR[p][t])
                            - 8. * ((NC2 - 1.) / NC2) * g3 *
                            (WC1(CuGR, s, r) * yddagI[p][t] + WC1(CuGI, s, r) * yddagR[p][t])

                            - 0.5 * ((3. * Yd2 + 2. * YuYd + 3. * Yu2) * g12 + 3. * g22
                            + 12. * (NC - (1. / NC)) * g32) * WC5(Cquqd1I, p, r, s, t)
                            -(1. / NC)*(
                            + (3. * Yd2 + 10 * YuYd + 3. * Yu2) * g12 - 3. * g22
                            + 8. * (NC - (1. / NC)) * g32
                            ) * WC5(Cquqd1I, s, r, p, t)
                            - 0.5 * (1 - (1. / NC2))*((3. * Yd2 + 10. * YuYd + 3. * Yu2) * g12
                            - 3. * g22 + 4. * (NC - (2. / NC)) * g32) * WC5(Cquqd8I, s, r, p, t)
                            + 2. * (1. - (1. / NC2)) * g32 * WC5(Cquqd8I, p, r, s, t)
                            //RGE 3
                            ;
                    count ++;


                    //Cquqd8R
                    f[c + count * NG * NG * NG * NG] = + 8. * g3 *
                            (WC1(CdGR, s, t) * yudagR[p][r] - WC1(CdGI, s, t) * yudagI[p][r])
                            - 16. * g1 * (Yq + Yu)*
                            (WC1(CdBR, p, t) * yudagR[s][r] - WC1(CdBI, p, t) * yudagI[s][r])
                            + 24. * g2 *
                            (WC1(CdWR, p, t) * yudagR[s][r] - WC1(CdWI, p, t) * yudagI[s][r])
                            +(16. / NC) * g3 *
                            (WC1(CdGR, p, t) * yudagR[s][r] - WC1(CdGI, p, t) * yudagI[s][r])

                            + 8. * g3 *
                            (WC1(CuGR, p, r) * yddagR[s][t] - WC1(CuGI, p, r) * yddagI[s][t])
                            - 16. * g1 * (Yq + Yd)*
                            (WC1(CuBR, s, r) * yddagR[p][t] - WC1(CuBI, s, r) * yddagI[p][t])
                            + 24. * g2 *
                            (WC1(CuWR, s, r) * yddagR[p][t] - WC1(CuWI, s, r) * yddagI[p][t])
                            +(16. / NC) * g3 *
                            (WC1(CuGR, s, r) * yddagR[p][t] - WC1(CuGI, s, r) * yddagI[p][t])

                            + 8. * g32 * WC5(Cquqd1R, p, r, s, t)
                            +(- 2. * (3. * Yd2 + 10. * YuYd + 3. * Yu2) * g12 + 6. * g22
                            + (16. / NC) * g32) * WC5(Cquqd1R, s, r, p, t)

                            +((- 1.5 * Yd2 - YuYd - 1.5 * Yu2) * g12 - 1.5 * g22
                            + 2. * (NC - (1. / NC)) * g32) * WC5(Cquqd8R, p, r, s, t)

                            +(1 / NC)*(
                            (3. * Yd2 + 10. * YuYd + 3. * Yu2) * g12 - 3. * g22
                            + 4. * (- NC - 2. / NC) * g32) * WC5(Cquqd8R, s, r, p, t)
                            ; //RGE 3
                    count ++;
                    //Cquqd8I
                    f[c + count * NG * NG * NG * NG] = + 8. * g3 *
                            (WC1(CdGR, s, t) * yudagI[p][r] + WC1(CdGI, s, t) * yudagR[p][r])
                            - 16. * g1 * (Yq + Yu)*
                            (WC1(CdBR, p, t) * yudagI[s][r] + WC1(CdBI, p, t) * yudagR[s][r])
                            + 24. * g2 *
                            (WC1(CdWR, p, t) * yudagI[s][r] + WC1(CdWI, p, t) * yudagR[s][r])
                            +(16. / NC) * g3 *
                            (WC1(CdGR, p, t) * yudagI[s][r] + WC1(CdGI, p, t) * yudagR[s][r])

                            + 8. * g3 *
                            (WC1(CuGR, p, r) * yddagI[s][t] + WC1(CuGI, p, r) * yddagR[s][t])
                            - 16. * g1 * (Yq + Yd)*
                            (WC1(CuBR, s, r) * yddagI[p][t] + WC1(CuBI, s, r) * yddagR[p][t])
                            + 24. * g2 *
                            (WC1(CuWR, s, r) * yddagI[p][t] + WC1(CuWI, s, r) * yddagR[p][t])
                            +(16. / NC) * g3 *
                            (WC1(CuGR, s, r) * yddagI[p][t] + WC1(CuGI, s, r) * yddagR[p][t])

                            + 8. * g32 * WC5(Cquqd1I, p, r, s, t)
                            +(- 2. * (3. * Yd2 + 10. * YuYd + 3. * Yu2) * g12 + 6. * g22
                            + (16. / NC) * g32) * WC5(Cquqd1I, s, r, p, t)
                            +((- 1.5 * Yd2 - YuYd - 1.5 * Yu2) * g12 - 1.5 * g22
                            + 2. * (NC - (1. / NC)) * g32) * WC5(Cquqd8I, p, r, s, t)

                            +(1. / NC)*((3. * Yd2 + 10. * YuYd + 3. * Yu2) * g12 - 3. * g22
                            + 4. * (- NC - 2. / NC) * g32) * WC5(Cquqd8I, s, r, p, t)
                            ; //RGE 3
                    count ++;
                    //Clequ1R
                    f[c + count * NG * NG * NG * NG] =
                            + 2. * (yudagR[s][t] * xieR[p][r] - yudagI[s][t] * xieI[p][r])
                            + 2. * (yedagR[p][r] * xiuR[s][t] - yedagI[p][r] * xiuI[s][t])
                            //RGE 2
                            -(6. * (Ye2 + Ye * (Yu - Yq) + YuYq) * g12
                            + 3. * (NC - (1. / NC)) * g32) * WC5(Clequ1R, p, r, s, t)
                            -(24. * (Yu + Yq)*(2. * Ye - Yq + Yu) * g12 - 18. * g22) * WC5(Clequ3R, p, r, s, t)
                            ; //RGE 3
                    count ++;
                    //Clequ1I
                    f[c + count * NG * NG * NG * NG] =
                            + 2. * (yudagI[s][t] * xieR[p][r] + yudagR[s][t] * xieI[p][r])
                            + 2. * (yedagI[p][r] * xiuR[s][t] + yedagR[p][r] * xiuI[s][t])
                            //RGE 2
                            -(6. * (Ye2 + Ye * (Yu - Yq) + YuYq) * g12
                            + 3. * (NC - (1 / NC)) * g32) * WC5(Clequ1I, p, r, s, t)
                            -(24. * (Yu + Yq)*(2. * Ye - Yq + Yu) * g12 - 18. * g22) * WC5(Clequ3I, p, r, s, t)
                            ; //RGE 3
                    count ++;
                    //Clequ3R
                    f[c + count * NG * NG * NG * NG] = g1 * (Yq + Yu)*
                            (WC1(CeBR, p, r) * yudagR[s][t] - WC1(CeBI, p, r) * yudagI[s][t])
                            - 1.5 * g2 *
                            (WC1(CuWR, s, t) * yedagR[p][r] - WC1(CuWI, s, t) * yedagI[p][r])
                            + g1 * (Yl + Ye)*
                            (WC1(CuBR, s, t) * yedagR[p][r] - WC1(CuBI, s, t) * yedagI[p][r])
                            - 1.5 * g2 *
                            (WC1(CeWR, p, r) * yudagR[s][t] - WC1(CeWI, p, r) * yudagI[s][t])

                            +((2. * (Ye2 - YeYq + YuYe - 2. * Yq2 + 5. * YuYq - 2. * Yu2) * g12 - 3. * g22)
                            +(NC - (1 / NC)) * g32) * WC5(Clequ3R, p, r, s, t)

                            + 0.125 * (- 4. * (Yq + Yu)*(2. * Ye - Yq + Yu) * g12 + 3. * g22) * WC5(Clequ1R, p, r, s, t)

                            ; //RGE 3
                    count ++;
                    //Clequ3I
                    f[c + count * NG * NG * NG * NG] = g1 * (Yq + Yu)*
                            (WC1(CeBR, p, r) * yudagI[s][t] + WC1(CeBI, p, r) * yudagR[s][t])
                            - 1.5 * g2 *
                            (WC1(CuWR, s, t) * yedagI[p][r] + WC1(CuWI, s, t) * yedagR[p][r])
                            + g1 * (Yl + Ye)*
                            (WC1(CuBR, s, t) * yedagI[p][r] + WC1(CuBI, s, t) * yedagR[p][r])
                            - 1.5 * g2 *
                            (WC1(CeWR, p, r) * yudagI[s][t] + WC1(CeWI, p, r) * yudagR[s][t])

                            +((2. * (Ye2 - YeYq + YuYe - 2. * Yq2 + 5. * YuYq - 2. * Yu2) * g12 - 3. * g22)
                            +(NC - (1 / NC)) * g32) * WC5(Clequ3I, p, r, s, t)

                            + 0.125 * (- 4. * (Yq + Yu)*(2. * Ye - Yq + Yu) * g12 + 3. * g22) * WC5(Clequ1I, p, r, s, t)

                            ; //RGE 3
                    count ++;
                    //Entries with 1 summed index
                    for (v = 0; v < NG; v ++) {
                        count = 0;
                        //Cquqd1R
                        f[c + count * NG * NG * NG * NG] +=
                                + (gqR[p][v] * WC5(Cquqd1R, v, r, s, t) - gqI[p][v] * WC5(Cquqd1I, v, r, s, t))
                                +(gqR[s][v] * WC5(Cquqd1R, p, r, v, t) - gqI[s][v] * WC5(Cquqd1I, p, r, v, t))
                                +(WC5(Cquqd1R, p, v, s, t) * guR[v][r] - WC5(Cquqd1I, p, v, s, t) * guI[v][r])
                                +(WC5(Cquqd1R, p, r, s, v) * gdR[v][t] - WC5(Cquqd1I, p, r, s, v) * gdI[v][t])
                                ; //RGE 2
                        count ++;
                        //Cquqd1I
                        f[c + count * NG * NG * NG * NG] +=
                                + (gqI[p][v] * WC5(Cquqd1R, v, r, s, t) + gqR[p][v] * WC5(Cquqd1I, v, r, s, t))
                                +(gqI[s][v] * WC5(Cquqd1R, p, r, v, t) + gqR[s][v] * WC5(Cquqd1I, p, r, v, t))
                                +(WC5(Cquqd1R, p, v, s, t) * guI[v][r] + WC5(Cquqd1I, p, v, s, t) * guR[v][r])
                                +(WC5(Cquqd1R, p, r, s, v) * gdI[v][t] + WC5(Cquqd1I, p, r, s, v) * gdR[v][t])
                                ; //RGE 2
                        count ++;
                        //Cquqd8R
                        f[c + count * NG * NG * NG * NG] +=
                                + (gqR[p][v] * WC5(Cquqd8R, v, r, s, t) - gqI[p][v] * WC5(Cquqd8I, v, r, s, t))
                                +(gqR[s][v] * WC5(Cquqd8R, p, r, v, t) - gqI[s][v] * WC5(Cquqd8I, p, r, v, t))
                                +(WC5(Cquqd8R, p, v, s, t) * guR[v][r] - WC5(Cquqd8I, p, v, s, t) * guI[v][r])
                                +(WC5(Cquqd8R, p, r, s, v) * gdR[v][t] - WC5(Cquqd8I, p, r, s, v) * gdI[v][t])
                                ; //RGE 2;
                        count ++;
                        //Cquqd8I
                        f[c + count * NG * NG * NG * NG] +=
                                + (gqI[p][v] * WC5(Cquqd8R, v, r, s, t) + gqR[p][v] * WC5(Cquqd8I, v, r, s, t))
                                +(gqI[s][v] * WC5(Cquqd8R, p, r, v, t) + gqR[s][v] * WC5(Cquqd8I, p, r, v, t))
                                +(WC5(Cquqd8R, p, v, s, t) * guI[v][r] + WC5(Cquqd8I, p, v, s, t) * guR[v][r])
                                +(WC5(Cquqd8R, p, r, s, v) * gdI[v][t] + WC5(Cquqd8I, p, r, s, v) * gdR[v][t])
                                ; //RGE 2
                        count ++;
                        //Clequ1R
                        f[c + count * NG * NG * NG * NG] +=
                                + (glR[p][v] * WC5(Clequ1R, v, r, s, t) - glI[p][v] * WC5(Clequ1I, v, r, s, t))
                                +(gqR[s][v] * WC5(Clequ1R, p, r, v, t) - gqI[s][v] * WC5(Clequ1I, p, r, v, t))
                                +(WC5(Clequ1R, p, v, s, t) * geR[v][r] - WC5(Clequ1I, p, v, s, t) * geI[v][r])
                                +(WC5(Clequ1R, p, r, s, v) * guR[v][t] - WC5(Clequ1I, p, r, s, v) * guI[v][t])
                                ; //RGE 2
                        count ++;
                        //Clequ1I
                        f[c + count * NG * NG * NG * NG] +=
                                + (glI[p][v] * WC5(Clequ1R, v, r, s, t) + glR[p][v] * WC5(Clequ1I, v, r, s, t))
                                +(gqI[s][v] * WC5(Clequ1R, p, r, v, t) + gqR[s][v] * WC5(Clequ1I, p, r, v, t))
                                +(WC5(Clequ1R, p, v, s, t) * geI[v][r] + WC5(Clequ1I, p, v, s, t) * geR[v][r])
                                +(WC5(Clequ1R, p, r, s, v) * guI[v][t] + WC5(Clequ1I, p, r, s, v) * guR[v][t])
                                ; //RGE 2
                        count ++;
                        //Clequ3R
                        f[c + count * NG * NG * NG * NG] +=
                                + (glR[p][v] * WC5(Clequ3R, v, r, s, t) - glI[p][v] * WC5(Clequ3I, v, r, s, t))
                                +(gqR[s][v] * WC5(Clequ3R, p, r, v, t) - gqI[s][v] * WC5(Clequ3I, p, r, v, t))
                                +(WC5(Clequ3R, p, v, s, t) * geR[v][r] - WC5(Clequ3I, p, v, s, t) * geI[v][r])
                                +(WC5(Clequ3R, p, r, s, v) * guR[v][t] - WC5(Clequ3I, p, r, s, v) * guI[v][t])
                                ; //RGE 2
                        count ++;
                        //Clequ1I
                        f[c + count * NG * NG * NG * NG] +=
                                + (glI[p][v] * WC5(Clequ3R, v, r, s, t) + glR[p][v] * WC5(Clequ3I, v, r, s, t))
                                +(gqI[s][v] * WC5(Clequ3R, p, r, v, t) + gqR[s][v] * WC5(Clequ3I, p, r, v, t))
                                +(WC5(Clequ3R, p, v, s, t) * geI[v][r] + WC5(Clequ3I, p, v, s, t) * geR[v][r])
                                +(WC5(Clequ3R, p, r, s, v) * guI[v][t] + WC5(Clequ3I, p, r, s, v) * guR[v][t])
                                ; //RGE 2
                        count ++;
                        //Entries with 2 summed indices
                        for (w = 0; w < NG; w ++) {
                            count = 0;
                            //Cquqd1R
                            f[c + count * NG * NG * NG * NG] += - (2. / NC2)*(
                                    + (yudagR[v][r] * yddagR[p][w] - yudagI[v][r] * yddagI[p][w]) * WC7R(Cqd8R, s, v, w, t)
                                    -(yudagI[v][r] * yddagR[p][w] + yudagR[v][r] * yddagI[p][w]) * WC7I(Cqd8I, s, v, w, t)
                                    +(yddagR[v][t] * yudagR[s][w] - yddagI[v][t] * yudagI[s][w]) * WC7R(Cqu8R, p, v, w, r)
                                    - (yddagI[v][t] * yudagR[s][w] + yddagR[v][t] * yudagI[s][w]) * WC7I(Cqu8I, p, v, w, r)
                                    +(yddagR[p][w] * yudagR[s][v] - yddagI[p][w] * yudagI[s][v]) * WC7R(Cud8R, v, r, w, t)
                                    - (yddagI[p][w] * yudagR[s][v] + yddagR[p][w] * yudagI[s][v]) * WC7I(Cud8I, v, r, w, t)
                                    )

                                    +(4. / NC)*(
                                    + (yddagR[w][t] * yudagR[v][r] - yddagI[w][t] * yudagI[v][r]) * WC6R(Cqq1R, s, v, p, w)
                                    -(yddagI[w][t] * yudagR[v][r] + yddagR[w][t] * yudagI[v][r]) * WC6I(Cqq1I, s, v, p, w)
                                    +(yddagR[v][t] * yudagR[w][r] - yddagI[v][t] * yudagI[w][r]) * WC6R(Cqq1R, p, v, s, w)
                                    -(yddagI[v][t] * yudagR[w][r] + yddagR[v][t] * yudagI[w][r]) * WC6I(Cqq1I, p, v, s, w)
                                    - 3. * ((yddagR[w][t] * yudagR[v][r] - yddagI[w][t] * yudagI[v][r]) * WC6R(Cqq3R, s, v, p, w)
                                    -(yddagI[w][t] * yudagR[v][r] + yddagR[w][t] * yudagI[v][r]) * WC6I(Cqq3I, s, v, p, w))
                                    - 3. * ((yddagR[v][t] * yudagR[w][r] - yddagI[v][t] * yudagI[w][r]) * WC6R(Cqq3R, p, v, s, w)
                                    -(yddagI[v][t] * yudagR[w][r] + yddagR[v][t] * yudagI[w][r]) * WC6I(Cqq3I, p, v, s, w))
                                    )

                                    +(4. / NC)*(
                                    + (yddagR[p][w] * yudagR[v][r] - yddagI[p][w] * yudagI[v][r]) * WC7R(Cqd1R, s, v, w, t)
                                    -(yddagI[p][w] * yudagR[v][r] + yddagR[p][w] * yudagI[v][r]) * WC7I(Cqd1I, s, v, w, t)
                                    +(yddagR[v][t] * yudagR[s][w] - yddagI[v][t] * yudagI[s][w]) * WC7R(Cqu1R, p, v, w, r)
                                    -(yddagI[v][t] * yudagR[s][w] + yddagR[v][t] * yudagI[s][w]) * WC7I(Cqu1I, p, v, w, r)
                                    +(yddagR[p][w] * yudagR[s][v] - yddagI[p][w] * yudagI[s][v]) * WC7R(Cud1R, v, r, w, t)
                                    -(yddagI[p][w] * yudagR[s][v] + yddagR[p][w] * yudagI[s][v]) * WC7I(Cud1I, v, r, w, t)
                                    )

                                    - 4. * (
                                    + (yddagR[w][t] * yudagR[v][r] - yddagI[w][t] * yudagI[v][r]) * WC6R(Cqq1R, p, v, s, w)
                                    -(yddagI[w][t] * yudagR[v][r] + yddagR[w][t] * yudagI[v][r]) * WC6I(Cqq1I, p, v, s, w)
                                    +(yddagR[v][t] * yudagR[w][r] - yddagI[v][t] * yudagI[w][r]) * WC6R(Cqq1R, s, v, p, w)
                                    -(yddagI[v][t] * yudagR[w][r] + yddagR[v][t] * yudagI[w][r]) * WC6I(Cqq1I, s, v, p, w)
                                    )

                                    + 12. * (
                                    + (yddagR[w][t] * yudagR[v][r] - yddagI[w][t] * yudagI[v][r]) * WC6R(Cqq3R, p, v, s, w)
                                    -(yddagI[w][t] * yudagR[v][r] + yddagR[w][t] * yudagI[v][r]) * WC6I(Cqq3I, p, v, s, w)
                                    +(yddagR[v][t] * yudagR[w][r] - yddagI[v][t] * yudagI[w][r]) * WC6R(Cqq3R, s, v, p, w)
                                    -(yddagI[v][t] * yudagR[w][r] + yddagR[v][t] * yudagI[w][r]) * WC6I(Cqq3I, s, v, p, w)
                                    )

                                    + 2. * (
                                    + (yddagR[p][w] * yudagR[v][r] - yddagI[p][w] * yudagI[v][r]) * WC7R(Cqd8R, s, v, w, t)
                                    - (yddagI[p][w] * yudagR[v][r] + yddagR[p][w] * yudagI[v][r]) * WC7I(Cqd8I, s, v, w, t)
                                    +(yddagR[v][t] * yudagR[s][w] - yddagI[v][t] * yudagI[s][w]) * WC7R(Cqu8R, p, v, w, r)
                                    - (yddagI[v][t] * yudagR[s][w] + yddagR[v][t] * yudagI[s][w]) * WC7I(Cqu8I, p, v, w, r)
                                    +(yddagR[p][w] * yudagR[s][v] - yddagI[p][w] * yudagI[s][v]) * WC7R(Cud8R, v, r, w, t)
                                    - (yddagI[p][w] * yudagR[s][v] + yddagR[p][w] * yudagI[s][v]) * WC7I(Cud8I, v, r, w, t)
                                    )
                                    - 4. * (
                                    + (yddagR[s][w] * yudagR[p][v] - yddagI[s][w] * yudagI[p][v]) * WC7R(Cud1R, v, r, w, t)
                                    -(yddagI[s][w] * yudagR[p][v] + yddagR[s][w] * yudagI[p][v]) * WC7I(Cud1I, v, r, w, t)
                                    )
                                    ; //RGE 2
                            count ++;
                            //Cquqd1I
                            f[c + count * NG * NG * NG * NG] += - (2. / NC2)*(
                                    + (yudagR[v][r] * yddagR[p][w] - yudagI[v][r] * yddagI[p][w]) * WC7I(Cqd8I, s, v, w, t)
                                    +(yudagI[v][r] * yddagR[p][w] + yudagR[v][r] * yddagI[p][w]) * WC7R(Cqd8R, s, v, w, t)
                                    +(yddagR[v][t] * yudagR[s][w] - yddagI[v][t] * yudagI[s][w]) * WC7I(Cqu8I, p, v, w, r)
                                    + (yddagI[v][t] * yudagR[s][w] + yddagR[v][t] * yudagI[s][w]) * WC7R(Cqu8R, p, v, w, r)
                                    +(yddagR[p][w] * yudagR[s][v] - yddagI[p][w] * yudagI[s][v]) * WC7I(Cud8I, v, r, w, t)
                                    + (yddagI[p][w] * yudagR[s][v] + yddagR[p][w] * yudagI[s][v]) * WC7R(Cud8R, v, r, w, t)
                                    )

                                    +(4. / NC)*(
                                    + (yddagR[w][t] * yudagR[v][r] - yddagI[w][t] * yudagI[v][r]) * WC6I(Cqq1I, s, v, p, w)
                                    +(yddagI[w][t] * yudagR[v][r] + yddagR[w][t] * yudagI[v][r]) * WC6R(Cqq1R, s, v, p, w)
                                    +(yddagR[v][t] * yudagR[w][r] - yddagI[v][t] * yudagI[w][r]) * WC6I(Cqq1I, p, v, s, w)
                                    +(yddagI[v][t] * yudagR[w][r] + yddagR[v][t] * yudagI[w][r]) * WC6R(Cqq1R, p, v, s, w)
                                    - 3. * ((yddagR[w][t] * yudagR[v][r] - yddagI[w][t] * yudagI[v][r]) * WC6I(Cqq3I, s, v, p, w)
                                    +(yddagI[w][t] * yudagR[v][r] + yddagR[w][t] * yudagI[v][r]) * WC6R(Cqq3R, s, v, p, w))
                                    - 3. * ((yddagR[v][t] * yudagR[w][r] - yddagI[v][t] * yudagI[w][r]) * WC6I(Cqq3I, p, v, s, w)
                                    +(yddagI[v][t] * yudagR[w][r] + yddagR[v][t] * yudagI[w][r]) * WC6R(Cqq3R, p, v, s, w))
                                    )

                                    +(4. / NC)*(
                                    + (yddagR[p][w] * yudagR[v][r] - yddagI[p][w] * yudagI[v][r]) * WC7I(Cqd1I, s, v, w, t)
                                    +(yddagI[p][w] * yudagR[v][r] + yddagR[p][w] * yudagI[v][r]) * WC7R(Cqd1R, s, v, w, t)
                                    +(yddagR[v][t] * yudagR[s][w] - yddagI[v][t] * yudagI[s][w]) * WC7I(Cqu1I, p, v, w, r)
                                    +(yddagI[v][t] * yudagR[s][w] + yddagR[v][t] * yudagI[s][w]) * WC7R(Cqu1R, p, v, w, r)
                                    +(yddagR[p][w] * yudagR[s][v] - yddagI[p][w] * yudagI[s][v]) * WC7I(Cud1I, v, r, w, t)
                                    +(yddagI[p][w] * yudagR[s][v] + yddagR[p][w] * yudagI[s][v]) * WC7R(Cud1R, v, r, w, t)
                                    )

                                    - 4. * (
                                    + (yddagR[w][t] * yudagR[v][r] - yddagI[w][t] * yudagI[v][r]) * WC6I(Cqq1I, p, v, s, w)
                                    +(yddagI[w][t] * yudagR[v][r] + yddagR[w][t] * yudagI[v][r]) * WC6R(Cqq1R, p, v, s, w)
                                    +(yddagR[v][t] * yudagR[w][r] - yddagI[v][t] * yudagI[w][r]) * WC6I(Cqq1I, s, v, p, w)
                                    +(yddagI[v][t] * yudagR[w][r] + yddagR[v][t] * yudagI[w][r]) * WC6R(Cqq1R, s, v, p, w)
                                    )

                                    + 12. * (
                                    + (yddagR[w][t] * yudagR[v][r] - yddagI[w][t] * yudagI[v][r]) * WC6I(Cqq3I, p, v, s, w)
                                    +(yddagI[w][t] * yudagR[v][r] + yddagR[w][t] * yudagI[v][r]) * WC6R(Cqq3R, p, v, s, w)
                                    +(yddagR[v][t] * yudagR[w][r] - yddagI[v][t] * yudagI[w][r]) * WC6I(Cqq3I, s, v, p, w)
                                    +(yddagI[v][t] * yudagR[w][r] + yddagR[v][t] * yudagI[w][r]) * WC6R(Cqq3R, s, v, p, w)
                                    )

                                    + 2. * (
                                    + (yddagR[p][w] * yudagR[v][r] - yddagI[p][w] * yudagI[v][r]) * WC7I(Cqd8I, s, v, w, t)
                                    + (yddagI[p][w] * yudagR[v][r] + yddagR[p][w] * yudagI[v][r]) * WC7R(Cqd8R, s, v, w, t)
                                    +(yddagR[v][t] * yudagR[s][w] - yddagI[v][t] * yudagI[s][w]) * WC7I(Cqu8I, p, v, w, r)
                                    +(yddagI[v][t] * yudagR[s][w] + yddagR[v][t] * yudagI[s][w]) * WC7R(Cqu8R, p, v, w, r)
                                    +(yddagR[p][w] * yudagR[s][v] - yddagI[p][w] * yudagI[s][v]) * WC7I(Cud8I, v, r, w, t)
                                    + (yddagI[p][w] * yudagR[s][v] + yddagR[p][w] * yudagI[s][v]) * WC7R(Cud8R, v, r, w, t)
                                    )

                                    - 4. * (
                                    + (yddagR[s][w] * yudagR[p][v] - yddagI[s][w] * yudagI[p][v]) * WC7I(Cud1I, v, r, w, t)
                                    +(yddagI[s][w] * yudagR[p][v] + yddagR[s][w] * yudagI[p][v]) * WC7R(Cud1R, v, r, w, t)
                                    )
                                    ; //RGE 2
                            count ++;
                            //Cquqd8R
                            f[c + count * NG * NG * NG * NG] += - (4. / NC) * (
                                    + (yddagR[p][w] * yudagR[v][r] - yddagI[p][w] * yudagI[v][r]) * WC7R(Cqd8R, s, v, w, t)
                                    -(yddagI[p][w] * yudagR[v][r] + yddagR[p][w] * yudagI[v][r]) * WC7I(Cqd8I, s, v, w, t)
                                    +(yddagR[v][t] * yudagR[s][w] - yddagI[v][t] * yudagI[s][w]) * WC7R(Cqu8R, p, v, w, r)
                                    -(yddagI[v][t] * yudagR[s][w] + yddagR[v][t] * yudagI[s][w]) * WC7I(Cqu8I, p, v, w, r)
                                    +(yddagR[p][w] * yudagR[s][v] - yddagI[p][w] * yudagI[s][v]) * WC7R(Cud8R, v, r, w, t)
                                    -(yddagI[p][w] * yudagR[s][v] + yddagR[p][w] * yudagI[s][v]) * WC7I(Cud8I, v, r, w, t)
                                    )

                                    + 8. * (
                                    + (yddagR[w][t] * yudagR[v][r] - yddagI[w][t] * yudagI[v][r]) * WC6R(Cqq1R, s, v, p, w)
                                    -(yddagI[w][t] * yudagR[v][r] + yddagR[w][t] * yudagI[v][r]) * WC6I(Cqq1I, s, v, p, w)
                                    +(yddagR[v][t] * yudagR[w][r] - yddagI[v][t] * yudagI[w][r]) * WC6R(Cqq1R, p, v, s, w)
                                    -(yddagI[v][t] * yudagR[w][r] + yddagR[v][t] * yudagI[w][r]) * WC6I(Cqq1I, p, v, s, w)
                                    )

                                    - 24. * (
                                    + (yddagR[w][t] * yudagR[v][r] - yddagI[w][t] * yudagI[v][r]) * WC6R(Cqq3R, s, v, p, w)
                                    -(yddagI[w][t] * yudagR[v][r] + yddagR[w][t] * yudagI[v][r]) * WC6I(Cqq3I, s, v, p, w)
                                    +(yddagR[v][t] * yudagR[w][r] - yddagI[v][t] * yudagI[w][r]) * WC6R(Cqq3R, p, v, s, w)
                                    -(yddagI[v][t] * yudagR[w][r] + yddagR[v][t] * yudagI[w][r]) * WC6I(Cqq3I, p, v, s, w)
                                    )

                                    + 8. * (
                                    + (yddagR[p][w] * yudagR[v][r] - yddagI[p][w] * yudagI[v][r]) * WC7R(Cqd1R, s, v, w, t)
                                    -(yddagI[p][w] * yudagR[v][r] + yddagR[p][w] * yudagI[v][r]) * WC7I(Cqd1I, s, v, w, t)
                                    +(yddagR[v][t] * yudagR[s][w] - yddagI[v][t] * yudagI[s][w]) * WC7R(Cqu1R, p, v, w, r)
                                    -(yddagI[v][t] * yudagR[s][w] + yddagR[v][t] * yudagI[s][w]) * WC7I(Cqu1I, p, v, w, r)
                                    +(yddagR[p][w] * yudagR[s][v] - yddagI[p][w] * yudagI[s][v]) * WC7R(Cud1R, v, r, w, t)
                                    -(yddagI[p][w] * yudagR[s][v] + yddagR[p][w] * yudagI[s][v]) * WC7I(Cud1I, v, r, w, t)
                                    )

                                    - 4. * (
                                    + (yddagR[s][w] * yudagR[p][v] - yddagI[s][w] * yudagI[p][v]) * WC7R(Cud8R, v, r, w, t)
                                    -(yddagI[s][w] * yudagR[p][v] + yddagR[s][w] * yudagI[p][v]) * WC7I(Cud8I, v, r, w, t)
                                    )
                                    ; //RGE 2
                            count ++;
                            //Cquqd8I
                            f[c + count * NG * NG * NG * NG] += - (4. / NC) * (
                                    + (yddagR[p][w] * yudagR[v][r] - yddagI[p][w] * yudagI[v][r]) * WC7I(Cqd8I, s, v, w, t)
                                    +(yddagI[p][w] * yudagR[v][r] + yddagR[p][w] * yudagI[v][r]) * WC7R(Cqd8R, s, v, w, t)
                                    +(yddagR[v][t] * yudagR[s][w] - yddagI[v][t] * yudagI[s][w]) * WC7I(Cqu8I, p, v, w, r)
                                    +(yddagI[v][t] * yudagR[s][w] + yddagR[v][t] * yudagI[s][w]) * WC7R(Cqu8R, p, v, w, r)
                                    +(yddagR[p][w] * yudagR[s][v] - yddagI[p][w] * yudagI[s][v]) * WC7I(Cud8I, v, r, w, t)
                                    +(yddagI[p][w] * yudagR[s][v] + yddagR[p][w] * yudagI[s][v]) * WC7R(Cud8R, v, r, w, t)
                                    )

                                    + 8. * (
                                    + (yddagR[w][t] * yudagR[v][r] - yddagI[w][t] * yudagI[v][r]) * WC6I(Cqq1I, s, v, p, w)
                                    +(yddagI[w][t] * yudagR[v][r] + yddagR[w][t] * yudagI[v][r]) * WC6R(Cqq1R, s, v, p, w)
                                    +(yddagR[v][t] * yudagR[w][r] - yddagI[v][t] * yudagI[w][r]) * WC6I(Cqq1I, p, v, s, w)
                                    +(yddagI[v][t] * yudagR[w][r] + yddagR[v][t] * yudagI[w][r]) * WC6R(Cqq1R, p, v, s, w)
                                    )

                                    - 24. * (
                                    + (yddagR[w][t] * yudagR[v][r] - yddagI[w][t] * yudagI[v][r]) * WC6I(Cqq3I, s, v, p, w)
                                    +(yddagI[w][t] * yudagR[v][r] + yddagR[w][t] * yudagI[v][r]) * WC6R(Cqq3R, s, v, p, w)
                                    +(yddagR[v][t] * yudagR[w][r] - yddagI[v][t] * yudagI[w][r]) * WC6I(Cqq3I, p, v, s, w)
                                    +(yddagI[v][t] * yudagR[w][r] + yddagR[v][t] * yudagI[w][r]) * WC6R(Cqq3R, p, v, s, w)
                                    )

                                    + 8. * (
                                    + (yddagR[p][w] * yudagR[v][r] - yddagI[p][w] * yudagI[v][r]) * WC7I(Cqd1I, s, v, w, t)
                                    +(yddagI[p][w] * yudagR[v][r] + yddagR[p][w] * yudagI[v][r]) * WC7R(Cqd1R, s, v, w, t)
                                    +(yddagR[v][t] * yudagR[s][w] - yddagI[v][t] * yudagI[s][w]) * WC7I(Cqu1I, p, v, w, r)
                                    +(yddagI[v][t] * yudagR[s][w] + yddagR[v][t] * yudagI[s][w]) * WC7R(Cqu1R, p, v, w, r)
                                    +(yddagR[p][w] * yudagR[s][v] - yddagI[p][w] * yudagI[s][v]) * WC7I(Cud1I, v, r, w, t)
                                    +(yddagI[p][w] * yudagR[s][v] + yddagR[p][w] * yudagI[s][v]) * WC7R(Cud1R, v, r, w, t)
                                    )

                                    - 4. * (
                                    + (yddagR[s][w] * yudagR[p][v] - yddagI[s][w] * yudagI[p][v]) * WC7I(Cud8I, v, r, w, t)
                                    +(yddagI[s][w] * yudagR[p][v] + yddagR[s][w] * yudagI[p][v]) * WC7R(Cud8R, v, r, w, t)
                                    )
                                    ; //RGE 2
                            count ++;
                            //Clequ1R
                            f[c + count * NG * NG * NG * NG] += 2. * (
                                    + (yddagR[s][v] * yudagR[w][t] - yddagI[s][v] * yudagI[w][t]) * WC5(CledqR, p, r, v, w)
                                    -(yddagI[s][v] * yudagR[w][t] + yddagR[s][v] * yudagI[w][t]) * WC5(CledqI, p, r, v, w)
                                    )

                                    + 2. * (
                                    + (yedagR[p][v] * yudagR[s][w] - yedagI[p][v] * yudagI[s][w]) * WC7R(CeuR, v, r, w, t)
                                    -(yedagI[p][v] * yudagR[s][w] + yedagR[p][v] * yudagI[s][w]) * WC7I(CeuI, v, r, w, t)
                                    )

                                    + 2. * (
                                    + (yedagR[v][r] * yudagR[w][t] - yedagI[v][r] * yudagI[w][t]) * WC7R(Clq1R, p, v, s, w)
                                    -(yedagI[v][r] * yudagR[w][t] + yedagR[v][r] * yudagI[w][t]) * WC7I(Clq1I, p, v, s, w)
                                    )

                                    - 6. * (
                                    + (yedagR[v][r] * yudagR[w][t] - yedagI[v][r] * yudagI[w][t]) * WC7R(Clq3R, p, v, s, w)
                                    -(yedagI[v][r] * yudagR[w][t] + yedagR[v][r] * yudagI[w][t]) * WC7I(Clq3I, p, v, s, w)
                                    )

                                    - 2. * (
                                    + (yedagR[v][r] * yudagR[s][w] - yedagI[v][r] * yudagI[s][w]) * WC7R(CluR, p, v, w, t)
                                    -(yedagI[v][r] * yudagR[s][w] + yedagR[v][r] * yudagI[s][w]) * WC7I(CluI, p, v, w, t)
                                    )

                                    - 2. * (
                                    + (yedagR[p][w] * yudagR[v][t] - yedagI[p][w] * yudagI[v][t]) * WC7R(CqeR, s, v, w, r)
                                    -(yedagI[p][w] * yudagR[v][t] + yedagR[p][w] * yudagI[v][t]) * WC7I(CqeI, s, v, w, r)
                                    )
                                    ; //RGE 2
                            count ++;
                            //Clequ1I
                            f[c + count * NG * NG * NG * NG] += 2. * (
                                    + (yddagR[s][v] * yudagR[w][t] - yddagI[s][v] * yudagI[w][t]) * WC5(CledqI, p, r, v, w)
                                    +(yddagI[s][v] * yudagR[w][t] + yddagR[s][v] * yudagI[w][t]) * WC5(CledqR, p, r, v, w)
                                    )

                                    + 2. * (
                                    + (yedagR[p][v] * yudagR[s][w] - yedagI[p][v] * yudagI[s][w]) * WC7I(CeuI, v, r, w, t)
                                    +(yedagI[p][v] * yudagR[s][w] + yedagR[p][v] * yudagI[s][w]) * WC7R(CeuR, v, r, w, t)
                                    )

                                    + 2. * (
                                    + (yedagR[v][r] * yudagR[w][t] - yedagI[v][r] * yudagI[w][t]) * WC7I(Clq1I, p, v, s, w)
                                    +(yedagI[v][r] * yudagR[w][t] + yedagR[v][r] * yudagI[w][t]) * WC7R(Clq1R, p, v, s, w)
                                    )

                                    - 6. * (
                                    + (yedagR[v][r] * yudagR[w][t] - yedagI[v][r] * yudagI[w][t]) * WC7I(Clq3I, p, v, s, w)
                                    +(yedagI[v][r] * yudagR[w][t] + yedagR[v][r] * yudagI[w][t]) * WC7R(Clq3R, p, v, s, w)
                                    )

                                    - 2. * (
                                    + (yedagR[v][r] * yudagR[s][w] - yedagI[v][r] * yudagI[s][w]) * WC7I(CluI, p, v, w, t)
                                    +(yedagI[v][r] * yudagR[s][w] + yedagR[v][r] * yudagI[s][w]) * WC7R(CluR, p, v, w, t)
                                    )

                                    - 2. * (
                                    + (yedagR[p][w] * yudagR[v][t] - yedagI[p][w] * yudagI[v][t]) * WC7I(CqeI, s, v, w, r)
                                    +(yedagI[p][w] * yudagR[v][t] + yedagR[p][w] * yudagI[v][t]) * WC7R(CqeR, s, v, w, r)
                                    )
                                    ; //RGE 2	
                            count ++;
                            //Clequ3R
                            f[c + count * NG * NG * NG * NG] += - 0.5 * (
                                    + (yudagR[s][w] * yedagR[p][v] - yudagI[s][w] * yedagI[p][v]) * WC7R(CeuR, v, r, w, t)
                                    -(yudagI[s][w] * yedagR[p][v] + yudagR[s][w] * yedagI[p][v]) * WC7I(CeuI, v, r, w, t)
                                    )
                                    - 0.5 * (
                                    + (yedagR[v][r] * yudagR[w][t] - yedagI[v][r] * yudagI[w][t]) * WC7R(Clq1R, p, v, s, w)
                                    -(yedagI[v][r] * yudagR[w][t] + yedagR[v][r] * yudagI[w][t]) * WC7I(Clq1I, p, v, s, w)
                                    )
                                    + 1.5 * (
                                    + (yedagR[v][r] * yudagR[w][t] - yedagI[v][r] * yudagI[w][t]) * WC7R(Clq3R, p, v, s, w)
                                    -(yedagI[v][r] * yudagR[w][t] + yedagR[v][r] * yudagI[w][t]) * WC7I(Clq3I, p, v, s, w)
                                    )
                                    - 0.5 * (
                                    + (yedagR[v][r] * yudagR[s][w] - yedagI[v][r] * yudagI[s][w]) * WC7R(CluR, p, v, w, t)
                                    -(yedagI[v][r] * yudagR[s][w] + yedagR[v][r] * yudagI[s][w]) * WC7I(CluI, p, v, w, t)
                                    )
                                    - 0.5 * (
                                    + (yedagR[p][w] * yudagR[v][t] - yedagI[p][w] * yudagI[v][t]) * WC7R(CqeR, s, v, w, r)
                                    -(yedagI[p][w] * yudagR[v][t] + yedagR[p][w] * yudagI[v][t]) * WC7I(CqeI, s, v, w, r)
                                    )
                                    ; //RGE 2
                            count ++;
                            //Clequ3I
                            f[c + count * NG * NG * NG * NG] += - 0.5 * (
                                    + (yudagR[s][w] * yedagR[p][v] - yudagI[s][w] * yedagI[p][v]) * WC7I(CeuI, v, r, w, t)
                                    +(yudagI[s][w] * yedagR[p][v] + yudagR[s][w] * yedagI[p][v]) * WC7R(CeuR, v, r, w, t)
                                    )
                                    - 0.5 * (
                                    + (yedagR[v][r] * yudagR[w][t] - yedagI[v][r] * yudagI[w][t]) * WC7I(Clq1I, p, v, s, w)
                                    +(yedagI[v][r] * yudagR[w][t] + yedagR[v][r] * yudagI[w][t]) * WC7R(Clq1R, p, v, s, w)
                                    )
                                    + 1.5 * (
                                    + (yedagR[v][r] * yudagR[w][t] - yedagI[v][r] * yudagI[w][t]) * WC7I(Clq3I, p, v, s, w)
                                    +(yedagI[v][r] * yudagR[w][t] + yedagR[v][r] * yudagI[w][t]) * WC7R(Clq3R, p, v, s, w)
                                    )
                                    - 0.5 * (
                                    + (yedagR[v][r] * yudagR[s][w] - yedagI[v][r] * yudagI[s][w]) * WC7I(CluI, p, v, w, t)
                                    +(yedagI[v][r] * yudagR[s][w] + yedagR[v][r] * yudagI[s][w]) * WC7R(CluR, p, v, w, t)
                                    )
                                    - 0.5 * (
                                    + (yedagR[p][w] * yudagR[v][t] - yedagI[p][w] * yudagI[v][t]) * WC7I(CqeI, s, v, w, r)
                                    +(yedagI[p][w] * yudagR[v][t] + yedagR[p][w] * yudagI[v][t]) * WC7R(CqeR, s, v, w, r)
                                    )
                                    ; //RGE 2
                            count ++;
                        }
                    }

                    for (count = 0; count < N8_LRLR * 2; count ++) {
                        f[c + count * NG * NG * NG * NG] *= loop_factor;
                    }
                    c ++;
                }
            }
        }
    }
    c += NG * NG * NG * NG * (2 * N8_LRLR - 1);


    //INDICES CHECK
    {

        //In fact decommenting this piece sets all the entries of the ADM to 0
        //so one should verify that all the initial condition
        //stay the same under the evolution.
        //If this happen, the switch from different structures 
        //in main (complex numbers, matrices etc) to an array of doubles
        //has been correctly managed  
        /*
                        for (i = 0; i < dim; i++) {
                                f[i]=0.;
                        }
         */
    }



    return GSL_SUCCESS;
}
