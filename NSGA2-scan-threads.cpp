#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "string.h" 
#include "direct.h" 
#include "NSGA2.hpp"
#include <omp.h>



double aerosol(int n, double wavenumber)
{
    double x, b0, b1;
    int i, j;

    x = 1.0e7 / wavenumber - 2000.0;
    if (x < 0.0)
        x = 0.0;
    if (x > 12000.0)
        x = 12000.0;

    i = (int)x;
    b1 = x - (i + 0.0);
    j = i + 1;
    if (j > 12000)
        j = 12000;

    b0 = (1.0 - b1) * aerosol_35[n][i] + b1 * aerosol_35[n][j];

    return(b0);
}

double EBwavenumber(double wavenumber1, double wavenumber2, double t0)
{
    double a, c1, c2;
    double b1, b2, b3, b4, b5;
    double g0[9][2];
    int i;

    c1 = 3.7418e-16;
    c2 = 1.4388e-2;

    g0[0][0] = -0.9681602395; g0[0][1] = 0.0812743884;
    g0[1][0] = -0.8360311073; g0[1][1] = 0.1806481697;
    g0[2][0] = -0.6133714327; g0[2][1] = 0.2606106964;
    g0[3][0] = -0.3242534234; g0[3][1] = 0.3123470770;
    g0[4][0] = 0.0; g0[4][1] = 0.3302393550;
    g0[5][0] = 0.3242534234; g0[5][1] = 0.3123470770;
    g0[6][0] = 0.6133714327; g0[6][1] = 0.2606106964;
    g0[7][0] = 0.8360311073; g0[7][1] = 0.1806481697;
    g0[8][0] = 0.9681602395; g0[8][1] = 0.0812743884;

    b1 = c2 * 100.0 * wavenumber1 / t0;
    b2 = c2 * 100.0 * wavenumber2 / t0;
    for (a = 0.0, i = 0; i < 9; i++)
    {
        b3 = (b1 + b2) / 2.0 + (b2 - b1) / 2.0 * g0[i][0];
        if (b3 < 200.0)
            a += g0[i][1] * b3 * b3 * b3 / (exp(b3) - 1.0);
    }
    b3 = t0 / c2;
    a = (b2 - b1) / 2.0 * a * c1 * b3 * b3 * b3 * b3;

    return(a);
}

double E3(double L)
{
    double b4;

    if (L > 200.0)
        L = 200.0;
    b4 = exp(-L);

    return(b4);
}

void test1(double gi[20][3], double gim[20][20], int NGL)
{
    int i, j, k, n;
    double b1, b2, b3;
    double gm[20][20][2];



    if (NGL == 20)
    {
        gi[19][0] = 0.998238; gi[19][1] = 0.004521;
        gi[18][0] = 0.990726; gi[18][1] = 0.010498;
        gi[17][0] = 0.977260; gi[17][1] = 0.016421;
        gi[16][0] = 0.957917; gi[16][1] = 0.022246;
        gi[15][0] = 0.932813; gi[15][1] = 0.027937;
        gi[14][0] = 0.902099; gi[14][1] = 0.033460;
        gi[13][0] = 0.865960; gi[13][1] = 0.038782;
        gi[12][0] = 0.824612; gi[12][1] = 0.043871;
        gi[11][0] = 0.778306; gi[11][1] = 0.048696;
        gi[10][0] = 0.727318; gi[10][1] = 0.053228;
        gi[9][0] = 0.671957;  gi[9][1] = 0.057440;
        gi[8][0] = 0.612554;  gi[8][1] = 0.061306;
        gi[7][0] = 0.549467;  gi[7][1] = 0.064804;
        gi[6][0] = 0.483076;  gi[6][1] = 0.067912;
        gi[5][0] = 0.413779;  gi[5][1] = 0.070612;
        gi[4][0] = 0.341994;  gi[4][1] = 0.072887;
        gi[3][0] = 0.268152;  gi[3][1] = 0.074723;
        gi[2][0] = 0.192698;  gi[2][1] = 0.076110;
        gi[1][0] = 0.116084;  gi[1][1] = 0.077040;
        gi[0][0] = 0.0387724; gi[0][1] = 0.077506;
    }
    if (NGL == 19)
    {
        gi[18][0] = 0.99805;	gi[18][1] = 0.005;
        gi[17][0] = 0.98974;	gi[17][1] = 0.011613;
        gi[16][0] = 0.97485;	gi[16][1] = 0.018157;
        gi[15][0] = 0.953466;	gi[15][1] = 0.02458;
        gi[14][0] = 0.925741;	gi[14][1] = 0.03084;
        gi[13][0] = 0.891856;	gi[13][1] = 0.0368941;
        gi[12][0] = 0.852035;	gi[12][1] = 0.0427;
        gi[11][0] = 0.806544;	gi[11][1] = 0.04823;
        gi[10][0] = 0.755686;	gi[10][1] = 0.05343;
        gi[9][0] = 0.699799;	gi[9][1] = 0.05828;
        gi[8][0] = 0.639254;	gi[8][1] = 0.062741;
        gi[7][0] = 0.574456;	gi[7][1] = 0.066784;
        gi[6][0] = 0.505835;	gi[6][1] = 0.070383;
        gi[5][0] = 0.43385;	  gi[5][1] = 0.073513;
        gi[4][0] = 0.358972;	gi[4][1] = 0.07615;
        gi[3][0] = 0.281709;	gi[3][1] = 0.07829;
        gi[2][0] = 0.202571;	gi[2][1] = 0.079901;
        gi[1][0] = 0.122084;	gi[1][1] = 0.080982;
        gi[0][0] = 0.040785;	gi[0][1] = 0.081525;
    }
    if (NGL == 18)
    {
        gi[17][0] = 0.997831;	gi[17][1] = 0.005566;
        gi[16][0] = 0.988587;	gi[16][1] = 0.012916;
        gi[15][0] = 0.972028;	gi[15][1] = 0.020182;
        gi[14][0] = 0.94827;	gi[14][1] = 0.027299;
        gi[13][0] = 0.917498;	gi[13][1] = 0.034214;
        gi[12][0] = 0.87993;	gi[12][1] = 0.040876;
        gi[11][0] = 0.835847;	gi[11][1] = 0.047235;
        gi[10][0] = 0.785576;	gi[10][1] = 0.053245;
        gi[9][0] = 0.729489;	gi[9][1] = 0.05886;
        gi[8][0] = 0.668001;	gi[8][1] = 0.06404;
        gi[7][0] = 0.60157;	  gi[7][1] = 0.068745;
        gi[6][0] = 0.53068;	  gi[6][1] = 0.072942;
        gi[5][0] = 0.455864;	gi[5][1] = 0.076598;
        gi[4][0] = 0.377673;	gi[4][1] = 0.079688;
        gi[3][0] = 0.296685;	gi[3][1] = 0.08219;
        gi[2][0] = 0.213501;	gi[2][1] = 0.08408;
        gi[1][0] = 0.128736;	gi[1][1] = 0.08535;
        gi[0][0] = 0.0430182;	gi[0][1] = 0.085983;
    }
    if (NGL == 17)
    {
        gi[16][0] = 0.997572;	gi[16][1] = 0.006229;
        gi[15][0] = 0.987228;	gi[15][1] = 0.01445;
        gi[14][0] = 0.968708;	gi[14][1] = 0.022564;
        gi[13][0] = 0.942162;	gi[13][1] = 0.030491;
        gi[12][0] = 0.90781;	gi[12][1] = 0.038167;
        gi[11][0] = 0.865935;	gi[11][1] = 0.045526;
        gi[10][0] = 0.816884;	gi[10][1] = 0.052507;
        gi[9][0] = 0.761065;	gi[9][1] = 0.05905;
        gi[8][0] = 0.698939;	gi[8][1] = 0.06511;
        gi[7][0] = 0.631022;	gi[7][1] = 0.070629;
        gi[6][0] = 0.557876;	gi[6][1] = 0.075562;
        gi[5][0] = 0.480107;	gi[5][1] = 0.079868;
        gi[4][0] = 0.398359;	gi[4][1] = 0.083513;
        gi[3][0] = 0.313311;	gi[3][1] = 0.0864657;
        gi[2][0] = 0.225667;	gi[2][1] = 0.088702;
        gi[1][0] = 0.136152;	gi[1][1] = 0.090203;
        gi[0][0] = 0.04551;	  gi[0][1] = 0.090957;
    }
    if (NGL == 16)
    {
        gi[15][0] = 0.997264;	gi[15][1] = 0.007019;
        gi[14][0] = 0.985612;	gi[14][1] = 0.016274;
        gi[13][0] = 0.964762;	gi[13][1] = 0.025392;
        gi[12][0] = 0.934906;	gi[12][1] = 0.034274;
        gi[11][0] = 0.896321;	gi[11][1] = 0.042836;
        gi[10][0] = 0.84937;	gi[10][1] = 0.050998;
        gi[9][0] = 0.794484;	gi[9][1] = 0.058684;
        gi[8][0] = 0.73218;	  gi[8][1] = 0.065822;
        gi[7][0] = 0.663044;	gi[7][1] = 0.072346;
        gi[6][0] = 0.587716;	gi[6][1] = 0.078194;
        gi[5][0] = 0.5069;	  gi[5][1] = 0.08331;
        gi[4][0] = 0.421351;	gi[4][1] = 0.08765;
        gi[3][0] = 0.331869;	gi[3][1] = 0.091174;
        gi[2][0] = 0.239287;	gi[2][1] = 0.093844;
        gi[1][0] = 0.144472;	gi[1][1] = 0.095639;
        gi[0][0] = 0.048308;	gi[0][1] = 0.09654;
    }
    if (NGL == 15)
    {
        gi[14][0] = 0.996894;	gi[14][1] = 0.00797;
        gi[13][0] = 0.983668;	gi[13][1] = 0.018466;
        gi[12][0] = 0.960022;	gi[12][1] = 0.028785;
        gi[11][0] = 0.9262;	  gi[11][1] = 0.038799;
        gi[10][0] = 0.882561;	gi[10][1] = 0.048403;
        gi[9][0] = 0.829566;	gi[9][1] = 0.057493;
        gi[8][0] = 0.767777;	gi[8][1] = 0.06597;
        gi[7][0] = 0.69785;	  gi[7][1] = 0.073756;
        gi[6][0] = 0.620526;	gi[6][1] = 0.08076;
        gi[5][0] = 0.536624;	gi[5][1] = 0.0869;
        gi[4][0] = 0.447034;	gi[4][1] = 0.092123;
        gi[3][0] = 0.352705;	gi[3][1] = 0.096369;
        gi[2][0] = 0.254637;	gi[2][1] = 0.099593;
        gi[1][0] = 0.15387;	  gi[1][1] = 0.101762;
        gi[0][0] = 0.051472;	gi[0][1] = 0.102853;
    }
    if (NGL == 14)
    {
        gi[13][0] = 0.996443;	gi[13][1] = 0.009124;
        gi[12][0] = 0.981303;	gi[12][1] = 0.021132;
        gi[11][0] = 0.954259;	gi[11][1] = 0.032901;
        gi[10][0] = 0.915633;	gi[10][1] = 0.04427;
        gi[9][0] = 0.865893;	gi[9][1] = 0.055107;
        gi[8][0] = 0.805641;	gi[8][1] = 0.065273;
        gi[7][0] = 0.73561;	  gi[7][1] = 0.074646;
        gi[6][0] = 0.656651;	gi[6][1] = 0.0831134;
        gi[5][0] = 0.569721;	gi[5][1] = 0.090572;
        gi[4][0] = 0.475874;	gi[4][1] = 0.096931;
        gi[3][0] = 0.376252;	gi[3][1] = 0.102113;
        gi[2][0] = 0.272062;	gi[2][1] = 0.106056;
        gi[1][0] = 0.164569;	gi[1][1] = 0.108711;
        gi[0][0] = 0.05508;	  gi[0][1] = 0.110047;
    }
    if (NGL == 13)
    {
        gi[12][0] = 0.995886;	gi[12][1] = 0.010551;
        gi[11][0] = 0.978385;	gi[11][1] = 0.0244179;
        gi[10][0] = 0.947159;	gi[10][1] = 0.0379624;
        gi[9][0] = 0.902638;	gi[9][1] = 0.050976;
        gi[8][0] = 0.845446;	gi[8][1] = 0.063274;
        gi[7][0] = 0.776386;	gi[7][1] = 0.074684;
        gi[6][0] = 0.696427;	gi[6][1] = 0.085046;
        gi[5][0] = 0.606692;	gi[5][1] = 0.094214;
        gi[4][0] = 0.508441;	gi[4][1] = 0.10206;
        gi[3][0] = 0.403052;	gi[3][1] = 0.108472;
        gi[2][0] = 0.292005;	gi[2][1] = 0.113362;
        gi[1][0] = 0.176859;	gi[1][1] = 0.11666;
        gi[0][0] = 0.05923;	  gi[0][1] = 0.11832;
    }
    if (NGL == 12)
    {
        gi[11][0] = 0.995187;	gi[11][1] = 0.012341;
        gi[10][0] = 0.974729;	gi[10][1] = 0.028531;
        gi[9][0] = 0.938275;	gi[9][1] = 0.044277;
        gi[8][0] = 0.886416;	gi[8][1] = 0.059299;
        gi[7][0] = 0.82;	    gi[7][1] = 0.07335;
        gi[6][0] = 0.740124;	gi[6][1] = 0.08619;
        gi[5][0] = 0.648094;	gi[5][1] = 0.097619;
        gi[4][0] = 0.54542;	  gi[4][1] = 0.107444;
        gi[3][0] = 0.433794;	gi[3][1] = 0.115506;
        gi[2][0] = 0.315043;	gi[2][1] = 0.121671;
        gi[1][0] = 0.191119;	gi[1][1] = 0.125838;
        gi[0][0] = 0.064057;	gi[0][1] = 0.127938;
    }
    if (NGL == 11)
    {
        gi[10][0] = 0.994295;	gi[10][1] = 0.014628;
        gi[9][0] = 0.970061;	gi[9][1] = 0.033775;
        gi[8][0] = 0.926957;	gi[8][1] = 0.052293;
        gi[7][0] = 0.865813;	gi[7][1] = 0.069796;
        gi[6][0] = 0.787817;	gi[6][1] = 0.0859416;
        gi[5][0] = 0.694487;	gi[5][1] = 0.100414;
        gi[4][0] = 0.58764;	  gi[4][1] = 0.112932;
        gi[3][0] = 0.469356;	gi[3][1] = 0.123252;
        gi[2][0] = 0.341936;	gi[2][1] = 0.131174;
        gi[1][0] = 0.20786;	  gi[1][1] = 0.13654;
        gi[0][0] = 0.06974;	  gi[0][1] = 0.13925;
    }
    if (NGL == 10)
    {
        gi[9][0] = 0.993129;	gi[9][1] = 0.017614;
        gi[8][0] = 0.963972;	gi[8][1] = 0.040601;
        gi[7][0] = 0.912234;	gi[7][1] = 0.062672;
        gi[6][0] = 0.839117;	gi[6][1] = 0.083277;
        gi[5][0] = 0.74633;	  gi[5][1] = 0.10193;
        gi[4][0] = 0.636054;	gi[4][1] = 0.118195;
        gi[3][0] = 0.510867;	gi[3][1] = 0.131689;
        gi[2][0] = 0.373706;	gi[2][1] = 0.142096;
        gi[1][0] = 0.227786;	gi[1][1] = 0.14917;
        gi[0][0] = 0.076527;	gi[0][1] = 0.152753;
    }
    if (NGL == 9)
    {
        gi[8][0] = 0.991565;	gi[8][1] = 0.021616;
        gi[7][0] = 0.955824;	gi[7][1] = 0.04971;
        gi[6][0] = 0.892603;	gi[6][1] = 0.0764257;
        gi[5][0] = 0.80371;	  gi[5][1] = 0.100942;
        gi[4][0] = 0.691687;	gi[4][1] = 0.122555;
        gi[3][0] = 0.559771;	gi[3][1] = 0.140643;
        gi[2][0] = 0.411751;	gi[2][1] = 0.154685;
        gi[1][0] = 0.251886;	gi[1][1] = 0.164277;
        gi[0][0] = 0.084775;	gi[0][1] = 0.169142;
    }
    if (NGL == 8)
    {
        gi[7][0] = 0.989401;	gi[7][1] = 0.027152;
        gi[6][0] = 0.944575;	gi[6][1] = 0.062254;
        gi[5][0] = 0.865631;	gi[5][1] = 0.095159;
        gi[4][0] = 0.755404;	gi[4][1] = 0.124629;
        gi[3][0] = 0.617876;	gi[3][1] = 0.149596;
        gi[2][0] = 0.458017;	gi[2][1] = 0.169157;
        gi[1][0] = 0.281604;	gi[1][1] = 0.182603;
        gi[0][0] = 0.095013;	gi[0][1] = 0.18945;
    }
    if (NGL == 7)
    {
        gi[6][0] = 0.98628;	  gi[6][1] = 0.0351195;
        gi[5][0] = 0.928435;	gi[5][1] = 0.080158;
        gi[4][0] = 0.827201;	gi[4][1] = 0.12152;
        gi[3][0] = 0.687293;	gi[3][1] = 0.1572;
        gi[2][0] = 0.515249;	gi[2][1] = 0.185538;
        gi[1][0] = 0.319112;	gi[1][1] = 0.205199;
        gi[0][0] = 0.108055;	gi[0][1] = 0.215264;
    }
    if (NGL == 6)
    {
        gi[5][0] = 0.981561;	gi[5][1] = 0.047175;
        gi[4][0] = 0.904117;	gi[4][1] = 0.106939;
        gi[3][0] = 0.769903;	gi[3][1] = 0.160078;
        gi[2][0] = 0.587318;	gi[2][1] = 0.203167;
        gi[1][0] = 0.367832;	gi[1][1] = 0.23349;
        gi[0][0] = 0.125233;	gi[0][1] = 0.24915;
    }
    if (NGL == 5)
    {
        gi[4][0] = 0.973907;	gi[4][1] = 0.066671;
        gi[3][0] = 0.865063;	gi[3][1] = 0.149451;
        gi[2][0] = 0.67941;	  gi[2][1] = 0.219086;
        gi[1][0] = 0.433395;	gi[1][1] = 0.26927;
        gi[0][0] = 0.148874;	gi[0][1] = 0.295524;
    }
    if (NGL == 4)
    {
        gi[3][0] = 0.96029;	  gi[3][1] = 0.101229;
        gi[2][0] = 0.796667;	gi[2][1] = 0.222381;
        gi[1][0] = 0.525532;	gi[1][1] = 0.313707;
        gi[0][0] = 0.183435;	gi[0][1] = 0.362684;
    }
    if (NGL == 3)
    {
        gi[2][0] = 0.93247;	  gi[2][1] = 0.171325;
        gi[1][0] = 0.661209;	gi[1][1] = 0.360762;
        gi[0][0] = 0.238619;	gi[0][1] = 0.467914;
    }



    for (i = 0; i < NGL; i++)
        for (j = 0; j < NGL; j++)
        {
            gm[i][j][0] = pow(gi[j][0], 2 * i);
            if (i == j)
                gm[i][j][1] = 1.0;
            else
                gm[i][j][1] = 0.0;
        }

    for (n = 0; n < NGL; n++)
    {
        b1 = gm[n][n][0];
        for (j = 0; j < NGL; j++)
        {
            gm[n][j][0] /= b1;
            gm[n][j][1] /= b1;
        }
        for (i = 0; i < NGL; i++)
            if (i != n)
            {
                b2 = gm[i][n][0];
                for (j = 0; j < NGL; j++)
                {
                    gm[i][j][0] -= b2 * gm[n][j][0];
                    gm[i][j][1] -= b2 * gm[n][j][1];
                }
            }
    }

    for (i = 0; i < NGL; i++)
        for (j = 0; j < NGL; j++)
            gim[i][j] = gm[i][j][1];
}

void gka0(double mole_concentration1, double mole_concentration2, double pressure, double temperature, int aerosol_type, double density_aerosol, double cutoff_wavenumber, double wavenumber_interval, FILE* fp1, int aa)
{
    int i, j, k, m, n, n1, n2, ns, mm;
    int k1, k2, m1, m2, m3, m4, m5;
    double b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11;
    int interval_number;
    double** interval_K, ** interval_K1, * interval_K2, ** interval_KT2;



    interval_number = (int)(upper_wave_number / wavenumber_interval) + 1;
    interval_K = (double**)calloc(interval_number, sizeof(double*));
    interval_K1 = (double**)calloc(interval_number, sizeof(double*));
    interval_K2 = (double*)calloc(interval_number, sizeof(double));
    for (i = 0; i < interval_number; i++)
    {
        interval_K[i] = (double*)calloc(7, sizeof(double));
        interval_K1[i] = (double*)calloc(parallel_threads_number, sizeof(double));
    }
    for (i = 0; i < interval_number; i++)
        interval_K[i][0] = wavenumber_interval * (i + 0.5);
    for (n = 0; n < 2 + 4 * aa; n++)
    {
        fread(interval_K2, sizeof(double), interval_number, fp1);
        for (k = 0; k < interval_number; k++)
        {

            interval_K[k][1 + n] = interval_K2[k];


            if (interval_K[k][1 + n] < 1.0e-10 && n < 2)
                interval_K[k][1 + n] = 1.0e-10;
        }
    }
    interval_KT2 = (double**)calloc(interval_number1, sizeof(double*));
    for (i = 0; i < interval_number1; i++)
        interval_KT2[i] = (double*)calloc(2, sizeof(double));



    for (ns = 0; ns < 2; ns++)
    {
#pragma omp parallel for schedule(guided) private(j,k,b1,b2,b4,b5,b6,b7)
        for (i = interval_number0; i < interval_number1; i++)
        {
            b1 = (lower_wave_number + dyita * (i + 0.5)) / wavenumber_interval - 0.5;
            if (b1 >= interval_number - 1.0)
            {
                interval_KT2[i][0] = interval_K[interval_number - 1][1 + ns];
                interval_KT2[i][1] = mole_concentration1 * interval_K[interval_number - 1][1] + mole_concentration2 * interval_K[interval_number - 1][2];
                for (k = 2; k < 2 + 4 * aa; k++)
                    interval_KT2[i][1] += mole_concentration[k] * interval_K[interval_number - 1][1 + k];
            }
            else
            {
                j = (int)b1;
                b2 = b1 - (j + 0.0);
                interval_KT2[i][0] = (1.0 - b2) * interval_K[j][1 + ns] + b2 * interval_K[j + 1][1 + ns];
                interval_KT2[i][1] = mole_concentration1 * ((1.0 - b2) * interval_K[j][1] + b2 * interval_K[j + 1][1]) + mole_concentration2 * ((1.0 - b2) * interval_K[j][2] + b2 * interval_K[j + 1][2]);
                for (k = 2; k < 2 + 4 * aa; k++)
                    interval_KT2[i][1] += mole_concentration[k] * ((1.0 - b2) * interval_K[j][1 + k] + b2 * interval_K[j + 1][1 + k]);
            }
            if (aerosol_type > 0)
                interval_KT2[i][1] += density_aerosol * aerosol(aerosol_type - 1, lower_wave_number + dyita * (i + 0.5)) / pressure;
            if (ns == 0)
                interval_KT2[i][1] = interval_KT2[i][1] / mole_concentration1;
            else
                interval_KT2[i][1] = interval_KT2[i][1] / mole_concentration2;

            for (j = 0; j < 2; j++)
                interval_KT1[NM][i][2 * ns + j] = interval_KT2[i][j];
        }
    }



    for (i = 0; i < interval_number; i++)
    {
        free(interval_K[i]);
        free(interval_K1[i]);
    }
    for (i = 0; i < interval_number1; i++)
        free(interval_KT2[i]);
    free(interval_K);
    free(interval_K1);
    free(interval_K2);
    free(interval_KT2);
    NM++;
}
void gka(double temperature, int nm, int ns, int NGL, double temperature0, double KA1[GN1][3][20], double KA2[GN2][3][20])
{
    int i, j, k, m, n, n1, n2, ng, mm;
    int k1, k2, m1, m2, m3, m4, m5;
    double b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11;
    double* ag, ** AG;
    double QKG[20], output_parameter[3][20];
    double gi[20][3], gim[20][20];
    double** interval_KT2;
    int thread_number;



    interval_KT2 = (double**)calloc(interval_number1, sizeof(double*));
    for (i = 0; i < interval_number1; i++)
        interval_KT2[i] = (double*)calloc(4, sizeof(double));
    ag = (double*)calloc(g_point_number2 + 1, sizeof(double));
    AG = (double**)calloc(g_point_number2 + 1, sizeof(double*));
    for (i = 0; i < g_point_number2 + 1; i++)
        AG[i] = (double*)calloc(2, sizeof(double));
    test1(gi, gim, NGL);



    for (i = interval_number0; i < interval_number1; i++)
    {
        interval_KT2[i][0] = interval_KT1[nm][i][2 * ns];
        interval_KT2[i][1] = interval_KT1[nm][i][2 * ns + 1];

        for (j = 0; j < 2; j++)
        {
            if (j == 0)
                b7 = temperature0;
            else
                b7 = temperature;

            b4 = EBwavenumber(lower_wave_number, upper_wave_number, b7);
            b5 = lower_wave_number + dyita * (i + 0.0);
            b6 = lower_wave_number + dyita * (i + 1.0);
            interval_KT2[i][2 + j] = EBwavenumber(b5, b6, b7) / b4;
        }
    }

    //*****************************************  

    for (ng = 0; ng < gn[ns]; ng++)
    {
        b1 = 1.0e12; b2 = 0.0;
        for (i = interval_number0; i < interval_number1; i++)
            if (interval_KG2[i][ns] == ng)
            {
                if (b1 > interval_KT2[i][0])
                    b1 = interval_KT2[i][0];
                if (b2 < interval_KT2[i][0])
                    b2 = interval_KT2[i][0];
            }

        for (i = 0; i < g_point_number2; i++)
        {
            b4 = (i + 1.0) / (g_point_number2 + 0.0);
            ag[i] = b1 * pow(b2 / b1, b4);
            for (j = 0; j < 2; j++)
                AG[i][j] = 0.0;
        }

        for (k = interval_number0; k < interval_number1; k++)
        {
            if (interval_KG2[k][ns] == ng)
            {
                i = 0;
                while ((interval_KT2[k][0] > ag[i + 10000]) && (i < g_point_number2 - 10000))
                    i += 10000;
                while ((interval_KT2[k][0] > ag[i + 1000]) && (i < g_point_number2 - 1000))
                    i += 1000;
                while ((interval_KT2[k][0] > ag[i + 100]) && (i < g_point_number2 - 100))
                    i += 100;
                while ((interval_KT2[k][0] > ag[i + 10]) && (i < g_point_number2 - 10))
                    i += 10;
                while ((interval_KT2[k][0] > ag[i]) && (i < g_point_number2 - 1))
                    i++;

                AG[i][0] += interval_KT2[k][2];
                AG[i][1] += interval_KT2[k][3];
            }
        }



        b10 = 0.0;
        for (i = 0; i < g_point_number2; i++)
            b10 += AG[i][0];
        for (i = 0; i < 20; i++)
            QKG[i] = 0.0;
        for (n = 0, i = 0; i < g_point_number2; i++)
        {
            if (AG[i][0] > 0.0)
            {
                ag[n] = ag[i];
                if (n == 0)
                {
                    AG[n][0] = AG[i][0];
                    AG[n][1] = AG[i][1];
                }
                else
                {
                    AG[n][0] = AG[n - 1][0] + AG[i][0];
                    AG[n][1] = AG[i][1];
                }
                b6 = AG[n][0] / b10;
                for (j = 0; j < NGL; j++)
                    QKG[j] += AG[n][1] * pow(b6, 2 * j);
                n++;
            }
        }



        for (j = 0; j < NGL; j++)
        {
            i = 0;
            while (AG[i + 10000][0] < gi[j][0] * b10 && i < n - 10000)
                i += 10000;
            while (AG[i + 1000][0] < gi[j][0] * b10 && i < n - 1000)
                i += 1000;
            while (AG[i + 100][0] < gi[j][0] * b10 && i < n - 100)
                i += 100;
            while (AG[i + 10][0] < gi[j][0] * b10 && i < n - 10)
                i += 10;
            while (AG[i][0] < gi[j][0] * b10 && i < n - 1)
                i++;

            if (i == 0)
            {
                b5 = (AG[i][0] - gi[j][0] * b10) / AG[i][0];
                output_parameter[0][j] = b1 * b5 + ag[0] * (1.0 - b5);
            }
            else
            {
                b5 = (AG[i][0] - gi[j][0] * b10) / (AG[i][0] - AG[i - 1][0]);
                if (b5 < 0.0)
                    b5 = 0.0;
                output_parameter[0][j] = ag[i - 1] * b5 + ag[i] * (1.0 - b5);
            }

            output_parameter[1][j] = 0.0;
            for (k = 0; k < NGL; k++)
                output_parameter[1][j] += QKG[k] * gim[j][k];
        }

        //*****************************************  

        b1 = 1.0e12; b2 = 0.0;
        for (i = interval_number0; i < interval_number1; i++)
            if (interval_KG2[i][ns] == ng)
            {
                if (b1 > interval_KT2[i][1])
                    b1 = interval_KT2[i][1];
                if (b2 < interval_KT2[i][1])
                    b2 = interval_KT2[i][1];
            }

        for (i = 0; i < g_point_number2; i++)
        {
            b4 = (i + 1.0) / (g_point_number2 + 0.0);
            ag[i] = b1 * pow(b2 / b1, b4);
            AG[i][0] = 0.0;
        }

        for (k = interval_number0; k < interval_number1; k++)
            if (interval_KG2[k][ns] == ng)
            {
                i = 0;
                while ((interval_KT2[k][1] > ag[i + 10000]) && (i < g_point_number2 - 10000))
                    i += 10000;
                while ((interval_KT2[k][1] > ag[i + 1000]) && (i < g_point_number2 - 1000))
                    i += 1000;
                while ((interval_KT2[k][1] > ag[i + 100]) && (i < g_point_number2 - 100))
                    i += 100;
                while ((interval_KT2[k][1] > ag[i + 10]) && (i < g_point_number2 - 10))
                    i += 10;
                while ((interval_KT2[k][1] > ag[i]) && (i < g_point_number2 - 1))
                    i++;

                AG[i][0] += interval_KT2[k][3] * interval_KT2[k][0];
            }

        n = 0;
        for (i = 0; i < g_point_number2; i++)
        {
            if (AG[i][0] > 0.0)
            {
                ag[n] = ag[i];
                if (n == 0)
                    AG[n][0] = AG[i][0];
                else
                    AG[n][0] = AG[n - 1][0] + AG[i][0];
                n++;
            }
        }

        for (j = 0; j < NGL; j++)
        {
            b1 = 0.0;
#pragma omp parallel for schedule(guided) reduction(+:b1)
            for (k = interval_number0; k < interval_number1; k++)
                if (interval_KG2[k][ns] == ng && interval_KT2[k][0] <= output_parameter[0][j])
                    b1 += interval_KT2[k][3] * interval_KT2[k][0];
            QKG[j] = b1;
        }



        for (m3 = 0; m3 < ((NGL + 1) / 2); m3++)
        {
            m5 = NGL - 2 * m3;
            if (m5 > 2)
                m5 = 2;
            for (m4 = 0; m4 < m5; m4++)
            {
                if (m4 == 0)
                    j = m3;
                else
                    j = NGL - 1 - m3;
                if (j == 0)
                {
                    m1 = 0;
                    m2 = n;
                }

                k1 = m1; k2 = m2;
                do
                {
                    k = (k1 + k2) / 2;
                    if (QKG[j] < AG[k][0])
                        k2 = k;
                    else
                        k1 = k;
                } while (k2 - k1 > 1);

                output_parameter[2][j] = ag[k2];
                if (m4 == 0)
                    m1 = k2;
                else
                    m2 = k2;
            }
        }



        if (ns == 0)
            for (i = 0; i < 3; i++)
                for (j = 0; j < NGL; j++)
                    KA1[ng][i][j] = output_parameter[i][j];
        if (ns == 1)
            for (i = 0; i < 3; i++)
                for (j = 0; j < NGL; j++)
                    KA2[ng][i][j] = output_parameter[i][j];
    }



    for (i = 0; i < g_point_number2 + 1; i++)
        free(AG[i]);
    free(AG);
    free(ag);
    for (i = 0; i < interval_number1; i++)
        free(interval_KT2[i]);
    free(interval_KT2);
}

//****************************************************************************************

#include "I-lh-lc-kal-msmgwb12.h"
#include "I-lh-lc-kal-msmgwb12w.h"

void test0()
{
	printf("test0\n");
    
    double b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10;
    int i, j, k, m, n, n1;
    double cutoff_wavenumber, wavenumber_interval;
    FILE* fp, * fp1;
    char c;
    char* c0;
    char c1[testN][300], c2[testN][300];

    printf("ppp0\n");
    
    double temperaturew;

    double Lh;
    double pressure1;
    double temperature1;
    double mole_concentration11;
    double mole_concentration12;
    double mole_concentration13;

    double Lc;
    double pressure2;
    double temperature2;
    double mole_concentration21;
    double mole_concentration22;
    double mole_concentration23;

    double dL;
    double pressure3;
    double temperature3;
    double mole_concentration31;
    double mole_concentration32;
    double mole_concentration33;
    double density_aerosol;

	printf("ppp1\n");

    interval_number0 = 0;
    interval_number1 = (int)((upper_wave_number - lower_wave_number) / dyita) + 1;
    interval_KG2 = (int**)calloc(interval_number1, sizeof(int*));
    for (i = 0; i < interval_number1; i++)
        interval_KG2[i] = (int*)calloc(2, sizeof(int));

	printf("ppp2\n");
    
    interval_KT1 = (double***)calloc(2 * testN + 2 * 14, sizeof(double**));
    for (i = 0; i < 2 * testN + 2 * 14; i++)
    {
        interval_KT1[i] = (double**)calloc(interval_number1, sizeof(double*));
        for (j = 0; j < interval_number1; j++)
            interval_KT1[i][j] = (double*)calloc(4, sizeof(double));
    }
    
	printf("ppp3\n");
    
    gn[0] = GN1;
    gn[1] = GN2;
    TP[0] = TP[9] = 100;
    TP[1] = TP[10] = 200;
    TP[2] = TP[11] = 300;
    TP[3] = TP[12] = 500;
    TP[4] = TP[13] = 800;
    TP[5] = TP[14] = 1100;
    TP[6] = TP[15] = 1500;
    TP[7] = TP[16] = 1900;
    TP[8] = TP[17] = 2500;

    printf("pp1\n");
    
    fp = fopen("atm_modes.txt", "r");
    for (j = 0; j < 4; j++)
        for (i = 0; i < 16; i++)
            fscanf(fp, "%lf%le%le%le%le%le", &PM[j][i][0], &PM[j][i][1], &PM[j][i][5], &PM[j][i][4], &PM[j][i][2], &PM[j][i][3]);
    fclose(fp);
    printf("大气模式参数读取完毕\n");

	printf("pp2\n");

    fp = fopen("sand_dust_0.4micron.txt", "r");
    for (i = 0; i <= 12000; i++)
        fscanf(fp, "%le%le", &b1, &aerosol_35[0][i]);
    fclose(fp);
    fp = fopen("sand_dust_4.0micron.txt", "r");
    for (i = 0; i <= 12000; i++)
        fscanf(fp, "%le%le", &b1, &aerosol_35[1][i]);
    fclose(fp);
    fp = fopen("sand_dust_10.0micron.txt", "r");
    for (i = 0; i <= 12000; i++)
        fscanf(fp, "%le%le", &b1, &aerosol_35[2][i]);
    fclose(fp);
    fp = fopen("soot_0.8micron.txt", "r");
    for (i = 0; i <= 12000; i++)
        fscanf(fp, "%le%le", &b1, &aerosol_35[3][i]);
    fclose(fp);
    printf("气溶胶消光系数读取完毕\n");

	printf("pp3\n");
    
    //******************
    fp = fopen("list.txt", "r");
    for (i = 0; i < testN; i++)
        fscanf(fp, "%s", c1[i]);
    fclose(fp);

    c0 = getcwd(NULL, 0);
    printf("当前文件夹为%s\n", c0);

    for (NM = 0, i = 0; i < testN; i++)
    {
        strcpy(c2[i], c0);
        strcat(c2[i], c1[i]);
        chdir(c2[i]);
        printf("任务%d文件夹为%s\n", i, c2[i]);
        NMn[i] = NM;


        if (i == 0)
        {
            pressure1 = 1.0;
            temperature1 = 800.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 1.0;
            temperature2 = 308.15;
            mole_concentration[0] = 0.05724;
            M5(2, pressure2);


        }
        if (i == 1)
        {
            pressure1 = 1.0;
            temperature1 = 600.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 1.0;
            temperature2 = 288.15;
            mole_concentration[0] = 0.00184;
            M5(0, pressure2);


        }
        if (i == 2)
        {
            pressure1 = 1.0;
            temperature1 = 800.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 1.0;
            temperature2 = 288.15;
            mole_concentration[0] = 0.00184;
            M5(0, pressure2);


        }
        if (i == 3)
        {
            pressure1 = 1.5;
            temperature1 = 1300.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 0.9;
            temperature2 = 298.15;
            mole_concentration[0] = 0.03226;
            M5(2, pressure2);


        }
        if (i == 4)
        {
            pressure1 = 1.0;
            temperature1 = 1300.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 1.0;
            temperature2 = 298.15;
            mole_concentration[0] = 0.03226;
            M5(3, pressure2);


        }
        if (i == 5)
        {
            pressure1 = 2.0;
            temperature1 = 1300.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 0.8;
            temperature2 = 285.2;
            mole_concentration[0] = 0.00959;
            M5(1, pressure2);


        }
        if (i == 6)
        {
            pressure1 = 1.0;
            temperature1 = 1300.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 1.0;
            temperature2 = 298.15;
            mole_concentration[0] = 0.003226;
            M5(0, pressure2);


        }
        if (i == 7)
        {
            pressure1 = 1.0;
            temperature1 = 400.0;
            mole_concentration11 = 0.11;
            mole_concentration12 = 0.11;


            pressure2 = 1.0;
            temperature2 = 300.0;
            mole_concentration[0] = 0.0323;
            M5(3, pressure2);


        }
        if (i == 8)
        {
            pressure1 = 1.0;
            temperature1 = 1600.0;
            mole_concentration11 = 0.11;
            mole_concentration12 = 0.11;
            pressure2 = 1.0;



            temperature2 = 300.0;
            mole_concentration[0] = 0.0323;
            M5(2, pressure2);


        }
        if (i == 9)
        {
            pressure1 = 2.0;
            temperature1 = 1600.0;
            mole_concentration11 = 0.11;
            mole_concentration12 = 0.11;


            pressure2 = 0.9;
            temperature2 = 300.0;
            mole_concentration[0] = 0.0323;
            M5(3, pressure2);


        }
        if (i == 10)
        {
            pressure1 = 1.0;
            temperature1 = 1900.0;
            mole_concentration11 = 0.12;
            mole_concentration12 = 0.12;


            pressure2 = 1.0;
            temperature2 = 300.0;
            mole_concentration[0] = 0.0323;
            M5(2, pressure2);


        }
        if (i == 11)
        {
            pressure1 = 1.5;
            temperature1 = 1900.0;
            mole_concentration11 = 0.12;
            mole_concentration12 = 0.12;


            pressure2 = 0.9;
            temperature2 = 300.0;
            mole_concentration[0] = 0.0323;
            M5(0, pressure2);


        }
        if (i == 12)
        {
            pressure1 = 1.5;
            temperature1 = 1900.0;
            mole_concentration11 = 0.12;
            mole_concentration12 = 0.12;


            pressure2 = 0.9;
            temperature2 = 294.2;
            mole_concentration[0] = 0.0184;
            M5(1, pressure2);


        }
        if (i == 13)
        {
            pressure1 = 1.6;
            temperature1 = 1800.0;
            mole_concentration11 = 0.14;
            mole_concentration12 = 0.12;


            pressure2 = 0.8;
            temperature2 = 300.0;
            mole_concentration[0] = 0.03;
            M5(1, pressure2);


        }
        if (i == 14)
        {
            pressure1 = 1.0;
            temperature1 = 1050.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 1.0;
            temperature2 = 294.2;
            mole_concentration[0] = 0.0184;
            M5(3, pressure2);


        }
        if (i == 15)
        {
            pressure1 = 2.0;
            temperature1 = 1050.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 0.9;
            temperature2 = 294.2;
            mole_concentration[0] = 0.0184;
            M5(0, pressure2);


        }
        if (i == 16)
        {
            pressure1 = 2.0;
            temperature1 = 1050.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 1.0;
            temperature2 = 294.2;
            mole_concentration[0] = 0.0184;
            M5(1, pressure2);


        }
        if (i == 17)
        {
            pressure1 = 2.0;
            temperature1 = 1050.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 0.7;
            temperature2 = 279.2;
            mole_concentration[0] = 0.00595;
            M5(0, pressure2);


        }
        if (i == 18)
        {
            pressure1 = 0.42;
            temperature1 = 1500.0;
            mole_concentration11 = 0.08;
            mole_concentration12 = 0.08;


            pressure2 = 0.42;
            temperature2 = 254.7;
            mole_concentration[0] = 0.00102;
            M5(2, pressure2);


        }
        if (i == 19)
        {
            pressure1 = 0.177;
            temperature1 = 1800.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 0.177;
            temperature2 = 215.8;
            mole_concentration[0] = 8.01e-6;
            M5(1, pressure2);


        }
        if (i == 20)
        {
            temperaturew = 900.0;

            pressure1 = 1.0;
            temperature1 = 1500.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 1.0;
            temperature2 = 550.0;
            mole_concentration21 = 0.05;
            mole_concentration22 = 0.05;


            pressure3 = 0.9;
            temperature3 = 300.0;
            mole_concentration[0] = 0.0323;
            M5(1, pressure3);


        }
        if (i == 21)
        {
            temperaturew = 450.0;

            pressure1 = 1.0;
            temperature1 = 1500.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 1.0;
            temperature2 = 550.0;
            mole_concentration21 = 0.05;
            mole_concentration22 = 0.05;


            pressure3 = 0.9;
            temperature3 = 300.0;
            mole_concentration[0] = 0.0323;
            M5(1, pressure3);


        }
        if (i == 22)
        {
            temperaturew = 450.0;

            pressure1 = 2.0;
            temperature1 = 1500.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 1.0;
            temperature2 = 550.0;
            mole_concentration21 = 0.05;
            mole_concentration22 = 0.05;


            pressure3 = 0.9;
            temperature3 = 300.0;
            mole_concentration[0] = 0.0323;
            M5(3, pressure3);


        }
        if (i == 23)
        {
            temperaturew = 650.0;

            pressure1 = 2.5;
            temperature1 = 1700.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 1.0;
            temperature2 = 550.0;
            mole_concentration21 = 0.05;
            mole_concentration22 = 0.05;


            pressure3 = 0.62;
            temperature3 = 273.2;
            mole_concentration[0] = 0.0038;
            M5(0, pressure3);


        }
        if (i == 24)
        {
            temperaturew = 650.0;

            pressure1 = 1.5;
            temperature1 = 1700.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 1.0;
            temperature2 = 550.0;
            mole_concentration21 = 0.05;
            mole_concentration22 = 0.05;


            pressure3 = 0.62;
            temperature3 = 273.2;
            mole_concentration[0] = 0.0038;
            M5(2, pressure3);


        }
        if (i == 25)
        {
            temperaturew = 650.0;

            pressure1 = 2.5;
            temperature1 = 1700.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 1.0;
            temperature2 = 550.0;
            mole_concentration21 = 0.05;
            mole_concentration22 = 0.05;


            pressure3 = 0.62;
            temperature3 = 273.2;
            mole_concentration[0] = 0.0038;
            M5(0, pressure3);


        }
        if (i == 26)
        {
            temperaturew = 450.0;

            pressure1 = 2.5;
            temperature1 = 1400.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 1.0;
            temperature2 = 550.0;
            mole_concentration21 = 0.05;
            mole_concentration22 = 0.05;


            pressure3 = 0.62;
            temperature3 = 273.2;
            mole_concentration[0] = 0.0038;
            M5(1, pressure3);


        }
        if (i == 27)
        {
            temperaturew = 500.0;

            pressure1 = 0.72;
            temperature1 = 1800.0;
            mole_concentration11 = 0.12;
            mole_concentration12 = 0.12;


            pressure2 = 0.48;
            temperature2 = 550.0;
            mole_concentration21 = 0.1;
            mole_concentration22 = 0.1;


            pressure3 = 0.32;
            temperature3 = 241.7;
            mole_concentration[0] = 4.13e-4;
            M5(2, pressure3);


        }
        if (i == 28)
        {
            pressure1 = 1.0;
            temperature1 = 800.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 1.0;
            temperature2 = 308.15;
            mole_concentration[0] = 0.05724;
            M5(2, pressure2);



            density_aerosol = 6.0e4;
        }
        if (i == 29)
        {
            pressure1 = 1.0;
            temperature1 = 600.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 1.0;
            temperature2 = 288.15;
            mole_concentration[0] = 0.00184;
            M5(3, pressure2);



            density_aerosol = 6.0e3;
        }
        if (i == 30)
        {
            pressure1 = 1.0;
            temperature1 = 800.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 1.0;
            temperature2 = 288.15;
            mole_concentration[0] = 0.00184;
            M5(0, pressure2);



            density_aerosol = 4.0e4;
        }
        if (i == 31)
        {
            pressure1 = 1.0;
            temperature1 = 1300.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 1.0;
            temperature2 = 298.15;
            mole_concentration[0] = 0.03226;
            M5(1, pressure2);



            density_aerosol = 0.3;
        }
        if (i == 32)
        {
            pressure1 = 1.0;
            temperature1 = 1300.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 1.0;
            temperature2 = 298.15;
            mole_concentration[0] = 0.003226;
            M5(2, pressure2);



            density_aerosol = 0.3;
        }
        if (i == 33)
        {
            pressure1 = 1.0;
            temperature1 = 1300.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 1.0;
            temperature2 = 298.15;
            mole_concentration[0] = 0.03226;
            M5(2, pressure2);



            density_aerosol = 3.0;
        }
        if (i == 34)
        {
            pressure1 = 1.0;
            temperature1 = 400.0;
            mole_concentration11 = 0.11;
            mole_concentration12 = 0.11;


            pressure2 = 1.0;
            temperature2 = 300.0;
            mole_concentration[0] = 0.0323;
            M5(0, pressure2);



            density_aerosol = 0.2;
        }
        if (i == 35)
        {
            pressure1 = 2.0;
            temperature1 = 1600.0;
            mole_concentration11 = 0.11;
            mole_concentration12 = 0.14;


            pressure2 = 0.9;
            temperature2 = 300.0;
            mole_concentration[0] = 0.0323;
            M5(1, pressure2);



            density_aerosol = 0.1;
        }
        if (i == 36)
        {
            pressure1 = 1.6;
            temperature1 = 1800.0;
            mole_concentration11 = 0.14;
            mole_concentration12 = 0.12;


            pressure2 = 0.8;
            temperature2 = 300.0;
            mole_concentration[0] = 0.024;
            M5(1, pressure2);



            density_aerosol = 0.4;
        }
        if (i == 37)
        {
            pressure1 = 1.0;
            temperature1 = 1050.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 1.0;
            temperature2 = 294.2;
            mole_concentration[0] = 0.0184;
            M5(0, pressure2);



            density_aerosol = 70.0;
        }
        if (i == 38)
        {
            pressure1 = 2.0;
            temperature1 = 1050.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 0.9;
            temperature2 = 294.2;
            mole_concentration[0] = 0.0184;
            M5(3, pressure2);



            density_aerosol = 70.0;
        }
        if (i == 39)
        {
            pressure1 = 0.42;
            temperature1 = 1500.0;
            mole_concentration11 = 0.08;
            mole_concentration12 = 0.08;


            pressure2 = 0.42;
            temperature2 = 254.7;
            mole_concentration[0] = 0.00102;
            M5(3, pressure2);



            density_aerosol = 10.0;
        }
        if (i == 40)
        {
            pressure1 = 0.177;
            temperature1 = 1800.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 0.177;
            temperature2 = 215.8;
            mole_concentration[0] = 8.01e-6;
            M5(0, pressure2);



            density_aerosol = 4.0;
        }
        if (i == 41)
        {
            temperaturew = 900.0;

            pressure1 = 1.0;
            temperature1 = 1500.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 1.0;
            temperature2 = 550.0;
            mole_concentration21 = 0.05;
            mole_concentration22 = 0.05;


            pressure3 = 0.9;
            temperature3 = 300.0;
            mole_concentration[0] = 0.0323;
            M5(0, pressure3);



            density_aerosol = 200.0;
        }
        if (i == 42)
        {
            pressure1 = 0.8;
            temperature1 = 650.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 0.8;
            temperature2 = 283.15;
            mole_concentration[0] = 0.005;
            M5(1, pressure2);


        }
        if (i == 43)
        {

            pressure1 = 1.0;
            temperature1 = 750.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 0.8;
            temperature2 = 288.15;
            mole_concentration[0] = 0.01;
            M5(2, pressure2);


        }
        if (i == 44)
        {

            pressure1 = 0.5;
            temperature1 = 900.0;
            mole_concentration11 = 0.12;
            mole_concentration12 = 0.1;


            pressure2 = 0.5;
            temperature2 = 263.15;
            mole_concentration[0] = 0.002;
            M5(3, pressure2);


        }
        if (i == 45)
        {
            pressure1 = 0.5;
            temperature1 = 500.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 1.0;
            temperature2 = 293.15;
            mole_concentration[0] = 0.015;
            M5(0, pressure2);


        }
        if (i == 46)
        {
            pressure1 = 1.0;
            temperature1 = 550.0;
            mole_concentration11 = 0.12;
            mole_concentration12 = 0.1;


            pressure2 = 0.6;
            temperature2 = 273.15;
            mole_concentration[0] = 0.004;
            M5(0, pressure2);


        }
        if (i == 47)
        {
            pressure1 = 1.0;
            temperature1 = 500.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.12;


            pressure2 = 1.0;
            temperature2 = 300.0;
            mole_concentration[0] = 0.012;
            M5(1, pressure2);


        }
        if (i == 48)
        {
            pressure1 = 2.5;
            temperature1 = 1600.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.1;


            pressure2 = 0.8;
            temperature2 = 288.15;
            mole_concentration[0] = 0.015;
            M5(3, pressure2);


        }
        if (i == 49)
        {
            pressure1 = 0.5;
            temperature1 = 1500.0;
            mole_concentration11 = 0.13;
            mole_concentration12 = 0.1;


            pressure2 = 0.9;
            temperature2 = 293.15;
            mole_concentration[0] = 0.02;
            M5(1, pressure2);


        }
        if (i == 50)
        {
            temperaturew = 1100.0;

            pressure1 = 1.0;
            temperature1 = 1400.0;
            mole_concentration11 = 0.12;
            mole_concentration12 = 0.12;


            pressure2 = 1.0;
            temperature2 = 600.0;
            mole_concentration21 = 0.06;
            mole_concentration22 = 0.07;


            pressure3 = 0.6;
            temperature3 = 260.0;
            mole_concentration[0] = 0.0015;
            M5(3, pressure3);


        }
        if (i == 51)
        {
            temperaturew = 1100.0;

            pressure1 = 0.6;
            temperature1 = 1700.0;
            mole_concentration11 = 0.1;
            mole_concentration12 = 0.12;


            pressure2 = 0.6;
            temperature2 = 700.0;
            mole_concentration21 = 0.06;
            mole_concentration22 = 0.05;


            pressure3 = 1.0;
            temperature3 = 300.0;
            mole_concentration[0] = 0.03;
            M5(3, pressure3);


        }
        if (i == 52)
        {
            pressure1 = 2.0;
            temperature1 = 1600.0;
            mole_concentration11 = 0.11;
            mole_concentration12 = 0.14;


            pressure2 = 0.9;
            temperature2 = 300.0;
            mole_concentration[0] = 0.0323;
            M5(0, pressure2);


        }
        if (i == 53)
        {
            temperaturew = 1000.0;

            pressure1 = 3.0;
            temperature1 = 1800.0;
            mole_concentration11 = 0.12;
            mole_concentration12 = 0.12;


            pressure2 = 1.5;
            temperature2 = 1400.0;
            mole_concentration21 = 0.1;
            mole_concentration22 = 0.1;


            pressure3 = 1.0;
            temperature3 = 300.0;
            mole_concentration[0] = 0.03;
            M5(3, pressure3);



            density_aerosol = 150.0;
        }
        if (i == 54)
        {
            temperaturew = 450.0;

            pressure1 = 0.8;
            temperature1 = 600.0;
            mole_concentration11 = 0.05;
            mole_concentration12 = 0.05;


            pressure2 = 0.8;
            temperature2 = 500.0;
            mole_concentration21 = 0.02;
            mole_concentration22 = 0.02;


            pressure3 = 1.0;
            temperature3 = 300.0;
            mole_concentration[0] = 0.02;
            M5(1, pressure3);



            density_aerosol = 2.0;
        }
        if (i == 55)
        {
            temperaturew = 450.0;

            pressure1 = 0.8;
            temperature1 = 350.0;
            mole_concentration11 = 0.05;
            mole_concentration12 = 0.05;


            pressure2 = 0.8;
            temperature2 = 400.0;
            mole_concentration21 = 0.02;
            mole_concentration22 = 0.02;


            pressure3 = 1.0;
            temperature3 = 288.0;
            mole_concentration[0] = 0.01;
            M5(0, pressure3);



            density_aerosol = 2.0e4;
        }
        NMn[i] = NM;



        if (i == 53)
            n = 4;
        else if (i == 54)
            n = 2;
        else if (i == 55)
            n = 1;
        else if (i < 28 || i>41)
            n = 0;
        else if (i < 31)
            n = 1;
        else if (i < 34)
            n = 2;
        else if (i < 37)
            n = 3;
        else
            n = 4;




        if ((i >= 20 && i <= 27) || i == 41 || i == 50 || i == 51 || i > 52)
        {
            fp1 = fopen("kk.txt", "rb");

            wavenumber_interval = 0.05 * pow(pressure1, 0.7);
            if (pressure1 < 0.177)
                wavenumber_interval = wavenumber_interval / 4.0;
            else if (pressure1 < 2.828)
                wavenumber_interval = wavenumber_interval / 2.0;
            cutoff_wavenumber = 7.5 * pow(pressure1, 0.65);
            gka0(mole_concentration11, mole_concentration12, pressure1, temperaturew, 0, density_aerosol, cutoff_wavenumber, wavenumber_interval, fp1, 0);
            gka0(mole_concentration11, mole_concentration12, pressure1, temperature1, 0, density_aerosol, cutoff_wavenumber, wavenumber_interval, fp1, 0);

            wavenumber_interval = 0.05 * pow(pressure2, 0.7);
            if (pressure2 < 0.177)
                wavenumber_interval = wavenumber_interval / 4.0;
            else if (pressure2 < 2.828)
                wavenumber_interval = wavenumber_interval / 2.0;
            cutoff_wavenumber = 7.5 * pow(pressure2, 0.65);
            gka0(mole_concentration21, mole_concentration22, pressure2, temperature2, 0, density_aerosol, cutoff_wavenumber, wavenumber_interval, fp1, 0);

            wavenumber_interval = 0.05 * pow(pressure3, 0.7);
            if (pressure3 < 0.177)
                wavenumber_interval = wavenumber_interval / 4.0;
            else if (pressure3 < 2.828)
                wavenumber_interval = wavenumber_interval / 2.0;
            cutoff_wavenumber = 7.5 * pow(pressure3, 0.65);
            gka0(mole_concentration[0], mole_concentration[1], pressure3, temperature3, n, density_aerosol, cutoff_wavenumber, wavenumber_interval, fp1, 1);

            fclose(fp1);
        }
        else
        {
            fp1 = fopen("kk.txt", "rb");

            wavenumber_interval = 0.05 * pow(pressure1, 0.7);
            if (pressure1 < 0.177)
                wavenumber_interval = wavenumber_interval / 4.0;
            else if (pressure1 < 2.828)
                wavenumber_interval = wavenumber_interval / 2.0;
            cutoff_wavenumber = 7.5 * pow(pressure1, 0.65);
            gka0(mole_concentration11, mole_concentration12, pressure1, temperature1, 0, density_aerosol, cutoff_wavenumber, wavenumber_interval, fp1, 0);

            wavenumber_interval = 0.05 * pow(pressure2, 0.7);
            if (pressure2 < 0.177)
                wavenumber_interval = wavenumber_interval / 4.0;
            else if (pressure2 < 2.828)
                wavenumber_interval = wavenumber_interval / 2.0;
            cutoff_wavenumber = 7.5 * pow(pressure2, 0.65);
            gka0(mole_concentration[0], mole_concentration[1], pressure2, temperature2, n, density_aerosol, cutoff_wavenumber, wavenumber_interval, fp1, 1);

            fclose(fp1);
        }



        fp = fopen("逐线计算结果-I.txt", "r");
        do
            c = fgetc(fp);
        while (c != '=');
        do
            c = fgetc(fp);
        while (c != '=');
        fscanf(fp, "%d", &k);
        for (j = 0; j < output_point_number; j++)
        {
            fscanf(fp, "%lf%lf", &b0, &output_parameter2[i][j]);
            if (b0 < LH0)
                startN[i][0] = j + 1;
            startN[i][1] = j;
        }
        fclose(fp);
    }
    chdir(c0);

	printf("pp4\n");
}

int test2(int m, int ns, double T1, double T2, double w1, double w2)
{
    int i, j, a[GN2 + 1], b;
    char c;
    char c1[300];
    char c2[300];

	sprintf(c1, ".\\groups\\1\\%.3f-%.3f-%.3f-%.3f\\group1-", T1, T2, w1, w2);
	sprintf(c2, ".\\groups\\2\\%.3f-%.3f-%.3f-%.3f\\group2-", T1, T2, w1, w2);
    
    char c0[300], cc[300];
    FILE* fp;

    fp = fopen("000.txt", "w");
    fprintf(fp, "%d%d%d%d.txt", m / 1000, (m % 1000) / 100, (m % 100) / 10, m % 10);
    fclose(fp);
    fp = fopen("000.txt", "r");
    fscanf(fp, "%s", cc);
    fclose(fp);
    if (ns == 0)
    {
        strcpy(c0, c1);
        strcat(c0, cc);
        fp = fopen(c0, "r");
        for (i = interval_number0; i < interval_number1; i++)
            fscanf(fp, "%d", &interval_KG2[i][0]);
        fclose(fp);
        b = 1;
    }
    else
    {
        strcpy(c0, c2);
        strcat(c0, cc);
        fp = fopen(c0, "r");
        for (j = 0; j <= GN2; j++)
            a[j] = 0;
        for (i = interval_number0; i < interval_number1; i++)
        {
            fscanf(fp, "%d", &interval_KG2[i][1]);
            a[interval_KG2[i][1]]++;
        }
        fclose(fp);
        for (b = 1, j = 0; j < GN2; j++)
            if (a[j] == 0)
                b = 0;
    }
    return(b);
}

int calc_err_and_write(const char* err_folder_name, int specie_type, int specie_num, double T1, double T2, double w1, double w2)
{
    int i, j, k, n;
    double output_parameter1[output_point_number][1 + GN1 + GN2],
        temperaturew,
        temperature0,
        Lh,
        pressure1,
        temperature1,
        mole_concentration11,
        mole_concentration12,
        mole_concentration13,
        Lc,
        pressure2,
        temperature2,
        mole_concentration21,
        mole_concentration22,
        mole_concentration23,
        dL,
        pressure3,
        temperature3,
        mole_concentration31,
        mole_concentration32,
        mole_concentration33,
        density_aerosol;



    int n1 = specie_type == 1 ? 0 : 1;
    int tnp = specie_type == 1 ? tnp1 : tnp2;

    for (int Tp_idx = 0; Tp_idx < tnp; ++Tp_idx) {
        FILE* fp_error;
        char error_file_name[256];
        snprintf(error_file_name, sizeof(error_file_name) - 1, "%s\\%d\\%d\\error%d-%d", err_folder_name, specie_type, Tp_idx, specie_type, specie_num);

        // 如果已经有现成存好的结果，直接读取
        if (fp_error = fopen(error_file_name, "rb")) {
            if (n1 == 0) { // H2O
                for (j = 0; j < GN1; j++)
                    for (i = 0; i < tng1; i++)
                        for (n = 0; n < testN; n++)
                            fread(error1[j][i][Tp_idx][n], sizeof(int), output_point_number, fp_error);
            }
            else if (n1 == 1) { // CO2
                for (j = 0; j < GN2; j++)
                    for (i = 0; i < tng2; i++)
                        for (n = 0; n < testN; n++)
                            fread(error2[j][i][Tp_idx][n], sizeof(int), output_point_number, fp_error);
            }
            else {
                printf("Error: specie_type = %d\n", specie_type);
                exit(-1);
            }
            fclose(fp_error);
        }
        // 如果没有现成存好的结果，计算结果并存储
        else {
            int NGL, NTP;

            int b = test2(specie_num, n1, T1, T2, w1, w2);
            if (b > 0)
            {
                int pit0, pit;
                if (n1 == 0)
                    pit0 = tng1 * testN;
                else
                    pit0 = tng2 * testN;
#pragma omp parallel for schedule(dynamic) private(i,j,n,NGL,NTP,temperature0,output_parameter1,temperaturew,Lh,pressure1,temperature1,mole_concentration11,mole_concentration12,Lc,pressure2,temperature2,mole_concentration21,mole_concentration22,dL,pressure3,temperature3,mole_concentration31,mole_concentration32,density_aerosol)                                          	
                for (pit = 0; pit < pit0; pit++)
                {
                    if (n1 == 0)
                    {
                        NGL = dng * (pit / (testN)) + ngls1;
                        NTP = Tp_idx;
                        temperature0 = TP[NTP];
                    }
                    else
                    {
                        NGL = dng * (pit / (testN)) + ngls2;
                        NTP = Tp_idx;
                        temperature0 = TP[tnp1 + NTP];
                    }
                    i = pit % testN;

                    // 很多个算例：不同波段时，这里的dL有可能不一样
                    {
                        if (i == 0)
                        {
                            Lh = 50.0;
                            pressure1 = 1.0;
                            temperature1 = 800.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            pressure2 = 1.0;
                            temperature2 = 308.15;
                            mole_concentration21 = 0.05724;

                            mole_concentration22 = 4.1e-4;
                            dL = 10000.0;
                        }
                        if (i == 1)
                        {
                            Lh = 50.0;
                            pressure1 = 1.0;
                            temperature1 = 600.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            pressure2 = 1.0;
                            temperature2 = 288.15;
                            mole_concentration21 = 0.00184;

                            mole_concentration22 = 4.1e-4;
                            dL = 10000.0;
                        }
                        if (i == 2)
                        {
                            Lh = 50.0;
                            pressure1 = 1.0;
                            temperature1 = 800.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            pressure2 = 1.0;
                            temperature2 = 288.15;
                            mole_concentration21 = 0.00184;

                            mole_concentration22 = 4.1e-4;
                            dL = 10000.0;
                        }
                        if (i == 3)
                        {
                            Lh = 80.0;
                            pressure1 = 1.5;
                            temperature1 = 1300.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            pressure2 = 0.9;
                            temperature2 = 298.15;
                            mole_concentration21 = 0.03226;

                            mole_concentration22 = 4.1e-4;
                            dL = 25000.0;
                        }
                        if (i == 4)
                        {
                            Lh = 80.0;
                            pressure1 = 1.0;
                            temperature1 = 1300.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            pressure2 = 1.0;
                            temperature2 = 298.15;
                            mole_concentration21 = 0.03226;

                            mole_concentration22 = 4.1e-4;
                            dL = 25000.0;
                        }
                        if (i == 5)
                        {
                            Lh = 80.0;
                            pressure1 = 2.0;
                            temperature1 = 1300.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            pressure2 = 0.8;
                            temperature2 = 285.2;
                            mole_concentration21 = 0.00959;

                            mole_concentration22 = 4.1e-4;
                            dL = 25000.0;
                        }
                        if (i == 6)
                        {
                            Lh = 80.0;
                            pressure1 = 1.0;
                            temperature1 = 1300.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            pressure2 = 1.0;
                            temperature2 = 298.15;
                            mole_concentration21 = 0.003226;

                            mole_concentration22 = 4.1e-4;
                            dL = 25000.0;
                        }
                        if (i == 7)
                        {
                            Lh = 80.0;
                            pressure1 = 1.0;
                            temperature1 = 400.0;
                            mole_concentration11 = 0.11;
                            mole_concentration12 = 0.11;


                            pressure2 = 1.0;
                            temperature2 = 300.0;
                            mole_concentration21 = 0.0323;

                            mole_concentration22 = 4.1e-4;
                            dL = 10000.0;
                        }
                        if (i == 8)
                        {
                            Lh = 80.0;
                            pressure1 = 1.0;
                            temperature1 = 1600.0;
                            mole_concentration11 = 0.11;
                            mole_concentration12 = 0.11;

                            pressure2 = 1.0;

                            temperature2 = 300.0;
                            mole_concentration21 = 0.0323;

                            mole_concentration22 = 4.1e-4;
                            dL = 25000.0;
                        }
                        if (i == 9)
                        {
                            Lh = 80.0;
                            pressure1 = 2.0;
                            temperature1 = 1600.0;
                            mole_concentration11 = 0.11;
                            mole_concentration12 = 0.11;


                            pressure2 = 0.9;
                            temperature2 = 300.0;
                            mole_concentration21 = 0.0323;

                            mole_concentration22 = 4.1e-4;
                            dL = 25000.0;
                        }
                        if (i == 10)
                        {
                            Lh = 80.0;
                            pressure1 = 1.0;
                            temperature1 = 1900.0;
                            mole_concentration11 = 0.12;
                            mole_concentration12 = 0.12;


                            pressure2 = 1.0;
                            temperature2 = 300.0;
                            mole_concentration21 = 0.0323;

                            mole_concentration22 = 4.1e-4;
                            dL = 25000.0;
                        }
                        if (i == 11)
                        {
                            Lh = 80.0;
                            pressure1 = 1.5;
                            temperature1 = 1900.0;
                            mole_concentration11 = 0.12;
                            mole_concentration12 = 0.12;


                            pressure2 = 0.9;
                            temperature2 = 300.0;
                            mole_concentration21 = 0.0323;

                            mole_concentration22 = 4.1e-4;
                            dL = 25000.0;
                        }
                        if (i == 12)
                        {
                            Lh = 80.0;
                            pressure1 = 1.5;
                            temperature1 = 1900.0;
                            mole_concentration11 = 0.12;
                            mole_concentration12 = 0.12;


                            pressure2 = 0.9;
                            temperature2 = 294.2;
                            mole_concentration21 = 0.0184;

                            mole_concentration22 = 4.1e-4;
                            dL = 25000.0;
                        }
                        if (i == 13)
                        {
                            Lh = 70.0;
                            pressure1 = 1.6;
                            temperature1 = 1800.0;
                            mole_concentration11 = 0.14;
                            mole_concentration12 = 0.12;


                            pressure2 = 0.8;
                            temperature2 = 300.0;
                            mole_concentration21 = 0.03;

                            mole_concentration22 = 4.1e-4;
                            dL = 25000.0;
                        }
                        if (i == 14)
                        {
                            Lh = 60.0;
                            pressure1 = 1.0;
                            temperature1 = 1050.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            pressure2 = 1.0;
                            temperature2 = 294.2;
                            mole_concentration21 = 0.0184;

                            mole_concentration22 = 4.1e-4;
                            dL = 20000.0;
                        }
                        if (i == 15)
                        {
                            Lh = 60.0;
                            pressure1 = 2.0;
                            temperature1 = 1050.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            pressure2 = 0.9;
                            temperature2 = 294.2;
                            mole_concentration21 = 0.0184;

                            mole_concentration22 = 4.1e-4;
                            dL = 20000.0;
                        }
                        if (i == 16)
                        {
                            Lh = 60.0;
                            pressure1 = 2.0;
                            temperature1 = 1050.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            pressure2 = 1.0;
                            temperature2 = 294.2;
                            mole_concentration21 = 0.0184;

                            mole_concentration22 = 4.1e-4;
                            dL = 15000.0;
                        }
                        if (i == 17)
                        {
                            Lh = 60.0;
                            pressure1 = 2.0;
                            temperature1 = 1050.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            pressure2 = 0.7;
                            temperature2 = 279.2;
                            mole_concentration21 = 0.00595;

                            mole_concentration22 = 4.1e-4;
                            dL = 20000.0;
                        }
                        if (i == 18)
                        {
                            Lh = 100.0;
                            pressure1 = 0.42;
                            temperature1 = 1500.0;
                            mole_concentration11 = 0.08;
                            mole_concentration12 = 0.08;


                            pressure2 = 0.42;
                            temperature2 = 254.7;
                            mole_concentration21 = 0.00102;
                            mole_concentration22 = 4.1e-4;

                            dL = 50000.0;
                        }
                        if (i == 19)
                        {
                            Lh = 100.0;
                            pressure1 = 0.177;
                            temperature1 = 1800.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            pressure2 = 0.177;
                            temperature2 = 215.8;
                            mole_concentration21 = 8.01e-6;

                            mole_concentration22 = 4.1e-4;
                            dL = 50000.0;
                        }
                        if (i == 20)
                        {
                            temperaturew = 900.0;

                            Lh = 50.0;
                            pressure1 = 1.0;
                            temperature1 = 1500.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            Lc = 150.0;
                            pressure2 = 1.0;
                            temperature2 = 550.0;
                            mole_concentration21 = 0.05;
                            mole_concentration22 = 0.05;


                            dL = 25000.0;
                            pressure3 = 0.9;
                            temperature3 = 300.0;
                            mole_concentration31 = 0.0323;

                            mole_concentration32 = 4.1e-4;
                        }
                        if (i == 21)
                        {
                            temperaturew = 450.0;

                            Lh = 50.0;
                            pressure1 = 1.0;
                            temperature1 = 1500.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            Lc = 150.0;
                            pressure2 = 1.0;
                            temperature2 = 550.0;
                            mole_concentration21 = 0.05;
                            mole_concentration22 = 0.05;


                            dL = 25000.0;
                            pressure3 = 0.9;
                            temperature3 = 300.0;
                            mole_concentration31 = 0.0323;

                            mole_concentration32 = 4.1e-4;
                        }
                        if (i == 22)
                        {
                            temperaturew = 450.0;

                            Lh = 50.0;
                            pressure1 = 2.0;
                            temperature1 = 1500.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            Lc = 150.0;
                            pressure2 = 1.0;
                            temperature2 = 550.0;
                            mole_concentration21 = 0.05;
                            mole_concentration22 = 0.05;


                            dL = 25000.0;
                            pressure3 = 0.9;
                            temperature3 = 300.0;
                            mole_concentration31 = 0.0323;

                            mole_concentration32 = 4.1e-4;
                        }
                        if (i == 23)
                        {
                            temperaturew = 650.0;

                            Lh = 50.0;
                            pressure1 = 2.5;
                            temperature1 = 1700.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            Lc = 300.0;
                            pressure2 = 1.0;
                            temperature2 = 550.0;
                            mole_concentration21 = 0.05;
                            mole_concentration22 = 0.05;


                            dL = 50000.0;
                            pressure3 = 0.62;
                            temperature3 = 273.2;
                            mole_concentration31 = 0.0038;

                            mole_concentration32 = 4.1e-4;
                        }
                        if (i == 24)
                        {
                            temperaturew = 650.0;

                            Lh = 50.0;
                            pressure1 = 1.5;
                            temperature1 = 1700.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            Lc = 300.0;
                            pressure2 = 1.0;
                            temperature2 = 550.0;
                            mole_concentration21 = 0.05;
                            mole_concentration22 = 0.05;


                            dL = 50000.0;
                            pressure3 = 0.62;
                            temperature3 = 273.2;
                            mole_concentration31 = 0.0038;

                            mole_concentration32 = 4.1e-4;
                        }
                        if (i == 25)
                        {
                            temperaturew = 650.0;

                            Lh = 50.0;
                            pressure1 = 2.5;
                            temperature1 = 1700.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            Lc = 100.0;
                            pressure2 = 1.0;
                            temperature2 = 550.0;
                            mole_concentration21 = 0.05;
                            mole_concentration22 = 0.05;


                            dL = 50000.0;
                            pressure3 = 0.62;
                            temperature3 = 273.2;
                            mole_concentration31 = 0.0038;

                            mole_concentration32 = 4.1e-4;
                        }
                        if (i == 26)
                        {
                            temperaturew = 450.0;

                            Lh = 50.0;
                            pressure1 = 2.5;
                            temperature1 = 1400.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            Lc = 300.0;
                            pressure2 = 1.0;
                            temperature2 = 550.0;
                            mole_concentration21 = 0.05;
                            mole_concentration22 = 0.05;


                            dL = 50000.0;
                            pressure3 = 0.62;
                            temperature3 = 273.2;
                            mole_concentration31 = 0.0038;

                            mole_concentration32 = 4.1e-4;
                        }
                        if (i == 27)
                        {
                            temperaturew = 500.0;

                            Lh = 50.0;
                            pressure1 = 0.72;
                            temperature1 = 1800.0;
                            mole_concentration11 = 0.12;
                            mole_concentration12 = 0.12;


                            Lc = 150.0;
                            pressure2 = 0.48;
                            temperature2 = 550.0;
                            mole_concentration21 = 0.1;
                            mole_concentration22 = 0.1;


                            dL = 50000.0;
                            pressure3 = 0.32;
                            temperature3 = 241.7;
                            mole_concentration31 = 4.13e-4;
                            mole_concentration32 = 4.1e-4;

                        }
                        if (i == 28)
                        {
                            Lh = 50.0;
                            pressure1 = 1.0;
                            temperature1 = 800.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            pressure2 = 1.0;
                            temperature2 = 308.15;
                            mole_concentration21 = 0.05724;

                            mole_concentration22 = 4.1e-4;
                            dL = 10000.0;

                            density_aerosol = 6.0e4;
                        }
                        if (i == 29)
                        {
                            Lh = 50.0;
                            pressure1 = 1.0;
                            temperature1 = 600.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            pressure2 = 1.0;
                            temperature2 = 288.15;
                            mole_concentration21 = 0.00184;

                            mole_concentration22 = 4.1e-4;
                            dL = 10000.0;

                            density_aerosol = 6.0e3;
                        }
                        if (i == 30)
                        {
                            Lh = 50.0;
                            pressure1 = 1.0;
                            temperature1 = 800.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            pressure2 = 1.0;
                            temperature2 = 288.15;
                            mole_concentration21 = 0.00184;

                            mole_concentration22 = 4.1e-4;
                            dL = 10000.0;

                            density_aerosol = 4.0e4;
                        }
                        if (i == 31)
                        {
                            Lh = 80.0;
                            pressure1 = 1.0;
                            temperature1 = 1300.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            pressure2 = 1.0;
                            temperature2 = 298.15;
                            mole_concentration21 = 0.03226;

                            mole_concentration22 = 4.1e-4;
                            dL = 25000.0;

                            density_aerosol = 0.3;
                        }
                        if (i == 32)
                        {
                            Lh = 80.0;
                            pressure1 = 1.0;
                            temperature1 = 1300.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            pressure2 = 1.0;
                            temperature2 = 298.15;
                            mole_concentration21 = 0.003226;

                            mole_concentration22 = 4.1e-4;
                            dL = 25000.0;

                            density_aerosol = 0.3;
                        }
                        if (i == 33)
                        {
                            Lh = 80.0;
                            pressure1 = 1.0;
                            temperature1 = 1300.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            pressure2 = 1.0;
                            temperature2 = 298.15;
                            mole_concentration21 = 0.03226;

                            mole_concentration22 = 4.1e-4;
                            dL = 7500.0;

                            density_aerosol = 3.0;
                        }
                        if (i == 34)
                        {
                            Lh = 80.0;
                            pressure1 = 1.0;
                            temperature1 = 400.0;
                            mole_concentration11 = 0.11;
                            mole_concentration12 = 0.11;


                            pressure2 = 1.0;
                            temperature2 = 300.0;
                            mole_concentration21 = 0.0323;

                            mole_concentration22 = 4.1e-4;
                            dL = 10000.0;

                            density_aerosol = 0.2;
                        }
                        if (i == 35)
                        {
                            Lh = 80.0;
                            pressure1 = 2.0;
                            temperature1 = 1600.0;
                            mole_concentration11 = 0.11;
                            mole_concentration12 = 0.14;


                            pressure2 = 0.9;
                            temperature2 = 300.0;
                            mole_concentration21 = 0.0323;

                            mole_concentration22 = 4.1e-4;
                            dL = 15000.0;

                            density_aerosol = 0.1;
                        }
                        if (i == 36)
                        {
                            Lh = 2.0;
                            pressure1 = 1.6;
                            temperature1 = 1800.0;
                            mole_concentration11 = 0.14;
                            mole_concentration12 = 0.12;


                            pressure2 = 0.8;
                            temperature2 = 300.0;
                            mole_concentration21 = 0.024;

                            mole_concentration22 = 4.1e-4;
                            dL = 5000.0;

                            density_aerosol = 0.4;
                        }
                        if (i == 37)
                        {
                            Lh = 60.0;
                            pressure1 = 1.0;
                            temperature1 = 1050.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            pressure2 = 1.0;
                            temperature2 = 294.2;
                            mole_concentration21 = 0.0184;

                            mole_concentration22 = 4.1e-4;
                            dL = 20000.0;

                            density_aerosol = 70.0;
                        }
                        if (i == 38)
                        {
                            Lh = 60.0;
                            pressure1 = 2.0;
                            temperature1 = 1050.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            pressure2 = 0.9;
                            temperature2 = 294.2;
                            mole_concentration21 = 0.0184;

                            mole_concentration22 = 4.1e-4;
                            dL = 15000.0;

                            density_aerosol = 70.0;
                        }
                        if (i == 39)
                        {
                            Lh = 100.0;
                            pressure1 = 0.42;
                            temperature1 = 1500.0;
                            mole_concentration11 = 0.08;
                            mole_concentration12 = 0.08;


                            pressure2 = 0.42;
                            temperature2 = 254.7;
                            mole_concentration21 = 0.00102;
                            mole_concentration22 = 4.1e-4;

                            dL = 50000.0;

                            density_aerosol = 10.0;
                        }
                        if (i == 40)
                        {
                            Lh = 100.0;
                            pressure1 = 0.177;
                            temperature1 = 1800.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            pressure2 = 0.177;
                            temperature2 = 215.8;
                            mole_concentration21 = 8.01e-6;

                            mole_concentration22 = 4.1e-4;
                            dL = 50000.0;

                            density_aerosol = 4.0;
                        }
                        if (i == 41)
                        {
                            temperaturew = 900.0;

                            Lh = 50.0;
                            pressure1 = 1.0;
                            temperature1 = 1500.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            Lc = 150.0;
                            pressure2 = 1.0;
                            temperature2 = 550.0;
                            mole_concentration21 = 0.05;
                            mole_concentration22 = 0.05;


                            dL = 25000.0;
                            pressure3 = 0.9;
                            temperature3 = 300.0;
                            mole_concentration31 = 0.0323;

                            mole_concentration32 = 4.1e-4;

                            density_aerosol = 200.0;
                        }
                        if (i == 42)
                        {
                            Lh = 15.0;
                            pressure1 = 0.8;
                            temperature1 = 650.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            pressure2 = 0.8;
                            temperature2 = 283.15;
                            mole_concentration21 = 0.005;

                            mole_concentration22 = 4.1e-4;
                            dL = 25000.0;
                        }
                        if (i == 43)
                        {
                            Lh = 5.0;
                            pressure1 = 1.0;
                            temperature1 = 750.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            pressure2 = 0.8;
                            temperature2 = 288.15;
                            mole_concentration21 = 0.01;

                            mole_concentration22 = 4.1e-4;
                            dL = 25000.0;
                        }
                        if (i == 44)
                        {
                            Lh = 5.0;
                            pressure1 = 0.5;
                            temperature1 = 900.0;
                            mole_concentration11 = 0.12;
                            mole_concentration12 = 0.1;


                            pressure2 = 0.5;
                            temperature2 = 263.15;
                            mole_concentration21 = 0.002;

                            mole_concentration22 = 4.1e-4;
                            dL = 30000.0;
                        }
                        if (i == 45)
                        {
                            Lh = 10.0;
                            pressure1 = 0.5;
                            temperature1 = 500.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            pressure2 = 1.0;
                            temperature2 = 293.15;
                            mole_concentration21 = 0.015;

                            mole_concentration22 = 4.1e-4;
                            dL = 10000.0;
                        }
                        if (i == 46)
                        {
                            Lh = 10.0;
                            pressure1 = 1.0;
                            temperature1 = 550.0;
                            mole_concentration11 = 0.12;
                            mole_concentration12 = 0.1;


                            pressure2 = 0.6;
                            temperature2 = 273.15;
                            mole_concentration21 = 0.004;

                            mole_concentration22 = 4.1e-4;
                            dL = 25000.0;
                        }
                        if (i == 47)
                        {
                            Lh = 10.0;
                            pressure1 = 1.0;
                            temperature1 = 500.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.12;


                            pressure2 = 1.0;
                            temperature2 = 300.0;
                            mole_concentration21 = 0.012;

                            mole_concentration22 = 4.1e-4;
                            dL = 20000.0;
                        }
                        if (i == 48)
                        {
                            Lh = 150.0;
                            pressure1 = 2.5;
                            temperature1 = 1600.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.1;


                            pressure2 = 0.8;
                            temperature2 = 288.15;
                            mole_concentration21 = 0.015;

                            mole_concentration22 = 4.1e-4;
                            dL = 25000.0;
                        }
                        if (i == 49)
                        {
                            Lh = 10.0;
                            pressure1 = 0.5;
                            temperature1 = 1500.0;
                            mole_concentration11 = 0.13;
                            mole_concentration12 = 0.1;


                            pressure2 = 0.9;
                            temperature2 = 293.15;
                            mole_concentration21 = 0.02;

                            mole_concentration22 = 4.1e-4;
                            dL = 20000.0;
                        }
                        if (i == 50)
                        {
                            temperaturew = 1100.0;

                            Lh = 20.0;
                            pressure1 = 1.0;
                            temperature1 = 1400.0;
                            mole_concentration11 = 0.12;
                            mole_concentration12 = 0.12;


                            Lc = 20.0;

                            pressure2 = 1.0;
                            temperature2 = 600.0;
                            mole_concentration21 = 0.06;
                            mole_concentration22 = 0.07;


                            dL = 30000.0;
                            pressure3 = 0.6;
                            temperature3 = 260.0;
                            mole_concentration31 = 0.0015;
                            mole_concentration32 = 4.1e-4;

                        }
                        if (i == 51)
                        {
                            temperaturew = 1100.0;

                            Lh = 50.0;
                            pressure1 = 0.6;
                            temperature1 = 1700.0;
                            mole_concentration11 = 0.1;
                            mole_concentration12 = 0.12;


                            Lc = 150.0;
                            pressure2 = 0.6;
                            temperature2 = 700.0;
                            mole_concentration21 = 0.06;
                            mole_concentration22 = 0.05;


                            dL = 10000.0;
                            pressure3 = 1.0;
                            temperature3 = 300.0;
                            mole_concentration31 = 0.03;
                            mole_concentration32 = 4.1e-4;

                        }
                        if (i == 52)
                        {
                            Lh = 80.0;
                            pressure1 = 2.0;
                            temperature1 = 1600.0;
                            mole_concentration11 = 0.11;
                            mole_concentration12 = 0.14;


                            pressure2 = 0.9;
                            temperature2 = 300.0;
                            mole_concentration21 = 0.0323;

                            mole_concentration22 = 4.1e-4;
                            dL = 25000.0;
                        }
                        if (i == 53)
                        {
                            temperaturew = 1000.0;

                            Lh = 20.0;
                            pressure1 = 3.0;
                            temperature1 = 1800.0;
                            mole_concentration11 = 0.12;
                            mole_concentration12 = 0.12;


                            Lc = 40.0;
                            pressure2 = 1.5;
                            temperature2 = 1400.0;
                            mole_concentration21 = 0.1;
                            mole_concentration22 = 0.1;


                            dL = 10000.0;
                            pressure3 = 1.0;
                            temperature3 = 300.0;
                            mole_concentration31 = 0.03;

                            mole_concentration32 = 4.1e-4;

                            density_aerosol = 150.0;
                        }
                        if (i == 54)
                        {
                            temperaturew = 450.0;

                            Lh = 3.0;
                            pressure1 = 0.8;
                            temperature1 = 600.0;
                            mole_concentration11 = 0.05;
                            mole_concentration12 = 0.05;


                            Lc = 5.0;
                            pressure2 = 0.8;
                            temperature2 = 500.0;
                            mole_concentration21 = 0.02;
                            mole_concentration22 = 0.02;


                            dL = 7500.0;
                            pressure3 = 1.0;
                            temperature3 = 300.0;
                            mole_concentration31 = 0.02;

                            mole_concentration32 = 4.1e-4;

                            density_aerosol = 2.0;
                        }
                        if (i == 55)
                        {
                            temperaturew = 450.0;

                            Lh = 3.0;
                            pressure1 = 0.8;
                            temperature1 = 350.0;
                            mole_concentration11 = 0.05;
                            mole_concentration12 = 0.05;


                            Lc = 5.0;
                            pressure2 = 0.8;
                            temperature2 = 400.0;
                            mole_concentration21 = 0.02;
                            mole_concentration22 = 0.02;



                            dL = 10000.0;
                            pressure3 = 1.0;
                            temperature3 = 288.0;
                            mole_concentration31 = 0.01;

                            mole_concentration32 = 4.1e-4;

                            density_aerosol = 2.0e4;
                        }
                    }

                    if ((i >= 20 && i <= 27) || i == 41 || i == 50 || i == 51 || i > 52)
                        msmgwb4(NMn[i], n1, NGL, temperature0, output_parameter1, temperaturew,
                            Lh, pressure1, temperature1, mole_concentration11, mole_concentration12,
                            Lc, pressure2, temperature2, mole_concentration21, mole_concentration22,
                            dL, pressure3, temperature3, mole_concentration31, mole_concentration32);
                    else
                        msmgwb2(NMn[i], n1, NGL, temperature0, output_parameter1,
                            Lh, pressure1, temperature1, mole_concentration11, mole_concentration12,
                            dL, pressure2, temperature2, mole_concentration21, mole_concentration22);

                    if (n1 == 0) {
                        for (int n = 0; n < GN1; n++)
                            for (int j = 0; j < output_point_number; j++)
                                error1[n][pit / (testN)][NTP][i][j] = (int)(10000.0 * output_parameter1[j][1 + n] / output_parameter2[i][j]);
                    }

                    if (n1 == 1) {

                        for (int n = 0; n < GN2; n++)
                            for (int j = 0; j < output_point_number; j++)
                                error2[n][pit / (testN)][NTP][i][j] = (int)(10000.0 * (output_parameter1[j][1 + n] - output_parameter2[i][j] / (GN2 + 0.0)) / output_parameter2[i][j]);
                    }

                    printf("\r%d-%d-%d ", NTP, NGL, i);
                }


                if (n1 == 0) {
                    fp_error = fopen(error_file_name, "wb");
                    for (int j = 0; j < GN1; j++)
                        for (int i = 0; i < tng1; i++)
                            for (int n = 0; n < testN; n++)
                                fwrite(error1[j][i][Tp_idx][n], sizeof(int), output_point_number, fp_error);
                    fclose(fp_error);
                }
                else if (n1 == 1) {
                    fp_error = fopen(error_file_name, "wb");
                    for (int j = 0; j < GN2; j++)
                        for (int i = 0; i < tng2; i++)
                            for (int n = 0; n < testN; n++)
                                fwrite(error2[j][i][Tp_idx][n], sizeof(int), output_point_number, fp_error);
                    fclose(fp_error);
                }
                else {
                    printf("ERROR: n1 = %d\n", n1);
                    exit(-1);
                }
            }

            else {
                printf("ERROR: b = %d\n", b);
                exit(-1);
            }

        }
    }
    printf("\n");

    return 0;
}



void main(int argc, char* argv[])
{

    //----------------------------------------------
    string mission_name = argv[1]; // 任务结果存放文件夹编号
    int step_num = atoi(argv[2]);  // 进行到第几步，从0开始
    int pop_size_small = atoi(argv[3]); // 粗算小种群规模
    int pop_size_large = atoi(argv[4]); // 精算大种群规模
    int total_num_processes = atoi(argv[5]);  // 总共分成几个程序算
    int present_process_num = atoi(argv[6]); // 当前程序编号，从0开始

    double T1 = atof(argv[7]);
    double T2 = atof(argv[8]);
    double w1 = atof(argv[9]);
    double w2 = atof(argv[10]);

    char method_name[256];
    snprintf(method_name, sizeof(method_name) - 1, "NSGA2-%s", mission_name.c_str());
    printf("method_name = %s, parallel_threads_number = %d, step_num = %d\n", method_name, parallel_threads_number, step_num);
    printf("H2O_num_groups = %d, CO2_num_groups = %d, pop_size_small = %d, pop_size_large = %d\n", M1, M2, pop_size_small, pop_size_large);
    printf("total_num_processes = %d, present_process_num = %d\n", total_num_processes, present_process_num);
    //----------------------------------------------


    double b5[GN1 * tng1 + GN2 * tng2][2], b6, b7;
    double bb[3][401];
    int i, j, k, m, m1, n, n1, n2, n3, NGL, NTP, CN[1 + 2 * GN1 + 2 * GN2];
    int nn[2], nn1[2];
    float b, b0, b1, b2, b3, b4;
    long long m2, n4;
    FILE* fp;
    char cc[300];
    int pit, pit0;

    printf("p1\n");
    
    for (j = 0; j < GN1 * tng1 + GN2 * tng2; j++)
    {
        b5[j][0] = 1000.0;
        b5[j][1] = 1000.0;
    }

    omp_set_num_threads(parallel_threads_number);
    test0();

    
	printf("p2\n");
    // error文件存放路径
    char error_data_dir[32] = ".\\error_data";
    mkdir(error_data_dir);

    // H2O
    char error_data_dir1[32] = ".\\error_data\\1";
    mkdir(error_data_dir1);
    for (int i = 0; i < tnp1; ++i) { // 不同Tp
        char error_data_dir1_Tp[32];
        snprintf(error_data_dir1_Tp, sizeof(error_data_dir1_Tp) - 1, "%s\\%d", error_data_dir1, i);
        mkdir(error_data_dir1_Tp);
    }

    // CO2
    char error_data_dir2[32] = ".\\error_data\\2";
    mkdir(error_data_dir2);
    for (int i = 0; i < tnp2; ++i) { // 不同Tp
        char error_data_dir2_Tp[32];
        snprintf(error_data_dir2_Tp, sizeof(error_data_dir2_Tp) - 1, "%s\\%d", error_data_dir2, i);
        mkdir(error_data_dir2_Tp);
    }

    // 计算结果文件存放路径
    char folder_name[256];
    snprintf(folder_name, sizeof(folder_name) - 1, ".\\%s", method_name);
    mkdir(folder_name);

    // 临时数据文件存放路径
    char temp_dir[256];
    snprintf(temp_dir, sizeof(temp_dir) - 1, "%s\\temp", folder_name);
    mkdir(temp_dir);
    int total_row_num = M1;
    int total_col_num = M2;

    // 扫描流程记录文件存放路径
    mkdir("scan_flow_records");
    char scan_flow_record_dir[256];
    snprintf(scan_flow_record_dir, sizeof(scan_flow_record_dir) - 1, "scan_flow_records\\scan_flow_record-%.3f-%.3f-%.3f-%.3f.txt", T1, T2, w1, w2);


    // 先用小种群数扫描
    Goal err_now;


    if (step_num < 9) {
        int pop_size = pop_size_small;

        // 计算出row的编号，写入文件
        if (step_num % 4 == 0) {
            vec1i ordered_rows_index(total_row_num, 0);

            // step0的话直接随机生成
            if (step_num == 0) {
                random_device rd;
                default_random_engine eng(rd());
                uniform_int_distribution<int> distr_int_row(0, total_row_num - 1);
                ordered_rows_index[0] = distr_int_row(eng);
            }
            else {
                vec1d err_rows(total_row_num, -1e5);
                for (int i = 0; i < total_num_processes; ++i) {
                    char res_dir[256];
                    snprintf(res_dir, sizeof(res_dir) - 1, "%s\\%d\\%d", temp_dir, step_num - 1, i);
                    vec1d temp(total_row_num);
                    read_from_file2(res_dir, temp);

                    int beg_group_num = int(total_row_num * (i + 0.0) / total_num_processes);
                    int end_group_num = int(total_row_num * (i + 1.0) / total_num_processes);
                    for (int it = beg_group_num; it < end_group_num; ++it)
                        err_rows[it] = temp[it];
                }
                writeVec2File(scan_flow_record_dir, 1, err_rows);
                ordered_rows_index = sort_indexes(err_rows);
            }
            char out_dir[256];
            snprintf(out_dir, sizeof(out_dir) - 1, "%s\\%d", temp_dir, step_num);
            write_to_file2(out_dir, ordered_rows_index);

            writeNum2File(scan_flow_record_dir, 1, ordered_rows_index[0]);
        }

        // 分布式计算err_cols，写入文件夹内多个文件
        if (step_num % 4 == 1) {

            // 读取上一步给出的scan_row_num
            vec1i ordered_rows_index(total_row_num);
            char last_dir[256];
            snprintf(last_dir, sizeof(last_dir) - 1, "%s\\%d", temp_dir, step_num - 1);
            read_from_file2(last_dir, ordered_rows_index);
            int scan_row_num = ordered_rows_index[0];


            // 计算err_cols
            vec1d err_cols(total_col_num, -1e5);
            int beg_group_num = int(total_col_num * (present_process_num + 0.0) / total_num_processes);
            int end_group_num = int(total_col_num * (present_process_num + 1.0) / total_num_processes);

            calc_err_and_write(error_data_dir, 1, scan_row_num, T1, T2, w1, w2);
            for (int it = beg_group_num; it < end_group_num; ++it) {
                calc_err_and_write(error_data_dir, 2, it, T1, T2, w1, w2);
                NSGA2_one(scan_row_num, it, pop_size, 10000, method_name, startN, err_now);
                err_cols[it] = err_now.error;
            }

            // 写入文件夹
            char out_folder_dir[256];
            snprintf(out_folder_dir, sizeof(out_folder_dir) - 1, "%s\\%d", temp_dir, step_num);
            mkdir(out_folder_dir);
            char out_dir[256];
            snprintf(out_dir, sizeof(out_dir) - 1, "%s\\%d", out_folder_dir, present_process_num);
            write_to_file2(out_dir, err_cols);
        }

        // 计算出col的编号，写入文件
        if (step_num % 4 == 2) {
            vec1i ordered_cols_index(total_col_num);
            vec1d err_cols(total_col_num, -1e5);
            for (int i = 0; i < total_num_processes; ++i) {
                char res_dir[256];
                snprintf(res_dir, sizeof(res_dir) - 1, "%s\\%d\\%d", temp_dir, step_num - 1, i);
                vec1d temp(total_col_num);
                read_from_file2(res_dir, temp);

                int beg_group_num = int(total_col_num * (i + 0.0) / total_num_processes);
                int end_group_num = int(total_col_num * (i + 1.0) / total_num_processes);
                for (int it = beg_group_num; it < end_group_num; ++it)
                    err_cols[it] = temp[it];
            }
            writeVec2File(scan_flow_record_dir, 1, err_cols);
            ordered_cols_index = sort_indexes(err_cols);

            char out_dir[256];
            snprintf(out_dir, sizeof(out_dir) - 1, "%s\\%d", temp_dir, step_num);
            write_to_file2(out_dir, ordered_cols_index);

            writeNum2File(scan_flow_record_dir, 1, ordered_cols_index[0]);
        }

        // 分布式计算err_rows，写入文件夹内多个文件
        if (step_num % 4 == 3) {

            // 读取上一步给出的scan_row_num

            vec1i ordered_cols_index(total_col_num);
            char last_dir[256];
            snprintf(last_dir, sizeof(last_dir) - 1, "%s\\%d", temp_dir, step_num - 1);
            read_from_file2(last_dir, ordered_cols_index);
            int scan_col_num = ordered_cols_index[0];

            // 计算err_rows
            vec1d err_rows(total_row_num, -1e5);
            int beg_group_num = int(total_row_num * (present_process_num + 0.0) / total_num_processes);
            int end_group_num = int(total_row_num * (present_process_num + 1.0) / total_num_processes);


            calc_err_and_write(error_data_dir, 2, scan_col_num, T1, T2, w1, w2);
            for (int it = beg_group_num; it < end_group_num; ++it) {
                calc_err_and_write(error_data_dir, 1, it, T1, T2, w1, w2);
                NSGA2_one(it, scan_col_num, pop_size, 10000, method_name, startN, err_now);
                err_rows[it] = err_now.error;
            }

            // 写入文件夹
            char out_folder_dir[256];
            snprintf(out_folder_dir, sizeof(out_folder_dir) - 1, "%s\\%d", temp_dir, step_num);
            mkdir(out_folder_dir);
            char out_dir[256];
            snprintf(out_dir, sizeof(out_dir) - 1, "%s\\%d", out_folder_dir, present_process_num);
            write_to_file2(out_dir, err_rows);
        }
    }
    else {

        int pop_size = pop_size_large; // 大种群数精算

        char ordered_rows_index_dir[256];
        snprintf(ordered_rows_index_dir, sizeof(ordered_rows_index_dir) - 1, "%s\\8", temp_dir);
        vec1i ordered_rows_index(total_row_num);
        read_from_file2(ordered_rows_index_dir, ordered_rows_index);

        char ordered_cols_index_dir[256];
        snprintf(ordered_cols_index_dir, sizeof(ordered_cols_index_dir) - 1, "%s\\6", temp_dir);
        vec1i ordered_cols_index(total_col_num);
        read_from_file2(ordered_cols_index_dir, ordered_cols_index);

        vec1i chosen_rows_index = vec1i(ordered_rows_index.begin(), ordered_rows_index.begin() + chosen_rows_num);
        vec1i chosen_cols_index = vec1i(ordered_cols_index.begin(), ordered_cols_index.begin() + chosen_cols_num);

        int total_num_cases = chosen_cols_num * chosen_rows_num;

        // 分布式精算
        if (step_num == 9) {



            vec1d err_final(total_num_cases);
            int beg_group_num = int(total_num_cases * (present_process_num + 0.0) / total_num_processes);
            int end_group_num = int(total_num_cases * (present_process_num + 1.0) / total_num_processes);

            // 先cols变，后rows变
            for (int i = beg_group_num; i < end_group_num; ++i) {
                int row_num = ordered_rows_index[i / chosen_cols_num];
                int col_num = ordered_cols_index[i % chosen_cols_num];
                calc_err_and_write(error_data_dir, 1, row_num, T1, T2, w1, w2);
                calc_err_and_write(error_data_dir, 2, col_num, T1, T2, w1, w2);
                NSGA2_one(row_num, col_num, pop_size, 10000, method_name, startN, err_now);
                err_final[i] = err_now.error;
            }

            // 写入文件夹
            char out_folder_dir[256];
            snprintf(out_folder_dir, sizeof(out_folder_dir) - 1, "%s\\%d", temp_dir, step_num);
            mkdir(out_folder_dir);
            char out_dir[256];
            snprintf(out_dir, sizeof(out_dir) - 1, "%s\\%d", out_folder_dir, present_process_num);
            write_to_file2(out_dir, err_final);
        }

        // 汇总精算结果
        if (step_num == 10) {

            writeVec2File(scan_flow_record_dir, 1, chosen_rows_index);
            writeVec2File(scan_flow_record_dir, 1, chosen_cols_index);

            vec1d err_final(total_num_cases);

            for (int i = 0; i < total_num_processes; ++i) {
                char res_dir[256];
                snprintf(res_dir, sizeof(res_dir) - 1, "%s\\%d\\%d", temp_dir, step_num - 1, i);
                vec1d temp(total_num_cases);
                read_from_file2(res_dir, temp);

                int beg_group_num = int(total_num_cases * (i + 0.0) / total_num_processes);
                int end_group_num = int(total_num_cases * (i + 1.0) / total_num_processes);
                for (int it = beg_group_num; it < end_group_num; ++it)
                    err_final[it] = temp[it];
            }

            // 写精算的最小err
            vec1d err_final_ordered(err_final);
            sort(err_final_ordered.begin(), err_final_ordered.end());
            writeNum2File(scan_flow_record_dir, 1, err_final_ordered[0]);

            for (int i = 0; i < chosen_rows_num; ++i) {
                vec1d temp = vec1d(err_final.begin() + i * chosen_cols_num, err_final.begin() + (i + 1) * chosen_cols_num);
                writeVec2File(scan_flow_record_dir, 1, temp);
            }

            char out_dir[256];
            snprintf(out_dir, sizeof(out_dir) - 1, "%s\\%d", temp_dir, step_num);
            write_to_file2(out_dir, err_final_ordered);
        }
    }







}
