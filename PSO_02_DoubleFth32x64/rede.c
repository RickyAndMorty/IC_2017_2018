#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void calculaRede(REDE* rede)
{
    int i;
    rede->Gamp = 100.00;
    rede->alfa = 0.0461;
    rede->a_star = 0.0016;
    rede->K = 32;
    rede->M = 64;
    rede->Pmax = 0.0100;
    rede->Pmin = 1.0 * pow(10,-10);
    rede->Vmax = 0.0020;
    rede->Vmin = - 0.0020;
    rede->SNR_target = 1.5259 * pow(10,-5);
    rede->sigma_cc2 = 0.01;
    rede->q = 256.00;

    rede->Rc = 1.00;
    rede->Rb = 1.00;
    rede->Pn = 2.00 * (rede->Gamp - 1.0) * 193.0 * pow(10,12) * 100.0* pow(10,9) * 6.63 * pow(10,-34);
    rede->g_t = alloca1DD(rede->K);
    rede->Ltx_i = alloca1DD(rede->K);
    rede->F = alloca1DD(rede->K);
    rede->H = alloca2DD(rede->K,rede->K);

    rede->Ltx_i[0] = 48.8779;
    rede->Ltx_i[1] = 26.9518;
    rede->Ltx_i[2] = 40.0057;
    rede->Ltx_i[3] = 46.3947;
    rede->Ltx_i[4] = 18.7903;
    rede->Ltx_i[5] = 22.8634;
    rede->Ltx_i[6] = 37.3182;
    rede->Ltx_i[7] = 5.2198;
    rede->Ltx_i[8] = 32.6824;
    rede->Ltx_i[9] = 48.5307;
    rede->Ltx_i[10] = 10.2101;
    rede->Ltx_i[11] = 40.0123;
    rede->Ltx_i[12] = 14.8661;
    rede->Ltx_i[13] = 26.3720;
    rede->Ltx_i[14] = 46.7481;
    rede->Ltx_i[15] = 29.7142;
    rede->Ltx_i[16] = 28.8779;
    rede->Ltx_i[17] = 36.9518;
    rede->Ltx_i[18] = 30.0057;
    rede->Ltx_i[19] = 16.3947;
    rede->Ltx_i[20] = 8.7903;
    rede->Ltx_i[21] = 32.8634;
    rede->Ltx_i[22] = 27.3182;
    rede->Ltx_i[23] = 4.2198;
    rede->Ltx_i[24] = 12.6824;
    rede->Ltx_i[25] = 28.5307;
    rede->Ltx_i[26] = 40.2101;
    rede->Ltx_i[27] = 30.0123;
    rede->Ltx_i[28] = 34.8661;
    rede->Ltx_i[29] = 26.3720;
    rede->Ltx_i[30] = 26.7481;
    rede->Ltx_i[31] = 19.7142;
    for(i = 0; i < rede->K; i ++)
    {
        rede->g_t[i] = rede->a_star * pow(2.7183,-rede->alfa * rede->Ltx_i[i]) * rede->Gamp * 0.0025;
    }
    calculaH(rede);
    calculaF(rede);
}
