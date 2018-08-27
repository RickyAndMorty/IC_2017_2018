#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void  calculaSNR(PSO* pso, REDE* rede)
{
    int i;
    for(i = 0; i < pso->K; i ++)
    {
        pso->SNR[i] = 10.00F * (logf(pso->SNIR[i][pso->M-1] * ((rede->q*rede->q) / rede->sigma_cc2)) / logf(10.00F));
    }

}
