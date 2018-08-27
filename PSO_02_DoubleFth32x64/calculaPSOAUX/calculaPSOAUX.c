#include <stdio.h>
#include <stdlib.h>

void calculaPSOAUX(PSOAUX* psoaux, PSO* pso)
{
    psoaux->K = pso->K;
    psoaux->M = pso->M;
    psoaux->SNIR = alloca2DD(pso->iteracoes,pso->K);
    psoaux->Pgbest = alloca2DD(pso->iteracoes,pso->K);
    psoaux->iteracoes = pso->iteracoes;
    psoaux->jP = alloca1DD(psoaux->iteracoes);
    psoaux->SNR = alloca1DD(psoaux->K);
}
