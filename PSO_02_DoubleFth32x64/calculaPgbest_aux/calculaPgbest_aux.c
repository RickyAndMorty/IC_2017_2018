#include <stdio.h>
#include <stdlib.h>



void calculaPgbest_aux(PSO* pso, PSOAUX* psoaux, int z)
{
    int i;
    for(i = 0; i < pso->K; i ++)
    {
        psoaux->Pgbest[z][i] = pso->Pgbest[i];
    }
}
