#include <stdio.h>
#include <stdlib.h>

void calculaPgbest(PSO * pso)
{
    int i;
    for(i = 0; i < pso->K; i ++)
    {
        pso->Pgbest[i] = (((double)(rand())) / RAND_MAX) / 100.00F;//pso->Pmax;
    }
}

