#include <stdio.h>
#include <stdlib.h>
void calculaPibest(PSO* pso)
{
    int i,j;
    for(i = 0; i < pso->K; i ++)
    {
        for(j = 0; j < pso->M; j ++)
        {
            pso->Pibest[i][j] = pso->P[i][j];
        }
    }
}
