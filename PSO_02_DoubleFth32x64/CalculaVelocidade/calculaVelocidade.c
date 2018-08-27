#include <stdio.h>
#include <stdlib.h>


void zerosV(PSO* pso)
{
    int i,j;
    for(i = 0; i < pso->K; i ++)
    {
        for(j = 0; j < pso->M; j ++)
        {
            pso->v[i][j] = 0.00F;
        }
    }
}
