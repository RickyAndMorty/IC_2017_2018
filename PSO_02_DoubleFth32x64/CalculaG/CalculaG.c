#include <stdio.h>
#include <stdlib.h>

void calculaG(PSO* pso, REDE* rede)
{
    int i,j;
    for(i = 0; i < pso->K; i ++)
    {
        for(j = 0; j < pso->K; j ++)
        {
            if(i == j)
            {
                pso->G[i][j] = 1.00F;
            }
            else
            {
                pso->G[i][j] = rede->H[i][j];
            }
        }
    }
}
