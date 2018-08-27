#include <stdio.h>
#include <stdlib.h>

void calculajPibest(PSO* pso)
{
    int i;
    for(i = 0; i < pso->M; i ++)
    {
        pso->jPibest[i] = pso->jP[i];
    }
}
