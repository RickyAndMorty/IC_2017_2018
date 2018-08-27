#include <stdio.h>
#include <stdlib.h>
#include <time.h>
void CalculaP(PSO* pso)
{
    int i,j;
    //float a;


    time_t t;
    srand((unsigned) time(&t));
    for(i = 0; i < pso->K; i ++)
    {
        for(j = 0; j < pso->M; j ++)
        {
            //a = 11.00F + (float)((rand() % 89));
            //pso->P[i][j] = a / 10000.00F;
            pso->P[i][j] = (((double)(rand())) / RAND_MAX) / 100.00F;
            //pso->P[i][j] = (((double)(rand())) / RAND_MAX) / 100.00F;
            //pso->P[i][j] = (((double)(rand())) / RAND_MAX) / 1000.00F;
            //pso->P[i][j] = (((double)(rand())) / RAND_MAX) / 10000.00F;
        }
    }

}
