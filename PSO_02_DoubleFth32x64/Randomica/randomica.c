#include <stdio.h>
#include <stdlib.h>
#include <time.h>


double** randomica(int linhas, int colunas)
{
    double** matriz = NULL;
    matriz = alloca2DD(linhas,colunas);
    int i,j;
    time_t t;
    srand((unsigned) time(&t));
    for(i = 0; i < linhas; i ++)
    {
        for(j = 0; j < colunas; j ++)
        {
            matriz[i][j] = ((double)(rand())) /RAND_MAX;
        }
    }
    return matriz;
}
