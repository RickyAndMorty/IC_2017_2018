#include <stdio.h>
#include <stdlib.h>

Fth* inserirFth(Fth* fth,int i,int j)
{
    Fth* aux1 = fth;
    Fth* aux2 = NULL;
    if(fth == NULL)
    {
        fth = allocaFth();
        fth->x = i;
        fth->y = j;
	fth->prox = NULL;
        return fth;
    }
    while(aux1->prox != NULL)
    {
        aux1 = aux1->prox;
    }
    aux2 = allocaFth();
    aux2->x = i;
    aux2->y = j;
    aux1->prox = aux2;
    return fth;
}


void calculaFth(PSO* pso, REDE* rede)
{
    int i,j;
    for(i = 0; i < pso->K; i ++)
    {
        for(j = 0; j < pso->M; j ++)
        {
            if(pso->SNIR[i][j] > rede->SNR_target)
            {
                pso->fth = inserirFth(pso->fth,i,j);
            }
        }
    }
}
