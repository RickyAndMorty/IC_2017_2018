#include <stdio.h>
#include <stdlib.h>


double* alloca1DD(int colunas)
{
    double* vetor = NULL;
    vetor = ((double*)calloc(colunas, sizeof(double)));
    if(vetor == NULL)
    {
        printf("Erro.\n");
        exit(1);
    }
    return vetor;
}
double** alloca2DD(int linhas, int colunas)
{
    int i;
    double** matriz = NULL;
    matriz = ((double**)calloc(linhas , sizeof(double*)));
    if(matriz == NULL)
    {
        printf("Erro.\n");
        exit(1);
    }
    for(i = 0; i < linhas; i ++)
    {
        matriz[i] = ((double*)calloc(colunas , sizeof(double)));
        if(matriz[i] == NULL)
        {
            printf("Erro.\n");
            exit(1);
        }
    }
    return matriz;
}
double* libera1DD(double* vetor)
{
    free(vetor);
    vetor = NULL;
    return vetor;
}
double** libera2DD(double** matriz,int linhas)
{
    int i;
    for(i = 0; i < linhas; i ++)
    {
        free(matriz[i]);
        matriz[i] = NULL;
    }
    free(matriz);
    matriz = NULL;
    return matriz;
}
int * alloca1DI(int colunas)
{
    int* vetor = NULL;
    vetor = ((int*) calloc(colunas , sizeof(int)));
    if(vetor == NULL)
    {
        printf("Erro.\n");
        exit(1);
    }
    return vetor;
}
int* libera1DI(int* vetor)
{
    free(vetor);
    vetor = NULL;
    return vetor;
}

REDE* allocaREDE()
{

    REDE* rede = ((REDE*)calloc(1, sizeof(REDE)));
    if(rede == NULL)
    {
        printf("Erro.\n");
        exit(1);
    }
    return rede;
}
PSO* allocaPSO()
{
    PSO* pso = ((PSO*)calloc(1, sizeof(PSO)));
    if(pso == NULL)
    {
        printf("Erro.\n");
        exit(1);
    }
    return pso;
}
PSOAUX* allocaPSOAUX()
{
    PSOAUX* psoaux = ((PSOAUX*)calloc(1, sizeof(PSOAUX)));
    if(psoaux == NULL)
    {
        printf("Erro.\n");
        exit(1);
    }
    return psoaux;
}
Fth* allocaFth()
{
    Fth* fth = NULL;
    fth = ((Fth*)calloc(1, sizeof(Fth)));
    if(fth == NULL)
    {
        printf("Erro.\n");
        exit(1);
    }
    fth->prox = NULL;
    return fth;
}

Fth* liberaFth(Fth* fth)
{
    Fth* aux = NULL;
    while(fth != NULL)
    {
        aux = fth;
        fth = fth->prox;
        free(aux);
    }
    fth = NULL;
    return fth;
}
















