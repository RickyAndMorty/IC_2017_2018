#include <stdio.h>
#include <stdlib.h>
#include <time.h>
typedef struct Fth// Lista dinâmica que armazena as posições onde a SNIR é maior que a SNIR target
{
    int x,y;
    struct Fth* prox;
}Fth;

typedef struct rede// Estrutura que armazena as variáveis da rede de fíbra óptica
{
    double** H;
    double* g_t, * Ltx_i, * F;
    double a_star,alfa,Gamp,Pn,Pmax,Pmin;
    double Vmax, Vmin,SNR_target, Rc, Rb, q, sigma_cc2;
    int K,M;
}REDE;
typedef struct pso//Estrutura que armazena dinãmicamente as variáveis do algoritmo PSO
{
    double** SNIR,** P, ** G, ** v, ** Pibest;
    double *Gii, * F, * jP, *Pgbest, * jPibest,*SNR;
    double Pn,Pmax, Pmin,jPgbest,c1,c2,Wadp,Vmin,Vmax;
    Fth* fth;
    int K, M,iteracoes;
}PSO;


typedef struct psoaux
{
    double Wadp,c1,c2;
    double** SNIR,** P, ** G,** Fth, ** v, ** Pibest, **Pgbest;
    double *Gii, * F, * jP, * jPibest,*SNR;
    int K,iteracoes,M;

}PSOAUX;
#include "aloca/aloca.c"
#include "InsertionSort/InsertSort.c"
#include "Randomica/randomica.c"
#include "imprimir.c"
#include "CalculaH/CalculaH.c"
#include "CalculaF/CalculaF.c"
#include "rede.c"
#include "pso.c"
#include "CalculaP/CalculaP.c"
#include "CalculaFth/CalculaFth.c"
#include "CalculaPibest/calculaPibest.c"
#include "CalculaPgbest/CalculaPgbest.c"
#include "CalculaG/CalculaG.c"
#include "calculajPibest/calculajPibest.c"
#include "CalculaVelocidade/calculaVelocidade.c"
#include "CalculaSNR/calculaSNR.c"
#include "calculaSNIR_aux/calculaSNIR_aux.c"
#include "calculaPgbest_aux/calculaPgbest_aux.c"
#include "calculaPSOAUX/calculaPSOAUX.c"
#include "gravarTxt/gravarTxt.c"



int verificaSNR(PSO* pso)
{
    int i,cont;
    cont = 0;
    for(i = 0; i < pso->K; i ++)
    {
        if(pso->SNR[i] > 21.00 || pso->SNR[i] < 20.00)
        {
            cont = 1;
        }
    }
    return cont;
}


int main()
{
    int (i);
    REDE* rede = allocaREDE();
    PSO* pso = allocaPSO();
    PSOAUX * psoaux = allocaPSOAUX();
    clock_t tempo;
	tempo = clock();
   // do
    //{
        calculaRede(rede);
        inserirPSO(pso,rede);
        calculaPSOAUX(psoaux,pso);
        CalculaP(pso);
        calculaPibest(pso);
        calculaPgbest(pso);
        calculaG(pso,rede);
        SNIR(pso);
        calculaFth(pso,rede);
        fitness(pso);

        for(i = 0; i < pso->iteracoes; i ++)
        {
            SNIR(pso);
            calculaFth(pso,rede);
            fitness(pso);
            bestLocal(pso);
            bestGlobal(pso);
            speed(pso);

            speedBounds(pso);
            populationUpdate(pso);
            powerBounds(pso);



            calculaSNIR_aux(pso,psoaux,i);
            calculaPgbest_aux(pso,psoaux,i);
            psoaux->jP[i] = pso->jPgbest;
        }
        calculaSNR(pso,rede);
        psoaux->SNR = pso->SNR;
   //}while(verificaSNR(pso));

    imprimir1DF(pso->SNR,pso->K);
	printf("Jpgbest: %f\n",pso->jPgbest);
	printf("Wadp: %f\n",pso->Wadp);
	printf("c1: %f\n",pso->c1);
	printf("c2: %f\n",pso->c2);
	psoaux->c1 = pso->c1;
	psoaux->c2 = pso->c2;
	psoaux->Wadp = pso->Wadp;
	printf("Tempo:%f\n ",(clock() - tempo) / (double)CLOCKS_PER_SEC);
    gravarTxt(psoaux);
    system("PAUSE");
    return 0;
}
