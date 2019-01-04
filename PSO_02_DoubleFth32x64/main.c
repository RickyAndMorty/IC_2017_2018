#include <stdio.h>
#include <stdlib.h>
#include <time.h>
typedef struct Fth// Dynamic list where the position of the SNIR is bigger than the SNIR target is stored. 
{
    int x,y;
    struct Fth* prox;
}Fth;


typedef struct rede// Structure that store the optical network variables
{
    double** H;
    double* g_t, * Ltx_i, * F;
    double a_star,alfa,Gamp,Pn,Pmax,Pmin;
    double Vmax, Vmin,SNR_target, Rc, Rb, q, sigma_cc2;
    int K,M;
}REDE;
typedef struct pso// Structure that store dynamically the PSO algorithm variables.
{
    double** SNIR,** P, ** G, ** v, ** Pibest;
    double *Gii, * F, * jP, *Pgbest, * jPibest,*SNR;
    double Pn,Pmax, Pmin,jPgbest,c1,c2,Wadp,Vmin,Vmax;
    Fth* fth;
    int K, M,iteracoes;
}PSO;
typedef struct psoaux// This structure is utilized to store the PSO algorithm variables in each round. 
{
    double Wadp,c1,c2;
    double** SNIR,** P, ** G,** Fth, ** v, ** Pibest, **Pgbest;
    double *Gii, * F, * jP, * jPibest,*SNR;
    int K,iteracoes,M;

}PSOAUX;
// This files are created separatelly, but belongs to the same folder.
#include "aloca/aloca.c"//This file contains the functions responsable for allocating the memory dynamically. 
#include "InsertionSort/InsertSort.c"// This file contains the function InsertionSort, that is, ordenation by insertion.
#include "Randomica/randomica.c"// This file constains a function responsable for generate a matrix which contains numbers randomically generated.
#include "imprimir.c"// This file contains the functions responsable for shown the results in the prompt windows.
#include "CalculaH/CalculaH.c"// This file contains the function responsable for generate the normalized matrix H of interference.
#include "CalculaF/CalculaF.c"// This file contais the function responsable for calculate the quocient of the bit rate by the chip rate.
#include "rede.c"// This file contains all the attributes related to the passive optical network, all the datas are generated mathematically.
#include "pso.c"// This file contains the functions related to the pso algorithm.
#include "CalculaP/CalculaP.c"// This file contains the function responsable for generate the surface of research of the pso algorithm. 
#include "CalculaFth/CalculaFth.c"// This file contains the function responsable for give guindace to the particle swarm of the algorithm. 
#include "CalculaPibest/calculaPibest.c"// This file contains the function responsable for, initially, generate the best local positions. 
#include "CalculaPgbest/CalculaPgbest.c"// This file contains the function responsable for, initially, generate the best global positions. 
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
