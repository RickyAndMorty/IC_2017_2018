#include <stdio.h>
#include <stdlib.h>

void calculaH(REDE* rede)
{
    int i,j;
    double** H1 = NULL;
    double** H2 = NULL;
    H1 = alloca2DD(rede->K,rede->K);
    H2 = alloca2DD(rede->K,rede->K);
    for(i = 0; i < rede->K; i ++)
    {
        for(j = 0; j < rede->K; j ++)
        {
            H1[i][j] = rede->g_t[j];
            H2[j][i] = rede->g_t[j];
        }
    }
    for(i = 0; i < rede->K; i ++)
    {
        for(j = 0; j < rede->K; j ++)
        {
            if(i == j)
            {
                rede->H[i][j] = 0.00F;
            }
            else
            {
                rede->H[i][j] = H1[i][j] / H2[i][j];
            }

        }
    }
    H1 = libera2DD(H1,rede->K);
    H2 = libera2DD(H2,rede->K);
}
