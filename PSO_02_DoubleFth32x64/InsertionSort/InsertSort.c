#include <stdio.h>
#include <stdlib.h>


void insertionSort(double* jP, double * MaxJ, int* idx, int M)
{
    int i,j;
    double chave;
    double * original = NULL;
    original = alloca1DD(M);
    for(i = 0; i < M; i ++)
    {
        MaxJ [i] = original[i] = jP[i];
    }
    for(j = 1; j < M; j ++)
    {
        chave = MaxJ[j];
        i  = j - 1;
        while((i > -1) && (MaxJ[i] > chave))
        {
            MaxJ[i + 1] = MaxJ[i];
            i --;
        }
        MaxJ[i + 1] = chave;
    }
    for(i = 0; i < M; i ++)
    {
        j = 0;
        while(MaxJ[i] != original[j])
        {
            j ++;
        }
        original[j] = -1;
        idx[i] = j;
    }
    original = libera1DD(original);
}
