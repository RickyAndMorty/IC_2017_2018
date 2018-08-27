#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void calculaF(REDE* rede)
{
    int i;
    for(i = 0; i < rede->K; i ++)
    {
        rede->F[i] = rede->Rc / rede->Rb;
    }
}
