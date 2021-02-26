#include <stdio.h>
#include <stdlib.h>
#include "cAdd.h"  

void cAdd(int m, int n, int blkMax, double* in1, int lda, int* jpvt, double* in2, int numel)
{
    int i;
    for (i=0; i<numel; i++) {
        in1[i] = in1[i] + in2[i];
    }
}
