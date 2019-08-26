//
//  decompose.h
//  Bifold
//
//  Created by Fenix Huang on 7/11/18.
//  Copyright © 2018 Fenix Huang. All rights reserved.
//

#ifndef decompose_h
#define decompose_h

#include <stdio.h>
#include "params.h"

#define MAXENG 100000

#define R GASCONST

extern paramT *P;

typedef struct Loop Loop;

struct Loop
{
    int type; // 0 unknown, 1 hairpin, 2 interior loop/helix, 3 multiloop
    int representative; 
    int *nucleotide;
    Loop *next;
};

typedef struct Bi_block Bi_block;

struct Bi_block
{
    int *vertics; // in 1, not in 0
    int *pair1; // 0 unpaired
    int *pair2;
    int *external; // [0]: # of nucleotides
    Bi_block *next;
};

#endif /* decompose_h */

extern void chain_decomposition(int *chain, int *chain_head, int *cycle_head, int *pair1, int *pair2);
extern patten *bicompatible_generator (char *str1, char *str2, int mul);

extern int est_max_index(Bi_block *head); 

extern Bi_block *sequential_removal (char *str1, char *str2);

extern int mfe_double_seq(Bi_block *head, char *seq, char *str1, char *str2);

extern double **PF_double_seq(Bi_block *head, char *str1, char *str2);

extern patten *Biseq_sampler (Bi_block *head, double **PF, int mul, char *str1, char *str2); 
