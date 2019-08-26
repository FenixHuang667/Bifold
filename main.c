//
//  main.c
//  Inv_seq
//
//  Created by Fenix Huang on 8/19/15.
//  Copyright (c) 2015 Fenix Huang. All rights reserved.
//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "head.h"
#include "utils.h"
#include "decompose.h"
#include "stat.h"
#include "pfunc.h"



char Inputfilename[500];
char **Infile;
char Outputfilename[500];
char *Outfile;
int sample, heat, mul;
FILE *Input;
FILE *Output;

float c1 = 0.5, c2 = 0.5;


void usage()
{
    printf("usage:\n"
           "Smapler [-OPTIONS}\n"
           "[-i inputfilename] Specify the input file (default input.in)"
           "[-o outputfilename] Specify the output file (default output.out)"
           "[-m number] Specify the sample size (defaout:1000)"
           "[-c] c1 coefficient"
           "[-v] c2 coefficient"
           "[-showS] show the sampled sequences"
           "[-showF] show the folded structures"
           "[-showE] show the energy"
           "[-showMI show the mutual information"
    );
    exit(0);
}

int main(int argc, const char * argv[]) {
    int i,j, k, r,distance, cnt = 0, length, DP_interation;
    char buff[1000];
    int *ptable, *seq_pat, *str_pat, *ref_seq, pair[1000], p_randA[500], p_randB[500];
    char structureA[1000], structureB[1000], seqA[1000], seqB[1000], name[500], randseq[500], *struc_randA, *struc_randB;
    patten *seq_list=NULL, *index;
    Bi_block *Block_index, *current_index;
    char *fold_seq;
    double **PF;
    Bi_block *head;
    int max_index, seq_ana = 0, mul = 1000;
    char struc_list[500][500], *prefix;
    int *s1, *s2;
    int sequence_length_scale = 1;
    
    //seq_ana flag
    //0: normal bifold
    
    srand48(time(NULL)); //random seed
    //strcpy(Inputfilename, "input.in");
    //strcpy(Outputfilename, "output.out");
    
    //strcpy(Inputfilename, "/Users/biocomplexity/Documents/Program/Bifold/Bifold/Input/rand_ensemble/rand_ensemble6.in");
    strcpy(Inputfilename, "/Users/biocomplexity/Documents/Program/Bifold/Riboswitch/Input/mgtE_Mg.in");
    //strcpy(Inputfilename, "/Users/biocomplexity/Documents/Program/Bifold/Bifold/Input/rand_structure/rand_struc83.in");
    //strcpy(Outputfilename, "/Users/biocomplexity/Documents/Program/Bifold/Output/BSUBT_yitI_SAM.out");
    
    strcpy(Outputfilename, "/Users/biocomplexity/Documents/Program/Bifold/Bifold/mgtE_Mg_on.out");
    //strcpy(Outputfilename, "/Users/biocomplexity/Documents/Program/Bifold/Bifold/Output/rand_ensemble/rand_ensemble1_off.out");
    //strcpy(Outputfilename, "/Users/biocomplexity/Documents/Program/Bifold/Bifold/Output/rand_ensemble/rand_ensemble6_off.out");
    //strcpy(Outputfilename, "./test.out");
    
    
    for (i=1; i<argc; i++) {
        if (argv[i][0]=='-') {
            switch (argv[i][1]) {
                case 'i':
                    if (i==argc-1) usage();
                    Infile = argv[++i];
                    strcpy(Inputfilename, Infile);
                    break;
                case 'm':
                    if (i==argc-1) usage();
                    mul = atoi(argv[++i]);
                    break;
                case 'c':
                    if (i==argc-1) usage();
                    c1 = atoi(argv[++i]);
                    break;
                case 'v':
                    if (i==argc-1) usage();
                    c2 = atoi(argv[++i]);
                    break;
                case 'o':
                    if (i==argc-1) usage();
                    Outfile = argv[++i];
                    strcpy(Outputfilename, Outfile);
                    break;
                default: usage();
            }
        }
    }
    
    Input=fopen(Inputfilename,"r");
    if (Input==NULL) {
        printf("Input file error!\n");
        exit(0);
    }
    
    Output=fopen(Outputfilename,"w");
    if (Output==NULL) {
        printf("Output file error!\n");
        exit(0);
    }
    
    if (seq_ana == 1) {
        seq_analysis(Input);
    } else if (seq_ana == 0){
        while(!feof(Input)) {
            fscanf(Input, "%s\n", name);
            fscanf(Input, "%s", structureA);
            fscanf(Input, "%s", structureB);
            fscanf(Input, "%s", seqA);
        }
    
        printf("%s\n", structureA);
        printf("%s\n", structureB);
        printf("%s\n", seqA);
        printf("Length: %d\n", strlen(structureA));
        
        if (sequence_length_scale) {
            c1 *=((float)  1)/((float) strlen(structureA));
            c2 *=((float)  1)/((float) strlen(structureA));
        }
        //scale energy by 1/n, n is the sequence length
        
        
    // printf("%s\n", seqB);
    
    /*
    seq_list = bicompatible_generator(structureA, structureB, 0);
    index = seq_list;
    while (index!=NULL) {
        printf("%s\n", index->pat);
        index = index->next;
    }
    */
    
        length = strlen(structureA);
        fold_seq = (char *) calloc (length+5, sizeof(char));
        for (i=0; i<length; i++) fold_seq[i] = '_';
        fold_seq[length] = 0;
    
        head = sequential_removal(structureA, structureB);
        max_index = est_max_index(head);
        printf("Max index:%d\n", max_index);

    
       // if (max_index < 20) {
       //     mfe_double_seq(head, fold_seq, structureA, structureB);
            // printf("%s\n", fold_seq);
       // } else {
       //     printf("Max index exceed limit!\n");
       // }
    
    
        if (max_index < 20) {
            
            PF = PF_double_seq(head, structureA, structureB);
            seq_list = Biseq_sampler(head, PF, mul, structureA, structureB);
            
            length = 0;
            Block_index = head;
            while (Block_index !=NULL) {
                length ++;
                Block_index = Block_index->next;
            }
            
            fprintf(Output, "PF: %e\n", PF[length][0]); 
        
            fprintf(Output, "Sample: \n");
            index = seq_list;
            while (index!=NULL) {
                fprintf(Output, "%s\n", index->pat);
                index = index->next;
            }
        } else {
            printf("Max index exceed limit!\n");
        }
    } else if (seq_ana == 2) {
        while(!feof(Input)) {
            fscanf(Input, "%s\n", name);
            fscanf(Input, "%s", structureA);
            fscanf(Input, "%s", structureB);
            fscanf(Input, "%s", seqA);
        }
        bicompatible_generator(structureA, structureB, 0);
    } else if (seq_ana == 3) {
        //fscanf(Input, "%s\n", name);
        //fscanf(Input, "%s", structureA);
        //fscanf(Input, "%s", structureB);
        //fscanf(Input, "%s", seqA);
        for (i=0; i<1000; i++) {
            randomseq(100, randseq);
            fold_sec(randseq, p_randA);
            randomseq(100, randseq);
            fold_sec(randseq, p_randB);
        
            struc_randA = pair2structure(p_randA);
            struc_randB = pair2structure(p_randB);
        
            distance = bp_distance(struc_randA, struc_randB);
            head = sequential_removal(struc_randA, struc_randB);
            max_index = est_max_index(head);
            // compute max index and bp distance of structures A and B
        
            if (max_index<=15 && distance<=30) {
                printf("%s\n", struc_randA);
                printf("%s\n", struc_randB);
                printf("Max_index:%d   bp distance: %d\n", max_index, distance);
            }
            
        
            Block_index = head;
            while (Block_index!=NULL) {
                current_index = Block_index;
                Block_index = Block_index ->next;
                free(current_index->external);
                free(current_index->pair1);
                free(current_index->pair2);
                free(current_index->vertics);
                free(current_index);
            }
            free(struc_randA);
            free(struc_randB);
            //free space
        }
    } else if (seq_ana == 4) {
        fclose(Input);
        fclose(Output);
        
        Input = fopen("rand_ensemble.in", "r");
        prefix = "rand_ensemble";
        cnt = 0;
        while (!feof(Input)) {
            fscanf(Input, "%s\n", struc_list[cnt]);
            cnt++;
        }
        r = 1;
        for (k=1; k<=300 && r<=100; k++) {
            i = rand() % cnt;
            j = rand() % cnt;
            //printf("%d %d %d\n", cnt, i,j);
            if (i==j) continue;
            s1 = structure2pair(struc_list[i]);
            s2 = structure2pair(struc_list[j]);
            distance = bp_distance_p(s1, s2);
            head = sequential_removal(struc_list[i], struc_list[j]);
            max_index = est_max_index(head);
                
            if (distance > 20 && max_index<15 && r<=100) {
                strcpy(Outputfilename, prefix);
                sprintf(buff, "%d", r);
                strcat(Outputfilename, buff);
                strcat(Outputfilename, ".in");
                    
                Output=fopen(Outputfilename,"w");
                if (Output==NULL) {
                    printf("Output file error!\n");
                    exit(0);
                }
                fprintf(Output, "%s%d\n", prefix, r);
                fprintf(Output, "%s\n", struc_list[i]);
                fprintf(Output, "%s\n", struc_list[j]);
                randomseq(strlen(struc_list[i]), randseq);
                fprintf(Output, "%s", randseq);
                fclose (Output);
                    
                r++;
                
                free(s1);
                free(s2);
            }
        }
        
    }
    
    fclose(Input);
    fclose(Output);
    
    printf("done!\n");
    
    return 0;
}
