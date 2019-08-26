//
//  decompose.c
//  Bifold
//
//  Created by Fenix Huang on 7/11/18.
//  Copyright © 2018 Fenix Huang. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "fold_vars.h"
#include "pair_mat.h"
#include "energy_par.h"
#include "energy_const.h"
#include "fold.h"

#include "head.h"
#include "decompose.h"
#include "utils.h"
#include "stat.h"

#define MAXINDEX 30

paramT *P; //energy parameters


long double EXP(int G)
{
    return (long double) exp(-1*(10*(G))/(R*(P->temperature+K0)));
}


Loop *new_loop (void)
{
    Loop *p;
    p = (Loop *) calloc (1, sizeof(Loop));
    p->type = 0;
    p->representative = 0;
    p->nucleotide = (int *) calloc (200, sizeof(int));
    p->next = NULL;
    return p;
}
// initiate a loop struc

Bi_block *new_Bi_block(int n)
{
    Bi_block *b;
    int i;
    b = (Bi_block *) calloc (1, sizeof(Bi_block));
    
    b->vertics = (int *) calloc (n+5, sizeof(int));
    b->vertics[0] = n;
    for (i=1; i<=n; i++) b->vertics[i] = 0;
    
    b->external = (int *) calloc (200, sizeof(int));
    b->external[0] = 0;
    for (i=1; i<50; i++) b->external[i] = 0;
    
    b->pair1 = (int *) calloc (n+5, sizeof(int));
    b->pair1[0] = n;
    for (i=1; i<=n; i++) b->pair1[i] = 0;
    
    b->pair2 = (int *) calloc (n+5, sizeof(int));
    b->pair2[0] = n;
    for (i=1; i<=n; i++) b->pair2[i] = 0;
    
    b->next = NULL;
    
    return b;
}
// initiate a bi-struc

void print_block(Bi_block *block)
{
    int i, n=block->pair1[0];
    for (i=1; i<=n; i++) printf("%d ", block->vertics[i]); printf("\n");
    for (i=1; i<=n; i++) printf("%d ", block->pair1[i]); printf("\n");
    for (i=1; i<=n; i++) printf("%d ", block->pair2[i]); printf("\n");
    for (i=1; i<=block->external[0]; i++) printf("%d ", block->external[i]); printf("\n");
}
// print block

void chain_decomposition(int *chain, int *chain_head, int *cycle_head, int *pair1, int *pair2) {
    int i, j, n;
    int num_chain = 0, num_cycle = 0;
    
    n = pair1[0];
    
    //chain decomposition
    for (i=0; i<=n; i++) chain[i] = 0;
    for (i=1; i<=n; i++) {
        if (chain[i] != 0) continue;
        if (pair1[i] == 0 && pair2[i] !=0) {
            j = i;
            num_chain++;
            chain_head[num_chain] = j;
            chain[j] = pair2[j];
            j = pair2[j];
            while (pair1[j]!=0 && pair2[j]!=0) {
                if (pair1[j]!=0) {
                    chain[j] = pair1[j];
                    j = pair1[j];
                }
                if (pair2[j]!=0) {
                    chain[j] = pair2[j];
                    j = pair2[j];
                }
            }
            chain[j] = -1;
        }
        
        if (pair1[i] != 0 && pair2[i] ==0) {
            j = i;
            num_chain++;
            chain_head[num_chain] = j;
            chain[j] = pair1[j];
            j = pair1[j];
            while (pair1[j]!=0 && pair2[j]!=0) {
                if (pair2[j]!=0) {
                    chain[j] = pair2[j];
                    j = pair2[j];
                }
                if (pair1[j]!=0) {
                    chain[j] = pair1[j];
                    j = pair1[j];
                }
            }
            chain[j] = -1;
        }
    }
    chain_head[0] = num_chain;
    
    // for (i=1; i<=n; i++) printf("%d ", chain[i]); printf("\n");
    
    for (i=1; i<=n; i++) {
        if (chain[i] != 0) continue;
        if (pair1[i] != 0 && pair2[i] !=0) {
            j=i;
            num_cycle++;
            cycle_head[num_cycle] = j;
            chain[j] = pair2[j];
            j = pair2[j];
            while (j!=cycle_head[num_cycle]) {
                if (j!=cycle_head[num_cycle]) {
                    chain[j] = pair1[j];
                    j = pair1[j];
                }
                if (j!=cycle_head[num_cycle]) {
                    chain[j] = pair2[j];
                    j = pair2[j];
                }
            }
        }
    }
    cycle_head[0] = num_cycle;
}
// decompose chain/cycle of bi-struc
// put in array chain, chain_head/cycle_head remember the representative
// chain_head[0]/cycle_head[0] is the number if cycle/chain


unsigned long long seq2index(char *str)
{
    int i, l = strlen(str), c;
    unsigned long long index = 0, base = 1;
    
    for (i = l-1; i>=0; i--) {
        c = nucleotide2code(str[i]);
        index += base * (c-1);
        base *= 4;
    }
    return index;
}

char * index2seq(unsigned long long index, int length)
{
    char *str;
    int l = 0, c;
    unsigned long long temp = index;
    
    str = (char *) calloc (length+2, sizeof(char));
    
    for (l = 0; l<length; l++) {
        c = temp % 4;
        str[length-l-1] = code2nucleotide(c+1);
        temp = (unsigned long long) ((temp - c) / 4);
    }
    str[length] = 0;
    
    return str;
}



Loop *loop_decomposition (int *ptable)
{
    int i, n = ptable[0], k, p, bp;
    Loop *head = NULL, *temp, *index;
    
    for (i=1; i<=n; i++) {
        bp = 0;
        if (i < ptable[i]) {
            temp = new_loop();
            temp->representative = i;
            k=1;
            temp->nucleotide[k++] = i;
            p = ptable[i];
            temp->nucleotide[k++] = p;
            p--;
            do {
                temp->nucleotide[k++] = p;
                if (ptable[p] !=0) {
                    p = ptable[p];
                    temp->nucleotide[k++] = p;
                    bp++;
                }
                p--;
            } while (p != i);
            
            temp->nucleotide[0] = k-1;
            if (bp == 0) temp->type = 1;
            else if (bp == 1) temp->type = 2;
            else temp->type = 3;
            
            temp ->next = head;
            head = temp;
        }
    }
    
    
    index = head;
    while (index!=NULL) {
        printf("Loop representative: %d\n", index->representative);
        printf("Loop type: %d\n", index->type);
        for (i=0; i<=index->nucleotide[0]; i++) printf("%d ", index->nucleotide[i]);
        printf("\n");
        index = index->next;
    }
    
    return head;
}
//decompose a struc into loop struc


int in_chain(int v, int *chain, int *chain_head, int *cycle_head)
{
    int i, num_chain = chain_head[0], j;
    for (i=1; i<=num_chain; i++) {
        j = chain_head[i];
        while (j!=-1) {
            if (v == j) return chain_head[i];
            j = chain[j];
        }
    }
    return 0;
}
// find the representative chain contains v



int in_cycle(int v, int *chain, int *chain_head, int *cycle_head)
{
    int i, num_cycle = cycle_head[0], j;
    for (i=1; i<=num_cycle; i++) {
        j = cycle_head[i];
        do {
            if (v == j) return cycle_head[i];
            j = chain[j];
            
        } while (j!=cycle_head[i]);
    }
    return 0;
}
// find the representative cycle contains v


int find_loop(int v, int *pair)
{
    int i, n= pair[0], p, q;
    p=0; q=n+1;
    for (i=1; i<=n;i++) {
        if (pair[i] && i < pair[i] && i<v && v<pair[i]) {
            if (i>p && pair[i]<q) {
                p = i; q = pair[i];
            }
        }
    }
    return p;
}
//Given a vertex v, find the loop contains it, presented by (p,q)

int looptype (int r, int *pair)
{
    int i, j, k=0, c;
    if (r>pair[r]) {
        i = pair[r]; j = r;
    } else {
        i = r; j = pair[r];
    }
    
    c = i+1;
    do {
        while (pair[c]==0) c++;
        if (c<j) {
            c = pair[c]+1;
            k++;
        }
    } while (c<j);
    
    return k;
}


int ** loop_v (int *pair1, int *pair2)
{
    int i, j, n=pair1[0];
    int **vertex;
    
    vertex = (int **) calloc(n+2, sizeof(int *));
    for (i=0; i<=n; i++) {
        vertex[i] = (int *) calloc (5, sizeof(int));
        for (j=0; j<=4; j++) vertex[i][j] = 0;
    }
    
    vertex[0][0] = n;

    for (i=1; i<=n; i++) {
        vertex[i][1] = find_loop(i, pair1);
        if (pair1[i]) {
            if (i<pair1[i]) vertex[i][2] = i;
            else vertex[i][2] = pair1[i];
        }
        
        vertex[i][3] = find_loop(i, pair2);
        if (pair2[i]) {
            if (i<pair2[i]) vertex[i][4] = i;
            else vertex[i][4] = pair2[i];
        }
    }
    return vertex;
}
// find loop of a vertex

int ** loop_v_truncated (int *pair1, int *pair2)
{
    int i, j, n=pair1[0], k, p, q;
    int **vertex, type;
    
    vertex = (int **) calloc(n+2, sizeof(int *));
    for (i=0; i<=n; i++) {
        vertex[i] = (int *) calloc (5, sizeof(int));
        for (j=0; j<=4; j++) vertex[i][j] = 0;
    }
    
    for (i=1; i<=n; i++) {
        if (pair1[i] == 0 || i>pair1[i]) continue;
        j = pair1[i];
        type = looptype(i, pair1);
        if (type == 0) {
            vertex[i][1] = i;
            vertex[j][1] = i;
            if (j-i-1 <= 4) {
                for (k=i+1; k<j; k++) {
                    vertex[k][1] = i;
                }
            } else {
                vertex[i+1][1] = i;
                vertex[j-1][1] = i;
            }
        } else if (type == 1) {
            vertex[i][1] = i;
            vertex[i+1][1] = i;
            vertex[j][1] = i;
            vertex[j-1][1] = i;
            p = i+1;
            while (pair1[p]==0) p++;
            q = pair1[p];
            vertex[p][1] = i;
            vertex[p-1][1] = i;
            vertex[q][1] = i;
            vertex[q+1][1] = i;
            
            vertex[p][2] = i;
            vertex[q][2] = i;
        } else if (type > 1) {
            vertex[i][1] = i;
            vertex[j][1] = i;
            
            /*
            p = i+1;
            do {
                while (pair1[p] == 0 && p<j) p++;
                if (p<j) {
                    q = pair1[p];
                    vertex[p][1] = i;
                    vertex[q][1] = i;
                    vertex[p][2] = i;
                    vertex[q][2] = i;
                    p = q+1;
                }
            } while (p<j);
            */
            // arcs nested in multiloop matters
        }
    }
    
    for (i=1; i<=n; i++) {
        if (pair2[i] == 0 || i>pair2[i]) continue;
        j = pair2[i];
        type = looptype(i, pair2);
        if (type == 0) {
            vertex[i][3] = i;
            vertex[j][3] = i;
            if (j-i-1 <= 4) {
                for (k=i+1; k<j; k++) {
                    vertex[k][3] = i;
                }
            } else {
                vertex[i+1][3] = i;
                vertex[j-1][3] = i;
            }
        } else if (type == 1) {
            vertex[i][3] = i;
            vertex[i+1][3] = i;
            vertex[j][3] = i;
            vertex[j-1][3] = i;
            p = i+1;
            while (pair2[p]==0) p++;
            q = pair2[p];
            vertex[p][3] = i;
            vertex[p-1][3] = i;
            vertex[q][3] = i;
            vertex[q+1][3] = i;
            
            vertex[p][4] = i;
            vertex[q][4] = i;
        } else if (type > 1) {
            vertex[i][3] = i;
            vertex[j][3] = i;
            
            /*
            p = i+1;
            do {
                while (pair2[p] == 0 && p<j) p++;
                if (p<j) {
                    q = pair2[p];
                    vertex[p][3] = i;
                    vertex[q][3] = i;
                    vertex[p][4] = i;
                    vertex[q][4] = i;
                    p = q+1;
                }
            } while (p<j);
            */
            // arcs nested in multiloop matters
        }
    }
    
    return vertex;
}
// practical model of loops



int ** update_v_loop(int *pair1, int *pair2, int **ref_v_loop)
{
    int i, j, n=pair1[0], **vertex;
    
    vertex = (int **) calloc(n+2, sizeof(int *));
    for (i=0; i<=n; i++) {
        vertex[i] = (int *) calloc (5, sizeof(int));
        for (j=0; j<=4; j++) vertex[i][j] = ref_v_loop[i][j];
    }
    
    for (i=1; i<=n; i++) {
        if (vertex[i][1]!=0 && pair1[vertex[i][1]] == 0) vertex[i][1] = 0;
        if (vertex[i][2]!=0 && pair1[vertex[i][2]] == 0) vertex[i][2] = 0;
        if (vertex[i][3]!=0 && pair2[vertex[i][3]] == 0) vertex[i][3] = 0;
        if (vertex[i][4]!=0 && pair2[vertex[i][4]] == 0) vertex[i][4] = 0;
    }
    return vertex;
}


void free_v_loop (int ** v_loop)
{
    int i, n=v_loop[0][0];
    for (i=0; i<=n; i++) {
        free(v_loop[i]);
    }
    free(v_loop);
}
//free the space of vertex loop table

void free_block(Bi_block *block)
{
    free(block->vertics);
    free(block->pair1);
    free(block->pair2);
    free(block->external);
    free(block);
}

int check_v_status(int v, int **v_loop, int ** ref_v_loop)
{
    int i;
    if (v_loop[v][1] == 0 && v_loop[v][2] == 0 && v_loop[v][3] == 0 && v_loop[v][4] == 0) {
        return -1;
    } else {
        for (i=1; i<=4; i++) {
            if (v_loop[v][i] != ref_v_loop[v][i]) return 0;
        }
    }
    return 1;
    // -1: not involved in any loop, to be removed
    // 0: expose
    // 1: complete
}

void remove_from_ext (int r, Bi_block *block)
{
    int i, n = block->external[0], j;
    for (i=1; i<=n; i++) {
        if (r == block->external[i]) {
            for (j=i; j<n; j++)
                block->external[j] = block->external[j+1];
            block->external[j] = 0;
            block->external[0]--;
        }
    }
}
//remove an ext vertex

int exist_ext (int r, int *ext)
{
    int i, n=ext[0];
    for (i=1; i<=n; i++) {
        if (r == ext[i]) return 1;
    }
    return 0;
}
//check if r exist in ext

void sort_ext(int *ext)
{
    int i, j, n=ext[0], swap;
    for (i=1; i<n; i++) {
        for (j=i; j<=n; j++) {
            if (ext[i]>ext[j]) {
                swap = ext[i];
                ext[i] = ext[j];
                ext[j] = swap;
            }
        }
    }
}

// remove an arc (p,q) update vertices set in Bi-struc
Bi_block *removal (int r, int side, int ** ref_v_loop, Bi_block *block)
{
    Bi_block *result;
    int ** new_v_loop, status;
    int n = block->vertics[0], i, j, s, k;
    
    result = new_Bi_block(n);
    for (i=0; i<=n; i++) {
        result ->vertics[i] = block->vertics[i];
        result ->pair1[i] = block->pair1[i];
        result ->pair2[i] = block->pair2[i];
    }
    for (i=0; i<=block->external[0]; i++) {
        result ->external[i] = block->external[i];
    }
    //copy block

    if (side == 1) {
        s = result->pair1[r];
        result->pair1[r] = 0;
        result->pair1[s] = 0;
    } else if (side ==2) {
        s = result->pair2[r];
        result->pair2[r] = 0;
        result->pair2[s] = 0;
    }
    //remove (r,s)
    
    new_v_loop = update_v_loop(result->pair1, result->pair2, ref_v_loop);
    
    /*
    for (i=1; i<=n; i++) {
        printf("%d: %d %d %d %d\n", i, new_v_loop[i][1], new_v_loop[i][2], new_v_loop[i][3], new_v_loop[i][4]);
    }
    */
    // debug display
    
    for (i=1; i<=n; i++) {
        status = check_v_status(i, new_v_loop, ref_v_loop);
        if (status == -1) {
            result->vertics[i] = 0;
            remove_from_ext(i, result);
        } else if (status == 0) {
            result->vertics[i] = 1;
            if (!exist_ext(i, result->external)) {
                k = result->external[0];
                result->external[0]++;
                result->external[k+1] = i;
            }
        } else {
            result->vertics[i] = 1;
        }
    }
    sort_ext(result->external);
    free_v_loop(new_v_loop);
    return result;
}
// remove arc (r,s) from the block, update vertex, pair1, pair2 and external

int differ_external (Bi_block *block, Bi_block *ref_block)
{
    int i, j, n = 0;
    n = block->external[0] + ref_block->external[0];
    for (i=1; i<=block->external[0]; i++) {
        for (j = 1; j<=ref_block->external[0]; j++) {
            if (block->external[i] == ref_block->external[j]) n--;
        }
    }
    
    return n;
}

int * Union_ext (Bi_block *B1, Bi_block *B2, int remove, int side, Bi_block *fullpair) //B2 remove an arc becomes B1, order important
{
    int l, flag =0, j, r_exist = 0, type, r,s, p, q;
    int *v_map;
    
    l = B1->external[0];
    if (side == 1) {
        if (remove < B2->pair1[remove]) {
            r = remove; s = B2->pair1[remove];
        } else {
            r = B2->pair1[remove]; s = remove;
        }
        type = looptype(r, fullpair->pair1);
    } else {
        if (remove < B2->pair2[remove]) {
            r = remove; s = B2->pair2[remove];
        } else {
            r = B2->pair2[remove]; s = remove;
        }
        type = looptype(r, fullpair->pair2);
    }
    //find base pair (remove, p) and its loop type
    v_map = (int *) calloc(B1->external[0]+B2->external[0]+10, sizeof(int));
    
    for (j=1; j<=B1->external[0]; j++) v_map[j] = B1->external[j];
    for (j=1; j<=B2->external[0]; j++) {
        if (!exist_ext(B2->external[j], B1->external)) {
            l++;
            v_map[l] = B2->external[j];
        }
    }
    //union external of B1 and B2
    
    v_map[0] = l;
    
    if (type <= 1) {
        if (!exist_ext(r, B1->external) && !exist_ext(r, B2->external)) {
            l++; v_map[l] = r; v_map[0]++;
        }
        if (!exist_ext(s, B1->external) && !exist_ext(s, B2->external)) {
            l++; v_map[l] = s; v_map[0]++;
        }
        if (!exist_ext(r+1, v_map) ) {
            l++; v_map[l] = r+1; v_map[0]++;
        }
        if (!exist_ext(s-1, v_map)) {
            l++; v_map[l] = s-1; v_map[0]++;
        }
        if (type == 0 && s-r-1 == 3) {
            l++; v_map[l] = r+2; v_map[0]++;
        }
        if (type == 0 && s-r-1 == 4) {
            l++; v_map[l] = r+2; v_map[0]++;
            l++; v_map[l] = r+3; v_map[0]++;
        }
        
        if (type == 1) {
            p  = r+1;
            if (side == 1) {
                while (fullpair->pair1[p]==0) p++;
                q = fullpair->pair1[p];
            } else {
                while (fullpair->pair2[p]==0) p++;
                q = fullpair->pair2[p];
            }
            if (!exist_ext(p-1, v_map) ) {
                l++; v_map[l] = p-1; v_map[0]++;
            }
            if (!exist_ext(q+1, v_map) ) {
                l++; v_map[l] = q+1; v_map[0]++;
            }
        }
    }
    sort_ext(v_map);
    return v_map;
}


Bi_block *next_block (int ** ref_v_loop, Bi_block *block, Bi_block *fullblock)
{
    int i, n = block->pair1[0], r, side = 0;
    int differ_ext, min_differ = n, det_r = 0;
    int *v_map;
    Bi_block *temp;
    
    for (r=1; r<=n; r++) {
        if (block->vertics[r] == 0) continue;
        
        if (block->pair1[r] && r<block->pair1[r]) {
            temp = removal(r, 1, ref_v_loop, block);
            // differ_ext = differ_external(temp, block);
        
            // v_map = Union_ext(temp, block, r, 1, fullblock);
            
            differ_ext = temp->external[0];
            //printf("Remove %d(1) differs %d \n", r, differ_ext);
            
            if (differ_ext < min_differ) {
                min_differ = differ_ext;
                det_r = r;
                side = 1;
            }
            
            free_block(temp);
        }
        
        if (block->pair2[r] && r<block->pair2[r]) {
            temp = removal(r, 2, ref_v_loop, block);
            // differ_ext = differ_external(temp, block);
            
            //v_map = Union_ext(temp, block, r, 2, fullblock);
            
            differ_ext = temp->external[0];
            //printf("Remove %d(2) differs %d \n", r, differ_ext);
            
            if (differ_ext < min_differ) {
                min_differ = differ_ext;
                det_r = r;
                side = 2;
            }
            
            free_block(temp);
        }
    }
    temp = removal(det_r, side, ref_v_loop, block);
    //printf("Remove %d(%d)\n", det_r, side);
    //for (i=1;i<=temp->external[0];i++) printf("%d ",temp->external[i]); printf("\n");
    return temp;
}
// given a block, decide which arc to be removed
                         

int empty_block(Bi_block *block)
{
    int i, n=block->vertics[0];
    for (i=1; i<=n; i++) {
        if (block->pair1[i] !=0 || block->pair2[i] !=0) return 0;
    }
    return 1;
}
//check whether block is empty

int removed_arc (Bi_block *block, Bi_block *removed_block, int *side)
{
    int i;
    for (i=1; i<=block->pair1[0]; i++) {
        if (block->pair1[i] != 0 && removed_block->pair1[i] == 0 ) {
            *side = 1;
            return i;
        }
    }
    
    for (i=1; i<=block->pair2[0]; i++) {
        if (block->pair2[i] != 0 && removed_block->pair2[i] == 0) {
            *side = 2;
            return i;
        }
    }

    return 0;
}




Bi_block *sequential_removal (char *str1, char *str2)
{
    int *pair1, *pair2, n, i;
    int ** ref_v_loop;
    Bi_block *head, *temp, *index, *fullblcok;
    
    n = strlen (str1);
    
    pair1 = structure2pair(str1);
    pair2 = structure2pair(str2);
    
    temp = new_Bi_block(n);
    for (i=0; i<=n; i++) {
        temp->vertics[i] = 1;
        temp->pair1[i] = pair1[i];
        temp->pair2[i] = pair2[i];
    }
    temp->vertics[0] = n;
    // ref_v_loop = loop_v(pair1, pair2);
    ref_v_loop = loop_v_truncated(pair1, pair2);

    
    /*
    for (i=1; i<=n; i++) {
        printf("%d: %d %d %d %d\n", i, ref_v_loop[i][1], ref_v_loop[i][2], ref_v_loop[i][3], ref_v_loop[i][4]);
    }
    */
    
    //debug display
    
    temp->external[0] = 1;
    
    temp->external[1] = 1;
    //temp->external[2] = n;
    // copy the origin structure to block, set ref v table
    
    head = temp;
    fullblcok = temp;
    
    // print_block(head);
    
    do {
        temp = next_block (ref_v_loop, head, fullblcok);
        temp->next = head;
        head = temp;
        
        // print_block(head);
        
    } while (!empty_block(temp));
    
    /*
    index = head;
    while (index!=NULL) {
        printf("%d, ", index->external[0]);
        //for (i=1; i<=index->external[0]; i++)
        //    printf("%d, ", index->external[i]);
        //printf("\n");
        index = index->next;
    }
    */
    //debug
    
    //printf("AAAAAGCUAACAUUC: %ld: %s\n", seq2index("AAAAAGCUAACAUUC"), index2seq(seq2index("AAAAAGCUAACAUUC"), 15));
    // debug index function
    return head;
    
}



int E_Loop (int r, Bi_block *block, int side, int *n_map, char * nucleotide)
{
    int energy = 0, k, i, j, p, l_type, a_type, a_type2, ni, nj, ni1, nj1, nt1, nu1, nt, nu, unpaired, ri, rj, ri1, rj1, rt1, ru1, rt, ru, t, u, n1, n2;
    char hairpin_str[30];
    int *rev_n_map;
    
    if (side == 1) {
        l_type = looptype(r, block->pair1);
    } else {
        l_type = looptype(r, block->pair2);
    }
    
    i = r;
    if (side == 1) {
        j = block->pair1[i];
    } else {
        j = block->pair2[i];
    }
    
    rev_n_map = (int *) calloc (block->pair1[0]+5, sizeof(int));
    rev_n_map[0] = block->pair1[0];
    for (k=1; k<=rev_n_map[0]; k++) rev_n_map[k] = 0;
    for (k=1; k<=n_map[0]; k++) rev_n_map[n_map[k]] = k;
    //reverse map
    
    ri = rev_n_map[i];
    rj = rev_n_map[j];
    ni = nucleotide2code(nucleotide[ri-1]);
    nj = nucleotide2code(nucleotide[rj-1]);
    a_type = BP_pair[ni][nj];
    if (a_type == 0) {
        free(rev_n_map);
        return MAXENG;
    }
    
    if (l_type == 0) { // hairpin
        ri1 = rev_n_map[i+1];
        rj1 = rev_n_map[j-1];
        ni1 = nucleotide2code(nucleotide[ri1-1]);
        nj1 = nucleotide2code(nucleotide[rj1-1]);
        
        if (j-i-1 <=4) {
            for (p=0; p<j-i+1; p++) {
                hairpin_str[p] = nucleotide[rev_n_map[i+p]-1];
            }
        } else {
            hairpin_str[0] = nucleotide[ri-1];
            hairpin_str[1] = nucleotide[ri1-1];
            hairpin_str[j-i] = nucleotide[rj-1];
            hairpin_str[j-i-1] = nucleotide[rj1-1];
            for (p=2; p<j-i-1; p++) hairpin_str[p] = 'C';
        }
        hairpin_str[j-i+1] = 0;
        energy = HairpinE(j-i-1, a_type, ni1, nj1, hairpin_str);
    } else if (l_type == 1) { // interior loop
        t = i+1; n1=0; n2=0;
        if (side == 1) {
            while (block->pair1[t]==0 && t<j) {
                t++;
                n1++;
            }
            u = block->pair1[t];
            n2 = j-u-1;
        } else {
            while (block->pair2[t]==0 && t<j) {
                t++;
                n1++;
            }
            u = block->pair2[t];
            n2 = j-u-1;
        }
        // identify arc (t,u)
        
        ri1 = rev_n_map[i+1];
        rj1 = rev_n_map[j-1];
        ni1 = nucleotide2code(nucleotide[ri1-1]);
        nj1 = nucleotide2code(nucleotide[rj1-1]);
        rt1 = rev_n_map[t-1];
        ru1 = rev_n_map[u+1];
        nt1 = nucleotide2code(nucleotide[rt1-1]);
        nu1 = nucleotide2code(nucleotide[ru1-1]);
        rt = rev_n_map[t];
        ru = rev_n_map[u];
        nt = nucleotide2code(nucleotide[rt-1]);
        nu = nucleotide2code(nucleotide[ru-1]);
        
        
        a_type2 = BP_pair[nu][nt];
        if (a_type2 == 0) {
            free(rev_n_map);
            return MAXENG;
        }
        
        energy = LoopEnergy(n1, n2, a_type, a_type2, ni1, nj1, nt1, nu1);
    } else if (l_type >= 2) { // multiloop

        unpaired = 0;
        if (side == 1) {
            p = i+1;
            do {
                while (block->pair1[p]==0) {
                    p++;
                    unpaired++;
                }
                if (p<j) {
                    if (rev_n_map[p] && rev_n_map[block->pair1[p]]) {
                        a_type2 = BP_pair[rev_n_map[p]][rev_n_map[block->pair1[p]]];
                        if (a_type2 == 0) {
                            free(rev_n_map);
                            return MAXENG;
                        }
                    }
                    p = block->pair1[p]+1;
                }
            } while (p<j);
        } else {
            p = i+1;
            do {
                while (block->pair2[p]==0) {
                    p++;
                    unpaired++;
                }
                if (p<j) {
                    if (rev_n_map[p] && rev_n_map[block->pair2[p]]) {
                        a_type2 = BP_pair[rev_n_map[p]][rev_n_map[block->pair2[p]]];
                        if (a_type2 == 0) {
                            free(rev_n_map);
                            return MAXENG;
                        }
                    }
                    p = block->pair2[p]+1;
                }
            } while (p<j);
        }
        
        energy = P->MLclosing + P->MLintern[a_type] + P->MLbase * unpaired;
    }
    free(rev_n_map);
    return energy;
}

int array_common(int *a, int *b)
{
    int i, j, cnt = 0;
    for (i=1; i<=a[0]; i++) {
        for (j=1; j<=b[0]; j++) {
            if (a[i] == b[j]) cnt++;
        }
    }
    return cnt;
}

int compatible (char *seq, int *pair) {
    int i, type;
    for (i=1; i<pair[0]; i++) {
        if (pair[i] && i<pair[i]) {
            type = BP_pair[nucleotide2code(seq[i-1])][nucleotide2code(seq[pair[i]-1])];
            if (type == 0) {
                printf("[%d %d]: [%c %c]\n", i, pair[i], seq[i-1], seq[pair[i]-1]);
                return 0;
            }
        }
    }
    return 1;
}



int est_max_index(Bi_block *head)
{
    Bi_block *index, *index2, *fullblock;
    int max_index = 0, *v_map, side, r;
    fullblock = head;
    while (fullblock->next!=NULL) fullblock = fullblock->next;
    
    index = head;
    while (index->next!=NULL) {
        index2 = index->next;
        r = removed_arc(index2, index, &side);
        //printf("%s\n%s\n\n", pair2structure(index->pair1), pair2structure(index->pair2));
        v_map = Union_ext(index, index2, r, side, fullblock);
        if (v_map[0] > max_index) max_index = v_map[0];
        index = index->next;
    }
    free(v_map);
    return max_index; 
}


int mfe_double_seq(Bi_block *head, char *seq, char *str1, char *str2)
{
    int **Q, length, i, side, r, c, l;
    //Q[i][code] is the result in DP. i is the step, code is the encoded index
    Bi_block *index, * block_arrray[1000], *index2, *fullblock;
    int *v_map,*rev_n_map;
    int energy, temp_energy, k, loop_energy, energy1 = 0, energy2 = 0;
    char *N, Nx[50], Ny[50], *N_temp, fold_struc[1000];
    unsigned long long loop_index, index_x, index_y, max_index = 0, m_index = 0;
    float fold_energy;
    
    //debug seq
    
    update_fold_params();
    P = scale_parameters();
    
    length = 0;
    index = head;
    while (index !=NULL) {
        block_arrray[length] = index;
        if (index->external[0] > max_index) max_index = index->external[0];
        length ++;
        index = index->next;
    }
    fullblock = block_arrray[length-1];
    //compute the length of iteration in DP
    
    Q = (int **) calloc (length+5, sizeof (int *));
    for (k=0; k<length; k++) {
        m_index = pow(4, block_arrray[k]->external[0]);
        Q[k] = (int *) calloc (m_index+5, sizeof(int));
        for (loop_index=0; loop_index<=m_index; loop_index++) Q[k][loop_index] = 10000;
    }
    Q[0][0] = 0;
    // initiate
    
    for (k=0; k<length-1; k++) {
        index = block_arrray[k];
        index2 = block_arrray[k+1];
        r = removed_arc(index2, index, &side);
        
        v_map = Union_ext(index, index2, r, side, fullblock);
        m_index = pow(4, v_map[0]);
        
        
        //for (i=1; i<=index->external[0];i++) printf("%d ", index->external[i]); printf("\n");
        //for (i=1; i<=index->next->external[0];i++) printf("%d ", index->next->external[i]); printf("\n");
        //for (i=1; i<=v_map[0]; i++) printf("%d ", v_map[i]); printf("\n");
        // debug show v_map results
        
        rev_n_map = (int *) calloc (index->pair1[0]+5, sizeof(int));
        rev_n_map[0] = index->pair1[0];
        for (i=1; i<=rev_n_map[0]; i++) rev_n_map[i] = 0;
        for (i=1; i<=v_map[0]; i++) rev_n_map[v_map[i]] = i;
        //construct reverse map
        
        
        // for (i=0; i<v_map[0]; i++) N[i] = test_str[v_map[i+1]-1]; //test string
        for (loop_index = 0; loop_index < m_index; loop_index++) {
            N = index2seq(loop_index, v_map[0]);
           
            
            // for (i=0; i<=index->pair1[0]; i++) printf("%d ", index->pair1[i]); printf("\n");
            // for (i=0; i<index->next->pair1[0]; i++) printf("%d ", index->next->pair1[i]); printf("\n");
            // for (i=0; i<=index->pair2[0]; i++) printf("%d ", index->pair2[i]); printf("\n");
            // for (i=0; i<index->next->pair2[0]; i++) printf("%d ", index->next->pair2[i]); printf("\n");
            // for (i=0; i<index->vertics[0]; i++) printf("%d ", index->vertics[i]); printf("\n");
            //printf("r=%d side=%d\n", r, side); // show removal arc
            //debug
            
            for (i=1; i<=index2->external[0]; i++) Nx[i-1] = N[rev_n_map[index2->external[i]]-1];
            for (i=1; i<=index->external[0]; i++) Ny[i-1] = N[rev_n_map[index->external[i]]-1];
            Nx[index2->external[0]] = 0;
            Ny[index->external[0]] = 0;
            //construct Nx and Ny
            
            index_x = seq2index(Nx);
            index_y = seq2index(Ny);
            
            loop_energy = E_Loop(r, block_arrray[length-1], side, v_map, N);
           
            if (side == 1) {
                if (loop_energy<1000) {
                    loop_energy *= c1;
                }
            }
            else if (side ==2) {
                if (loop_energy<1000) {
                    loop_energy *= c2;
                }
            }
            //scale
            
            
            if (loop_energy + Q[k][index_y] < Q[k+1][index_x]) {
                Q[k+1][index_x] = loop_energy + Q[k][index_y];
                // fprintf(Output, "%d %d %llu %d\n", k, loop_energy, index_x, Q[k][index_y]);
            }
            
            free(N);
        }
        
        free(rev_n_map);
        free(v_map);
    }
    energy = Q[length-1][0];
    index_x = 0;
    
    //show energy of mfe
    
    for (i=0; i<4; i++) {
        if (Q[length-1][i]<energy) {
            index_x = i;
            energy = Q[length-1][i];
        }
    }
    
    
    fprintf(Output, "Total Energy: %d\n", energy);
    
    N = index2seq(index_x, 1);
    Nx[0] = N[0]; Nx[1] = 0;
    free(N);
    
    for (k=length-1; k>0; k--) {
        index2 = block_arrray[k];
        index = block_arrray[k-1];
        r = removed_arc(index2, index, &side);
        v_map = Union_ext(index, index2, r, side, fullblock);
        //c = array_common(index->external, index2->external);
        m_index = pow(4, v_map[0]-index2->external[0]);
        
        
        rev_n_map = (int *) calloc (index->pair1[0]+5, sizeof(int));
        rev_n_map[0] = index->pair1[0];
        for (i=1; i<=v_map[0]; i++) rev_n_map[i] = 0;
        for (i=1; i<=v_map[0]; i++) rev_n_map[v_map[i]] = i;
        //construct reverse map
        
        // for (i=0; i<v_map[0]; i++) N[i] = test_str[v_map[i+1]-1]; //test string
        for (loop_index = 0; loop_index < m_index; loop_index++) {
            N_temp = index2seq(loop_index, v_map[0]-index2->external[0]);
            
            
            // for (i=0; i<=index->pair1[0]; i++) printf("%d ", index->pair1[i]); printf("\n");
            // for (i=0; i<index->next->pair1[0]; i++) printf("%d ", index->next->pair1[i]); printf("\n");
            // for (i=0; i<=index->pair2[0]; i++) printf("%d ", index->pair2[i]); printf("\n");
            // for (i=0; i<index->next->pair2[0]; i++) printf("%d ", index->next->pair2[i]); printf("\n");
            // for (i=0; i<index->vertics[0]; i++) printf("%d ", index->vertics[i]); printf("\n");
            //printf("r=%d side=%d\n", r, side); // show removal arc
            //debug
    
            N = (char *) calloc (v_map[0]+5, sizeof(char));
            for (i=0; i<v_map[0]; i++) N[i] = '_';
            N[v_map[0]] = 0;
            for (i=1; i<=index2->external[0]; i++) N[rev_n_map[index2->external[i]]-1] = Nx[i-1];
            l = 0;
            for (i=1; i<=v_map[0]; i++) {
                if (N[i-1] == '_') {
                    N[i-1] = N_temp[l];
                    l++;
                }
            }
            for (i=1; i<=index->external[0]; i++) {
                Ny[i-1] = N[rev_n_map[index->external[i]]-1];
            }
            free(N_temp);
            Ny[index->external[0]] = 0;
            //construct N and Ny
            
            // index_x = seq2index(Nx);
            index_y = seq2index(Ny);
            
            loop_energy = E_Loop(r, block_arrray[length-1], side, v_map, N);
            
            if (side == 1) {
                if (loop_energy<1000) {
                    loop_energy *= c1;
                }
            }
            else if (side ==2) {
                if (loop_energy<1000) {
                    loop_energy *= c2;
                }
            }
            //scale
            
            if (loop_energy + Q[k-1][index_y] == Q[k][index_x]) {
                for (i=1; i<=v_map[0]; i++) {
                    seq[v_map[i]-1] = N[rev_n_map[v_map[i]]-1];
                }
                temp_energy = Q[k-1][index_y];
                free(N);
                
                index_x = index_y;
                for (i=0; i<50; i++) Nx[i] = Ny[i];
                
                if (side == 1) energy1 += loop_energy;
                else if (side == 2) energy2 += loop_energy;
                break;
            }
            free(N);
        }
        
        if (loop_index>=m_index) {
            printf("tracing back error!\n");
        }

    
        free(rev_n_map);
        free(v_map);
    }
    
    for (i=0;i<strlen(str1);i++) {
        if (seq[i] == '_') seq[i] = randomnucleotide();
    }
    
    // printf("Energy Side1: %d    Energy Side1: %d\n", energy1, energy2);
    
    //fold_energy = fold(seq, fold_struc);
    // printf("%s %f\n", fold_struc, fold_energy);
    //fprintf(Output, "MFE seq: %s\n", seq);
    //fprintf(Output, "Fold: %s %f\n", fold_struc, fold_energy);
    return energy;
}


double **PF_double_seq(Bi_block *head, char *str1, char *str2)
{
    int length, i, side, r;
    double **PF;
    //PF[i][code] is the partition function in DP. i is the step, code is the encoded index
    Bi_block *index, * block_arrray[1000], *index2, *fullblock;
    int *v_map,*rev_n_map;
    int k, loop_energy;
    char *N, Nx[50], Ny[50];
    unsigned long long loop_index, index_x, index_y, max_index = 0, m_index = 0;
    FILE *fp;

    
    update_fold_params();
    P = scale_parameters();
    
    length = 0;
    index = head;
    while (index !=NULL) {
        block_arrray[length] = index;
        if (index->external[0] > max_index) max_index = index->external[0];
        length ++;
        index = index->next;
    }
    fullblock = block_arrray[length-1];
    //compute the length of iteration in DP
    
    PF = (double **) calloc (length+5, sizeof (double *));
    for (k=0; k<length; k++) {
        m_index = pow(4, block_arrray[k]->external[0]);
        PF[k] = (double *) calloc (m_index+5, sizeof(double));
        for (loop_index=0; loop_index<=m_index; loop_index++) PF[k][loop_index] = 0;
    }
    PF[0][0] = 1;
    PF[length] = (double *) calloc (10, sizeof(double));
    // initiate
    
    for (k=0; k<length-1; k++) {
        index = block_arrray[k];
        index2 = block_arrray[k+1];
        r = removed_arc(index2, index, &side);
        v_map = Union_ext(index, index2, r, side,fullblock);
        m_index = pow(4, v_map[0]);
        
        //for (i=1; i<=index->external[0];i++) printf("%d ", index->external[i]); printf("\n");
        //for (i=1; i<=index->next->external[0];i++) printf("%d ", index->next->external[i]); printf("\n");
        //for (i=1; i<=v_map[0]; i++) printf("%d ", v_map[i]); printf("\n");
        // debug show v_map results
        
        rev_n_map = (int *) calloc (index->pair1[0]+5, sizeof(int));
        rev_n_map[0] = index->pair1[0];
        for (i=1; i<=rev_n_map[0]; i++) rev_n_map[i] = 0;
        for (i=1; i<=v_map[0]; i++) rev_n_map[v_map[i]] = i;
        //construct reverse map
        
        // for (i=0; i<v_map[0]; i++) N[i] = test_str[v_map[i+1]-1]; //test string
        for (loop_index = 0; loop_index < m_index; loop_index++) {
            N = index2seq(loop_index, v_map[0]);
            
            
            // for (i=0; i<=index->pair1[0]; i++) printf("%d ", index->pair1[i]); printf("\n");
            // for (i=0; i<index->next->pair1[0]; i++) printf("%d ", index->next->pair1[i]); printf("\n");
            // for (i=0; i<=index->pair2[0]; i++) printf("%d ", index->pair2[i]); printf("\n");
            // for (i=0; i<index->next->pair2[0]; i++) printf("%d ", index->next->pair2[i]); printf("\n");
            // for (i=0; i<index->vertics[0]; i++) printf("%d ", index->vertics[i]); printf("\n");
            //printf("r=%d side=%d\n", r, side); // show removal arc
            //debug
            
            for (i=1; i<=index2->external[0]; i++) Nx[i-1] = N[rev_n_map[index2->external[i]]-1];
            for (i=1; i<=index->external[0]; i++) Ny[i-1] = N[rev_n_map[index->external[i]]-1];
            Nx[index2->external[0]] = 0;
            Ny[index->external[0]] = 0;
            //construct Nx and Ny
            
            index_x = seq2index(Nx);
            index_y = seq2index(Ny);
            
            loop_energy = E_Loop(r, block_arrray[length-1], side, v_map, N);
            
            if (side == 1) {
                if (loop_energy<1000) {
                    loop_energy *= c1;
                }
            }
            else if (side ==2) {
                if (loop_energy<1000) {
                    loop_energy *= c2;
                }
            }
            //scale
            
            if (loop_energy < 0.1*MAXENG) {
                PF[k+1][index_x] += EXP(loop_energy) * PF[k][index_y];
            }
            
            //fprintf(Output, "%d %d %llu %e\n", k, loop_energy, index_x, PF[k][index_y]);
            
            free(N);
        }
        
        free(rev_n_map);
        free(v_map);
    }
    index_x = 0;
    
    //show energy of mfe
    
    // printf("Total Energy: %d\n", energy);
    
    return PF;
}




patten *Biseq_sampler (Bi_block *head, double **PF, int mul, char *str1, char *str2)
{
    
    int seq_length = strlen(str1), length, i, side, r, c, l, m, ud = 0;
    float seed;
    Bi_block *index, * block_arrray[1000], *index2, *fullblock;
    int *v_map, *rev_n_map;
    int energy, temp_energy, k, loop_energy, energy1 = 0, energy2 = 0;
    char *N, Nx[50], Ny[50], *N_temp, *sample_seq;
    unsigned long long loop_index, index_x, index_y, max_index = 0, m_index = 0;
    double PFtemp, PFpro;
    patten *List = NULL, *seq_list_index;
    FILE *fp;
    
    length = 0;
    index = head;
    while (index !=NULL) {
        block_arrray[length] = index;
        if (index->external[0] > max_index) max_index = index->external[0];
        length ++;
        index = index->next;
    }
    fullblock = block_arrray[length-1];
    //compute the length of iteration in DP
    
    for (m=0; m<mul; m++) {
        
        sample_seq = (char *) calloc (seq_length+5, sizeof(char));
        for (i=0; i<seq_length; i++) sample_seq[i] = '_';
        sample_seq[seq_length] = 0;
        seed = (1-drand48());
        PFtemp = seed * (PF[length-1][0] + PF[length-1][1] + PF[length-1][2] + PF[length-1][3]);
        if (PFtemp < PF[length-1][0]) {
            index_x = 0;
        } else if (PFtemp < PF[length-1][0] + PF[length-1][1]) {
            index_x = 1;
        } else if (PFtemp < PF[length-1][0] + PF[length-1][1] + PF[length-1][2]) {
            index_x = 2;
        } else {
            index_x = 3;
        }
        N = index2seq(index_x, 1);
        Nx[0] = N[0]; Nx[1] = 0;
        sample_seq[0] = Nx[0];
        free(N);

        for (k=length-1; k>0; k--) {
            index2 = block_arrray[k];
            index = block_arrray[k-1];
            r = removed_arc(index2, index, &side);
            v_map = Union_ext(index, index2, r, side, fullblock);
            // c = array_common(index->external, index2->external);
            m_index = pow(4, v_map[0]-index2->external[0]);
        
        
            rev_n_map = (int *) calloc (seq_length+5, sizeof(int));
            rev_n_map[0] = seq_length;
            for (i=1; i<=v_map[0]; i++) rev_n_map[i] = 0;
            for (i=1; i<=v_map[0]; i++) rev_n_map[v_map[i]] = i;
            //construct reverse map
        
            PFtemp = 0;
            seed = (1-drand48());
            PFpro = seed * PF[k][index_x];
            
            // for (i=0; i<v_map[0]; i++) N[i] = test_str[v_map[i+1]-1]; //test string
            for (loop_index = 0; loop_index < m_index; loop_index++) {
                N_temp = index2seq(loop_index, v_map[0]-index2->external[0]);
                
            
                // for (i=0; i<=index->pair1[0]; i++) printf("%d ", index->pair1[i]); printf("\n");
                // for (i=0; i<index->next->pair1[0]; i++) printf("%d ", index->next->pair1[i]); printf("\n");
                // for (i=0; i<=index->pair2[0]; i++) printf("%d ", index->pair2[i]); printf("\n");
                // for (i=0; i<index->next->pair2[0]; i++) printf("%d ", index->next->pair2[i]); printf("\n");
                // for (i=0; i<index->vertics[0]; i++) printf("%d ", index->vertics[i]); printf("\n");
                //printf("r=%d side=%d\n", r, side); // show removal arc
                //debug
            
                N = (char *) calloc (v_map[0]+5, sizeof(char));
                for (i=0; i<v_map[0]; i++) N[i] = '_';
                N[v_map[0]] = 0;
                for (i=1; i<=index2->external[0]; i++) N[rev_n_map[index2->external[i]]-1] = Nx[i-1];
                l = 0;
                for (i=1; i<=v_map[0]; i++) {
                    if (N[i-1] == '_') {
                        N[i-1] = N_temp[l];
                        l++;
                    }
                }
                for (i=1; i<=index->external[0]; i++) {
                    Ny[i-1] = N[rev_n_map[index->external[i]]-1];
                }
                free(N_temp);
                Ny[index->external[0]] = 0;
                //construct N and Ny
            
                // index_x = seq2index(Nx);
                index_y = seq2index(Ny);
            
                loop_energy = E_Loop(r, block_arrray[length-1], side, v_map, N);
                
                if (side == 1) {
                    if (loop_energy<1000) {
                        loop_energy *= c1;
                    }
                }
                else if (side ==2) {
                    if (loop_energy<1000) {
                        loop_energy *= c2;
                    }
                }
                //scale
                
                if (loop_energy < 0.1*MAXENG) {
                    PFtemp += EXP(loop_energy) * PF[k-1][index_y];
                }
                
                if (PFtemp > PFpro) {
                    for (i=1; i<=v_map[0]; i++) {
                        sample_seq[v_map[i]-1] = N[rev_n_map[v_map[i]]-1];
                    }
                    free(N);
                    
                    index_x = index_y;
                    for (i=0; i<50; i++) Nx[i] = Ny[i];
                    
                    break;
                }
                
                free(N);
            }
            if (loop_index>=m_index) {
                printf("Error! Can't match!\n");
            }
            
            free(rev_n_map);
            free(v_map);
        }
        ud = 0;
        for (i=0;i<seq_length;i++) {
            if (sample_seq[i] == '_') {
                sample_seq[i] = randomnucleotide();
                ud++;
            }
        }
        
        List = patten_detector(sample_seq, 0, seq_length-1, 1, List);
        // printf("Biseq: %s  %d %d\n", sample_seq, compatible(sample_seq, fullblock->pair1),  compatible(sample_seq, fullblock->pair2));
        
    
        free(sample_seq);
    }
    
    
    //printf("%d\n", ud);
    PFtemp =  PF[length-1][0] + PF[length-1][1] + PF[length-1][2] + PF[length-1][3];
    //printf("PF: %e\n", PFtemp);
    PFtemp *= pow(4, ud);
    
    PF[length][0] = PFtemp;
    
    return List;

}







patten *bicompatible_generator (char *str1, char *str2, int mul)
{
    char *seq;
    int *pair1, *pair2, i, n, j, m, total, length, p, k, ub=0;
    int F[200], C[200], *chain, num_chain = 0, num_cycle = 0, chain_head[200], cycle_head[200];
    float seed;
    patten *seq_list = NULL, *temp, *index;
    double num_compatible = 1;
    
    n = strlen(str1);
    seq = (char *) calloc (n+5, sizeof(char));
    chain = (int *) calloc (n+5, sizeof(int));
    for (i=0; i<n; i++) seq[i] = '_';
    seq[n] = 0;
    pair1 = structure2pair(str1);
    pair2 = structure2pair(str2);
    
    
    F[0] = n;
    F[1] = 1;
    F[2] = 1;
    for (i=3; i<=n; i++) F[i] = F[i-1] + F[i-2];
    // prnitf("%d\n", F[n]);
    
    for (i=0; i<100; i++) {
        chain_head[i] = 0;
    }
    
    
    chain_decomposition(chain, chain_head, cycle_head, pair1, pair2);
    num_chain = chain_head[0];
    num_cycle = cycle_head[0];
    
    printf("%d\n", num_chain);
      for (i=1; i<=num_chain; i++) {
          j = chain_head[i];
          while (j!=-1) {
              printf("%d ", j);
              j = chain[j];
              
          }
          printf("\n");
     }
    
    
    
    printf("%d\n", num_cycle);
    for (i=1; i<=num_cycle; i++) {
        j = cycle_head[i];
        do {
            printf("%d ", j);
            j = chain[j];
            
        } while (j!=cycle_head[i]);
        printf("\n");
    }
    // for (i=1; i<=n; i++) printf("%d ", chain[i]); printf("\n");
    
    
    for (i=1; i<=n; i++) {
        if (chain[i] == 0) {
            ub++;
        }
    }
    for (i=1; i<=num_chain; i++) {
        j = chain_head[i];
        length = 0;
        do {
            length ++;
            j = chain[j];
        } while (j!=-1);
        num_compatible *= 2*(F[length] + F[length+1]);
    }
    for (i=1; i<=num_cycle; i++) {
        j = cycle_head[i];
        length = 0;
        do {
            length ++;
            j = chain[j];
        } while (j!=cycle_head[i]);
        num_compatible *= 2*(F[length+1] + F[length-1]);
    }
    num_compatible *= pow(4,ub);
    //compute the compatible space
    printf("Total Compatible: %e\n", num_compatible);
    
    //generate sequence
    for (m=1; m<=mul; m++) {
        for (i=1; i<=n; i++) {
            if (chain[i] == 0) {
                seed=(1-drand48())*4;
                if (seed < 1) seq[i-1] = 'A';
                else if (seed < 2) seq[i-1] = 'C';
                else if (seed < 3) seq[i-1] = 'G';
                else seq[i-1] = 'U';
            }
        }
        for (i=1; i<=num_chain; i++) {
            j = chain_head[i];
            length = 0;
            do {
                length ++;
                j = chain[j];
            } while (j!=-1);
            
            j = chain_head[i];
            k = 0;
            total = 2*(F[length] + F[length+1]);
            seed=(1-drand48())*total;
            
            if (seed < F[length]) seq[j-1] = 'A';
            else if (seed < 2*F[length]) seq[j-1] = 'C';
            else if (seed < 2*F[length]+F[length+1]) seq[j-1] = 'G';
            else seq[j-1] = 'U';
            
            while (j!=-1) {
                p = j;
                j = chain[j];
                k++;
                
                if (seq[p-1] == 'A') seq[j-1] = 'U';
                if (seq[p-1] == 'C') seq[j-1] = 'G';
                if (seq[p-1] == 'G') {
                    total = F[length+1-k] + F[length-k];
                    seed=(1-drand48())*total;
                    if (seed < F[length+1-k]) seq[j-1] = 'U';
                    else seq[j-1] = 'C';
                }
                if (seq[p-1] == 'U') {
                    total = F[length+1-k] + F[length-k];
                    seed=(1-drand48())*total;
                    if (seed < F[length+1-k]) seq[j-1] = 'G';
                    else seq[j-1] = 'A';
                }
            }
        }
        
        for (i=1; i<=num_cycle; i++) {
            j = cycle_head[i];
            length = 0;
            do {
                length ++;
                j = chain[j];
            } while (j!=cycle_head[i]);
            
            j = cycle_head[i];
            k = 0;
            total = 2*(F[length+1] + F[length-1]);
            seed=(1-drand48())*total;
            
            if (seed < F[length-1]) seq[j-1] = 'A';
            else if (seed < 2*F[length-1]) seq[j-1] = 'C';
            else if (seed < 2*F[length-1]+F[length+1]) seq[j-1] = 'G';
            else seq[j-1] = 'U';
            
            if (length < 4) {
                if (seq[j-1] == 'A') seq[chain[j]-1] = 'U';
                if (seq[j-1] == 'C') seq[chain[j]-1] = 'G';
                if (seq[j-1] == 'G') seq[chain[j]-1] = 'C';
                if (seq[j-1] == 'U') seq[chain[j]-1] = 'A';
                j = chain[j];
            }
            
            while (j!=cycle_head[i]) {
                p = j;
                j = chain[j];
                k++;
                
                if (seq[p-1] == 'A') seq[j-1] = 'U';
                if (seq[p-1] == 'C') seq[j-1] = 'G';
                if ((seq[p-1] == 'G' || seq[p-1] == 'U') && (seq[cycle_head[i]-1] == 'A' || seq[cycle_head[i]-1] == 'C')) {
                    total = F[length-k-1];
                    seed=(1-drand48())*total;
                    if (seed < F[length-k-3]) {
                        if (seq[p-1] == 'G') seq[j-1] = 'C';
                        else if (seq[p-1] =='U') seq[j-1] = 'A';
                    } else {
                        if (seq[p-1] == 'G') seq[j-1] = 'U';
                        else if (seq[p-1] == 'U') seq[j-1] = 'G';
                    }
                }
                if ((seq[p-1] == 'G' || seq[p-1] == 'U') && (seq[cycle_head[i]-1] == 'U' || seq[cycle_head[i]-1] == 'G')) {
                    total = F[length-k];
                    seed=(1-drand48())*total;
                    if (seed < F[length-k-2]) {
                        if (seq[p-1] == 'G') seq[j-1] = 'C';
                        else if (seq[p-1] =='U') seq[j-1] = 'A';
                    } else {
                        if (seq[p-1] == 'G') seq[j-1] = 'U';
                        else if (seq[p-1] == 'U') seq[j-1] = 'G';
                    }
                }
            }
        }
        
        seq_list = patten_detector(seq, 0, n-1, 1, seq_list);
    }
    /*
     index = seq_list;
     while (index!=NULL) {
     printf("%s  %d\n", index->pat, index->cnt);
     index = index->next;
     }
     */
    
    
    int ** vertex;
    vertex = loop_v(pair1,pair2);
    //for (i=1; i<=n; i++) {
    //    printf("%d: %d %d %d %d\n", i, vertex[i][1], vertex[i][2], vertex[i][3], vertex[i][4]);
    //}
    
    free(pair1);
    free(pair2);
    free(chain);
    return seq_list;
}
//uniform generate bicompatible sequence from two structrus










