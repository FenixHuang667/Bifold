//
//  head.h
//  Bifold
//
//  Created by Fenix Huang on 7/11/18.
//  Copyright © 2018 Fenix Huang. All rights reserved.
//

#ifndef head_h
#define head_h

typedef struct patten patten;

struct patten
{
    char *pat;
    unsigned long cnt;
    patten *next;
};

#endif /* head_h */

extern char Outputfilename[500];
extern char *Outfile;
extern FILE *Output;
extern float c1;
extern float c2;
