

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include "graph.h"

#include "mmio.h"

//#define FAST_SEQUENTIAL

/*
 * Load graph market file
 */

typedef struct {
    size_t row, col;
    double val;
} mmelement_t;

int
compareMMElement( const void* a, const void* b ) {
    mmelement_t* elemA = (mmelement_t*)a;
    mmelement_t* elemB = (mmelement_t*)b;

    int i = elemA->row - elemB->row;

    return i ? i : elemA->col - elemB->col;
}

static bool
read_matrix_market2(const char *filename, graph_t* g )
{
    const int bs = 4;
    FILE *fh = fopen(filename, "r");
    if (!fh)
    {
        perror("fopen");
        return false;
    }

    //printf("Reading graph '%s'...\n", filename);

    MM_typecode matcode;
    if (mm_read_banner(fh, &matcode) != 0)
    {
        fprintf(stderr, "(e) mm_read_banner failed\n");
        return false;
    }

    int M, N, size;
    int ret_code;
    ret_code = mm_read_mtx_crd_size(fh, &M, &N, &size);
    if (ret_code != 0)
    {
        fprintf(stderr, "(e) mm_read_mtx_crd_size failed\n");
        return false;
    }

    if( M != N || !mm_is_symmetric(matcode) ) {
        fprintf(stderr, "(e) Matrix is not an undirected graph!\n");
        return false;
    }

    g->size =0;
    g->m =M;
    g->nodePtr = calloc( M, sizeof( node_t ) );


    for( int i =0; i < M; i++ ) {
        node_t* n =&g->nodePtr[i];
        // We use the fact that calloc lazily allocates this memory
        n->size =bs; n->degree =0;
        n->index =i;
        n->edges = calloc( bs, sizeof(idx_t) );
        n->weights=calloc( bs, sizeof(float) );
    }

    for (int i = 0; i < size; i++)
    {
        int row, col;
        double val;

        if (mm_is_pattern(matcode))
        {
            fscanf(fh, "%d %d\n", &row, &col);
            val = 1.0;
        }
        else
            fscanf(fh, "%d %d %lg\n", &row, &col, &val);

        row--; /* adjust from 1-based to 0-based */
        col--;

        //if( val == 0.0 ) continue;

#ifdef FAST_SEQUENTIAL
        bool swap =true;
#else
        bool swap =false;
#endif
        node_t* n;
L0:
        if( row == col ) continue;

        n =&g->nodePtr[ swap ? col : row ];
        if( n->degree >= n->size ) {
            n->size += bs;
            n->edges =realloc( n->edges, n->size * sizeof(idx_t) ); 
            n->weights =realloc( n->weights, n->size * sizeof(float) ); 
        }
        n->edges[n->degree]   = swap ? row : col;
        n->weights[n->degree] = val;
        n->degree++;
        g->size++;

#ifndef FAST_SEQUENTIAL
        if (mm_is_symmetric(matcode) && row != col && !swap ) {
            swap =true;
            goto L0;
        }
#endif
    }

    fclose(fh);
    for( int i =0; i < M; i++ ) {
        node_t* n =&g->nodePtr[i];
        n->size =n->degree;
        n->edges = realloc( n->edges,   n->degree * sizeof(idx_t) );
        n->weights=realloc( n->weights, n->degree * sizeof(float) );
    }

    return true;
}

void
graph_dump( const graph_t* mat ) {

    for (int i = 0; i < mat->m; ++i)
    {
        if( mat->nodePtr[i].index != i ) continue; // contracted node
        node_t* n =&mat->nodePtr[i];
        for (idx_t idx =0; idx < n->size; ++idx)
        {
            printf("(%d)->(%d) weight %f\n", i, n->edges[idx], n->weights[idx] );
        }
    }
}


bool
graph_loadMM(const char *filename, graph_t *mat ) {

    if( !read_matrix_market2( filename, mat ) )
        return false;


    return true;
}

void
graph_free( graph_t* mat ) {
    for( size_t i =0; i < mat->m; i++ ) {
        free( mat->nodePtr[i].weights );
        free( mat->nodePtr[i].edges );
    }
    free( mat->nodePtr );
}

