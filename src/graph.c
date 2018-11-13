

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include "graph.h"

#include "mmio.h"


/*
 * Load graph market file
 */

typedef struct {
    size_t row, col;
    float val;
} mmelement_t;

int
compareMMElement( const void* a, const void* b ) {
    mmelement_t* elemA = (mmelement_t*)a;
    mmelement_t* elemB = (mmelement_t*)b;

    int i = elemA->row - elemB->row;

    return i ? i : elemA->col - elemB->col;
}

static bool
read_matrix_market(const char *filename,
        mmelement_t **elements_ptr, uint32_t* size_ptr,
        uint32_t* n_nodes_ptr )
{
    FILE *fh = fopen(filename, "r");
    if (!fh)
    {
        perror("fopen");
        return false;
    }

    printf("Reading graph '%s'...\n", filename);

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

    *n_nodes_ptr = M;

    // Reserve space for all elements

    const size_t nmemb = mm_is_symmetric(matcode) ? 2 * size : size;
    mmelement_t* elements = calloc( nmemb, sizeof( mmelement_t ) );
    size_t off =0;

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

        mmelement_t e = { row, col, val };
        elements[off++] = e;
        if (mm_is_symmetric(matcode) && row != col) {
            e.row = col; e.col =row;
            elements[off++] = e;
        }
    }

    fclose(fh);
    qsort( elements, off, sizeof( mmelement_t ), compareMMElement );

    *elements_ptr = elements;
    *size_ptr =off;

    return true;
}

/*
 * Transfer graph elements into Compressed Row Storage structure
 */

static void
load_elements(const mmelement_t* elements, size_t n, size_t size, node_t* nodePtr[] )
{
    idx_t idx = 0;
    idx_t node = elements[0].row;
    

    idx_t edgesTemp[n];
    float  weightsTemp[n];

    for( size_t i =0; i < size; i++ ) {

        const mmelement_t *it = &elements[i];

        weightsTemp[idx] = it->val;
        edgesTemp[idx] = it->col;
        idx++;
        
        // We have complete one node, write it to the datastructure
        if ( i+1 == size || elements[i+1].row != node)
        {
            nodePtr[node] = calloc( 1, sizeof(node_t) );
            assert( nodePtr[node] != NULL );

            nodePtr[node]->size = idx;
            nodePtr[node]->degree =idx;
            nodePtr[node]->weights = calloc( idx, sizeof(float) );
            assert( nodePtr[node]->weights != NULL );
            nodePtr[node]->edges = calloc( idx, sizeof(idx_t) );
            assert( nodePtr[node]->edges != NULL );

            memcpy( nodePtr[node]->weights, weightsTemp, idx * sizeof(float) );
            memcpy( nodePtr[node]->edges, edgesTemp, idx * sizeof( idx_t) );
            nodePtr[node]->index =node;
            nodePtr[node]->first =node;
            nodePtr[node]->next  =-1;

            if( i+1 != size )
                node = elements[i+1].row;
            idx =0;
        }
    }
}

void
graph_dump( const graph_t* mat ) {

    for (int i = 0; i < mat->m; ++i)
    {
        if( mat->nodePtr[i] == NULL ) continue;
        if( mat->nodePtr[i]->index != i ) continue; // contracted node
        for (idx_t idx =0; idx < mat->nodePtr[i]->size; ++idx)
        {
            printf("(%d)->(%d) weight %f\n", i, mat->nodePtr[i]->edges[idx], mat->nodePtr[i]->weights[idx] );
        }
    }
}


bool
graph_loadMM(const char *filename, graph_t *mat ) {

    mmelement_t *elements = NULL;

    if (!read_matrix_market(filename, &elements, &mat->size, &mat->m ))
        return false;

    mat->nodePtr = calloc( mat->m, sizeof(node_t*) );

    load_elements(elements, mat->m, mat->size, mat->nodePtr);
    free( elements );

    printf("(i) Import ok: graph with %d nodes, %d data points\n",
            mat->m, mat->size);

    return true;
}

void
graph_free( graph_t* mat ) {
    for( size_t i =0; i < mat->m; i++ ) {
        if( mat->nodePtr[i] != NULL && mat->nodePtr[i]->index == i ) {
            free( mat->nodePtr[i]->weights );
            free( mat->nodePtr[i]->edges );
            free( mat->nodePtr[i] );
        }
    }
    free( mat->nodePtr );
}

