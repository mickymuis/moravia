

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

    g->m =M;
    g->nodePtr = calloc( M, sizeof( node_t ) );
    g->dataPtr = calloc( M, sizeof( nodedata_t ) );

    for( int i =0; i < M; i++ ) {
        node_t* n =&g->nodePtr[i];
        nodedata_t* d =&g->dataPtr[i];
        n->index =i;
        n->first =i; n->next =-1;
        n->level =0;
        // We use the fact that calloc lazily allocates this memory
        d->size =bs; d->degree =0;
        d->edges = calloc( bs, sizeof(idx_t) );
        d->weights=calloc( bs, sizeof(double) );
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

        bool swap =false;
        nodedata_t* d;
L0:
        if( row == col ) continue;

//        if( (swap ? col : row ) < (swap ? row : col) ) {
            d =&g->dataPtr[ swap ? col : row ];
            if( d->degree >= d->size ) {
                d->size += bs;
                d->edges =realloc( d->edges, d->size * sizeof(idx_t) ); 
                d->weights =realloc( d->weights, d->size * sizeof(double) ); 
            }
            d->edges[d->degree]   = swap ? row : col;
            d->weights[d->degree] = val;
            d->degree++;
            g->size++;
//        }

        if (mm_is_symmetric(matcode) && row != col && !swap ) {
            swap =true;
            goto L0;
        }
    }

    fclose(fh);
    for( int i =0; i < M; i++ ) {
        nodedata_t* d =&g->dataPtr[i];
        d->size =d->degree;
        d->edges = realloc( d->edges,   d->degree * sizeof(idx_t) );
        d->weights=realloc( d->weights, d->degree * sizeof(double) );
    }

    return true;
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

        //if( val == 0.0 ) continue;

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

/*static void
load_elements(const mmelement_t* elements, size_t n, size_t size, node_t* nodePtr[] )
{
    idx_t idx = 0;
    idx_t node = elements[0].row;
    

    idx_t *edgesTemp = calloc(n, sizeof(idx_t));
    double *weightsTemp = calloc(n, sizeof(double));

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
            nodePtr[node]->weights = calloc( idx, sizeof(double) );
            assert( nodePtr[node]->weights != NULL );
            nodePtr[node]->edges = calloc( idx, sizeof(idx_t) );
            assert( nodePtr[node]->edges != NULL );

            memcpy( nodePtr[node]->weights, weightsTemp, idx * sizeof(double) );
            memcpy( nodePtr[node]->edges, edgesTemp, idx * sizeof( idx_t) );
            nodePtr[node]->index =node;
            nodePtr[node]->first =node;
            nodePtr[node]->next  =-1;

            if( i+1 != size ) {
                assert( node < elements[i+1].row ); // List of elements is sorted by row
                node = elements[i+1].row;
            }
            idx =0;
        }
    }

    free( edgesTemp ); free( weightsTemp );
}*/

void
graph_dump( const graph_t* mat ) {

    for (int i = 0; i < mat->m; ++i)
    {
        if( mat->nodePtr[i].index != i ) continue; // contracted node
        nodedata_t* d =&mat->dataPtr[i];
        for (idx_t idx =0; idx < d->size; ++idx)
        {
            printf("(%d)->(%d) weight %f\n", i, d->edges[idx], d->weights[idx] );
        }
    }
}


bool
graph_loadMM(const char *filename, graph_t *mat ) {

    /*mmelement_t *elements = NULL;

    if (!read_matrix_market(filename, &elements, &mat->size, &mat->m ))
        return false;

    mat->nodePtr = calloc( mat->m, sizeof(node_t*) );

    load_elements(elements, mat->m, mat->size, mat->nodePtr);
    free( elements );*/

    if( !read_matrix_market2( filename, mat ) )
        return false;

    printf("(i) Import ok: graph with %d nodes, %d data points\n",
            mat->m, mat->size);

    return true;
}

void
graph_free( graph_t* mat ) {
    for( size_t i =0; i < mat->m; i++ ) {
        free( mat->dataPtr[i].weights );
        free( mat->dataPtr[i].edges );
    }
    free( mat->nodePtr );
    free( mat->dataPtr );
}

