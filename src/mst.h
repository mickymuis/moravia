#ifndef __MST_H__
#define __MST_H__

#include "graph.h"
#include <mpi.h>

typedef struct {
    idx_t from, to, hindex;
    float weight;
} cedge_t;

typedef struct {
    idx_t    index;       // Index of the hypergraph
    edge_t*  edges;       // Pointer to an array of edges
    double   totalweight; // Total weight of the hypergraph
    uint32_t n;           // Number of edges
    uint32_t capacity;    // Capacity of the array
} stree_t;

typedef struct {
    stree_t* treePtr;     // Pointer to an array of mst nodes
    uint32_t n;           // Total number of subtrees
    uint32_t m;           // Number of nodes in the graph
    idx_t* hindex;        // Hypernode index for each graph node
    cedge_t* edgePtr;     // Buffer that stores all currently processed edges
    uint32_t nedges;      // Number of edges so far
    uint32_t capacity;    // Capacity of the arrays
} mst_t;

void mst_init( mst_t *m, int nnodes );

void mst_free( mst_t *m );

/** Given graph node index @nidx, compute the hypergraph index it belongs to.
 */
inline static idx_t
mst_computeHIndex( const mst_t *m, idx_t nidx ) {
    idx_t i =nidx;
    while( m->hindex[i] != i )
        i = m->hindex[i];
    return i;
}

/** Takes the hypergraph index from @m and follows each redirection and 
 *  replaces each element with its actual index, such that 
 *  all indices can be computed in constant time. */
inline static void
mst_flattenHIndex( mst_t *m ) {
    idx_t *queue = calloc( m->m, sizeof(idx_t) );
    bool *flat = calloc( m->m, sizeof(bool) );
    idx_t n=0;

    for( idx_t i=0; i < m->m; i++ ) {
        if( flat[i] ) continue;
        n =0;
        idx_t j =i;
        while( j != m->hindex[j] ) {
            flat[j] =true;
            queue[n++] =j;
            j = m->hindex[j];
        }
        for( idx_t k =0; k < n; k++ ) {
            m->hindex[queue[k]] =j;
        }
    }
    free( flat );
    free( queue );
}

inline static void
mst_updateHIndex( mst_t *m, const cedge_t *edgePtr, int n ) {
    for( int i =0; i < n; i++ ) {
        cedge_t *edge = &edgePtr[i];
        m->hindex[edge->from] = edge->hindex;
        m->hindex[edge->to]   = edge->hindex;
        /*m->hindex[mst_computeHIndex(m,edge->from)] = edge->hindex;
        m->hindex[mst_computeHIndex(m,edge->to)]   = edge->hindex;*/
    }
    mst_flattenHIndex( m );
}

/** Returns the MPI datatype to cover the cedge_t type used by the mst_t structure
 */
MPI_Datatype
mst_mpiEdgeType();

#endif
