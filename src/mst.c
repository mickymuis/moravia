#include "mst.h"

void 
mst_init( mst_t *m, int nnodes ) {
    m->n = m->m =nnodes;
    m->hindex = calloc( m->m, sizeof(idx_t) );
    m->edgePtr = calloc( m->m, sizeof(cedge_t) );
    m->treePtr = calloc( m->m, sizeof(stree_t) );
    for( int i=0; i < m->m; i++ ) {
        m->hindex[i] =i;
        stree_t* t =&m->treePtr[i];
        t->index =i;
        t->n =0;
        t->capacity =1;
        t->edges = calloc( 1, sizeof(edge_t) );
    }
}

void 
mst_free( mst_t *m ) {
    for( int i=0; i < m->m; i++ ) {
        stree_t* t =&m->treePtr[i];
        free( t->edges );
    }
    free( m->hindex );
    free( m->edgePtr );
    free( m->treePtr );
}

