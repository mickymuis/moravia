#include "mst.h"
#include <stddef.h>

void 
mst_init( mst_t *m, int nnodes ) {
    m->n = m->m =nnodes;
    m->nedges =0;
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

MPI_Datatype
mst_mpiEdgeType() {
    static                       MPI_Datatype MPI_CEDGE = MPI_INT;
/*    static const int             blocklengths[4] = { 1, 1, 1, 1 };
    static const MPI_Datatype    types[4]        = { MPI_INT, MPI_INT, MPI_INT, MPI_FLOAT };
    static const MPI_Aint        offsets[4]      = { 
                                                    offsetof(cedge_t,from), 
                                                    offsetof(cedge_t,to), 
                                                    offsetof(cedge_t,hindex), 
                                                    offsetof(cedge_t,weight) };
*/
    if( MPI_CEDGE == MPI_INT ) {

//        MPI_Type_create_struct( 4, blocklengths, offsets, types, &MPI_CEDGE );
        MPI_Type_contiguous(4, MPI_INT, &MPI_CEDGE);
        MPI_Type_commit( &MPI_CEDGE );

    }
    return MPI_CEDGE;
}
