#include "graph.h"
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <assert.h>


typedef struct {
    idx_t* idxPtr;
    unsigned int capacity;
    unsigned int count;
} nodeset_t;

void
nodeset_make( nodeset_t* ns, unsigned int capacity ) {
    ns->idxPtr =calloc( capacity, sizeof(idx_t) );
    ns->capacity =capacity;
    ns->count =0;
}

void
nodeset_free( nodeset_t* ns ) { 
    free( ns->idxPtr );
}

void
nodeset_append( nodeset_t* ns, idx_t i ) {
    if( ns->capacity == ns->count ) {
        ns->idxPtr =realloc( ns->idxPtr, (ns->capacity+128) * sizeof(idx_t) );
        ns->capacity +=128;
    }
    ns->idxPtr[ns->count++] = i;
}

void
nodeset_print( nodeset_t* ns ) {
    for( int i =0; i < ns->count; i++ )
        printf( "(%d) ", ns->idxPtr[i] );
    printf( "\n" );
}

typedef struct {
    cedge_t *edgePtr;
    int capacity;
    int count;
} edgelist_t;

void
edgelist_make( edgelist_t* el, unsigned int capacity ) {
    el->edgePtr =calloc( capacity, sizeof(cedge_t) );
    el->capacity =capacity;
    el->count =0;
}

void
edgelist_free( edgelist_t* el ) { 
    free( el->edgePtr );
}

void
edgelist_clear( edgelist_t* el ) {
    el->count =0;
}

cedge_t*
edgelist_append( edgelist_t* el, cedge_t e ) {
    if( el->capacity == el->count ) {
        el->edgePtr =realloc( el->edgePtr, (el->capacity+128) * sizeof(cedge_t) );
        el->capacity +=128;
    }
    el->edgePtr[el->count++] = e;
    return &el->edgePtr[el->count-1];
}

cedge_t*
edgelist_geth( edgelist_t* el, idx_t hindex ) {
    for( int i =0; i < el->count; i++ ) {
        if( el->edgePtr[i].hindex == hindex ) {
            return &el->edgePtr[i];
        }
    }
    cedge_t e= { -1, -1, hindex, INFINITY };
    return edgelist_append( el, e );
}

/** Compute the contracted node from @n1 and @n2 by 
 *  merging their supernode pointer chains together.
 *  The level of all nodes in the chain is set to the highest level + 1
 */
int
contractNodes( mst_t* m, stree_t* t1, stree_t* t2, cedge_t* edge ) {

    if( !( t1->edges && t2->edges ) ) return -1;

    //printf( "Merging component [%d] with [%d]\n", t1->index, t2->index );
    if( t1->index == t2->index ) return -1;

    if( t1->capacity < (t1->n + t2->n + 1) ) {
        t1->edges = realloc( t1->edges, (t1->capacity + t2->n + 1) * sizeof(edge_t) );
        t1->capacity += t2->n + 1;
        assert( t1->edges != NULL );
    }

    if( t2->n ) {
        memcpy( &t1->edges[t1->n], t2->edges, t2->n * sizeof(edge_t) );
        t1->n += t2->n;
    }

    t1->edges[t1->n].from = edge->from;
    t1->edges[t1->n++].to = edge->to;
    t1->totalweight += t2->totalweight + edge->weight;

//#pragma omp parallel for
 /*   for( int i =0; i < m->m; i++ ) {
        if( m->hnode[i] == t2->index ) m->hnode[i] = t1->index;
      //  printf( "%d ", m->hnode[i] );
    }*/
   
    m->hnode[edge->from] = t1->index;
    m->hnode[edge->to] = t1->index;
    for( int i =0; i < t2->n; i++ ) {
        m->hnode[t2->edges[i].from] = t1->index;
        m->hnode[t2->edges[i].to] = t1->index;
    }
    //printf( "\n" );
    
    free( t2->edges );
    t2->edges =NULL;
    t2->capacity =0;
    t2->n =0;
    //printf( "." );
    return 0;
}

bool
computeBestEdge( graph_t* g, mst_t* m, node_t* n, cedge_t bestedge[] ) {
    
    //node_t *ptr =&g->nodePtr[node->first];
    idx_t hindex =m->hnode[n->index];
    cedge_t *be = &bestedge[hindex];
    if( be->from == -1 ) be->weight = INFINITY;

    for( int i =0; i < n->size; i++ ) {
        idx_t hindex2 =m->hnode[n->edges[i]];
        if( hindex2 == hindex ) continue;
        
        if( n->weights[i] < be->weight ) {
            // END
            be->weight = n->weights[i];
            be->from = n->index;
            be->to = n->edges[i];
        }

        cedge_t *be2 =&bestedge[hindex2];
        if( be2->from == -1 ) be2->weight = INFINITY;
        if( n->weights[i] < be2->weight ) {
            be2->weight = n->weights[i];
            be2->from = n->index;
            be2->to = n->edges[i];
        }
    
    }

    return true;
}

idx_t
selectStartNode( const graph_t* g, bool nodes[] ) {
    for( int i=0; i < g->m; i++ ) {
        if( nodes[i] ) {
            nodes[i] =false;
            return i;
        }
    }
    return -1;
}

void
makePartition( const graph_t* g, bool nodes[], nodeset_t* pre, nodeset_t* candidate, nodeset_t* post ) {
    nodeset_t* inset[2] = { pre, candidate };
    nodeset_t* outset[2]= { candidate, post };
    candidate->count = post->count =0;

    for( int s =0; s < 2; s++ ) {
        outset[s]->count =0;
        for( int i =0; i < inset[s]->count; i++ ) {
            node_t* n =&g->nodePtr[inset[s]->idxPtr[i]];
            for( int j =0; j < n->size; j++ ) {
                idx_t idx =n->edges[j];
                if( nodes[idx] == false ) continue;
                
                nodeset_append( outset[s], idx );
                nodes[idx] = false;
            }
        }
    }
}

int
boruvka3( graph_t* g, mst_t *m  ) {

    /* 
     * (1) Setup the partitions
     */
    const int bs =128;
    const int np =16;
    const int psize =g->m / 16;

    /*edgelist_t cedges[16];
    for( int i=0; i < np; i++ )
        edgelist_make( &cedges[i], bs );
*/
    /*
     * (2) Setup the mst_t data structure
     */
    m->n = m->m =g->m;
    m->hnode = calloc( m->m, sizeof(idx_t) );
    m->treePtr = calloc( m->m, sizeof(stree_t) );
    for( int i=0; i < m->m; i++ ) {
        m->hnode[i] =i;
        stree_t* t =&m->treePtr[i];
        t->index =i;
        t->n =0;
        t->capacity =1;
        t->edges = calloc( 1, sizeof(edge_t) );
    }

    /*
     * (3) Do stuff
     */

    cedge_t *bestedge = calloc( g->m, sizeof(cedge_t) );
    
    /*edgelist_t bestedge;
    edgelist_make( &bestedge, 128 );*/

    while( 1 ) {
        memset( bestedge, -1, g->m * sizeof(cedge_t) );
        edgelist_clear( &bestedge );
        
        printf( "Computing best edges...\n" );

        // Process the nodes from each partition by finding their minimum edges
//#pragma omp parallel for

    /*    for( int i =0; i < g->m; i++ ) {
            computeBestEdge( g, m, &g->nodePtr[i], bestedge );
        }*/

        for( int p =0; p < np; p++ ) {
            int ps =psize;
            if( p == 0 ) ps + g->m % np;
            for( int i= p*psize; i < p*psize + ps; i++ ) {
                computeBestEdge( g, m, &g->nodePtr[i], bestedge );
            }
        }

        int mops =0, calls =0;

        printf( "Merging nodes ... \n" );

        // Merge the nodes that were selected by the previous step
        for( int i =0; i < g->m; i++ ) {
            if( bestedge[i].from == -1 ) continue;
            // This edge is part of the minimum spanning tree
            cedge_t *e =&bestedge[i];
     //       printf( "Best edge for component [%d]: (%d)-(%d) with weight %f\n", m->hnode[e.from], e.from, e.to, e.weight );
            
            stree_t* t1 =&m->treePtr[m->hnode[e->from]];
            stree_t* t2 =&m->treePtr[m->hnode[e->to]];
            
            // Contract the source node with the destination node of the edge
            int l = contractNodes( m, t1, t2, e );
            if( l != -1 ) {
                mops++;
            }
            calls++;
        }
        printf( "Merged %d nodes with %d calls.\n", mops, calls );
        // No nodes were contracted?
        if( !mops ) break;
    }

    return 0;
}

#if 0
int
boruvka2( graph_t* g, mst_t *m  ) {

    /* 
     * (1) Setup the partitions
     */
    const int bs =128;
    int np =bs;  // Number of allocated partitions
    int ns =0;          // Actual number of partitions
    nodeset_t* separator = calloc( np, sizeof(nodeset_t) );
    nodeset_t* candidate = calloc( np, sizeof(nodeset_t) );
    for( int i=0; i < np; i++ ) {
        nodeset_make( &separator[i], bs );
        nodeset_make( &candidate[i], bs );
    }
    nodeset_t s0;
    nodeset_make( &s0, 1 );
    s0.count =1;

    bool *nodes = calloc( g->m, sizeof(bool) );
    memset( nodes, true, g->m*sizeof(bool) );
    
    // Create up to np partitions
    nodeset_t *pre, *cand, *post;
    
    // Pick the first node from which we build the set of partitions
PICK:
    s0.count =1;
    idx_t start =selectStartNode( g, nodes );
//    printf( "Picking %d\n", start );
    if( start != -1 ) {
        s0.idxPtr[0] =start;
        pre =&s0;


        while( 1 ) {
            cand = &candidate[ns];
            post = &separator[ns];

            makePartition( g, nodes, pre, cand, post );
            if( pre == &s0 ) {
                candidate[ns].idxPtr[candidate[ns].count++] =s0.idxPtr[0];
            }
            if( cand->count == 0 ) {
                goto PICK;
            }
            /*printf( "PARTITION #%d\n------------\n", ns );
            printf( "candidates " ); nodeset_print( cand );
            printf( "separator  " ); nodeset_print( post );*/
            ns++;

            if( ns == np ) {
                candidate =realloc( candidate, (np+bs) * sizeof(nodeset_t) );
                separator =realloc( separator, (np+bs) * sizeof(nodeset_t) );
                for( int i =ns; i < ns+bs; i++ ) {
                    nodeset_make( &candidate[i], bs );
                    nodeset_make( &separator[i], bs );
                }
                np+=bs;
            }
            pre =&separator[ns-1];
        }
    }

    assert( ns < np );

    // No partitions could be made for the given s0
    if( ns == 0 ) { 
        printf( "Could not make any partition for the given start node\n" );
        return;
    }
    printf( "Created %d partitions \n", ns );
    free( nodes );

    /*
     * (2) Setup the mst_t data structure
     */
    m->n = m->m =g->m;
    m->hnode = calloc( m->m, sizeof(idx_t) );
    m->treePtr = calloc( m->m, sizeof(stree_t) );
    for( int i=0; i < m->m; i++ ) {
        m->hnode[i] =i;
        stree_t* t =&m->treePtr[i];
        t->index =i;
        t->n =0;
        t->capacity =1;
        t->edges = calloc( 1, sizeof(edge_t) );
    }

    /*
     * (3) Do stuff
     */

    cedge_t *bestedge = calloc( g->m, sizeof(cedge_t) );
    
    while( 1 ) {
        memset( bestedge, -1, g->m * sizeof(cedge_t) );
        
        printf( "Computing best edges...\n" );

        // Process the nodes from each partition by finding their minimum edges
//#pragma omp parallel for
        for( int p =0; p < ns; p++ ) {
           // printf( "%d ", p );
            for( int i =0; i < candidate[p].count; i++ ) {
                node_t *n = &g->nodePtr[candidate[p].idxPtr[i]];
                
                computeBestEdge( g, m, n, bestedge );

               // printf( "Partition %d: (%d)->(%d)\n", p, from, to );
            }
            for( int i =0; i < separator[p].count; i++ ) {
                node_t *n = &g->nodePtr[separator[p].idxPtr[i]];
                
                computeBestEdge( g, m, n, bestedge );

               // printf( "Partition %d: (%d)->(%d)\n", p, from, to );
            }
        }

        /*for( int i =0; i < g->m; i++ ) {
            computeBestEdge( g, m, &g->nodePtr[i], bestedge );
        }*/

        int mops =0, calls =0;

        printf( "Merging nodes ... \n" );

        // Merge the nodes that were selected by the previous step
        for( int i =0; i < g->m; i++ ) {
            if( bestedge[i].from == -1 ) continue;
            // This edge is part of the minimum spanning tree
            cedge_t e =bestedge[i];
     //       printf( "Best edge for component [%d]: (%d)-(%d) with weight %f\n", m->hnode[e.from], e.from, e.to, e.weight );
            
            stree_t* t1 =&m->treePtr[m->hnode[e.from]];
            stree_t* t2 =&m->treePtr[m->hnode[e.to]];
            
            // Contract the source node with the destination node of the edge
            int l = contractNodes( m, t1, t2, &e );
            if( l != -1 ) {
                mops++;
            }
            calls++;
        }
        printf( "Merged %d nodes with %d calls.\n", mops, calls );
        // No nodes were contracted?
        if( !mops ) break;
    }

    for( int i=0; i < np; i++ ) {
        nodeset_free( &separator[i] );
        nodeset_free( &candidate[i] );
//        free( edges[i] );
    }
    nodeset_free( &s0 );
    free( candidate );
    free( separator );
    //free( nodes );
    free( bestedge );
    return 0;
}
#endif

int
main( int argc, const char** argv ) {

    if( argc != 2 ) {
        fprintf( stderr, "(i) Usage: %s [matrix market file]\n", argv[0] );
        return -1;
    }
    graph_t g;

    if( !graph_loadMM( argv[1], &g ) ) {
        fprintf( stderr, "(e) Could not load matrix market file `%s'\n", argv[1] );
        return -1;
    }

    mst_t mst;

    printf( "Computing MST ... " );
    int numedges =boruvka3( &g, &mst );
    printf( "done" );

   // graph_dump( &g );
  
    int sgraphs =0;

    for( int i=0; i < mst.m; i++ ) {
        stree_t *t =&mst.treePtr[i];
        if( t->edges != NULL ) {
                for( int j =0; j < t->n; j++ ) {
                }
                printf( "[%d] %d edges, total weight %f\n", ++sgraphs, t->n, t->totalweight );
        }
    }

    graph_free( &g );

    return 0;
}
