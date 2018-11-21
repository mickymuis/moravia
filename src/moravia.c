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
        t1->edges = realloc( t1->edges, (t1->capacity + t2->n + 128) * sizeof(edge_t) );
        t1->capacity += t2->n + 128;
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
//#pragma omp parallel for
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
computeBestEdge( const graph_t* g, const mst_t* m, const node_t* n, cedge_t bestedge[] ) {
    
    //node_t *ptr =&g->nodePtr[node->first];
    idx_t hindex =m->hnode[n->index];
    cedge_t *be = &bestedge[hindex];
    if( be->from == -2 ) return true;
    if( be->from == -1 ) be->weight = INFINITY;

    for( int i =0; i < n->size; i++ ) {
        idx_t hindex2 =m->hnode[n->edges[i]];
        if( hindex2 == hindex ) continue;
        
        if( n->weights[i] < be->weight ) {
            be->weight = n->weights[i];
            be->from = n->index;
            be->to = n->edges[i];
        }

        cedge_t *be2 =&bestedge[hindex2];
        if( be2->from == -2 ) continue;
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

inline void
makePartition( const graph_t* g, bool nodes[], nodeset_t* in, nodeset_t* out ) {
    out->count =0;

    for( int i =0; i < in->count; i++ ) {
        node_t* n =&g->nodePtr[in->idxPtr[i]];
        for( int j =0; j < n->size; j++ ) {
            idx_t idx =n->edges[j];
            if( nodes[idx] == false ) continue;
            
            nodeset_append( out, idx );
            nodes[idx] = false;
        }
    }
}

void
computeMST( const graph_t* g, mst_t* m, nodeset_t* partition, int nparts, bool all ) {
    
    cedge_t *bestedge = calloc( g->m, sizeof(cedge_t) );
    memset( bestedge, -1, g->m * sizeof(cedge_t) );
    int hbegin =0, hend = g->m;
    // Given the partitions which we need to process, 
    // we need to exclude certain components from being computed
    if( !all ) {
        for( int i=0; i < partition[0].count; i++ ) {
            bestedge[m->hnode[partition[0].idxPtr[i]]].from =-2;
            hbegin =fminf( hbegin, m->hnode[partition[0].idxPtr[i]] );
        }
        for( int i=0; i < partition[nparts-1].count; i++ ) {
            bestedge[m->hnode[partition[nparts-1].idxPtr[i]]].from =-2;
            hend =fmaxf( hend, m->hnode[partition[nparts-1].idxPtr[i]] );
        }
    }

    while( 1 ) {
       
        if( all ) printf( "Computing best edges...\n" );

        // Process the nodes from each partition by finding their minimum edges
        for( int p =0; p < nparts; p++ ) {
           // printf( "%d ", p );
            for( int i =0; i < partition[p].count ; i++ ) {
                node_t *n = &g->nodePtr[partition[p].idxPtr[i]];
                
                computeBestEdge( g, m, n, bestedge );

               // printf( "Partition %d: (%d)->(%d)\n", p, from, to );
            }
        }

        int mops =0, calls =0;

        if( all ) printf( "Merging nodes ... \n" );

        // Merge the nodes that were selected by the previous step
        for( int i =hbegin; i < hend; i++ ) {
            if( bestedge[i].from < 0 ) continue;
            // This edge is part of the minimum spanning tree
            cedge_t *e =&bestedge[i];

            idx_t h1 = m->hnode[e->from];
            idx_t h2 = m->hnode[e->to];

            stree_t* t1 =&m->treePtr[h1];
            stree_t* t2 =&m->treePtr[h2];
            
            // Contract the source node with the destination node of the edge
            int l = contractNodes( m, t1, t2, e );
            if( l != -1 ) {
                mops++;
            }
            if( bestedge[h2].from == -2 ) {
                bestedge[h1].from = -2;
                e->from =-2;
            }
            else e->from =-1;
            calls++;
        }
        /*if( all )*/ printf( "Merged %d nodes with %d calls.\n", mops, calls );
        // No nodes were contracted?
        if( !mops ) break;
    }
    free( bestedge );

}

int
boruvka4( graph_t* g, mst_t *m  ) {

    /* 
     * (1) Setup the partitions
     */
    const int procs =16;// Number of processes
    const int bs =1024; // Block size for pre-allocating
    int cparts =bs;     // Number of allocated partitions
    int nparts =0;      // Actual number of partitions
    nodeset_t* partition = calloc( cparts, sizeof(nodeset_t) );
    for( int i=0; i < cparts; i++ ) {
        nodeset_make( &partition[i], bs );
    }
    nodeset_t s0;
    nodeset_make( &s0, 1 );
    s0.count =1;

    bool *nodes = calloc( g->m, sizeof(bool) );
    memset( nodes, true, g->m*sizeof(bool) );
    
    // Create up to np partitions
    nodeset_t *in, *out;
    
    // Pick the first node from which we build the set of partitions
PICK:
    s0.count =1;
    idx_t start =selectStartNode( g, nodes );
//    printf( "Picking %d\n", start );
    if( start != -1 ) {
        s0.idxPtr[0] =start;
        in =&s0;

        while( 1 ) {
            out = &partition[nparts];

            makePartition( g, nodes, in, out );
            if( in == &s0 ) {
                nodeset_append( out, s0.idxPtr[0] );
            }
            if( out->count == 0 ) {
                goto PICK;
            }
            /*printf( "PARTITION #%d\n------------\n", nparts );
            nodeset_print( out );*/
            nparts++;

            if( nparts == cparts ) {
                partition =realloc( partition, (cparts+bs) * sizeof(nodeset_t) );
                for( int i =nparts; i < nparts+bs; i++ ) {
                    nodeset_make( &partition[i], bs );
                }
                cparts+=bs;
            }
            in =&partition[nparts-1];
        }
    }

    // No partitions could be made for the given s0
    if( nparts == 0 ) { 
        printf( "Could not make any partition for the given start node\n" );
        return -1;
    }
    printf( "Created %d partitions \n", nparts );
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
     * (3) Distribute the partitions to different nodes
     */

    int nparts_div = (nparts/procs) * procs;
    for( int np =procs; np > 0; np = np/2 ) {
    int ppn = nparts_div / np; // Part per node
//#pragma omp parallel for
        for( int p=0; p < np; p++ ) {
            int begin = p*ppn;
            int end = (p+1)*ppn;
            if( p==np-1 ) end = nparts;
            int n = end - begin;
            if( n<=0 ) continue;
            printf( "[np=%d] Computing parts %d - %d...\n", np, begin, p*ppn + n );
            computeMST( g, m, &partition[begin], n, np==1?true:false );
        }
    }

    /*
     * (4) Finish the MST on a single node
     */

    //printf( "Computing from all parts...\n" );
    //computeMST( g, m, partition, nparts, true );

    /*
     * (5) Clean up
     */
    for( int i=0; i < cparts; i++ ) {
        nodeset_free( &partition[i] );
    }
    nodeset_free( &s0 );
    free( partition );
    return 0;
}


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
    int numedges =boruvka4( &g, &mst );
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
