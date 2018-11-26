#include "graph.h"
#include "mst.h"
#include "nodeset.h"
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <stdarg.h>
#include <mpi.h>

void
rprintf( const char* fmt, ... ) {
    va_list ap;
    va_start(ap,fmt);
    int pid;
    MPI_Comm_rank( MPI_COMM_WORLD, &pid );

    if( pid == 0 )
        vprintf( fmt, ap ); 
    va_end(ap);
}

/* Functions related to sending and received MST parts over MPI */

void
sendEdge( const cedge_t* edge, int pid, int dst_pid ) {

}

/* Core functions for MST computation */

/** Compute the contracted subtree from @t1 and @t2 by adding @edge to the list of MST edges. 
 */
int
contractNodes( mst_t* m, stree_t* t1, stree_t* t2, cedge_t* edge ) {

    if( t1->index == t2->index ) return -1;

    edge->hindex =t1->index;
    m->edgePtr[m->nedges++] =*edge;
   
    m->hindex[edge->from] = t1->index;
    m->hindex[edge->to] = t1->index;
    m->hindex[t2->index] = t1->index;
    return 0;
}

/** Given graph @g, mst @m and node @n, computes the shortest edge from @n that isn't a self edge.
 *  The resulting shortest edge (if any) is stored in @bestedge[i] where i is @n's hypergraph index.*/
bool
computeBestEdge( const graph_t* g, const mst_t* m, const node_t* n, cedge_t bestedge[] ) {
    
    //node_t *ptr =&g->nodePtr[node->first];
    //idx_t hindex =m->hindex[n->index];
    idx_t hindex =mst_computeHIndex( m, n->index );
    cedge_t *be = &bestedge[hindex];
    if( be->from == -2 ) return true;
    if( be->from == -1 ) be->weight = INFINITY;

    for( int i =0; i < n->size; i++ ) {
        //idx_t hindex2 =m->hindex[n->edges[i]];
        idx_t hindex2 =mst_computeHIndex( m, n->edges[i] );
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

/** Return the first node index i in @g for which @nodes[i] is true */ 
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

/** For a set of nodes @in, add all direct descendants of @in for which @nodes[] is true, to @out. */
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
computeMST( const graph_t* g, mst_t* m, nodeset_t* partition, int nparts, bool root, int pid, int dst_pid ) {
    
    cedge_t *bestedge = calloc( g->m, sizeof(cedge_t) );
    memset( bestedge, -1, g->m * sizeof(cedge_t) );
    int hbegin =0, hend = g->m;
    // Given the partitions which we need to process, 
    // we need to exclude certain components from being computed
    if( !root ) {
        for( int i=0; i < partition[0].count; i++ ) {
            bestedge[m->hindex[partition[0].idxPtr[i]]].from =-2;
            hbegin =fminf( hbegin, m->hindex[partition[0].idxPtr[i]] );
        }
        for( int i=0; i < partition[nparts-1].count; i++ ) {
            bestedge[m->hindex[partition[nparts-1].idxPtr[i]]].from =-2;
            hend =fmaxf( hend, m->hindex[partition[nparts-1].idxPtr[i]] );
        }
    }

    while( 1 ) {
       
        //if( all ) printf( "Computing best edges...\n" );

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

        //if( all ) printf( "Merging nodes ... \n" );

        // Merge the nodes that were selected by the previous step
        for( int i =hbegin; i < hend; i++ ) {
            if( bestedge[i].from < 0 ) continue;
            // This edge is part of the minimum spanning tree
            cedge_t *e =&bestedge[i];

            /*idx_t h1 = m->hindex[e->from];
            idx_t h2 = m->hindex[e->to];*/
            idx_t h1 = mst_computeHIndex( m, e->from );
            idx_t h2 = mst_computeHIndex( m, e->to );

            stree_t* t1 =&m->treePtr[h1];
            stree_t* t2 =&m->treePtr[h2];
            
            // Contract the source node with the destination node of the edge
            int l = contractNodes( m, t1, t2, e );
            if( l != -1 ) {
                mops++;
            }

            //e->hindex =h1;
            sendEdge( e, pid, dst_pid );

            if( bestedge[h2].from == -2 ) {
                bestedge[h1].from = -2;
                e->from =-2;
            }
            else e->from =-1;
            calls++;
        }
        /*if( all )*/ //printf( "Merged %d nodes with %d calls.\n", mops, calls );
        
        mst_flattenHIndex( m );

        // No nodes were contracted?
        if( !mops ) break;
    }
    free( bestedge );

}

/** If @m only contains a list of raw edges @m->edgePtr, compute all MST subtrees by iterating over @m->edgePtr
 *  and adding the edges to their corresponding hypergraph in @m->treePtr
 */
void
postprocessMST( mst_t* m ) {
    // Optimze the indices to save some time
    mst_flattenHIndex( m );

    // Iterate over all raw edges
    for( int i=0; i < m->nedges; i++ ) {
        cedge_t* edge =&m->edgePtr[i];
        idx_t hindex =m->hindex[edge->hindex];
        // Hypergraph subtree which this edge belongs to
        stree_t *t =&m->treePtr[hindex];
        // Increse the capacity of the array if neccesary
        if( t->capacity < t->n + 1 ) {
            t->edges = realloc( t->edges, (t->capacity + 128) * sizeof(edge_t) );
            t->capacity += 128;
            assert( t->edges != NULL );
        }

        // Add the edge to the subtree
        edge_t e = { edge->from, edge->to };
        t->edges[t->n++] =e;
        t->totalweight +=edge->weight;
    }
}

int
boruvka5( graph_t* g, mst_t *m, int procs, int pid  ) {

    /* 
     * (1) Setup the partitions
     */
    //const int procs =16;// Number of processes
    const int bs =1024; // Block size for pre-allocating
    int cparts =bs;     // Number of allocated partitions
    int nparts =0;      // Actual number of partitions
    nodeset_t* partition = calloc( cparts, sizeof(nodeset_t) );
    for( int i=0; i < cparts; i++ ) {
        nodeset_init( &partition[i], bs );
    }
    nodeset_t s0;
    nodeset_init( &s0, 1 );
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
                    nodeset_init( &partition[i], bs );
                }
                cparts+=bs;
            }
            in =&partition[nparts-1];
        }
    }

    // No partitions could be made for the given s0
    if( nparts == 0 ) { 
        rprintf( "Could not make any partition for the given start node\n" );
        return -1;
    }
    rprintf( "Created %d partitions \n", nparts );
    free( nodes );

    /*
     * (2) Setup the mst_t data structure
     */

    mst_init( m, g->m );

    /*
     * (3) Distribute the partitions to different nodes
     */

    int nparts_div = (nparts/procs) * procs;
    for( int np =procs; np > 0; np = np/2 ) {
    int ppn = nparts_div / np; // Part per node
        for( int p=0; p < np; p++ ) {
            int begin = p*ppn;
            int end = (p+1)*ppn;
            if( p==np-1 ) end = nparts;
            int size = end - begin;
            if( size<=0 ) continue;
            
            int src_pid = p * procs/np;
            int dst_pid = (p/2) * procs/np * 2;
            
            rprintf( "[level %d] Rank %d, sending to %d. Computing parts %d - %d...\n", 
                    procs/np, src_pid, dst_pid, begin, p*ppn + size );
            
            if( pid == src_pid )
                computeMST( g, m, &partition[begin], size, np==1?true:false, pid, dst_pid );
        }
    }

    /*
     * (4) Finish the MST by processing the raw edge list into subtrees
     */


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
main( int argc, char** argv ) {

    int procs, pid;

    MPI_Init( &argc, &argv );

    MPI_Comm_rank( MPI_COMM_WORLD, &pid );
    MPI_Comm_size( MPI_COMM_WORLD, &procs );

    if( argc != 2 ) {
        if( pid == 0 )
            fprintf( stderr, "(i) Usage: %s [matrix market file]\n", argv[0] );
        MPI_Finalize();
        return -1;
    }
    graph_t g;

    double start_time, end_time;
    start_time = MPI_Wtime();

    if( !graph_loadMM( argv[1], &g ) ) {
        fprintf( stderr, "(e) [%d] Could not load matrix market file `%s'\n", pid, argv[1] );
        MPI_Finalize();
        return -1;
    }

    MPI_Barrier( MPI_COMM_WORLD );
    
    rprintf("(i) Import ok: graph with %d nodes, %d data points\n", g.m, g.size);

    end_time = MPI_Wtime();

    rprintf( "(t) Loaded matrix market in %gs\n", end_time-start_time );

    mst_t mst;

    start_time = MPI_Wtime();

    rprintf( "Computing MST ... " );
    boruvka5( &g, &mst, procs, pid );
    
    end_time = MPI_Wtime();

    rprintf( "(t) Computed MST in %gs\n", end_time-start_time );
    
    postprocessMST( &mst );

   // graph_dump( &g );
  
    if( pid == 0 ) {
        int sgraphs =0;

        for( int i=0; i < mst.m; i++ ) {
            stree_t *t =&mst.treePtr[i];
            if( t->n != 0 ) {
                    for( int j =0; j < t->n; j++ ) {
                    }
                    printf( "[%d] %d edges, total weight %f\n", ++sgraphs, t->n, t->totalweight );
            }
        }
    }

    graph_free( &g );
    mst_free( &mst );

    MPI_Finalize();
    return 0;
}
