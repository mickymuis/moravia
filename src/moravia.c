#include "graph.h"
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <assert.h>

/* The edge_t type is only used to store the edges of the resulting MSP */
typedef struct {
    idx_t from, to;
    double weight;
} edge_t;

static const edge_t NOEDGE = { -1, -1, INFINITY };

void
printSuperNodes( graph_t* g ) {

    for( idx_t i=0; i < g->m; i++ ) {
        printf( "%d: ", i );
        node_t* n =&g->nodePtr[g->nodePtr[i].first];
        while( 1 ) {
            printf( "%d -> ", n->index );  
            if( n->next == -1 ) break;
            n = &g->nodePtr[n->next];
        }

        printf( "X\n" );
    }
}

/** Compute the contracted node from @n1 and @n2 by 
 *  merging their supernode pointer chains together.
 *  The level of all nodes in the chain is set to the highest level + 1
 */
int
contractNodes( graph_t* g, node_t* n1, node_t* n2 ) {

    if( n1->first == n2->first ) return -1;

    node_t* n1_last = &g->nodePtr[n1->first];
    int level = (int)fmaxf( n1->level, n2->level) + 1;
    
    // Find the last node in n1's chain and update all nodes' level
    while( 1 ) {
        n1_last->level = level;
        if( n1_last->next == -1 ) break;
        n1_last =&g->nodePtr[n1_last->next];
    }
     
    n1_last->next = n2->first;
    node_t* n =&g->nodePtr[n2->first];

    // Also update the levels and first pointers in n2's chain
    while( 1 ) {
        n->level =level;
        n->first =n1->first;
        if( n->next == -1 ) break;
        n = &g->nodePtr[n->next];
    }

    return level;
}

bool
computeBestEdge( graph_t* g, node_t* node, edge_t *bestedge ) {
    
    node_t *ptr =&g->nodePtr[node->first];
    edge_t *be = &bestedge[node->first];
    if( be->from == -1 ) be->weight = INFINITY;

    while( 1 ) {
        nodedata_t *d =&g->dataPtr[ptr->index];

        for( idx_t i =0; i < d->size; i++ ) {
            //if( ptr->weights[i] == 0.f ) continue; // Skip explicit zeroes
            if( d->weights[i] < be->weight ) {
                if( g->nodePtr[d->edges[i]].first == node->first ) continue; // Selfedge

                // END
                be->weight = d->weights[i];
                be->from = ptr->index;
                be->to = d->edges[i];
            }
        }
        if( ptr->next == -1 ) break;
        ptr =&g->nodePtr[ptr->next];
    }

    return true;
}

void
removeSuperNode( graph_t* g, bool nodes[], node_t* n ) {
    if( !n ) return;

    node_t* ptr =&g->nodePtr[n->first];
    nodes[ptr->index] = false;

    // Remove one or more nodes depending on the size of the super node n belongs to
    /*while( 1 ){
        nodes[ptr->index] = NULL;
        if( ptr->next == -1 ) break;
        ptr = g->nodePtr[ptr->next];
    }*/
}

void
addSuperNode( graph_t* g, bool nodes[], node_t* n ) {
    if( !n ) return;

    node_t* ptr =&g->nodePtr[n->first];

    while( 1 ){
        nodes[ptr->index] = true;
        if( ptr->next == -1 ) break;
        ptr = &g->nodePtr[ptr->next];
    }
}

typedef struct {
    node_t** nodePtr;
    unsigned int capacity;
    unsigned int count;
} nodeset_t;

void
nodeset_make( nodeset_t* ns, unsigned int capacity ) {
    ns->nodePtr =calloc( capacity, sizeof(node_t) );
    ns->capacity =capacity;
    ns->count =0;
}

void
nodeset_free( nodeset_t* ns ) { 
    free( ns->nodePtr );
}

void
nodeset_print( nodeset_t* ns ) {
    for( int i =0; i < ns->count; i++ )
        printf( "(%d) ", ns->nodePtr[i]->index );
    printf( "\n" );
}

node_t*
selectStartNode( const graph_t* g, bool nodes[], int level ) {

    for( int i =0; i < g->m; i++ ) {
        idx_t idx=g->nodePtr[i].first;
        if( nodes[idx] != false /*&& g->nodePtr[i].level == level*/ ) {
            node_t* n=&g->nodePtr[idx];
            //if( nodes[n->first] == false ) continue; // supernode already selected
            nodes[idx] = false;
            //nodes[n->first] = NULL;
            //removeSuperNode( g, nodes, &g->nodePtr[i] );
            return n;
        }
    }
    return NULL;
}

void
makePartition( const graph_t* g, bool nodes[], nodeset_t* pre, nodeset_t* candidate, nodeset_t* post ) {
    nodeset_t* inset[2] = { pre, candidate };
    nodeset_t* outset[2]= { candidate, post };
    candidate->count = post->count =0;

    for( int s =0; s < 2; s++ ) {
        outset[s]->count =0;
        for( int i =0; i < inset[s]->count; i++ ) {
            idx_t first =inset[s]->nodePtr[i]->first;
            if( first != inset[s]->nodePtr[i]->index ) continue; //?
            node_t* n =&g->nodePtr[first];
            while(1) {
                nodedata_t* d =&g->dataPtr[n->index];
                for( int j =0; j < d->size; j++ ) {
                    idx_t idx =g->nodePtr[d->edges[j]].first;
                    if( nodes[idx] == false ) continue;
                    //if( nodes[g->nodePtr[idx].first] == false ) continue;
                    //if( n->weights[j] != 0.f ) // Skip explicit zeroes
                    
                    outset[s]->nodePtr[outset[s]->count++] = &g->nodePtr[idx];
                    nodes[idx] = false;
                    //removeSuperNode( g, nodes, &g->nodePtr[idx] );
                }
                if( n->next == -1 ) break;
                n =&g->nodePtr[n->next];
            }
        }
    }
}

/*void 
makePartition2( const graph_t* g, node_t* nodes[], nodeset_t* pre, nodeset_t* candidate, nodeset_t* post ) {
    candidate->count =0;
    post->count =0;

    for( int i =0; i < pre->count; i++ ) {
        idx_t first =pre->nodePtr[i]->first;
        //if( nodes[first] == NULL ) continue;
        node_t* n =g->nodePtr[first];
        while( 1 ) {
            for( int j =0; j < n->size; j++ ) {
                idx_t idx =n->edges[j];
                if( nodes[idx] == NULL ) continue;
                node_t *n2 = nodes[g->nodePtr[idx]->first];
                if( n2 == NULL ) continue;
                removeSuperNode( g, nodes, nodes[idx] );

                //if( n->weights[j] != 0.f ) // Skip explicit zeroes
                candidate->nodePtr[candidate->count++] = g->nodePtr[idx];

                while( 1 ) {
                    for( int k =0; k < n2->size; k++ ) {
                        idx_t idx2 =n2->edges[k];
                        if( nodes[idx2] == NULL ) continue;
                        if( nodes[g->nodePtr[idx2]->first] == NULL ) continue;
                        removeSuperNode( g, nodes, nodes[idx2] );

                        //if( n->weights[j] != 0.f ) // Skip explicit zeroes
                        post->nodePtr[post->count++] = g->nodePtr[idx2];
                    }

                    if( n2->next == -1 ) break;
                    n2 =g->nodePtr[n2->next];
                }
            }
            if( n->next == -1 ) break;
            n =g->nodePtr[n->next];
        }
    }
}*/

int
boruvka2( graph_t* g, edge_t *edgelist  ) {
    const int np =128;

    nodeset_t separator[np];
    nodeset_t candidate[np];
    //edge_t *edges[np];
    for( int i=0; i < np; i++ ) {
        nodeset_make( &separator[i], g->m );
        nodeset_make( &candidate[i], g->m );
      //  edges[i] = calloc( g->m, sizeof(edge_t) );
    }
    nodeset_t s0;
    nodeset_make( &s0, 1 );
    s0.count =1;

    bool *nodes = calloc( g->m, sizeof(bool) );

    edge_t *bestedge = calloc( g->m, sizeof(edge_t) );

    int level =0, maxlevel =0, numedges =0, ns =0;
    
    while( 1 ) {
        memset( nodes, true, g->m*sizeof(bool) );
        memset( bestedge, -1, g->m * sizeof(edge_t) );
PICK:
        // Pick the first node from which we build the set of partitions
        s0.count =1;
        s0.nodePtr[0] =selectStartNode( g, nodes, level );

        // There are no nodes of this level left
        if( s0.nodePtr[0] == NULL && !ns ) {
           /* printf( "No nodes left in current level (%d)\n", level );
            memset( nodes, true, g->m*sizeof(bool) );
            level++;
            continue;*/
            break;
        }

        // Create up to np partitions
        nodeset_t *pre, *cand, *post;
        pre =&s0;

        if( s0.nodePtr[0] != NULL ) {
            //printf( "Picked s0=(%d)\n", s0.nodePtr[0]->index );
            while( ns < np ) {
                cand = &candidate[ns];
                post = &separator[ns];

                makePartition( g, nodes, pre, cand, post );
                if( pre == &s0 ) {
                    candidate[ns].nodePtr[candidate[ns].count++] =s0.nodePtr[0];
                }
                if( cand->count == 0 ) {
                    goto PICK; // Pick another start node
                }
/*                printf( "PARTITION #%d\n------------\n", ns );
                printf( "candidates " ); nodeset_print( cand );
                printf( "separator  " ); nodeset_print( post );*/
                pre =&separator[ns];
                ns++;
            }
        }

        // No partitions could be made for the given s0
        if( ns == 0 ) { 
            printf( "Could not make any partition for the given start node\n" );
            continue;
        }

        printf( "Created %d partitions, computing minmum edges... ", ns );
        

        // Process the nodes from each partition by finding their minimum edges
#pragma omp parallel for
        for( int p =0; p < ns; p++ ) {
           // printf( "%d ", p );
            for( int i =0; i < candidate[p].count; i++ ) {
                node_t *n = candidate[p].nodePtr[i];
                
                computeBestEdge( g, n, bestedge );

               // printf( "Partition %d: (%d)->(%d)\n", p, from, to );
              //  printf( "(%d)-(%d) with weight %f\n", e.from, e.to, e.weight );
            }
        }

        int mops =0;

        printf( "done\nMerging nodes ... " );

        // Merge the nodes that were selected by the previous step
        for( int i =0; i < g->m; i++ ) {
            if( bestedge[i].from == -1 ) continue;
            // This edge is part of the minimum spanning tree
            edge_t e =bestedge[i];
            node_t* dnode =&g->nodePtr[e.to];
            node_t* snode =&g->nodePtr[e.from];
            
            // Contract the source node with the destination node of the edge
            int l = contractNodes( g, snode, dnode );
            if( l != -1 ) {
                maxlevel =fmaxf( maxlevel, l );
                edgelist[numedges++] =e;
                mops++;
                //node_t* first_node =&g->nodePtr[snode->first];
                //addSuperNode( g, nodes, snode );
                //candidate[first_node->index] =first_node;
            }
        }
        printf( " merged %d nodes.\n", mops );
        ns =0;
        // No nodes were contracted?
        if( !mops ) break;
    }

    for( int i=0; i < np; i++ ) {
        nodeset_free( &separator[i] );
        nodeset_free( &candidate[i] );
//        free( edges[i] );
    }
    nodeset_free( &s0 );
    free( nodes );
    free( bestedge );
    return numedges;
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

    edge_t *edgelist = calloc( g.m, sizeof(edge_t) );
  
    printf( "Computing MST ... " );
    int numedges =boruvka2( &g, edgelist );
    printf( "done" );

//    int numedges = boruvka( &g, (idx_t *)edgelist );
   // graph_dump( &g );
  
    // We now need to spit the raw list of edges into different subgraphs, if any
    // Unfortunately this requires us to traverse through the list of edges repeatedly

    int sgraphs =0;
    for( int i=0; i < g.m; i++ ) {
        node_t* first_node =&g.nodePtr[i];
        if( first_node->first != i ) continue; // Not a first node 

//        printf( "MSP for subgraph #%d:\n", sgraphs );

        double totalweight =0.f;
        int totaledges =0;

        for( int j=0; j < numedges; j++ ) {
            node_t* n =&g.nodePtr[edgelist[j].from];
            if( n->first != i ) continue; // This edge is not part of the subgraph
    //        printf( "(%d)->(%d) ", edgelist[j].from, edgelist[j].to );
            totalweight += edgelist[j].weight;
            totaledges++;
        }
        if( totaledges != 0 )
            printf( "[%d] %d edges, total weight: %f\n", ++sgraphs, totaledges, totalweight );
    }

    graph_free( &g );
    free( edgelist );

    return 0;
}
