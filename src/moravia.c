#include "graph.h"
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <assert.h>

void
printSuperNodes( graph_t* g ) {

    for( size_t i=0; i < g->m; i++ ) {
        printf( "%d: ", i );
        node_t* n =g->nodePtr[g->nodePtr[i]->first];
        while( 1 ) {
            printf( "%d -> ", n->index );  
            if( n->next == -1 ) break;
            n = g->nodePtr[n->next];
        }

        printf( "X\n" );
    }
}

void
contractNodes( graph_t* g, node_t* n1, node_t* n2 ) {
    node_t* n1_last = g->nodePtr[n1->first];
    int level = (int)fmaxf( n1->level, n2->level) + 1;
    
    // Find the last node in n1's chain and update all nodes' level
    while( 1 ) {
        n1_last->level = level;
        if( n1_last->next == -1 ) break;
        n1_last =g->nodePtr[n1_last->next];
    }
     
    n1_last->next = n2->first;
    node_t* n =g->nodePtr[n2->first];

    // Also update the levels and first pointers in n2's chain
    while( 1 ) {
        n->level =level;
        n->first =n1->first;
        if( n->next == -1 ) break;
        n = g->nodePtr[n->next];
    }
}

/** Compute the contracted node from @n1 and @n2 by 
 *  computing the union of their edges.
 *  For edges with identical end nodes, the minimum edge is taken.
 */
node_t*
contractNodes2( graph_t* g, node_t* n1, node_t* n2 ) {
    node_t* contracted =malloc( sizeof(node_t) );
    // Allocate the maximum amount of memory we could need, we will free it partially later
    contracted->weights =malloc( (n1->size+n2->size) * sizeof(float) );
    contracted->edges =malloc( (n1->size+n2->size) * sizeof(idx_t) );

    contracted->index =n1->index;

    // Compute the union of edges
    size_t i1 =0, i2 =0, off =0;
    while( i1 < n1->size || i2 < n2->size ) {
        idx_t idx =-1;
        float weight1 =INFINITY, weight2 =INFINITY;

        if( i1 < n1->size && ( i2 == n2->size || n1->edges[i1] <= n2->edges[i2] )) {
            idx =n1->edges[i1];
            weight1 = n1->weights[i1];
        }
        if( i2 < n2->size && ( i1 == n1->size || n2->edges[i2] <= n1->edges[i1] )) {
            idx =n2->edges[i2];
            weight2 = n2->weights[i2];
        }

        if( !isinf( weight1 ) ) i1++;
        if( !isinf( weight2 ) ) i2++;
        
        // We don't want to introduce self-edges
        if( idx == n1->index || idx == n2->index || idx < 0 ) continue;


        contracted->edges[off] =idx;
        contracted->weights[off] =fminf( weight1, weight2 );
        
        //printf( "-%ld: %d\n", off, contracted->edges[off] );
        if( off ) assert( contracted->edges[off-1] < contracted->edges[off] );

        // Both n1 and n2 have an edge to the same node
        if( !isinf( weight1 ) && !isinf( weight2 ) ) {
            node_t* adjnode = g->nodePtr[idx];
            // Reduce the degree of this node by one. 
            // This is only important to optimize the search process
            adjnode->degree--;
        }
        
        // Now, unfortunately, we have to fix all edges from adjacent nodes to point to this supernode
#if 0
        node_t* adjnode = g->nodePtr[idx];
        idx_t n1_off =-1, n2_off =-1, lt_n1_off =-1;
        assert( adjnode != NULL ); // Assert graph is consistent
        for( size_t i =0; i < adjnode->size; i++ ) {
            if( adjnode->edges[i] == n1->index ) {
                assert( n1_off == -1 );
                n1_off =i;
            }
            else if( adjnode->edges[i] == n2->index ) {
                assert( n2_off == -1 );
                n2_off =i;
            }
            // For reordering the array later, lest it be neccessary
            if( adjnode->edges[i] != -1 && adjnode->edges[i] < n1->index ) lt_n1_off =i;
        }
        // Fix edges pointing to n2 by pointing them to n1 instead
        if( n2_off != -1 ) {
            // Take the min weight
            //printf( "debug: removing (%d)->(%d)\n", adjnode->index, n2->index );
            if( n1_off != -1 ) {
                adjnode->weights[n1_off] = fminf( adjnode->weights[n1_off], adjnode->weights[n2_off] );
                adjnode->edges[n2_off] = -1;
                        
            } else {
                /*printf( "Swapping n1=%d for n2=%d in (%d), before: ", n1->index, n2->index, adjnode->index );
                for( int i =0; i < adjnode->size; i++ ) printf( "%d ", adjnode->edges[i] );
                printf( "\n" );*/
                // We cannot just swap the index of n2 with the index of n1 as the array becomes unsorted
                // adjnode->edges[n2_off] = n1->index;
                float weight =adjnode->weights[n2_off];
                // Handle the edge cases first
                if( lt_n1_off == n2_off 
                 || lt_n1_off == n2_off-1 
                 || (n2_off == adjnode->size-1 && n1->index > n2->index )
                 || (n2_off == 0 && n1->index < n2->index ) ) {
                    lt_n1_off = n2_off;
                }
                else if( n1->index < n2->index ) {
                    lt_n1_off++;
                    memmove( &adjnode->edges[lt_n1_off+1], &adjnode->edges[lt_n1_off], (n2_off-lt_n1_off) * sizeof(idx_t) );
                    memmove( &adjnode->weights[lt_n1_off+1], &adjnode->weights[lt_n1_off], (n2_off-lt_n1_off) * sizeof(float) );

                } else { // n2 < n1
                    memmove( &adjnode->edges[n2_off], &adjnode->edges[n2_off+1], (lt_n1_off-n2_off) * sizeof(idx_t) );
                    memmove( &adjnode->weights[n2_off], &adjnode->weights[n2_off+1], (lt_n1_off-n2_off) * sizeof(float) );

                }
                adjnode->edges[lt_n1_off] =n1->index;
                adjnode->weights[lt_n1_off] =weight;

                adjnode->degree--; // We have actually removed an edge
                /*printf( "after: " );
                for( int i =0; i < adjnode->size; i++ ) printf( "%d ", adjnode->edges[i] );
                printf( "\n" );*/
            }
        }
#endif

        off++;
    }
    // Resize the arrays to free unused space
    contracted->weights =realloc( contracted->weights, off * sizeof(float) );
    contracted->edges =realloc( contracted->edges, off * sizeof(idx_t) );
    contracted->size =off;
    contracted->degree =off;


    return contracted;
}


bool
computeMinEdge( graph_t* g, node_t* node, idx_t* from, idx_t* to, float* weight ) {
    
    double minWeight = INFINITY;
    idx_t minTo =0, minFrom =0;

    node_t *ptr =g->nodePtr[node->first];

    while( 1 ) {
        for( idx_t i =0; i < ptr->size; i++ ) {
            if( ptr->weights[i] < minWeight ) {
                // TODO OPTIMIZE 
                
                node_t *ptr2 =g->nodePtr[node->first];
                while( 1 ) {
                    if( ptr2->index == ptr->edges[i] ) goto SELFEDGE;
                    if( ptr2->next == -1 ) break;
                    ptr2 = g->nodePtr[ptr2->next];
                }

                // END
                minWeight = ptr->weights[i];
                minFrom = ptr->index;
                minTo = ptr->edges[i];
SELFEDGE:;
            }
        }
        if( ptr->next == -1 ) break;
        ptr =g->nodePtr[ptr->next];
    }

    if( isinf( minWeight ) ) return false;

    *weight = minWeight;
    *from = minFrom;
    *to = minTo;

    return true;
}

bool
minDegree( node_t* nodes[], size_t m, int level, idx_t* idx ) {
    idx_t minDeg = INT_MAX;
    idx_t minRow = 0;
    for( size_t i =0; i < m; i++ ) {
        // This node has been removed
        if( nodes[i] == NULL ) continue;
        // This node has been contracted with another node 
        //if( nodes[i]->first != i ) continue;
        //assert( nodes[i]->index == i );
        if( level != 0 && nodes[i]->level == level ) {
            minDeg =nodes[i]->degree;
            minRow =i;
            break;
        }
        if( level == 0 && nodes[i]->level == level && nodes[i]->degree < minDeg && nodes[i]->size > 0 ) {
            minDeg =nodes[i]->degree;
            minRow =i;
        }
    }
    if( minDeg == INT_MAX ) return false;

    *idx =minRow;
    return true;
}

void
removeSuperNode( graph_t* g, node_t* nodes[], node_t* n ) {
    if( !n ) return;

    node_t* ptr =g->nodePtr[n->first];

    // Remove one or more nodes depending on the size of the super node n belongs to
    while( 1 ) {
        nodes[ptr->index] = NULL;
        if( ptr->next == -1 ) break;
        ptr = g->nodePtr[ptr->next];
    }
}

void
removeSubGraph( graph_t* g, node_t* nodes[], node_t* n ) {
    if( !n ) return;

    removeSuperNode( g, nodes, n );

    for( size_t i=0; i < n->size; i++ ) {
        node_t *ptr =g->nodePtr[n->edges[i]];
        removeSuperNode( g, nodes, ptr );
        for( size_t j=0; j < ptr->size; j++ ) {
            node_t *ptr2 =nodes[ptr->edges[j]];
            removeSuperNode( g, nodes, ptr2 );
        }
    }
}

int
boruvka( graph_t* g ) {

    const int np =16;
    const size_t m =g->m;
    node_t* candidate[m];
    idx_t selected[np];
    float totalWeight =0.f; // Total weight of the MSP
    int level =0;

    while(1) {

        memcpy( candidate, g->nodePtr, m * sizeof(node_t*) );
        size_t ns =0;
        while( ns < np ) {
            idx_t sidx;
            // Take the candidate node with the lowest degree at the current level, if any
            if( !minDegree( candidate, m, level, &sidx ) ) {
                if( ns > 0 ) break;
                //printf( "up!\n" );
                if( !minDegree( candidate, m, ++level, &sidx ) ) break;
            }

            node_t* selectedNode =candidate[sidx];
            assert( selectedNode != NULL );

            //printf( "debug: selected node %d=%d degree=%d\n", sidx, selectedNode->index, selectedNode->degree );

            removeSubGraph( g, candidate, selectedNode );

            // Add the selected node to the list
            selected[ns++] = sidx;

        }

        if( !ns ) break;

        printf( "Merging %ld selected nodes:\n", ns );

        for( idx_t i =0; i < ns; i++ ) {
            idx_t from, to;
            float weight;
            node_t* snode =g->nodePtr[selected[i]];

            if( !computeMinEdge( g, snode, &from, &to, &weight ) ) continue;

            // This edge is part of the minimum spanning tree
            totalWeight += weight;
            printf( "(%d)-(%d) with weight %f\n", from, to, weight );

            node_t* dnode =g->nodePtr[to];
            assert( dnode->index == to );

            // Contract the source node with the destination node of the edge
            contractNodes( g, snode, dnode );
            //cnode->level =level+1;
        }

        /*printf ("DEBUG\n-----\n");
        graph_dump( g );*/

        /*printSuperNodes( g );
        char c;
        scanf ( "%c", &c );*/

    }

    printf( "Total weight: %f\n", totalWeight );
    return 0;
}

int
main( int argc, const char** argv ) {

    if( argc != 2 ) {
        fprintf( stderr, "(i) Usage: %s [graph market file]\n", argv[0] );
        return -1;
    }
    graph_t g;

    if( !graph_loadMM( argv[1], &g ) ) {
        fprintf( stderr, "(e) Could not load graph market file `%s'\n", argv[1] );
        return -1;
    }
    
   // graph_dump( &g );

    boruvka( &g );
   // graph_dump( &g );
    

    graph_free( &g );

    return 0;
}
