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
    float weight;
} edge_t;

void
printSuperNodes( graph_t* g ) {

    for( idx_t i=0; i < g->m; i++ ) {
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

/** Compute the contracted node from @n1 and @n2 by 
 *  merging their supernode pointer chains together.
 *  The level of all nodes in the chain is set to the highest level + 1
 */
int
contractNodes( graph_t* g, node_t* n1, node_t* n2 ) {

    if( n1->first == n2->first ) return -1;

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

    return level;
}

bool
computeMinEdge( graph_t* g, node_t* node, idx_t* from, idx_t* to, float* weight ) {
    
    double minWeight = INFINITY;
    idx_t minTo =0, minFrom =0;
/*    float minWeights[m];
    memset( minWeights, INFINITY, g->m * sizeof(float) );*/

    node_t *ptr =g->nodePtr[node->first];

    while( 1 ) {
        for( idx_t i =0; i < ptr->size; i++ ) {
            if( ptr->weights[i] < minWeight ) {
                // TODO OPTIMIZE 
                
    /*            node_t *ptr2 =g->nodePtr[node->first];
                while( 1 ) {
                    if( ptr2->index == ptr->edges[i] ) goto SELFEDGE;
                    if( ptr2->next == -1 ) break;
                    ptr2 = g->nodePtr[ptr2->next];
                }*/
                if( g->nodePtr[ptr->edges[i]]->first == node->first ) continue; // Selfedge

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
minDegree( const node_t* nodes[], size_t m, int level, idx_t* idx ) {
    idx_t minDeg = INT_MAX;
    idx_t minRow = 0;


    for( size_t i =0; i < m; i++ ) {
        // This node has been removed
        if( nodes[i] == NULL ) continue;
        // The supernode to which this node belongs has been removed
        if( nodes[nodes[i]->first] == NULL ) continue; 
        if( /*level != 0 && */nodes[i]->level == level ) {
            minDeg =nodes[i]->degree;
            minRow =i;
            break;
        }
       /* if( level == 0 && nodes[i]->level == level && nodes[i]->degree < minDeg && nodes[i]->size > 0 ) {
            minDeg =nodes[i]->degree;
            minRow =i;
        }*/
    }
    if( minDeg == INT_MAX ) return false;

    *idx =minRow;
    return true;
}

void
removeSuperNode( graph_t* g, node_t* nodes[], node_t* n ) {
    if( !n ) return;

    node_t* ptr =g->nodePtr[n->first];
    nodes[ptr->index] = NULL;

    // Remove one or more nodes depending on the size of the super node n belongs to
   /* while( 1 ){
        nodes[ptr->index] = NULL;
        if( ptr->next == -1 ) break;
        ptr = g->nodePtr[ptr->next];
    }*/
}

void
removeSubGraph( graph_t* g, node_t* nodes[], node_t* n ) {
    if( !n ) return;

    removeSuperNode( g, nodes, n );

    /*for( size_t i=0; i < n->size; i++ ) {
        node_t *ptr =g->nodePtr[n->edges[i]];
        removeSuperNode( g, nodes, ptr );
        for( size_t j=0; j < ptr->size; j++ ) {
            node_t *ptr2 =nodes[ptr->edges[j]];
            removeSuperNode( g, nodes, ptr2 );
        }
    }*/
}

int
boruvka( graph_t* g, edge_t *edgelist ) {

    int numedges =0;
    const int np =8;
    const size_t m =g->m;
    node_t* candidate[m];
    idx_t selected[np];
    edge_t bestedges[np];
    int level =0, maxlevel =0;
    memcpy( candidate, g->nodePtr, m * sizeof(node_t*) );

    while(1) {

    //    memcpy( candidate, g->nodePtr, m * sizeof(node_t*) );
        size_t ns =0;
        while( ns < np ) {
            idx_t sidx;
            // Take the candidate node with the lowest degree at the current level, if any
            if( !minDegree( candidate, m, level, &sidx ) ) {
                if( ns > 0 ) break;
                if( ++level > maxlevel ) break;
                continue;
            }

            node_t* selectedNode =candidate[sidx];
            assert( selectedNode != NULL );

            //printf( "debug: selected node %d=%d degree=%d\n", sidx, selectedNode->index, selectedNode->degree );

            removeSubGraph( g, candidate, selectedNode );

            // Add the selected node to the list
            selected[ns++] = sidx;

        }

        if( !ns ) break;

        //printf( "Merging %ld selected nodes\n", ns );
        
//#pragma omp parallel for
        for( idx_t i =0; i < ns; i++ ) {
            idx_t from, to;
            float weight;
            node_t* snode =g->nodePtr[selected[i]];

            bestedges[i].from =-1;

            //printf( "from %d ", selected[i] );

            if( !computeMinEdge( g, snode, &from, &to, &weight ) ) continue;
            
            //printf( "to %d\n", to );

            edge_t e = { .from = from, .to = to, .weight =weight };
            bestedges[i] =e;

        }

        int mops =0; // Number of successfull merges operations

        for( idx_t i =0; i < ns; i++ ) {

            if( bestedges[i].from == -1 ) continue;
            // This edge is part of the minimum spanning tree
            //totalWeight += weight;
            //printf( "(%d)-(%d) with weight %f\n", from, to, weight );
            edge_t e =bestedges[i];

            node_t* dnode =g->nodePtr[e.to];
            assert( dnode->index == e.to );

            node_t* snode =g->nodePtr[e.from];
            assert( snode->index == e.from );

            // Contract the source node with the destination node of the edge
            int l = contractNodes( g, snode, dnode );
            if( l != -1 ) {
                edgelist[numedges++] =e;
                mops++;
                maxlevel =fmaxf( maxlevel, l );
                node_t* first_node =g->nodePtr[snode->first];
                candidate[first_node->index] =first_node;
            }
            //printf( "." );
        }

        if( !mops ) break; // There is nothing left to do

//        printSuperNodes( g );
        /*char c;
        scanf ( "%c", &c );*/

    }

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
    
    //graph_dump( &g );

    int numedges = boruvka( &g, (idx_t *)edgelist );
   // graph_dump( &g );
  
    // We now need to spit the raw list of edges into different subgraphs, if any
    // Unfortunately this requires us to traverse through the list of edges repeatedly

    int sgraphs =0;
    for( int i=0; i < g.m; i++ ) {
        node_t* first_node =g.nodePtr[i];
        if( first_node->first != i ) continue; // Not a first node 

//        printf( "MSP for subgraph #%d:\n", sgraphs );

        float totalweight =0.f;

        for( int j=0; j < numedges; j++ ) {
            node_t* n =g.nodePtr[edgelist[j].from];
            if( n->first != i ) continue; // This edge is not part of the subgraph
    //        printf( "(%d)->(%d) ", edgelist[j].from, edgelist[j].to );
            totalweight += edgelist[j].weight;
        }
        printf( "[%d] Total weight: %f\n", ++sgraphs, totalweight );
    }

    graph_free( &g );
    free( edgelist );

    return 0;
}
