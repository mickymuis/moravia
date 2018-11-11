#include "matrix.h"
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <assert.h>

/** Compute the contracted node from @n1 and @n2 by 
 *  computing the union of their edges.
 *  For edges with identical end nodes, the minimum edge is taken.
 */
row_t*
contractNodes( matrix_t* g, row_t* n1, row_t* n2 ) {
    row_t* contracted =malloc( sizeof(row_t) );
    // Allocate the maximum amount of memory we could need, we will free it partially later
    contracted->values =malloc( (n1->nz+n2->nz) * sizeof(float) );
    contracted->cols =malloc( (n1->nz+n2->nz) * sizeof(idx_t) );

    contracted->index =n1->index;

    // Compute the union of edges
    size_t i1 =0, i2 =0, off =0;
    while( i1 < n1->nz || i2 < n2->nz ) {
        idx_t idx =0;
        float weight1 =INFINITY, weight2 =INFINITY;

        if( i1 < n1->nz && ( i2 == n2->nz || n1->cols[i1] <= n2->cols[i2] )) {
            idx =n1->cols[i1];
            weight1 = n1->values[i1];
        }
        if( i2 < n2->nz && ( i1 == n1->nz || n2->cols[i2] <= n1->cols[i1] )) {
            idx =n2->cols[i2];
            weight2 = n2->values[i2];
        }

        if( !isinf( weight1 ) ) i1++;
        if( !isinf( weight2 ) ) i2++;
        
        // We don't want to introduce self-edges
        if( idx == n1->index || idx == n2->index || idx < 0 ) continue;

        contracted->cols[off] =idx;
        contracted->values[off] =fminf( weight1, weight2 );
        // Now, unfortunately, we have to fix all edges from adjacent nodes to point to this supernode
        row_t* adjnode = g->rowPtr[idx];
        idx_t n1_offs =-1, n2_offs =-1;
        assert( adjnode != NULL ); // Assert graph is consistent
        for( size_t i =0; i < adjnode->nz; i++ ) {
            if( adjnode->cols[i] == n1->index )
                n1_offs =i;
            else if( adjnode->cols[i] == n2->index )
                n2_offs =i;
        }
        // Fix edges pointing to n2 by pointing them to n1 instead
        if( n2_offs != -1 ) {
            // Take the min weight
            if( n1_offs != -1 ) {
                adjnode->values[n1_offs] = fminf( adjnode->values[n1_offs], adjnode->values[n2_offs] );
                adjnode->cols[n2_offs] = -1;
                        
            } else
                adjnode->cols[n2_offs] = n1->index;
        }

        off++;
    }
    // Resize the arrays to free unused space
    contracted->values =realloc( contracted->values, off * sizeof(float) );
    contracted->cols =realloc( contracted->cols, off * sizeof(idx_t) );
    contracted->nz =off;


    return contracted;
}


bool
computeMinEdge( matrix_t* mat, row_t* node, idx_t* to, float* weight ) {
    
    double minWeight = INFINITY;
    idx_t minTo =0;
    for( idx_t i =0; i < node->nz; i++ ) {
        if( node->values[i] < minWeight && node->cols[i] > 0 ) {
            minWeight = node->values[i];
            minTo = node->cols[i];
        }
    }
    if( isinf( minWeight ) ) return false;

    *weight = minWeight;
    *to = minTo;

    return true;
}

bool
minDegree( row_t* nodes[], size_t m, idx_t* idx ) {
    idx_t minDeg = INT_MAX;
    idx_t minRow = 0;
    for( size_t i =0; i < m; i++ ) {
        // This node has been removed
        if( nodes[i] == NULL ) continue;
        // This node has been contracted with another node 
        if( nodes[i]->index != i ) continue;
        if( nodes[i]->nz < minDeg && nodes[i]->nz > 0 ) {
            minDeg =nodes[i]->nz;
            minRow =i;
        }
    }
    if( minDeg == INT_MAX ) return false;

    *idx =minRow;
    return true;
}

int
boruvka( matrix_t* mat ) {

    const int np =4;
    const size_t m =mat->m;
    row_t* candidate[m];
    row_t* selected[np];
    float totalWeight =0.f; // Total weight of the MSP


    while(1) {

        memcpy( candidate, mat->rowPtr, m * sizeof(row_t*) );
        size_t ns =0;
        while( ns < np ) {
            idx_t sidx;
            // Take the candidate node with the lowest degree, if any
            if( !minDegree( candidate, m, &sidx ) ) break;

            row_t* selectedRow =candidate[sidx];
            candidate[sidx] = NULL;

            printf( "debug: selected row %d=%d degree=%d\n", sidx, selectedRow->index, selectedRow->nz );

            // Traverse all the nodes connected to the selected node
            for( int i =0; i < selectedRow->nz; i++ ) {
                // Obtain the ith node that is connected to the selected node
                if( selectedRow->cols[i] < 0 ) continue; // Removed

                idx_t idx1 =selectedRow->cols[i];
                row_t* row1 =candidate[idx1];
                if( row1 == NULL ) continue;

                //idx1 =row1->index; // Take the actual index
                // Remove it from the candidate list
                candidate[idx1] = NULL;

                // Traverse the second order connected nodes to the selected node
                for(int j =0; j < row1->nz; j++ ) {
                    if( row1->cols[j] < 0 ) continue; // Removed edge

                    idx_t idx2 =row1->cols[j];
                    row_t* row2 =candidate[idx2];
                    if( row2 == NULL ) continue;

                    //idx2 =row2->index; // Take the actual index
                    // Remove it from the candidate list
                    candidate[idx2] = NULL;
                }

            }
            // Add the selected node to the list
            selected[ns++] = selectedRow;

        }

        if( !ns ) break;

        printf( "Merging %ld selected nodes:\n", ns );

        for( idx_t i =0; i < ns; i++ ) {
            idx_t j;
            float weight;
            row_t* snode =selected[i];

            if( !computeMinEdge( mat, snode, &j, &weight ) ) continue;

            // This edge is part of the minimum spanning tree
            totalWeight += weight;
            printf( "(%d)-(%d) with weight %f\n", snode->index,j,weight );

            row_t* dnode =mat->rowPtr[j];

            // Contract the source node with the destination node of the edge
            row_t* cnode =contractNodes( mat, snode, dnode );
            // Replace the source node with the contracted node
            mat->rowPtr[cnode->index] =cnode;
            candidate[cnode->index] =cnode;
            mat->rowPtr[dnode->index] =NULL;
            candidate[dnode->index] =NULL;
            // Destroy the constituents
            free( snode ); free( dnode );
        }

        printf ("DEBUG\n-----\n");
        matrix_dump( mat );

    }

    printf( "Total weight: %f\n", totalWeight );
    return 0;
}

int
main( int argc, const char** argv ) {

    if( argc != 2 ) {
        fprintf( stderr, "(i) Usage: %s [matrix market file]\n", argv[0] );
        return -1;
    }
    matrix_t mat;

    if( !matrix_loadMM( argv[1], &mat ) ) {
        fprintf( stderr, "(e) Could not load matrix market file `%s'\n", argv[1] );
        return -1;
    }
    
    matrix_dump( &mat );

    boruvka( &mat );
    matrix_dump( &mat );
    

    matrix_free( &mat );

    return 0;
}
