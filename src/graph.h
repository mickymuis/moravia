#ifndef __MATRIX_H__
#define __MATRIX_H__
#include <stdlib.h>
#include <stddef.h>
#include <inttypes.h>
#include <stdbool.h>

#if defined(_MSC_VER)
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#else
#include <sys/types.h>
#endif

typedef int32_t idx_t;  // Type large enough to hold an arbitrary index 

typedef struct {
    idx_t    index;     // Index of the node
    uint32_t degree;    // Number of edges
    uint32_t size;      // Size of the array
    idx_t*   edges;     // Array of node indices denoting an edge destination
    float*   weights;   // Matching array of weights
    int      level;     // Used internally by the Boruvka algorithm
    //uint32_t group;     // Group ID if the node belongs to a hypergraph, 0 if it doesn't
    idx_t    first;
    idx_t    next;
} node_t;

typedef struct {
    node_t** nodePtr;
    uint32_t m;         // Number of nodes
    uint32_t size;      // Total number of stored weights/indices
} graph_t;

void graph_dump( const graph_t* mat );

bool graph_loadMM( const char *filename,  graph_t *mat );

void graph_free( graph_t* mat );

#endif /* __MATRIX_H__ */
