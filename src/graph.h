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
    idx_t from, to;
} edge_t;

typedef struct {
    idx_t from, to, hindex;
    float weight;
} cedge_t;

typedef struct {
    idx_t    index;       // Index of the hypergraph
    edge_t*  edges;       // Pointer to an array of edges
    double   totalweight; // Total weight of the hypergraph
    uint32_t n;           // Number of edges
    uint32_t capacity;    // Capacity of the array
} stree_t;

typedef struct {
    uint32_t size;        // Size of the array
    uint32_t degree;      // Number of outgoing edges
    idx_t*   edges;       // Array of node indices denoting an edge destination
    float*  weights;     // Matching array of weights
    idx_t    index;       // Index in the array
} node_t;

typedef struct {
    node_t* nodePtr;      // Pointer to an array of nodes
    uint32_t m;           // Number of nodes
    uint32_t size;        // Total number of stored weights/indices
} graph_t;

typedef struct {
    stree_t* treePtr;     // Pointer to an array of mst nodes
    uint32_t n;           // Total number of subtrees
    uint32_t m;           // Number of nodes in the graph
    idx_t* hindex;        // Hypernode index for each graph node
    edge_t* edgePtr;      // Buffer that stores all currently processed edges
    uint32_t nedges;      // Number of edges so far
    uint32_t capacity;    // Capacity of the arrays
} mst_t;

void graph_dump( const graph_t* mat );

bool graph_loadMM( const char *filename,  graph_t *mat );

void graph_free( graph_t* mat );

#endif /* __MATRIX_H__ */
