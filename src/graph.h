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
    int      level;     // Used internally by the Boruvka algorithm
    idx_t    first;
    idx_t    next;
} node_t;

typedef struct {
    uint32_t size;      // Size of the array
    uint32_t degree;
    idx_t*   edges;     // Array of node indices denoting an edge destination
    double*   weights;   // Matching array of weights

} nodedata_t;

typedef struct {
    node_t* nodePtr;
    nodedata_t* dataPtr;
    uint32_t m;         // Number of nodes
    uint32_t size;      // Total number of stored weights/indices
} graph_t;

void graph_dump( const graph_t* mat );

bool graph_loadMM( const char *filename,  graph_t *mat );

void graph_free( graph_t* mat );

#endif /* __MATRIX_H__ */
