#ifndef __NODESET_H__
#define __NODESET_H__

#include "graph.h"

typedef struct {
    idx_t* idxPtr;
    unsigned int capacity;
    unsigned int count;
} nodeset_t;

void nodeset_init( nodeset_t* ns, unsigned int capacity );

void nodeset_free( nodeset_t* ns );

void nodeset_append( nodeset_t* ns, idx_t i );

void nodeset_print( nodeset_t* ns );

#endif
