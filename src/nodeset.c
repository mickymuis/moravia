#include "nodeset.h"
#include <stdio.h>

void
nodeset_init( nodeset_t* ns, unsigned int capacity ) {
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
