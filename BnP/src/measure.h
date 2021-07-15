#ifndef _MEASURE_H
#define _MEASURE_H

#include <cstddef>
#include <cstdlib>
#include <cstring>

typedef struct {
    double total;
    int index;
} MEASURE_BLOCK;

MEASURE_BLOCK ** block_arr = NULL;
int blocks_size = 0;

void block_store(MEASURE_BLOCK *);
void measure_block_init(MEASURE_BLOCK *);
void measure_before(MEASURE_BLOCK block);
void measure_after(MEASURE_BLOCK block);
double measure_total(MEASURE_BLOCK block);
void measure_end();


void block_store(MEASURE_BLOCK * block){
    MEASURE_BLOCK ** cpy = (MEASURE_BLOCK **) calloc(blocks_size, sizeof(MEASURE_BLOCK *));
    if(block_arr != NULL){
        memcpy(cpy, block_arr, sizeof(MEASURE_BLOCK *)*blocks_size);
        free(block_arr);
    }

    block_arr = (MEASURE_BLOCK **) calloc(++blocks_size, sizeof(MEASURE_BLOCK *));
    memcpy(block_arr, cpy, sizeof(MEASURE_BLOCK *)*(blocks_size-1));
    block_arr[blocks_size - 1] = block;

    free(cpy);
}

void measure_block_init(MEASURE_BLOCK * block){
    block->total = 0;

    block_store(block);
    
    block->index = blocks_size - 1;
}

double now()
{
    struct timespec now;
    clock_gettime(CLOCK_REALTIME, &now);
    return now.tv_sec + now.tv_nsec*1e-9;
}

void measure_before(MEASURE_BLOCK block){
    int i = block.index;

    block_arr[i]->total -= now();
}

void measure_after(MEASURE_BLOCK block){
    int i = block.index;

    block_arr[i]->total += now();
}

double measure_total(MEASURE_BLOCK block){
    return block.total;
}

void measure_end(){
    free(block_arr);
}

#endif
