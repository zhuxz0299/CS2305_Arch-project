// Xiaozhi Zhu 521021910299
#define _GNU_SOURCE
#include "cachelab.h"
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#define HELP_STRING "Usage: ./csim [-hv] -s <num> -E <num> -b <num> -t <file>\
\n Options:\
\n\t  -h         Print this help message.\
\n\t  -v         Optional verbose flag.\
\n\t  -s <num>   Number of set index bits.\
\n\t  -E <num>   Number of lines per set.\
\n\t  -b <num>   Number of block offset bits.\
\n\t  -t <file>  Trace file.\
\n\
\nExamples:\
\n\t  linux>  ./csim-ref -s 4 -E 1 -b 4 -t traces/yi.trace\
\n\t  linux>  ./csim-ref -v -s 8 -E 2 -b 4 -t traces/yi.trace"

typedef struct cmd_opts
{
    bool verbose;
    unsigned int s;
    unsigned int E;
    unsigned int b;
    char *tracefile;
} cmd_opts;

typedef struct cache_line
{
    bool valid_bit;
    unsigned int tag;
    unsigned int lru_counter;
    char *block;
} cache_line;

typedef struct cache_set
{
    unsigned int E;
    cache_line **lines;
} cache_set;

typedef struct cache
{
    unsigned int hit_count;
    unsigned int miss_count;
    unsigned int eviction_count;

    unsigned int s;
    unsigned int S;
    unsigned int b;
    unsigned int B;
    unsigned int E;

    cache_set **sets;
} cache;

typedef struct addr_seg
{
    unsigned int tag;
    unsigned int index;
    unsigned int block_offset;
} addr_seg;

cache_line *create_cache_line(unsigned int B);
cache_set *create_cache_set(unsigned int E, unsigned int B);
cache *create_cache(cmd_opts *opts);
cmd_opts *parse_cmd(int argc, char *argv[]);
addr_seg *parse_addr(unsigned long long addr, cmd_opts *opts);
void run_sim(cache *cache, cmd_opts *opts);
void load(unsigned long long addr, cache *cache, cmd_opts *opts);
void store(unsigned long long addr, cache *cache, cmd_opts *opts);
void free_cache(cache *cache);

int main(int argc, char *argv[])
{
    cmd_opts *cmd_opts = parse_cmd(argc, argv);
    cache *cache = create_cache(cmd_opts);
    run_sim(cache, cmd_opts);
    printSummary(cache->hit_count, cache->miss_count, cache->eviction_count);
    free_cache(cache);
    return 0;
}

cache_line *create_cache_line(unsigned int B)
{
    cache_line *line = (cache_line *)malloc(sizeof(cache_line));
    line->valid_bit = false;
    line->tag = 0;
    line->lru_counter = 0;
    line->block = (char *)malloc(sizeof(char) * B);
    return line;
}

cache_set *create_cache_set(unsigned int E, unsigned int B)
{
    cache_set *set = (cache_set *)malloc(sizeof(cache_set));
    set->E = E;
    set->lines = (cache_line **)malloc(sizeof(cache_line *) * E);
    for (int i = 0; i < E; i++)
        set->lines[i] = create_cache_line(B);
    return set;
}

cache *create_cache(cmd_opts *opts)
{
    cache *c = (cache *)malloc(sizeof(cache));
    c->hit_count = 0;
    c->miss_count = 0;
    c->eviction_count = 0;

    c->s = opts->s;
    c->S = 1 << opts->s;
    c->b = opts->b;
    c->B = 1 << opts->b;
    c->E = opts->E;

    c->sets = (cache_set **)malloc(sizeof(cache_set *) * c->S);
    for (int i = 0; i < c->S; i++)
        c->sets[i] = create_cache_set(c->E, c->B);
    return c;
}

cmd_opts *parse_cmd(int argc, char *argv[])
{
    cmd_opts *opts = (cmd_opts *)malloc(sizeof(cmd_opts));
    opts->verbose = false;
    opts->s = 0;
    opts->E = 0;
    opts->b = 0;
    opts->tracefile = NULL;

    extern char *optarg;
    char opt;
    while ((opt = getopt(argc, argv, "hvs:E:b:t:")) != -1)
    {
        switch (opt)
        {
        case 'h':
            printf(HELP_STRING);
            exit(0);
        case 'v':
            opts->verbose = true;
            break;
        case 's':
            opts->s = atoi(optarg);
            break;
        case 'E':
            opts->E = atoi(optarg);
            break;
        case 'b':
            opts->b = atoi(optarg);
            break;
        case 't':
            opts->tracefile = optarg;
            break;
        default:
            printf(HELP_STRING);
            exit(1);
        }
    }
    return opts;
}

addr_seg *parse_addr(unsigned long long addr, cmd_opts *opts)
{
    addr_seg *seg = (addr_seg *)malloc(sizeof(addr_seg));
    seg->tag = addr >> (opts->s + opts->b);
    seg->index = (addr >> opts->b) & ((1 << opts->s) - 1);
    seg->block_offset = addr & ((1 << opts->b) - 1);
    return seg;
}

void run_sim(cache *cache, cmd_opts *opts)
{
    FILE *fp;
    int read_len = 0;
    size_t buffer_size = 0;
    char *line = 0;
    fp = fopen(opts->tracefile, "r");
    while ((read_len = getline(&line, &buffer_size, fp)) != -1)
    {
        if (read_len > 0 && line[read_len - 1] == '\n')
            line[read_len - 1] = '\0'; // remove trailing newline
        if (line[0] == 'I')
            continue;
        if (opts->verbose)
            printf("%s", line);

        char op = line[1];
        char *addr_str = strtok(line + 3, ",");
        unsigned long long addr = strtoull(addr_str, NULL, 16);

        // printf("\nDEBUG:addr=%llx\n", addr);

        switch (op)
        {
        case 'L':
            load(addr, cache, opts);
            break;
        case 'S':
            store(addr, cache, opts);
            break;
        case 'M':
            load(addr, cache, opts);
            store(addr, cache, opts);
            break;
        }

        if (opts->verbose)
            printf("\n");
    }
    fclose(fp);
}

void load(unsigned long long addr, cache *cache, cmd_opts *opts)
{
    addr_seg *addr_seg = parse_addr(addr, opts);
    cache_set *current_set = cache->sets[addr_seg->index];

    // for (int i = 0; i < cache->E; i++)
    //     printf("\nDEBUG:current_set->lines[%d]->valid_bit=%d\n", i, current_set->lines[i]->valid_bit);

    // check if hit, and choose the line to evict
    int hit_line = -1, eviction_line = 0;
    for (int i = 0; i < cache->E; ++i)
    {
        if (current_set->lines[i]->valid_bit == false)
        {
            eviction_line = i;
            break;
        }
        if (current_set->lines[i]->tag == addr_seg->tag)
        {
            hit_line = i;
            break;
        }
        if (current_set->lines[i]->lru_counter > current_set->lines[eviction_line]->lru_counter)
            eviction_line = i;
    }

    // printf("\nDEBUG:hit_line=%d,eviction_line=%d\n", hit_line, eviction_line);
    // for (int i = 0; i < cache->E; i++)
    //     printf("\nDEBUG:current_set->lines[%d]->lru_counter=%d\n", i, current_set->lines[i]->lru_counter);

    if (hit_line != -1) // hit
    {
        cache->hit_count++;
        // update lru counter
        current_set->lines[hit_line]->lru_counter = 0;
        for (int i = 0; i < cache->E; ++i)
            if (i != hit_line)
                current_set->lines[i]->lru_counter++;

        if (opts->verbose)
            printf(" hit");
    }
    else
    {
        cache->miss_count++;
        if (opts->verbose)
            printf(" miss");
        if (current_set->lines[eviction_line]->valid_bit)
        {
            cache->eviction_count++;
            if (opts->verbose)
                printf(" eviction");
        }
        current_set->lines[eviction_line]->valid_bit = true;
        current_set->lines[eviction_line]->tag = addr_seg->tag;
        current_set->lines[eviction_line]->lru_counter = 0;
        for (int i = 0; i < cache->E; ++i)
            if (i != eviction_line)
                current_set->lines[i]->lru_counter++;
    }

    free(addr_seg);
}

void store(unsigned long long addr, cache *cache, cmd_opts *opts)
{
    load(addr, cache, opts);
}

void free_cache(cache *cache)
{
    for (int i = 0; i < cache->S; ++i)
    {
        for (int j = 0; j < cache->E; ++j)
        {
            free(cache->sets[i]->lines[j]->block);
            free(cache->sets[i]->lines[j]);
        }
        free(cache->sets[i]->lines);
        free(cache->sets[i]);
    }
    free(cache->sets);
    free(cache);
}