#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include "htslib/sam.h"
#include "htslib/thread_pool.h"


int main(int argc, char *argv[])
{
    int c, ret = 0;
    samFile *in = 0, *out = 0;
    sam_hdr_t *header = NULL;
    char *fn_out = 0;
    int threads = 1;
    htsThreadPool p = {NULL, 0};

    while ((c = getopt(argc, argv, "o:@:")) >= 1){
        if (c=='o') fn_out = strdup(optarg);
        if (c=='@') threads =  atoi(optarg);
    }

    if (!fn_out){
        fprintf(stderr, "No output file, please specify output path\n");
        return 1;
    }

    if ((in = hts_open_format("-", "r", NULL)) == 0) {
        fprintf(stderr, "alignment failed, check options\n");
        ret = 1;
        goto view_end;
    }
    if ((out = hts_open_format(fn_out, "wb", NULL)) == 0) {
        fprintf(stderr, "failed to open \"%s\" for writing output path not valid\n", fn_out);
        hts_close(in);
        return 1;
    }
    if ((header = sam_hdr_read(in)) == 0) {
        fprintf(stderr, "alignment failed, check options\n");
        ret = 1;
        goto view_end;
    }
    if (sam_hdr_write(out, header) != 0) {
        fprintf(stderr, "failed to write the SAM header\n");
        ret = 1;
        goto view_end;
    }
    if (threads > 1) {
        if (!(p.pool = hts_tpool_init(threads))) {
            fprintf(stderr, "Error creating thread pool\n");
            ret = 1;
            goto view_end;
        }
        hts_set_opt(in,  HTS_OPT_THREAD_POOL, &p);
        if (out) hts_set_opt(out, HTS_OPT_THREAD_POOL, &p);
    } 
    bam1_t *b = bam_init1();
    int r;
    while ((r = sam_read1(in, header, b)) >= 0) { // read one alignment from `in'
        if (( r = sam_write1(out, header, b)) <0) break;}
    if (r < -1) {
        fprintf(stderr, "truncated file.\n");
        ret = 1;
    }
    bam_destroy1(b);
    view_end:

    // close files, free and return
    if (in){
        if (hts_close(in)){
            fprintf(stderr, "failed to close stdin\n");
            }
        }
    if (out){
        if (hts_close(out)){
            fprintf(stderr, "failed to close alignment output\n");
        }
    }
    free(fn_out);
    sam_hdr_destroy(header);
    if (p.pool) hts_tpool_destroy(p.pool);
    return ret;
}