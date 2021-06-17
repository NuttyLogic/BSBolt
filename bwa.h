#ifndef BWA_H_
#define BWA_H_

#include <stdint.h>
#include "bntseq.h"
#include "bwt.h"

#define BWA_IDX_BWT 0x1
#define BWA_IDX_BNS 0x2
#define BWA_IDX_PAC 0x4
#define BWA_IDX_ALL 0x7

#define BWA_CTL_SIZE 0x10000

#define BWTALGO_AUTO  0
#define BWTALGO_RB2   1
#define BWTALGO_BWTSW 2
#define BWTALGO_IS    3

typedef struct {
	bwt_t    *bwt; // FM-index
	bntseq_t *bns; // information on the reference sequences
	uint8_t  *pac, *opac; // the actual 2-bit encoded reference sequences with 'N' converted to a random base

	int    is_shm;
	int64_t l_mem;
	uint8_t  *mem;
} bwaidx_t;

typedef struct {
	int l_seq, id, first, read_group, subsitution_pattern, alignment_score, mapped;
	int conversion, bs_conflict, crick, paired;
	char *name, *comment, *seq, *qual, *sam, *oseq;
} bseq1_t;

extern int bwa_verbose;
extern char bwa_rg_id[256];

#ifdef __cplusplus
extern "C" {
#endif

	bseq1_t *bseq_read(int64_t chunk_size, int *n_, void *ks1_, void *ks2_, int64_t *s, int conversion, 
                       int undirectional, float substitution_proportion);
	void bseq_classify(int n, bseq1_t *seqs, int m[2], bseq1_t *sep[2]);

	void bwa_fill_scmat(int a, int b, int8_t mat[25]);
	uint32_t *bwa_gen_cigar(const int8_t mat[25], int q, int r, int w_, int64_t l_pac, const uint8_t *pac, const uint8_t *opac, int l_query, uint8_t *query, int64_t rb, int64_t re, int *score, int *n_cigar, int *NM, int64_t crick_l, uint8_t *oquery, int *cg_meth, int *cg_unmeth, int *ch_meth, int *ch_unmeth);
	uint32_t *bwa_gen_cigar2(const int8_t mat[25], int o_del, int e_del, int o_ins, int e_ins, int w_, int64_t l_pac, const uint8_t *pac, const uint8_t *opac, int l_query, uint8_t *query, int64_t rb, int64_t re, int *score, int *n_cigar, int *NM, int64_t crick_l, uint8_t *oquery, int *cg_meth, int *cg_unmeth, int *ch_meth, int *ch_unmeth);

	int bwa_idx_build(const char *fa, const char *prefix, int algo_type, int block_size);

	char *bwa_idx_infer_prefix(const char *hint);
	bwt_t *bwa_idx_load_bwt(const char *hint);

	bwaidx_t *bwa_idx_load_from_shm(const char *hint);
	bwaidx_t *bwa_idx_load_from_disk(const char *hint, int which);
	bwaidx_t *bwa_idx_load(const char *hint, int which);
	void bwa_idx_destroy(bwaidx_t *idx);
	int bwa_idx2mem(bwaidx_t *idx);
	int bwa_mem2idx(int64_t l_mem, uint8_t *mem, bwaidx_t *idx);

	void bwa_print_sam_hdr(const bntseq_t *bns, const char *hdr_line);
	char *bwa_set_rg(const char *s);
	char *bwa_insert_header(const char *s, char *hdr);

#ifdef __cplusplus
}
#endif

#endif
