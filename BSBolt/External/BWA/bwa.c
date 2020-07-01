#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include <assert.h>
#include "bntseq.h"
#include "bs_helpers.h"
#include "bwa.h"
#include "ksw.h"
#include "utils.h"
#include "kstring.h"
#include "kvec.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

int bwa_verbose = 3;
char bwa_rg_id[256];
char *bwa_pg;

/************************
 * Batch FASTA/Q reader *
 ************************/

#include "kseq.h"
KSEQ_DECLARE(gzFile)

static inline void trim_readno(kstring_t *s)
{
	if (s->l > 2 && s->s[s->l-2] == '/' && isdigit(s->s[s->l-1]))
		s->l -= 2, s->s[s->l] = 0;
}

static inline char *dupkstring(const kstring_t *str, int dupempty)
{
	char *s = (str->l > 0 || dupempty)? malloc(str->l + 1) : NULL;
	if (!s) return NULL;

	memcpy(s, str->s, str->l);
	s[str->l] = '\0';
	return s;
}

static inline void kseq2bseq1(const kseq_t *ks, bseq1_t *s, int conversion, int first,
                              int read_group, int substitution_pattern)
{ // TODO: it would be better to allocate one chunk of memory, but probably it does not matter in practice
	s->name = strdup(ks->name.s);
	s->comment = ks->comment.l? strdup(ks->comment.s) : 0;
	s->seq = strdup(ks->seq.s);
	s->qual = ks->qual.l? strdup(ks->qual.s) : 0;
	s->l_seq = strlen(s->seq);
	s->first = first;
	s->read_group = read_group;
	s->subsitution_pattern = substitution_pattern;
	s->conversion = conversion;
	s->bs_conflict = 0;
	s->crick = 0;
	s->mapped = 0;
	s->paired = 0;
	s->alignment_score = 0;
	if (conversion){
		s->oseq = strdup(ks->seq.s);
		bsConversion(s->seq, s->subsitution_pattern, s->l_seq);
		for (int i = 0; i < s->l_seq; ++i){
			s->oseq[i] = s->oseq[i] < 4 ? s->oseq[i] : nst_nt4_table[(int)s->oseq[i]];
		}
	}
	else{
		s->oseq = NULL;
	}
}

bseq1_t *bseq_read(int64_t chunk_size, int *n_, void *ks1_, void *ks2_, int64_t *s, int conversion, 
                        int undirectional, float substitution_proportion)
{
	kseq_t *ks = (kseq_t*)ks1_, *ks2 = (kseq_t*)ks2_;
	int64_t size = 0, m, n;
	bseq1_t *seqs;
	m = n = 0; seqs = 0;
	while (kseq_read(ks) >= 0)
	{
		if (ks2 && kseq_read(ks2) < 0) { // the 2nd file has fewer reads
			fprintf(stderr, "[W::%s] the 2nd file has fewer sequences.\n", __func__);
			break;
		}
		// reallocate memory m -2 to account for undirectional library bseq read number uncertainty 
		if (n >= m - 2) {
			m = m? m<<1 : 256;
			seqs = (bseq1_t*) realloc(seqs, m * sizeof(bseq1_t));
		}
		trim_readno(&ks->name);
		int substitution_pattern = 0;
		int compare_reads = 0;
		// assess read conversion type for undirectional reads
		if (undirectional){
			int un_type = 0;
			if (ks2) un_type = assessConversion(ks->seq.s, ks2->seq.s, 1, substitution_proportion);
			else un_type = assessConversion(ks->seq.s, ks->seq.s, 0, substitution_proportion);
			if (un_type == 2) compare_reads = 1;
			else  substitution_pattern = un_type;
			if (substitution_pattern >1) fprintf(stderr, "conversion error\n");
		}
		kseq2bseq1(ks, &seqs[n], conversion, 0, 0, substitution_pattern);
		seqs[n].id = n;
		//{
		//	size += strlen(seqs[n].name);
		//	size += strlen(seqs[n].comment);
		//	size += strlen(seqs[n].qual);
		//	// fprintf(stderr, "qual len: %d %d\n", strlen(seqs[n].qual), seqs[n].l_seq);
		//	size += 7; // non accounted chars
		//}
		size += seqs[n++].l_seq;

		if (ks2) {
			trim_readno(&ks2->name);
			kseq2bseq1(ks2, &seqs[n], conversion, 1, 0, substitution_pattern ? 0: 1);
			seqs[n].id = n;
			size += seqs[n++].l_seq;
		}
		if (compare_reads){
			trim_readno(&ks->name);
			kseq2bseq1(ks, &seqs[n], conversion, 0, 1, 1);
			seqs[n].id = n;
			
			size += seqs[n++].l_seq;

			if (ks2) {
				trim_readno(&ks2->name);
				kseq2bseq1(ks2, &seqs[n], conversion, 1, 1, 0);
				seqs[n].id = n;
				size += seqs[n++].l_seq;
			}
		}
		if (size >= chunk_size && (n&1) == 0) break;
		// if (size >= chunk_size) {
		//	break;
		// }
	}	
	if (size == 0) { // test if the 2nd file is finished
		if (ks2 && kseq_read(ks2) >= 0)
			fprintf(stderr, "[W::%s] the 1st file has fewer sequences.\n", __func__);
	}
	*n_ = n;
	return seqs;
}

void bseq_classify(int n, bseq1_t *seqs, int m[2], bseq1_t *sep[2])
{
	int i, has_last;
	kvec_t(bseq1_t) a[2] = {{0,0,0}, {0,0,0}};
	for (i = 1, has_last = 1; i < n; ++i) {
		if (has_last) {
			if (strcmp(seqs[i].name, seqs[i-1].name) == 0) {
				kv_push(bseq1_t, a[1], seqs[i-1]);
				kv_push(bseq1_t, a[1], seqs[i]);
				has_last = 0;
			} else kv_push(bseq1_t, a[0], seqs[i-1]);
		} else has_last = 1;
	}
	if (has_last) kv_push(bseq1_t, a[0], seqs[i-1]);
	sep[0] = a[0].a, m[0] = a[0].n;
	sep[1] = a[1].a, m[1] = a[1].n;
}

/*****************
 * CIGAR related *
 *****************/

void bwa_fill_scmat(int a, int b, int8_t mat[25])
{
	int i, j, k;
	for (i = k = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j)
			mat[k++] = i == j? a : -b;
		mat[k++] = -1; // ambiguous base
	}
	for (j = 0; j < 5; ++j) mat[k++] = -1;
}

char getMethylationContext(int base1, int base2, int base3, int meth_status, int *cg_count, int *ch_count){
	if (base1 == 1 && base2 == 2){
		*cg_count += 1;
		if (meth_status) return 'X';
		return 'x';
	}
	else if (base1 == 1 && base3 == 2){
		*ch_count += 1;
		if (meth_status) return 'Y';
		return 'y';
	}
	else{
		*ch_count += 1;
		if (meth_status) return 'Z';
    	return 'z';
	}
}

// Generate CIGAR when the alignment end points are known
uint32_t *bwa_gen_cigar2(const int8_t mat[25], int o_del, int e_del, int o_ins, int e_ins, int w_, int64_t l_pac, const uint8_t *pac, const uint8_t *opac, int l_query, 
                         uint8_t *query, int64_t rb, int64_t re, int *score, int *n_cigar, int *NM, 
                         int64_t crick_l, uint8_t *oquery, int *cg_meth, int *cg_unmeth, int *ch_meth, int *ch_unmeth)
{
	uint32_t *cigar = 0;
	uint8_t tmp, *rseq;
	int i, reverse;
	int64_t rlen;
	kstring_t str;
	const char *int2base;
	if (n_cigar) *n_cigar = 0;
	if (NM) *NM = -1;
	if (l_query <= 0 || rb >= re || (rb < l_pac && re > l_pac)) return 0; // reject if negative length or bridging the forward and reverse strand
	reverse = rb > l_pac ? 1: 0; 

	rseq = bns_get_seq(l_pac, pac, rb, re, &rlen);
	if (re - rb != rlen) goto ret_gen_cigar; // possible if out of range
	if (rb >= l_pac) { // then reverse both query and rseq; this is to ensure indels to be placed at the leftmost position
		for (i = 0; i < l_query>>1; ++i)
			tmp = query[i], query[i] = query[l_query - 1 - i], query[l_query - 1 - i] = tmp;
		for (i = 0; i < rlen>>1; ++i)
			tmp = rseq[i], rseq[i] = rseq[rlen - 1 - i], rseq[rlen - 1 - i] = tmp;
	}
	if (l_query == re - rb && w_ == 0) { // no gap; no need to do DP
		// UPDATE: we come to this block now... FIXME: due to an issue in mem_reg2aln(), we never come to this block. This does not affect accuracy, but it hurts performance.
		if (n_cigar) {
			cigar = malloc(4);
			cigar[0] = l_query<<4 | 0;
			*n_cigar = 1;
		}
		for (i = 0, *score = 0; i < l_query; ++i)
			*score += mat[rseq[i]*5 + query[i]];
	} else {
		int w, max_gap, max_ins, max_del, min_w;
		// set the band-width
		max_ins = (int)((double)(((l_query+1)>>1) * mat[0] - o_ins) / e_ins + 1.);
		max_del = (int)((double)(((l_query+1)>>1) * mat[0] - o_del) / e_del + 1.);
		max_gap = max_ins > max_del? max_ins : max_del;
		max_gap = max_gap > 1? max_gap : 1;
		w = (max_gap + abs((int)rlen - l_query) + 1) >> 1;
		w = w < w_? w : w_;
		min_w = abs((int)rlen - l_query) + 3;
		w = w > min_w? w : min_w;
		// NW alignment
		if (bwa_verbose >= 4) {
			printf("* Global bandwidth: %d\n", w);
			printf("* Global ref:   "); for (i = 0; i < rlen; ++i) putchar("ACGTN"[(int)rseq[i]]); putchar('\n');
			printf("* Global query: "); for (i = 0; i < l_query; ++i) putchar("ACGTN"[(int)query[i]]); putchar('\n');
		}
		*score = ksw_global2(l_query, query, rlen, rseq, 5, mat, o_del, e_del, o_ins, e_ins, w, n_cigar, &cigar);
	}
	if (n_cigar) {// compute NM and MD
		int is_crick, ref_base, methyl_base;
		char meth_state;
		uint8_t *cseq;
		int64_t clen;
		// check crick status, crick pac len is 2/4 and 4/4 of pac 
		is_crick = 0;
		if (rb >= crick_l && rb < l_pac) is_crick = 1;
		else if (rb >= crick_l + l_pac) is_crick = 1;

		// set reference base by mapping strand
		ref_base = is_crick? 2: 1;
		methyl_base = is_crick? 0: 3; 

		if (rb > re) cseq = bns_get_seq(l_pac, opac, rb + 2, re - 2, &clen);
		else cseq = bns_get_seq(l_pac, opac, rb - 2, re + 2, &clen);
		if (rb >= l_pac) { // then reverse both query and rseq; this is to ensure indels to be placed at the leftmost position
			for (i = 0; i < clen>>1; ++i)
				tmp = cseq[i], cseq[i] = cseq[clen - 1 - i], cseq[clen - 1 - i] = tmp;
			for (i = 0; i < l_query>>1; ++i)
				tmp = oquery[i], oquery[i] = oquery[l_query - 1 - i], oquery[l_query - 1 - i] = tmp;
		}
		int k, x, y, u, n_mm = 0, n_gap = 0;
		int base1, base2, base3, meth_pos = 0;
		int ch_meth_count = 0, ch_unmeth_count = 0, cg_meth_count = 0, cg_unmeth_count = 0;
		kstring_t meth_str;
		meth_str.l = meth_str.m = 0; meth_str.s = 0;
		str.l = str.m = *n_cigar * 4; str.s = (char*)cigar; // append MD to CIGAR
		int2base = rb < l_pac? "ACGTN" : "TGCAN";
		for (k = 0, x = y = u = 0; k < *n_cigar; ++k) {
			int op, len;
			cigar = (uint32_t*)str.s;
			op  = cigar[k]&0xf, len = cigar[k]>>4;
			if (op == 0) { // match
				for (i = 0; i < len; ++i) {
					if (!cseq) break;
					if (cseq[y + i + 2] == ref_base){
						// select bases to get relative methylation context
						if (!reverse && !is_crick) base1 = cseq[y + i + 2], base2 = cseq[y +i + 3], base3 = cseq[y + i + 4]; // wc2t
						else if (!reverse && is_crick) base1 = cseq[y + i + 1], base2 = cseq[y + i + 2], base3 = cseq[y + i + 3]; // cg2a
						else if (reverse && !is_crick) base1 = cseq[y + i + 2], base2 = cseq[y +i + 1], base3 = cseq[y + i]; // wg2a
						else base1 = cseq[y + i + 3], base2 = cseq[y + i + 2], base3 = cseq[y + i + 1]; // cc2t
						// check methylation status
						if (oquery[x + i] == ref_base) {
							meth_state = getMethylationContext(base1, base2, base3, 1, &cg_meth_count, &ch_meth_count);
							if (meth_pos > 0) kputw(meth_pos, &meth_str);
							kputc(meth_state, &meth_str);
							meth_pos = 0, ++u;
						}
						else if (oquery[x + i] == methyl_base) {
							meth_state = getMethylationContext(base1, base2, base3, 0, &cg_unmeth_count, &ch_unmeth_count);
							if (meth_pos > 0) kputw(meth_pos, &meth_str);
							kputc(meth_state, &meth_str);
							meth_pos = 0, ++u;
						}
						else {
							kputw(u, &str);
							kputc(int2base[cseq[ y+ i +2]], &str);
							++n_mm; ++meth_pos, u = 0;
						}
					}
					else if (oquery[x + i] != cseq[y + i +2]) {
						kputw(u, &str);
						kputc(int2base[cseq[ y+ i +2]], &str);
						++n_mm; ++meth_pos, u = 0;
					} else ++u, ++meth_pos;
				}
				x += len; y += len;
			} else if (op == 2) { // deletion
				if (k > 0 && k < *n_cigar - 1) { // don't do the following if D is the first or the last CIGAR
					kputw(u, &str); kputc('^', &str);
					for (i = 0; i < len; ++i)
						kputc(int2base[cseq[y+i+2]], &str);
					u = 0; n_gap += len;
				}
				y += len;
			} else if (op == 1) x += len, n_gap += len, meth_pos += len; // insertion
		}
		kputw(u, &str); 
		if (meth_pos > 0) kputw(meth_pos, &meth_str); 
		kputsn("\tXB:Z:", 6, &str);
		kputsn(meth_str.s, meth_str.l, &str);
		kputc(0, &str);
		*NM = n_mm + n_gap;
		cigar = (uint32_t*)str.s;
		if (rb >= l_pac) // reverse back query
			for (i = 0; i < l_query>>1; ++i)
				tmp = oquery[i], oquery[i] = oquery[l_query - 1 - i], oquery[l_query - 1 - i] = tmp;
		*cg_meth = cg_meth_count; *cg_unmeth = cg_unmeth_count;
		*ch_meth = ch_meth_count; *ch_unmeth = ch_unmeth_count;
		free(meth_str.s);
		free(cseq);
	}
	if (rb >= l_pac) // reverse back query
		for (i = 0; i < l_query>>1; ++i)
			tmp = query[i], query[i] = query[l_query - 1 - i], query[l_query - 1 - i] = tmp;

ret_gen_cigar:
	free(rseq);
	return cigar;
}

uint32_t *bwa_gen_cigar(const int8_t mat[25], int q, int r, int w_, int64_t l_pac, const uint8_t *pac, const uint8_t *opac, int l_query, uint8_t *query, int64_t rb, int64_t re, int *score, int *n_cigar, int *NM, 
                        int64_t crick_l, uint8_t *oquery, int *cg_meth, int *cg_unmeth, int *ch_meth, int *ch_unmeth)
{
	return bwa_gen_cigar2(mat, q, r, q, r, w_, l_pac, pac, opac, l_query, query, rb, re, score, n_cigar, NM, crick_l, oquery, cg_meth, cg_unmeth, ch_meth, ch_unmeth);
}

/*********************
 * Full index reader *
 *********************/

char *bwa_idx_infer_prefix(const char *hint)
{
	char *prefix;
	int l_hint;
	FILE *fp;
	l_hint = strlen(hint);
	prefix = malloc(l_hint + 3 + 4 + 1);
	strcpy(prefix, hint);
	strcpy(prefix + l_hint, ".64.bwt");
	if ((fp = fopen(prefix, "rb")) != 0) {
		fclose(fp);
		prefix[l_hint + 3] = 0;
		return prefix;
	} else {
		strcpy(prefix + l_hint, ".bwt");
		if ((fp = fopen(prefix, "rb")) == 0) {
			free(prefix);
			return 0;
		} else {
			fclose(fp);
			prefix[l_hint] = 0;
			return prefix;
		}
	}
}

bwt_t *bwa_idx_load_bwt(const char *hint)
{
	char *tmp, *prefix;
	bwt_t *bwt;
	prefix = bwa_idx_infer_prefix(hint);
	if (prefix == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to locate the index files\n", __func__);
		return 0;
	}
	tmp = calloc(strlen(prefix) + 5, 1);
	strcat(strcpy(tmp, prefix), ".bwt"); // FM-index
	bwt = bwt_restore_bwt(tmp);
	strcat(strcpy(tmp, prefix), ".sa");  // partial suffix array (SA)
	bwt_restore_sa(tmp, bwt);
	free(tmp); free(prefix);
	return bwt;
}

bwaidx_t *bwa_idx_load_from_disk(const char *hint, int which)
{
	bwaidx_t *idx;
	char *prefix;
	prefix = bwa_idx_infer_prefix(hint);
	if (prefix == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to locate the index files\n", __func__);
		return 0;
	}
	idx = calloc(1, sizeof(bwaidx_t));
	if (which & BWA_IDX_BWT) idx->bwt = bwa_idx_load_bwt(hint);
	if (which & BWA_IDX_BNS) {
		int i, c;
		idx->bns = bns_restore(prefix);
		for (i = c = 0; i < idx->bns->n_seqs; ++i)
			if (idx->bns->anns[i].is_alt) ++c;
		if (bwa_verbose >= 3)
			fprintf(stderr, "[M::%s] read %d ALT contigs\n", __func__, c);
		if (which & BWA_IDX_PAC) {
			idx->pac = calloc(idx->bns->l_pac/4+1, 1);
			err_fread_noeof(idx->pac, 1, idx->bns->l_pac/4+1, idx->bns->fp_pac); // concatenated 2-bit encoded sequence
			err_fclose(idx->bns->fp_pac);
			idx->bns->fp_pac = 0;
			idx->opac = calloc(idx->bns->l_pac/4+1, 1);
			err_fread_noeof(idx->opac, 1, idx->bns->l_pac/4+1, idx->bns->fp_opac); // concatenated 2-bit encoded sequence
			err_fclose(idx->bns->fp_opac);
			idx->bns->fp_opac = 0;
		}
	}
	free(prefix);
	return idx;
}

bwaidx_t *bwa_idx_load(const char *hint, int which)
{
	return bwa_idx_load_from_disk(hint, which);
}

void bwa_idx_destroy(bwaidx_t *idx)
{
	if (idx == 0) return;
	if (idx->mem == 0) {
		if (idx->bwt) bwt_destroy(idx->bwt);
		if (idx->bns) bns_destroy(idx->bns);
		if (idx->pac) free(idx->pac);
		if (idx->opac) free(idx->opac);
	} else {
		free(idx->bwt); free(idx->bns->anns); free(idx->bns);
		if (!idx->is_shm) free(idx->mem);
	}
	free(idx);
}

int bwa_mem2idx(int64_t l_mem, uint8_t *mem, bwaidx_t *idx)
{
	int64_t k = 0, x;
	int i;

	// generate idx->bwt
	x = sizeof(bwt_t); idx->bwt = malloc(x); memcpy(idx->bwt, mem + k, x); k += x;
	x = idx->bwt->bwt_size * 4; idx->bwt->bwt = (uint32_t*)(mem + k); k += x;
	x = idx->bwt->n_sa * sizeof(bwtint_t); idx->bwt->sa = (bwtint_t*)(mem + k); k += x;

	// generate idx->bns and idx->pac
	x = sizeof(bntseq_t); idx->bns = malloc(x); memcpy(idx->bns, mem + k, x); k += x;
	x = idx->bns->n_holes * sizeof(bntamb1_t); idx->bns->ambs = (bntamb1_t*)(mem + k); k += x;
	x = idx->bns->n_seqs  * sizeof(bntann1_t); idx->bns->anns = malloc(x); memcpy(idx->bns->anns, mem + k, x); k += x;
	for (i = 0; i < idx->bns->n_seqs; ++i) {
		idx->bns->anns[i].name = (char*)(mem + k); k += strlen(idx->bns->anns[i].name) + 1;
		idx->bns->anns[i].anno = (char*)(mem + k); k += strlen(idx->bns->anns[i].anno) + 1;
	}
	idx->pac = (uint8_t*)(mem + k); k += idx->bns->l_pac/4+1;
	assert(k == l_mem);

	idx->l_mem = k; idx->mem = mem;
	return 0;
}

int bwa_idx2mem(bwaidx_t *idx)
{
	int i;
	int64_t k, x, tmp;
	uint8_t *mem;

	// copy idx->bwt
	x = idx->bwt->bwt_size * 4;
	mem = realloc(idx->bwt->bwt, sizeof(bwt_t) + x); idx->bwt->bwt = 0;
	memmove(mem + sizeof(bwt_t), mem, x);
	memcpy(mem, idx->bwt, sizeof(bwt_t)); k = sizeof(bwt_t) + x;
	x = idx->bwt->n_sa * sizeof(bwtint_t); mem = realloc(mem, k + x); memcpy(mem + k, idx->bwt->sa, x); k += x;
	free(idx->bwt->sa);
	free(idx->bwt); idx->bwt = 0;

	// copy idx->bns
	tmp = idx->bns->n_seqs * sizeof(bntann1_t) + idx->bns->n_holes * sizeof(bntamb1_t);
	for (i = 0; i < idx->bns->n_seqs; ++i) // compute the size of heap-allocated memory
		tmp += strlen(idx->bns->anns[i].name) + strlen(idx->bns->anns[i].anno) + 2;
	mem = realloc(mem, k + sizeof(bntseq_t) + tmp);
	x = sizeof(bntseq_t); memcpy(mem + k, idx->bns, x); k += x;
	x = idx->bns->n_holes * sizeof(bntamb1_t); memcpy(mem + k, idx->bns->ambs, x); k += x;
	free(idx->bns->ambs);
	x = idx->bns->n_seqs * sizeof(bntann1_t); memcpy(mem + k, idx->bns->anns, x); k += x;
	for (i = 0; i < idx->bns->n_seqs; ++i) {
		x = strlen(idx->bns->anns[i].name) + 1; memcpy(mem + k, idx->bns->anns[i].name, x); k += x;
		x = strlen(idx->bns->anns[i].anno) + 1; memcpy(mem + k, idx->bns->anns[i].anno, x); k += x;
		free(idx->bns->anns[i].name); free(idx->bns->anns[i].anno);
	}
	free(idx->bns->anns);

	// copy idx->pac
	x = idx->bns->l_pac/4+1;
	mem = realloc(mem, k + x);
	memcpy(mem + k, idx->pac, x); k += x;
	free(idx->bns); idx->bns = 0;
	free(idx->pac); idx->pac = 0;

	return bwa_mem2idx(k, mem, idx);
}

/***********************
 * SAM header routines *
 ***********************/

void bwa_print_sam_hdr(const bntseq_t *bns, const char *hdr_line)
{
	int i, n_SQ = 0;
	extern char *bwa_pg;
	if (hdr_line) {
		const char *p = hdr_line;
		while ((p = strstr(p, "@SQ\t")) != 0) {
			if (p == hdr_line || *(p-1) == '\n') ++n_SQ;
			p += 4;
		}
	}
	if (n_SQ == 0) {
		for (i = 0; i < bns->n_seqs; ++i) {
			if (! bns->anns[i].is_crick){
				err_printf("@SQ\tSN:%s\tLN:%d", bns->anns[i].name, bns->anns[i].len);
				if (bns->anns[i].is_alt) err_printf("\tAH:*\n");
				else err_fputc('\n', stdout);
			}
		}
	} else if (n_SQ != bns->n_seqs && bwa_verbose >= 2)
		fprintf(stderr, "[W::%s] %d @SQ lines provided with -H; %d sequences in the index. Continue anyway.\n", __func__, n_SQ, bns->n_seqs);
	if (hdr_line) err_printf("%s\n", hdr_line);
	if (bwa_pg) err_printf("%s\n", bwa_pg);
}

static char *bwa_escape(char *s)
{
	char *p, *q;
	for (p = q = s; *p; ++p) {
		if (*p == '\\') {
			++p;
			if (*p == 't') *q++ = '\t';
			else if (*p == 'n') *q++ = '\n';
			else if (*p == 'r') *q++ = '\r';
			else if (*p == '\\') *q++ = '\\';
		} else *q++ = *p;
	}
	*q = '\0';
	return s;
}

char *bwa_set_rg(const char *s)
{
	char *p, *q, *r, *rg_line = 0;
	memset(bwa_rg_id, 0, 256);
	if (strstr(s, "@RG") != s) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] the read group line is not started with @RG\n", __func__);
		goto err_set_rg;
	}
	if (strstr(s, "\t") != NULL) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] the read group line contained literal <tab> characters -- replace with escaped tabs: \\t\n", __func__);
		goto err_set_rg;
	}
	rg_line = strdup(s);
	bwa_escape(rg_line);
	if ((p = strstr(rg_line, "\tID:")) == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] no ID within the read group line\n", __func__);
		goto err_set_rg;
	}
	p += 4;
	for (q = p; *q && *q != '\t' && *q != '\n'; ++q);
	if (q - p + 1 > 256) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] @RG:ID is longer than 255 characters\n", __func__);
		goto err_set_rg;
	}
	for (q = p, r = bwa_rg_id; *q && *q != '\t' && *q != '\n'; ++q)
		*r++ = *q;
	return rg_line;

err_set_rg:
	free(rg_line);
	return 0;
}

char *bwa_insert_header(const char *s, char *hdr)
{
	int len = 0;
	if (s == 0 || s[0] != '@') return hdr;
	if (hdr) {
		len = strlen(hdr);
		hdr = realloc(hdr, len + strlen(s) + 2);
		hdr[len++] = '\n';
		strcpy(hdr + len, s);
	} else hdr = strdup(s);
	bwa_escape(hdr + len);
	return hdr;
}
