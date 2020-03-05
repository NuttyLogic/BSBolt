#include <iostream>
#include <queue>
#include <string.h>
#include "bs_sorter.h"
#include "bwa.h"
#include "fastmap.h"
#include "kstring.h"
#include "kseq.h"





samSorter::samSorter(bseq1_t *first_read, ktp_aux_t *aux){
        group_1_score = 0;
        group_2_score = 0;
        wg2a = 0;
        wc2t = 0;
        cg2a = 0;
        cc2t = 0;
        unaligned = 0;
        bs_ambiguous = 0;
        reads_observed = 0;
        alignments_observed = 0;
        current_read_name = strdup(first_read->name);
        aux_opts = aux;
        samSorter::bankRead(first_read);
}

samSorter::~samSorter(){
    int read_group = samSorter::pickReadGroup();
    if (read_group == 2) samSorter::outputReads(0, true);
    else if (read_group == 1) samSorter::outputReads(1);
    else samSorter::outputReads(0);
    free(samSorter::current_read_name);
    samSorter::resetBuffer();
}

void samSorter::processRead(bseq1_t *read){
    char *current_read = samSorter::current_read_name;
    if ( strcmp(read->name, samSorter::current_read_name) == 0){
        samSorter::bankRead(read);
    } else{
        int read_group = samSorter::pickReadGroup();
        if (read_group == 2) samSorter::outputReads(0, true);
        else if (read_group == 1) samSorter::outputReads(1);
        else if (read_group == 0) samSorter::outputReads(0);
        free(current_read_name);
        samSorter::current_read_name = strdup(read->name);
        samSorter::resetBuffer();
        samSorter::bankRead(read);
    }
}

void samSorter::setUnmapped(bseq1_t *read){
    free(read->sam);
    kstring_t *str;
    str = (kstring_t*)calloc(1, sizeof(kstring_t));
    int l_name = strlen(read->name);
    ks_resize(str, str->l + read->l_seq + l_name + (read->qual? read->l_seq : 0) + 20);
    kputsn(read->name, l_name, str); 
    kputc('\t', str);
    
    if (read->paired){
        if (read->first) kputsn("77\t", 3, str);
        else kputsn("141\t", 3, str);
    }
    else kputsn("4\t", 2, str);
    kputsn("*\t0\t0\t*\t*\t0\t0\t", 14, str);
    
    int i, qb = 0, qe = read->l_seq;
    ks_resize(str, str->l + (qe - qb) + 1);
    for (i = qb; i < qe; ++i){
        str->s[str->l++] = read->oseq[i];
    }
    kputc('\t', str);
    if (read->qual) { 
        ks_resize(str, str->l + (qe - qb) + 1);
        for (i = qb; i < qe; ++i) str->s[str->l++] = read->qual[i];
        str->s[str->l] = 0;
    }
    kputsn("\tAS:i:0", 7, str); 
    kputsn("\tXO:Z:WC\n", 9, str);
    read->sam = strdup(str->s);
    free(str->s);
}

void samSorter::updateMappingStats(int mapped, bool bs_conflict){
    ++samSorter::alignments_observed;
    if (mapped == 0){ 
        ++samSorter::unaligned;
        if (bs_conflict) ++samSorter::bs_ambiguous;
    }
    else if (mapped == 1) ++samSorter::cg2a;
    else if (mapped == 2) ++samSorter::cc2t;
    else if (mapped == 3) ++samSorter::wg2a;
    else if (mapped == 4) ++samSorter::wc2t;
}

void samSorter::outputReads(int read_group, bool output_unmapped){
    ++samSorter::reads_observed;
    if (read_group){
        while (not samSorter::group_2.empty()){
            bseq1_t *read = samSorter::group_2.front();
            samSorter::updateMappingStats(read->mapped, read->bs_conflict);
            fputs(read->sam, samSorter::aux_opts->fp);
            samSorter::group_2.pop();
            samSorter::freeRead(read);
        }
    } else{
        while (not samSorter::group_1.empty()){
            bseq1_t *read = samSorter::group_1.front();
            if (output_unmapped){
                if (read->mapped) samSorter::setUnmapped(read);
                read->mapped = 0;
                read->bs_conflict = true;
            }
            samSorter::updateMappingStats(read->mapped, read->bs_conflict);
            fputs(read->sam, samSorter::aux_opts->fp);
            samSorter::group_1.pop();
            samSorter::freeRead(read);
        }
    }
}

void samSorter::bankRead(bseq1_t *read){
    if (read->read_group == 0){
        if (not read->bs_conflict) samSorter::group_1_score += read->alignment_score;
        samSorter::group_1.push(read);
    }
    else{
         if (not read->bs_conflict) samSorter::group_2_score += read->alignment_score;
          samSorter::group_2.push(read);
    }
}

int samSorter::pickReadGroup(){
    if (samSorter::group_1_score > samSorter::group_2_score) return 0;
    else if (samSorter::group_1_score < samSorter::group_2_score) return 1;
    return 2;
}

void samSorter::resetBuffer(){
    samSorter::group_1_score = 0;
    samSorter::group_2_score = 0;
    while (not samSorter::group_1.empty()){
        bseq1_t *read = samSorter::group_1.front();
        samSorter::group_1.pop();
        samSorter::freeRead(read);
    }
    while (not samSorter::group_2.empty()){
        bseq1_t *read = samSorter::group_2.front();
        samSorter::group_2.pop();
        samSorter::freeRead(read);
    }
}

void samSorter::freeRead(bseq1_t *read){
    free(read->name); free(read->comment);
    free(read->seq); free(read->qual);
	free(read->sam); free(read->oseq);
}

/*
int main(){
    bseq1_t *test;
    test = (bseq1_t*)calloc(1, sizeof(bseq1_t));
    test->name = strdup("name");
	test->comment = strdup("++++");
	test->seq = strdup("atcg");
    test->oseq = strdup("atcg");
	test->qual = strdup("atcg");
    test->sam = strdup("my read");
	test->l_seq = 4;
	test->first = true;
	test->read_group = 0;
	test->subsitution_pattern = 0;
	test->conversion = 1;
	test->bs_conflict = false;
	test->crick = false;
	test->mapped = true;
	test->paired = false;
	test->alignment_score = 0;
    ktp_aux_t *test_aux; 
    samSorter test_sort = samSorter(test, test_aux);
    std::cout << test_sort.group_1_score << std::endl;
    return 0;
}
*/