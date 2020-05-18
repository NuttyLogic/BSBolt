#ifndef __bs_sorter_wrapper_H
#define __bs_sorter_wrapper_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct samSorter samSorter;

samSorter* new_samSorter(bseq1_t *first_read);
void samSorter_processRead(samSorter* ss, bseq1_t *read);
void samSorter_processEnd(samSorter* ss);
void samSorter_delete(samSorter* ss);

#ifdef __cplusplus
}
#endif
#endif