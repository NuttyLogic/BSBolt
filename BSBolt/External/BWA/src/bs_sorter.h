#include <queue>
#include <stdlib.h>
#include <string.h>
#include "bwa.h"
#include "fastmap.h"
#include "kseq.h"
#include "kstring.h"


using namespace std;

class samSorter
{
    public:
    std::queue<bseq1_t*> group_1;
    std::queue<bseq1_t*> group_2;
    int group_1_score, group_2_score, wg2a, wc2t, cg2a, cc2t; 
    int unaligned, bs_ambiguous, reads_observed, alignments_observed;
    char *current_read_name;
    ktp_aux_t *aux_opts;

    samSorter(bseq1_t *first_read, ktp_aux_t *aux);
    ~samSorter();
    void processRead(bseq1_t *read);
    void processEnd();
    void bankRead(bseq1_t *read);
    void outputReads(int read_group, bool output_unmapped=false);
    int pickReadGroup();
    void resetBuffer();
    void freeRead(bseq1_t *read);
    void setUnmapped(bseq1_t *read);
    void updateMappingStats(int mapping_group, bool bs_conflict);
};