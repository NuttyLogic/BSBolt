#include "bs_sorter.h"
#include "bs_sorter_wrapper.h"

extern "C" {
        samSorter* new_samSorter(bseq1_t *first_read) {
                return new samSorter(first_read);
        }

        void samSorter_processRead(samSorter* ss, bseq1_t *read) {
            ss->processRead(read);
        }

        void samSorter_processEnd(samSorter* ss) {
                ss->processEnd();
                fprintf(stderr, "BSStat TotalReads: %d\n", ss->reads_observed);
		        fprintf(stderr, "BSStat TotalAlignments: %d\n", ss->alignments_observed);
		        fprintf(stderr, "BSStat W_C2T: %d\n", ss->wc2t);
		        fprintf(stderr, "BSStat W_G2A: %d\n", ss->wg2a);
		        fprintf(stderr, "BSStat C_C2T: %d\n", ss->cc2t);
		        fprintf(stderr, "BSStat C_G2A: %d\n", ss->cg2a);
		        fprintf(stderr, "BSStat Unaligned: %d\n", ss->unaligned);
		        fprintf(stderr, "BSStat BSAmbiguous: %d\n", ss->bs_ambiguous);
        }

        void samSorter_delete(samSorter* ss) {
                delete ss;
        }
}