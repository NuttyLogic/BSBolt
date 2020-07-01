
//using namespace std;

#ifdef __cplusplus
  extern "C" {
#endif

void bsConversion(char sequence[], int substitution_pattern, int iter_length);
int checkRname(char rname[]);
void formatCrickRname(char rname[]);
int formatSeqQual(int is_crick, char sequence[], char qual[], int is_reverse);
int countAlts(char *tag, int is_crick);
int assessConversion(char sequence1[], char sequence2[], int paired_end, float substitution_proportion);

#ifdef __cplusplus
  }
#endif