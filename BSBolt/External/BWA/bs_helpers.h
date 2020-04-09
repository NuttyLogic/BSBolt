
//using namespace std;

#ifdef __cplusplus
  extern "C" {
#endif
void reverseSeq(char sequence[], int sequence_length);
void complement(char sequence[], int sequence_length);


void bsConversion(char sequence[], int substitution_pattern, int iter_length);
void reverseComplement(char sequence[], int sequence_length);
int checkRname(char rname[]);
void formatCrickRname(char rname[]);
void reverseCigar(unsigned int cigar[], int cigar_length);
int getCrickMappingPos(unsigned cigar[], int cigar_length, int reference_length, int mapping_position);
int formatSeqQual(int is_crick, char sequence[], char qual[], int is_reverse);
int countAlts(char *tag, int is_crick);
void swapSeq(char const ** str1, char const ** str2);
void swapCigar(const unsigned int ** str1, const unsigned int ** str2);

#ifdef __cplusplus
  }
#endif