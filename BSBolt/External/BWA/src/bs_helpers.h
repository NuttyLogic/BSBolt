#include <string>

using namespace std;

void bsConversion(char sequence[], int substitution_pattern, int iter_length);
void reverseSeq(char sequence[], int sequence_length);
void complement(char sequence[], int sequence_length);
void reverseComplement(char sequence[], int sequence_length);
bool checkRname(char rname[]);
void formatCrickRname(char rname[]);
void reverseCigar(unsigned int cigar[], int cigar_length);
int getCrickMappingPos(unsigned cigar[], int cigar_length, int reference_length, int mapping_position);
bool formatSeqQual(bool is_crick, char sequence[], char qual[], bool is_reverse);
int countAlts(char *tag, bool is_crick);