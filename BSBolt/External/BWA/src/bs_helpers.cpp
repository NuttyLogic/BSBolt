#include <string>
#include <cstring>
#include "bs_helpers.h"

using namespace std;
 


int countAlts(char *tag, bool is_crick){
    char strand = '-';
    if (is_crick) strand = '+';
    for(int i=0; i< strlen(tag); i++){
        if (tag[i] == strand) return 1;
    }
    return 0;
}

void bsConversion(char sequence[], int substitution_pattern, int iter_length){
    char ref_base = 'C';
    char sub_base = 'T';
    if (substitution_pattern){
        ref_base = 'G';
        sub_base = 'A';
    }
    for (int counter = 0; counter < iter_length; ++counter){
        if (sequence[counter] == ref_base){
            sequence[counter] = sub_base;
        }
    }
}

void reverseSeq(char sequence[], int sequence_length){
    for(int i=0; i< sequence_length/2; i++)
    {
        swap(sequence[i], sequence[sequence_length-i-1]);
    }
}

void complement(char sequence[], int sequence_length){
    for(int i=0; i< sequence_length; i++){
        if (sequence[i] == 'A') sequence[i] = 'T';
        else if (sequence[i] == 'T') sequence[i] = 'A';
        else if (sequence[i] == 'C') sequence[i] = 'G';
        else if (sequence[i] == 'G') sequence[i] = 'C';
        else sequence[i] = 'N';
    }
}

void reverseComplement(char sequence[], int sequence_length){
    reverseSeq(sequence, sequence_length);
    complement(sequence, sequence_length);
}

bool checkRname(char rname[]){
    if (strstr(rname, "_crick_bs")){
        return true;
    }
    else{
        return false;
    }
}

void formatCrickRname(char rname[]){
    rname[strlen(rname) - 9] = '\0';
}

void reverseCigar(unsigned int cigar[], int cigar_length){
    for(int i=0; i< cigar_length/2; i++)
    {
        swap(cigar[i], cigar[cigar_length-i-1]);
    }
}

int getCrickMappingPos(unsigned int cigar[], int cigar_length, int reference_length, int mapping_position){
    int mapping_length = 0;
    for (int i = 0; i < cigar_length; ++i) {
		int c = cigar[i]&0xf;
        //reference_consumers M and D (0, 2)
        if (c == 0 or c == 2){
            mapping_length += cigar[i]>>4;
        }
    }
    int watson_mapping_pos = reference_length - mapping_position - mapping_length;
    return watson_mapping_pos;
}

bool formatSeqQual(bool is_crick, char sequence[], char qual[], bool is_reverse){
    if (is_crick){
        if (not is_reverse){ 
            reverseComplement(sequence, strlen(sequence));
            reverseSeq(qual, strlen(qual));
        }
        return true;
    }
    else{
        if (is_reverse){
            reverseComplement(sequence, strlen(sequence));
            reverseSeq(qual, strlen(qual));
        }
    return false;
    }
}