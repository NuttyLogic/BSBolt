#include <string>
#include <cstring>
#include <iostream>
#include "bs_helpers.h"
#include "malloc_wrap.h"

using namespace std;

int countAlts(char *tag, int is_crick){
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

int countBase(char sequence[], char ref_base){
    float base_count = 0;
    for(int i=0; i< strlen(sequence); i++){
        if (sequence[i] == ref_base) ++base_count;
    }
    return base_count;
}

int assessConversion(char sequence1[], char sequence2[], int paired_end, float substitution_proportion){
    // check to see which substitution pattern will result in the fewest number of changes
    // if a high number of bases have to changed with both substitution patterns align both and sort dowstream
    float observed_bases = strlen(sequence1);
    float c_count = countBase(sequence1, 'C');
    float g_count = countBase(sequence1, 'G');
    if (paired_end){
        observed_bases += strlen(sequence2);
        // add count based on expected Illumina paired end library prep
        g_count += countBase(sequence2, 'C');
        c_count += countBase(sequence2, 'G');
    }
    float c_proportion = c_count / observed_bases;
    float g_proportion = g_count / observed_bases;
    float c_g_difference = c_proportion - g_proportion;
    if (c_g_difference < 0) c_g_difference = c_g_difference * -1.0;
    if (c_proportion > substitution_proportion && g_proportion > substitution_proportion) return 2;
    else if (c_proportion == g_proportion) return 2;
    else if (c_g_difference < 0.02) return 2;
    else if (c_proportion < g_proportion) return 0;
    else return 1;
}


int checkRname(char rname[]){
    if (strstr(rname, "_crick_bs")){
        return 1;
    }
    else{
        return 0;
    }
}

void formatCrickRname(char rname[]){
    rname[strlen(rname) - 9] = '\0';
}

int formatSeqQual(int is_crick, char sequence[], char qual[], int is_reverse){
    if (is_crick){
        return 1;
    }
    else{
        return 0;
    }
}
