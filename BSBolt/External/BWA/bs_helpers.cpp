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
