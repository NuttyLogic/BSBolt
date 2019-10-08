#!/bin/bash

depbase=`echo art_illumina_src/art_illumina.o | sed 's|[^/]*$|.deps/&|;s|\.o$||'`;\
g++ -DHAVE_CONFIG_H -I.  -O3 -I/opt/local/include  -g -O2 -MT art_illumina_src/art_illumina.o -MD -MP -MF $depbase.Tpo -c -o art_illumina_src/art_illumina.o art_illumina_src/art_illumina.cpp &&\
mv -f $depbase.Tpo $depbase.Po
depbase=`echo art_illumina_src/art_qual_scale.o | sed 's|[^/]*$|.deps/&|;s|\.o$||'`;\
g++ -DHAVE_CONFIG_H -I.  -O3 -I/opt/local/include  -g -O2 -MT art_illumina_src/art_qual_scale.o -MD -MP -MF $depbase.Tpo -c -o art_illumina_src/art_qual_scale.o art_illumina_src/art_qual_scale.cpp &&\
mv -f $depbase.Tpo $depbase.Po
depbase=`echo art_illumina_src/empdist.o | sed 's|[^/]*$|.deps/&|;s|\.o$||'`;\
g++ -DHAVE_CONFIG_H -I.  -O3 -I/opt/local/include  -g -O2 -MT art_illumina_src/empdist.o -MD -MP -MF $depbase.Tpo -c -o art_illumina_src/empdist.o art_illumina_src/empdist.cpp &&\
mv -f $depbase.Tpo $depbase.Po
depbase=`echo art_illumina_src/readSeqFile.o | sed 's|[^/]*$|.deps/&|;s|\.o$||'`;\
g++ -DHAVE_CONFIG_H -I.  -O3 -I/opt/local/include  -g -O2 -MT art_illumina_src/readSeqFile.o -MD -MP -MF $depbase.Tpo -c -o art_illumina_src/readSeqFile.o art_illumina_src/readSeqFile.cpp &&\
mv -f $depbase.Tpo $depbase.Po
depbase=`echo art_illumina_src/seqRead.o | sed 's|[^/]*$|.deps/&|;s|\.o$||'`;\
g++ -DHAVE_CONFIG_H -I.  -O3 -I/opt/local/include  -g -O2 -MT art_illumina_src/seqRead.o -MD -MP -MF $depbase.Tpo -c -o art_illumina_src/seqRead.o art_illumina_src/seqRead.cpp &&\
mv -f $depbase.Tpo $depbase.Po
depbase=`echo art_illumina_src/samRead.o | sed 's|[^/]*$|.deps/&|;s|\.o$||'`;\
g++ -DHAVE_CONFIG_H -I.  -O3 -I/opt/local/include  -g -O2 -MT art_illumina_src/samRead.o -MD -MP -MF $depbase.Tpo -c -o art_illumina_src/samRead.o art_illumina_src/samRead.cpp &&\
mv -f $depbase.Tpo $depbase.Po
g++  -g -O2  -L/opt/local/lib -o art_illumina art_illumina_src/art_illumina.o art_illumina_src/art_qual_scale.o art_illumina_src/empdist.o art_illumina_src/readSeqFile.o art_illumina_src/seqRead.o art_illumina_src/samRead.o  -lgsl -lgslcblas -lm 
