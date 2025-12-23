#!/usr/bin/env bash

## $1   EntityID
export PATH=$(dirname $0):$PATH
cd $PWD

cat $1.DataQC.tsv | perl -ne 'if(1..1){next}else{@l=split/\s+/;$l[3]=~s/,//g;$l[4]=~s/,//g;$num+=$l[3];$len+=$l[4];}END{print"Total\tFASTQ\tDNA\t$num\t$len\n"}' | digit_formatter.pl -k > tmp1.seqkit.QC.stat.tsv


cat $1.DataQC.tsv | sed 's/ \+/     /g' | cut -f 1,4,5,7,13,18,17 > tmp2.seqkit.QC.stat.tsv
cat tmp2.seqkit.QC.stat.tsv tmp1.seqkit.QC.stat.tsv > $1.seqkit.QC.stat.tsv
rm tmp2.seqkit.QC.stat.tsv tmp2.seqkit.QC.stat.tsv

export data=$(tail -1 $1.seqkit.QC.stat.tsv | cut -f5 | sed 's/,//g')

cat $1.sort.digestion.pairs.gz.stat | perl -ne 'BEGIN{print "type\tpairwise contact\tyield pairwise contacts per Gb\n"}{chomp;@l=split;$ratio=($l[1]/$ENV{data})*1000000000;$ratio=~s/\.\d+$//;print"$l[0]\t$l[1]\t$ratio\n"}' | digit_formatter.pl -k > $1.pairwise.stat.tsv

awk -v OFS="\t" 'BEGIN{print "Resolution(bp)\tRatio of bins over 1k contact(%)"}{if(/^Resolution/){print $2,$NF*100}}' $1.sort.valid.pairs.gz.resolution > $1.resolution.stat.tsv

cat $1.map.stat.tsv | digit_formatter.pl -k > $1.align.stat.tsv

cat $1.map.stat.tsv |awk -v OFS="\t" 'BEGIN{print "type\tread_cout\tread_ratio(%)\tbase_count\tbase_ratio(%)"}{if(NR==2){r=$2;b=$3}else if(NR>2){rr=$2/r*100;br=$3/b*100;print $1,$3,br,$2,rr}}' | digit_formatter.pl -k -f 2 > $1.align.stat.tsv
