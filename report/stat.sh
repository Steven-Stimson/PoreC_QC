# To generate statistic for HTML report

# $1    EntityID
# $2    raw seqkit stat
# $3    raw pair wise stat
# $4    raw resolution stat

cat $2 | perl -ne 'if(1..1){print}else{@l=split/\s+/;$l[3]=~s/,//g;$l[4]=~s/,//g;$num+=$l[3];$len+=$l[4];print}END{print"Total\tFASTQ\tDNA\t$num\t$len\n"}' > $1.seqkit.QC.stat.tsv

export data=$(tail -1 $1.seqkit.QC.stat.tsv |cut -f5)

cat $3 | perl -ne 'BEGIN{print "type\tpairwise contact (M)\tyield pairwise contacts per Gb (M)\n"}{chomp;@l=split;$l[1]/=1000000;$ratio=($l[1]/$ENV{data})*1000000000;$l[1]=sprintf("%.2f",$l[1]);$ratio=sprintf("%.2f",$ratio);print"$l[0]\t$l[1]\t$ratio\n"}' > $1.pairwise.stat.tsv

awk -v OFS="\t" 'BEGIN{print "Resolution(bp)\tRatio of bins over 1k contact(%)"}{if(/^Resolution/){print $2,$NF*100}}' $4 > $1.resolution.stat.tsv
