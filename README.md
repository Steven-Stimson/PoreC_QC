# This is a Quality Control workflow for Pore-C Alignmet

## Pore-C data statistic by seqkit
```
seqkit stat --all --type dna --out-file ~{EntityID}.DataQC.tsv pore-c.data1.fq.gz pore-c.data2.fq.gz pore-c.data3.fq.gz
```
## Align Pore-C data to reference genome
```
falign \
-enzyme ~{enzyme} \
-outfmt frag-bam \
-out ~{EntityID}.frag.bam \
~{refernece} \
pore-c.data1.fq.gz \
pore-c.data2.fq.gz \
pore-c.data3.fq.gz
```
## Quality Control (all scripts below can be found in `script` directory)
```
Falign.QC.py --frag_bam ~{EntityID}.frag.bam --seqkit_stat ~{EntityID}.DataQC.tsv --ID ~{EntityID}
order_read_type.py --input ~{EntityID}.order_read_type_read.stat --ID ~{EntityID}.order_read_type_read
order_read_type.py --input ~{EntityID}.order_read_type_base.stat --ID ~{EntityID}.order_read_type_base
order_ref_num.py --input ~{EntityID}.order_ref_num_read.stat --ID ~{EntityID}.order_ref_num_read
order_ref_num.py --input ~{EntityID}.order_ref_num_base.stat --ID ~{EntityID}.order_ref_num_base
read_len_distribution.py --input ~{EntityID}.read_len_distribution.tsv --ID ~{EntityID}.read_len_distribution
simulate_digestion.py ~{reference} -r $enzyme -o ~{EntityID}.digestion.frag.length -e 100 75 50 25 20
samtools view ~{Align}/~{EntityID}.frag.bam | perl -ne '/qs:i:(\d+)/ && (\$s=\$1); /qe:i:(\d+)/ && (\$e=\$1); \$l=\$e-\$s;print(\"~{EntityID}\t\$l\\n\")' >> ~{EntityID}.digestion.frag.length
length_by_window.py ~{EntityID}.digestion.frag.length -o ~{EntityID}.digestion.frag.window.stat -w 100
plot.digestion.py ~{EntityID}.digestion.frag.window.stat -o ~{EntityID}.digestion.frag.window --x_celling 10000
digest.genome.py ~{reference} -r ^GATC -o ~{EntityID}.digestion.bed
bam2pairs.py ~{Align}/~{EntityID}.frag.bam > ~{EntityID}.pairs
pairtools sort --nproc ~{req_cpu} --memory ~{req_memory} ~{EntityID}.pairs --output ~{EntityID}.sort.pairs.gz
check_valid_pairs.py ~{EntityID}.sort.pairs.gz ~{EntityID}.digestion.bed ~{EntityID}.sort.digestion.pairs.gz
count_valid_pairs.py ~{EntityID}.sort.add.pairs.gz > ~{EntityID}.sort.digestion.pairs.gz.stat
gunzip -cd ~{EntityID}.sort.digestion.pairs.gz | awk '/^#/ || \$9!=\$10' | cut -f1-8 |gzip > ~{EntityID}.sort.valid.pairs.gz
calculate_resolution_for_pairs.py ~{EntityID}.sort.valid.pairs.gz > ~{EntityID}.sort.valid.pairs.gz.resolution
split_inter_intra.py ~{EntityID}.sort.valid.pairs.gz ./
```
