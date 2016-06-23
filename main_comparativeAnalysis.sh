#!/bin/sh

## construction of the basic diretory structure
git clone https://github.com/drtamermansour/p_asteroides_CompAna.git
cd p_asteroides_CompAna
p_asteroides=$(pwd)

## create working directory and define paths for raw data and scripts
mkdir -p $p_asteroides/{c_abundFilter,compAnalysis,resources,blast_out}
script_path=$p_asteroides/scripts
abundFilter=$p_asteroides/c_abundFilter   ## copy the filtered reads from the p_asteroides assembly project
compAnalysis=$p_asteroides/compAnalysis   
## copy the necessary files from the p_asteroides assembly project
FCS_ann_exp_tran=$p_asteroides/compAnalysis/Trinity.clean.201.exp.FCS.fasta
swiss_ann=$p_asteroides/compAnalysis/uniprot_sprot.blastx.outfmt6.sig.best.exp2.FCS
LongOrfs=$p_asteroides/compAnalysis/longest_orfs.pep.exp2.FCS
gene_transcript_map=$p_asteroides/compAnalysis/gene_trans_map.3
########################
## blast out transcriptome aganist other known coral and Cnidarian sequences
bash prepResources.sh "${p_asteroides}"

module load BLAST+/2.2.30
cd ${p_asteroides}/resources/symb_trans
makeblastdb -in symb_trans.fasta -input_type fasta -dbtype nucl
blastn -query $FCS_ann_exp_tran \
       -db ${p_asteroides}/resources/symb_trans/symb_trans.fasta \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp" \
       -dust 'yes' -best_hit_overhang 0.25 -best_hit_score_edge 0.25 \
       -max_target_seqs 10  -out ${p_asteroides}/blast_out/trinityVsSymb ## -perc_identity 90 -qcov_hsp_perc 50
sort -k1,1 -k11,11g ${p_asteroides}/blast_out/trinityVsSymb | sort -u -k1,1 --merge > ${p_asteroides}/blast_out/trinityVsSymb_best

cd ${p_asteroides}/blast_out
grep "^>" $FCS_ann_exp_tran | wc -l        ## 867255
wc -l trinityVsSymb_best ## 188089
cat trinityVsSymb_best | awk '$11 <= 1e-5' > trinityVsSymb_best.sig
wc -l trinityVsSymb_best.sig ## 186177
cat trinityVsSymb_best.sig | awk -F '[\t.]' '{A[$2]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > symb.mapping_rates

module load QIIME/1.8.0
filter_fasta.py --input_fasta_fp $FCS_ann_exp_tran --output_fasta_fp trinityVsSymb_best.sig.fasta --seq_id_fp trinityVsSymb_best.sig
grep "^>" trinityVsSymb_best.sig.fasta | wc -l  ## 186177
zoox_trans=${p_asteroides}/blast_out/trinityVsSymb_best.sig.fasta

filter_fasta.py --input_fasta_fp $FCS_ann_exp_tran --output_fasta_fp noSymb.fasta --seq_id_fp trinityVsSymb_best.sig --negate
grep "^>" noSymb.fasta | wc -l  ## 681078
noSymb_transcriptome=${p_asteroides}/blast_out/noSymb.fasta

## select annotation for significant hits
cat trinityVsSymb_best.sig | awk '{print ">"$1"|"}' | grep -F -f - $LongOrfs | sed 's/>//' > LongOrfs.trinityVsSymb.key
filter_fasta.py --input_fasta_fp $LongOrfs --output_fasta_fp $LongOrfs.trinityVsSymb --seq_id_fp LongOrfs.trinityVsSymb.key
filter_fasta.py --input_fasta_fp $LongOrfs.complete --output_fasta_fp $LongOrfs.trinityVsSymb.complete --seq_id_fp LongOrfs.trinityVsSymb.key
grep "^>" $LongOrfs.trinityVsSymb | wc -l ## 221299 ## no of all possible ORF
grep "^>" $LongOrfs.trinityVsSymb | awk -F '[>|]' '{print $2}' | sort | uniq | wc -l ## 89425 ## no of transcripts with ORF
grep "^>" $LongOrfs.trinityVsSymb.complete | wc -l ## 124810 ## no of all possible complete ORF
grep "^>" $LongOrfs.trinityVsSymb.complete | awk -F '[>|]' '{print $2}' | sort | uniq | wc -l ## 59126 ## no of transcripts with complete ORF
cat trinityVsSymb_best.sig | awk '{print $1}' | grep -w -F -f - $swiss_ann > $swiss_ann.trinityVsSymb
cat $swiss_ann.trinityVsSymb | awk '{print $1}' | sort | uniq | wc -l
######################################
cd ${p_asteroides}/resources/coral_trans
makeblastdb -in coral_trans.fasta -input_type fasta -dbtype nucl
blastn -query $noSymb_transcriptome \
   -db ${p_asteroides}/resources/coral_trans/coral_trans.fasta \
   -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp" \
   -dust 'yes' -best_hit_overhang 0.25 -best_hit_score_edge 0.25 \
   -max_target_seqs 10  -out ${p_asteroides}/blast_out/trinityVscoralTrans ## -perc_identity 90 -qcov_hsp_perc 50
sort -k1,1 -k11,11g ${p_asteroides}/blast_out/trinityVscoralTrans | sort -u -k1,1 --merge > ${p_asteroides}/blast_out/trinityVscoralTrans_best

cd ${p_asteroides}/blast_out
wc -l trinityVscoralTrans_best ## 92944
cat trinityVscoralTrans_best | awk '$11 <= 1e-5' > trinityVscoralTrans_best.sig
wc -l trinityVscoralTrans_best.sig ## 92332
cat trinityVscoralTrans_best.sig | awk -F '[\t.]' '{A[$2]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > coralTrans.mapping_rates

filter_fasta.py --input_fasta_fp $noSymb_transcriptome --output_fasta_fp trinityVscoralTrans_best.sig.fasta --seq_id_fp trinityVscoralTrans_best.sig
grep "^>" trinityVscoralTrans_best.sig.fasta | wc -l  ## 92332

filter_fasta.py --input_fasta_fp $noSymb_transcriptome --output_fasta_fp noSymb_noCoralTrans.fasta --seq_id_fp trinityVscoralTrans_best.sig --negate
grep "^>" noSymb_noCoralTrans.fasta | wc -l  ## 588746
noSymb_noCoralTrans_transcriptome=${p_asteroides}/blast_out/noSymb_noCoralTrans.fasta

## select annotation for significant hits
#cat trinityVscoralTrans_best.sig | awk '{print ">"$1"|"}' | grep -F -f - $LongOrfs | sed 's/>//' > LongOrfs.trinityVscoralTrans.key
#filter_fasta.py --input_fasta_fp $LongOrfs --output_fasta_fp $LongOrfs.trinityVscoralTrans --seq_id_fp LongOrfs.trinityVscoralTrans.key
#filter_fasta.py --input_fasta_fp $LongOrfs.complete --output_fasta_fp $LongOrfs.trinityVscoralTrans.complete --seq_id_fp LongOrfs.trinityVscoralTrans.key
#cat trinityVscoralTrans_best.sig | awk '{print $1}' | grep -w -F -f - $swiss_ann > $swiss_ann.trinityVscoralTrans
######################################
cd ${p_asteroides}/resources/coral_genomic
makeblastdb -in coral_genomic.fasta -input_type fasta -dbtype nucl
blastn -query $noSymb_noCoralTrans_transcriptome \
   -db ${p_asteroides}/resources/coral_genomic/coral_genomic.fasta \
   -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp" \
   -dust 'yes' -best_hit_overhang 0.25 -best_hit_score_edge 0.25 \
   -max_target_seqs 10  -out ${p_asteroides}/blast_out/trinityVscoralGenomic ## -perc_identity 90 -qcov_hsp_perc 50
sort -k1,1 -k11,11g ${p_asteroides}/blast_out/trinityVscoralGenomic | sort -u -k1,1 --merge > ${p_asteroides}/blast_out/trinityVscoralGenomic_best

cd ${p_asteroides}/blast_out
wc -l trinityVscoralGenomic_best ## 38665
cat trinityVscoralGenomic_best | awk '$11 <= 1e-5' > trinityVscoralGenomic_best.sig
wc -l trinityVscoralGenomic_best.sig ## 37414
cat trinityVscoralGenomic_best.sig | awk -F '[\t.]' '{A[$2]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > coralGenomics.mapping_rates

filter_fasta.py --input_fasta_fp $noSymb_noCoralTrans_transcriptome --output_fasta_fp trinityVscoralGenomic_best.sig.fasta --seq_id_fp trinityVscoralGenomic_best.sig
grep "^>" trinityVscoralGenomic_best.sig.fasta | wc -l  ## 37414

filter_fasta.py --input_fasta_fp $noSymb_noCoralTrans_transcriptome --output_fasta_fp unrecognized.fasta --seq_id_fp trinityVscoralGenomic_best.sig --negate
grep "^>" unrecognized.fasta | wc -l  ## 551332
unrecognized=${p_asteroides}/blast_out/unrecognized.fasta

cat trinityVscoralTrans_best.sig trinityVscoralGenomic_best.sig > trinityVscoral_best.sig ## 129746
cat trinityVsSymb_best.sig trinityVscoral_best.sig > recognized_best.sig
cat trinityVscoralTrans_best.sig.fasta trinityVscoralGenomic_best.sig.fasta > coral_transcriptome.fasta
cat trinityVsSymb_best.sig.fasta coral_transcriptome.fasta > recognized.fasta

coral_transcriptome=${p_asteroides}/blast_out/coral_transcriptome.fasta
recognized_transcriptome=${p_asteroides}/blast_out/recognized.fasta

## select annotation for significant hits
#cat trinityVscoralGenomic_best.sig | awk '{print ">"$1"|"}' | grep -F -f - $LongOrfs | sed 's/>//' > LongOrfs.trinityVscoralGenomic.key
#filter_fasta.py --input_fasta_fp $LongOrfs --output_fasta_fp $LongOrfs.trinityVscoralGenomic --seq_id_fp LongOrfs.trinityVscoralGenomic.key
#filter_fasta.py --input_fasta_fp $LongOrfs.complete --output_fasta_fp $LongOrfs.trinityVscoralGenomic.complete --seq_id_fp LongOrfs.trinityVscoralGenomic.key
#cat trinityVscoralGenomic_best.sig | awk '{print $1}' | grep -w -F -f - $swiss_ann > $swiss_ann.trinityVscoralGenomic

#cat $LongOrfs.trinityVscoralTrans $LongOrfs.trinityVscoralGenomic > $LongOrfs.trinityVscoral
#cat $LongOrfs.trinityVsSymb $LongOrfs.trinityVscoral > $LongOrfs.recognized
#cat $LongOrfs.trinityVscoralTrans.complete $LongOrfs.trinityVscoralGenomic.complete > $LongOrfs.trinityVscoral.complete
#cat $LongOrfs.trinityVsSymb.complete $LongOrfs.trinityVscoral.complete > $LongOrfs.recognized.complete
#cat $swiss_ann.trinityVscoralTrans $swiss_ann.trinityVscoralGenomic > $swiss_ann.trinityVscoral
#cat $swiss_ann.trinityVsSymb $swiss_ann.trinityVscoral > $swiss_ann.recognized
##################
## Abundance estimation (This general expression is not used anymore)
cd $compAnalysis
qsub -v index="salmon_index",transcriptome="$FCS_ann_exp_tran" ${script_path}/salmonIndex.sh

cd $abundFilter
for f in $p_asteroides/c_abundFilter/*.s_pe.fq.1; do if [ -f $f ]; then
 identifier=$(basename ${f%_L00*.s_pe.fq.1}); echo $identifier;
fi;done | uniq > identifiers.txt
identifiers=$p_asteroides/c_abundFilter/identifiers.txt

while read identifier;do
 ls ${identifier}_L00*.s_pe.fq.1 ${identifier}_L00*.s_pe2.fq ${identifier}_L00*.s_se.fq
 qsub -v index="$compAnalysis/salmon_index",identifier=$identifier ${script_path}/salmonQuant_PE.sh
 qsub -v index="$compAnalysis/salmon_index",identifier=$identifier ${script_path}/salmonQuant_SE.sh
done < $identifiers
find ./*.quant -name *.sf -exec grep -H "mapping rate" {} \; | sort > salmonQuant_summary.txt
python $script_path/gather-counts.py -i "$(pwd)"
echo "transcript"$'\t'"length" > transcripts.lengthes
sf=$(find ./*.quant -name \*.sf | head -n1)
cat $sf | grep -v "^#" | awk -F "\t" -v OFS='\t' '{print $1,$2}' >> transcripts.lengthes

#module load R/3.0.1
while read identifier;do
  echo $(pwd) $identifier
  bash $script_path/run_calcTPM.sh "$(pwd)" "$identifier" "transcripts.lengthes" "$gene_transcript_map" ${script_path}/calcTPM2.R
  #Rscript ${script_path}/calcTPM2.R "$(pwd)" "$identifier" "transcripts.lengthes" "$gene_transcript_map" >> targets_list
done < $identifiers
bash $script_path/abund_est.sh  ## This produce (for gene & isoform) allTissues_TPM, adultOnly_TPM, larva_TPM, larvaOnly_TPM, unexp_TPM, exp_TPM

## expression stastics 
cd ${p_asteroides}/blast_out
## isoforms expressed in larva (568353) & has sig hit against symb. (186177)
comm -12 <(cat $abundFilter/larva_isoformTPM | tail -n+2 | awk '{print $1}' |sort) <(cat $p_asteroides/blast_out/trinityVsSymb_best.sig | awk '{print $1}' | sort) > larva_isoformTPM_vs_trinityVsSymb  ## 179515
comm -12 <(cat $swiss_ann | tail -n+2 | awk '{print $1}' |sort) larva_isoformTPM_vs_trinityVsSymb > larva_isoformTPM_vs_trinityVsSymb.ann  ## 49780

## isoforms expressed in larva (568353) & has sig hit against coral seq (129746)
comm -12 <(cat $abundFilter/larva_isoformTPM | tail -n+2 | awk '{print $1}' |sort) <(cat $p_asteroides/blast_out/trinityVscoral_best.sig | awk '{print $1}' | sort) > larva_isoformTPM_vs_trinityVscoral  ## 125873
comm -12 <(cat $swiss_ann | tail -n+2 | awk '{print $1}' |sort) larva_isoformTPM_vs_trinityVscoral > larva_isoformTPM_vs_trinityVscoral.ann  ## 24248

## isoforms expressed in adult only (298617) & has sig hit against symb. (186177)
comm -12 <(cat $abundFilter/adultOnly_isoformTPM | tail -n+2 | awk '{print $1}' |sort) <(cat $p_asteroides/blast_out/trinityVsSymb_best.sig | awk '{print $1}' | sort) > adultOnly_isoformTPM_vs_trinityVsSymb  ## 6615
comm -12 <(cat $swiss_ann | tail -n+2 | awk '{print $1}' |sort) adultOnly_isoformTPM_vs_trinityVsSymb > adultOnly_isoformTPM_vs_trinityVsSymb.ann  ## 2518

## isoforms expressed in adult only (298617)  & has sig hit against coral seq (129746)
comm -12 <(cat $abundFilter/adultOnly_isoformTPM | tail -n+2 | awk '{print $1}' |sort) <(cat $p_asteroides/blast_out/trinityVscoral_best.sig | awk '{print $1}' | sort) > adultOnly_isoformTPM_vs_trinityVscoral  ## 3845
comm -12 <(cat $swiss_ann | tail -n+2 | awk '{print $1}' |sort) adultOnly_isoformTPM_vs_trinityVscoral > adultOnly_isoformTPM_vs_trinityVscoral.ann  ## 1136
#####
## strain specific transcriptome abundance and annotation
cd $abundFilter
mkdir coral
cat ${p_asteroides}/blast_out/trinityVscoral_best.sig | awk '{print $1}' > coral_transIDs
head -n1 transcripts.lengthes > coral/transcripts.lengthes
grep -w -F -f coral_transIDs transcripts.lengthes >> coral/transcripts.lengthes
grep -w -F -f coral_transIDs $gene_transcript_map > coral/gene_transcript_map
for f in *.quant.counts;do
  head -n1 $f > coral/$f
  grep -w -F -f coral_transIDs $f >> coral/$f
done
cd coral

#module load R/3.0.1
while read identifier;do
  echo $(pwd) $identifier
  bash $script_path/run_calcTPM.sh "$(pwd)" "$identifier" "transcripts.lengthes" "gene_transcript_map" ${script_path}/calcTPM2.R
  #Rscript ${script_path}/calcTPM2.R "$(pwd)" "$identifier" "transcripts.lengthes" "gene_transcript_map" >> targets_list
done < $identifiers
bash $script_path/abund_est.sh

## define the species specific transcriptome and related annotations 
filter_fasta.py --input_fasta_fp $FCS_ann_exp_tran --output_fasta_fp p_ast2016.fasta --seq_id_fp exp_isoformTPM
p_ast2016_trans=$abundFilter/coral/p_ast2016.fasta ## This file should be matching to "coral_transcriptome.fasta" but missing few transcripts with TPM=0
Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1], header=T,row.names=NULL);head(data1);'\
'data2=read.table(args[2], header=T,row.names=NULL,sep="\t",quote="");head(data2);'\
'data3=read.table(args[3], header=F,row.names=NULL,sep="\t");colnames(data3)[2]="coral_blast_hit";head(data3);'\
'ann_isoExp=merge(data1,data2[,c(1,2,17,11)],by.x="geneName",by.y="qseqid",all.x = TRUE);'\
'ann_isoExp2=merge(ann_isoExp,data3[,c(1,2)],by.x="geneName",by.y="V1",all.x = TRUE);'\
'write.table(ann_isoExp2,"ann_isoExp", sep="\t", quote=F, row.names=F, col.names=T);' exp_isoformTPM $swiss_ann ${p_asteroides}/blast_out/trinityVscoral_best.sig
tail -n+2 exp_isoformTPM | awk '{print $1"|"}' | grep -F -f - $LongOrfs | sed 's/>//' > LongOrfs.coral.key
filter_fasta.py --input_fasta_fp $LongOrfs --output_fasta_fp LongOrfs.coral.pep --seq_id_fp LongOrfs.coral.key
filter_fasta.py --input_fasta_fp $LongOrfs.complete --output_fasta_fp LongOrfs.coral.pep.complete --seq_id_fp LongOrfs.coral.key
grep "^>" LongOrfs.coral.pep | wc -l ## 49348 ## no of all possible ORF
grep "^>" LongOrfs.coral.pep | awk -F '[>|]' '{print $2}' | sort | uniq | wc -l ## 33331 ## no of transcripts with ORF
grep "^>" LongOrfs.coral.pep.complete | wc -l ## 26802 ## no of all possible complete ORF
grep "^>" LongOrfs.coral.pep.complete | awk -F '[>|]' '{print $2}' | sort | uniq | wc -l ## 17277 ## no of transcripts with complete ORF
cat exp_isoformTPM | awk '{print $1}' | grep -w -F -f - $swiss_ann > $swiss_ann.p_ast2016
cat $swiss_ann.p_ast2016 | awk '{print $1}' | sort | uniq | wc -l
##############
cd $abundFilter
mkdir spC15
grep "S_spC15.est" ${p_asteroides}/blast_out/trinityVsSymb_best.sig > spC15/trinityVsSymb_best.sig.spC15
cat spC15/trinityVsSymb_best.sig.spC15 | awk '{print $1}' > spC15_transIDs
head -n1 transcripts.lengthes > spC15/transcripts.lengthes
grep -w -F -f spC15_transIDs transcripts.lengthes >> spC15/transcripts.lengthes
grep -w -F -f spC15_transIDs $gene_transcript_map > spC15/gene_transcript_map
for f in *.quant.counts;do
  head -n1 $f > spC15/$f
  grep -w -F -f spC15_transIDs $f >> spC15/$f
done
cd spC15

#module load R/3.0.1
while read identifier;do
  echo $(pwd) $identifier
  bash $script_path/run_calcTPM.sh "$(pwd)" "$identifier" "transcripts.lengthes" "gene_transcript_map" ${script_path}/calcTPM2.R
  #Rscript ${script_path}/calcTPM2.R "$(pwd)" "$identifier" "transcripts.lengthes" "gene_transcript_map" >> targets_list
done < $identifiers
bash $script_path/abund_est.sh

## define the species specific transcriptome and related annotations
filter_fasta.py --input_fasta_fp $FCS_ann_exp_tran --output_fasta_fp spC15_2016.fasta --seq_id_fp exp_isoformTPM
spC15_2016_trans=$abundFilter/spC15/spC15_2016.fasta
Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1], header=T,row.names=NULL);head(data1);'\
'data2=read.table(args[2], header=T,row.names=NULL,sep="\t",quote="");head(data2);'\
'data3=read.table(args[3], header=F,row.names=NULL,sep="\t");colnames(data3)[2]="coral_blast_hit";head(data3);'\
'ann_isoExp=merge(data1,data2[,c(1,2,17,11)],by.x="geneName",by.y="qseqid",all.x = TRUE);'\
'ann_isoExp2=merge(ann_isoExp,data3[,c(1,2)],by.x="geneName",by.y="V1",all.x = TRUE);'\
'write.table(ann_isoExp2,"ann_isoExp", sep="\t", quote=F, row.names=F, col.names=T);' exp_isoformTPM $swiss_ann trinityVsSymb_best.sig.spC15
tail -n+2 exp_isoformTPM | awk '{print $1"|"}' | grep -F -f - $LongOrfs | sed 's/>//' > LongOrfs.spC15.key
filter_fasta.py --input_fasta_fp $LongOrfs --output_fasta_fp LongOrfs.spC15.pep --seq_id_fp LongOrfs.spC15.key
filter_fasta.py --input_fasta_fp $LongOrfs.complete --output_fasta_fp LongOrfs.spC15.pep.complete --seq_id_fp LongOrfs.spC15.key
grep "^>" LongOrfs.spC15.pep | wc -l ## 47310 ## no of all possible ORF
grep "^>" LongOrfs.spC15.pep | awk -F '[>|]' '{print $2}' | sort | uniq | wc -l ## 29166 ## no of transcripts with ORF
grep "^>" LongOrfs.spC15.pep.complete | wc -l ## 34480 ## no of all possible complete ORF
grep "^>" LongOrfs.spC15.pep.complete | awk -F '[>|]' '{print $2}' | sort | uniq | wc -l ## 20741 ## no of transcripts with complete ORF
cat exp_isoformTPM | awk '{print $1}' | grep -w -F -f - $swiss_ann > $swiss_ann.spC15
cat $swiss_ann.spC15 | awk '{print $1}' | sort | uniq | wc -l
##############
cd $abundFilter
mkdir cladeA
grep "S_cladeA.est" ${p_asteroides}/blast_out/trinityVsSymb_best.sig > cladeA/trinityVsSymb_best.sig.cladeA
cat cladeA/trinityVsSymb_best.sig.cladeA | awk '{print $1}' > cladeA_transIDs
head -n1 transcripts.lengthes > cladeA/transcripts.lengthes
grep -w -F -f cladeA_transIDs transcripts.lengthes >> cladeA/transcripts.lengthes
grep -w -F -f cladeA_transIDs $gene_transcript_map > cladeA/gene_transcript_map
for f in *.quant.counts;do
  head -n1 $f > cladeA/$f
  grep -w -F -f cladeA_transIDs $f >> cladeA/$f
done
cd cladeA

#module load R/3.0.1
while read identifier;do
  echo $(pwd) $identifier
  bash $script_path/run_calcTPM.sh "$(pwd)" "$identifier" "transcripts.lengthes" "gene_transcript_map" ${script_path}/calcTPM2.R
  #Rscript ${script_path}/calcTPM2.R "$(pwd)" "$identifier" "transcripts.lengthes" "gene_transcript_map" >> targets_list
done < $identifiers
bash $script_path/abund_est.sh

## define the species specific transcriptome and related annotations
filter_fasta.py --input_fasta_fp $FCS_ann_exp_tran --output_fasta_fp cladeA_2016.fasta --seq_id_fp exp_isoformTPM
cladeA_2016_trans=$abundFilter/cladeA/cladeA_2016.fasta
Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1], header=T,row.names=NULL);head(data1);'\
'data2=read.table(args[2], header=T,row.names=NULL,sep="\t",quote="");head(data2);'\
'data3=read.table(args[3], header=F,row.names=NULL,sep="\t");colnames(data3)[2]="coral_blast_hit";head(data3);'\
'ann_isoExp=merge(data1,data2[,c(1,2,17,11)],by.x="geneName",by.y="qseqid",all.x = TRUE);'\
'ann_isoExp2=merge(ann_isoExp,data3[,c(1,2)],by.x="geneName",by.y="V1",all.x = TRUE);'\
'write.table(ann_isoExp2,"ann_isoExp", sep="\t", quote=F, row.names=F, col.names=T);' exp_isoformTPM $swiss_ann trinityVsSymb_best.sig.cladeA
tail -n+2 exp_isoformTPM | awk '{print $1"|"}' | grep -F -f - $LongOrfs | sed 's/>//' > LongOrfs.cladeA.key
filter_fasta.py --input_fasta_fp $LongOrfs --output_fasta_fp LongOrfs.cladeA.pep --seq_id_fp LongOrfs.cladeA.key
filter_fasta.py --input_fasta_fp $LongOrfs.complete --output_fasta_fp LongOrfs.cladeA.pep.complete --seq_id_fp LongOrfs.cladeA.key
grep "^>" LongOrfs.cladeA.pep | wc -l ## 147700 ## no of all possible ORF
grep "^>" LongOrfs.cladeA.pep | awk -F '[>|]' '{print $2}' | sort | uniq | wc -l ## 48939 ## no of transcripts with ORF
grep "^>" LongOrfs.cladeA.pep.complete | wc -l ## 81218 ## no of all possible complete ORF
grep "^>" LongOrfs.cladeA.pep.complete | awk -F '[>|]' '{print $2}' | sort | uniq | wc -l ## 33553 ## no of transcripts with complete ORF
cat exp_isoformTPM | awk '{print $1}' | grep -w -F -f - $swiss_ann > $swiss_ann.cladeA
cat $swiss_ann.cladeA | awk '{print $1}' | sort | uniq | wc -l
##############
cd $abundFilter
mkdir S_spCCMP2430
grep "S_spCCMP2430.est" ${p_asteroides}/blast_out/trinityVsSymb_best.sig > S_spCCMP2430/trinityVsSymb_best.sig.S_spCCMP2430
cat S_spCCMP2430/trinityVsSymb_best.sig.S_spCCMP2430 | awk '{print $1}' > S_spCCMP2430_transIDs
head -n1 transcripts.lengthes > S_spCCMP2430/transcripts.lengthes
grep -w -F -f S_spCCMP2430_transIDs transcripts.lengthes >> S_spCCMP2430/transcripts.lengthes
grep -w -F -f S_spCCMP2430_transIDs $gene_transcript_map > S_spCCMP2430/gene_transcript_map
for f in *.quant.counts;do
  head -n1 $f > S_spCCMP2430/$f
  grep -w -F -f S_spCCMP2430_transIDs $f >> S_spCCMP2430/$f
done
cd S_spCCMP2430

#module load R/3.0.1
while read identifier;do
  echo $(pwd) $identifier
  bash $script_path/run_calcTPM.sh "$(pwd)" "$identifier" "transcripts.lengthes" "gene_transcript_map" ${script_path}/calcTPM2.R
  #Rscript ${script_path}/calcTPM2.R "$(pwd)" "$identifier" "transcripts.lengthes" "gene_transcript_map" >> targets_list
done < $identifiers
bash $script_path/abund_est.sh

## define the species specific transcriptome and related annotations
filter_fasta.py --input_fasta_fp $FCS_ann_exp_tran --output_fasta_fp S_spCCMP2430_2016.fasta --seq_id_fp exp_isoformTPM
S_spCCMP2430_2016_trans=$abundFilter/S_spCCMP2430/S_spCCMP2430_2016.fasta
Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1], header=T,row.names=NULL);head(data1);'\
'data2=read.table(args[2], header=T,row.names=NULL,sep="\t",quote="");head(data2);'\
'data3=read.table(args[3], header=F,row.names=NULL,sep="\t");colnames(data3)[2]="coral_blast_hit";head(data3);'\
'ann_isoExp=merge(data1,data2[,c(1,2,17,11)],by.x="geneName",by.y="qseqid",all.x = TRUE);'\
'ann_isoExp2=merge(ann_isoExp,data3[,c(1,2)],by.x="geneName",by.y="V1",all.x = TRUE);'\
'write.table(ann_isoExp2,"ann_isoExp", sep="\t", quote=F, row.names=F, col.names=T);' exp_isoformTPM $swiss_ann trinityVsSymb_best.sig.S_spCCMP2430
tail -n+2 exp_isoformTPM | awk '{print $1"|"}' | grep -F -f - $LongOrfs | sed 's/>//' > LongOrfs.S_spCCMP2430.key
filter_fasta.py --input_fasta_fp $LongOrfs --output_fasta_fp LongOrfs.S_spCCMP2430.pep --seq_id_fp LongOrfs.S_spCCMP2430.key
filter_fasta.py --input_fasta_fp $LongOrfs.complete --output_fasta_fp LongOrfs.S_spCCMP2430.pep.complete --seq_id_fp LongOrfs.S_spCCMP2430.key
grep "^>" LongOrfs.S_spCCMP2430.pep | wc -l ## 21039 ## no of all possible ORF
grep "^>" LongOrfs.S_spCCMP2430.pep | awk -F '[>|]' '{print $2}' | sort | uniq | wc -l ## 8840 ## no of transcripts with ORF
grep "^>" LongOrfs.S_spCCMP2430.pep.complete | wc -l ## 7261 ## no of all possible complete ORF
grep "^>" LongOrfs.S_spCCMP2430.pep.complete | awk -F '[>|]' '{print $2}' | sort | uniq | wc -l ## 3884 ## no of transcripts with complete ORF
cat exp_isoformTPM | awk '{print $1}' | grep -w -F -f - $swiss_ann > $swiss_ann.S_spCCMP2430
cat $swiss_ann.S_spCCMP2430 | awk '{print $1}' | sort | uniq | wc -l
##############
## Assessement of the transcriptome
cd $compAnalysis
module load Bioperl/1.6.923
perl ${script_path}/seq_stats.pl $p_ast2016_trans > $p_ast2016_trans.MatzStat
perl ${script_path}/seq_stats.pl $spC15_2016_trans > $spC15_2016_trans.MatzStat
perl ${script_path}/seq_stats.pl $cladeA_2016_trans > $cladeA_2016_trans.MatzStat
perl ${script_path}/seq_stats.pl $S_spCCMP2430_2016_trans > $S_spCCMP2430_2016_trans.MatzStat
perl ${script_path}/seq_stats.pl $zoox_trans > $zoox_trans.MatzStat

module load trinity/2.2.0
TrinityStats.pl $p_ast2016_trans > $p_ast2016_trans.TrinityStat
TrinityStats.pl $spC15_2016_trans > $spC15_2016_trans.TrinityStat
TrinityStats.pl $cladeA_2016_trans > $cladeA_2016_trans.TrinityStat
TrinityStats.pl $S_spCCMP2430_2016_trans > $S_spCCMP2430_2016_trans.TrinityStat
TrinityStats.pl $zoox_trans > $zoox_trans.TrinityStat
#####################
## blast against the published pastreoids transcriptome
# a Matz lab paper (http://onlinelibrary.wiley.com/doi/10.1111/mec.12390/abstract)
# The annotated transcriptome (http://www.bio.utexas.edu/research/matz_lab/matzlab/Data.html).
# Porites astreoides (adult, Symbiodinium-specific reads excluded)
cd ${p_asteroides}/resources/coral_trans/P_astreoides
module load Bioperl/1.6.923
perl ${script_path}/seq_stats.pl past.fasta > past.fasta.MatzStat
# change to Trinity format and calc Trinity stats
awk '{print $1,$2}' past.fasta > TRN1_past.fasta
sed 's/>.*=/>/' TRN1_past.fasta > TRN2_past.fasta
grep "^>" TRN2_past.fasta | sort | uniq -c > TRN_count
gene=1
while read x id;do
 isoform=1
 while [  $isoform -le $x ]; do
  echo $isoform
  sed -i "0,/$id/ s//>TR${gene}_c0_g1_i${isoform}/" TRN2_past.fasta;
  let isoform=isoform+1
 done
 let gene=gene+1
done < TRN_count
module load trinity/2.2.0
TrinityStats.pl TRN2_past.fasta > TRN2_past.fasta.TrinityStat
# identify Long ORFs
module load TransDecoder/2.0.1
TransDecoder.LongOrfs -t past.fasta
matzLongOrfs=$(pwd)/past.fasta.transdecoder_dir/longest_orfs.pep
grep -A1 "type:complete" $matzLongOrfs | grep -v "^--" > $matzLongOrfs.complete 
grep "^>" $matzLongOrfs | wc -l ## 18351 ## no of all possible ORF
grep "^>" $matzLongOrfs | awk -F '[>|]' '{print $2}' | sort | uniq | wc -l ## 15183 ## no of transcripts with ORF
grep "^>" $matzLongOrfs.complete | wc -l ## 3666 ## no of all possible complete ORF
grep "^>" $matzLongOrfs.complete | awk -F '[>|]' '{print $2}' | sort | uniq | wc -l ## 2932 ## no of transcripts with complete ORF

## Run blast
module load BLAST+/2.2.30
## header of blast output= qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
mkdir ${p_asteroides}/uniprot
cd ${p_asteroides}/uniprot
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
makeblastdb -in uniprot_sprot.fasta -dbtype prot

mkdir ${p_asteroides}/resources/coral_trans/P_astreoides/blastx_dir
cd  ${p_asteroides}/resources/coral_trans/P_astreoides/blastx_dir
cp ../past.fasta .
perl ${script_path}/splitFasta.pl past.fasta 30  ## 500 for uniprot_uniref90

# header of blast output= qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
for f in subset*_past.fasta; do
  qsub -v input=$f,DB=${p_asteroides}/uniprot/uniprot_sprot.fasta,label="uniprot_sprot" ${script_path}/blastx_targetDB.sh;
done
cat subset*_past.fasta.bx > ../uniprot_sprot.blastx.outfmt6             
cd ../
cat uniprot_sprot.blastx.outfmt6 | awk '$11 <= 1e-3' > uniprot_sprot.blastx.outfmt6.sig    ## 15492
sort -k1,1 -k12,12nr -k11,11n  uniprot_sprot.blastx.outfmt6.sig | sort -u -k1,1 --merge > uniprot_sprot.blastx.outfmt6.sig.best             ## 206780
header='qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tqcovs\tqcovhsp\tstitle'
sed -i -e 1i"${header}" uniprot_sprot.blastx.outfmt6.sig.best

#mkdir ${p_asteroides}/resources/coral_trans/P_astreoides/P_ast.BlastDB
#cd ${p_asteroides}/resources/coral_trans/P_astreoides/P_ast.BlastDB
#cp ../past.fasta .
#module load BLAST+/2.2.30
#makeblastdb -in past.fasta -input_type fasta -dbtype nucl
#blastn -query $p_ast2016_trans \
#       -db past.fasta \
#       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp" \
#       -dust 'yes' -best_hit_overhang 0.25 -best_hit_score_edge 0.25 -evalue 1e-5 \
#       -max_target_seqs 20  -out p_ast2016VsMatz2013 ## -perc_identity 90 -qcov_hsp_perc 50
#sort -k1,1 -k12,12nr -k11,11n p_ast2016VsMatz2013 | sort -u -k1,1 --merge > p_ast2016VsMatz2013_best
#wc -l p_ast2016VsMatz2013_best ## 82384  (out of 129718)
#cat p_ast2016VsMatz2013_best | awk '$11 <= 1e-5' | wc -l        ##  82384
#cat p_ast2016VsMatz2013 | awk -F'\t' '{print $2}' | sort | uniq | wc -l  ## 22212
#cat p_ast2016VsMatz2013_best | awk -F'\t' '{print $2}' | sort | uniq | wc -l ## 18944
 
#mkdir ${p_asteroides}/resources/coral_trans/P_astreoides/newAsm.BlastDB
#cd ${p_asteroides}/resources/coral_trans/P_astreoides/newAsm.BlastDB
#cp $p_ast2016_trans p_ast2016.fasta
#module load BLAST+/2.2.30
#makeblastdb -in p_ast2016.fasta -input_type fasta -dbtype nucl
#blastn -query ../past.fasta \
#       -db p_ast2016.fasta \
#       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp" \
#       -dust 'yes' -best_hit_overhang 0.25 -best_hit_score_edge 0.25 -evalue 1e-5 \
#       -max_target_seqs 20  -out Matz2013Vsp_ast2016 ## -perc_identity 90 -qcov_hsp_perc 50
#sort -k1,1 -k12,12nr -k11,11n Matz2013Vsp_ast2016 | sort -u -k1,1 --merge > Matz2013Vsp_ast2016_best
#wc -l Matz2013Vsp_ast2016_best ## 22204  (out of 30740)
#cat Matz2013Vsp_ast2016_best | awk '$11 <= 1e-5' | wc -l        ##  22204
#cat Matz2013Vsp_ast2016 | awk -F'\t' '{print $2}' | sort | uniq | wc -l  ## 69154
#cat Matz2013Vsp_ast2016_best | awk -F'\t' '{print $2}' | sort | uniq | wc -l ## 17477

mkdir ${p_asteroides}/compAsm
cd ${p_asteroides}/compAsm
cp $p_asteroides/compAnalysis/Trinity.clean.201.exp.FCS.fasta Trinity.fasta ## Trinity assembly (the submiited version to ncbi)
cp ${p_asteroides}/blast_out/trinityVsSymb_best.sig.fasta zoox.fasta ## symbiotic transcriptome
cp $abundFilter/spC15/spC15_2016.fasta .
cp $abundFilter/cladeA/cladeA_2016.fasta .
cp $abundFilter/S_spCCMP2430/S_spCCMP2430_2016.fasta .
cp $p_ast2016_trans p_ast2016.fasta  ## our new clean coral assembly
cp ${p_asteroides}/resources/coral_trans/P_astreoides/past.fasta matz2013.fasta
#python ${script_path}/blast_rbh.py -a nucl -t blastn -o RBH p_ast2016.fasta matz2013.fasta
#cat RBH | awk '$3 > $4' > oldAssVsnewAss_best_better
#wc -l oldAssVsnewAss_best_better ## 1382
#cat RBH | awk '$3 == $4' > oldAssVsnewAss_best_equal
#wc -l oldAssVsnewAss_best_equal ## 6
#cat RBH | awk '$3 < $4' > oldAssVsnewAss_best_less
#wc -l oldAssVsnewAss_best_less ## 1085
#cat oldAssVsnewAss_best_better | awk '{ sum+=$4} END {print sum}' ## 824062
#cat oldAssVsnewAss_best_better | awk '{ sum+=$3} END {print sum}' ## 1164786 (i.e. difference of 340724 =~0.34Mb)
#cat oldAssVsnewAss_best_less | awk '{ sum+=$4} END {print sum}' ## 639843
#cat oldAssVsnewAss_best_less | awk '{ sum+=$3} END {print sum}' ## 484854 (i.e. difference of 154989 =~0.15Mb)

module load transrate/1.0.1
transrate --assembly Trinity.fasta,zoox.fasta,spC15_2016.fasta,cladeA_2016.fasta,S_spCCMP2430_2016.fasta,p_ast2016.fasta,matz2013.fasta --output transrate_results &> trans.log
transrate --assembly matz2013.fasta --reference p_ast2016.fasta --output matz2013VSp_ast2016_transrate &> matz2013VSp_ast2016_trans.log
awk -F',' 'BEGIN{OFS=",";} {if($10=="true")print $1,$2,$11,$12;}' matz2013VSp_ast2016_transrate/matz2013/contigs.csv | sort -t "," -k 4b,4 > CRBB_hits.csv
tail -n+2 transrate_results/p_ast2016/contigs.csv | awk -F',' 'BEGIN{OFS=",";} {print $1,$2;}' | sort -t "," -k 1b,1 > p_ast2016_len.csv
echo "qname qlen rcoverage rname rlen" > CRBB_withReflen
join -1 4 -2 1 -t "," CRBB_hits.csv p_ast2016_len.csv | awk -F',' '{print $2,$3,$4,$1,$5;}' >> CRBB_withReflen
tail -n+2 CRBB_withReflen | awk '$5 > $2' > CRBB_oldAssVsnewAss_best_better
better=$(cat CRBB_oldAssVsnewAss_best_better | wc -l) ## 15115
tail -n+2 CRBB_withReflen | awk '$5 == $2' > CRBB_oldAssVsnewAss_best_equal
equal=$(cat CRBB_oldAssVsnewAss_best_equal | wc -l) ## 14
tail -n+2 CRBB_withReflen | awk '$5 < $2' > CRBB_oldAssVsnewAss_best_less
less=$(cat CRBB_oldAssVsnewAss_best_less | wc -l) ## 6104
newBetter=$(cat CRBB_oldAssVsnewAss_best_better | awk '{ sum+=$5} END {print sum}') ## 36863704
oldBetter=$(cat CRBB_oldAssVsnewAss_best_better | awk '{ sum+=$2} END {print sum}') ## 7022034 (i.e. difference of 29841670 =~29.5Mb)
newLess=$(cat CRBB_oldAssVsnewAss_best_less | awk '{ sum+=$5} END {print sum}') ## 2539349
oldLess=$(cat CRBB_oldAssVsnewAss_best_less | awk '{ sum+=$2} END {print sum}') ## 5032274 (i.e. difference of 2492925 =~2.5Mb)
cat trans.log matz2013VSp_ast2016_trans.log > compAsm.log ## edit manually to remove reduandancy
echo "no of transcripts that gain length in the new assembly = $better transcript" >> compAsm.log
echo "total length of longer transcripts in the new assembly = $newBetter bp" >> compAsm.log
echo "total gain of length = $(($newBetter-$oldBetter)) bp" >> compAsm.log
echo "=============================================" >> compAsm.log
echo "no of transcripts that did not change = $equal transcript" >> compAsm.log
echo "=============================================" >> compAsm.log
echo "no of transcripts that lost length in the new assembly = $less transcript" >> compAsm.log
echo "total length of shorter transcripts in the new assembly = $newLess bp" >> compAsm.log
echo "total loss of length = $(($oldLess-$newLess)) bp" >> compAsm.log
#####################
scp $p_asteroides/compAnalysis/Trinity.clean.201.exp.FCS.fasta tmansour@loretta.hpcf.upr.edu:/storage/prcen/coral/p_ast.assemblies.2016/. ## final all transcriptome with initial annotation
scp ${p_asteroides}/blast_out/trinityVsSymb_best.sig.fasta tmansour@loretta.hpcf.upr.edu:/storage/prcen/coral/p_ast.assemblies.2016/. ## the total zooxanrthellae transcriptome with initial annotation
scp $abundFilter/coral/{p_ast2016.fasta,ann_isoExp,LongOrfs.coral.key,*.*Stat} tmansour@loretta.hpcf.upr.edu:/storage/prcen/coral/p_ast.assemblies.2016/coral/. ## clean porites asteroids transcriptome with annotation files
scp $abundFilter/spC15/{spC15_2016.fasta,ann_isoExp,LongOrfs.spC15.key,*.*Stat} tmansour@loretta.hpcf.upr.edu:/storage/prcen/coral/p_ast.assemblies.2016/spC15/. ## clean porites asteroids transcriptome with annotation files
scp $abundFilter/cladeA/{cladeA_2016.fasta,ann_isoExp,LongOrfs.cladeA.key,*.*Stat} tmansour@loretta.hpcf.upr.edu:/storage/prcen/coral/p_ast.assemblies.2016/cladeA/. ## clean porites asteroids transcriptome with annotation files
scp $abundFilter/S_spCCMP2430/{S_spCCMP2430_2016.fasta,ann_isoExp,LongOrfs.S_spCCMP2430.key,*.*Stat} tmansour@loretta.hpcf.upr.edu:/storage/prcen/coral/p_ast.assemblies.2016/S_spCCMP2430/. ## clean porites asteroids transcriptome with annotation files
#####################
#####################
