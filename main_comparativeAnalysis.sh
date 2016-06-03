#!/bin/sh

## construction of the basic diretory structure
git clone https://github.com/drtamermansour/p_asteroides_CompAna.git
cd p_asteroides_CompAna
p_asteroides=$(pwd)

## create working directory and define paths for raw data and scripts
mkdir -p $p_asteroides/{c_abundFilter,compAnalysis,resources,blast_out}
script_path=$p_asteroides/scripts
abundFilter=$p_asteroides/c_abundFilter   ## copy the filtered reads from the p_asteroides assembly project
compAnalysis=$p_asteroides/compAnalysis   ## copy the annotated transcripome from the p_asteroides assembly project
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

filter_fasta.py --input_fasta_fp $FCS_ann_exp_tran --output_fasta_fp noSymb.fasta --seq_id_fp trinityVsSymb_best.sig --negate
grep "^>" noSymb.fasta | wc -l  ## 681078
noSymb_transcriptome=${p_asteroides}/blast_out/noSymb.fasta

## select annotation for significant hits
#cat trinityVsSymb_best.sig | awk '{print ">"$1"|"}' | grep -F -f - $LongOrfs | sed 's/>//' > LongOrfs.trinityVsSymb.key
#filter_fasta.py --input_fasta_fp $LongOrfs --output_fasta_fp $LongOrfs.trinityVsSymb --seq_id_fp LongOrfs.trinityVsSymb.key
#filter_fasta.py --input_fasta_fp $LongOrfs.complete --output_fasta_fp $LongOrfs.trinityVsSymb.complete --seq_id_fp LongOrfs.trinityVsSymb.key
#cat trinityVsSymb_best.sig | awk '{print $1}' | grep -w -F -f - $swiss_ann > $swiss_ann.trinityVsSymb
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
## Abundance estimation
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
###################
## Annotation files (has header)
# $compAnalysis/uniprot_sprot.blastx.outfmt6.sig.best.exp

## TPM (adult or larva specific) (has header)
# $p_asteroides/c_abundFilter/larva_isoformTPM
# $p_asteroides/c_abundFilter/adultOnly_isoformTPM
# $p_asteroides/c_abundFilter/larvaOnly_isoformTPM
# $p_asteroides/c_abundFilter/unexp_isoformTPM
# $p_asteroides/c_abundFilter/exp_isoformTPM

## hits in symb or coral
# ${p_asteroides}/blast_out/trinityVsSymb_best.sig
# ${p_asteroides}/blast_out/trinityVscoral_best.sig
# ${p_asteroides}/blast_out/recognized_best.sig

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
grep "^>" LongOrfs.coral.pep.complete | awk -F '[>|]' '{print $2}' | sort | uniq | wc -l ## 17277
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
grep "^>" LongOrfs.spC15.pep.complete | awk -F '[>|]' '{print $2}' | sort | uniq | wc -l ## 20741
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
grep "^>" LongOrfs.cladeA.pep.complete | awk -F '[>|]' '{print $2}' | sort | uniq | wc -l ## 33553
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
grep "^>" LongOrfs.S_spCCMP2430.pep.complete | awk -F '[>|]' '{print $2}' | sort | uniq | wc -l ## 3884
##############
## Assessement of the transcriptome
cd $compAnalysis
module load Bioperl/1.6.923
perl ${script_path}/seq_stats.pl $p_ast2016_trans > $p_ast2016_trans.MatzStat
perl ${script_path}/seq_stats.pl $spC15_2016_trans > $spC15_2016_trans.MatzStat
perl ${script_path}/seq_stats.pl $cladeA_2016_trans > $cladeA_2016_trans.MatzStat
perl ${script_path}/seq_stats.pl $S_spCCMP2430_2016_trans > $S_spCCMP2430_2016_trans.MatzStat

module load trinity/2.2.0
TrinityStats.pl $p_ast2016_trans > $p_ast2016_trans.TrinityStat
TrinityStats.pl $spC15_2016_trans > $spC15_2016_trans.TrinityStat
TrinityStats.pl $cladeA_2016_trans > $cladeA_2016_trans.TrinityStat
TrinityStats.pl $S_spCCMP2430_2016_trans > $S_spCCMP2430_2016_trans.TrinityStat
#####################
## blast against the published pastreoids transcriptome
# a Matz lab paper (http://onlinelibrary.wiley.com/doi/10.1111/mec.12390/abstract)
# The annotated transcriptome (http://www.bio.utexas.edu/research/matz_lab/matzlab/Data.html).
# Porites astreoides (adult, Symbiodinium-specific reads excluded)
cd ${p_asteroides}/resources/coral_trans/P_astreoides
module load Bioperl/1.6.923
perl ${script_path}/seq_stats.pl past.fasta > past.fasta.MatzStat

mkdir ${p_asteroides}/resources/coral_trans/P_astreoides/P_ast.BlastDB
cd ${p_asteroides}/resources/coral_trans/P_astreoides/P_ast.BlastDB
cp ../past.fasta .
module load BLAST+/2.2.30
makeblastdb -in past.fasta -input_type fasta -dbtype nucl
blastn -query $p_ast2016_trans \
       -db past.fasta \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp" \
       -dust 'yes' -best_hit_overhang 0.25 -best_hit_score_edge 0.25 -evalue 1e-5 \
       -max_target_seqs 20  -out p_ast2016VsMatz2013 ## -perc_identity 90 -qcov_hsp_perc 50
sort -k1,1 -k12,12nr -k11,11n p_ast2016VsMatz2013 | sort -u -k1,1 --merge > p_ast2016VsMatz2013_best
wc -l p_ast2016VsMatz2013_best ## 82384  (out of 129718)
cat p_ast2016VsMatz2013_best | awk '$11 <= 1e-5' | wc -l        ##  82384
cat p_ast2016VsMatz2013 | awk -F'\t' '{print $2}' | sort | uniq | wc -l  ## 22212
cat p_ast2016VsMatz2013_best | awk -F'\t' '{print $2}' | sort | uniq | wc -l ## 18944
 
mkdir ${p_asteroides}/resources/coral_trans/P_astreoides/newAsm.BlastDB
cd ${p_asteroides}/resources/coral_trans/P_astreoides/newAsm.BlastDB
cp $p_ast2016_trans p_ast2016.fasta
module load BLAST+/2.2.30
makeblastdb -in p_ast2016.fasta -input_type fasta -dbtype nucl
blastn -query ../past.fasta \
       -db p_ast2016.fasta \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp" \
       -dust 'yes' -best_hit_overhang 0.25 -best_hit_score_edge 0.25 -evalue 1e-5 \
       -max_target_seqs 20  -out Matz2013Vsp_ast2016 ## -perc_identity 90 -qcov_hsp_perc 50
sort -k1,1 -k12,12nr -k11,11n Matz2013Vsp_ast2016 | sort -u -k1,1 --merge > Matz2013Vsp_ast2016_best
wc -l Matz2013Vsp_ast2016_best ## 22204  (out of 30740)
cat Matz2013Vsp_ast2016_best | awk '$11 <= 1e-5' | wc -l        ##  22204
cat Matz2013Vsp_ast2016 | awk -F'\t' '{print $2}' | sort | uniq | wc -l  ## 69154
cat Matz2013Vsp_ast2016_best | awk -F'\t' '{print $2}' | sort | uniq | wc -l ## 17477

mkdir ${p_asteroides}/compAsm
cd ${p_asteroides}/compAsm
cp $p_ast2016_trans p_ast2016.fasta
cp ${p_asteroides}/resources/coral_trans/P_astreoides/past.fasta matz2013.fasta
python ${script_path}/blast_rbh.py -a nucl -t blastn -o RBH p_ast2016.fasta matz2013.fasta
cat RBH | awk '$3 > $4' > oldAssVsnewAss_best_better
wc -l oldAssVsnewAss_best_better ## 1382
cat RBH | awk '$3 == $4' > oldAssVsnewAss_best_equal
wc -l oldAssVsnewAss_best_equal ## 6
cat RBH | awk '$3 < $4' > oldAssVsnewAss_best_less
wc -l oldAssVsnewAss_best_less ## 1085
cat oldAssVsnewAss_best_better | awk '{ sum+=$4} END {print sum}' ## 824062
cat oldAssVsnewAss_best_better | awk '{ sum+=$3} END {print sum}' ## 1164786 (i.e. difference of 340724 =~0.34Mb)
cat oldAssVsnewAss_best_less | awk '{ sum+=$4} END {print sum}' ## 639843
cat oldAssVsnewAss_best_less | awk '{ sum+=$3} END {print sum}' ## 484854 (i.e. difference of 154989 =~0.15Mb)

module load transrate/1.0.1
transrate --assembly p_ast2016.fasta,matz2013.fasta --output transrate_results &> trans.log
transrate --assembly matz2013.fasta --reference p_ast2016.fasta --output matz2013VSp_ast2016_transrate &> matz2013VSp_ast2016_trans.log
awk -F',' 'BEGIN{OFS=",";} {if($10=="true")print $1,$2,$11,$12;}' matz2013VSp_ast2016_transrate/matz2013/contigs.csv | sort -t "," -k 4b,4 > CRBB_hits.csv
tail -n+2 transrate_results/p_ast2016/contigs.csv | awk -F',' 'BEGIN{OFS=",";} {print $1,$2;}' | sort -t "," -k 1b,1 > p_ast2016_len.csv
echo "qname qlen rcoverage rname rlen" > CRBB_withReflen
join -1 4 -2 1 -t "," CRBB_hits.csv p_ast2016_len.csv | awk -F',' '{print $2,$3,$4,$1,$5;}' >> CRBB_withReflen
cat CRBB_withReflen | awk '$5 > $2' > CRBB_oldAssVsnewAss_best_better
better=$(cat CRBB_oldAssVsnewAss_best_better | wc -l) ## 15115
cat CRBB_withReflen | awk '$5 == $2' > CRBB_oldAssVsnewAss_best_equal
equal=$(cat CRBB_oldAssVsnewAss_best_equal | wc -l) ## 14
cat CRBB_withReflen | awk '$5 < $2' > CRBB_oldAssVsnewAss_best_less
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
## blast unrecog transcripts aganist NCBI nuclutide database
module load BLAST+/2.2.30
cd $p_asteroides/resources
mkdir nt_DB && cd nt_DB
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
gunzip *.gz
tar -xvf *.tar
export BLASTDB= $p_asteroides/resources/nt_DB:$BLASTDB
## alternatively we can use the pre-made HPC databse
export BLASTDB=/mnt/research/common-data/Bio/blastdb:$BLASTDB

# splitting the transcriptome into 40 chunks
mkdir ${p_asteroides}/blast_out/unrec_chunks && cd ${p_asteroides}/blast_out/unrec_chunks
cp ../unrecognizied_est.fasta .
perl ${script_path}/splitFasta.pl unrecognizied_est.fasta 40

# blasting all chunks to nt database
for f in subset*.fasta; do
    qsub -v input=$f ${script_path}/blastn.sh; done
cat subset*.fasta.br > unrecognizied.br
rm subset*.fasta subset*.fasta.br
wc -l unrecognizied.br ## 204266

sort -k1,1 -k11,11g  unrecognizied.br | sort -u -k1,1 --merge >  unrecognizied.br_best
wc -l unrecognizied.br_best ## 35830

## filter the hits to exclude transcripts recognizied by blast to Matz lab assembly (the blast output  ${p_asteroides}/resources/P_ast.transcriptome/P_ast.BlastDB/unrecogniziedVspublishedAss_best)
## I got important unix syntex help from: http://unix.stackexchange.com/questions/110645/select-lines-from-text-file-which-have-ids-listed-in-another-file
cd  ${p_asteroides}/blast_out
grep "^>" recog.uncontaminated2.fasta | awk '{print $1}' | cut -d ">" -f2 > recog.uncontaminated2.headers
grep "^>" unrecognizied_est2.fasta | awk '{print $1}' | cut -d ">" -f2 > unrecognizied_est2.headers

cd ${p_asteroides}/blast_out/unrec_chunks
#grep -Fwf unrecognizied.br_best ${p_asteroides}/blast_out/recog.uncontaminated2.headers > unrecognizied.br_best.rec
#grep -Fwf unrecognizied.br_best ${p_asteroides}/blast_out/unrecognizied_est2.headers > unrecognizied.br_best.unrec
while read pat; do grep -w "^$pat" unrecognizied.br_best; done < ${p_asteroides}/blast_out/recog.uncontaminated2.headers > unrecognizied.br_best.rec
wc -l unrecognizied.br_best.rec  ## 7792
cat unrecognizied.br_best.rec | awk '$11 <= 0.001' >  unrecognizied.br_best.rec_sig
wc -l unrecognizied.br_best.rec_sig ## 6059

while read pat; do grep -w "^$pat" unrecognizied.br_best; done < ${p_asteroides}/blast_out/unrecognizied_est2.headers > unrecognizied.br_best.unrec
wc -l unrecognizied.br_best.unrec  ## 28038
cat unrecognizied.br_best.unrec | awk '$11 <= 0.001' >  unrecognizied.br_best.unrec_sig
wc -l unrecognizied.br_best.unrec_sig ## 23593

cat unrecognizied.br_best.rec_sig | awk -F $'\t' '{A[$20]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > taxids.rec.count
cat unrecognizied.br_best.rec_sig | awk -F $'\t' '{A[$22]++}END{for(i in A)print i,"\t",A[i]}' | sort -k2,2nr -t "Ctr+v then tab" > species.rec.count
cat unrecognizied.br_best.rec_sig | awk -F $'\t' '{A[$21]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > kingdom.rec.count

cat unrecognizied.br_best.unrec_sig | awk -F $'\t' '{A[$20]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > taxids.unrec.count
cat unrecognizied.br_best.unrec_sig | awk -F $'\t' '{A[$22]++}END{for(i in A)print i,"\t",A[i]}' | sort -k2,2nr -t "Ctr+v then tab" > species.unrec.count
cat unrecognizied.br_best.unrec_sig | awk -F $'\t' '{A[$21]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > kingdom.unrec.count

#cat unrecognizied.br | awk -F $'\t' '$20 == 7227' > unrecognizied.br.D_melanogaster
#sort -k1,1 -k11,11g  unrecognizied.br.D_melanogaster | sort -u -k1,1 --merge >  unrecognizied.br.D_melanogaster_best
#wc -l unrecognizied.br.D_melanogaster_best ## 2596
grep -w "Porites lobata" unrecognizied.br_best.rec > P_lobata.rec
grep -w "Porites lobata" unrecognizied.br_best.unrec > P_lobata.unrec 

###
# splitting the transcriptome into 20 chunks
mkdir ${p_asteroides}/blast_out/rec1_chunks && cd ${p_asteroides}/blast_out/rec1_chunks
cp ../recog.uncontaminated.fasta .
perl ${script_path}/splitFasta.pl recog.uncontaminated.fasta 20

# blasting all chunks to nt database
module load BLAST+/2.2.30
export BLASTDB=/mnt/research/common-data/Bio/blastdb:$BLASTDB
for f in subset11_*.fasta; do blastn -query $f -db nt -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp stitle salltitles staxids sskingdoms sscinames" -dust 'yes' -num_threads 4 -max_target_seqs 10 -out $f.br; done
cat subset*.fasta.br > recog.br
rm subset*.fasta subset*.fasta.br
sort -k1,1 -k11,11g  recog.br | sort -u -k1,1 --merge >  recog.br_best
wc -l recog.br_best ## 10228
cat recog.br_best | awk '$11 <= 0.001' >  recog.br_best_sig
wc -l recog.br_best_sig ## 9127

cat recog.br_best_sig | awk -F $'\t' '{A[$20]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > taxids.rec.count
cat recog.br_best_sig | awk -F $'\t' '{A[$22]++}END{for(i in A)print i,"\t",A[i]}' | sort -k2,2nr -t "Ctr+v then tab" > species.rec.count
cat recog.br_best_sig | awk -F $'\t' '{A[$21]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > kingdom.rec.count

###
mkdir ${p_asteroides}/blast_out/symp && cd ${p_asteroides}/blast_out/symp
cp ../cladeAandB.fasta .
perl ${script_path}/splitFasta.pl cladeAandB.fasta 20

# blasting all chunks to nt database
module load BLAST+/2.2.30
export BLASTDB=/mnt/research/common-data/Bio/blastdb:$BLASTDB
for f in subset4_*.fasta; do blastn -query $f -db nt -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp stitle salltitles staxids sskingdoms sscinames" -dust 'yes' -num_threads 4 -max_target_seqs 10 -out $f.br; done
cat subset*.fasta.br > symp.br
rm subset*.fasta subset*.fasta.br
sort -k1,1 -k11,11g  symp.br | sort -u -k1,1 --merge >  symp.br_best
wc -l symp.br_best ## 5647
cat symp.br_best | awk '$11 <= 0.001' >  symp.br_best_sig
wc -l symp.br_best_sig ## 4691

cat symp.br_best_sig | awk -F $'\t' '{A[$20]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > taxids.symp.count
cat symp.br_best_sig | awk -F $'\t' '{A[$22]++}END{for(i in A)print i,"\t",A[i]}' | sort -k2,2nr -t "Ctr+v then tab" > species.symp.count
cat symp.br_best_sig | awk -F $'\t' '{A[$21]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > kingdom.symp.count

###
mkdir ${p_asteroides}/blast_out/matz && cd ${p_asteroides}/blast_out/matz
cp ${p_asteroides}/resources/P_ast.transcriptome/past.fasta .
perl ${script_path}/splitFasta.pl past.fasta 20

# blasting all chunks to nt database
module load BLAST+/2.2.30
export BLASTDB=/mnt/research/common-data/Bio/blastdb:$BLASTDB
for f in subset20_*.fasta; do blastn -query $f -db nt -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp stitle salltitles staxids sskingdoms sscinames" -dust 'yes' -num_threads 4 -max_target_seqs 10 -out $f.br; done
cat subset*.fasta.br > matz.br
rm subset*.fasta subset*.fasta.br
sort -k1,1 -k11,11g  matz.br | sort -u -k1,1 --merge >  matz.br_best
wc -l matz.br_best ## 1820
cat matz.br_best | awk '$11 <= 0.001' >  matz.br_best_sig
wc -l matz.br_best_sig ## 1650

cat matz.br_best_sig | awk -F $'\t' '{A[$20]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > taxids.matz.count
cat matz.br_best_sig | awk -F $'\t' '{A[$22]++}END{for(i in A)print i,"\t",A[i]}' | sort -k2,2nr -t "Ctr+v then tab" > species.matz.count
cat matz.br_best_sig | awk -F $'\t' '{A[$21]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > kingdom.matz.count

###
R CMD BATCH ${script_path}/mergetaxids.R
####################################################################################
# splitting the transcriptome into 500 chunks
mkdir ${p_asteroides}/blast_out/unrec_chunks && cd ${p_asteroides}/blast_out/unrec_chunks
cp ../unrecognizied_est.fasta .
perl ${script_path}/splitFasta.pl unrecognizied_est.fasta 500

for f in subset140_*.fasta; do
    echo $f
    mkdir $f.temp && cd $f.temp
    cp ../$f .
    perl ${script_path}/splitFasta.pl $f 100
    for s in subset1*_$f; do
	echo $s
	blastx -query $s -db nr -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp stitle salltitles staxids sskingdoms sscinames" -num_threads 4 -max_target_seqs 4 -out $s.bx
	cat $s.bx >> $f.bx; done
    cd ../; done

# blasting all chunks to nr database
for f in subset2*.fasta; do
    qsub -v input=$f ${script_path}/blastx.sh; done
#cat subset*.fasta.bx > unrecognizied.bx
#rm subset*.fasta subset*.fasta.bx
#wc -l unrecognizied.bx ## 
for f in subset*.fasta.bx; do
    sort -k1,1 -k11,11g  $f | sort -u -k1,1 --merge > $f"_best"; done
cat subset*.fasta.bx_best > unrecognizied.bx_best
while read pat; do grep -w "^$pat" unrecognizied.bx_best; done < ${p_asteroides}/blast_out/recog.uncontaminated2.headers > unrecognizied.bx_best.rec
wc -l unrecognizied.bx_best.rec  ## 
while read pat; do grep -w "^$pat" unrecognizied.bx_best; done < ${p_asteroides}/blast_out/unrecognizied_est2.headers > unrecognizied.bx_best.unrec
wc -l unrecognizied.bx_best.unrec  ## 

wc -l unrecognizied.bx_best ##
cat unrecognizied.bx_best | awk -F $'\t' '{A[$20]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > taxids.bx.count
cat unrecognizied.bx_best | awk -F $'\t' '{A[$22]++}END{for(i in A)print i,"\t",A[i]}' | sort -k2,2nr -t "Ctr+v then tab" > species.bx.count
cat unrecognizied.bx_best | awk -F $'\t' '{A[$21]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > kingdom.bx.count

cat unrecognizied.bx_best | awk '$11 <= 1e-3' >  unrecognizied.bx_best_sig

wc -l unrecognizied.bx_best_sig ##
cat unrecognizied.bx_best_sig | awk -F $'\t' '{A[$20]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > taxids.bx.sig.count

################
cd $p_asteroides/resources
mkdir uniprot_DB && cd uniprot_DB
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip -c ../uniprot_sprot.fasta.gz > uniprot_sprot.fasta
module load BLAST+/2.2.30
makeblastdb -in uniprot_sprot.fasta -parse_seqids -dbtype prot
blastx -query $f -db ${p_asteroides}/resources/uniprot_DB/uniprot_sprot.fasta -evalue 0.0001 \
       -num_threads 4 -max_target_seqs 4 -out $f.br \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp stitle salltitles staxids sskingdoms sscinames"

################
grep "^>" transcriptome.fa.annot | wc -l
grep "h=" transcriptome.fa.annot | wc -l
grep "ortho:" transcriptome.fa.annot | wc -l
grep "transcript family ortho to:" transcriptome.fa.annot | wc -l
grep "transcript family homol to:" transcriptome.fa.annot | wc -l

###############
cat $f | awk '$2 < 100'


#############

source env/bin/activate
cd khmer
p_asteroides=$"/mnt/lustre_scratch_2012/Tamer/p_asteroides"
data_path=$"/mnt/lustre_scratch_2012/Tamer/p_asteroides/data"
script_path=$"/mnt/home/mansourt/p_asteroides/scripts"
trimmed_data=${p_asteroides}/b_adap_remove 
exp_transcriptome=${p_asteroides}/b_diginormC25k20_2/trinity_out_dir/seqclean/Trinity.fasta.clean.201
compAnalysis=$"/mnt/lustre_scratch_2012/Tamer/p_asteroides/b_diginormC25k20_2/trinity_out_dir/seqclean"
#############

##R
#setwd("/Users/drtamermansour/Desktop/FileZilla_client/p_asteroides/histograms")
#data = read.table("allsamples.keep.hist", header=FALSE)
#plot(data$V1,log(data$V2), xlim=range(0,10000))
#plot(data$V1,data$V2, log = "xy", type="p", xlab="Kmer multiplicity" , ylab="Frequency of appearance", cex=1, lwd=1)
#plot(data$V1,data$V2, log = "xy", type="p", xlab="Kmer multiplicity" , ylab="Frequency of appearance", cex=1, lwd=1, xlim=range(1,100), ylim=range(1e+05,1e+10))




#########################################################################################################
## blast unrecog transcripts aganist NCBI est  database
module load BLAST+/2.2.30
cd $p_asteroides/resources
mkdir est_DB && cd est_DB
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/est.tar.gz
#wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/est_human.*.tar.gz 
#wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/est_mouse.*.tar.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/est_others.*.tar.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
gunzip *.gz
for f in *.tar; do tar -xvf $f; done
export BLASTDB=$p_asteroides/resources/est_DB:$BLASTDB

# splitting the transcriptome into 20 chunks
cd ${p_asteroides}/blast_out/rec1_chunks

# blasting all chunks to est database
for f in subset*_*.fasta; do blastn -query $f -db est_others -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp stitle salltitles staxids sskingdoms sscinames" -dust 'yes' -num_threads 4 -max_target_seqs 10 -out $f.est; done
cat subset*.fasta.est > recog.est
sort -k1,1 -k11,11g  recog.est | sort -u -k1,1 --merge >  recog.est_best
wc -l recog.est_best ## 41203
cat recog.est_best | awk '$11 <= 0.001' >  recog.est_best_sig
wc -l recog.est_best_sig ## 40806

cat recog.est_best_sig | awk -F $'\t' '{A[$20]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > taxids.rec_est.count
cat recog.est_best_sig | awk -F $'\t' '{A[$22]++}END{for(i in A)print i,"\t",A[i]}' | sort -k2,2nr -t "Ctr+v then tab" > species.rec_est.count
cat recog.est_best_sig | awk -F $'\t' '{A[$21]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > kingdom.rec_est.count
##
cd ${p_asteroides}/blast_out/symp

# blasting all chunks to est database
for f in subset*3_*.fasta; do blastn -query $f -db est_others -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp stitle salltitles staxids sskingdoms sscinames" -dust 'yes' -num_threads 4 -max_target_seqs 10 -out $f.est; done
cat subset*.fasta.est > symp.est
sort -k1,1 -k11,11g  symp.est | sort -u -k1,1 --merge >  symp.est_best
wc -l symp.est_best ## 8223
cat symp.est_best | awk '$11 <= 0.001' >  symp.est_best_sig
wc -l symp.est_best_sig ## 7820

cat symp.est_best_sig | awk -F $'\t' '{A[$20]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > taxids.symp_est.count
cat symp.est_best_sig | awk -F $'\t' '{A[$22]++}END{for(i in A)print i,"\t",A[i]}' | sort -k2,2nr -t "Ctr+v then tab" > species.symp_est.count
cat symp.est_best_sig | awk -F $'\t' '{A[$21]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > kingdom.symp_est.count
##
mkdir ${p_asteroides}/blast_out/unrec_chunks_unrec && cd ${p_asteroides}/blast_out/unrec_chunks_unrec
cp ../unrecognizied_est2.fasta .
perl ${script_path}/splitFasta.pl unrecognizied_est2.fasta 40

# blasting all chunks to est database
for f in subset*1_*.fasta; do blastn -query $f -db est_others -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp stitle salltitles staxids sskingdoms sscinames" -dust 'yes' -num_threads 4 -max_target_seqs 10 -out $f.est; done

cat subset*.fasta.est > unrecognizied.est
sort -k1,1 -k11,11g  unrecognizied.est | sort -u -k1,1 --merge >  unrecognizied.est_best
wc -l unrecognizied.est_best ## 40357
cat unrecognizied.est_best | awk '$11 <= 0.001' >  unrecognizied.est_best_sig
wc -l unrecognizied.est_best_sig ## 38831

cat unrecognizied.est_best_sig | awk -F $'\t' '{A[$20]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > taxids.bx.count
cat unrecognizied.est_best_sig | awk -F $'\t' '{A[$22]++}END{for(i in A)print i,"\t",A[i]}' | sort -k2,2nr -t "Ctr+v then tab" > species.bx.count
cat unrecognizied.est_best_sig | awk -F $'\t' '{A[$21]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > kingdom.bx.count
##
mkdir ${p_asteroides}/blast_out/unrec_chunks_rec && cd ${p_asteroides}/blast_out/unrec_chunks_rec
cp ../recog.uncontaminated2.fasta .
perl ${script_path}/splitFasta.pl recog.uncontaminated2.fasta 30

# blasting all chunks to est database
for f in subset*1_*.fasta; do blastn -query $f -db est_others -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp stitle salltitles staxids sskingdoms sscinames" -dust 'yes' -num_threads 4 -max_target_seqs 10 -out $f.est; done

cat subset*.fasta.est > recog2.est
sort -k1,1 -k11,11g  recog2.est | sort -u -k1,1 --merge >  recog2.est_best
wc -l recog2.est_best ## 69024
cat recog2.est_best | awk '$11 <= 0.001' >  recog2.est_best_sig
wc -l recog2.est_best_sig ## 68282

cat recog2.est_best_sig | awk -F $'\t' '{A[$20]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > taxids.rec2_est.count
cat recog2.est_best_sig | awk -F $'\t' '{A[$22]++}END{for(i in A)print i,"\t",A[i]}' | sort -k2,2nr -t "Ctr+v then tab" > species.rec2_est.count
cat recog2.est_best_sig | awk -F $'\t' '{A[$21]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > kingdom.rec2_est.count
##
cd ${p_asteroides}/blast_out/matz

# blasting all chunks to est database
for f in subset*6_*.fasta; do blastn -query $f -db est_others -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp stitle salltitles staxids sskingdoms sscinames" -dust 'yes' -num_threads 4 -max_target_seqs 10 -out $f.est; done
cat subset*.fasta.est > matz.est
sort -k1,1 -k11,11g  matz.est | sort -u -k1,1 --merge >  matz.est_best
wc -l matz.est_best ## 7050
cat matz.est_best | awk '$11 <= 0.001' >  matz.est_best_sig
wc -l matz.est_best_sig ## 7003

cat matz.est_best_sig | awk -F $'\t' '{A[$20]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > taxids.matz_est.count
cat matz.est_best_sig | awk -F $'\t' '{A[$22]++}END{for(i in A)print i,"\t",A[i]}' | sort -k2,2nr -t "Ctr+v then tab" > species.matz_est.count
cat matz.est_best_sig | awk -F $'\t' '{A[$21]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > kingdom.matz_est.count

###
R CMD BATCH ${script_path}/mergetaxids2.R
###
grep -w "7227" unrec_chunks/unrecognizied.br_best.unrec_sig > Drosophila.unrec_nt
grep -w "7227" unrec_chunks_unrec/unrecognizied.est_best_sig > Drosophila.unrec_est


###
recog.br_best_sig
unrecognizied.br_best.rec_sig
unrecognizied.br_best.unrec_sig
symp.br_best_sig
matz.br_best_sig

setwd("/mnt/lustre_scratch_2012/Tamer/p_asteroides/blast_out")
d1=read.table("rec1_chunks/taxids.rec.count")
colnames(d1)=c("taxids","rec1")
d2=read.table("unrec_chunks/taxids.rec.count")
colnames(d2)=c("taxids","rec2")
d3=read.table("unrec_chunks/taxids.unrec.count")
colnames(d3)=c("taxids","unrec")
d4=read.table("symp/taxids.symp.count")
colnames(d4)=c("taxids","symp")
d5=read.table("matz/taxids.matz.count")
colnames(d5)=c("taxids","matz")

recog.est_best_sig
recog2.est_best_sig
unrecognizied.est_best_sig
symp.est_best_sig
matz.est_best_sig


d1=read.table("rec1_chunks/taxids.rec_est.count")
colnames(d1)=c("taxids","rec1")
d2=read.table("unrec_chunks_rec/taxids.rec2_est.count")
colnames(d2)=c("taxids","rec2")
d3=read.table("unrec_chunks_unrec/taxids.bx.count")
colnames(d3)=c("taxids","unrec")
d4=read.table("symp/taxids.symp_est.count")
colnames(d4)=c("taxids","symp")
d5=read.table("matz/taxids.matz_est.count")
colnames(d5)=c("taxids","matz")

####################
cd ${p_asteroides}/resources
for coral in "A_millepora" "A_tenuis" "A_hyacinthus" "A_digitifera" "N_vectensis" "H.vulgaris";do cd ${p_asteroides}/resources/$coral; ref=$(find . -name "*.fasta"); echo $ref; makeblastdb -in $ref -input_type fasta -dbtype nucl;done
mkdir ${p_asteroides}/symbTrans
cd ${p_asteroides}/pool
cat S_*.est.fasta K_*.est.fasta P_*.est.fasta A_fundyense.est.fasta A_monilatum.est.fasta A_temarense.est.fasta > ${p_asteroides}/symbTrans/symb.fasta
cd ${p_asteroides}/symbTrans

for coral in "A_millepora" "A_tenuis" "A_hyacinthus" "A_digitifera" "N_vectensis" "H.vulgaris";do
cd ${p_asteroides}/resources/$coral;
ref=$(find . -name "*.fasta"); echo $ref;
blastn -query ${p_asteroides}/symbTrans/symb.fasta \
   -db $ref \
   -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp" \
   -dust 'yes' -best_hit_overhang 0.25 -best_hit_score_edge 0.25 \
   -max_target_seqs 10  -out ${p_asteroides}/symbTrans/$coral ## -perc_identity 90 -qcov_hsp_perc 50
done


for coral in "A_millepora" "A_tenuis" "A_hyacinthus" "A_digitifera" "N_vectensis" "H.vulgaris";do

sort -k1,1 -k11,11g ${p_asteroides}/symbTrans/$coral | sort -u -k1,1 --merge > ${p_asteroides}/symbTrans/$coral.best
cat ${p_asteroides}/symbTrans/$coral.best | awk '$11 <= 1e-5' > ${p_asteroides}/symbTrans/$coral.best.sig

done


#####################
#####################
scp $p_asteroides/compAnalysis/Trinity.clean.201.exp.FCS.fasta tmansour@loretta.hpcf.upr.edu:/storage/prcen/coral/p_ast.assemblies.2016/. ## final all transcriptome with initial annotation
scp ${p_asteroides}/blast_out/trinityVsSymb_best.sig.fasta tmansour@loretta.hpcf.upr.edu:/storage/prcen/coral/p_ast.assemblies.2016/. ## the total zooxanrthellae transcriptome with initial annotation
scp $abundFilter/coral/{p_ast2016.fasta,ann_isoExp,LongOrfs.coral.key,*.*Stat} tmansour@loretta.hpcf.upr.edu:/storage/prcen/coral/p_ast.assemblies.2016/coral/. ## clean porites asteroids transcriptome with annotation files
scp $abundFilter/spC15/{spC15_2016.fasta,ann_isoExp,LongOrfs.spC15.key,*.*Stat} tmansour@loretta.hpcf.upr.edu:/storage/prcen/coral/p_ast.assemblies.2016/spC15/. ## clean porites asteroids transcriptome with annotation files
scp $abundFilter/cladeA/{cladeA_2016.fasta,ann_isoExp,LongOrfs.cladeA.key,*.*Stat} tmansour@loretta.hpcf.upr.edu:/storage/prcen/coral/p_ast.assemblies.2016/cladeA/. ## clean porites asteroids transcriptome with annotation files
scp $abundFilter/S_spCCMP2430/{S_spCCMP2430_2016.fasta,ann_isoExp,LongOrfs.S_spCCMP2430.key,*.*Stat} tmansour@loretta.hpcf.upr.edu:/storage/prcen/coral/p_ast.assemblies.2016/S_spCCMP2430/. ## clean porites asteroids transcriptome with annotation files




