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
univec_ann_exp_tran=$p_asteroides/compAnalysis/Trinity.clean.201.exp.UniVec.fasta
swiss_ann=$p_asteroides/compAnalysis/uniprot_sprot.blastx.outfmt6.sig.best.exp2.univec
LongORFs=$p_asteroides/compAnalysis/longest_orfs.pep.exp2.univec
########################
## blast out transcriptome aganist other known coral and Cnidarian sequences
bash prepResources.sh "${p_asteroides}"

module load BLAST+/2.2.30
cd ${p_asteroides}/resources/symb_trans
makeblastdb -in symb_trans.fasta -input_type fasta -dbtype nucl
blastn -query $univec_ann_exp_tran \
       -db ${p_asteroides}/resources/symb_trans/symb_trans.fasta \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp" \
       -dust 'yes' -best_hit_overhang 0.25 -best_hit_score_edge 0.25 \
       -max_target_seqs 10  -out ${p_asteroides}/blast_out/trinityVsSymb ## -perc_identity 90 -qcov_hsp_perc 50
sort -k1,1 -k11,11g ${p_asteroides}/blast_out/trinityVsSymb | sort -u -k1,1 --merge > ${p_asteroides}/blast_out/trinityVsSymb_best

cd ${p_asteroides}/blast_out
grep "^>" $univec_ann_exp_tran | wc -l        ## 868905
wc -l trinityVsSymb_best ## 188592
cat trinityVsSymb_best | awk '$11 <= 1e-5' > trinityVsSymb_best.sig
wc -l trinityVsSymb_best.sig ## 186670
cat trinityVsSymb_best.sig | awk -F '[\t.]' '{A[$2]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > symb.mapping_rates

module load QIIME/1.8.0
filter_fasta.py --input_fasta_fp $univec_ann_exp_tran --output_fasta_fp trinityVsSymb_best.sig.fasta --seq_id_fp trinityVsSymb_best.sig
grep "^>" trinityVsSymb_best.sig.fasta | wc -l  ## 186670

filter_fasta.py --input_fasta_fp $univec_ann_exp_tran --output_fasta_fp noSymb.fasta --seq_id_fp trinityVsSymb_best.sig --negate
grep "^>" noSymb.fasta | wc -l  ## 682235
noSymb_transcriptome=${p_asteroides}/blast_out/noSymb.fasta
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
wc -l trinityVscoralTrans_best ## 93078
cat trinityVscoralTrans_best | awk '$11 <= 1e-5' > trinityVscoralTrans_best.sig
wc -l trinityVscoralTrans_best.sig ## 92454
cat trinityVscoralTrans_best.sig | awk -F '[\t.]' '{A[$2]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > coralTrans.mapping_rates

filter_fasta.py --input_fasta_fp $noSymb_transcriptome --output_fasta_fp trinityVscoralTrans_best.sig.fasta --seq_id_fp trinityVscoralTrans_best.sig
grep "^>" trinityVscoralTrans_best.sig.fasta | wc -l  ## 92454

filter_fasta.py --input_fasta_fp $noSymb_transcriptome --output_fasta_fp noSymb_noCoralTrans.fasta --seq_id_fp trinityVscoralTrans_best.sig --negate
grep "^>" noSymb_noCoralTrans.fasta | wc -l  ## 589781
noSymb_noCoralTrans_transcriptome=${p_asteroides}/blast_out/noSymb_noCoralTrans.fasta
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
wc -l trinityVscoralGenomic_best ## 38701
cat trinityVscoralGenomic_best | awk '$11 <= 1e-5' > trinityVscoralGenomic_best.sig
wc -l trinityVscoralGenomic_best.sig ## 37448
cat trinityVscoralGenomic_best.sig | awk -F '[\t.]' '{A[$2]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > coralGenomics.mapping_rates

filter_fasta.py --input_fasta_fp $noSymb_noCoralTrans_transcriptome --output_fasta_fp trinityVscoralGenomic_best.sig.fasta --seq_id_fp trinityVscoralGenomic_best.sig
grep "^>" trinityVscoralGenomic_best.sig.fasta | wc -l  ## 37448

filter_fasta.py --input_fasta_fp $noSymb_noCoralTrans_transcriptome --output_fasta_fp unrecognized.fasta --seq_id_fp trinityVscoralGenomic_best.sig --negate
grep "^>" unrecognized.fasta | wc -l  ## 552333
unrecognized=${p_asteroides}/blast_out/unrecognized.fasta

cat trinityVscoralTrans_best.sig trinityVscoralGenomic_best.sig > trinityVscoral_best.sig
cat trinityVsSymb_best.sig trinityVscoral_best.sig > recognized_best.sig
cat trinityVscoralTrans_best.sig.fasta trinityVscoralGenomic_best.sig.fasta > coral_transcriptome.fasta
cat trinityVsSymb_best.sig.fasta coral_transcriptome.fasta > recognized.fasta
coral_transcriptome=${p_asteroides}/blast_out/coral_transcriptome.fasta
recognized_transcriptome=${p_asteroides}/blast_out/recognized.fasta
##################
## Abundance estimation
cd $compAnalysis
qsub -v index="salmon_index",transcriptome="$univec_ann_exp_tran" ${script_path}/salmonIndex.sh

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

module load R/3.0.1
while read identifier;do
  echo $(pwd) $identifier
  Rscript ${script_path}/calcTPM_tis.R "$(pwd)" "$identifier" "transcripts.lengthes" "$gene_transcript_map" >> targets_list
done < $identifiers
bash $script_path/abund_est.sh
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

## isoforms expressed in larva (569961) & has sig hit against symb. (186670)
comm -12 <(cat $abundFilter/larva_isoformTPM | tail -n+2 | awk '{print $1}' |sort) <(cat $p_asteroides/blast_out/trinityVsSymb_best.sig | awk '{print $1}' | sort) > larva_isoformTPM_vs_trinityVsSymb  ## 179948
comm -12 <(cat $swiss_ann | tail -n+2 | awk '{print $1}' |sort) larva_isoformTPM_vs_trinityVsSymb > larva_isoformTPM_vs_trinityVsSymb.ann  ## 49900

## isoforms expressed in larva (569961) & has sig hit against coral seq (129902)
comm -12 <(cat $abundFilter/larva_isoformTPM | tail -n+2 | awk '{print $1}' |sort) <(cat $p_asteroides/blast_out/trinityVscoral_best.sig | awk '{print $1}' | sort) > larva_isoformTPM_vs_trinityVscoral  ## 126038
comm -12 <(cat $swiss_ann | tail -n+2 | awk '{print $1}' |sort) larva_isoformTPM_vs_trinityVscoral > larva_isoformTPM_vs_trinityVscoral.ann  ## 24292

## isoforms expressed in adult only (298944) & has sig hit against symb. (186670)
comm -12 <(cat $abundFilter/adultOnly_isoformTPM | tail -n+2 | awk '{print $1}' |sort) <(cat $p_asteroides/blast_out/trinityVsSymb_best.sig | awk '{print $1}' | sort) > adultOnly_isoformTPM_vs_trinityVsSymb  ## 6722
comm -12 <(cat $swiss_ann | tail -n+2 | awk '{print $1}' |sort) adultOnly_isoformTPM_vs_trinityVsSymb > adultOnly_isoformTPM_vs_trinityVsSymb.ann  ## 2534

## isoforms expressed in adult only (298944)  & has sig hit against coral seq (129902)
comm -12 <(cat $abundFilter/adultOnly_isoformTPM | tail -n+2 | awk '{print $1}' |sort) <(cat $p_asteroides/blast_out/trinityVscoral_best.sig | awk '{print $1}' | sort) > adultOnly_isoformTPM_vs_trinityVscoral  ## 3864
comm -12 <(cat $swiss_ann | tail -n+2 | awk '{print $1}' |sort) adultOnly_isoformTPM_vs_trinityVscoral > adultOnly_isoformTPM_vs_trinityVscoral.ann  ## 1146

#####
cd $abundFilter
mkdir coral
cat ${p_asteroides}/blast_out/trinityVscoral_best.sig | awk '{print $1}' > coral_transIDs
head -n1 transcripts.lengthes > coral/transcripts.lengthes
grep -w -F -f coral_transIDs transcripts.lengthes >> coral/transcripts.lengthes
grep -w -F -f coral_transIDs gene_transcript_map > coral/gene_transcript_map
for f in *.quant.counts;do
  head -n1 $f > coral/$f
  grep -w -F -f coral_transIDs $f >> coral/$f
done
cd coral

module load R/3.0.1
while read identifier;do
  echo $(pwd) $identifier
  Rscript ${script_path}/calcTPM_tis.R "$(pwd)" "$identifier" "transcripts.lengthes" "gene_transcript_map" >> targets_list
done < $identifiers
bash $script_path/abund_est.sh

Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1], header=T,row.names=NULL);head(data1); data2=read.table(args[2], header=T,row.names=NULL,sep="\t",quote="");head(data2); data3=read.table(args[3], header=F,row.names=NULL,sep="\t");head(data3);ann_isoExp=merge(data1,data2[,c(1,2,17)],by.x="geneName",by.y="qseqid");ann_isoExp2=merge(ann_isoExp,data3[,c(1,2)],by.x="geneName",by.y="V1");write.table(ann_isoExp2,"ann_isoExp", sep="\t", quote=F, row.names=F, col.names=T);' $abundFilter/coral/exp_isoformTPM $swiss_ann ${p_asteroides}/blast_out/trinityVscoral_best.sig

filter_fasta.py --input_fasta_fp $univec_ann_exp_tran --output_fasta_fp p_ast2016.fasta --seq_id_fp ../coral_transIDs
p_ast2016_trans=$abundFilter/coral/p_ast2016.fasta

filter_fasta.py --input_fasta_fp $LongOrfs --output_fasta_fp longest_orfs.pep --seq_id_fp ../coral_transIDs

while read gene;do
  grep -A1 $gene":" $LongOrfs.complete
done < ../coral_transIDs > longest_orfs.pep.complete  ## 26806
grep "^>" longest_orfs.pep.complete | awk -F '[>|]' '{print $2"|"$3}' | sort | uniq | wc -l ## 17282

##############
cd $abundFilter
mkdir spC15
grep "S_spC15.est" ${p_asteroides}/blast_out/trinityVsSymb_best.sig | awk '{print $1}' > spC15_transIDs
head -n1 transcripts.lengthes > spC15/transcripts.lengthes
grep -w -F -f spC15_transIDs transcripts.lengthes >> spC15/transcripts.lengthes
grep -w -F -f spC15_transIDs gene_transcript_map > spC15/gene_transcript_map
for f in *.quant.counts;do
  head -n1 $f > spC15/$f
  grep -w -F -f spC15_transIDs $f >> spC15/$f
done
cd spC15

module load R/3.0.1
while read identifier;do
  echo $(pwd) $identifier
  Rscript ${script_path}/calcTPM_tis.R "$(pwd)" "$identifier" "transcripts.lengthes" "gene_transcript_map" >> targets_list
done < $identifiers
bash $script_path/abund_est.sh


filter_fasta.py --input_fasta_fp $univec_ann_exp_tran --output_fasta_fp spC15_2016.fasta --seq_id_fp ../spC15_transIDs
##############
cd $abundFilter
mkdir cladeA
grep "S_cladeA.est" ${p_asteroides}/blast_out/trinityVsSymb_best.sig | awk '{print $1}' > cladeA_transIDs
head -n1 transcripts.lengthes > cladeA/transcripts.lengthes
grep -w -F -f cladeA_transIDs transcripts.lengthes >> cladeA/transcripts.lengthes
grep -w -F -f cladeA_transIDs gene_transcript_map > cladeA/gene_transcript_map
for f in *.quant.counts;do
  head -n1 $f > cladeA/$f
  grep -w -F -f cladeA_transIDs $f >> cladeA/$f
done
cd cladeA

module load R/3.0.1
while read identifier;do
  echo $(pwd) $identifier
  Rscript ${script_path}/calcTPM_tis.R "$(pwd)" "$identifier" "transcripts.lengthes" "gene_transcript_map" >> targets_list
done < $identifiers
bash $script_path/abund_est.sh

filter_fasta.py --input_fasta_fp $univec_ann_exp_tran --output_fasta_fp cladeA_2016.fasta --seq_id_fp ../cladeA_transIDs

##############
## Assessement of the transcriptome
cd $compAnalysis
module load Bioperl/1.6.923
perl ${script_path}/seq_stats.pl $p_ast2016_trans > $p_ast2016_trans.MatzStat

module load trinity/6.0.2
TrinityStats.pl $p_ast2016_trans > $p_ast2016_trans.TrinityStat

## calc the the no of Complete ORFs
grep "type:complete" $LongOrfs | wc -l  ##225219
grep -A1 "type:complete" $LongOrfs | grep -v "^--" > $LongOrfs.complete ## 225219
grep "^>" $LongOrfs.complete | awk -F '[>|]' '{print $2"|"$3}' | sort | uniq | wc -l ## 124271
while read gene;do grep $gene"|" $LongOrfs; done < <(cat $compAnalysis/unexpIDs $p_asteroides/UniVec/excludeIDs_orig | sort |uniq) > $LongOrfs.exclude
grep "type:complete" $LongOrfs.exclude | wc -l ## 810
grep "type:complete" $LongOrfs.exclude > $LongOrfs.exclude.complete ## 810
grep "^>" $LongOrfs.exclude.complete | awk -F '[>|]' '{print $2"|"$3}' | sort | uniq | wc -l ## 449
## Total no of complte ORFs (one transcript may have many ORFs): 225219 - 810 = 224409
## no of uniqe transcripts with complete ORFS: 124271 - 449 = 123822


## instaling and running assemblathon2
#cd ${script_path}
#git clone https://github.com/ucdavis-bioinformatics/assemblathon2-analysis.git
#cd assemblathon2-analysis ## you have to be in this folder so that the file can use the FAlite.pm script
#perl ./assemblathon_stats.pl ${compAnalysis}/Trinity.fasta.clean > ${compAnalysis}/trinity_assemblathon2.stat
#perl ./assemblathon_stats.pl ${p_asteroides}/data/Porites.Astreoides.uniprot2013.fa > ${p_asteroides}/data/trinity_assemblathon2.stat

## http://deweylab.biostat.wisc.edu/detonate/vignette.html
#cd $ass_dir
#module load DETONATE/1.8.1
#lf="${lf_files[*]}"
#rt="${rt_files[*]}"
#rsem-eval-calculate-score --paired-end $(echo ${lf[*]} | tr ' ' ',') $(echo ${rt[*]} | tr ' ' ',') \
#			  trinity_out_dir/seqclean/Trinity.fasta.clean \
#			  Trinity.fasta.clean.detonate \
#			  400 \
#			  --transcript-length-parameters rsem-eval/true_transcript_length_distribution/mouse.txt \
#			  -p 16

#####################
## TSA Submission Guide
## http://www.ncbi.nlm.nih.gov/genbank/tsaguide
## create ASN file
## https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2#tbl
module load tbl2asn/20150331
mkdir ${p_asteroides}/ASN
cd ${p_asteroides}/ASN
## download the template.sbt & assembly.cmt
cat $univec_ann_exp_tran | awk '{print $1}' > allTrans.fsa
qsub $script_path/createASN.sh
##############
## copy the assemblies to the other server
cd ${p_asteroides}/ASN
scp allTrans.fsa tmansour@loretta.hpcf.upr.edu:/storage/prcen/coral/NCBI_submission/allTrans.fsa
scp allTrans.sqn tmansour@loretta.hpcf.upr.edu:/storage/prcen/coral/NCBI_submission/allTrans.sqn

scp $p_asteroides/c_abundFilter/coral/p_ast2016.fasta tmansour@loretta.hpcf.upr.edu:/export/home/tmansour/p_ast.assemblies.2016/.
scp $p_asteroides/c_abundFilter/spC15/spC15_2016.fasta tmansour@loretta.hpcf.upr.edu:/export/home/tmansour/p_ast.assemblies.2016/.
scp $p_asteroides/c_abundFilter/cladeA/cladeA_2016.fasta tmansour@loretta.hpcf.upr.edu:/export/home/tmansour/p_ast.assemblies.2016/.
##############
## compare the new Assembly with the older assembly
bash $script_path/compareVSoldTrans.sh
##################
## blast against the published pastreoids transcriptome
# a Matz lab paper (http://onlinelibrary.wiley.com/doi/10.1111/mec.12390/abstract)
# The annotated transcriptome (http://www.bio.utexas.edu/research/matz_lab/matzlab/Data.html).
# Porites astreoides (adult, Symbiodinium-specific reads excluded)
cd ${p_asteroides}/resources
wget https://dl.dropboxusercontent.com/u/37523721/pastreoides_transcriptome_july2014.zip
mkdir ${p_asteroides}/resources/P_ast.transcriptome
unzip pastreoides_transcriptome_july2014.zip -d P_ast.transcriptome
mkdir ${p_asteroides}/resources/P_ast.transcriptome/P_ast.BlastDB
cd ${p_asteroides}/resources/P_ast.transcriptome/P_ast.BlastDB
cp ../past.fasta .
module load BLAST+/2.2.30
makeblastdb -in past.fasta -input_type fasta -dbtype nucl
blastn -query ${p_asteroides}/blast_out/unrecognizied_est.fasta \
       -db past.fasta \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp" \
       -dust 'yes' -best_hit_overhang 0.25 -best_hit_score_edge 0.25 -evalue 1e-5 \
       -max_target_seqs 10  -out unrecogniziedVspublishedAss ## -perc_identity 90 -qcov_hsp_perc 50
sort -k1,1 -k12,12nr -k11,11n unrecogniziedVspublishedAss | sort -u -k1,1 --merge > unrecogniziedVspublishedAss_best
wc -l unrecogniziedVspublishedAss_best ## 153020  (out of 809187)
cat unrecogniziedVspublishedAss_best | awk '$11 <= 1e-5' | wc -l        ##  153020

module load QIIME/1.8.0  ## it might be better to start new screen to avoid module conflict on HPC
filter_fasta.py --input_fasta_fp ${p_asteroides}/blast_out/unrecognizied_est.fasta --output_fasta_fp ${p_asteroides}/blast_out/recog.uncontaminated2.fasta --seq_id_fp unrecogniziedVspublishedAss_best

filter_fasta.py --input_fasta_fp ${p_asteroides}/blast_out/unrecognizied_est.fasta --output_fasta_fp ${p_asteroides}/blast_out/unrecognizied_est2.fasta --seq_id_fp unrecogniziedVspublishedAss_best --negate

cd ${p_asteroides}/blast_out/
cat recog.uncontaminated.fasta recog.uncontaminated2.fasta > recog.uncontaminated.total.fasta
mkdir HiconfTrans.total_BlastDB && cd HiconfTrans.total_BlastDB
cp ../recog.uncontaminated.total.fasta .
module load BLAST+/2.2.30
makeblastdb -in recog.uncontaminated.total.fasta -input_type fasta -dbtype nucl

####################
grep "^>" ${p_asteroides}/resources/P_ast.transcriptome/P_ast.BlastDB/past.fasta | wc -l ## 30740

cd ${p_asteroides}/blast_out/HiconfTrans_BlastDB/
blastn -query ${p_asteroides}/resources/P_ast.transcriptome/P_ast.BlastDB/past.fasta \
       -db recog.uncontaminated.fasta \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp" \
       -dust 'yes' -best_hit_overhang 0.25 -best_hit_score_edge 0.25 -evalue 1e-5 \
       -max_target_seqs 10  -out publishedAssVsrecog ## -perc_identity 90 -qcov_hsp_perc 50
sort -k1,1 -k12,12nr -k11,11n  publishedAssVsrecog | sort -u -k1,1 --merge > publishedAssVsrecog_best 
wc -l publishedAssVsrecog_best ## 16503 (out of 30740)
cat publishedAssVsrecog_best | awk '$11 <= 1e-5' | wc -l        ## 16503 


cd ${p_asteroides}/blast_out/HiconfTrans.total_BlastDB/
blastn -query ${p_asteroides}/resources/P_ast.transcriptome/P_ast.BlastDB/past.fasta \
       -db recog.uncontaminated.total.fasta \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp" \
       -dust 'yes' -best_hit_overhang 0.25 -best_hit_score_edge 0.25 -evalue 1e-5 \
       -max_target_seqs 10  -out publishedAssVsrecog.total ## -perc_identity 90 -qcov_hsp_perc 50
sort -k1,1 -k12,12nr -k11,11n  publishedAssVsrecog.total | sort -u -k1,1 --merge > publishedAssVsrecog.total_best 
wc -l publishedAssVsrecog.total_best ## 29929 (out of 30740)
cat publishedAssVsrecog.total_best | awk '$11 <= 1e-5' | wc -l        ## 29929 

cat publishedAssVsrecog.total_best | awk '$13 < $14' > publishedAssVsrecog.total_best_better        
wc -l publishedAssVsrecog.total_best_better ## 27564
cat publishedAssVsrecog.total_best | awk '$13 == $14' > publishedAssVsrecog.total_best_equal        
wc -l publishedAssVsrecog.total_best_equal ## 7
cat publishedAssVsrecog.total_best | awk '$13 > $14' > publishedAssVsrecog.total_best_less        
wc -l publishedAssVsrecog.total_best_less ## 2354
cat publishedAssVsrecog.total_best_better | awk '{ sum+=$14} END {print sum}' ## 92886988
cat publishedAssVsrecog.total_best_better | awk '{ sum+=$13} END {print sum}' ## 14688495 (i.e. difference of 78198493 =~78Mb)
cat publishedAssVsrecog.total_best_less | awk '{ sum+=$14} END {print sum}' ## 1300368
cat publishedAssVsrecog.total_best_less | awk '{ sum+=$13} END {print sum}' ## 1888017 (i.e. difference of 587649 =~0.6Mb)

module load QIIME/1.8.0  ## it might be better to start new screen to avoid module conflict on HPC
filter_fasta.py --input_fasta_fp ${p_asteroides}/resources/P_ast.transcriptome/P_ast.BlastDB/past.fasta --output_fasta_fp publishedunrecognizied_est.fasta --seq_id_fp publishedAssVsrecog.total_best --negate

module load Bioperl/1.6.923
perl ${script_path}/seq_stats.pl recog.uncontaminated.total.fasta >  recog.uncontaminated.total.fasta.MatzStat
perl ${script_path}/seq_stats.pl ${p_asteroides}/resources/P_ast.transcriptome/P_ast.BlastDB/past.fasta > ${p_asteroides}/resources/P_ast.transcriptome/P_ast.BlastDB/past.fasta.MatzStat

module load trinity/6.0.2
TrinityStats.pl recog.uncontaminated.total.fasta > recog.uncontaminated.total.fasta.TrinityStat

#cd ${p_asteroides}/blast_out/uncontaminated_BlastDB/
#blastn -query ${p_asteroides}/resources/P_ast.transcriptome/P_ast.BlastDB/past.fasta \
#       -db uncontaminated.fasta \
#       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp" \
#       -dust 'yes' -best_hit_overhang 0.25 -best_hit_score_edge 0.25 -evalue 1e-5 \
#       -max_target_seqs 10  -out publishedAssVsuncontaminated ## -perc_identity 90 -qcov_hsp_perc 50
#sort -k1,1 -k12,12nr -k11,11n  publishedAssVsuncontaminated | sort -u -k1,1 --merge > publishedAssVsuncontaminated_best 
#wc -l publishedAssVsuncontaminated_best ## 29932 (out of 30740)



##############
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







