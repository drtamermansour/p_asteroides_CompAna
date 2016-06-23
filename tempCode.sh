## Temp code
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
