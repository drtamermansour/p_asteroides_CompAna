p_asteroides="$1"

mkdir ${p_asteroides}/resources/coral_trans
coral_trans="${p_asteroides}/resources/coral_trans"
mkdir ${p_asteroides}/resources/coral_genomic
coral_genomic="${p_asteroides}/resources/coral_genomic"
mkdir ${p_asteroides}/resources/symb_trans
symb_trans="${p_asteroides}/resources/symb_trans"

## Acropora digitifera
## http://www.nature.com/nature/journal/v476/n7360/full/nature10249.html
## http://marinegenomics.oist.jp/genomes/downloads?project_id=3 ## also for genome download
mkdir -p ${p_asteroides}/resources/coral_trans/A_digitifera
cd ${p_asteroides}/resources/coral_trans/A_digitifera
wget http://marinegenomics.oist.jp/coral/download/adi_transcriptome_assembly.v1.fa.gz
gunzip -c adi_transcriptome_assembly.v1.fa.gz > A_digit.est.fasta
sed 's/^>/>A_digitifera.est./' A_digit.est.fasta > $coral_trans/A_digit.est.fasta

## Acropora millepora (larva and adult) initally publised by Moya Mol Ecol 2012, 21: 2440
## http://onlinelibrary.wiley.com/doi/10.1111/j.1365-294X.2012.05554.x/abstract
## The download is avialble from Matz lab who re-annotated the transcriptome
## The new download links can be obtained from: http://www.bio.utexas.edu/research/matz_lab/matzlab/Data.html
mkdir ${p_asteroides}/resources/coral_trans/A_millepora
cd ${p_asteroides}/resources/coral_trans/A_millepora
wget https://www.dropbox.com/s/w5nb1siz7uu7v5b/Amillepora_transcriptome_July2014.zip
unzip Amillepora_transcriptome_July2014.zip -d .
sed 's/^>/>A_millepora.est./' amil_apr2014/amil.fasta > $coral_trans/A_millepora.est.fasta

## Acropora hyacinthus (larvae) by matzlab
mkdir ${p_asteroides}/resources/coral_trans/A_hyacinthus
cd ${p_asteroides}/resources/coral_trans/A_hyacinthus
wget https://dl.dropboxusercontent.com/u/37523721/ahyacinthus_transcriptome_july2014.zip
unzip ahyacinthus_transcriptome_july2014.zip -d .
sed 's/^>/>A_hyacinthus.est./' ahya.fasta > $coral_trans/A_hyacinthus.est.fasta

## Acropora tenuis (larvae and aposymbiotic recruits)
mkdir ${p_asteroides}/resources/coral_trans/A_tenuis
cd ${p_asteroides}/resources/coral_trans/A_tenuis
wget https://dl.dropboxusercontent.com/u/37523721/atenuis_transcriptome_july2014.zip
unzip atenuis_transcriptome_july2014.zip -d .
sed 's/^>/>A_tenuis.est./' aten.fasta > $coral_trans/A_tenuis.est.fasta

## Porites astreoides (adult, Symbiodinium specific reads excluded)
mkdir ${p_asteroides}/resources/coral_trans/P_astreoides
cd ${p_asteroides}/resources/coral_trans/P_astreoides
wget https://dl.dropboxusercontent.com/u/37523721/pastreoides_transcriptome_july2014.zip
unzip pastreoides_transcriptome_july2014.zip -d .
sed 's/^>/>P_astreoides.est./' past.fasta > $coral_trans/P_astreoides.est.fasta

## Nematostella vectensis by JGI
## http://cnidarians.bu.edu/stellabase/cgi-bin/libraries.cgi ## also for genome download
## http://genome.jgi.doe.gov/Nemve1/Nemve1.home.html
mkdir ${p_asteroides}/resources/coral_trans/N_vectensis
cd ${p_asteroides}/resources/coral_trans/N_vectensis
wget http://cnidarians.bu.edu/stellabase/assembly/NvT1.fasta
sed 's/^>/>N_vectensis.est./' NvT1.fasta > $coral_trans/N_vectensis.est.fasta

## downloads from dbEST
## https://www.biostars.org/p/518/
## http://www.ncbi.nlm.nih.gov/books/NBK25501/
## my lines to retrive ests of hydra vulgaris from dbest
## http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucest&term="Hydra vulgaris"[Organism]&&usehistory=y
## http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucest&query_key=1&WebEnv=NCID_1_1188506488_130.14.18.34_9001_1454075766_1185569046_0MetA0_S_MegaStore_F_1&rettype=fasta&retmode=text
## running E-utilities at command line
## http://www.ncbi.nlm.nih.gov/books/NBK179288/
mkdir ${p_asteroides}/resources/coral_trans/H_vulgaris
cd ${p_asteroides}/resources/coral_trans/H_vulgaris
esearch -db nucest -query "Hydra vulgaris[Organism]" | efetch -format fasta -mode text > H_vulgaris.est.fasta
sed 's/^>/>H_vulgaris.est./' H_vulgaris.est.fasta > $coral_trans/H_vulgaris.est.fasta

cd $coral_trans
cat *.fasta > coral_trans.fasta
 
################################################################################################
## Acropora millepora from Trace archive
mkdir ${p_asteroides}/resources/coral_genomic/A_milleporaSRA
cd ${p_asteroides}/resources/coral_genomic/A_milleporaSRA
wget ftp://ftp.ncbi.nlm.nih.gov/pub/TraceDB/acropora_millepora/fasta.acropora_millepora.001.gz
gunzip fasta.acropora_millepora.001.gz
sed 's/^>/>A_milleporaSRA./' fasta.acropora_millepora.001 > $coral_genomic/A_milleporaSRA.genomic.fasta

## Acropora palmata from Trace archive
mkdir ${p_asteroides}/resources/coral_genomic/A_palmataSRA
cd ${p_asteroides}/resources/coral_genomic/A_palmataSRA
wget ftp://ftp.ncbi.nlm.nih.gov/pub/TraceDB/acropora_palmata/fasta.acropora_palmata.001.gz
gunzip fasta.acropora_palmata.001.gz
sed 's/^>/>A_palmataSRA.est./' fasta.acropora_palmata.001 > $coral_genomic/A_palmataSRA.genomic.fasta

## Porites lobata from Trace archive
mkdir ${p_asteroides}/resources/coral_genomic/P_lobataSRA
cd ${p_asteroides}/resources/coral_genomic/P_lobataSRA
wget ftp://ftp.ncbi.nlm.nih.gov/pub/TraceDB/porites_lobata/fasta.porites_lobata.001.gz
gunzip fasta.porites_lobata.001.gz
sed 's/^>/>P_lobataSRA./' fasta.porites_lobata.001 > $coral_genomic/P_lobataSRA.genomic.fasta

## Montastraea faveolata from Trace archive
mkdir ${p_asteroides}/resources/coral_genomic/M_faveolataSRA
cd ${p_asteroides}/resources/coral_genomic/M_faveolataSRA
wget ftp://ftp.ncbi.nlm.nih.gov/pub/TraceDB/montastraea_faveolata/fasta.montastraea_faveolata.001.gz
gunzip fasta.montastraea_faveolata.001.gz
sed 's/^>/>M_faveolataSRA./' fasta.montastraea_faveolata.001 > $coral_genomic/M_faveolataSRA.genomic.fasta

## Acropora digitifera genome
mkdir ${p_asteroides}/resources/coral_genomic/A_digitiferaGenome
cd ${p_asteroides}/resources/coral_genomic/A_digitiferaGenome
wget http://marinegenomics.oist.jp/coral/download/adi_v1.0.scaffold.fa.gz
gunzip adi_v1.0.scaffold.fa.gz
sed 's/^>/>A_digitiferaGenome./' adi_v1.0.scaffold.fa > $coral_genomic/A_digitifera.genomic.fasta

## Nematostella vectensis genome by JGI
mkdir ${p_asteroides}/resources/coral_genomic/N_vectensisGenome
cd ${p_asteroides}/resources/coral_genomic/N_vectensisGenome
wget ftp://ftp.jgi-psf.org/pub/JGI_data/Nematostella_vectensis/v1.0/assembly/Nemve1.allmasked.gz
gunzip Nemve1.allmasked.gz
sed 's/^>/>N_vectensisGenome./' Nemve1.allmasked > $coral_genomic/N_vectensis.genomic.fasta

cd ${p_asteroides}/resources/coral_genomic
cat *.fasta > coral_genomic.fasta

################################################################################################
## blast against the published Symbiodinium transcriptomes
# Symbiodinium transcriptomes: http://www-ncbi-nlm-gov.proxy1.cl.msu.edu/pubmed/22529998
# download page: http://medinalab.org/zoox/
mkdir -p ${p_asteroides}/resources/Symbiodinium/cladeA
cd ${p_asteroides}/resources/Symbiodinium/cladeA
wget http://medinalab.org/zoox/kb8_assembly.fasta.bz2
bzip2 -d *.bz2
sed 's/^>/>S_cladeA.est./' kb8_assembly.fasta > ${p_asteroides}/resources/symb_trans/S_cladeA.est.fasta

mkdir -p ${p_asteroides}/resources/Symbiodinium/cladeB
cd ${p_asteroides}/resources/Symbiodinium/cladeB
wget http://medinalab.org/zoox/mf105_assembly.fasta.bz2
bzip2 -d *.bz2
sed 's/^>/>S_cladeB.est./' mf105_assembly.fasta > ${p_asteroides}/resources/symb_trans/S_cladeB.est.fasta

# download page: http://marinegenomics.oist.jp/genomes/gallery
mkdir -p ${p_asteroides}/resources/Symbiodinium/minutum
cd ${p_asteroides}/resources/Symbiodinium/minutum
wget http://marinegenomics.oist.jp/symb/download/symbB1_v1.0.transcriptome_assembly.fa.gz
gunzip symbB1_v1.0.transcriptome_assembly.fa.gz
sed 's/^>/>S_minutum.est./' symbB1_v1.0.transcriptome_assembly.fa > ${p_asteroides}/resources/symb_trans/S_minutum.est.fasta

## Project: Marine Microbial Eukaryote Transcriptome Sequencing Project (MMETSP)
# http://data.imicrobe.us/project/view/104
# download page: ftp://ftp.imicrobe.us/camera/combined_assemblies (genome & Annotations available)
mkdir ${p_asteroides}/resources/Symbiodinium/kawagutii
cd ${p_asteroides}/resources/Symbiodinium/kawagutii
wget ftp://ftp.imicrobe.us/camera/combined_assemblies/Symbiodinium-kawagutii-CCMP2468.nt.fa.gz
gunzip Symbiodinium-kawagutii-CCMP2468.nt.fa.gz
sed 's/^>/>S_kawagutii.est./' Symbiodinium-kawagutii-CCMP2468.nt.fa > ${p_asteroides}/resources/symb_trans/S_kawagutii.est.fasta

mkdir ${p_asteroides}/resources/Symbiodinium/sp_C1
cd ${p_asteroides}/resources/Symbiodinium/sp_C1
wget ftp://ftp.imicrobe.us/camera/combined_assemblies/Symbiodinium-sp-C1.nt.fa.gz
gunzip Symbiodinium-sp-C1.nt.fa.gz
sed 's/^>/>S_spC1.est./' Symbiodinium-sp-C1.nt.fa > ${p_asteroides}/resources/symb_trans/S_spC1.est.fasta

mkdir ${p_asteroides}/resources/Symbiodinium/sp_C15
cd ${p_asteroides}/resources/Symbiodinium/sp_C15
wget ftp://ftp.imicrobe.us/camera/combined_assemblies/Symbiodinium-sp-C15.nt.fa.gz
gunzip Symbiodinium-sp-C15.nt.fa.gz
sed 's/^>/>S_spC15.est./' Symbiodinium-sp-C15.nt.fa > ${p_asteroides}/resources/symb_trans/S_spC15.est.fasta

mkdir ${p_asteroides}/resources/Symbiodinium/sp_CCMP2430
cd ${p_asteroides}/resources/Symbiodinium/sp_CCMP2430
wget ftp://ftp.imicrobe.us/camera/combined_assemblies/Symbiodinium-sp-CCMP2430.nt.fa.gz
gunzip Symbiodinium-sp-CCMP2430.nt.fa.gz
sed 's/^>/>S_spCCMP2430.est./' Symbiodinium-sp-CCMP2430.nt.fa > ${p_asteroides}/resources/symb_trans/S_spCCMP2430.est.fasta

mkdir ${p_asteroides}/resources/Symbiodinium/sp_Mp
cd ${p_asteroides}/resources/Symbiodinium/sp_Mp
wget ftp://ftp.imicrobe.us/camera/combined_assemblies/Symbiodinium-sp-Mp.nt.fa.gz
gunzip Symbiodinium-sp-Mp.nt.fa.gz
sed 's/^>/>S_spMp.est./' Symbiodinium-sp-Mp.nt.fa > ${p_asteroides}/resources/symb_trans/S_spMp.est.fasta

mkdir -p ${p_asteroides}/resources/Alexandrium/fundyense
cd ${p_asteroides}/resources/Alexandrium/fundyense
wget ftp://ftp.imicrobe.us/camera/combined_assemblies/Alexandrium-fundyense-CCMP1719.nt.fa.gz
gunzip Alexandrium-fundyense-CCMP1719.nt.fa.gz
sed 's/^>/>A_fundyense.est./' Alexandrium-fundyense-CCMP1719.nt.fa > ${p_asteroides}/resources/symb_trans/A_fundyense.est.fasta

mkdir ${p_asteroides}/resources/Alexandrium/monilatum
cd ${p_asteroides}/resources/Alexandrium/monilatum
wget ftp://ftp.imicrobe.us/camera/combined_assemblies/Alexandrium-monilatum-CCMP3105.nt.fa.gz
gunzip Alexandrium-monilatum-CCMP3105.nt.fa.gz
sed 's/^>/>A_monilatum.est./' Alexandrium-monilatum-CCMP3105.nt.fa > ${p_asteroides}/resources/symb_trans/A_monilatum.est.fasta

mkdir -p ${p_asteroides}/resources/Alexandrium/temarense
cd ${p_asteroides}/resources/Alexandrium/temarense
wget ftp://ftp.imicrobe.us/camera/combined_assemblies/Alexandrium-temarense-CCMP1771.nt.fa.gz
gunzip Alexandrium-temarense-CCMP1771.nt.fa.gz
sed 's/^>/>A_temarense.est./' Alexandrium-temarense-CCMP1771.nt.fa > ${p_asteroides}/resources/symb_trans/A_temarense.est.fasta

mkdir -p ${p_asteroides}/resources/Peridinium/aciculiferum
cd ${p_asteroides}/resources/Peridinium/aciculiferum
wget ftp://ftp.imicrobe.us/camera/combined_assemblies/Peridinium-aciculiferum-PAER_2.nt.fa.gz
gunzip Peridinium-aciculiferum-PAER_2.nt.fa.gz
sed 's/^>/>P_aciculiferum.est./' Peridinium-aciculiferum-PAER_2.nt.fa > ${p_asteroides}/resources/symb_trans/P_aciculiferum.est.fasta

mkdir -p ${p_asteroides}/resources/Karenia/brevis_CCMP2229
cd ${p_asteroides}/resources/Karenia/brevis_CCMP2229
wget ftp://ftp.imicrobe.us/camera/combined_assemblies/Karenia-brevis-CCMP2229.nt.fa.gz
gunzip Karenia-brevis-CCMP2229.nt.fa.gz
sed 's/^>/>K_brevisCCMP2229.est./' Karenia-brevis-CCMP2229.nt.fa > ${p_asteroides}/resources/symb_trans/K_brevisCCMP2229.est.fasta

mkdir -p ${p_asteroides}/resources/Karenia/brevis_SP1
cd ${p_asteroides}/resources/Karenia/brevis_SP1
wget ftp://ftp.imicrobe.us/camera/combined_assemblies/Karenia-brevis-SP1.nt.fa.gz
gunzip Karenia-brevis-SP1.nt.fa.gz
sed 's/^>/>K_brevisSP1.est./' Karenia-brevis-SP1.nt.fa > ${p_asteroides}/resources/symb_trans/K_brevisSP1.est.fasta

mkdir -p ${p_asteroides}/resources/Karenia/brevis_SP3
cd ${p_asteroides}/resources/Karenia/brevis_SP3
wget ftp://ftp.imicrobe.us/camera/combined_assemblies/Karenia-brevis-SP3.nt.fa.gz
gunzip Karenia-brevis-SP3.nt.fa.gz
sed 's/^>/>K_brevisSP3.est./' Karenia-brevis-SP3.nt.fa > ${p_asteroides}/resources/symb_trans/K_brevisSP3.est.fasta

mkdir -p ${p_asteroides}/resources/Karenia/brevis_Wilson
cd ${p_asteroides}/resources/Karenia/brevis_Wilson
wget ftp://ftp.imicrobe.us/camera/combined_assemblies/Karenia-brevis-Wilson.nt.fa.gz
gunzip Karenia-brevis-Wilson.nt.fa.gz
sed 's/^>/>K_brevisWilson.est./' Karenia-brevis-Wilson.nt.fa > ${p_asteroides}/resources/symb_trans/K_brevisWilson.est.fasta

mkdir -p ${p_asteroides}/resources/Prorocentrum/minimum_CCMP1329
cd ${p_asteroides}/resources/Prorocentrum/minimum_CCMP1329
wget ftp://ftp.imicrobe.us/camera/combined_assemblies/Prorocentrum-minimum-CCMP1329.nt.fa.gz
gunzip Prorocentrum-minimum-CCMP1329.nt.fa.gz
sed 's/^>/>P_minimumCCMP1329.est./' Prorocentrum-minimum-CCMP1329.nt.fa > ${p_asteroides}/resources/symb_trans/P_minimumCCMP1329.est.fasta

mkdir -p ${p_asteroides}/resources/Prorocentrum/minimum_CCMP2233
cd ${p_asteroides}/resources/Prorocentrum/minimum_CCMP2233
wget ftp://ftp.imicrobe.us/camera/combined_assemblies/Prorocentrum-minimum-CCMP2233.nt.fa.gz
gunzip Prorocentrum-minimum-CCMP2233.nt.fa.gz
sed 's/^>/>P_minimumCCMP2233.est./' Prorocentrum-minimum-CCMP2233.nt.fa > ${p_asteroides}/resources/symb_trans/P_minimumCCMP2233.est.fasta

cd ${p_asteroides}/resources/symb_trans
cat *.fasta > symb_trans.fasta




