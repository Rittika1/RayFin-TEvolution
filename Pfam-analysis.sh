# ##--To install hmmer before this, use 
# sudo apt install hmmer
# 
# gunzip uniprot_sprot.pep.gz
# gunzip Pfam-A.hmm.gz
# makeblastdb -in uniprot_sprot.pep -dbtype prot
# hmmpress Pfam-A.hmm

blastx -query trinity_out_dir.Trinity.fasta \
  -db uniprot_sprot.pep \
  -num_threads 8 \
  -max_target_seqs 1 \
  -outfmt 6 > gar_blastx.outfmt6
  
blastp -query Lep_oss_0048.pep \
  -db uniprot_sprot.pep \
  -num_threads 8 \
  -max_target_seqs 1 \
  -outfmt 6 > gar_blastp.outfmt6
  
hmmscan --cpu 8 \
  --domtblout TrinotatePFAM.out \
  Pfam-A.hmm Lep_oss_0048.pep > pfam.log

signalp -format short -batch 2000000 -fasta Lep_oss_0048.pep -stdout

tmhmm --short < Lep_oss_0048.pep > tmhmm.out

RnammerTranscriptome.pl --transcriptome ttrinity.fasta \
  --path_to_rnammer /usr/bin/software/rnammer_v1.2/rnammer

wget "https://data.broadinstitute.org/Trinity/Trinotate_v3_RESOURCES/Trinotate_v3.sqlite.gz" -O Trinotate.sqlite.gz

get_Trinity_gene_to_trans_map.pl trinity.fasta >  Trinity.fasta.gene_trans_map


##--- before this you need to install a bunch of stuffs, like the Trinotate and the perl dependencies
#
TRINOTATE_HOME/admin/Build_Trinotate_Boilerplate_SQLite_db.pl  Trinotate
# 
Trinotate-Trinotate-v3.2.2/Trinotate Trinotate.sqlite init --gene_trans_map Trinity.fasta.gene_trans_map --transcript_fasta trinity_out_dir.Trinity.fasta --transdecoder_pep Lep_oss_0048.pep

Trinotate-Trinotate-v3.2.2/Trinotate Trinotate.sqlite LOAD_swissprot_blastp gar_blastp.outfmt6

Trinotate-Trinotate-v3.2.2/Trinotate Trinotate.sqlite LOAD_swissprot_blastx gar_blastx.outfmt6

Trinotate-Trinotate-v3.2.2/Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
Trinotate-Trinotate-v3.2.2/Trinotate Trinotate.sqlite LOAD_tmhmm gar_tmhmm.out
Trinotate-Trinotate-v3.2.2/Trinotate Trinotate.sqlite LOAD_signalp signalp_gar.out


Trinotate-Trinotate-v3.2.2/Trinotate Trinotate.sqlite report > trinotate_annotation_report.xls


##--Go go_annotations
Trinotate-Trinotate-v3.2.2/util/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls trinotate_annotation_report.xls -G --include_ancestral_terms  > go_annotations.txt

###---subsetdef: goslim_generic "Generic GO slim"
