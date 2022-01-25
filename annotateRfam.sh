#!/usr/bin/env bash

# for f in /home/petarp/scripts/Borgs/bins/fasta/*.fna; do 
# 	g=${f##*/};
# 	cmscan --cpu $SLURM_CPUS_ON_NODE --rfam --cut_ga --nohmmonly --tblout /home/petarp/scripts/Borgs/RfamAnnotations/${g%.fna}.tblout --fmt 2 --clanin /home/petarp/scripts/RPhmm/Rfam/Rfam.clanin /home/petarp/scripts/RPhmm/Rfam/Rfam.cm $f > /home/petarp/scripts/Borgs/RfamAnnotations/${g%.fna}.cmscan
# done

for f in /groups/banfield/ggkbase/exports/petar-public-contigs/*; do 
	g=${f##*/};
	errs=$(zcat $f/$g.contigs.fa.gz | grep -v ">" | grep -n "read_length_" | wc -l)
	if [ $errs -gt 0 ]; then
		continue
	fi
	dbSize=$(esl-seqstat $f/$g.contigs.fa.gz | grep "^Total # residues:" | awk '{ printf ($4*2/1000000); }')
	cmscan -Z $dbSize --cpu $SLURM_CPUS_ON_NODE --noali --tblout ./rRNAOUT/RfamOUT/$g.tblout --fmt 2 --oskip --oclan --cut_ga --rfam --nohmmonly --clanin /home/petarp/scripts/RPhmm/Rfam/Rfam.clanin /home/petarp/scripts/RPhmm/Rfam/Rfam.cm $f/$g.contigs.fa.gz;
done