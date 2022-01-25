#!/usr/bin/env bash
#for f in RP_HMMs/2/*; do g=${f##*/}; hmmsearch --tblout ./hmmOUT/Bacteria/${g%.hmm}_RifSed_csp1_19ft_2_Methanoperedens_44_10.out --cpu 12 $f ./proteinSeqs/RifSed_csp1_19ft_2_Methanoperedens_44_10.proteins.faa ;done

# for protSeq in ./BorgProtSeqs/*.faa; do
# 	outSeqName=${protSeq##*/};
# 	outDir=./BorgProtSeqs/hmmOUT/${outSeqName%.proteins.faa}
# 	mkdir $outDir
# 	for f in RP_HMMs/2157/*; do 
# 		g=${f##*/}; 
# 		hmmsearch --tblout $outDir/${g%.hmm}.out -T 50 --cpu 12 $f $protSeq ; 
# 	done;
# done

#for f in RP_HMMs/2759/e*.hmm; do g=${f##*/}; hmmsearch --tblout ./hmmOUT/Eukarya/${g%.hmm}_RifSed_csp1_19ft_2_Methanoperedens_44_10.out --cpu 12 $f ./proteinSeqs/RifSed_csp1_19ft_2_Methanoperedens_44_10.proteins.faa ;done

for f in /groups/banfield/ggkbase/exports/petar-public-contigs/*; do 
	for hmm in rRNA/*.hmm; do
		outFolder=${hmm##*/};
		g=${f##*/};
		hmmsearch --noali --pfamtblout ./rRNAOUT/${outFolder%.hmm}/$g.out --cpu $SLURM_CPUS_ON_NODE $hmm $f/$g.contigs.fa.gz;
	done;
done