for f in RP_HMMs/2/*; do g=${f##*/}; hmmsearch --tblout ./hmmOUT/Bacteria/${g%.hmm}_RifSed_csp1_19ft_2_Methanoperedens_44_10.out --cpu 12 $f ./proteinSeqs/RifSed_csp1_19ft_2_Methanoperedens_44_10.proteins.faa ;done

for f in RP_HMMs/2157/*; do g=${f##*/}; hmmsearch --tblout ./hmmOUT/Archaea/${g%.hmm}_RifSed_csp1_19ft_2_Methanoperedens_44_10.out --cpu 12 $f ./proteinSeqs/RifSed_csp1_19ft_2_Methanoperedens_44_10.proteins.faa ;done

for f in RP_HMMs/2759/e*.hmm; do g=${f##*/}; hmmsearch --tblout ./hmmOUT/Eukarya/${g%.hmm}_RifSed_csp1_19ft_2_Methanoperedens_44_10.out --cpu 12 $f ./proteinSeqs/RifSed_csp1_19ft_2_Methanoperedens_44_10.proteins.faa ;done
