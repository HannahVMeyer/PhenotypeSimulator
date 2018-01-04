cd ~/PhenotypeSimulator/inst/extdata/genotypes/plink

echo -e "1000\tSNP\t0.00\t1.00\t1.00\t1.00" > plink_sim_par.txt

plink --simulate plink_sim_par.txt \
      --simulate-ncases 0 \
      --simulate-ncontrols 100 \
      --out genotypes_plink
