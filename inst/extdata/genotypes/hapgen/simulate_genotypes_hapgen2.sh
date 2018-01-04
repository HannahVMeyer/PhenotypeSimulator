cd ~/PhenotypeSimulator/inst/extdata/genotypes/hapgen

# directory with data downloaded from https://mathgen.stats.ox.ac.uk/impute/impute_v1.html#Using_IMPUTE_with_the_HapMap_Data
hapdir=~/data/hmeyer/Supplementary/CEU.0908.impute.files

hapgen2 -m $hapdir/genetic_map_chr22_combined_b36.txt \
        -l legend_file \
        -h $hapdir/CEU.0908.chr22.hap \
        -o genotypes_hapgen \
        -n 100 0 \
        -dl 45162 0 0 0 \
        -no_haps_output
