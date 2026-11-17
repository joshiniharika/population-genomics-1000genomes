# Methods

## Data Source
**1000 Genomes Project Phase 3** (Auton et al., 2015, *Nature*)  
- 2,504 individuals from 26 populations across 5 superpopulations (AFR, AMR, EAS, EUR, SAS)  
- Chromosome 22 used as a computationally tractable demonstration dataset  
- Downloaded from: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

## Quality Control
Standard QC filters applied using PLINK 1.9:
- Minor allele frequency (MAF) > 0.05
- Per-SNP missingness < 5%
- Per-individual missingness < 10%
- Hardy-Weinberg equilibrium p > 1×10⁻⁶
- LD pruning: window 50 SNPs, step 10, r² < 0.2

## Principal Component Analysis
- PCA computed using PLINK `--pca` on LD-pruned SNP set
- First 10 PCs retained
- Visualised using R/ggplot2, coloured by superpopulation

## ADMIXTURE
- Ancestry inference performed using ADMIXTURE v1.3 (Alexander et al., 2009)
- K = 2 to 7 tested with 10-fold cross-validation
- Optimal K selected by minimum cross-validation error

## Runs of Homozygosity (ROH)
- ROH detected using PLINK `--homozyg`
- Minimum ROH length: 500 kb
- Minimum SNPs per ROH: 50
- Maximum gap within ROH: 1,000 kb
- FROH calculated as total ROH length / chromosome 22 length (~51 Mb)

### Connection to Ancient DNA Methods
ROH-based inbreeding inference is the methodological foundation for **hapROH** (Ringbauer et al., 2021), which extends this approach to low-coverage ancient DNA. hapROH accounts for the specific noise characteristics of ancient genomes (post-mortem damage, high missingness) while retaining the same underlying principle: long ROH indicate recent shared ancestry between the two haplotypes of an individual, revealing inbreeding or close parental relatedness.

Similarly, **ancIBD** (Ringbauer et al., 2024) extends IBD segment detection to ancient DNA, enabling pedigree reconstruction from prehistoric genomic data — as applied in this repository author's MSc research.

## References
- Auton A et al. (2015). A global reference for human genetic variation. *Nature* 526, 68–74.
- Alexander DH et al. (2009). Fast model-based estimation of ancestry in unrelated individuals. *Genome Research* 19, 1655–1664.
- Purcell S et al. (2007). PLINK: a toolset for whole-genome association and population-based linkage analyses. *AJHG* 81, 559–575.
- Ringbauer H et al. (2021). Insights into human genetic diversity using ancient DNA. *bioRxiv*.
- Ringbauer H et al. (2024). Accurate detection of identity-by-descent segments in human ancient DNA. *Nature Genetics* 56, 143–151.
