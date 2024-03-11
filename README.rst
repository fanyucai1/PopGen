Get started with PopGen
++++++++++++++++++++++++++++++
`Illumina DRAGEN Resources for PopGen:https://developer.illumina.com/dragen/dragen-popgen <https://developer.illumina.com/dragen/dragen-popgen>`_

`Accurate and efficient calling of small and large variants from popgen datasets using the DRAGEN Bio-IT Platform <https://sapac.illumina.com/science/genomics-research/articles/popgen-variant-calling-with-dragen.html>`_

`DRAGEN reanalysis of the 1000 Genomes Dataset now available on the Registry of Open Data <https://aws.amazon.com/cn/blogs/industries/dragen-reanalysis-of-the-1000-genomes-dataset-now-available-on-the-registry-of-open-data/>`_

`biobank主页:https://www.ukbiobank.ac.uk <https://www.ukbiobank.ac.uk>`_

.. image:: png/genetic-data-sept2022.jpg

`UKBB Command Line for DRAGEN <https://developer.illumina.com/dragen/dragen-popgen>`_
::

    dragen -r <hg38-ref-dir> \
    --bam-input <input BAM file> \
    --output-directory <out-dir> \
    --output-file-prefix <prefix> \
    --enable-variant-caller=true \
    --vc-emit-ref-confidence=GVCF \
    --vc-enable-vcf-output=true \
    --enable-duplicate-marking=true \
    --enable-map-align=true \
    --enable-map-align-output=true \
    --output-format=CRAM \
    --vc-hard-filter 'DRAGENHardQUAL:all:QUAL<5.0;LowDepth:all:DP<=1' \
    --vc-frd-max-effective-depth=40 \
    --qc-cross-cont-vcf <path-to>/SNP_NCBI_GRCh38.vcf \
    --qc-coverage-region-1 <path-to>/wgs_coverage_regions.hg38_minus_N.interval_list.bed \
    --qc-coverage-reports-1 cov_report \
    --qc-coverage-region-2 <path-to>/acmg59_allofus_19dec2019.GRC38.wGenes.NEW.bed \
    --qc-coverage-reports-2 cov_report \
    --qc-coverage-ignore-overlaps=true \
    --qc-coverage-count-soft-clipped-bases=true \
    --read-trimmers polyg \
    --soft-read-trimmers none \
    --intermediate-results-dir=/ephemeral/staging/tmp/ \
    --repeat-genotype-enable=true \
    --enable-cyp2d6=true \
    --enable-sv=true \
    --enable-cnv=true \
    --cnv-enable-self-normalization=true \
    --vc-enable-joint-detection=true

PopGen data processing and analysis workflows using the DRAGEN Platform (left) and GATK best practices (right) workflows
########################################################################################################################################

.. image:: png/dragen-popgen.png

Global genomic biobanks and studies
########################################################################################################################################
.. image:: png/cohort-studies.png

`Carress H, Lawson D J, Elhaik E. Population genetic considerations for using biobanks as international resources in the pandemic era and beyond[J]. BMC genomics, 2021, 22: 1-19. <https://link.springer.com/article/10.1186/s12864-021-07618-x>`_

Large-scale cohort studies with genomic information
########################################################################################################################################
+--------------------------------+-------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------+----------------+
| Project                        | Description                                                                                           | Website                                                                                     | Country        |
+================================+=======================================================================================================+=============================================================================================+================+
| Human Genome Project (HGP)     | The Initial sequencing program of the human genome                                                    | https://www.genome.gov/human-genome-project                                                 | International  |
+--------------------------------+-------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------+----------------+
| International HapMap Project   | Study of the common pattern of human genetic variation using SNP array                                | https://www.genome.gov/10001688/international-hapmap-project                                | International  |
+--------------------------------+-------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------+----------------+
| 1000 Genomes Project           | Determining the human genetic variation by means of whole-genome sequencing in population scale       | https://www.internationalgenome.org                                                         | International  |
+--------------------------------+-------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------+----------------+
| Human Genome Diversity Project | Biological samples and genetic data collection from different population groups throughout the world  | https://www.hagsc.org/hgdp/                                                                 | International  |
+--------------------------------+-------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------+----------------+
| Simon Genome Diversity Project | Whole-genome sequencing project of diverse human populations                                          | https://docs.cancergenomicscloud.org/v1.0/docs/simons-genome-diversity-project-sgdp-dataset | International  |
+--------------------------------+-------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------+----------------+
| Genome Asia 100k               | WGS-based genome study of people in South and East Asia                                               | https://genomeasia100k.org/                                                                 | International  |
+--------------------------------+-------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------+----------------+
| UK Biobank                     | Biobank study involving 500,000 residents in the UK                                                   | https://www.ukbiobank.ac.uk                                                                 | UK             |
+--------------------------------+-------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------+----------------+
| Genomics England               | WGS-based genome study of patient with rare disease and their families and cancer patients in England | https://www.genomicsengland.co.uk/                                                          | UK             |
+--------------------------------+-------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------+----------------+
| FinnGen                        | Nationwide biobank and genome cohort study in Finland                                                 | https://www.finngen.fi/en                                                                   | Finnland       |
+--------------------------------+-------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------+----------------+
`Tanjo T, Kawai Y, Tokunaga K, et al. Practical guide for managing large-scale human genome data in research[J]. Journal of Human Genetics, 2021, 66(1): 39-52. <https://www.nature.com/articles/s10038-020-00862-1>`_

`WGS-biobank </Biobank/>`_
####################################################################

`WES-biobank </WES/>`_
####################################################################

`Bioinformatics Documents <./bioinformatics/>`_
#####################################################################

`PGx_STR(药物基因组以及短重复序列) <./PGx_STR/>`_
####################################################################

`Medical_genes(医学临床相关基因) <./Medical_genes/>`_
####################################################################

`contamination(样本污染) <./contamination/>`_
####################################################################

`genotyping <./genotyping/>`_
####################################################################

`全基因组关联分析(GWAS)与多基因组风险评估polygenic risk scores (PRS) <./GWAS_PRS/>`_
###################################################################################

`pangenome reference(人泛基因组研究) <./pan-genome/>`_
####################################################################

Imputation(基因型填充)
####################################################################
`GLIMPSE2 is a set of tools for low-coverage whole genome sequencing imputation.  <https://odelaneau.github.io/GLIMPSE/>`_

`Rubinacci S, Hofmeister R J, Sousa da Mota B, et al. Imputation of low-coverage sequencing data from 150,119 UK Biobank genomes[J]. Nature Genetics, 2023, 55(7): 1088-1090. <https://www.nature.com/articles/s41588-023-01438-3>`_

`phasing_genotyping <./phasing_genotyping/>`_
####################################################################

`rare disease and cancer 专病队列研究 <./Genomics_England/>`_
#########################################################################################################################################

`The effect of sequencing coverage on structural variation detection sensitivity测序深度 <./coverage_depth/>`_
#########################################################################################################################################

long-read sequencing for All of Us
####################################################################
`Mahmoud M, Huang Y, Garimella K, et al. Utility of long-read sequencing for All of Us[J]. bioRxiv, 2023: 2023.01. 23.525236. <https://www.biorxiv.org/content/10.1101/2023.01.23.525236v1.abstract>`_

`UK Biobank Allele Frequency Browser <https://afb.ukbiobank.ac.uk/>`_
############################################################################################
