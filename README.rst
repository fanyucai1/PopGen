Get started with PopGen
++++++++++++++++++++++++++++++

`biobank主页:https://www.ukbiobank.ac.uk <https://www.ukbiobank.ac.uk>`_

.. image:: genetic-data-sept2022.jpg


`UKBB Command Line for DRAGEN <https://developer.illumina.com/dragen/dragen-popgen>`_
::

    dragen \
    -r <hg38-ref-dir> \
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

关于WES
####################################################################

`Van Hout, C. V. et al. Exome sequencing and characterization of 49,960 individuals in the UK Biobank. Nature 586, 749–756 (2020). <https://www.nature.com/articles/s41586-020-2853-0>`_

`Wang Q, Dhindsa R S, Carss K, et al. Rare variant contribution to human disease in 281,104 UK Biobank exomes[J]. Nature, 2021, 597(7877): 527-532. <https://www.nature.com/articles/s41586-021-03855-y>`_

`Kar S P, Quiros P M, Gu M, et al. Genome-wide analyses of 200,453 individuals yield new insights into the causes and consequences of clonal hematopoiesis[J]. Nature Genetics, 2022, 54(8): 1155-1166. <https://www.nature.com/articles/s41588-022-01121-z>`_

`Protocol for Processing UKB Whole Exome Sequencing Data Sets:https://dnanexus.gitbook.io/uk-biobank-rap/science-corner/whole-exome-sequencing-oqfe-protocol/protocol-for-processing-ukb-whole-exome-sequencing-data-sets <https://dnanexus.gitbook.io/uk-biobank-rap/science-corner/whole-exome-sequencing-oqfe-protocol/protocol-for-processing-ukb-whole-exome-sequencing-data-sets>`_

`Original Quality Functional Equivalent (OQFE) Pipeline - Docker Image <https://hub.docker.com/r/dnanexus/oqfe>`_

`Krasheninina O, Hwang Y C, Bai X, et al. Open-source mapping and variant calling for large-scale NGS data from original base-quality scores[J]. bioRxiv, 2020: 2020.12. 15.356360. <https://www.biorxiv.org/content/10.1101/2020.12.15.356360v1>`_

`2022-uk-biobank-exome-release-faq.pdf <https://github.com/fanyucai1/PopGen/blob/main/2022-uk-biobank-exome-release-faq.pdf>`_

关于WGS
####################################################################

`Halldorsson B V, Eggertsson H P, Moore K H S, et al. The sequences of 150,119 genomes in the UK Biobank[J]. Nature, 2022, 607(7920): 732-740. <https://www.nature.com/articles/s41586-022-04965-x>`_

`VerifyBamID2 <https://github.com/Griffan/VerifyBamID>`_
######################################################################################################

`Zhang F, Flickinger M, Taliun S A G, et al. Ancestry-agnostic estimation of DNA sample contamination from sequence reads[J]. Genome research, 2020, 30(2): 185-194. <https://genome.cshlp.org/content/30/2/185.short>`_

`graphtyper <https://github.com/DecodeGenetics/graphtyper>`_
####################################################################
::

    graphtyper is a graph-based variant caller capable of genotyping population-scale short read data sets.
    It represents a reference genome and known variants of a genomic region using an acyclic graph structure (a "pangenome reference"), which high-throughput sequence reads are re-aligned to for the purpose of discovering and genotyping SNPs, small indels, and structural variants.

`weCall <https://github.com/Genomicsplc/wecall>`_
####################################################################
::

    weCall is a fast, accurate and simple to use command line tool for variant detection in Next Generation Sequencing (NGS) data.
