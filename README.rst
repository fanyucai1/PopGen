Get started with DRAGEN Germline for PopGen
################################################################################################

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
        --vc-enable-joint-detection=true \

WES
#################
::

    Van Hout, C. V. et al. Exome sequencing and characterization of 49,960 individuals in the UK Biobank. Nature 586, 749–756 (2020).

WGS
#################
::

   Halldorsson B V, Eggertsson H P, Moore K H S, et al. The sequences of 150,119 genomes in the UK Biobank[J]. Nature, 2022, 607(7920): 732-740.

