# 1.VCF操作

[Liao W W, Asri M, Ebler J, et al. A draft human pangenome reference[J]. Nature, 2023, 617(7960): 312-324.](https://www.nature.com/articles/s41586-023-05896-x)

a.multi-sample VCF files were then converted to per-sample VCF files using bcftools view -a -I -s <sample name>
```{.cs}
bcftools view -a -I -s <sample name>
```

b.Left-align and normalize indels, check if REF alleles match the reference, split multiallelic sites into multiple rows; Left-alignment and normalization will only be applied if the --fasta-ref option is supplied.
```{.cs}
bcftools norm -m -any -o out.vcf in.vcf
```

c.The multi-nucleotide polymorphisms and complex indels were further decomposed into SNPs and simple indels,https://github.com/RealTimeGenomics/rtg-tools
```{.cs}
vcfdecompose --break-mnps --no-header --no-gzip --break-indels -i out.vcf -o out1.vcf 
```

d.select snp(-v, --types snps|indels|mnps|other)
```{.cs}
bcftools view -v snps out1.format.vcf >out2.snp.vcf
```

e.converts from \[start, end\] to [start-1, start),BEDOPS:https://bedops.readthedocs.io/en/latest/index.html
```{.cs}
# three custom options for filtering input for each of the three types of variants listed: --snvs, --insertions and --deletions. 
vcf2bed < input.vcf > output.bed
```

# 2.Haplotype Comparison Tools

[llumina/hap.py https://github.com/Illumina/hap.py](https://github.com/Illumina/hap.py)

```
yum install -y devtoolset-8-gcc devtoolset-3-binutils devtoolset-8-gcc-c++ devtoolset-8-gcc-gfortran centos-release-scl boost169-devel autoconf automake glibc-static libstdc++-static python-devel ant cmake bzip2 bzip2-devel ncurses-devel zlib-devel
pip install -y cython numpy scipy biopython matplotlib pandas pysam bx-python pyvcf cyvcf2 nose
wget https://github.com/Illumina/hap.py/archive/refs/tags/v0.3.15.tar.gz
tar xzvf v0.3.15.tar.gz
cd hap.py-0.3.15/
python install.py ~/hap.py-v0.3.15/ --with-rtgtools --no-tests
```

[Platinum Genomes https://github.com/Illumina/PlatinumGenomes/blob/master/files/2017-1.0.files](https://github.com/Illumina/PlatinumGenomes/blob/master/files/2017-1.0.files)

```
hap.py NA12878.vcf.gz test.vcf.gz -o test -r hg19.fa --threads 40 -f ConfidentRegions.bed.gz
```

## 3.参考链接

[GA4GH benchmarking-tools:Germline Small Variant Benchmarking Tools and Standards](https://github.com/ga4gh/benchmarking-tools/tree/master)

[The International Sample Genome Resource (IGSR)](https://github.com/igsr)

[UCSC hg19:giab+platinumGenomes](https://hgdownload.soe.ucsc.edu/gbdb/hg19/)

## 4.BAM 文件数据抽取

```
sambamba-1.0.1-linux-amd64-static view -t 36 -f bam -s 0.1313 -o out.bam in.bam
-s, --subsample=FRACTION
                    subsample reads (read pairs)
```

## 5. SAM文件

[SAM format specification:https://www.samformat.info/sam-format-flag](https://www.samformat.info/sam-format-flag)

![SAM file](./sam.jpg)


## 6.参考文献

[Deng L, Xie B, Wang Y, et al. A protocol for applying a population-specific reference genome assembly to population genetics and medical studies[J]. STAR protocols, 2022, 3(2): 101440.](https://www.sciencedirect.com/science/article/pii/S2666166722003203)

[Gao Y, Yang X, Chen H, et al. A pangenome reference of 36 Chinese populations[J]. Nature, 2023: 1-10.](https://www.nature.com/articles/s41586-023-06173-7)

[Liao W W, Asri M, Ebler J, et al. A draft human pangenome reference[J]. Nature, 2023, 617(7960): 312-324.](https://www.nature.com/articles/s41586-023-05896-x)
