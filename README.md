###Walkthrough to determine parent-of-origin effects on gene expression (allele bias)

This analysis is tailored for pooled sequencing experiments in which DNA-seq is available for each parental population and RNA-Seq for the hybrid offspring. Parental sequence processing is only shown for one parent, but must be done separately. Not shown: standard preprocessing for sequencing reads such as error correction and quality filtering. 

Prerequisites:  
[bowtie1](https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/)
[MMSEQ](https://github.com/eturro/mmseq)
[GATK](https://software.broadinstitute.org/gatk/download/)
[FreeBayes](https://github.com/ekg/freebayes)
[VCFfilter](https://github.com/vcflib/vcflib)
[bedtools](https://github.com/arq5x/bedtools2/releases)
[seqmagick](https://github.com/fhcrc/seqmagick)
[GNU Parallel](https://www.gnu.org/software/parallel/)

####Add read groups to DNA-seq parental samples

Brief overview given [here](https://www.biostars.org/p/47487/). Here are the commands I used for four samples from one parental population (two technical replicates for each sex). 

`java -jar /home/hdd/3/picard/dist/picard.jar AddOrReplaceReadGroups I=C-ma.bam O=C-ma_RG.bam RGID=1 RGLB=1 RGPL=illumina RGPU=1 RGSM=C-F`
`java -jar /home/hdd/3/picard/dist/picard.jar AddOrReplaceReadGroups I=C-mb.bam O=C-mb_RG.bam RGID=2 RGLB=1 RGPL=illumina RGPU=2 RGSM=C-F`
`java -jar /home/hdd/3/picard/dist/picard.jar AddOrReplaceReadGroups I=C-fa.bam O=C-fa_RG.bam RGID=3 RGLB=2 RGPL=illumina RGPU=3 RGSM=C-M`
`java -jar /home/hdd/3/picard/dist/picard.jar AddOrReplaceReadGroups I=C-fb.bam O=C-fb_RG.bam RGID=4 RGLB=2 RGPL=illumina RGPU=4 RGSM=C-M`

####Mark and remove PCR duplicates

See [here](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates) for usage. 

`java -jar picard.jar MarkDuplicates I=C-fa_RG.bam O=C-fa_RG_DUO.bam REMOVE_DUPLICATES=TRUE`

####Find fixed SNPs in each parental population with FreeBayes

Note that I am only calling SNPs for annotated gene regions (transcript_all.bed). The parameters are set a little loosely to minimize false negatives. Stringent filtering will be applied later to avoid false positives. Since heterozygosity is impossible to determine in pooled samples, only those SNPs segregating near 100% frequency are kept. 

`freebayes -f tcas331.fa -N -n 2 --ploidy 100 -t transcript_all.bed --pooled-discrete --min-alternate-fraction 0.95 --strict-vcf --min-coverage 50 -b C-fa_RG_DUO.bam -b C-fb_RG_DUO.bam -b C-ma_RG_DUO.bam -b C-mb_RG_DUO.bam > csm_all_transcripts.vcf`

####Apply stringent filtering to VCF output

`vcffilter -f "DP < 140" -f "QUAL > 20" -f "AF > 0.99" csm_all_transcripts.vcf > csm.vcf`

* DP -- filter hits with extreme read depth (140 is two std dev. above average sequencing coverage)
* AF -- allele frequency of 100% (fixed in population)
* QUAL -- threshold for Freebayes' variant quality score

####Make alternate reference genome for each parental population

`java -jar GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R tcas331.fa -o csm.fa -V csm.vcf`

Remove additional header information added by GATK (this command edits in place)

`sed -i 's/>.* />/g;s/:.*//g' csm.fa`

####Extract gene regions from new reference, add parental identifier, and concatenate

Extract gene regions
`bedtools getfasta -fi csm.fa -fo csm_transcript.fa -bed transcript_all.bed -s -name`

Here I am adding '_csm' as my maternal parent identifier
`sed 's/>.*/&_csm/' csm_transcript.fa > csm_header.fa`

All previous steps must be done separately for each parent. Now we concatenate the files
`cat csm_header.fa des_header.fa > hybrid_vigor.fa`

Add header that MMSEQ requires
`perl -pe 's/\>(\w+)/\>$1 gene:\1 transcript:\1/' hybrid_vigor.fa`

Deduplicate as safety measure. My transcript reference contained duplicates for some reason, and this will throw off read mapping big time!

`seqmagick --deduplicate-sequences hybrid_vigor.fa`

####Map hybrid offspring RNA-seq reads to new merged transcript reference

Now that we have a bi-allelic reference transcriptome, we can determine biased expression if F1 hybrid RNA-seq reads uniquely map to one parent's copy in a dominant manner. 

Quick note: some raw read error-correction tools (such as [BFC](https://github.com/lh3/bfc) add tags to the FASTQ header which will throw an error with bowtie such as "[main_samview] truncated file." This can be alleviated by 

`awk '{print $1}' sample.1.fq > sample_.1.fq`

Read mapping for each F1 hybrid sample (it is important you use bowtie1 per MMSEQ author) :

`bowtie --wrapper basic-0 -a --best --strata -S -m 100 -X 500 --norc --fullref --chunkmbs 256 -p 24 hybrid_vigor -1 ../mapping/bfc_1p1n.1_val_1.fq -2 ../mapping/bfc_1p1n.2_val_2.fq | samtools view -F 0xC -bS - | samtools sort -n - 1p1n`

####Count hits from each alignment file produced in previous steps

`for i in *.bam; do bam2hits -m "(\S+).*gene:(\S+).*" 1 2 hybrid_vigor.fa $i > ${i/.bam/.hits}; done`

`find . -name "*.hits" -exec mmseq {} {} \;`

Here you must replace "csm" and "des with your own unique parental identifiers added to the reference transcriptome:

`for i in *.gene.mmseq; do grep 'csm' $i > $i.csm && grep 'des' $i > $i.des; done`

Since we've separated the parents from the output, we must add the header back to each file:

`find . -name '*.csm' -exec sed -i '1 i\feature_id\tlog_mu\tsd\tmcse\tiact\teffective_length\ttrue_length\tunique_hits\tmean_proportion\tmean_probit_proportion\tsd_probit_proportion\tlog_mu_em\tobserved\tntranscripts\tpercentiles5,25,50,75,95\tpercentiles_proportion5,25,50,75,95' {} \;`

`find . -name '*.des' -exec sed -i '1 i\feature_id\tlog_mu\tsd\tmcse\tiact\teffective_length\ttrue_length\tunique_hits\tmean_proportion\tmean_probit_proportion\tsd_probit_proportion\tlog_mu_em\tobserved\tntranscripts\tpercentiles5,25,50,75,95\tpercentiles_proportion5,25,50,75,95' {} \;`

`find . -name '*mmseq.csm' -exec sed -i 's/_csm//g' {} \; ;  find . -name '*mmseq.des' -exec sed -i 's/_des//g' {} \;`

####Set up the experimental contrasts between F1 hybrid samples

Here you must perform a contrast between parental read counts for each experimental condition. Refer to [MMSEQ](https://github.com/eturro/mmseq) website for additional design matrices. 

`mmdiff -de 2 2 1p1n.hits.gene.mmseq.csm 3p1n.hits.gene.mmseq.csm 1p1n.hits.gene.mmseq.des 3p1n.hits.gene.mmseq.des > NORM_P1_csm_VS_des.mmdiff`

When interpreting the output, here is the author's recommendation

`Regarding the threshold for mmdiff, I suggest looking at the distribution of the "posterior_probability" values. You should see a small number of features that rise up well above twice the prior (i.e. > 0.2, by default). That should give reasonable results. In theory the FDR is the mean of 1-posterior_probability of the selected transcripts above a particular threshold, but that assumes that the posterior probabilities are well calibrated, so I wouldn't rely on this notion too much.`
