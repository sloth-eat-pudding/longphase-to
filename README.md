# LongPhase
LongPhase is an ultra-fast program for simultaneously co-phasing SNPs, small indels, large SVs, and modifications for Nanopore and PacBio platforms. It can produce nearly chromosome-scale haplotype blocks by using Nanpore ultra-long reads without the need for additional trios, chromosome conformation, and strand-seq data. On an 8-core machine, LongPhase can finish phasing a human genome in 10-20 minutes.

---
- [Installation](#installation)
- [Usage](#usage)
	- [Phase commad](#phase-command)
		- [SNP-only phasing](#snp-only-phasing)
		- [SNP and indel co-phasing](#snp-and-indel-co-phasing)
		- [SNP and SV co-phasing](#snp-and-sv-co-phasing)
		- [SNP and modification co-phasing](#snp-and-modification-co-phasing)
		- [The complete list of phase parameters](#the-complete-list-of-phase-parameters)
		- [Output files of SNP-only phasing](#output-files-of-snp-only-phasing)
		- [Output files of SNP and SV co-phasing](#output-files-of-snp-and-sv-co-phasing)
	- [Haplotag command](#haplotag-command)
		- [The complete list of haplotag parameters](#the-complete-list-of-haplotag-parameters)
- [Input Preparation](#input-preparation)
	- [Generate reference index](#generate-reference-index)
	- [Generate alignment and index files](#generate-alignment-and-index-files)
	- [Generate single nucleotide polymorphism (SNP) file](#generate-single-nucleotide-polymorphism-snp-file)
	- [Generate Structural variation (SV) file](#generate-structural-variation-sv-file)
 	- [Modified bam process](#modified-bam-process)
- [Comparison with other SNP-phasing programs](#comparison-with-other-snp-phasing-programs)
- [Citation](#citation)
- [Contact](#contact)
---
## Installation
You are recommended to download a [linux 64bit binary release](https://github.com/twolinin/longphase/releases/download/v1.5/longphase_linux-x64.tar.xz) without compilation. 

```
wget https://github.com/twolinin/longphase/releases/download/v1.5/longphase_linux-x64.tar.xz
tar -xJf longphase_linux-x64.tar.xz
```

An executable file, longphase_linux-x64, can be executed directly. If you need to compile a local version, you can clone and compile using the following commands. 

```
git clone https://github.com/twolinin/longphase.git
cd longphase
autoreconf -i
./configure
make -j 4
```

---
## Usage
### Phase command
#### SNP-only phasing
For SNP-only phasing, the input of LongPhase consists of SNPs in VCF (e.g., SNP.vcf), an indexed reference in Fasta (e.g., reference.fasta, reference.fasta.fai), and one (or multiple) indexed read-to-reference alignment in BAM (e.g., alignment1.bam, alignment1.bai, alignment2.bam, ...) (see [Input Preparation](#input-preparation)). The users should specify the sequencing platform (--ont for Nanopore and --pb for PacBio). An example of SNP phasing usage is shown below.
```
longphase phase \
-s SNP.vcf \
-b alignment1.bam \
-b alignment2.bam \
-r reference.fasta \
-t 8 \
-o phased_prefix \
--ont # or --pb for PacBio Hifi
```
#### SNP and indel co-phasing
When the SNP.vcf file contains both SNP and small insertion/deletion variations (indels), use the `--indels` parameter for co-phasing SNPs and indels.
```
longphase phase \
-s SNP.vcf \
-b alignment.bam \
-r reference.fasta \
-t 8 \
-o phased_prefix \
--indels \
--ont # or --pb for PacBio Hifi
```
#### SNP and SV co-phasing
When co-phasing SNPs and SVs, except for the same input (i.e., SNPs, reference, and alignments), LongPhase takes an extra input of called SVs in VCF (e.g., [SV_sniffles.vcf](http://bioinfo.cs.ccu.edu.tw/bioinfo/HG002_60x/)), which should be generated by Sniffles (with --num_reads_report in sniffles1 and --output-rnames in sniffles2) or CuteSV (with --report_readid--genotype)(see [Input Preparation](#input-preparation)).
```
longphase phase \
-s SNP.vcf \
--sv-file SV.vcf \
-b alignment.bam \
-r reference.fasta \
-t 8 \
-o phased_prefix \
--ont # or --pb for PacBio Hifi
```

#### SNP and modification co-phasing
When co-phasing SNPs and modifications (5mC supported at this moment), it is necessary to first call the modified loci (in VCF) using [modcall](#modcall-command) and then activate SNP-modification co-phasing via specifying the `--mod-file` with the modcall-generated VCF. The phased SNPs and modifications will be written into new VCFs accordingly.
```
longphase phase \
-s SNP.vcf \
--mod-file modcall.vcf \
-b alignment.bam \
-r reference.fasta \
-t 8 \
-o phased_prefix \
--ont # or --pb for PacBio Hifi
```

#### The complete list of phase parameters
```
Usage:  phase [OPTION] ... READSFILE
   --help                                 display this help and exit.

require arguments:
   -s, --snp-file=NAME                    input SNP vcf file.
   -b, --bam-file=NAME                    input bam file.
   -r, --reference=NAME                   reference fasta.
   --ont, --pb                            ont: Oxford Nanopore genomic reads.
                                          pb: PacBio HiFi/CCS genomic reads.

optional arguments:
   --sv-file=NAME                         input SV vcf file.
   --mod-file=NAME                        input modified vcf file.(produce by longphase modcall)
   -t, --threads=Num                      number of thread. default:1
   -o, --out-prefix=NAME                  prefix of phasing result. default: result
   --indels                               phase small indel. default: False
   --dot                                  each contig/chromosome will generate dot file.

parse alignment arguments:
   -q, --mappingQuality=Num               filter alignment if mapping quality is lower than threshold. default:1

phasing graph arguments:
   -a, --connectAdjacent=Num              connect adjacent N SNPs. default:6
   -d, --distance=Num                     phasing two variant if distance less than threshold. default:300000
   -1, --readsThreshold=[0~1]             give up SNP-SNP phasing pair if the number of reads of the
                                          two combinations are similar. default:0.05
   -v, --confidentHaplotype=[0~1]         the haplotype of the current SNP is judged by the haplotype of the previous N SNPs.
                                          if the threshold is higher, the consistency of SNP needs to be higher. default:0.5
   -j, --judgeInconsistent=[0~1]          the proportion of inconsistent haplotypes among the haplotypes of the previous N SNPs.
                                          inconsistent SNPs are tagged if the proportion is below the threshold. default:0.4
   -i, --inconsistentThreshold=Num        phased genotype correction is performed when a SNP is tagged multiple times. default:5

haplotag read correction arguments:
   -m, --readConfidence=[0.5~1]           The confidence of a read being assigned to any haplotype. default:0.65
   -n, --snpConfidence=[0.5~1]            The confidence of assigning two alleles of a SNP to different haplotypes. default:0.75

```
---
### Output files of SNP-only phasing
When phasing SNPs alone, LongPhase outputs the results into a VCF file. The alleles of the two haplotypes are stored in the GT field (e.g., 1|0), whereas the left and right alleles of the vertical bar represents the paternal or maternal haplotypes. The last PS field (e.g., 16809) represents the identifier of the block. For instance, the following example illustrates two haplotypes of five phased SNPs, CCCCC and GATGT, in the same block 16809.

```
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  default
1       16809   .       C       G       8.4     PASS    .       GT:GQ:DP:AD:VAF:PL:PS   1|0:8:79:51,28:0.35443:7,0,16:16809
1       16949   .       A       C       8.4     PASS    .       GT:GQ:DP:AD:VAF:PL:PS   0|1:8:67:43,21:0.313433:7,0,25:16809
1       21580   .       C       T       13.9    PASS    .       GT:GQ:DP:AD:VAF:PL:PS   1|0:14:75:50,24:0.32:13,0,30:16809
1       23359   .       C       G       13.2    PASS    .       GT:GQ:DP:AD:VAF:PL:PS   1|0:13:52:24,18:0.346154:13,0,42:16809
1       24132   .       C       T       11.1    PASS    .       GT:GQ:DP:AD:VAF:PL:PS   1|0:11:63:41,17:0.269841:10,0,29:16809
```

---
### Output files of SNP and SV co-phasing
When co-phasing SNPs and SVs, two VCFs (one for SNPs and one for SVs) are outputted. Similarly, the phased SVs are stored in the GT field and the block ID is in the PS field. For instance, the following example illustrates two haplotypes of five SNPs and two SVs, A\<INS\>G\<noSV\>TCC and G\<noSV\>A\<INS\>ATT, which are co-phased in the same block 382189.

An example of SNP VCF file
```
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  default
1       465289  .       A       G       34.2    PASS    .       GT:GQ:DP:AD:VAF:PL:PS   1|0:4:16:2,11:0.6875:31,0,1:382189
1       544890  .       G       A       3.1     PASS    .       GT:GQ:DP:AD:VAF:PL:PS   1|0:3:31:30,0:0:0,0,31:382189
1       545612  .       T       A       6.8     PASS    .       GT:GQ:DP:AD:VAF:PL:PS   1|0:7:37:28,8:0.216216:5,0,32:382189
1       545653  .       C       T       14      PASS    .       GT:GQ:DP:AD:VAF:PL:PS   1|0:14:44:30,14:0.318182:13,0,26:382189
1       561458  .       C       T       5.1     PASS    .       GT:GQ:DP:AD:VAF:PL:PS   1|0:5:17:14,0:0:3,0,17:382189
```
An example of SV VCF file
```
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  default
1       534057  7       N       AAATTCGCGCATATCACGGGTGCCGCCTCTGTGCAGCTCACGAAACGCCATACTACGGTGCTGCTCAGCAGCTACGGAATCGCTATACCTACGCGAGCTGCCTCAGCAGCCAC       .       PASS    IMPRECISE;SVMETHOD=Snifflesv1.0.11;CHR2=1;END=534128;STD_quant_start=25.369041;STD_quant_stop=29.191054;Kurtosis_quant_start=-0.239012;Kurtosis_quant_stop=-1.333733;SVTYPE=INS;RNAMES=0120d560-50f0-4298-8b03-7bd30f3cf139,030ac5d4-e616-4ce9-8ad3-243835335085,0cf1b0d9-2b4d-463d-a658-01b4b040dc63,1c9e982a-8af7-4ba0-8cc2-154a679c72e2,22e11f79-0067-4735-8b69-97d951ca702f,2ca8a6f4-be9d-4df5-80d2-dc1743f97a84,35dff960-22b6-4216-af69-8878b8860362,390d2fb4-9224-41a1-a9fe-6cb3bbe4273a,3e333422-12ca-4f16-afb8-ed7611dcbc2c,3e8ed78a-b857-4941-bbc1-52ca51e26c08,4191371c-49ea-466d-aadc-06f27cdf1050,4aaae789-54fe-4fa5-84b3-5524dc2b3796,581e0cfb-2491-44d7-a2e1-ba1516ba0f2f,59749531-9abf-4ff4-a4a1-31484ba3d32d,5c97b0a9-925e-4153-952d-0f437171d3dc,6067590e-956c-442f-bbb7-cae597d616ad,623804bb-e2fe-415d-96ae-3d06aec63e5d,672244ce-2d5d-45cf-beb2-ddeddae917e8,6b79aa23-7c9c-49dc-9b88-8419c88c7a36,6e60d235-6654-4ef8-9feb-70f12a397721,6fbb55c5-57fc-43bc-8a24-b0058778054c,8e10bf13-9674-489c-924e-182a42e08a34,aa6ba092-4221-4d54-8819-811448c34983,af2169b3-b308-4db5-9675-15ff5f68d8dd,b214fcbd-77de-4dd5-84db-6d2b7e1f3158,c140eaba-e0e7-44e7-9f16-c8c67fd4a2f2,c7835cf7-44c0-44da-b10e-b2468fc8caab,ca4aa84d-34d1-4639-8634-b6a5540129ca,caba4bde-cdc5-4344-9803-a3c158525b0c,e0747feb-60bb-40db-a144-a9b43dd13256,e6992c7d-c00e-40e7-b80b-562094a9b60f,e8bb376c-20e0-4bed-a61f-b82b5c37ef6f,f3242a61-deec-49e7-b99f-335a1ba13791,f87dfdf7-7b68-421b-b395-3769a5fa3ac1,f91a7627-7fdb-4f03-8f33-0ed1649d96fe;SUPTYPE=AL;SVLEN=43;STRANDS=+-;RE=35;REF_strand=44;AF=0.795455    GT:DR:DV:PS     0|1:9:35:382189
1       545892  8       N       ACACGCGGGCCGTGGCCAGCAGGCGGCGCTGCAGGAGAGGAGATGCCCAGGCCTGGCGGCC   .       PASS    IMPRECISE;SVMETHOD=Snifflesv1.0.11;CHR2=1;END=545893;STD_quant_start=28.919840;STD_quant_stop=28.543200;Kurtosis_quant_start=-0.382251;Kurtosis_quant_stop=-0.130808;SVTYPE=INS;RNAMES=0120d560-50f0-4298-8b03-7bd30f3cf139,030ac5d4-e616-4ce9-8ad3-243835335085,0cf1b0d9-2b4d-463d-a658-01b4b040dc63,22e11f79-0067-4735-8b69-97d951ca702f,2ca8a6f4-be9d-4df5-80d2-dc1743f97a84,3977c988-9901-4e5b-9f9c-b8ebfcce8e93,3e333422-12ca-4f16-afb8-ed7611dcbc2c,4191371c-49ea-466d-aadc-06f27cdf1050,4aaae789-54fe-4fa5-84b3-5524dc2b3796,5933e1b7-1aeb-4437-a875-3befbf703420,623804bb-e2fe-415d-96ae-3d06aec63e5d,672244ce-2d5d-45cf-beb2-ddeddae917e8,6b79aa23-7c9c-49dc-9b88-8419c88c7a36,7842d9f1-9a77-4c9a-ab5b-5a644ed2d355,7ba26d64-d9b0-475f-8d5f-1fa73fc42d93,8e10bf13-9674-489c-924e-182a42e08a34,a2b1b2ef-1e28-465e-8b3f-c44e15990d8b,a45514f1-4aae-40eb-94eb-2969722a7b05,b8181546-6839-49cd-b64f-b65c96369a2b,c140eaba-e0e7-44e7-9f16-c8c67fd4a2f2,c7835cf7-44c0-44da-b10e-b2468fc8caab,ca4aa84d-34d1-4639-8634-b6a5540129ca,d56f0abe-4389-4197-a151-0eb567fb99f0,e6992c7d-c00e-40e7-b80b-562094a9b60f,e8bb376c-20e0-4bed-a61f-b82b5c37ef6f,ec325153-0c55-4ece-8f3c-c432701e6750,f3242a61-deec-49e7-b99f-335a1ba13791,f91a7627-7fdb-4f03-8f33-0ed1649d96fe;SUPTYPE=AL;SVLEN=62;STRANDS=+-;RE=28;REF_strand=51;AF=0.54902        GT:DR:DV:PS     1|0:23:28:382189
```

### Haplotag command
This command tags (assigns) each read (in BAM) to one haplotype in the phased SNP/SV VCF. i.e., reads will be tagged as HP:i:1 or HP:i:2. In addition, the haplotype block of each read is stored in the PS tag. The phased VCF can be also generated by other programs as long as the PS or HP tags are encoded. The author can specify ```--log``` for additionally output a plain-text file containing haplotype tags of each read without parsing BAM.
```
longphase haplotag \
-s phased_snp.vcf \
--sv-file phased_sv.vcf \
-b alignment.bam \
-t 8 \
-o tagged_bam_prefix
```

#### The complete list of haplotag parameters

```
Usage:  haplotag [OPTION] ... READSFILE
      --help                          display this help and exit.

require arguments:
      -s, --snp-file=NAME             input SNP vcf file.
      -b, --bam-file=NAME             input bam file.
      -r, --reference=NAME            reference fasta.
optional arguments:
      --tagSupplementary              tag supplementary alignment. default:false
      --sv-file=NAME                  input phased SV vcf file.
      -q, --qualityThreshold=Num      not tag alignment if the mapping quality less than threshold. default:0
      -p, --percentageThreshold=Num   the alignment will be tagged according to the haplotype corresponding to most alleles.
                                      if the alignment has no obvious corresponding haplotype, it will not be tagged. default:0.6
      -t, --threads=Num               number of thread. default:1
      -o, --out-prefix=NAME           prefix of phasing result. default:result
      --log                           an additional log file records the result of each read. default:false
```

### ModCall command
Modcall implements a module for calling allele-specific modifications in latest Nanopore and PacBio basecalled reads, assuming the modifications are stored in the raw BAM files using `MM` and `ML` tags. Unaligned BAM files can be converted into alignment BAM files through [preprocessing](#modified-bam-process). Modcall identifies allele-specific modifications from the BAM files and store them in a VCF file. An example of this command is shown below.
```
longphase modcall \
-b alignment.bam \
-r reference.fasta \
-t 8 \
-o modcall
```

An example of ModCall VCF file
```
##INFO=<ID=RS,Number=.,Type=String,Description="Read Strand">
##INFO=<ID=MR,Number=.,Type=String,Description="Read Name of Modified position">
##INFO=<ID=NR,Number=.,Type=String,Description="Read Name of nonModified position">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=MD,Number=1,Type=Integer,Description="Modified Depth">
##FORMAT=<ID=UD,Number=1,Type=Integer,Description="Unmodified Depth">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth"
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE
chr1    11200   .       C       .       .       .       RS=P;MR=eb459876-8c81-4714-a496-a90ea8be94d2,6ca3a71f-62fd-416e-8c6e-8c4a9c054e1a,0d8b7c68-d98d-4045-a572-82fedac62da5,4ee6d422-e677-41fc-9ef3-8b387c05973a,71b5dfe9-7cd1-4959-b634-b9d162468edb,6e205f04-560a-4975-932d-dbd60ead695d,ae2238f8-b622-4cb5-8f02-4c3f54ab8ca3,d7ef87d0-bffe-404d-9f10-62eceb0c5121,c0e0c526-e193-4ee2-81d6-1fbe0c970dc1;NR=8249c7c7-fd04-4a4c-985d-2bbcb2030bc4,5db6fcd3-5780-494b-bfdd-4d9d7282d012,e3b03dc0-8399-4e1e-bd0d-e5e4e3e2911a,b8dc7af8-9f78-4dac-b8cc-35e157c51621,0b2638c1-8380-48b4-b08c-a6161495ad9d,7b1e5c16-f0a6-47f6-9726-247c762a10ca,ac7ae685-082e-413e-b5a6-b4c12d49a1c2;     GT:MD:UD:DP     0/1:9:7:16
chr1    11201   .       G       .       .       .       RS=N;MR=7f505676-704b-41b7-a292-666559b40159,7b070045-225e-49ed-9fff-e7cd1bffdfd5,e7556960-d21d-4859-931c-41c6bc196141,69c925f2-ae2c-4cf1-a8c6-6b39868f0d86,7e6a4472-5f8d-4d20-a70b-ec9db403121b,5f879a79-5b28-4a58-b12a-8445c44a37aa;NR=123518e1-4dbe-4a72-9b31-3450b6b8d3da,1783e61c-546a-4414-a8ef-92297ee1dc1f,c9772456-bb97-41c8-b561-001cbcc1233f,f40a87f4-0240-4137-b3ed-b98ba5a479cb,3c6d0a5e-23a5-4a25-9b3b-ed90b19188a0,d584347d-6db5-42bf-8cd8-95ac8a2c51d1; GT:MD:UD:DP     0/1:6:6:12
chr1    11434   .       C       .       .       .       RS=P;MR=eb459876-8c81-4714-a496-a90ea8be94d2,6ca3a71f-62fd-416e-8c6e-8c4a9c054e1a,0d8b7c68-d98d-4045-a572-82fedac62da5,71b5dfe9-7cd1-4959-b634-b9d162468edb,5db6fcd3-5780-494b-bfdd-4d9d7282d012,6e205f04-560a-4975-932d-dbd60ead695d,ae2238f8-b622-4cb5-8f02-4c3f54ab8ca3,d7ef87d0-bffe-404d-9f10-62eceb0c5121;NR=8249c7c7-fd04-4a4c-985d-2bbcb2030bc4,e3b03dc0-8399-4e1e-bd0d-e5e4e3e2911a,b8dc7af8-9f78-4dac-b8cc-35e157c51621,0b2638c1-8380-48b4-b08c-a6161495ad9d,7b1e5c16-f0a6-47f6-9726-247c762a10ca,ac7ae685-082e-413e-b5a6-b4c12d49a1c2,c0e0c526-e193-4ee2-81d6-1fbe0c970dc1;  GT:MD:UD:DP     0/1:8:7:16
chr1    11435   .       G       .       .       .       RS=N;MR=7f505676-704b-41b7-a292-666559b40159,7b070045-225e-49ed-9fff-e7cd1bffdfd5,e7556960-d21d-4859-931c-41c6bc196141,7e6a4472-5f8d-4d20-a70b-ec9db403121b,5f879a79-5b28-4a58-b12a-8445c44a37aa;NR=123518e1-4dbe-4a72-9b31-3450b6b8d3da,1783e61c-546a-4414-a8ef-92297ee1dc1f,c9772456-bb97-41c8-b561-001cbcc1233f,69c925f2-ae2c-4cf1-a8c6-6b39868f0d86,f40a87f4-0240-4137-b3ed-b98ba5a479cb,d584347d-6db5-42bf-8cd8-95ac8a2c51d1;      GT:MD:UD:DP     0/1:5:6:12
```

#### The complete list of modcall parameters

```
Usage:  modcall [OPTION] ... READSFILE
      --help                        display this help and exit.
require arguments:
      -b, --bam-file=NAME           modified sorted bam file.
      -r, --reference=NAME          reference fasta.
optional arguments:
      -o, --out-prefix=NAME         prefix of phasing result. default: modcall_result
      -t, --threads=Num             number of thread. default:1
phasing arguments:
      -m, --modThreshold=[0~1]      value extracted from MM tag and ML tag.
                                    above the threshold means modification occurred. default: 0.8
      -u, --unModThreshold=[0~1]    value extracted from MM tag and ML tag.
                                    above the threshold means no modification occurred. default: 0.2
      -e, --heterRatio=[0~1]        modified and unmodified scales.
                                    a higher threshold means that the two quantities need to be closer. default: 0.7
      -i, --noiseRatio=[0~1]        not being judged as modified and unmodified is noise.
                                    higher threshold means lower noise needs. default: 0.2
      -a, --connectAdjacent=Num     connect adjacent N METHs. default:6
      -c, --connectConfidence=[0~1] determine the confidence of phasing two ASMs.
                                    higher threshold requires greater consistency in the reads. default: 0.9
```

---
## Input Preparation
#### Generate reference index
Index the reference genome with [samtools](https://github.com/samtools/samtools).
```
samtools faidx reference.fasta
```
#### Generate alignment and index files
Produce read-to-reference alignment via [minimap2](https://github.com/lh3/minimap2) and sort/index the bam by [samtools](https://github.com/samtools/samtools).
```
# generate alignment flie with minimap2 according to the sequencing platform e.g. map-pb/map-ont/map-hifi
# Note that the MD-tag is required by sniffles (–MD).
minimap2 --MD -ax map-ont -t 10 reference.fasta reads.fastq -o alignment.sam

# sort alignment file
samtools sort -@ 10 alignment.sam -o alignment.bam

# index alignment file
samtools index -@ 10 alignment.bam
```
#### Generate single nucleotide polymorphism (SNP) file
e.g. [PEPPER-Margin-DeepVariant](https://github.com/kishwarshafin/pepper) or [Clair3](https://github.com/HKU-BAL/Clair3) pipeline.
```
INPUT_DIR={input data path}
OUTPUT_DIR={output data path}
BAM=alignment.bam
REF=reference.fasta
THREADS=10

sudo docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/pepper_deepvariant:r0.7 \
run_pepper_margin_deepvariant call_variant \
-b "${INPUT_DIR}/${BAM}" \
-f "${INPUT_DIR}/${REF}" \
-o "${OUTPUT_DIR}" \
-t "${THREADS}" \
--ont_r9_guppy5_sup

# --ont_r9_guppy5_sup is preset for ONT R9.4.1 Guppy 5 "Sup" basecaller
# for ONT R10.4 Q20 reads: --ont_r10_q20
# for PacBio-HiFi reads: --hifi
```
#### Generate Structural variation (SV) file
e.g. [sniffles](https://github.com/fritzsedlazeck/Sniffles) or [CuteSV](https://github.com/tjiangHIT/cuteSV).
```
# In sniffles1 please specofic --num_reads_report -1. For sniffles2 please specify --output-rnames instead.
sniffles -t 10 --num_reads_report -1 -m alignment.bam -v SV.vcf # for sniffles1
sniffles -t 10 --output-rnames -1 -m alignment.bam -v SV.vcf # for sniffles2

# cuteSV command for PacBio CLR data:
cuteSV alignment.bam reference.fasta SV.vcf work_dir --report_readid --genotype

# additional platform-specific parameters suggested by cuteSV
# PacBio CLR data: 
--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5
# PacBio CCS(HIFI) data: 
--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5
# ONT data: 
--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3
```

#### [Modified bam process](https://github.com/nanoporetech/dorado/issues/145)
The `-T` parameter in samtools fastq allows for the extraction of tags from the BAM file and stores them in the header of the FASTQ file. Please ensure that the BAM file includes both `MM` and `ML` tags.
```
samtools fastq -T '*' methylcall.raw.bam > methylcall.raw.fastq
```

The `-y` option in minimap2 allows appending tags stored in the FASTQ header to the alignment.
```
minimap2 -ax map-ont -y reference.fasta methylcall.raw.fastq 
```

---
## Comparison with other SNP-phasing programs
LongPhase is 10x faster than WhatsHap and Margin and produces much larger blocks when tested on HG002, HG003,and HG004.
![btac058f3](https://github.com/twolinin/longphase/assets/6086073/af3a75a1-6268-4700-9dcc-4a6f34e86f7a)

---
## Citation
Jyun-Hong Lin, Liang-Chi Chen, Shu-Qi Yu and Yao-Ting Huang, [LongPhase: an ultra-fast chromosome-scale phasing algorithm for small and large variants](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btac058/6519151), Bioinformatics, 2022.

---
## Contact
Yao-Ting Huang, ythuang at cs.ccu.edu.tw


