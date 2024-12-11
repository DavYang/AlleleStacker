### Date: 12/2/24

Note:
Due to how the current CpG tools does methylation pielups at heterozygous variants (CpG decay variants) the methylation at the variant is reported as a CpG site with 0% methylation probability. When MethBat Segment is used with haplotype segmentation enabled, this results in regions with CpG decay to be called as unmethylated, when in actuality, these alleles lack CpGs. 

Additionally, I will re-segment using the previous parameters (min-cpg=5, max-gap=500) and re-work the bash script strucutre so that the variables for each script can be set within the python notebook

### Date: 11/26/24

Re-working variant mapper using a small subset of the data to optimize the analysis and output structure: 


##### Variant Mapper Results Documentation

##### Column Headers Description

#### Region Information
| Column | Description |
|--------|-------------|
| chrom | Chromosome containing the methylation region |
| start | Start coordinate of the methylation region |
| end | End coordinate of the methylation region |

#### Variant Information
| Column | Description |
|--------|-------------|
| variant_chr | Chromosome containing the variant |
| variant_start | Start position of the variant |
| variant_end | End position of the variant (important for structural variants) |
| variant_id | Unique identifier for the variant (from VCF or constructed as [type]_[chr]_[pos]) |
| type | Type of variant: SNP, CNV (Copy Number Variant), or SV (Structural Variant) |

#### Sample Fractions
| Column | Description |
|--------|-------------|
| meth_fraction | Ratio of methylated samples containing variant (format: "n/N" where n = samples with variant, N = total methylated samples) |
| unmeth_fraction | Ratio of unmethylated samples containing variant (format: "n/N" where n = samples with variant, N = total unmethylated samples) |

#### Sample Details
| Column | Description |
|--------|-------------|
| meth_samples | List of methylated samples containing the variant with their genotypes (format: "SAMPLE_ID:GENOTYPE", comma-separated) |
| unmeth_samples | List of unmethylated samples containing the variant with their genotypes (format: "SAMPLE_ID:GENOTYPE", comma-separated) |

#### Genotype Notation
- For SNPs/SVs: "SAMPLE:A|B" where A and B are alleles (0=reference, 1=alternate)
- For CNVs: "SAMPLE:CN=X" where X is the copy number
- Multiple samples are separated by commas
- "." indicates no samples in that category

#### Example Entry
```
chr1  100000  120000  chr1  105678  105679  SNP_chr1_105678  SNP  2/8  5/12  SAMPLE1:1|0,SAMPLE4:1|0  SAMPLE8:1|0,SAMPLE9:1|0
```
This shows a SNP in a methylation region on chr1, present in 2 out of 8 methylated samples and 5 out of 12 unmethylated samples.



### Date : 11/21/24

Updates : 
- Variant mapping runs for all VCF file types, not production ready (not sure if mapping is 100% accurate, variant IDs need to be formatted better)
- Merged variants need to be QCed to remove low quality variants.

### Date : 11/20/24

Goals:

- Map variants using phased data from pb-WGS-wdl outputs (these are all phased)
- This will require merging variant files for each type and applying QC.
- Map variants for all types available.
- Include count of input regins, input variants, and outputs.


#### Variant file preparation

1) Compile sample data from sample directories

- sample_phased_small_variant_vcfs
- sample_phased_sv_vcfs
- hificnv_vcfs
- trgt_repeat_vcf
* See "/gs/gsfs0/shared-lab/greally-lab/David/6_base-seq_SC1/WGS-analysis/outputs_compiled/20240501-results" for compiled data

2) Merge VCFs and apply soft QC
3) Map to candidate regions




### Date : 11/19/24

This test run will use a new set of intial segmentation paramters to enhance resolution (purposeful over segmentation) and boost the generation of query regions.

- MethBatv0.13.2 segment : 
    --min-cpgs 2
    --max-gap 200

###### Variant mapping details:

- Inital attempt was to map all variant types available (SNVs, indels, SVs, repeats, CNVs) but ran into issues:
        - See "./bash/4a_map_variants.sh"

- Proceeded with using just small variants (deepvariant) cohort VCF, see "./bash/4b_map_small_variants.sh"
-  For this analysis, I am using the outputs generated as part of the ColoRSDb pipeline, which omit's SPM232 (low geno, ambiguous sex) and also to note is that SPM276-1 is written as SPM276_1 in VCFs (this is taken into account in the parsing step for mapping variants)
* These variants are not phased, 

- Took ~53 minutes to complete.
- 244,379 small variants identified.
- * The script only looked at haplotype 1, needs correction. 

### Next Steps:

- Map variants using phased data from pb-WGS-wdl outputs (these are all phased)
- This will require merging variant files for each type and applying QC.
- Map variants for all types available.
- Include count of input regins, input variants, and outputs.



