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



