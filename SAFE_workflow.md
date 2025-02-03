SAFE workflow

Generate communal NYCRD github with shared credentials

1)	Edit files on github (scripts, config files, paths and commands)
2)	Clone into SAFE HPC
3)	Copy and paste from commands list

AlleleStacker Configuration To Dos:

Step 1: CpG Pileup QC (requires downloading hg38 reference genome) : https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/
Step 2: Segmentation
Step 3: Candidate Region Generation 
Step 4: Allele Stack Generation
Step 5: Variant Mapping


Step 1: CpG Pileup QC
├── pileup_QC
│   ├── config.yaml
│   ├── sample_list.txt
│   └── submit_pileup_qc.sh

Step 2: Segmentation
├── segmentation
│   ├── a_segment_samples.sh *Update MethBat path 
│   ├── b_extract-regions.sh
│   ├── plot_region_distribution_counts.sh
│   ├── plot_region_distribution_size.sh
│   ├── plot_segmentation_group.sh
│   └── plot_segmentation_sample.sh

Step 3: Candidate Region Generation 
├── candidate_regions
│   ├── a_consensus_regions_meth.sh
│   ├── a_consensus_regions.sh
│   └── b_filter-consensus.sh

Step 4: Allele Stack Generation *revisit and compare with originals
├── igv_viewing
│   ├── a_IGV_all-samples.sh
│   ├── b_consensus-IGV.sh
│   └── c_IGV-segmentation.sh

Step 5: Variant Mapping *needs an update to the file_paths folder
└── variant_map
    ├── a_merge_small_vars.sh
    ├── b_merge_sv_vars.sh
    ├── c_merge_CNVs_vars.sh
    ├── d_merge_repeats_vars.sh
    ├── e_map_variants.sh

