```
# Variant Methylation Mapper

## Purpose

This script analyzes the association between genetic variants and methylation patterns in a set of samples. It takes as input a BED file defining genomic regions with methylation information and VCF files containing variant calls. The script maps variants to these regions, assesses their association with methylation status, and generates scored and prioritized output files.

## Functions

### 1. Data Loading and Preprocessing

* Loads methylation regions and sample information from a BED file.
* Loads variant calls and sample genotypes from VCF files (supports various variant types: SNPs, indels, CNVs, SVs, TRs).
* Normalizes sample names to ensure consistency across files.
* Handles phased and unphased variants.

### 2. Variant Mapping and Overlap Analysis

* Maps variants to methylation regions using `tabix` for efficient retrieval.
* Determines the overlap between variants and methylation states for each sample.
* Tracks 8 possible states for each variant-sample combination:
    - `M+V+`: Methylated with variant
    - `M-V+`: Unmethylated with variant
    - `M+V-`: Methylated without variant
    - `M-V-`: Unmethylated without variant
    - `X(M-)V+`: No methylation data, has variant
    - `X(V-)M+`: No variant data, methylated
    - `X(V-)M-`: No variant data, unmethylated
    - `X(M-)X(V-)`: No data for either methylation or variant

### 3. Association Scoring and Prioritization

* Calculates association metrics for each variant:
    - `data_coverage`: Proportion of samples with methylation data.
    - `meth_association`: Proportion of methylated samples with the variant.
    - `unmeth_association`: Proportion of unmethylated samples with the variant.
* Performs Fisher's exact test to assess the statistical significance of the association (p-value).
* Calculates two separate scores for each variant:
    - `methylation_score`: Reflects the association with methylation.
    - `unmethylation_score`: Reflects the association with unmethylation.
* Applies confidence thresholds to adjust the weight of association strength in the scores based on data coverage.

### 4. Output Generation

* Generates separate output files for methylated and unmethylated variants.
* Includes both methylation and unmethylation scores, p-value, and detailed sample information for each variant in the output.
* Creates a summary statistics file with variant pattern information.

## Scoring Equations

### Methylation Score

```
methylation_score = confidence_weight * (
                     WEIGHT_METH_ASSOCIATION * meth_ratio + 
                     WEIGHT_UNMETH_EXCLUSIVITY * (1 - unmeth_ratio)
                     ) - WEIGHT_P_VALUE * np.log10(p_value)
```

### Unmethylation Score

```
unmethylation_score = confidence_weight * (
                       WEIGHT_METH_ASSOCIATION * unmeth_ratio + 
                       WEIGHT_UNMETH_EXCLUSIVITY * (1 - meth_ratio)
                       ) - WEIGHT_P_VALUE * np.log10(p_value)
```

Where:

* `confidence_weight`: 1.0 (high), 0.7 (medium), 0.4 (low), or 0.1 (very low) based on `data_coverage` thresholds.
* `meth_ratio`: Proportion of methylated samples with the variant.
* `unmeth_ratio`: Proportion of unmethylated samples with the variant.
* `p_value`: P-value from Fisher's exact test.
* `WEIGHT_METH_ASSOCIATION`, `WEIGHT_UNMETH_EXCLUSIVITY`, `WEIGHT_P_VALUE`: Tunable weights for scoring components.

## Rationale for Tunable Parameters

* **Confidence Thresholds:**  Adjust the importance of association strength based on data coverage. Higher coverage = higher confidence = greater weight to association strength.
* **Weights:** Allow customization of the scoring function to prioritize different aspects of the association (methylation association, unmethylation exclusivity, statistical significance).

By adjusting these parameters, users can fine-tune the analysis to their specific needs and data characteristics, enabling more effective identification of candidate variants associated with methylation patterns.
```

**Filename:** `variant_methylation_mapper_summary.txt`
