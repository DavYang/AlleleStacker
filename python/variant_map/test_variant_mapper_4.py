import pytest
from pathlib import Path
import pandas as pd
from .variant_mapper_4 import VariantMethylationMapper, VariantCounts

@pytest.fixture
def test_mapper():
    return VariantMethylationMapper(
        output_prefix="test_output",
        haplotype="H1"
    )

@pytest.fixture
def sample_variant_counts():
    counts = VariantCounts()
    counts.total_meth = 100
    counts.total_unmeth = 100
    counts.meth_with_var = ["sample1:1|0", "sample2:1|1"]
    counts.unmeth_with_var = ["sample3:1|0"]
    counts.no_var_meth = 98
    counts.no_var_unmeth = 99
    return counts

def test_variant_counts_initialization():
    counts = VariantCounts()
    assert counts.total_meth == 0
    assert counts.total_unmeth == 0
    assert counts.meth_with_var == []
    assert counts.unmeth_with_var == []
    assert counts.no_meth_with_var == []

def test_normalize_sample_name(test_mapper):
    assert test_mapper._normalize_sample_name("sample-1") == "sample1"
    assert test_mapper._normalize_sample_name("sample_2") == "sample2"
    assert test_mapper._normalize_sample_name("sample-3_test") == "sample3test"

def test_get_samples(test_mapper):
    sample_str = "sample-1,sample_2;sample-3"
    expected = {"sample1", "sample2", "sample3"}
    assert test_mapper._get_samples(sample_str) == expected

def test_calculate_state_enrichment_score(test_mapper, sample_variant_counts):
    # Test methylated state
    meth_score = test_mapper.calculate_state_enrichment_score(sample_variant_counts, 'methylated')
    assert meth_score is not None
    assert 'enrichment_score' in meth_score
    assert 'confidence' in meth_score
    assert 'quality' in meth_score
    
    # Test unmethylated state
    unmeth_score = test_mapper.calculate_state_enrichment_score(sample_variant_counts, 'unmethylated')
    assert unmeth_score is not None
    assert 'enrichment_score' in unmeth_score
    assert 'confidence' in unmeth_score
    assert 'quality' in unmeth_score

def test_determine_methylation_association(test_mapper, sample_variant_counts):
    association = test_mapper.determine_methylation_association(sample_variant_counts)
    assert association is not None
    assert 'classification' in association
    assert 'methylated_metrics' in association
    assert 'unmethylated_metrics' in association
    assert 'primary_metrics' in association

def test_insufficient_data_handling(test_mapper):
    small_counts = VariantCounts()
    small_counts.total_meth = 5
    small_counts.total_unmeth = 5
    
    association = test_mapper.determine_methylation_association(small_counts)
    assert association['classification'] == 'insufficient_data'
    assert association['reason'] == 'insufficient_samples'
