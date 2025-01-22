import pytest
import pandas as pd
from pathlib import Path
from .variant_mapper_4 import VariantMethylationMapper, VariantCounts

def test_normalize_sample_name():
    """Test sample name normalization"""
    mapper = VariantMethylationMapper("test", "H1")
    assert mapper._normalize_sample_name("SAMPLE-01") == "SAMPLE01"
    assert mapper._normalize_sample_name("SAMPLE_01") == "SAMPLE01"
    assert mapper._normalize_sample_name("SAMPLE01") == "SAMPLE01"

def test_get_samples():
    """Test sample string parsing"""
    mapper = VariantMethylationMapper("test", "H1")
    
    # Test various formats
    assert mapper._get_samples("SAMPLE-01,SAMPLE-02") == {"SAMPLE01", "SAMPLE02"}
    assert mapper._get_samples("SAMPLE-01;SAMPLE-02") == {"SAMPLE01", "SAMPLE02"}
    assert mapper._get_samples(".") == set()
    assert mapper._get_samples("") == set()

def test_variant_counts():
    """Test VariantCounts class"""
    counts = VariantCounts()
    
    # Test initialization
    assert counts.total_meth == 0
    assert counts.total_unmeth == 0
    assert counts.meth_with_var == []
    
    # Test property
    counts.total_meth = 5
    counts.total_unmeth = 3
    counts.meth_with_var = ["sample1:1|0", "sample2:1|1"]
    assert counts.total_samples == 8

def test_methylation_scoring():
    """Test methylation association scoring"""
    mapper = VariantMethylationMapper("test", "H1")
    
    # Create test counts
    counts = VariantCounts()
    counts.total_meth = 20
    counts.total_unmeth = 20
    counts.meth_with_var = ["sample1:1|0"] * 15  # 75% methylated with variant
    counts.unmeth_with_var = ["sample2:1|0"] * 5  # 25% unmethylated with variant
    counts.no_var_meth = 5
    counts.no_var_unmeth = 15
    
    # Test association calculation
    result = mapper.determine_methylation_association(counts)
    assert result["classification"] == "methylated_associated"
    assert "enrichment_score" in result["methylated_metrics"]
    assert result["methylated_metrics"]["confidence"] > 0

@pytest.fixture
def sample_bed_content():
    """Create a sample BED file for testing"""
    return pd.DataFrame({
        'chrom': ['chr1', 'chr1'],
        'start': [1000, 2000],
        'end': [1500, 2500],
        'methylated_samples': ['SAMPLE-01,SAMPLE-02', 'SAMPLE-03'],
        'unmethylated_samples': ['SAMPLE-04', 'SAMPLE-05,SAMPLE-06']
    })

def test_process_region(tmp_path, sample_bed_content):
    """Test region processing"""
    # Save sample content to temp file
    bed_file = tmp_path / "test.bed"
    sample_bed_content.to_csv(bed_file, sep='\t', index=False)
    
    mapper = VariantMethylationMapper("test", "H1")
    results = mapper.process_region(sample_bed_content.iloc[0])
    
    # Basic structure tests - actual results will depend on VCF content
    assert isinstance(results, list)
