from dataclasses import dataclass
from scipy.stats import fisher_exact
import numpy as np

# Constants for confidence thresholds
HIGH_CONFIDENCE = 0.8
MEDIUM_CONFIDENCE = 0.6
LOW_CONFIDENCE = 0.4  # Added low confidence threshold

@dataclass
class VariantCounts:
    """Track sample counts and overlap patterns"""
    total_meth: int = 0
    total_unmeth: int = 0
    total_no_data: int = 0
    meth_with_var: int = 0
    unmeth_with_var: int = 0
    no_data_with_var: int = 0  # No methylation data, has variant
    no_meth_with_var: int = 0   # No methylation data, has variant (duplicate for clarity)
    no_var_meth: int = 0      # No variant data, methylated
    no_var_unmeth: int = 0    # No variant data, unmethylated
    no_data_both: int = 0     # No data for both methylation and variant

    @property
    def total_samples(self) -> int:
        return (self.total_meth + self.total_unmeth + self.total_no_data +
                self.no_meth_with_var + self.no_var_meth + self.no_var_unmeth + self.no_data_both)

    @property
    def total_with_variant(self) -> int:
        return self.meth_with_var + self.unmeth_with_var + self.no_meth_with_var

    @property
    def total_methylated(self) -> int:
        return self.total_meth + self.no_var_meth

    @property
    def total_unmethylated(self) -> int:
        return self.total_unmeth + self.no_var_unmeth

    def get_meth_ratio(self, include_no_variant_data=False) -> float:
        """
        Calculate methylated sample ratio.

        Args:
            include_no_variant_data (bool): Whether to include samples with no variant data in the calculation.
        """
        total_meth = self.total_methylated if include_no_variant_data else self.total_meth
        return self.meth_with_var / total_meth if total_meth > 0 else 0.0

    def get_unmeth_ratio(self, include_no_variant_data=False) -> float:
        """
        Calculate unmethylated sample ratio.

        Args:
            include_no_variant_data (bool): Whether to include samples with no variant data in the calculation.
        """
        total_unmeth = self.total_unmethylated if include_no_variant_data else self.total_unmeth
        return self.unmeth_with_var / total_unmeth if total_unmeth > 0 else 0.0
