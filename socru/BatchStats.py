"""Batch-level statistics from multiple Socru analysis results.

This module computes aggregate statistics across a collection of AnalysisResult
dictionaries, including type distributions, quality summaries, confidence
scoring, QC flag aggregation, and outlier detection.

Classes:
    BatchStats: Compute batch-level statistics from multiple AnalysisResult objects.
"""

import math
from collections import Counter
from typing import Any, Dict, List, Optional


class BatchStats:
    """Compute batch-level statistics from multiple AnalysisResult objects.

    Attributes:
        results: List of AnalysisResult.to_dict() dictionaries.
    """

    def __init__(self, results: List[Dict[str, Any]]) -> None:
        """Initialize with a list of result dictionaries.

        Args:
            results: List of AnalysisResult.to_dict() dicts, each containing
                keys such as gs_type, quality, confidence_score, qc_flags,
                genome_file, and genome_length.
        """
        self.results = results

    def type_distribution(self) -> Dict[str, int]:
        """Return dict of {gs_type: count} sorted by count descending.

        Returns:
            Dictionary mapping GS type strings to their occurrence counts,
            ordered from most frequent to least frequent.
        """
        counts = Counter(r.get("gs_type", "Unknown") for r in self.results)
        return dict(counts.most_common())

    def quality_summary(self) -> Dict[str, int]:
        """Return dict of {quality_level: count}.

        Returns:
            Dictionary with keys GREEN, AMBER, RED and their respective counts.
            Keys with zero counts are included.
        """
        counts: Dict[str, int] = {"GREEN": 0, "AMBER": 0, "RED": 0}
        for r in self.results:
            quality = r.get("quality", "").upper()
            if quality in counts:
                counts[quality] += 1
        return counts

    def mean_confidence(self) -> Optional[float]:
        """Return mean confidence score across all results.

        Returns:
            Mean confidence score as a float, or None if no results contain
            a confidence_score value.
        """
        scores = [
            r["confidence_score"]
            for r in self.results
            if r.get("confidence_score") is not None
        ]
        if not scores:
            return None
        return sum(scores) / len(scores)

    def flag_summary(self) -> Dict[str, int]:
        """Return dict of {flag_code: count} for all QC flags.

        Aggregates QC flag codes across all results.

        Returns:
            Dictionary mapping flag code strings to their total occurrence
            counts, sorted by count descending.
        """
        counts: Counter = Counter()
        for r in self.results:
            for flag in r.get("qc_flags", []):
                code = flag.get("code", "UNKNOWN")
                counts[code] += 1
        return dict(counts.most_common())

    def outlier_assemblies(self, z_threshold: float = 2.0) -> List[str]:
        """Return list of assembly names that are statistical outliers.

        An assembly is flagged as an outlier if its genome_length or
        confidence_score deviates from the batch mean by more than
        z_threshold standard deviations.

        Args:
            z_threshold: Number of standard deviations from the mean to
                consider a value an outlier. Defaults to 2.0.

        Returns:
            List of genome_file strings for outlier assemblies.
        """
        if len(self.results) < 2:
            return []

        outliers: set = set()

        for key in ("genome_length", "confidence_score"):
            values = [
                (r.get("genome_file", ""), r.get(key))
                for r in self.results
                if r.get(key) is not None
            ]
            if len(values) < 2:
                continue

            nums = [v for _, v in values]
            mean = sum(nums) / len(nums)
            variance = sum((x - mean) ** 2 for x in nums) / len(nums)
            std = math.sqrt(variance)

            if std == 0:
                continue

            for name, val in values:
                z = abs(val - mean) / std
                if z > z_threshold:
                    outliers.add(name)

        return sorted(outliers)
