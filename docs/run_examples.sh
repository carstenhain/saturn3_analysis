#!/bin/bash
# Example usage of SNV merger tool

echo "==== SNV Merger Usage Examples ===="
echo ""

# Example 1: Basic merge with all variants from both callers
echo "Example 1: Basic merge (all variants from both callers)"
python scripts/snv_merger.py \
    caller1:docs/example_caller1.vcf \
    caller2:docs/example_caller2.vcf \
    -o /tmp/merged_all.vcf

echo ""
echo "Result: merged all unique variants from both callers"
echo "Output saved to /tmp/merged_all.vcf"
echo ""

# Example 2: Require at least 2 callers to agree
echo "Example 2: High confidence merge (min 2 callers)"
python scripts/snv_merger.py \
    caller1:docs/example_caller1.vcf \
    caller2:docs/example_caller2.vcf \
    -o /tmp/merged_high_confidence.vcf \
    --min-callers 2

echo ""
echo "Result: only variants detected by both callers"
echo "Output saved to /tmp/merged_high_confidence.vcf"
echo ""

# Example 3: Quality filtering
echo "Example 3: Quality-filtered merge (min quality 40)"
python scripts/snv_merger.py \
    caller1:docs/example_caller1.vcf \
    caller2:docs/example_caller2.vcf \
    -o /tmp/merged_high_quality.vcf \
    --min-quality 40.0

echo ""
echo "Result: only high-quality variants (QUAL >= 40)"
echo "Output saved to /tmp/merged_high_quality.vcf"
echo ""

echo "==== Examples Complete ===="
echo "Check the output files in /tmp/ to see the results"
