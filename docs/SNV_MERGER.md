# SNV Merger

A tool to merge Single Nucleotide Variant (SNV) calls from multiple variant callers.

## Overview

The SNV Merger combines variant calls from different variant callers (e.g., GATK, MuTect2, VarScan) into a unified VCF file. This is useful for:

- Increasing confidence in variant calls by requiring multiple callers to agree
- Creating a comprehensive set of variants from multiple sources
- Tracking which callers detected each variant

## Features

- **VCF Format Support**: Parses standard VCF files from variant callers
- **Flexible Merging**: Configure minimum number of callers and quality thresholds
- **Caller Tracking**: Output includes which callers detected each variant
- **Quality-based Selection**: When multiple callers detect the same variant, uses the highest quality score
- **Standard Output**: Generates standard VCF format output

## Installation

No special installation required. Just ensure Python 3.7+ is installed.

## Usage

### Basic Usage

Merge variants from two callers with default settings:

```bash
python scripts/snv_merger.py \
    caller1:/path/to/caller1.vcf \
    caller2:/path/to/caller2.vcf \
    -o merged_variants.vcf
```

### Advanced Usage

Require at least 2 callers to agree and minimum quality of 30:

```bash
python scripts/snv_merger.py \
    gatk:/data/gatk_variants.vcf \
    mutect2:/data/mutect2_variants.vcf \
    varscan:/data/varscan_variants.vcf \
    -o merged_variants.vcf \
    --min-callers 2 \
    --min-quality 30.0
```

### Command-line Options

- `vcf_files`: One or more VCF files in format `caller_name:/path/to/file.vcf`
- `-o, --output`: Output VCF file for merged variants (required)
- `--min-callers`: Minimum number of callers required to support a variant (default: 1)
- `--min-quality`: Minimum quality score for a variant (default: 0.0)

## Input Format

Input VCF files should follow standard VCF format. Each input must be prefixed with a caller name:

```
caller_name:/path/to/file.vcf
```

Example:
```
gatk:/data/gatk_variants.vcf
mutect2:/data/mutect2_variants.vcf
```

## Output Format

The output is a standard VCF file with additional INFO fields:

- `CALLERS`: Comma-separated list of callers that detected this variant
- `CALLER_COUNT`: Number of callers that detected this variant

Example output line:
```
chr1    12345   .   A   T   45.00   PASS    CALLERS=gatk,mutect2;CALLER_COUNT=2;DP=100
```

## Running Tests

Run the unit tests:

```bash
python -m pytest tests/test_snv_merger.py
```

Or using unittest:

```bash
python -m unittest tests/test_snv_merger.py
```

## Example Workflow

1. Run multiple variant callers on your sample:
   ```bash
   # Run GATK
   gatk HaplotypeCaller -I sample.bam -O gatk_variants.vcf
   
   # Run MuTect2
   gatk Mutect2 -I sample.bam -O mutect2_variants.vcf
   
   # Run VarScan
   varscan mpileup2snp sample.pileup > varscan_variants.vcf
   ```

2. Merge the results:
   ```bash
   python scripts/snv_merger.py \
       gatk:gatk_variants.vcf \
       mutect2:mutect2_variants.vcf \
       varscan:varscan_variants.vcf \
       -o merged_variants.vcf \
       --min-callers 2
   ```

3. Use the merged variants for downstream analysis

## Future Extensions

This tool is designed to be extended to support:

- Structural Variant (SV) merging
- Additional variant types (indels, CNVs, etc.)
- Custom merging strategies
- Integration with workflow managers (Snakemake, Nextflow)

## Contributing

Contributions are welcome! Please submit issues or pull requests.

## License

See LICENSE file for details.
