# saturn3_analysis
Scripts, workflows and visualisation for the PDAC SATURN3 project

## Overview

This repository contains tools and workflows for analyzing variant calls in the SATURN3 project.

## Features

### SNV Merger

Merge Single Nucleotide Variant (SNV) calls from multiple variant callers into a unified output.

**Quick Start:**
```bash
python scripts/snv_merger.py \
    caller1:path/to/caller1.vcf \
    caller2:path/to/caller2.vcf \
    -o merged_variants.vcf \
    --min-callers 2
```

See [docs/SNV_MERGER.md](docs/SNV_MERGER.md) for detailed documentation.

### Future Features

- Structural Variant (SV) merging
- Workflow integration
- Visualization tools

## Project Structure

```
.
├── scripts/           # Analysis scripts
│   └── snv_merger.py  # SNV call merging tool
├── tests/             # Unit tests
├── workflows/         # Workflow definitions (future)
└── docs/              # Documentation
```

## Testing

Run the test suite:
```bash
python -m unittest discover tests
```

## Contributing

Please see individual tool documentation for usage and contribution guidelines.

## License

See LICENSE file for details.
