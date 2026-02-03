#!/usr/bin/env python3
"""
SNV Merger - Merge SNV calls from different variant callers

This module provides functionality to merge Single Nucleotide Variant (SNV)
calls from multiple variant callers into a unified output.
"""

import argparse
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict, Set, Optional


@dataclass
class Variant:
    """Represents a single nucleotide variant"""
    chrom: str
    pos: int
    ref: str
    alt: str
    caller: str
    quality: float = 0.0
    filters: Set[str] = None
    info: Dict[str, str] = None
    
    def __post_init__(self):
        if self.filters is None:
            self.filters = set()
        if self.info is None:
            self.info = {}
    
    def key(self):
        """Returns a unique key for this variant position"""
        return (self.chrom, self.pos, self.ref, self.alt)
    
    def __hash__(self):
        return hash(self.key())
    
    def __eq__(self, other):
        return isinstance(other, Variant) and self.key() == other.key()


class VCFParser:
    """Simple VCF parser for variant calls"""
    
    @staticmethod
    def parse_vcf(vcf_file: Path, caller_name: str) -> List[Variant]:
        """
        Parse a VCF file and return a list of Variant objects
        
        Args:
            vcf_file: Path to the VCF file
            caller_name: Name of the variant caller
            
        Returns:
            List of Variant objects
        """
        variants = []
        
        try:
            with open(vcf_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    
                    # Skip header lines
                    if line.startswith('#'):
                        continue
                    
                    # Skip empty lines
                    if not line:
                        continue
                    
                    # Parse variant line
                    fields = line.split('\t')
                    if len(fields) < 5:
                        continue
                    
                    chrom = fields[0]
                    pos = int(fields[1])
                    ref = fields[3]
                    alt = fields[4]
                    
                    # Get quality if available
                    quality = 0.0
                    if len(fields) > 5 and fields[5] != '.':
                        try:
                            quality = float(fields[5])
                        except ValueError:
                            pass
                    
                    # Get filter information
                    filters = set()
                    if len(fields) > 6 and fields[6] != '.' and fields[6] != 'PASS':
                        filters = set(fields[6].split(';'))
                    
                    # Get INFO field
                    info = {}
                    if len(fields) > 7 and fields[7] != '.':
                        for item in fields[7].split(';'):
                            if '=' in item:
                                key, value = item.split('=', 1)
                                info[key] = value
                            else:
                                info[item] = 'True'
                    
                    variant = Variant(
                        chrom=chrom,
                        pos=pos,
                        ref=ref,
                        alt=alt,
                        caller=caller_name,
                        quality=quality,
                        filters=filters,
                        info=info
                    )
                    variants.append(variant)
        
        except FileNotFoundError:
            print(f"Error: VCF file not found: {vcf_file}", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"Error parsing VCF file {vcf_file}: {e}", file=sys.stderr)
            sys.exit(1)
        
        return variants


class SNVMerger:
    """Merge SNV calls from multiple variant callers"""
    
    def __init__(self, min_callers: int = 1, min_quality: float = 0.0):
        """
        Initialize the SNV merger
        
        Args:
            min_callers: Minimum number of callers that must support a variant
            min_quality: Minimum quality score for a variant to be included
        """
        self.min_callers = min_callers
        self.min_quality = min_quality
        self.variants_by_position = defaultdict(list)
    
    def add_variants(self, variants: List[Variant]):
        """Add variants from a caller to the merger"""
        for variant in variants:
            self.variants_by_position[variant.key()].append(variant)
    
    def merge(self) -> List[Variant]:
        """
        Merge variants from all callers
        
        Returns:
            List of merged Variant objects
        """
        merged_variants = []
        
        for var_key, variant_list in self.variants_by_position.items():
            # Count how many callers support this variant
            caller_count = len(variant_list)
            
            if caller_count < self.min_callers:
                continue
            
            # Get the variant with the highest quality score
            best_variant = max(variant_list, key=lambda v: v.quality)
            
            if best_variant.quality < self.min_quality:
                continue
            
            # Create merged variant with caller information
            callers = [v.caller for v in variant_list]
            
            merged_info = best_variant.info.copy()
            merged_info['CALLERS'] = ','.join(sorted(callers))
            merged_info['CALLER_COUNT'] = str(caller_count)
            
            merged_variant = Variant(
                chrom=best_variant.chrom,
                pos=best_variant.pos,
                ref=best_variant.ref,
                alt=best_variant.alt,
                caller='MERGED',
                quality=best_variant.quality,
                filters=best_variant.filters,
                info=merged_info
            )
            
            merged_variants.append(merged_variant)
        
        # Sort by chromosome and position
        merged_variants.sort(key=lambda v: (v.chrom, v.pos))
        
        return merged_variants
    
    def write_vcf(self, output_file: Path, variants: List[Variant]):
        """
        Write merged variants to a VCF file
        
        Args:
            output_file: Path to output VCF file
            variants: List of variants to write
        """
        with open(output_file, 'w') as f:
            # Write VCF header
            f.write("##fileformat=VCFv4.2\n")
            f.write("##source=SNVMerger\n")
            f.write("##INFO=<ID=CALLERS,Number=.,Type=String,Description=\"Variant callers that detected this variant\">\n")
            f.write("##INFO=<ID=CALLER_COUNT,Number=1,Type=Integer,Description=\"Number of callers that detected this variant\">\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            
            # Write variants
            for variant in variants:
                filter_str = 'PASS' if not variant.filters else ';'.join(sorted(variant.filters))
                
                info_items = []
                for key, value in sorted(variant.info.items()):
                    if value == 'True':
                        info_items.append(key)
                    else:
                        info_items.append(f"{key}={value}")
                info_str = ';'.join(info_items) if info_items else '.'
                
                f.write(f"{variant.chrom}\t{variant.pos}\t.\t{variant.ref}\t{variant.alt}\t"
                       f"{variant.quality:.2f}\t{filter_str}\t{info_str}\n")


def main():
    """Main entry point for SNV merger"""
    parser = argparse.ArgumentParser(
        description='Merge SNV calls from multiple variant callers'
    )
    parser.add_argument(
        'vcf_files',
        nargs='+',
        help='VCF files from different variant callers (format: caller_name:path/to/file.vcf)'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output VCF file for merged variants'
    )
    parser.add_argument(
        '--min-callers',
        type=int,
        default=1,
        help='Minimum number of callers required to support a variant (default: 1)'
    )
    parser.add_argument(
        '--min-quality',
        type=float,
        default=0.0,
        help='Minimum quality score for a variant (default: 0.0)'
    )
    
    args = parser.parse_args()
    
    # Initialize merger
    merger = SNVMerger(min_callers=args.min_callers, min_quality=args.min_quality)
    
    # Parse and add variants from each caller
    for vcf_input in args.vcf_files:
        if ':' not in vcf_input:
            print(f"Error: VCF file must be in format 'caller_name:path/to/file.vcf'", file=sys.stderr)
            print(f"Got: {vcf_input}", file=sys.stderr)
            sys.exit(1)
        
        caller_name, vcf_path = vcf_input.split(':', 1)
        vcf_file = Path(vcf_path)
        
        print(f"Reading variants from {caller_name}: {vcf_file}")
        variants = VCFParser.parse_vcf(vcf_file, caller_name)
        print(f"  Found {len(variants)} variants")
        
        merger.add_variants(variants)
    
    # Merge variants
    print(f"\nMerging variants (min_callers={args.min_callers}, min_quality={args.min_quality})")
    merged_variants = merger.merge()
    print(f"  Merged to {len(merged_variants)} variants")
    
    # Write output
    output_path = Path(args.output)
    print(f"\nWriting merged variants to {output_path}")
    merger.write_vcf(output_path, merged_variants)
    print("Done!")


if __name__ == '__main__':
    main()
