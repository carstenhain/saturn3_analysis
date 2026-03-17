#!/usr/bin/env python3
"""
Convert BND (breakend) calls representing DEL, INV, DUP into proper SV type calls.

This script reads a VCF file where long SV events are represented as single 
breakend calls and converts them to proper DEL, INV, and DUP calls.
"""

import cyvcf2 # type: ignore
import argparse
import sys
import re


def parse_bnd_alt(alt_string:str) -> tuple:
    """
    Parses the ALT string of a BND record to extract mate chromosome, position, and orientation.

    Args:
        alt_string (str): ALT field string from a BND record.

    Returns:
        tuple: A 3-tuple containing:
            - mate_chr (str | None): Chromosome name of the breakend mate, or None if parsing fails
            - mate_pos (int | None): Position of the breakend mate, or None if parsing fails
            - orientation (str | None): Orientation type ('right_after', 'left_after', 
                                       'left_before', 'right_before'), or None if parsing fails
    """

    # Pattern for BND notation
    patterns = [
        (r't\[(.+):(\d+)\[', 'right_after'),      # t[p[
        (r't\](.+):(\d+)\]', 'left_after'),       # t]p]
        (r'\](.+):(\d+)\]t', 'left_before'),      # ]p]t
        (r'\[(.+):(\d+)\[t', 'right_before')      # [p[t
    ]
    
    for pattern, orientation in patterns:
        match = re.search(pattern, alt_string.replace('N', 't'))
        if match:
            mate_chr = match.group(1)
            mate_pos = int(match.group(2))
            return mate_chr, mate_pos, orientation
    
    return None, None, None


def classify_bnd_as_sv(record):
    """
    Classifies a BND record as DEL, INV, or DUP based on the ALT field and orientation of the breakends.

    Args:
        record (_type_): _description_

    Returns:
        _type_: _description_
    """
    if record.INFO.get("SVTYPE") != "BND":
        return None, None
    
    alt_string = record.ALT[0] if record.ALT else ""
    mate_chr, mate_pos, orientation = parse_bnd_alt(alt_string)
    
    if mate_chr is None or mate_chr != record.CHROM:
        # Not on same chromosome, can't be DEL/INV/DUP
        return None, None
    
    pos = record.POS
    
    # Ensure pos < mate_pos for consistency
    if pos > mate_pos:
        pos, mate_pos = mate_pos, pos
    
    # Determine SV type based on orientation
    # DEL: ].p] or [p.[  - breakends point away from each other
    # DUP: [.p[ or ]p.]  - breakends point toward each other  
    # INV: mix of orientations
    
    if orientation in ['left_after', 'right_before']:
        # Deletion-like pattern
        return "DEL", mate_pos
    elif orientation in ['right_after', 'left_before']:
        # Could be DUP or INV, need more sophisticated logic
        # For now, classify as DUP
        return "DUP", mate_pos
    
    return None, None


def convert_bnd_to_sv(input_vcf, output_vcf, min_size=50):
    """
    Convert BND calls to DEL/INV/DUP calls.
    
    Args:
        input_vcf: Path to input VCF file
        output_vcf: Path to output VCF file
        min_size: Minimum SV size to consider (default: 50bp)
    """
    vcf = cyvcf2.VCF(input_vcf)
    
    # Modify header to ensure SVTYPE and END are defined
    if 'SVTYPE' not in vcf.contains:
        vcf.add_info_to_header({
            'ID': 'SVTYPE',
            'Description': 'Type of structural variant',
            'Type': 'String',
            'Number': '1'
        })
    
    if 'END' not in vcf.contains:
        vcf.add_info_to_header({
            'ID': 'END',
            'Description': 'End position of the variant',
            'Type': 'Integer',
            'Number': '1'
        })
    
    if 'SVLEN' not in vcf.contains:
        vcf.add_info_to_header({
            'ID': 'SVLEN',
            'Description': 'Length of structural variant',
            'Type': 'Integer',
            'Number': '1'
        })
    
    writer = cyvcf2.Writer(output_vcf, vcf)
    
    converted_count = 0
    kept_count = 0
    skipped_count = 0
    
    for record in vcf:
        sv_type, end_pos = classify_bnd_as_sv(record)
        
        if sv_type and end_pos:
            sv_len = abs(end_pos - record.POS)
            
            if sv_len >= min_size:
                # Convert BND to proper SV type
                record.INFO["SVTYPE"] = sv_type
                record.INFO["END"] = end_pos
                record.INFO["SVLEN"] = sv_len if sv_type != "INV" else -sv_len
                
                # Update ALT field to symbolic notation
                record.ALT = [f"<{sv_type}>"]
                
                converted_count += 1
                writer.write_record(record)
            else:
                skipped_count += 1
        else:
            # Keep non-BND records or BND that can't be classified
            kept_count += 1
            writer.write_record(record)
    
    writer.close()
    vcf.close()
    
    print(f"Conversion complete:")
    print(f"  Converted BND to SV: {converted_count}")
    print(f"  Kept unchanged: {kept_count}")
    print(f"  Skipped (too small): {skipped_count}")
    print(f"Output written to: {output_vcf}")


def main():
    parser = argparse.ArgumentParser(
        description='Convert BND calls to DEL/INV/DUP calls in VCF files'
    )
    parser.add_argument(
        'input_vcf',
        help='Input VCF file with BND calls'
    )
    parser.add_argument(
        'output_vcf',
        help='Output VCF file with converted SV calls'
    )
    parser.add_argument(
        '--min-size',
        type=int,
        default=50,
        help='Minimum SV size in bp (default: 50)'
    )
    
    args = parser.parse_args()
    
    try:
        convert_bnd_to_sv(args.input_vcf, args.output_vcf, args.min_size)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
