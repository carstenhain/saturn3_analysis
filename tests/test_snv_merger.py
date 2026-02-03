#!/usr/bin/env python3
"""
Unit tests for SNV merger
"""

import unittest
import tempfile
import os
from pathlib import Path
import sys

# Add parent directory to path to import snv_merger
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from scripts.snv_merger import Variant, VCFParser, SNVMerger


class TestVariant(unittest.TestCase):
    """Test Variant class"""
    
    def test_variant_creation(self):
        """Test creating a variant"""
        v = Variant(chrom='chr1', pos=100, ref='A', alt='T', caller='test')
        self.assertEqual(v.chrom, 'chr1')
        self.assertEqual(v.pos, 100)
        self.assertEqual(v.ref, 'A')
        self.assertEqual(v.alt, 'T')
        self.assertEqual(v.caller, 'test')
    
    def test_variant_key(self):
        """Test variant key generation"""
        v1 = Variant(chrom='chr1', pos=100, ref='A', alt='T', caller='caller1')
        v2 = Variant(chrom='chr1', pos=100, ref='A', alt='T', caller='caller2')
        self.assertEqual(v1.key(), v2.key())
    
    def test_variant_equality(self):
        """Test variant equality"""
        v1 = Variant(chrom='chr1', pos=100, ref='A', alt='T', caller='caller1')
        v2 = Variant(chrom='chr1', pos=100, ref='A', alt='T', caller='caller2')
        self.assertEqual(v1, v2)


class TestVCFParser(unittest.TestCase):
    """Test VCF parser"""
    
    def test_parse_vcf(self):
        """Test parsing a VCF file"""
        # Create a temporary VCF file
        vcf_content = """##fileformat=VCFv4.2
##source=TestCaller
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	T	30.0	PASS	DP=50
chr1	200	.	G	C	40.0	PASS	DP=60
chr2	300	.	C	G	20.0	LowQual	DP=30
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write(vcf_content)
            vcf_file = f.name
        
        try:
            variants = VCFParser.parse_vcf(Path(vcf_file), 'TestCaller')
            
            self.assertEqual(len(variants), 3)
            
            # Check first variant
            v1 = variants[0]
            self.assertEqual(v1.chrom, 'chr1')
            self.assertEqual(v1.pos, 100)
            self.assertEqual(v1.ref, 'A')
            self.assertEqual(v1.alt, 'T')
            self.assertEqual(v1.quality, 30.0)
            self.assertEqual(v1.caller, 'TestCaller')
            
            # Check second variant
            v2 = variants[1]
            self.assertEqual(v2.chrom, 'chr1')
            self.assertEqual(v2.pos, 200)
            
            # Check third variant with filter
            v3 = variants[2]
            self.assertEqual(v3.chrom, 'chr2')
            self.assertIn('LowQual', v3.filters)
        
        finally:
            os.unlink(vcf_file)


class TestSNVMerger(unittest.TestCase):
    """Test SNV merger"""
    
    def test_merge_single_caller(self):
        """Test merging variants from a single caller"""
        merger = SNVMerger(min_callers=1)
        
        variants = [
            Variant(chrom='chr1', pos=100, ref='A', alt='T', caller='caller1', quality=30.0),
            Variant(chrom='chr1', pos=200, ref='G', alt='C', caller='caller1', quality=40.0),
        ]
        
        merger.add_variants(variants)
        merged = merger.merge()
        
        self.assertEqual(len(merged), 2)
    
    def test_merge_multiple_callers(self):
        """Test merging variants from multiple callers"""
        merger = SNVMerger(min_callers=1)
        
        # Variants from caller1
        variants1 = [
            Variant(chrom='chr1', pos=100, ref='A', alt='T', caller='caller1', quality=30.0),
            Variant(chrom='chr1', pos=200, ref='G', alt='C', caller='caller1', quality=40.0),
        ]
        
        # Variants from caller2
        variants2 = [
            Variant(chrom='chr1', pos=100, ref='A', alt='T', caller='caller2', quality=35.0),
            Variant(chrom='chr1', pos=300, ref='C', alt='G', caller='caller2', quality=25.0),
        ]
        
        merger.add_variants(variants1)
        merger.add_variants(variants2)
        merged = merger.merge()
        
        self.assertEqual(len(merged), 3)
        
        # Check that overlapping variant has both callers
        overlapping = [v for v in merged if v.pos == 100][0]
        self.assertIn('caller1', overlapping.info['CALLERS'])
        self.assertIn('caller2', overlapping.info['CALLERS'])
        self.assertEqual(overlapping.info['CALLER_COUNT'], '2')
    
    def test_min_callers_filter(self):
        """Test filtering by minimum callers"""
        merger = SNVMerger(min_callers=2)
        
        # Variants from caller1
        variants1 = [
            Variant(chrom='chr1', pos=100, ref='A', alt='T', caller='caller1', quality=30.0),
            Variant(chrom='chr1', pos=200, ref='G', alt='C', caller='caller1', quality=40.0),
        ]
        
        # Variants from caller2
        variants2 = [
            Variant(chrom='chr1', pos=100, ref='A', alt='T', caller='caller2', quality=35.0),
        ]
        
        merger.add_variants(variants1)
        merger.add_variants(variants2)
        merged = merger.merge()
        
        # Only the overlapping variant should remain
        self.assertEqual(len(merged), 1)
        self.assertEqual(merged[0].pos, 100)
    
    def test_min_quality_filter(self):
        """Test filtering by minimum quality"""
        merger = SNVMerger(min_callers=1, min_quality=35.0)
        
        variants = [
            Variant(chrom='chr1', pos=100, ref='A', alt='T', caller='caller1', quality=30.0),
            Variant(chrom='chr1', pos=200, ref='G', alt='C', caller='caller1', quality=40.0),
        ]
        
        merger.add_variants(variants)
        merged = merger.merge()
        
        # Only the high-quality variant should remain
        self.assertEqual(len(merged), 1)
        self.assertEqual(merged[0].pos, 200)
    
    def test_write_vcf(self):
        """Test writing merged variants to VCF"""
        merger = SNVMerger(min_callers=1)
        
        variants = [
            Variant(chrom='chr1', pos=100, ref='A', alt='T', caller='caller1', quality=30.0),
        ]
        
        merger.add_variants(variants)
        merged = merger.merge()
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            output_file = f.name
        
        try:
            merger.write_vcf(Path(output_file), merged)
            
            # Read back the file
            with open(output_file, 'r') as f:
                content = f.read()
            
            self.assertIn('##fileformat=VCFv4.2', content)
            self.assertIn('##source=SNVMerger', content)
            self.assertIn('chr1\t100', content)
        
        finally:
            os.unlink(output_file)


if __name__ == '__main__':
    unittest.main()
