import cyvcf2 # type: ignore
import numpy as np # type: ignore
import os
import pandas as pd # type: ignore
import cyvcf2 # type: ignore


def modify_mutect2_vcf (vcf_file_path:str, normal_sample:str, tumor_sample:str) -> str:
    """
    Modifies Mutect2 VCF by adding the information from the FILTER field into the FORMAT field.
    PASS equals a 1, while any other FILTER value leads to a 0
    This writes a new vcf file with the suffix ".modified.vcf" or ".modified.vcf.gz"

    Args:
        vcf_file_path (str): Path to the Mutect2 VCF file
        normal_sample (str): Name of the normal sample in the VCF
        tumor_sample (str): Name of the tumor sample in the VCF

    Raises:
        ValueError: Raises Error if the input VCF file does not end with .vcf or .vcf.gz
    
    Returns:
        str: Path to the modified VCF file
    """

    ### open vcf
    vcf = cyvcf2.VCF(vcf_file_path)

    ### get sample indices
    normal_index = vcf.samples.index(normal_sample)
    tumor_index = vcf.samples.index(tumor_sample)

    ### add FILTER to FORMAT header
    vcf.add_format_to_header({
        "ID": "FILTER",
        "Number": "1",
        "Type": "Integer",
        "Description": "FILTER state"
    })

    ### open modified vcf for writing
    if os.path.basename(vcf_file_path).endswith(".vcf"):
        ending = ".vcf"
    elif os.path.basename(vcf_file_path).endswith(".vcf.gz"):
        ending = ".vcf.gz"
    else:
        raise ValueError("Input VCF file must end with .vcf or .vcf.gz")
    
    ### open vcf writer
    writer = cyvcf2.Writer(f"{vcf_file_path[:-len(ending)]}.modified{ending}", vcf)

    for record in vcf:

        ### skip multiallelic sites
        if len(record.ALT) > 1:
            continue

        ### get filter state
        # filter_state: 1 = PASS, 0 = filtered
        if record.FILTER == None:
            filter_state = 1
        else:
            filter_state = 0

        ### set FILTER format
        filter_states = np.array([[-1], [-1]])
        filter_states[normal_index] = 0
        filter_states[tumor_index] = filter_state
        record.set_format("FILTER", filter_states)

        ### write record
        writer.write_record(record)

    ### close files
    vcf.close()
    writer.close()

    ### index modified vcf
    os.system(f"bcftools index {vcf_file_path[:-len(ending)]}.modified{ending}")

    return f"{vcf_file_path[:-len(ending)]}.modified{ending}"

def modify_strelka_vcf (vcf_file_path:str, normal_sample:str, tumor_sample:str) -> str:
    """
    Modifies Strelka VCF by adding the information from the FILTER field into the FORMAT field and adding GT, AD, AF fields.
    For the filter field PASS equals a 1, while any other FILTER value leads to a 0
    GT is set to 0/0 for normal and 0/1 for tumor
    AD and AF are calculated from the U allele depth fields provided by Strelka
    This writes a new vcf file with the suffix ".modified.vcf" or ".modified.vcf.gz"

    Args:
        vcf_file_path (str): Path to the Strelka VCF file
        normal_sample (str): Name of the normal sample in the VCF
        tumor_sample (str): Name of the tumor sample in the VCF

    Raises:
        ValueError: Raises Error if the input VCF file does not end with .vcf or .vcf.gz

    Returns:
        str: Path to the modified VCF file
    """

    ### open vcf
    vcf = cyvcf2.VCF(vcf_file_path)

    ### get sample indices
    normal_index = vcf.samples.index(normal_sample)
    tumor_index = vcf.samples.index(tumor_sample)

    ### add FILTER to FORMAT header
    vcf.add_format_to_header({
        "ID": "FILTER",
        "Number": "1",
        "Type": "Integer",
        "Description": "FILTER state"
    })
    ### add GT to FORMAT header
    vcf.add_format_to_header({
        "ID": "GT",
        "Number": "1",
        "Type": "String",
        "Description": "Genotype"
    })
    ### add AD to FORMAT header
    vcf.add_format_to_header({
        "ID": "AD",
        "Number": "2",
        "Type": "Integer",
        "Description": "Allelic depths for the ref and alt alleles in the order listed"
    })
    ### add AF to FORMAT header
    vcf.add_format_to_header({
        "ID": "AF",
        "Number": "A",
        "Type": "Float",
        "Description": "Allele fractions of alternate alleles in the tumor"
    })

    ### open modified vcf for writing
    if os.path.basename(vcf_file_path).endswith(".vcf"):
        ending = ".vcf"
    elif os.path.basename(vcf_file_path).endswith(".vcf.gz"):
        ending = ".vcf.gz"
    else:
        raise ValueError("Input VCF file must end with .vcf or .vcf.gz")
    
    ### open vcf writer
    writer = cyvcf2.Writer(f"{vcf_file_path[:-len(ending)]}.modified{ending}", vcf)

    for record in vcf:

        ### skip multiallelic sites
        if len(record.ALT) > 1:
            continue

        ### get filter state
        # filter_state: 1 = PASS, 0 = filtered
        if record.FILTER == None:
            filter_state = 1
        else:
            filter_state = 0

        ### set FILTER format
        filter_states = np.array([[-1], [-1]])
        filter_states[normal_index] = 0
        filter_states[tumor_index] = filter_state
        record.set_format("FILTER", filter_states)

        ### set GT format
        gt_array = np.array([[-1, -1], [-1, -1]])
        gt_array[normal_index] = [2, 2]
        gt_array[tumor_index] = [2, 4]
        record.set_format("GT", gt_array)

        ### calculate AD and AF
        ### check if SNV
        if len(record.ALT[0]) == 1 and len(record.REF) == 1:
            ad_left = [record.format(f"{record.REF}U").flatten()[0], record.format(f"{record.ALT[0]}U").flatten()[0]]
            dp_left = np.sum(ad_left)
            af_left = record.format(f"{record.ALT[0]}U").flatten()[0] / dp_left

            ad_right = [record.format(f"{record.REF}U").flatten()[2], record.format(f"{record.ALT[0]}U").flatten()[2]]
            dp_right = np.sum(ad_right)
            af_right = record.format(f"{record.ALT[0]}U").flatten()[2] / dp_right
        else:
            ad_left = [record.format("TAR").flatten()[0], record.format("TIR").flatten()[0]]
            dp_left = np.sum(ad_left)
            af_left = record.format("TIR").flatten()[0] / dp_left

            ad_right = [record.format("TAR").flatten()[2], record.format("TIR").flatten()[2]]
            dp_right = np.sum(ad_right)
            af_right = record.format("TIR").flatten()[2] / dp_right


        ### set AF format
        af_array = np.array([[af_left], [af_right]])
        record.set_format("AF", af_array)

        ### set AD format
        ad_array = np.array([ad_left, ad_right])
        record.set_format("AD", ad_array)

        ### write record
        writer.write_record(record)

    ### close files
    vcf.close()
    writer.close()

    ### index modified vcf
    os.system(f"bcftools index {vcf_file_path[:-len(ending)]}.modified{ending}")

    return f"{vcf_file_path[:-len(ending)]}.modified{ending}"

def load_vcf_into_dataframe (vcf_file_path:str, mutect2_sample:str, strelka_sample:str) -> pd.DataFrame:
    """
    Load merged vcf short variant calls from Mutect2 and Strelka into a DataFrame, keep important metrics for downstream analysis

    Args:
        vcf_file_path (str): Path to merged VCF file
        mutect2_sample (str): Mutect2 tumor sample name
        strelka_sample (str): Strelka tumor sample name

    Raises:
        ValueError: If variant not found in any caller

    Returns:
        pd.DataFrame: DataFrame containing variant information
    """

    ### load merged VCF
    vcf = cyvcf2.VCF(vcf_file_path)

    ### get index of mutect2 and strelka tumor samples
    mutect2_idx = vcf.samples.index(mutect2_sample)
    strelka_idx = vcf.samples.index(strelka_sample)

    ### store data here
    data = []

    ### query all variants
    for variant in vcf:
    
        ### get mutect genotype
        gt = variant.genotypes[mutect2_idx]
        if gt[0] == -1 and gt[1] == -1:
            gt_mutect2 = False
        else:
            gt_mutect2 = True
        
        ### get strelka genotype
        gt = variant.genotypes[strelka_idx]
        if gt[0] == -1 and gt[1] == -1:
            gt_strelka = False
        else:
            gt_strelka = True

        ### get number of supporting callers
        # nc is number of callers finding this variant (can be PASS or FAIL)
        nc = np.sum([gt_mutect2, gt_strelka])
        # nc_filtered is number of callers finding this variant with PASS filter
        nc_filtered = np.sum(variant.format("FILTER") > 0)

        ### get variant metrics
        dp = -1
        ad = [-1, -1]
        af = -1.0

        ### use Mutect2 metrics if possible
        if gt_mutect2:
            dp = variant.format("DP")[mutect2_idx][0]
            ad = variant.format("AD")[mutect2_idx]
            af = variant.format("AF")[mutect2_idx][0]
        ### otherwise use Strelka metrics
        elif gt_strelka:
            dp = variant.format("DP")[strelka_idx][0]
            ad = variant.format("AD")[strelka_idx]
            af = variant.format("AF")[strelka_idx][0]
        else:
            raise ValueError("Variant not found in any caller")

        ### store record
        data.append({
            "CHROM":variant.CHROM,
            "POS":variant.POS,
            "REF":variant.REF,
            "ALT":variant.ALT[0],        
            "NC":nc,
            "NC_PASS":nc_filtered,
            "AF":af,
            "AD":ad,
            "DP":dp,
            "M2":np.sum(gt_mutect2),
            "ST":np.sum(gt_strelka),
            "M2_PASS":1 if gt_mutect2 and (variant.format("FILTER")[mutect2_idx] == 1) else 0,
            "ST_PASS":1 if gt_strelka and (variant.format("FILTER")[strelka_idx] == 1) else 0
        })
    
    ### return DataFrame
    return pd.DataFrame(data)

def write_clean_merged_vcf (variant_data:pd.DataFrame, output_vcf_path:str, vcf_header_path:str, sample_name:str):
    """
    Write cleaned merged VCF file with all variants from either Mutect2 or Strelka2. The information regarding PASS/FAIL is included in the FORMAT fields for each sample.

    Args:
        variant_data (pd.DataFrame): DataFrame with variant data. Output of vcf_methods.load_vcf_into_dataframe
        output_vcf_path (str): Path to output VCF file.
        vcf_header_path (str): Path to VCF header file to include in output VCF.
        sample_name (str): Name of the sample to include in the VCF.
    """

    ### open clean VCF file for writing
    with open(output_vcf_path, "w") as f:

        ### write header
        for line in open(vcf_header_path, "r"):
            f.write(line)
    
        ### write column names
        f.write(f"\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")

        ### write variants
        for _, row in variant_data.iterrows():

            format_string = ":".join(
                [
                    "0/1",
                    str(int(row["DP"])),
                    f"{int(row['AD'][0])},{int(row['AD'][1])}",
                    f"{row['AF']:.4f}",
                    str(int(row['NC'])),
                    str(int(row['NC_PASS'])),
                    str(int(row['M2'])),
                    str(int(row['ST'])),
                    str(int(row['M2_PASS'])),
                    str(int(row['ST_PASS']))
                ])
            
            f.write("\t".join(
                [
                    row["CHROM"],
                    str(row["POS"]),
                    ".",
                    row["REF"],
                    row["ALT"],
                    ".",
                    "PASS",
                    ".",
                    "GT:DP:AD:AF:NC:NC_PASS:M2:ST:M2_PASS:ST_PASS",
                    format_string
                ]
            ) + "\n")