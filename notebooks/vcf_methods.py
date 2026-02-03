import cyvcf2 # type: ignore
import numpy as np # type: ignore
import os


def modify_mutect2_vcf (vcf_file_path:str, normal_sample:str, tumor_sample:str):
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


def modify_strelka_vcf (vcf_file_path:str, normal_sample:str, tumor_sample:str):
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
        "ID": "FILTER",
        "Number": "1",
        "Type": "String",
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
        ad_left = [record.format(f"{record.REF}U").flatten()[0], record.format(f"{record.ALT[0]}U").flatten()[0]]
        dp_left = np.sum(ad_left)
        af_left = record.format(f"{record.ALT[0]}U").flatten()[0] / dp_left

        ad_right = [record.format(f"{record.REF}U").flatten()[2], record.format(f"{record.ALT[0]}U").flatten()[2]]
        dp_right = np.sum(ad_right)
        af_right = record.format(f"{record.ALT[0]}U").flatten()[2] / dp_right

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
