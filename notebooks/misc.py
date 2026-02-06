def parse_format (format_string:str) -> tuple[int, int, int]:
    """
    Parses format from harmoized vcf files

    Args:
        format_string (str): AD:PASS

    Raises:
        ValueError: Unexpected format string in VCF

    Returns:
        tuple[int, int, int]: AD_REF, AD_ALT, PASS
    """
    ### return -1 for empty calls
    if format_string == ".:.":
        return -1, -1, -1
    ### return -1 for multiallelic calls
    if len(format_string.split(":")[0].split(",")) > 2:
        return -1, -1, -1
    ad = format_string.split(":")[0]
    pass_value = format_string.split(":")[1]
    try:
        return int(ad.split(",")[0]), int(ad.split(",")[1]), int(pass_value)
    except ValueError:
        print(format_string)
        raise ValueError("Unexpected format string in VCF")