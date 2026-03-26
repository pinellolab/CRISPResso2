import gzip
import argparse


def get_dict_from_plus_line(plus_line):
    """Create a dictionary from the 3rd line (plus line) of a CRISPResso-annotated fastq file

    Args:
        plus_line (str): The 3rd line of a CRISPResso-annotated fastq file

    Returns:
        dict: dictionary of CRISPResso annotations > values

    """
    equals_line_els = plus_line.split(" ")
    equals_dict = {}
    for el in equals_line_els:
        el_els = el.split("=")
        if len(el_els) > 1:
            equals_dict[el_els[0]] = el_els[1]
        else:
            equals_dict[el_els[0]] = None
    return equals_dict


def get_ref_details_from_aln_details(aln_details):
    """Get the reference details from the 'ALN_DETAILS' element of the 3rd line (plus line) of a CRISPResso-annotated fastq file

    Args:
        aln_details (str): The ALN_DETAILS element consisting of the reference name, read sequence alignment, reference sequence alignment, and alignment score for each reference

    Returns:
        dict: dictionary of ref_name > read_seq_alignment, ref_seq_alignment, aln_score

    """
    ref_details = {}
    ref_els = aln_details.split("&")
    for el in ref_els:
        el_els = el.split(",")
        if len(el_els) != 4:
            raise Exception("got unexpected number of elements in ALN_DETAILS: " + str(el_els))
        (ref_name, read_seq, ref_seq, aln_score) = el.split(",")
        ref_details[ref_name] = (read_seq, ref_seq, aln_score)
    return ref_details


def get_ref_coordinates_from_aln(read_aln, ref_aln):
    """Get the reference coordinates from the read and reference alignments
    So if we want to check the base at the 5th base in the reference, we would go to the ref_pos_in_aln[5]th position in the alignment
    If we want to get the base in the read corresponding to the 5th base in the reference, we would go to the ref_pos_in_read[5]th position in the read

    Args:
        read_aln (str): read alignment
        ref_aln (str): reference alignment

    Returns:
        arr(int): array of reference positions in the alignment
        arr(int): array of reference positions in the read

    """
    read_pos = -1
    ref_pos_in_aln = []  # list of reference positions in the alignment
    ref_pos_in_read = []  # list of reference positions in the read
    for i in range(len(read_aln)):
        if read_aln[i] != "-":
            read_pos += 1
        if ref_aln[i] != "-":
            ref_pos_in_aln.append(i)
            ref_pos_in_read.append(read_pos)

    return (ref_pos_in_aln, ref_pos_in_read)


def count_high_quality_bases(file, ref_name, base_pos_in_ref, original_base, target_base, quality_cutoff):
    """Counts high-quality bases at a certain position in the reference

    Args:
        file (str): CRISPResso_output.fastq.gz file to analyze
        ref (str): Reference name
        base_pos_in_ref (int): 0-based position in the reference to check
        original_base (char): original base to check
        target_base (char): target base to check
        quality_cutoff (int): quality cutoff for high-quality bases

    """
    # file = "CRISPResso_output.fastq.gz"
    # ref_name = "Reference"
    # base_pos_in_ref = 100
    # original_base = "A"
    # target_base = "T"
    # quality_cutoff = 30

    count_total_reads = 0
    count_unaligned = 0  # also includes reads not aligned to specified ref_name

    count_original_base_highQual = 0
    count_original_base_lowQual = 0

    count_target_base_highQual = 0
    count_target_base_lowQual = 0

    count_base_deletion = 0

    count_other_base_highQual = 0
    count_other_base_lowQual = 0

    with gzip.open(file, 'rt') as f:
        while True:
            id_line = f.readline()
            seq_line = f.readline()
            plus_line = f.readline()
            qual_line = f.readline()

            if not id_line:
                break

            count_total_reads += 1

            plus_line_dict = get_dict_from_plus_line(plus_line)
            if 'ALN' not in plus_line_dict:
                raise Exception("ALN not in plus line: " + plus_line)
            if plus_line_dict['ALN'] != ref_name:
                count_unaligned += 1
            else:
                ref_details = get_ref_details_from_aln_details(plus_line_dict['ALN_DETAILS'])
                if ref_name in ref_details:
                    read_aln = ref_details[ref_name][0]  # sequence of read aligned to ref (includes gaps)
                    ref_aln = ref_details[ref_name][1]  # sequence of ref aligned to read (includes gaps)
                    ref_coordinates_in_aln, ref_coordinates_in_read = get_ref_coordinates_from_aln(read_aln, ref_aln)
                    base_in_read = read_aln[ref_coordinates_in_aln[base_pos_in_ref]]
                    quality_in_read = qual_line[ref_coordinates_in_read[base_pos_in_ref]]

                    if base_in_read == original_base:
                        if ord(quality_in_read) - 33 >= quality_cutoff:
                            count_original_base_highQual += 1
                        else:
                            count_original_base_lowQual += 1
                    elif base_in_read == target_base:
                        if ord(quality_in_read) - 33 >= quality_cutoff:
                            count_target_base_highQual += 1
                        else:
                            count_target_base_lowQual += 1
                    elif base_in_read == "-":
                        count_base_deletion += 1
                    elif ord(quality_in_read) - 33 >= quality_cutoff:
                        count_other_base_highQual += 1
                    else:
                        count_other_base_lowQual += 1

    print("Total reads read: " + str(count_total_reads))
    print("Unaligned reads: " + str(count_unaligned))
    print("Original base high quality: " + str(count_original_base_highQual))
    print("Original base low quality: " + str(count_original_base_lowQual))
    print("Target base high quality: " + str(count_target_base_highQual))
    print("Target base low quality: " + str(count_target_base_lowQual))
    print("Base deletion: " + str(count_base_deletion))
    print("Other base high quality: " + str(count_other_base_highQual))
    print("Other base low quality: " + str(count_other_base_lowQual))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="CRISPResso_output.fastq.gz file to analyze", required=True)
    parser.add_argument("-r", "--ref_name", help="Reference name", default="Reference")
    parser.add_argument("-p", "--pos", type=int, help="0-based position in the reference to check", required=True)
    parser.add_argument("-o", "--original", help="original base to check", required=True)
    parser.add_argument("-t", "--target", help="target base to check", required=True)
    parser.add_argument("-q", "--quality", type=int, help="quality cutoff for high-quality bases", default=30)
    args = parser.parse_args()

    count_high_quality_bases(args.file, args.ref_name, args.pos, args.original, args.target, args.quality)
