import argparse
import os
import pandas as pd
import zipfile

from collections import defaultdict
from CRISPResso2 import CRISPRessoShared

def reconstitute_reads(crispresso_output_folder, fastq_output_file):
    print('Processing ' + crispresso_output_folder)
    crispresso2_info = CRISPRessoShared.load_crispresso_info(crispresso_output_folder)

    total_alinged_reads = crispresso2_info['running_info']['alignment_stats']['N_CACHED_ALN'] + crispresso2_info['running_info']['alignment_stats']['N_COMPUTED_ALN']

    z = zipfile.ZipFile(os.path.join(crispresso_output_folder, crispresso2_info['running_info']['allele_frequency_table_zip_filename']))
    zf = z.open(crispresso2_info['running_info']['allele_frequency_table_filename'])
    df_alleles = pd.read_csv(zf,sep="\t") # use pandas to open because sometimes there's a quoted newline for long lines in the file that pandas writes??

    if df_alleles.shape[0] == 0:
        raise Exception('No alleles found in allele frequency table.')
    if df_alleles.columns[0] != '#Reads':
        raise Exception('Cannot parse allele frequency table, unexpected first column: ' + df_alleles.columns[0])
    if df_alleles.columns[1] != 'Aligned_Sequence':
        raise Exception('Cannot parse allele frequency table, unexpected second column: ' + df_alleles.columns[1])

    with open(fastq_output_file, "wt") as fout:
        allele_count = 0
        read_count = 0
        for _, allele_row in df_alleles.iterrows():
            seq_count = allele_row['#Reads']

            for i in range(seq_count):
                seq = allele_row['Aligned_Sequence'].replace('-', '')
                fout.write(f"@SEQ_ID_{allele_count}_{i}\n")
                fout.write(f"{seq}\n")
                fout.write("+\n")
                fout.write('A' * len(seq) + '\n')
                read_count += 1
            allele_count += 1
    print(f"Wrote {allele_count} alleles in {read_count} reads to {fastq_output_file}")
    if read_count != total_alinged_reads:
        print(f"Warning: total reads written ({read_count}) does not match total aligned reads ({total_alinged_reads})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Reconstitute reads from CRISPResso2 allele frequency table, producing a fastq file for reanalysis.')
    parser.add_argument('crispresso_output_folder', type=str, help='Path to the CRISPResso2 output folder.')
    parser.add_argument('fastq_output_file', type=str, help='Path to the output FASTQ file.')
    args = parser.parse_args()

    reconstitute_reads(args.crispresso_output_folder, args.fastq_output_file)
