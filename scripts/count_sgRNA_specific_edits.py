import argparse
import os
import pandas as pd
import zipfile

from collections import defaultdict
from CRISPResso2 import CRISPRessoShared

def count_sgRNA_specific_edits(crispresso_output_folder):
    crispresso2_info = CRISPRessoShared.load_crispresso_info(crispresso_output_folder)

    run_version = crispresso2_info['running_info']['version']
    version_parts = run_version.split('.')
    if int(version_parts[0]) != 2 or int(version_parts[1]) < 2 or int(version_parts[2]) < 15:
        raise Exception('CRISPResso run must be run with CRISPResso2 v2.2.15 or later (this run was run with version v' + str(run_version) + ')')

    if not crispresso2_info['running_info']['args'].write_detailed_allele_table:
        raise Exception('CRISPResso run must be run with the parameter --write_detailed_allele_table')

    include_idxs = {} # ref_name > guide_name > idxs
    all_reference_guides = {} # order of guides for each reference
    modified_counts_by_amplicon = {} # ref_name > category > count
    for reference_name in crispresso2_info['results']['ref_names']:
        include_idxs[reference_name] = {}
        all_reference_guides[reference_name] = []
        modified_counts_by_amplicon[reference_name] = defaultdict(int)
        for idx, guide_name in enumerate(crispresso2_info['results']['refs'][reference_name]['sgRNA_names']):
            this_guide_name = guide_name
            if this_guide_name == '':
                this_guide_name = crispresso2_info['results']['refs'][reference_name]['sgRNA_orig_sequences'][idx]
            if this_guide_name in include_idxs[reference_name]:
                this_guide_name = this_guide_name + '_' + str(crispresso2_info['results']['refs'][reference_name]['sgRNA_intervals'][idx][0])
            include_idxs[reference_name][this_guide_name] = crispresso2_info['results']['refs'][reference_name]['sgRNA_include_idxs'][idx]
            all_reference_guides[reference_name].append(this_guide_name)


    z = zipfile.ZipFile(os.path.join(crispresso_output_folder, crispresso2_info['running_info']['allele_frequency_table_zip_filename']))
    zf = z.open(crispresso2_info['running_info']['allele_frequency_table_filename'])
    df_alleles = pd.read_csv(zf,sep="\t")

    read_alleles_count = 0 # all alleles read
    reference_count = defaultdict(int) # counts of modification at reference level
    reference_modified_count = defaultdict(int)
    total_counts = defaultdict(int)
    modified_counts = defaultdict(int)
    for idx, row in df_alleles.iterrows():
        read_alleles_count += row['#Reads']
        this_reference_name = row['Reference_Name']
        if this_reference_name in include_idxs: # sometimes (e.g. AMBIGUOUS) reads aren't assigned to a reference, so exclude them here..
            reference_count[this_reference_name] += row['#Reads']

            this_allele_modified_guide_names = []
            ref_is_modified = False
            row_insertion_positions = [int(x) for x in row['insertion_positions'][1:-1].split(",")] if row['insertion_positions'] != '[]' else []
            row_deletion_positions = [int(x) for x in row['deletion_positions'][1:-1].split(",")] if row['deletion_positions'] != '[]' else []
            row_substitution_positions = [int(x) for x in row['substitution_positions'][1:-1].split(",")] if row['substitution_positions'] != '[]' else []
            for guide_name in all_reference_guides[this_reference_name]:
                is_modified = False
                for ind in include_idxs[this_reference_name][guide_name]:
                    if ind in row_insertion_positions:
                        is_modified = True
                    if ind in row_deletion_positions:
                        is_modified = True
                    if ind in row_substitution_positions:
                        is_modified = True
                if is_modified:
                    this_allele_modified_guide_names.append(guide_name)
                    ref_is_modified = True
            
            this_category = 'UNMODIFIED'
            if ref_is_modified:
                reference_modified_count[this_reference_name] += row['#Reads']
                this_category = 'MODIFIED ' + ' + '.join(this_allele_modified_guide_names)
            modified_counts_by_amplicon[this_reference_name][this_category] += row['#Reads']


    print('Processed ' + str(read_alleles_count) + ' alleles')
    for reference_name in crispresso2_info['results']['ref_names']:
        print ('Reference: ' + reference_name + ' (' + str(reference_modified_count[reference_name]) + '/' + str(reference_count[reference_name]) + ' modified reads)')
        for category in modified_counts_by_amplicon[reference_name]:
            print("\t" + category + ": " + str(modified_counts_by_amplicon[reference_name][category]))



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Count sgRNA-specific edits')
    parser.add_argument("-f", "--folder", help="CRISPResso output folder", required=True)
    args = parser.parse_args()

    count_sgRNA_specific_edits(args.folder)
    