set -e
echo Running CRISPResso
CRISPResso -r1 FANC.Cas9.fastq -a CGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAGCTTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGGACCCCGCCACCGTGCGCCGGGCCTTGCAGTGGGCGCGCTACCTGCGCCACATCCATCGGCGCTTTGGTCGG -g GGAATCCCTTCTGCAGCACC --debug &> CRISPResso_on_FANC.Cas9.log

echo Running CRISPResso with parameters
CRISPResso -r1 FANC.Cas9.fastq -a CGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAGCTTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGGACCCCGCCACCGTGCGCCGGGCCTTGCAGTGGGCGCGCTACCTGCGCCACATCCATCGGCGCTTTGGTCGG -g GGAATCCCTTCTGCAGCACC -e CGGCCGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCTGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGCTTTTCCGAGCTTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGGACCCCGCCACCGTGCGCCGGGCCTTGCAGTGGGCGCGCTACCTGCGCCACATCCATCGGCGCTTTGGTCGG -c GGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGCACCTGGATCGCTTTT --dump -qwc 20-30_45-50 -q 30 --default_min_aln_score 80 -an FANC -n params --base_edit -fg AGCCTTGCAGTGGGCGCGCTA,CCCACTGAAGGCCC --dsODN GCTAGATTTCCCAAGAAGA -gn hi -fgn dear --debug &> CRISPResso_on_params.log

echo Running CRISPRessoBatch
CRISPRessoBatch -bs FANC.local.batch -a CGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAGCTTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGGACCCCGCCACCGTGCGCCGGGCCTTGCAGTGGGCGCGCTACCTGCGCCACATCCATCGGCGCTTTGGTCGG -g GGAATCCCTTCTGCAGCACC -p 2 --debug --base_editor -n FANC --debug &> CRISPRessoBatch_on_FANC.log

echo Running CRISPRessoPooled
CRISPRessoPooled -r1 Both.Cas9.fastq -f Cas9.amplicons.txt -p 2 --keep_intermediate --min_reads_to_use_region 100 --debug &> CRISPRessoPooled_on_Both.Cas9.log

echo Running CRISPRessoWGS
CRISPRessoWGS -b Both.Cas9.fastq.smallGenome.bam -r smallGenome/smallGenome.fa -f Cas9.regions.txt --debug &> CRISPRessoWGS_on_Both.Cas9.fastq.smallGenome.log

echo TESTING CRISPRESSO2
diff CRISPResso_on_FANC.Cas9/Nucleotide_frequency_table.txt expectedResults/CRISPResso_on_FANC.Cas9/Nucleotide_frequency_table.txt
diff CRISPResso_on_FANC.Cas9/CRISPResso_quantification_of_editing_frequency.txt expectedResults/CRISPResso_on_FANC.Cas9/CRISPResso_quantification_of_editing_frequency.txt

echo TESTING CRISPRESSO2 PARAMS
diff CRISPResso_on_params/FANC.Nucleotide_frequency_table.txt expectedResults/CRISPResso_on_params/FANC.Nucleotide_frequency_table.txt
diff CRISPResso_on_params/CRISPResso_quantification_of_editing_frequency.txt expectedResults/CRISPResso_on_params/CRISPResso_quantification_of_editing_frequency.txt

echo TESTING BATCH
diff CRISPRessoBatch_on_FANC/MODIFICATION_FREQUENCY_SUMMARY.txt expectedResults/CRISPRessoBatch_on_FANC/MODIFICATION_FREQUENCY_SUMMARY.txt

echo TESTING POOLED
diff CRISPRessoPooled_on_Both.Cas9/SAMPLES_QUANTIFICATION_SUMMARY.txt expectedResults/CRISPRessoPooled_on_Both.Cas9/SAMPLES_QUANTIFICATION_SUMMARY.txt

echo TESTING WGS
diff CRISPRessoWGS_on_Both.Cas9.fastq.smallGenome/SAMPLES_QUANTIFICATION_SUMMARY.txt expectedResults/CRISPRessoWGS_on_Both.Cas9.fastq.smallGenome/SAMPLES_QUANTIFICATION_SUMMARY.txt

echo Finished
