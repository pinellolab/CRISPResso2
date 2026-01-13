# Changelog

## Unreleased
### ADDED

- Add an amino acid nucleotide quilt plot by [@mbowcut2](https://github.com/mbowcut2) in [#552](https://github.com/pinellolab/CRISPResso2/pull/552)

- Add `scripts/reconstituteReads.py` to generate FASTQ from CRISPResso2 output by [@kclem](https://github.com/kclem) in [`a800762`](https://github.com/pinellolab/CRISPResso2/commit/a800762712e692a9fc7005b3510d013924c843a5) and [`cd79dcc`](https://github.com/pinellolab/CRISPResso2/commit/cd79dcc3ca9c06736dc6cd409d43dee8e4d1ae69)

- Add an [UpSet plot](https://upset.app/) to represent bystander edits for Base Editing analyses by [@mbowcut2](https://github.com/mbowcut2) in [#554](https://github.com/pinellolab/CRISPResso2/pull/554)

- Allow for messages to be served via CRISPResso reports by [@Colelyman](https://github.com/Colelyman) in [#583](https://github.com/pinellolab/CRISPResso2/pull/583)

- Add a plot that shows the distribution of homology scores for reads by [@mbowcut2](https://github.com/mbowcut2) in [#600](https://github.com/pinellolab/CRISPResso2/pull/600)

### FIXED

- Fix the `x_lim` settings on plot 3b by [@kclem](https://github.com/kclem) in [`56bd430`](https://github.com/pinellolab/CRISPResso2/commit/56bd4306292136ed95d7032b0582c6ad370dd79b)

- Fix parsing the `CRISPResso2_info.json` in CRISPRessoPooled by [@kclem](https://github.com/kclem) in [#558](https://github.com/pinellolab/CRISPResso2/pull/558)

- Forced cloned `include_idxs` to be `np.array`s by [@kclem](https://github.com/kclem) in [`da4badb`](https://github.com/pinellolab/CRISPResso2/commit/da4badb34282928fe62cbe56dabf8b0d9b1acdee)

- Fix the link to the CRISPResso cup in reports (so that SSL works correctly) by [@Colelyman](https://github.com/Colelyman) in [#571](https://github.com/pinellolab/CRISPResso2/pull/571)

- Fix the quantification of deletions at the second position of the sequence by [@Colelyman](https://github.com/Colelyman) in [#574](https://github.com/pinellolab/CRISPResso2/pull/574)

- Fix an issue with unaligned reads not being reported correctly when writing BAM output by [@trevormartinj7](https://github.com/trevormartinj7) in [#578](https://github.com/pinellolab/CRISPResso2/pull/578)

- Fix an issue where quantification window coordinates we not being correctly inferred by [@Colelyman](https://github.com/Colelyman) in [#598](https://github.com/pinellolab/CRISPResso2/pull/598)
  - This issue is present when there is a single quantifcation window coordinate provided and multiple amplicons. What happens is CRISPResso aligns the second amplicon to the first and then infers what the quantification window coordinates should be based on the alignment. A regression was introduced where the inference of the quantification window coordinates for the second amplicon was no longer correct. This change fixes the regression and brings the behavior back to match that of v2.2.9.
  - If you don't set quantification window coordinates and don't use multiple amplicons, there is no need for this fix and therefore no change in behavior.

- Fix a `SyntaxWarning` for an unescaped sequence in a matplotlib function by [@Colelyman](https://github.com/Colelyman) in [#600](https://github.com/pinellolab/CRISPResso2/pull/600)

- Fix a bug during `--bam_output` when there is an unaligned read, the remainder of the reads will not bu output by [@Colelyman](https://github.com/Colelyman) in [#602](https://github.com/pinellolab/CRISPResso2/pull/602)

### CHANGED

- Update the base Docker image to `mambaorg/micromamba:2.3.3` and remove dependency on Anaconda `defaults` channel by [@Colelyman](https://github.com/Colelyman) in [#575](https://github.com/pinellolab/CRISPResso2/pull/575)

### REMOVED


## v2.3.3 - Activity Fulton - 07/01/2025
### ADDED
- Asymmetrical Allele Plots by [@mbowcut2](https://github.com/mbowcut2) and [@Colelyman](https://github.com/Colelyman) in [#527](https://github.com/pinellolab/CRISPResso2/pull/527)

- Native CRISPResso Paired End Read Merging by [@Colelyman](https://github.com/Colelyman) and [@Snicker7](https://github.com/Snicker7) in [#537](https://github.com/pinellolab/CRISPResso2/pull/537)

- Support for amplicon names with emojis 🎉🤩 and other non-standard characters by [@kclem](https://github.com/kclem) in [`1455609`](https://github.com/pinellolab/CRISPResso2/commit/14556097fe23ef114071a7b3361c5a5c59b2919b)

- Multiplexing for CRISPRessoPooled subruns in [`1455609`](https://github.com/pinellolab/CRISPResso2/commit/14556097fe23ef114071a7b3361c5a5c59b2919b)

### FIXED
- Fix setting of 99%ile in negative direction for deletions plot in [`90ac42a`](https://github.com/pinellolab/CRISPResso2/commit/90ac42a7a74b88e1e249e72b5e87c9c3db21e1c9)

### CHANGED
- Make fig_filename_root default to None, in which case the figure is shown interactively (e.g. in a jupyter notebook) in [`c2a10c4`](https://github.com/pinellolab/CRISPResso2/commit/c2a10c4fd296e49a987176d2af7ce1a7507465c2)

- Don't rerun if --no_rerun is set but --verbosity has changed in [`6562a08`](https://github.com/pinellolab/CRISPResso2/commit/6562a08ac69b1ff7928c9cb3cb325466a4f17f07)

### REMOVED
- Remove warning for zipping nonexistant files in [`3784ea5`](https://github.com/pinellolab/CRISPResso2/commit/3784ea5e4bd35190c83b838f20b0bbd4955fee39)


## v2.3.2 - Junction Salt - 01/16/2025
### ADDED
- New parameters, `--flexiguide_gap_open_penalty` and `--flexiguide_gap_extend_penalty`, to customize flexiguide alignment in [#491](https://github.com/pinellolab/CRISPResso2/pull/491)

- New parameter `--halt_on_plot_fail` so that errors and exceptions in plots don't fail silently in [#494](https://github.com/pinellolab/CRISPResso2/pull/494)

- New parameter `--samtools_exclude_flag` to customize the filtering of reads in [#503](https://github.com/pinellolab/CRISPResso2/pull/503)

- New documentation website at <docs.crispresso.com>.

- Asymmetrical cut point by [@kclem](https://github.com/kclem) in [#457](https://github.com/pinellolab/CRISPResso2/pull/457)

- d3 plot enhancements by [@trevormartinj7](https://github.com/trevormartinj7) in [#459](https://github.com/pinellolab/CRISPResso2/pull/459)

- Add flexiguide alignment parameters by [@Colelyman](https://github.com/Colelyman) in [#491](https://github.com/pinellolab/CRISPResso2/pull/491)

- Add pyproject.toml and support numpy v2 by [@Snicker7](https://github.com/Snicker7) in [#496](https://github.com/pinellolab/CRISPResso2/pull/496)

- Add customizable samtools exclude flag by [@Colelyman](https://github.com/Colelyman) in [#503](https://github.com/pinellolab/CRISPResso2/pull/503)

- Add support for octal and comma separated samtools exclude flags ([#113](https://github.com/pinellolab/CRISPResso2/pull/113)) by [@kclem](https://github.com/kclem) in [#507](https://github.com/pinellolab/CRISPResso2/pull/507)

### FIXED
- Fix typo and move flexiguide to debug ([#77](https://github.com/pinellolab/CRISPResso2/issues/77)) by [@Colelyman](https://github.com/Colelyman) in [#438](https://github.com/pinellolab/CRISPResso2/pull/438)

- Matplotlib Compatibility Fix by [@mbowcut2](https://github.com/mbowcut2) and [@Snicker7](https://github.com/Snicker7) in [#464](https://github.com/pinellolab/CRISPResso2/pull/464)

- Fix CRISPRessoAggregate bug and other improvements ([#95](https://github.com/pinellolab/CRISPResso2/issues/95)) by [@Colelyman](https://github.com/Colelyman) in [#470](https://github.com/pinellolab/CRISPResso2/pull/470)

- Fix missing substitution in name of WGS, Compare and Meta reports by [@Colelyman](https://github.com/Colelyman) in [#498](https://github.com/pinellolab/CRISPResso2/pull/498)

- Fix `get_n_fastq` function by [@trevormartinj7](https://github.com/trevormartinj7) in [#508](https://github.com/pinellolab/CRISPResso2/pull/508)

### CHANGED
- Improvement of processing speed by
    - parallelization of alignments in [#480](https://github.com/pinellolab/CRISPResso2/pull/480)

    - memory usage reduction in [#478](https://github.com/pinellolab/CRISPResso2/pull/478) and [#509](https://github.com/pinellolab/CRISPResso2/pull/509)

- Progress percentages are displayed in the CLI output.

- Prefix the release Docker tag with a `v` by [@Colelyman](https://github.com/Colelyman) in [#434](https://github.com/pinellolab/CRISPResso2/pull/434)

- Pin versions of numpy and matplotlib in CI environment by [@Snicker7](https://github.com/Snicker7) in [#452](https://github.com/pinellolab/CRISPResso2/pull/452)

- Implement new pooled mixed-mode default behavior by [@mbowcut2](https://github.com/mbowcut2) in [#454](https://github.com/pinellolab/CRISPResso2/pull/454)

- Update README by [@Snicker7](https://github.com/Snicker7), [@mbowcut2](https://github.com/mbowcut2), [@trevormartinj7](https://github.com/trevormartinj7), [@Colelyman](https://github.com/Colelyman) and [@kclem](https://github.com/kclem) in [#456](https://github.com/pinellolab/CRISPResso2/pull/456)

- Cache conda packages in GIthub Actions by [@Colelyman](https://github.com/Colelyman) in [#466](https://github.com/pinellolab/CRISPResso2/pull/466)

- Replace zcat by [@Colelyman](https://github.com/Colelyman) in [#468](https://github.com/pinellolab/CRISPResso2/pull/468)

- Cache read merging step in CRISPRessoPooled on no_rerun by [@kclem](https://github.com/kclem) in [#467](https://github.com/pinellolab/CRISPResso2/pull/467)

- Display percentages in the CLI output by [@Colelyman](https://github.com/Colelyman) in [#473](https://github.com/pinellolab/CRISPResso2/pull/473)

- No processor pool when running in single thread by [@Snicker7](https://github.com/Snicker7) in [#474](https://github.com/pinellolab/CRISPResso2/pull/474)

- Round percentage complete in CLI and add initial 0% complete by [@Colelyman](https://github.com/Colelyman) in [#477](https://github.com/pinellolab/CRISPResso2/pull/477)

- Reduce memory usage for allele plots by [@Colelyman](https://github.com/Colelyman) in [#478](https://github.com/pinellolab/CRISPResso2/pull/478)

- Sync reports by [@Colelyman](https://github.com/Colelyman) in [#479](https://github.com/pinellolab/CRISPResso2/pull/479)

- Read Alignment Parallelization (#98) by [@trevormartinj7](https://github.com/trevormartinj7) in [#480](https://github.com/pinellolab/CRISPResso2/pull/480)

- Add `all_deletion_coordinates` to be returned by `find_indels_substitutions_legacy` function by [@Colelyman](https://github.com/Colelyman) in [#486](https://github.com/pinellolab/CRISPResso2/pull/486)

- Mckay/halt on plot fail by [@mbowcut2](https://github.com/mbowcut2) in [#494](https://github.com/pinellolab/CRISPResso2/pull/494)

- Update jinja_partials and bring Reports into sync by [@Colelyman](https://github.com/Colelyman) in [#500](https://github.com/pinellolab/CRISPResso2/pull/500)

- Update CRISPRessoPooledCORE.py by [@kclem](https://github.com/kclem) in [#502](https://github.com/pinellolab/CRISPResso2/pull/502)

- Pooled multi region map fix by [@kclem](https://github.com/kclem) in [#505](https://github.com/pinellolab/CRISPResso2/pull/505)

- Slots implementation by [@mbowcut2](https://github.com/mbowcut2) in [#509](https://github.com/pinellolab/CRISPResso2/pull/509)

- Suppress printing by [@kclem](https://github.com/kclem) in [#511](https://github.com/pinellolab/CRISPResso2/pull/511)

- Update detailed alleles table help option by [@Colelyman](https://github.com/Colelyman) in [#513](https://github.com/pinellolab/CRISPResso2/pull/513)


## v2.3.1 - Screen King - 05/13/2024
### FIXED
- Fix CRISPRessoPooled error reporting by [@kclem](https://github.com/kclem) in [#423](https://github.com/pinellolab/CRISPResso2/pull/423)

- Extract `jinja_partials` and fix CRISPRessoPooled fastp errors by [@Colelyman](https://github.com/Colelyman) and [@trevormartinj7](https://github.com/trevormartinj7) in [#425](https://github.com/pinellolab/CRISPResso2/pull/425)

- Fix batch mode pandas warning. ([#70](https://github.com/pinellolab/CRISPResso2/issues/70)) by [@mbowcut2](https://github.com/mbowcut2) and [@Colelyman](https://github.com/Colelyman) in [#429](https://github.com/pinellolab/CRISPResso2/pull/429)

- Fix issues with `file_prefix` by [@Colelyman](https://github.com/Colelyman) and [@Snicker7](https://github.com/Snicker7) in [#430](https://github.com/pinellolab/CRISPResso2/pull/430)

- Fix plots and improve plot error handling by [@Snicker7](https://github.com/Snicker7) and [@mbowcut2](https://github.com/mbowcut2) in [#431](https://github.com/pinellolab/CRISPResso2/pull/431)

### CHANGED
- Cole/refactor jinja undefined ([#66](https://github.com/pinellolab/CRISPResso2/issues/66)) by [@Colelyman](https://github.com/Colelyman) and @Skicker7 in [#421](https://github.com/pinellolab/CRISPResso2/pull/421)

- Update README by [@Colelyman](https://github.com/Colelyman) in [#424](https://github.com/pinellolab/CRISPResso2/pull/424)

- Bump version to 2.3.1 and change default CRISPRessoPooled behavior to change in 2.3.2 by [@Colelyman](https://github.com/Colelyman) in [#428](https://github.com/pinellolab/CRISPResso2/pull/428)

- Showing sgRNA sequences on hover in CRISPRessoPro by [@Colelyman](https://github.com/Colelyman) in [#432](https://github.com/pinellolab/CRISPResso2/pull/432)

### REMOVED
- Remove extra imports from CRISPRessoCore by [@Colelyman](https://github.com/Colelyman) in [#422](https://github.com/pinellolab/CRISPResso2/pull/422)


## v2.3.0 - Targeting Minato - 04/10/2024
### ADDED
- Guardrails (checking experimental conditions and raising warnings) by [@Snicker7](https://github.com/Snicker7)

- Enable quantification by sgRNA by [@kclem](https://github.com/kclem) in [#348](https://github.com/pinellolab/CRISPResso2/pull/348)

### FIXED
- Fix samtools piping by [@Colelyman](https://github.com/Colelyman) in [#325](https://github.com/pinellolab/CRISPResso2/pull/325)

- Fix for recent Matplotlib v3.8 by [@mbowcut2](https://github.com/mbowcut2) in [#346](https://github.com/pinellolab/CRISPResso2/pull/346)

- Fix Matplotlib breaking change issue by [@mbowcut2](https://github.com/mbowcut2) in [#352](https://github.com/pinellolab/CRISPResso2/pull/352)

- Fix assigning multiple qwc by [@Snicker7](https://github.com/Snicker7) in [#375](https://github.com/pinellolab/CRISPResso2/pull/375)

- Fix interleaved fastq input in CRISPRessoPooled and suppress CRISPRessoWGS params by [@Colelyman](https://github.com/Colelyman) in [#392](https://github.com/pinellolab/CRISPResso2/pull/392)

- Fix [#367](https://github.com/pinellolab/CRISPResso2/issues/367), reads only align to prime edited amplicon, not to reference by [@mbowcut2](https://github.com/mbowcut2) in [#393](https://github.com/pinellolab/CRISPResso2/pull/393)

- Fix the assignment of multiple quantification window coordinates by [@Snicker7](https://github.com/Snicker7) in [#403](https://github.com/pinellolab/CRISPResso2/pull/403)

- Fix `space` character in README by [@DennyDai](https://github.com/DennyDai) in [#400](https://github.com/pinellolab/CRISPResso2/pull/400)

- Fix Jinja2 undefined variables by [@Colelyman](https://github.com/Colelyman) in [#417](https://github.com/pinellolab/CRISPResso2/pull/417)

### CHANGED
- Flash and Trimmomatic are replaced with Fastp by [@trevormartinj7](https://github.com/trevormartinj7), [@Snicker7](https://github.com/Snicker7), and [@Colelyman](https://github.com/Colelyman)

- Failed runs are displayed with the error by [@trevormartinj7](https://github.com/trevormartinj7)

- Replace link to CRISPResso schematic with raw URL in README by [@Colelyman](https://github.com/Colelyman) in [#329](https://github.com/pinellolab/CRISPResso2/pull/329)

- Prime editing alignment params by [@kclem](https://github.com/kclem) in [#336](https://github.com/pinellolab/CRISPResso2/pull/336)

- Run unit tests via Github Actions and fix matplotlib character issue by [@Snicker7](https://github.com/Snicker7) and [@mbowcut2](https://github.com/mbowcut2) in [#386](https://github.com/pinellolab/CRISPResso2/pull/386)

- Remove future Pandas warnings and sort CRISPRessoCompare tables by [@mbowcut2](https://github.com/mbowcut2) and [@Snicker7](https://github.com/Snicker7) in [#389](https://github.com/pinellolab/CRISPResso2/pull/389)

- Run integration tests on every push by [@Snicker7](https://github.com/Snicker7) in [#394](https://github.com/pinellolab/CRISPResso2/pull/394)

- Move read filtering to after merging in CRISPResso by [@Colelyman](https://github.com/Colelyman) in [#397](https://github.com/pinellolab/CRISPResso2/pull/397)

- Decrease Docker image size and fix PE naming and parameter behavior by [@Colelyman](https://github.com/Colelyman), [@Snicker7](https://github.com/Snicker7) and [@mbowcut2](https://github.com/mbowcut2) in [#404](https://github.com/pinellolab/CRISPResso2/pull/404)


## v2.2.14 - Specific São Paulo - 08/10/2023
### FIXED
- [Fix missing ' in CRISPRessoPooled --demultiplex_only_at_amplicons](https://github.com/pinellolab/CRISPResso2/commit/3e04d1d402bcf95edd39fc7c8c9af61bb380f9db)

- [CRISPRessoPooled --compile_postrun_references bug fixes](https://github.com/pinellolab/CRISPResso2/commit/2b36a1a5c35e8a93516ce8baf464595615e0f402)

- [Fix bug in prime-editing scaffold-incorporation plotting](https://github.com/pinellolab/CRISPResso2/commit/4d9c71ecf2248c9bb1e10430178dc318b6621c8b)


## v2.2.13 - With Montgomery - 07/28/2023
### ADDED
- Add verbosity argument to CRISPRessoAggregate ([#18](https://github.com/pinellolab/CRISPResso2/issues/18)) fixes [#306](https://github.com/pinellolab/CRISPResso2/issues/306) by [@Colelyman](https://github.com/Colelyman) in [#307](https://github.com/pinellolab/CRISPResso2/pull/307)

### FIXED
- Parallel plotting fix by [@Colelyman](https://github.com/Colelyman) and [@kclem](https://github.com/kclem) in [`546446e`](https://github.com/pinellolab/CRISPResso2/commit/546446e36e7e68b527767d6c31ec341a49df2059) and [#286](https://github.com/pinellolab/CRISPResso2/pull/286)

- Fix multiprocessing lambda pickling by [@Colelyman](https://github.com/Colelyman) in [#311](https://github.com/pinellolab/CRISPResso2/pull/311)

### CHANGED
- Don't start pool when only using single thread by [@Colelyman](https://github.com/Colelyman) in [#302](https://github.com/pinellolab/CRISPResso2/pull/302)

- Raise exceptions from within futures in plot_pool in [`a439f09`](https://github.com/pinellolab/CRISPResso2/commit/a439f094745b2b5e7f032f0777d4c67e6d6f93c5)

- Enable CRISPRessoPooled multiprocessing when os allows multi-thread file append in [`ebb016d`](https://github.com/pinellolab/CRISPResso2/commit/ebb016dff46c280dce8c3c09e8ac0e0cc25d4d74)

- Allow multiple overlapping sgRNA matches in reference (previous behavior was to only search for non-overlapping sgRNA sites in the reference sequence in [`32e1e97`](https://github.com/pinellolab/CRISPResso2/commit/32e1e9797da5c3033cdc588e92f06b8813961953)

- Assert correct input fastq file format in [`7248ba8`](https://github.com/pinellolab/CRISPResso2/commit/7248ba8c4deee125ad1ec12fdf1294a84d5f6f93)

- Update plotCustomAllelePlot.py script for [#292](https://github.com/pinellolab/CRISPResso2/issues/292) by [@kclem](https://github.com/kclem) in [#293](https://github.com/pinellolab/CRISPResso2/pull/293)

- Clarify CRISPRessoWGS intended use by [@Colelyman](https://github.com/Colelyman) in [#303](https://github.com/pinellolab/CRISPResso2/pull/303)

- Case-insensitive headers accepted in CRISPRessoPooled [`e577318`](https://github.com/pinellolab/CRISPResso2/commit/e577318006cd17b2725bd028e5e56634c6eb829a)

- Allow dashes in filenames in [`712eb2a`](https://github.com/pinellolab/CRISPResso2/commit/712eb2a11825e8d36f2870deb12b35486bd633fb)

- Sort pandas dataframes by # of reads and sequences so that the order is consistent for testing by [@Snicker7](https://github.com/Snicker7) and [@Colelyman](https://github.com/Colelyman) in [#316](https://github.com/pinellolab/CRISPResso2/pull/316)

- Update `base_editor` parameters in README and add Plot Harness by [@Colelyman](https://github.com/Colelyman) in [#301](https://github.com/pinellolab/CRISPResso2/pull/301)


## v2.2.12 - Protospace Utah - 02/01/2023
### ADDED
- Add deprecation notice in [#260](https://github.com/pinellolab/CRISPResso2/pull/260)

- Add snippet about installing CRISPResso2 via bioconda on Apple silicon in [#274](https://github.com/pinellolab/CRISPResso2/pull/274)

### FIXED
- Fix CRISPRessoPooled bam input in [#265](https://github.com/pinellolab/CRISPResso2/pull/265)

- Fix deprecated numpy type names (fixes [#269](https://github.com/pinellolab/CRISPResso2/issues/269)) in [#270](https://github.com/pinellolab/CRISPResso2/pull/270)

- CRISPRessoPooled custom header fix in [#278](https://github.com/pinellolab/CRISPResso2/pull/278)

### CHANGED
- Status Updates + Pooled Mixed Mode Update in [#279](https://github.com/pinellolab/CRISPResso2/pull/279)


## v2.2.11 - Of Weber - 10/11/2022
### FIXED
- Fix batch quilt plot name by [@Colelyman](https://github.com/Colelyman) in [#249](https://github.com/pinellolab/CRISPResso2/pull/249)

- Batch amplicon plots by [@kclem](https://github.com/kclem) in [#251](https://github.com/pinellolab/CRISPResso2/pull/251)

- Fix typo of CRISPResssoPlot when plotting nucleotide quilt by [@Colelyman](https://github.com/Colelyman) in [#250](https://github.com/pinellolab/CRISPResso2/pull/250)


## v2.2.10 - Overhangs Alameda - 09/15/2022
### ADDED
- Add `--zip_output` parameter to produce a zipped file report by [@Colelyman](https://github.com/Colelyman) and [@Snicker7](https://github.com/Snicker7) in [`c80f828`](https://github.com/pinellolab/CRISPResso2/commit/c80f82838c5a228b79ad4484092877cfee08e02c)

- Allow N's in bam output by [@kclem](https://github.com/kclem) in [`b0b7d41`](https://github.com/pinellolab/CRISPResso2/commit/b0b7d41d697304d0d5fc93e3346c9de1b98ba41d)

- Autodetect reference amplicons from interleaved fastq input

### FIXED
- Fix bug when comparing two samples with the same name in [#228](https://github.com/pinellolab/CRISPResso2/pull/228)

- Fix bug when name is provided instead of amplicon_name in pooled input file in [#229](https://github.com/pinellolab/CRISPResso2/pull/229)

- Fix for aggregate plots in Batch mode in [#237](https://github.com/pinellolab/CRISPResso2/pull/237)

- Fix loading of crispressoInfo from WGS and pooled in [`49740ba`](https://github.com/pinellolab/CRISPResso2/commit/49740ba1d66ed6d13a9e154b8b17bc8b5186581d)

### CHANGED
- Parallel plot refactor in [#247](https://github.com/pinellolab/CRISPResso2/pull/247)

- Add plotly to dockerfile in [`b68a432`](https://github.com/pinellolab/CRISPResso2/commit/b68a43271115251b18e8955e285ccc18f549e8cd)


## v2.2.9 - Long Surrey - 06/23/2022
### ADDED
- fastq_to_bam implementation in [#219](https://github.com/pinellolab/CRISPResso2/pull/219)\
    - If the parameter --bam_output is provided, CRISPResso alignments will be written to a file called 'CRISPResso_output.bam' with the alignments in bam format. If the `bowtie2_index` is provided, alignments will be reported in reference to that genome. If the `bowtie2_index` is not provided, alignments will be reported in reference to a custom reference created by the amplicon sequence(s) and written to the file 'CRISPResso_output.fa'.\
    - This enables the viewing of CRISPResso alignments in other browsers (e.g., IGV). If no `bowtie2_index` is provided, the reference genome should be set to the produced 'CRISPResso_output.fa' file, and then the alignment bam can be loaded into IGV.

### FIXED
- Don't run global frameshift plot when there are no modified reads by [@Colelyman](https://github.com/Colelyman) in [#226](https://github.com/pinellolab/CRISPResso2/pull/226)

### CHANGED
- CRISPRessobatch: put directory in quotes by [@sshen8](https://github.com/sshen8) in [#222](https://github.com/pinellolab/CRISPResso2/pull/222)


## v2.2.8 - Welcome to High Waikato - 05/13/2022
### ADDED
- Interactive plotly summary plots in CRISPRessoAggregate and CRISPRessoBatch for visualizing and comparisons

- CRISPRessoPooled enhancement that allows the amplicons file to have a header and additional columns to be provided

- CRISPRessoCompare generates a report of the number of significant reads at each base

### CHANGED
- Minor bug fixes for plotCustomAllelePlot.py to work with Python3 by [@dharjanto](https://github.com/dharjanto) in [#212](https://github.com/pinellolab/CRISPResso2/pull/212)

- Coerce ints in batch file checking by [@Snicker7](https://github.com/Snicker7) in [#200](https://github.com/pinellolab/CRISPResso2/pull/200)

- Large aggregation by [@Colelyman](https://github.com/Colelyman) in [#192](https://github.com/pinellolab/CRISPResso2/pull/192)

- Flexible pooled input by [@Snicker7](https://github.com/Snicker7) and [@kclem](https://github.com/kclem) in [#217](https://github.com/pinellolab/CRISPResso2/pull/217)


## v2.2.7 - Literature and Los Angeles - 02/11/2022
### ADDED
- Adds features for providing aligned bams as input to CRISPRessoPooled and for a faster demultiplexing when amplicons and genome are provided. The added parameters are:
    - `--aligned_pooled_bam`: Path to aligned input for CRISPRessoPooled processing. If this parameter is specified, the alignments in the given bam will be used to demultiplex reads. If this parameter is not set (default), input reads provided by `--fastq_r1` (and optionally `--fastq_r2`) will be aligned to the reference genome using bowtie2. If the input bam is given, the corresponding reference fasta must also be given to extract reference genomic sequences via the parameter `--bowtie2_index`. Note that the aligned reads are paired-end seqenced, they should already be merged into 1 read (e.g. via Flash) before alignment.
    - `--demultiplex_only_at_amplicons`: If set, and an amplicon file (`--amplicons_file`) and reference sequence (`--bowtie2_index`) are provided, reads overlapping alignment positions of amplicons will be demultiplexed and assigned to that amplicon. If this flag is not set, the entire genome will be demultiplexed and reads with the same start and stop coordinates as an amplicon will be assigned to that amplicon.

### FIXED
- Fix int bug for CRISPRessoPooled n_reads ([`ef15cae`](https://github.com/pinellolab/CRISPResso2/commit/ef15caee29380f58aaae392c897fabe47587486e))

- Fix deprecated pandas indexers ([`eea442a`](https://github.com/pinellolab/CRISPResso2/commit/eea442a763e0c6c41da16a15d6e11ddf6d222dc8), [`f4b6cfc`](https://github.com/pinellolab/CRISPResso2/commit/f4b6cfc03951215c3b9019dc47beb2913a5448ab))

### CHANGED
- Improve performance by removing regex from indel location analysis by [@Colelyman](https://github.com/Colelyman) [#182](https://github.com/pinellolab/CRISPResso2/pull/182)

- Fastq output produced by `--fastq_output` now includes the inserted bases. Previously, a string like "DEL= INS=78(1) SUB= " would indicate a 1bp insertion at site 78. This update outputs strings like "DEL= INS=78(1+G) SUB= " with the insertion described as a plus character followed by the inserted bases. ([`2f84dd0`](https://github.com/pinellolab/CRISPResso2/commit/2f84dd02787abffa6d39efbc50c82c92d1c87528))

- Update ylabel_values -> y_label_values by [@swrosati](https://github.com/swrosati) in [#174](https://github.com/pinellolab/CRISPResso2/pull/174)

- Allow mixed-case prime-editing input ([`e999079`](https://github.com/pinellolab/CRISPResso2/commit/e9990790a0081b765c1f54f4a9b18db522ab4693))


## v2.2.6 - Basepairing Bern - 10/21/2021
### ADDED
- Add param --plot_center to allow custom plots centered at a given point in plotCustomAllelePlot script [ecf23ef](https://github.com/pinellolab/CRISPResso2/commit/ecf23ef23e5701b232bba547a6d7d4b96f085f26)

- Add unit tests [3e6c281](https://github.com/pinellolab/CRISPResso2/commit/3e6c281c931756acd26acca96249b2fa1ad1db31)

### FIXED
- Fix allele plotting error for plot 5 referring to uninitialized y_max variable [53197e6](https://github.com/pinellolab/CRISPResso2/commit/53197e62706e37db54f7ed50c94f38a938955e59)

- Fix unicode errors for bam read/write [8196b6a](https://github.com/pinellolab/CRISPResso2/commit/8196b6a81f477ddcb0e34d61dfb54085de20c1a0)

### CHANGED
- All sub-CRISPResso runs are run with 1 process in Batch, WGS, Pooled, etc. Because we added multiprocessing capabilities to CRISPResso (the plotting part) we thought it would be slick for CRISPRessoPooled to run sqrt(n_processes) CRISPResso processes with sqrt(n_processes) processes each. Unfortunately, sqrt(n) isn't a really useful number for the number of processes people usually run (e.g. 2 or 3 or even 8), and a lot of the CRISPResso processing isn't enabled to take advantage of multiprocessing (e.g. the alignment step), so processes were being wasted. So we reverted back to having n_processes CRISPResso processes, each with 1 process. [a923a7c](https://github.com/pinellolab/CRISPResso2/commit/a923a7c2ef182238bd6b8aa77289bac487b7679b)

- Convert columns in nucleotide count and modification tables to numeric for PE analysis [cabebbe](https://github.com/pinellolab/CRISPResso2/commit/cabebbefb2646967dbeee80af08ac14156b1b53c)

- Make loggers module-specific so matplotlib debug doesn't get spewn out in the CRISPResso log [c2bdd96](https://github.com/pinellolab/CRISPResso2/commit/c2bdd96651eef5af38fb7bbc11d257a827ac080d)

### REMOVED
- Remove version checks for numpy and seaborn [90b43ea](https://github.com/pinellolab/CRISPResso2/commit/90b43eaa03c0ea0fdee62a7b244204cad50056cc)


## v2.2.5 - Immunity from Bonneville - 09/22/2021
### FIXED
- Fixes bug when sequencing reads are much longer than the given reference sequence.


## v2.2.4 - Mutant Maricopa - 09/09/2021
### ADDED
- This release adds an additional parameter --assign_ambiguous_alignements_to_first_allele. For ambiguous alignments, setting this flag will force them to be assigned to the first (as provided by the references -a first and then -e second) amplicon. Thus, no reads will be discarded as 'ambiguous' and all reads will be counted once in the analysis.

### CHANGED
- Batch summaries are produced for amplicons present in only one sample.


## v2.2.3 - Collateral Cardston - 08/30/2021
### FIXED
- Fixes database_id bug


## v2.2.2 - Large Honolulu - 08/20/2021
### FIXED
- For some reason some of the previous commits resolving problems with filterFastqs weren't picked up in v2.2.1. So I'm hoping they'll be included here.


## v2.2.1 - Sequence Length Salt Lake - 08/20/2021
### FIXED
- More unicode bug fixes for filtering fastqs

### CHANGED
- CRISPRessBatch now outputs summary of splicing/frameshift mutation status


## v2.2.0 - Matches Sanpete - 08/13/2021
### CHANGED
- Python 3 release

    - Incorporates updates and changes from python2 up to this point

    - Adds multiprocessing to CRISPRessoBase to parallelize image generation to speed up results
        - CRISPRessoBatch, CRISPRessoPooled, and CRISPRessoWGS allocate processes to sub-CRISPResso commands so that sqrt(n_processes) sub-commands are run, each with sqrt(n_processes) unless plotting is turned off (via --suppress_report or --suppress_plots in which case n_processes are run, each with 1 process.

    - The crispresso_info dictionary containing run information is saved as json so it can be read across versions of python and by other programs (e.g. R). The dict structure of the object has also been changed to be more navigable and hierarchical.

    - Because of changes to crispresso_info, this version of CRISPResso will not be able to finish incomplete runs (e.g. checkpointing of CRISPRessoPooled) that were started by previous python2 versions of CRISPResso, or aggregate or access information from previous runs (e.g. using CRISPRessoAggregate, CRISPRessoCompare, or the custom python plotting scripts in scripts). However, even if we were to have stuck with pickle format for crispresso_info, the python2 and python3 versions were incompatible anyway. So we figured this was a good time to move toward a better format.


## v2.1.3 - Lentiviral Cache - 06/29/2021
## FIXED
- Fixes bug for CRISPRessoPooled analyzes of many amplicons where samtools sort writes status updates that can't be parsed as reads.


## v2.1.2 - Single Guide Washington - 06/23/2021
### ADDED
- Addition of a script for custom allele plotting [947fbab](https://github.com/pinellolab/CRISPResso2/commit/947fbab70e0f4fa78cc43273b4fe2be5225043cc)

### CHANGED
- Updates to CRISPRessoPooledWGSCompare, used for comparing multiple amplicons in CRISPRessoWGS or CRISPRessoPooled experiments [48d6c87](https://github.com/pinellolab/CRISPResso2/commit/48d6c8724f541b285e5d6889723cc37fc99e5cfa)
    - CRISPRessoPooledWGSCompare now produces html report linking to sub-CRISPRessoCompare reports


## v2.1.1 - Nicking Hancock - 05/22/2021
### ADDED
- This release incorporates changes to make bowtie2 alignment in CRISPRessoPooled more permissive [4dc9e7](https://github.com/pinellolab/CRISPResso2/commit/44dc9e75c660f4b3b683f1c80f5a964aa55e75bd), and remove duplicate rows in the Alleles_frequency_table.txt due to reads being in the forward or reverse direction [0e08cd0](https://github.com/pinellolab/CRISPResso2/commit/0e08cd05c2f3279ac95a068be76f1d36a4b0224d).

### CHANGED
- When given a genome file, CRISPRessoPooled aligns reads to the genome using the Bowtie2 aligner. The legacy parameters were somewhat strict. The new parameters reflect the 'default_min_aln_score' parameter in allowing for substantially more indels and mismatches than previous.
    - The parameter --use_legacy_bowtie2_options_string has been added to use the legacy settings. Otherwise, the bowtie2 alignment settings will be calculated as follows:
    - --end-to-end - no clipping, match bonus -ma is set to 0
    - -N 0 number of mismatches allowed in seed alignment
    - --np 0 where read (or ref have ambiguous character (N)) penalty is 0
    - --mp 3,2 mismatch penalty - set max mismatch to -3 to coincide with the gap extension penalty (2 is the default min mismatch penalty)
    - --score-min L,-5,-3*(1-H) For a given homology score, we allow up to (1-H) mismatches (-3) or gap extensions (-3) and one gap open (-5). This score translates to -5 + -3(1-H)L where L is the sequence length


## v2.1.0 - Knockout Lake - 03/23/2021
### CHANGED
- Starting in version 2.1.0, insertion quantification has been changed to only include insertions completely contained by the quantification window. To use the legacy quantification method (i.e. include insertions directly adjacent to the quantification window) please use the parameter --use_legacy_insertion_quantification

- In Prime Editing mode pegRNA spacer sequences given in the incorrect orientation are no longer tolerated

- HDR: Ambiguous alignments don't contribute to the plot 4g (except when --expand_ambiguous_alignments is provided) --fastq_output now also writes alignment scores and alignments for every read


## v2.0.45 - 12/30/2020
### ADDED
- CRISPRessoAggregate can be used to aggregate multiple completed CRISPResso runs.


## v2.0.44 - 11/17/2020
### CHANGED
- Improvements in inferring quantification windows across amplicons/alleles.


## v2.0.43 - 11/07/2020
### ADDED
- Add ticks to appropriate plots

- New parameter --plot_histogram_outlier to plot 100% of data

### FIXED
- Update the function of histograms

- By default 99% of data is shown in plots, now 100% of data is written to data files.


## v2.0.42 - 09/30/2020
### FIXED
- Fixed % character in CRISPRessoPooled arg string

## v2.0.41 - 09/30/2020
### ADDED
- Added --fastq-out parameter to report the CRISPResso analysis separately for each read. Note that this should be used with caution. I'm still trying to figure out what information should be reported for each read, and what format it should be in. Open to feedback on this issue!

### FIXED
- WGS parallelization mode bug fixed

### CHANGED
- WGS and Pooled summary figures scale height based on the number of entries so that they are legible in html reports.


## v2.0.40 - 07/09/2020
### CHANGED
- Prime editing updates - scaffold parameter is now called --prime_editing_pegRNA_scaffold_seq.
Guide names with spaces produce file names with hyphens instead of spaces

## v2.0.39 - 07/07/2020
### ADDED
- Batch mode supports bam and multiple quantification windows

## v2.0.38 - 07/01/2020
### ADDED
- Add new paramter --annotate_wildtype_allele to annotate wildtype alleles on the allele plots

- Input can now be read from bam using the parameter --bam_input and (optionally) --bam_chr_loc to use the reads in the bam at this location as input.

- An output bam is produced with an additional soace-separated field prefixed by c2 (e.g. c2:Z:ALN=Inferred CLASS=Inferred_MODIFIED MODS=D47;I0;S0 DEL=56(47) INS= SUB= ALN_REF=TTGGCGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAGCTTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGGACCCCGCCACCGTGCGCCGGGCCTTGCAGTGGGCGCGCTACCTGCGCCACATCCATCGGCGCTTTGGTCGGCATGGCCCCATTCGCACGGCTCT----------------------------------------------- ALN_SEQ=ACACCGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGA-----------------------------------------------TCGCTTTTCCGAGCTTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGGACCCCGCCACCGTGCGCCGGGCCTTGCAGTGGGCGCGCTACCTGCGCCACATCCATCGGCGCTTTGGTCGGCATGGCCCCATTCGCACGGCTCTGGAGCGGCGGCTGCACAACCAGTGGAGGCAAGAGGGCGGCTTTGGGC). Note that the alignment details (location, cigar string, etc) are not modified.. this may be done in the future). Bam file input cannot be trimmed or pre-processed with quality filtering.

### CHANGED
- Prime editing scaffold incorporation is now more accurate (looks for the scaffold sequence at the expected position directly after the extension sequence). A plot showing the number of bases matching the scaffold, as well as insertions after the extension sequence, and a data file with these numbers is produced. Added parameter --prime_editing_pegRNA_scaffold_min_match_length to define the minimum length required to classify a read as 'Scaffold-incorporated'

- Renamed split_paired_end parameter to --split_interleaved_input for interleaved input

- Auto mode now considers 5000 reads to detect amplicon sequences

- Update output when reporting missing files -- only lists first 15 files in the current directory and directory of input parameter

- --reference https instead of http


## v2.0.37 - 05/09/2020
### ADDED
- Max processors can be used in WGS and Pooled modes by setting -p max

- Prime editing analysis can be performed by specifying the parameters:
--prime_editing_pegRNA_spacer_seq
--prime_editing_pegRNA_extension_seq
and optionally
--prime_editing_pegRNA_extension_quantification_window_size
--prime_editing_pegRNA_scaffold_sequence
--prime_editing_nicking_guide_seq
with a summary shown in the report

- Extended read analysis data available with --write_detailed_allele_table flag

### CHANGED
- CRISPRessoPooled demultiplexing is performed in parallel and with reduced filesystem demand

- N's don't count as substitutions

- Nucleotide plots are shaded when the nucleotide matches the reference sequence

- sgRNA improvements:
    sgRNA annotations are plotted on multiple lines if they overlap
    sgRNAs can have their own cut site and quantification window size


## v2.0.34 - 04/06/2020
### ADDED
- Pooled Set flag to skip reporting problematic regions


## v2.0.33 - 04/03/2020
### CHANGED
- Plotting computation window is shaded

- Parallelization and checkpointing of CRISPRessoWGS and Pooled

- Increase of alignment efficiency of CRISPRessoPooled amplicons in genome +amplicons mode


## v2.0.32 - 02/25/2020
### CHANGED
- Plotting updates, dsODN detection, and general improvements and bug fixes.


## v2.0.31 - 09/26/2019
### ADDED
- Add custom post-processing plot functions for allele tables

### FIXED
- Fix CRISPRessoPooled handling of chromosomes with underscores

### CHANGED
- Update dependency requirements


## v2.0.30 - 07/02/2019
### ADDED
- Add nucleotide summary for batch mode

### FIXED
- Fix bug for reporting amplicons with no reads

### CHANGED
- Case-insensitive checking for guides


## v2.0.29 - 05/30/2019
### CHANGED
- By default, the html report is created on the outside of the output folder, so if the output is:
CRISPResso_on_SAMPLE/
the html report will be at
CRISPResso_on_SAMPLE.html

- This functionality can be reverted to place the report inside of the output folder using the parameter --place_report_in_output_folder
which will place the html report at:
CRISPResso_on_SAMPLE/CRISPResso2_report.html


## v2.0.28 - 05/24/2019
### ADDED
- Standardize file names

- Add CRISPREssoCompare output html

### CHANGED
- CRISPRessoBatch guide-specific output are plotted as separate plots

- Standardize window definitions (plot window and quantification window specify the distance from the cut site to the edge of the window, so the entire window is 2*plot window)


## v2.0.27 - 04/05/2019
### Added
- Add reports for pooled and WGS

- Add Batch pickle info

### CHANGED
- More precise plotting of cleavage cut site and quantification window

- Bioconda updates


## v2.0.26 - 03/06/2019
### Added
- Add report display name, remove paths from stored files, fix sgRNA plot, CRISPRessoPooled report HTML, add citation to report


## v2.0.25 - 02/21/2019
### Added
- Add inferring of guides


## v2.0.24 - 02/13/2019
### Changed
- Update docker, setup.py


## v2.0.23 - 01/24/2019
### Added
- Add manifest.in


## v2.0.22 - 01/23/2019
### Changed
- Change license location, license update


## v2.0.21 - 01/22/2019
### Changed
- Detangled root location dependency from params


## v2.0.20b - 01/22/2019
### Changed
- Prepare for bioconda integration
