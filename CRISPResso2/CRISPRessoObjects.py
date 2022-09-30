from CRISPResso2 import CRISPRessoShared

class Guide:
    """
    mismatches : a list of positions where a guide may have mismatches (e.g. for flexiguides)
    names : guide name
    quantification_window_center : position where quantification is centered for this guide
    quantification_window_sizes : size of quantification windows extending from quantification_window_center
    quantification_window_sizes_5prime: size of quantification windows extending in the sgRNA 5' direction from quantification_window_center. Overrides quantification_window_size if set, otherwise if None, is set to quantification_window_sizes.
    quantification_window_sizes_3prime: size of lengths of quantification windows extending in the sgRNA 3' direction from quantification_window_center. Overrides quantification_window_size if set, otherwise if None, is set to quantification_window_sizes.
    pw_size : size of window extending from quantification_window_center to plot
    pw_size_5prime: size of window extending from quantification_window_center to plot upstream (sgRNA 5' side). Overrides plot_window_size if set, otherwise if None, is set to plot_window_size.
    pw_size_3prime: size of window extending from quantification_window_center to plot downstream (sgRNA 3' side). Overrides plot_window_size if set, otherwise if None, is set to plot_window_size.
    show_cut_point : whether or not to add cut point to plot (prime editing flaps don't have cut points)
    """
    def __init__(self,seq='',name='',orig_seq='',mismatches=[],
            qw_center=None,qw_size=None,qw_size_5prime=None,qw_size_3prime=None,
            show_cut_point=False,pw_size=None,pw_size_5prime=None,pw_size_3prime=None):
        self.seq = seq
        self.name = name
        self.orig_seq=orig_seq
        self.mismatches = mismatches
        self.qw_center = qw_center
        self.qw_size = qw_size
        self.qw_size_5prime = qw_size_5prime
        self.qw_size_3prime = qw_size_3prime
        self.show_cut_point = show_cut_point
        self.pw_size = pw_size
        self.pw_size_5prime = pw_size_5prime
        self.pw_size_3prime = pw_size_3prime

        if seq != '' and orig_seq == '':
            self.orig_seq = seq

    @classmethod
    def default_from_args(cls,default_args):
        """
        Creates a new guide object using CRISPResso args parsed from the command line
        By default, if a value is a comma-separated string, it takes the first value. Eg. for --guide_names one,two,three it would set the default to 'one'

        params:
        default_args: arg object from argparser

        returs:
        Guide object with default values
        """
        self = cls()

        self.qw_center = _get_first_int_value_from_comma_sep_array(default_args.quantification_window_center, 'quantification_window_center')
        self.qw_size = _get_first_int_value_from_comma_sep_array(default_args.quantification_window_size, 'quantification_window_size')

        self.qw_size_5prime = _get_first_int_value_from_comma_sep_array(default_args.quantification_window_size_5prime, 'quantification_window_size_5prime')
        if self.qw_size_5prime is None:
            self.qw_size_5prime = self.qw_size

        self.qw_size_3prime = _get_first_int_value_from_comma_sep_array(default_args.quantification_window_size_3prime, 'quantification_window_size_3prime')
        if self.qw_size_3prime is None:
            self.qw_size_3prime = self.qw_size

        self.pw_size = _get_first_int_value_from_comma_sep_array(default_args.plot_window_size, 'plot_window_size')

        self.pw_size_5prime = _get_first_int_value_from_comma_sep_array(default_args.plot_window_size_5prime, 'plot_window_size')
        if self.pw_size_5prime is None:
            self.pw_size_5prime = self.pw_size

        self.pw_size_3prime = _get_first_int_value_from_comma_sep_array(default_args.plot_window_size_3prime, 'plot_window_size')
        if self.pw_size_3prime is None:
            self.pw_size_3prime = self.pw_size

        return self

    def __copy__(self):
        """
        Creates a copy of a guide object

        returs:
            Guide object with copies of values
        """
        copy = Guide()
        copy.seq = self.seq
        copy.name = self.name
        copy.orig_seq=self.rig_seq
        copy.mismatches = self.mismatches[:]
        copy.qw_center = self.qw_center
        copy.qw_size = self.qw_size
        copy.qw_size_5prime = self.qw_size_5prime
        copy.qw_size_3prime = self.qw_size_3prime
        copy.show_cut_point = self.show_cut_point
        copy.pw_size = self.pw_size
        copy.pw_size_5prime = self.pw_size_5prime
        copy.pw_size_3prime = self.pw_size_3prime

        return copy

    def set_sequence(self,seq):
        wrong_nt=CRISPRessoShared.find_wrong_nt(seq)
        if wrong_nt:
            raise CRISPRessoShared.NTException('The sgRNA sequence ' + seq + ' contains bad characters:%s'  % ' '.join(wrong_nt))
        self.seq = seq

    def set_name(self,name):
        if name is not None:
            self.name = name
        else:
            self.name = ''

    def set_original_sequence(self,orig_seq):
        wrong_nt=CRISPRessoShared.find_wrong_nt(orig_seq)
        if wrong_nt:
            raise CRISPRessoShared.NTException('The sgRNA sequence ' + orig_seq + ' contains bad characters:%s'  % ' '.join(wrong_nt))
        self.orig_seq = orig_seq

    def set_qw_center(self,qw_center):
        self.qw_center = _get_assert_int(qw_center,'quantification_window_center')

    def set_qw_sizes(self,qw_size):
        """
        set all the quantification window sizes (general, 5' and 3')
        """
        val = _get_assert_int(qw_size,'quantification_window_size')
        self.qw_site = val
        self.qw_size_5prime = val
        self.qw_size_3prime = val

    def set_qw_size(self,qw_size):
        self.qw_size = _get_assert_int(qw_size,'quantification_window_size')

    def set_qw_size_5prime(self,qw_size_5prime):
        self.qw_size_5prime = _get_assert_int(qw_size_5prime,'quantification_window_size_5prime')

    def set_qw_size_3prime(self,qw_size_3prime):
        self.qw_size_3prime = _get_assert_int(qw_size_3prime,'quantification_window_size_3prime')

    def set_show_cut_point(self,show_cut_point):
        self.show_cut_point = show_cut_point

    def set_pw_sizes(self,pw_size):
        """
        set all the plot window sizes (general, 5' and 3')
        """
        val = _get_assert_int(pw_size,'plot_window_size')
        self.pw_size = val
        self.pw_size_5prime = val
        self.pw_size_3prime = val

    def set_pw_size(self,pw_size):
        self.pw_size = _get_assert_int(pw_size,'plot_window_size')

    def set_pw_size_5prime(self,pw_size_5prime):
        self.pw_size_5prime = _get_assert_int(pw_size_5prime,'plot_window_size_5prime')

    def set_pw_size_3prime(self,pw_size_3prime):
        self.pw_size_3prime = _get_assert_int(pw_size_3prime,'plot_window_size_3prime')

    def get_seq(self):
        return self.seq
    def get_name(self):
        return self.name
    def get_rig_seq(self):
        return self.rig_seq
    def get_mismatches(self):
        return self.mismatches
    def get_qw_center(self):
        return self.qw_center
    def get_qw_size(self):
        return self.qw_size
    def get_qw_size_5prime(self):
        return self.qw_size_5prime
    def get_qw_size_3prime(self):
        return self.qw_size_3prime
    def get_show_cut_point(self):
        return self.show_cut_point
    def get_pw_size(self):
        return self.pw_size
    def get_pw_size_5prime(self):
        return self.pw_size_5prime
    def get_pw_size_3prime(self):
        return self.pw_size_3prime

class Amplicon:
    def __init__(self, seq=None, name=None, sequence=None, sequence_length=None, min_aln_score=None, gap_incentive=None,
            guides=[],
            sgRNA_cut_points=None, sgRNA_intervals=None, sgRNA_sequences=None, sgRNA_plot_cut_points=None,
            sgRNA_plot_idxs=None, sgRNA_names=None, sgRNA_mismatches=None, sgRNA_orig_sequences=None, contains_guide=None,
            contains_coding_seq=None, exon_positions=None, exon_len_mods=None, exon_intervals=None, splicing_positions=None,
            include_idxs=None, exclude_idxs=None, idx_cloned_from=None, fw_seeds=None, rc_seeds=None,
            aln_genome=None, aln_chr=None, aln_start=None, aln_end=None, aln_strand=None, default_args=None):

        self.seq=seq
        self.name=name
        self.sequence=sequence
        self.sequence_length=sequence_length
        self.min_aln_score=min_aln_score
        self.gap_incentive=gap_incentive
        self.guides=guides
        self.sgRNA_cut_points=sgRNA_cut_points
        self.sgRNA_intervals=sgRNA_intervals
        self.sgRNA_sequences=sgRNA_sequences
        self.sgRNA_plot_cut_points=sgRNA_plot_cut_points
        self.sgRNA_plot_idxs=sgRNA_plot_idxs
        self.sgRNA_names=sgRNA_names
        self.sgRNA_mismatches=sgRNA_mismatches
        self.sgRNA_orig_sequences=sgRNA_orig_sequences
        self.contains_guide=contains_guide
        self.contains_coding_seq=contains_coding_seq
        self.exon_positions=exon_positions
        self.exon_len_mods=exon_len_mods
        self.exon_intervals=exon_intervals
        self.splicing_positions=splicing_positions
        self.include_idxs=include_idxs
        self.exclude_idxs=exclude_idxs
        self.idx_cloned_from=idx_cloned_from
        self.fw_seeds=fw_seeds
        self.rc_seeds=rc_seeds
        self.aln_genome=aln_genome
        self.aln_chr=aln_chr
        self.aln_start=aln_start
        self.aln_end=aln_end
        self.aln_strand=aln_strand

    def set_seq(self,seq):
        self.seq = seq
    def set_name(self,name):
        self.name = name
    def set_sequence(self,sequence):
        self.sequence = sequence
    def set_sequence_length(self,sequence_length):
        self.sequence_length = _get_assert_int(sequence_length,'sequence_length')
    def set_min_aln_score(self,min_aln_score):
        self.min_aln_score = _get_assert_int(min_aln_score,'min_aln_score')
    def set_gap_incentive(self,gap_incentive):
        self.gap_incentive = _get_assert_int(gap_incentive,'gap_incentive')
    def set_sgRNA_cut_points(self,sgRNA_cut_points):
        self.sgRNA_cut_points = sgRNA_cut_points
    def set_sgRNA_intervals(self,sgRNA_intervals):
        self.sgRNA_intervals = sgRNA_intervals
    def set_sgRNA_sequences(self,sgRNA_sequences):
        self.sgRNA_sequences = sgRNA_sequences
    def set_sgRNA_plot_cut_points(self,sgRNA_plot_cut_points):
        self.sgRNA_plot_cut_points = sgRNA_plot_cut_points
    def set_sgRNA_plot_idxs(self,sgRNA_plot_idxs):
        self.sgRNA_plot_idxs = sgRNA_plot_idxs
    def set_sgRNA_names(self,sgRNA_names):
        self.sgRNA_names = sgRNA_names
    def set_sgRNA_mismatches(self,sgRNA_mismatches):
        self.sgRNA_mismatches = sgRNA_mismatches
    def set_sgRNA_orig_sequences(self,sgRNA_orig_sequences):
        self.sgRNA_orig_sequences = sgRNA_orig_sequences
    def set_contains_guide(self,contains_guide):
        self.contains_guide = contains_guide
    def set_contains_coding_seq(self,contains_coding_seq):
        self.contains_coding_seq = contains_coding_seq
    def set_exon_positions(self,exon_positions):
        self.exon_positions = exon_positions
    def set_exon_len_mods(self,exon_len_mods):
        self.exon_len_mods = exon_len_mods
    def set_exon_intervals(self,exon_intervals):
        self.exon_intervals = exon_intervals
    def set_splicing_positions(self,splicing_positions):
        self.splicing_positions = splicing_positions
    def set_include_idxs(self,include_idxs):
        self.include_idxs = include_idxs
    def set_exclude_idxs(self,exclude_idxs):
        self.exclude_idxs = exclude_idxs
    def set_idx_cloned_from(self,idx_cloned_from):
        self.idx_cloned_from = idx_cloned_from
    def set_fw_seeds(self,fw_seeds):
        self.fw_seeds = fw_seeds
    def set_rc_seeds(self,rc_seeds):
        self.rc_seeds = rc_seeds
    def set_aln_genome(self,aln_genome):
        self.aln_genome = aln_genome
    def set_aln_chr(self,aln_chr):
        self.aln_chr = aln_chr
    def set_aln_start(self,aln_start):
        self.aln_start = aln_start
    def set_aln_end(self,aln_end):
        self.aln_end = aln_end
    def set_aln_strand(self,aln_strand):
        self.aln_strand = aln_strand

def set_guides_from_args(args):
    """
    Set guide information from command line args

    params:
    args: argparser object with args passed in from command line

    returns:
    guides: array of Guide objects
    """
    guides = []
    if 'guide_seq' in args:
        for current_guide_seq in args.guide_seq.split(','):
            guide = Guide.default_from_args(default_args=args)
            guide.set_sequence(current_guide_seq)
            guides.append(guide)

    #each guide has a name, quantification window centers (relative to the end of the guide), quantification window sizes, and plotting window sizes
    if 'guide_name' in args:
        if args.guide_name:
            for idx, guide_name in enumerate(args.guide_name.split(",")):
                if idx > len(guides):
                    raise CRISPRessoShared.BadParameterException("More guide names were given than guides. Guides: %s Guide names: %s"%(args.guide_seq, args.guide_name))
                guides[idx].set_name(guide_name)

    if 'quantification_window_center' in args:
        for idx, qw_center in enumerate(args.quantification_window_center.split(",")):
            if idx > len(guides):
                raise CRISPRessoShared.BadParameterException("More quantification_window_center values were given than guides. Guides: %s quantication_window_center: %s"%(args.guide_seq, args.quantification_window_center))
            print('setting center to ' + str(qw_center))
            guides[idx].set_qw_center(qw_center)

    if 'quantification_window_size' in args:
        for idx, qw_size in enumerate(args.quantification_window_size.split(",")):
            if idx > len(guides):
                raise CRISPRessoShared.BadParameterException("More quantification_window_size values were given than guides. Guides: %s quantication_window_size: %s"%(args.guide_seq, args.quantification_window_size))
            guides[idx].set_qw_size(qw_size)

    if 'quantification_window_size_5prime' in args and args.quantification_window_size_5prime is not None:
        for idx, qw_size_5prime in enumerate(args.quantification_window_size_5prime.split(",")):
            if idx > len(guides):
                raise CRISPRessoShared.BadParameterException("More quantification_window_size_5prime values were given than guides. Guides: %s quantication_window_size_5prime: %s"%(args.guide_seq, args.quantification_window_size_5prime))
            guides[idx].set_qw_size_5prime(qw_size_5prime)

    if 'quantification_window_size_3prime' in args and args.quantification_window_size_3prime is not None:
        for idx, qw_size_3prime in enumerate(args.quantification_window_size_3prime.split(",")):
            if idx > len(guides):
                raise CRISPRessoShared.BadParameterException("More quantification_window_size_3prime values were given than guides. Guides: %s quantication_window_size_3prime: %s"%(args.guide_seq, args.quantification_window_size_3prime))
            guides[idx].set_qw_size_5prime(qw_size_5prime)

    if 'plot_window_size' in args:
        for idx, pw_size in enumerate(args.plot_window_size.split(",")):
            if idx > len(guides):
                raise CRISPRessoShared.BadParameterException("More plot_window_size values were given than guides. Guides: %s plot_window_size: %s"%(args.guide_seq, args.plot_window_size))
            guides[idx].set_pw_size(pw_size)

    if 'plot_window_size_5prime' in args and args.plot_window_size_5prime is not None:
        for idx, pw_size_5prime in enumerate(args.plot_window_size_5prime.split(",")):
            if idx > len(guides):
                raise CRISPRessoShared.BadParameterException("More plot_window_size_5prime values were given than guides. Guides: %s plot_window_size_5prime: %s"%(args.guide_seq, args.plot_window_size_5prime))
            guides[idx].set_pw_size_5prime(pw_size_5prime)

    if 'plot_window_size_3prime' in args and args.plot_window_size_3prime is not None:
        for idx, pw_size_3prime in enumerate(args.plot_window_size_3prime.split(",")):
            if idx > len(guides):
                raise CRISPRessoShared.BadParameterException("More plot_window_size_3prime values were given than guides. Guides: %s plot_window_size_3prime: %s"%(args.guide_seq, args.plot_window_size_3prime))
            guides[idx].set_pw_size_3prime(pw_size_3prime)


    for guide in guides:
        guide.show_cut_point = True
    if 'base_editor_output' in args:
        if args.base_editor_output:
            for guide in guides:
                guide.show_cut_point = False #whether to plot cut point -- base editor and prime editing flap cut points aren't plotted

    for guide in guides:
        if guide.qw_size_5prime is None:
            guide.qw_size_5prime = guide.qw_size
        if guide.qw_size_3prime is None:
            guide.qw_size_3prime = guide.qw_size
        if guide.pw_size_5prime is None:
            guide.pw_size_5prime = guide.pw_size
        if guide.pw_size_3prime is None:
            guide.pw_size_3prime = guide.pw_size
    return guides

def _get_first_int_value_from_comma_sep_array(vals,property_name_for_display):
    """
    Given a comma-separated string, get the first value as an int
    Useful for when users provided mutliple values for multiple guides as a parameter
    Raises a BadParameterException if value can't be parsed
    Returns None if vals is None
    e.g given '3,5,6' this function would return 3

    params:
        vals: string, comma-separated string to parse
        property_name_for_display: property name to display to user in case of error

    returns:
        ret_val: either None or an int val

    """
    if str(vals) == 'None':
        return None
    val_array = str(vals).split(",")
    try:
        ret_val = int(val_array[0])
    except (NameError, TypeError):
        raise CRISPRessoShared.BadParameterException("%s values must be provided as integers. In '%s' got unexpected value: '%s'" % (property_name_for_display,vals,vals_array[0]))
    return ret_val

def _get_assert_int(val,property_name_for_display):
    """
    Assert a character is an int, then return it as an int
    Raises a BadParameterException if value can't be parsed
    Returns None if val is None

    params:
        vals: string
        property_name_for_display: property name to display to user in case of error
    returns:
        ret_val: int
    """
    if str(val) == 'None':
        return None
    try:
        ret_val = int(val)
    except (NameError, TypeError):
        raise CRISPRessoShared.BadParameterException("%s value must be provided as integers. In '%s' got unexpected value: '%s'" % (property_name_for_display,vals,vals_array[0]))
    return ret_val
