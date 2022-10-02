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
            show_cut_point=True,pw_size=None,pw_size_5prime=None,pw_size_3prime=None):
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

        returns:
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

    def copy(self):
        """
        Creates a copy of a guide object

        returns:
            Guide object with copies of values
        """
        copy = Guide()
        copy.seq = self.seq
        copy.name = self.name
        copy.orig_seq=self.orig_seq
        copy.mismatches = self.mismatches[:]
        copy.qw_center = self.qw_center
        copy.qw_size_5prime = self.qw_size_5prime
        copy.qw_size_3prime = self.qw_size_3prime
        copy.show_cut_point = self.show_cut_point
        copy.pw_size_5prime = self.pw_size_5prime
        copy.pw_size_3prime = self.pw_size_3prime

        return copy

    def __str__(self):
        return "Guide '" +str(self.name) + "' (" + str(self.seq) + ")"

    def print(self):
        print("\tseq: " + str(self.seq))
        print("\tname: " + str(self.name))
        print("\torig_seq: " + str(self.orig_seq))
        print("\tmismatches: " + str(self.mismatches))
        print("\tqw_center: " + str(self.qw_center))
        print("\tqw_size: " + str(self.qw_size))
        print("\tqw_size_5prime: " + str(self.qw_size_5prime))
        print("\tqw_size_3prime: " + str(self.qw_size_3prime))
        print("\tshow_cut_point: " + str(self.show_cut_point))
        print("\tpw_size: " + str(self.pw_size))
        print("\tpw_size_5prime: " + str(self.pw_size_5prime))
        print("\tpw_size_3prime: " + str(self.pw_size_3prime))

    def set_sequence(self,seq):
        """
        Set the seq for this guide

        Args:
            seq (str): The seq for this guide
        """
        wrong_nt=CRISPRessoShared.find_wrong_nt(seq)
        if wrong_nt:
            raise CRISPRessoShared.NTException('The sgRNA sequence ' + seq + ' contains bad characters:%s'  % ' '.join(wrong_nt))
        self.seq = seq
        if self.orig_seq == '':
            self.orig_seq = seq

    def set_name(self,name):
        """
        Set the name for this guide

        Args:
            name (str): The name for this guide
        """
        if name is not None:
            self.name = name
        else:
            self.name = ''

    def set_original_sequence(self,orig_seq):
        """
        Set the orig_seq for this guide

        Args:
            orig_seq (str): The orig_seq for this guide
        """
        wrong_nt=CRISPRessoShared.find_wrong_nt(orig_seq)
        if wrong_nt:
            raise CRISPRessoShared.NTException('The sgRNA sequence ' + orig_seq + ' contains bad characters:%s'  % ' '.join(wrong_nt))
        self.orig_seq = orig_seq

    def set_mismatches(self,mismatches):
        """

        Args:
            mismatches (array of ints): The mismatches for this guide
        """
        self.mismatches = mismatches

    def set_qw_center(self,qw_center):
        """
        Set the quantification window center for this guide
        Args:
            
            qw_center (str): The quantification window center for this guide
        """
        self.qw_center = _get_assert_int(qw_center,'quantification_window_center')

    def set_qw_size(self,qw_size):
        """
        Set both the quantification window sizes (5' and 3')

        Args:
            qw_size (int): The quantification window to set for this guide
        """
        val = _get_assert_int(qw_size,'quantification_window_size')
        self.qw_size_5prime = val
        self.qw_size_3prime = val

    def set_qw_size_5prime(self,qw_size_5prime):
        """
        Set the qw_size_5prime for this guide
        
        Args:
            qw_size_5prime (int): The qw_size_5prime for this guide
        """
        if qw_size_5prime is not None:
            self.qw_size_5prime = _get_assert_int(qw_size_5prime,'quantification_window_size_5prime')

    def set_qw_size_3prime(self,qw_size_3prime):
        """
        Set the 3' quantification window size for this guide
        
        Args:
            qw_size_3prime (int): The 3' quantification window for this guide
        """
        if qw_size_3prime is not None:
            self.qw_size_3prime = _get_assert_int(qw_size_3prime,'quantification_window_size_3prime')

    def set_show_cut_point(self,show_cut_point):
        """
        Set the show_cut_point for this guide - whether to show the cut point on plots
        
        Args:
            show_cut_point (boolean): Whether to show the cut point on plots for this guide
        """
        self.show_cut_point = show_cut_point

    def set_pw_size(self,pw_size):
        """
        set both 5' and 5' plot window size
        """
        val = _get_assert_int(pw_size,'plot_window_size')
        self.pw_size_5prime = val
        self.pw_size_3prime = val

    def set_pw_size_5prime(self,pw_size_5prime):
        """
        Set the 5' plotting window for this guide
        
        Args:
            pw_size_5prime (int): The 5' plotting window for this guide
        """
        if pw_size_5prime is not None:
            self.pw_size_5prime = _get_assert_int(pw_size_5prime,'plot_window_size_5prime')

    def set_pw_size_3prime(self,pw_size_3prime):
        """
        Set the 3' plotting window for this guide
        
        Args:
            pw_size_3prime (int): The 3' plotting window for this guide
        """
        if pw_size_3prime is not None:
            self.pw_size_3prime = _get_assert_int(pw_size_3prime,'plot_window_size_3prime')

    def get_sequence(self):
        """
        Get the sequence for this guide
        
        Returns: 
            sequence (str): the sequence for this guide
        """
        return self.seq

    def get_name(self):
        """
        Get the name for this guide
        
        Returns: 
            name (str): the name for this guide
        """
        return self.name

    def get_original_seq(self):
        """
        Get the original_seq for this guide
        
        Returns: 
            original_seq (str): the original_seq for this guide
        """
        return self.orig_seq

    def get_mismatches(self):
        """
        Get the mismatches for this guide
        
        Returns: 
            mismatches (str): the mismatches for this guide
        """
        return self.mismatches

    def get_qw_center(self):
        """
        Get the qw_center for this guide
        
        Returns: 
            qw_center (str): the qw_center for this guide
        """
        return self.qw_center

    #can't get quantification window - get either 5' or 3'
#    def get_qw_size(self):

    def get_qw_size_5prime(self):
        """
        Get the 5' quantification window size for this guide
        
        Returns: 
            qw_size_5prime (int): the 5' quantification window size for this guide
        """
        return self.qw_size_5prime

    def get_qw_size_3prime(self):
        """
        Get the 3' quantification window for this guide
        
        Returns: 
            qw_size_3prime (int): the 3' quantification window size for this guide
        """
        return self.qw_size_3prime

    def get_show_cut_point(self):
        """
        Get the show_cut_point for this guide
        
        Returns: 
            show_cut_point (bool): whether to show the cut point for this guide
        """
        return self.show_cut_point

#    can't get plotting window size for this guide - get either 3' or 5'
#    def get_pw_size(self):

    def get_pw_size_5prime(self):
        """
        Get the 5' plotting window size for this guide
        
        Returns: 
            pw_size_5prime (int): the 5' plotting window size for this guide
        """
        return self.pw_size_5prime

    def get_pw_size_3prime(self):
        """
        Get the 3' plotting window size for this guide
        
        Returns: 
            pw_size_3prime (int): the 3' plotting window size for this guide
        """
        return self.pw_size_3prime


class Guide_instance:
    """
    Information for a guide in an amplicon including locations where the guide cuts, ect.

    guide: Guide object for which this Guide_instance contains data
    sequence : guide sequence (matches amplicon sequence)
    interval : indexes in the amplicon where the guide aligns
    orientation: whether guide is aligned in the forward ('F') or reverse ('R') direction relative to the amplicon sequence
    cut_point : cut point index in the amplicon
    quantification_idx: indices along the reference sequence for quantification
    plot_idx : indices along reference sequence for which to plot the allele plot (allele frequency plot around sgRNA)
    name_in_ref : names of sgRNAs (in case there are multiple matches for a single sgRNA)
    mismatches_in_ref : indices along guide that are mismatched against a 'flexiguide_seq'

    """
    def __init__(self, guide=None, sequence='', interval=None, orientation='F', cut_point=None,
            quantification_idx=None, plot_idx=None, name_in_ref='',mismatches_in_ref=[]):
        self.guide = guide
        self.sequence = sequence
        self.interval = interval
        self.orientation = orientation
        self.cut_point = cut_point
        self.quantification_idx = quantification_idx
        self.plot_idx = plot_idx
        self.name_in_ref = name_in_ref
        self.mismatches_in_ref = mismatches_in_ref

    def set_guide(self, guide):
        """
        Set the Guide object this Guide_in_amplicon represents

        Args:
            guide (Guide): the Guide this Guide_in_amplicon represents
        """
        self.guide= guide

    def set_sequence(self, sequence):
        """
        Set the sequence for this guide in this amplicon

        Args:
            sequence (str): the sequence for this guide in this amplicon
        """
        self.sequence = sequence

    def set_interval(self, interval):
        """
        Set the alignment interval for this guide in this amplicon
        
        Args:
            interval (tuple of (int,int)): the alignment interval for this guide in this amplicon (start, stop)
        """
        if len(interval) != 2:
            raise CRISPRessoShared.RuntimeException('Guide_instance interval must be a tuple of length 2 (got "'+",".join(interval)+'")')
        self.interval = interval

    def set_orientation(self, orientation):
        """
        Set the orientation for the alignment guide in this amplicon

        Args:
            cut_point (str): the orientation for the alignment of this guide in this amplicon, 'F' for forward, or 'R' for reverse
        """
        if orientation not in ['F','R']:
            raise CRISPRessoShared.RuntimeException('Guide_instance orientation must be "F" for forward or "R" for reverse (got "'+str(orientation)+'")')
        self.orientation = orientation

    def set_cut_point(self, cut_point):
        """
        Set the cut_point for this guide in this amplicon

        Args:
            cut_point (int): the cut_point for this guide in this amplicon
        """
        self.cut_point = _get_assert_int(cut_point, 'Guide cut point')

    def set_quantification_idx(self, quantification_idx):
        """
        Set the quantification_idx for this guide in this amplicon

        Args:
            quantification_idx (array of int): the indices for quantification for this guide in this amplicon
        """
        self.quantification_idx = quantification_idx

    def set_plot_idx(self, plot_idx):
        """
        Set the plot_idx for this guide in this amplicon

        Args:
            plot_idx (array of int): the indices to plot for this guide in this amplicon
        """
        self.plot_idx = plot_idx

    def set_name_in_ref(self, name_in_ref):
        """
        Set the name_in_ref for this guide in this amplicon

        Args:
            name_in_ref (str): the name_in_ref for this guide in this amplicon
        """
        self.name_in_ref = name_in_ref

    def set_mismatches_in_ref(self, mismatches_in_ref):
        """
        Set the mismatches_in_ref for this guide in this amplicon

        Args:
            mismatches_in_ref (array of int): the mismatches_in_ref for this guide in this amplicon
        """
        self.mismatches_in_ref = mismatches_in_ref

    def get_guide(self):
        """
        Get the Guide object this Guide_instance represents

        Returns:
            guide (Guide): the Guide this Guide_instance represents
        """
        return self.guide

    def get_sequence(self):
        """
        Get the sequence for this guide in this amplicon

        Returns:
            sequence (str): the sequence for this guide in this amplicon
        """
        return self.sequence

    def get_interval(self):
        """
        Get the alignment interval for this guide in this amplicon
        
        Returns:
            interval (tuple of (int,int)): the alignment interval for this guide in this amplicon
        """
        return self.interval

    def get_orientation(self):
        """
        Get the alignment orientation for this guide in this amplicon
        
        Returns:
        orientation (str): the alignment orientation for this guide in this amplicon, 'F' for forward, 'R' for reverse
        """
        return self.orientation
        

    def get_cut_point(self):
        """
        Get the cut_point for this guide in this amplicon

        Returns:
            cut_point (int): the cut_point for this guide in this amplicon
        """
        return self.cut_point
    
    def get_quantification_idx(self):
        """
        Get the quantification_idx for this guide in this amplicon

        Returns:
            plot_idx (array of int): the quantification_idx for this guide in this amplicon
        """
        return self.quantification_idx

    def get_plot_idx(self):
        """
        Get the plot_idx for this guide in this amplicon

        Returns:
            plot_idx (array of int): the plot_idx for this guide in this amplicon
        """
        return self.plot_idx

    def get_name_in_ref(self):
        """
        Get the name_in_ref for this guide in this amplicon

        Returns:
            name_in_ref (str): the name_in_ref for this guide in this amplicon
        """
        return self.name_in_ref

    def get_mismatches_in_ref(self):
        """
        Get the mismatches_in_ref for this guide in this amplicon

        Returns:
            mismatches_in_ref (array of int): the mismatches_in_ref for this guide in this amplicon
        """
        return self.mismatches_in_ref


class Amplicon:
    def __init__(self, name=None, sequence=None, sequence_length=None, min_aln_score=None, gap_incentive=None,
            guide_instances=[], sgRNA_cut_points=None, sgRNA_intervals=None, sgRNA_sequences=None, sgRNA_plot_cut_points=None,
            sgRNA_plot_idxs=None, sgRNA_names=None, sgRNA_mismatches=None, sgRNA_orig_sequences=None, contains_guide=None,
            contains_coding_seq=None, exon_positions=None, exon_len_mods=None, exon_intervals=None, splicing_positions=None,
            include_idxs=None, exclude_idxs=None, idx_cloned_from=None, fw_seeds=None, rc_seeds=None,
            aln_genome=None, aln_chr=None, aln_start=None, aln_end=None, aln_strand=None):

        self.name=name
        self.sequence=sequence
        self.sequence_length=sequence_length
        self.min_aln_score=min_aln_score
        self.gap_incentive=gap_incentive
        self.guide_instances=guide_instances
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


    def set_name(self,name):
        """
        Set name for this amplicon

        Args:
            name (str): Name for this amplicon
        """
        self.name = name

    def set_seq(self,seq):
        """
        Set sequence for this amplicon

        Args:
            seq (str): sequence of this amplicon
        """
        self.seq = seq

    def set_sequence(self,sequence):
        """
        Set sequence for this amplicon

        Args:
            sequence (str): sequence for this guide
        """
        self.sequence = sequence

    def set_sequence_length(self,sequence_length):
        """
        Set sequence_length for this guide

        Args:
            sequence_length (int): length of this amplicon sequence
        """
        self.sequence_length = _get_assert_int(sequence_length,'sequence_length')

    def set_min_aln_score(self,min_aln_score):
        """
        Set min_aln_score for this amplicon

        Args:
            min_aln_score (int): min_aln_score for this amplicon
        """
        self.min_aln_score = _get_assert_int(min_aln_score,'min_aln_score')

    def set_gap_incentive(self,gap_incentive):
        """
        Set gap_incentive for this amplicon

        Args:
            gap_incentive (array of int): gap_incentive for this amplicon with incentives where there are cut sites
        """
        self.gap_incentive = _get_assert_int(gap_incentive,'gap_incentive')

    def set_sgRNA_cut_points(self,sgRNA_cut_points):
        """
        Set sgRNA_cut_points for this amplicon

        Args:
            sgRNA_cut_points (array of int): sgRNA_cut_points for this amplicon
        """
        self.sgRNA_cut_points = sgRNA_cut_points

    def set_sgRNA_intervals(self,sgRNA_intervals):
        """
        Set sgRNA_intervals for this amplicon

        Args:
            sgRNA_intervals (array of indices): sgRNA_intervals for this amplicon
        """
        self.sgRNA_intervals = sgRNA_intervals

    def set_contains_coding_seq(self,contains_coding_seq):
        """
        Set contains_coding_seq for this amplicon

        Args:
            contains_coding_seq (boolean): boolean for if this amplicon contains the coding sequence
        """
        self.contains_coding_seq = contains_coding_seq

    def set_exon_positions(self,exon_positions):
        """
        Set exon_positions for this amplicon

        Args:
            exon_positions (array of int): exon_positions for this amplicon
        """
        self.exon_positions = exon_positions

    def set_exon_len_mods(self,exon_len_mods):
        """
        Set exon_len_mods for this amplicon

        Args:
            exon_len_mods (int): modification lengths for exons for this amplicon
        """
        self.exon_len_mods = exon_len_mods

    def set_exon_intervals(self,exon_intervals):
        """
        Set exon_intervals for this amplicon

        Args:
            exon_intervals (array of int): exon_intervals for this amplicon
        """
        self.exon_intervals = exon_intervals

    def set_splicing_positions(self,splicing_positions):
        """
        Set splicing_positions for this amplicon

        Args:
            splicing_positions (range of int): splicing_positions for this amplicon
        """
        self.splicing_positions = splicing_positions

    def set_include_idxs(self,include_idxs):
        """
        Set include_idxs for this amplicon

        Args:
            include_idxs (array of int): indices to be included in quantification for this amplicon
        """
        self.include_idxs = include_idxs

    def set_exclude_idxs(self,exclude_idxs):
        """
        Set exclude_idxs for this amplicon

        Args:
            exclude_idxs (array of int): indices to be excluded in quantification for this amplicon
        """
        self.exclude_idxs = exclude_idxs

    def set_idx_cloned_from(self,idx_cloned_from):
        """
        Set idx_cloned_from for this amplicon

        Args:
            idx_cloned_from (int): index of the amplicon this amplicon was cloned from
        """
        self.idx_cloned_from = idx_cloned_from

    def set_fw_seeds(self,fw_seeds):
        """
        Set fw_seeds for this amplicon

        Args:
            fw_seeds (int): alignment seeds for forward alignment to this amplicon
        """
        self.fw_seeds = fw_seeds

    def set_rc_seeds(self,rc_seeds):
        """
        Set rc_seeds for this amplicon

        Args:
            rc_seeds (int): alignment seeds for reverse-complement alignment to this amplicon
        """
        self.rc_seeds = rc_seeds

    def set_aln_genome(self,aln_genome):
        """
        Set aln_genome for this amplicon

        Args:
            aln_genome (str): genome this amplicon is aligned to
        """
        self.aln_genome = aln_genome

    def set_aln_chr(self,aln_chr):
        """
        Set aln_chr for this amplicon

        Args:
            aln_chr (str): chromosome this amplicon is aligned to
        """
        self.aln_chr = aln_chr

    def set_aln_start(self,aln_start):
        """
        Set aln_start for this amplicon

        Args:
            aln_start (int): index of alignment start for this amplicon
        """
        self.aln_start = aln_start

    def set_aln_end(self,aln_end):
        """
        Set aln_end for this amplicon

        Args:
            aln_end (int): index of alignment end for this amplicon
        """
        self.aln_end = aln_end

    def set_aln_strand(self,aln_strand):
        """
        Set aln_strand for this amplicon

        Args:
            aln_strand (int): which strand this amplicon aligns to
        """
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
            for idx, guide_seq in enumerate(args.guide_seq.split(",")):
                guide = Guide.default_from_args(default_args=args)
                guide.set_sequence(guide_seq)
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
                raise CRISPRessoShared.BadParameterException("More quantification_window_center values were given than guides. Guides: %s quantification_window_center: %s"%(args.guide_seq, args.quantification_window_center))
            guides[idx].set_qw_center(qw_center)

    if 'quantification_window_size' in args:
        for idx, qw_size in enumerate(args.quantification_window_size.split(",")):
            if idx > len(guides):
                raise CRISPRessoShared.BadParameterException("More quantification_window_size values were given than guides. Guides: %s quantification_window_size: %s"%(args.guide_seq, args.quantification_window_size))
            guides[idx].set_qw_size(qw_size)

    if 'quantification_window_size_5prime' in args and args.quantification_window_size_5prime is not None:
        for idx, qw_size_5prime in enumerate(args.quantification_window_size_5prime.split(",")):
            if idx > len(guides):
                raise CRISPRessoShared.BadParameterException("More quantification_window_size_5prime values were given than guides. Guides: %s quantification_window_size_5prime: %s"%(args.guide_seq, args.quantification_window_size_5prime))
            guides[idx].set_qw_size_5prime(qw_size_5prime)

    if 'quantification_window_size_3prime' in args and args.quantification_window_size_3prime is not None:
        for idx, qw_size_3prime in enumerate(args.quantification_window_size_3prime.split(",")):
            if idx > len(guides):
                raise CRISPRessoShared.BadParameterException("More quantification_window_size_3prime values were given than guides. Guides: %s quantification_window_size_3prime: %s"%(args.guide_seq, args.quantification_window_size_3prime))
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
    Useful for when users provided multiple values for multiple guides as a parameter
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
