import gzip
import argparse
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(description="Filter reads based on sequence presence")
    parser.add_argument("--fastq_r1", type=str, help="Fastq R1 to filter", required=True)
    parser.add_argument("--fastq_r2", type=str, help="Fastq R2 to filter", default=None)
    parser.add_argument("--fastq_r1_out", type=str, help="Fastq R1 output file to contain filtered reads", default=None)
    parser.add_argument("--fastq_r2_out", type=str, help="Fastq R2 output file to contain filtered reads", default=None)
    parser.add_argument("--include_seq", type=str, help="Sequence that must be present in every read. Reads without this sequence will be filtered out. " +
        "This argument may be specified multiple times, and reads must contain at least one of these sequences", default=[], action='append')
    parser.add_argument("--exclude_seq", type=str, help="Sequence that must NOT be present in every read. Reads with this sequence will be filtered out. " +
        "This argument may be specified multiple times, and read must not contain any of these sequences", default=[], action='append')

    args = parser.parse_args()

    if not args.include_seq and not args.exclude_seq:
        raise ValueError("You must specify either --include_seq or --exclude_seq")

    fastq_r1_out = args.fastq_r1_out
    if fastq_r1_out is None:
        fastq_r1_out = args.fastq_r1.replace('.fastq', '').replace('.fq', '').replace('.gz', '') + ".filtered.fq"
        if args.fastq_r1.endswith('.gz'):
            fastq_r1_out += ".gz"

    fastq_r2_out = args.fastq_r2_out
    if args.fastq_r2 is not None and fastq_r2_out is None:
        fastq_r2_out = args.fastq_r2.replace('.fastq', '').replace('.fq', '').replace('.gz', '') + ".filtered.fq"
        if args.fastq_r2.endswith('.gz'):
            fastq_r2_out += ".gz"

    # CREATION OF FILEHANDLES##
    if args.fastq_r1.endswith('.gz'):
        f1_in = gzip.open(args.fastq_r1, 'rt')
    else:
        f1_in = open(args.fastq_r1, 'rt')
    if fastq_r1_out.endswith('.gz'):
        f1_out = gzip.open(fastq_r1_out, 'wt')
    else:
        f1_out = open(fastq_r1_out, 'w')

    if args.fastq_r2:
        if args.fastq_r2.endswith('.gz'):
            f2_in = gzip.open(args.fastq_r2, 'rt')
        else:
            f2_in = open(args.fastq_r2, 'rt')
        if fastq_r2_out.endswith('.gz'):
            f2_out = gzip.open(fastq_r2_out, 'wt')
        else:
            f2_out = open(fastq_r2_out, 'w')
    # END CREATION OF FILEHANDLES##

    print('Fastq R1: %s' % args.fastq_r1)
    print('Fastq R2: %s' % args.fastq_r2)
    print('Fastq R1 out: %s' % fastq_r1_out)
    print('Fastq R2 out: %s' % fastq_r2_out)

    print('Include seq: %s' % args.include_seq)
    print('Exclude seq: %s' % args.exclude_seq)

    read_read_count = 0
    read_written_count = 0
    read_include_count = defaultdict(int)
    read_exclude_count = defaultdict(int)
    if args.fastq_r2:
        while True:
            r1_id_line = f1_in.readline().rstrip()
            r1_seq_line = f1_in.readline().rstrip()
            r1_plus_line = f1_in.readline().rstrip()
            r1_qual_line = f1_in.readline().rstrip()

            r2_id_line = f2_in.readline().rstrip()
            r2_seq_line = f2_in.readline().rstrip()
            r2_plus_line = f2_in.readline().rstrip()
            r2_qual_line = f2_in.readline().rstrip()

            if not r1_id_line:
                break
            read_read_count += 1

            include_read = True
            if args.include_seq != []:
                include_read = False
                for include_seq in args.include_seq:
                    if include_seq in r1_seq_line or include_seq in r2_seq_line:
                        include_read = True
                        read_include_count[include_seq] += 1

            exclude_read = False
            if args.exclude_seq != []:
                for exclude_seq in args.exclude_seq:
                    if exclude_seq in r1_seq_line or exclude_seq in r2_seq_line:
                        exclude_read = True
                        read_exclude_count[exclude_seq] += 1

            if include_read and not exclude_read:
                read_written_count += 1
                f1_out.write("%s\n%s\n%s\n%s\n" % (r1_id_line, r1_seq_line, r1_plus_line, r1_qual_line))
                f2_out.write("%s\n%s\n%s\n%s\n" % (r2_id_line, r2_seq_line, r2_plus_line, r2_qual_line))

        print('Printed %s/%s reads to %s and %s' % (read_written_count, read_read_count, fastq_r1_out, fastq_r2_out))

    else:
        while True:
            r1_id_line = f1_in.readline().rstrip()
            r1_seq_line = f1_in.readline().rstrip()
            r1_plus_line = f1_in.readline().rstrip()
            r1_qual_line = f1_in.readline().rstrip()

            if not r1_id_line:
                break
            read_read_count += 1

            include_read = True
            if args.include_seq != []:
                include_read = False
                for include_seq in args.include_seq:
                    if include_seq in r1_seq_line:
                        include_read = True
                        read_include_count[include_seq] += 1

            exclude_read = False
            if args.exclude_seq != []:
                for exclude_seq in args.exclude_seq:
                    if exclude_seq in r1_seq_line:
                        exclude_read = True
                        read_exclude_count[exclude_seq] += 1

            if include_read and not exclude_read:
                read_written_count += 1
                f1_out.write("%s\n%s\n%s\n%s" % (r1_id_line, r1_seq_line, r1_plus_line, r1_qual_line))
        print('Printed %s/%s reads to %s' % (read_written_count, read_read_count, fastq_r1_out))

    for seq in args.include_seq:
        print('Reads containing included sequence %s: %s' % (seq, read_include_count[seq]))
    for seq in args.exclude_seq:
        print('Reads containing excluded sequence %s: %s' % (seq, read_exclude_count[seq]))


if __name__ == "__main__":
    main()
