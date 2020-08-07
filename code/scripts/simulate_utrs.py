#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Add simulated 5' UTR regions to GFF records.
    Wrapper script for workaround provided by Joshua Griffin Dunn
    <https://www.linkedin.com/in/joshuagriffindunn>
    https://github.com/joshuagryphon/plastid/issues/3#issuecomment-202572593
"""

import argparse
import logging

from plastid import GFF3_TranscriptAssembler
from plastid import GenomicSegment
from plastid import Transcript

__version__ = "0.1"
__author__ = "Lance Parsons"
__author_email__ = "lparsons@princeton.edu"
__copyright__ = "Copyright 2019, Lance Parsons, Princeton University"
__license__ = \
    "BSD 2-Clause License http://www.opensource.org/licenses/BSD-2-Clause"


def main(args, loglevel):

    logging.basicConfig(format="%(levelname)s: %(message)s", level=loglevel)

    ORFS = GFF3_TranscriptAssembler(open(args.gff_file))
    with open(args.output_file, 'w') as output_file:
        for my_orf in ORFS:
            span = my_orf.spanning_segment
            if my_orf.strand == "+":
                new_region = GenomicSegment(
                        span.chrom,
                        span.start - args.upstream_utr_length,
                        span.end + args.downstream_utr_length,
                        span.strand)
            else:
                new_region = GenomicSegment(
                        span.chrom,
                        span.start - args.downstream_utr_length,
                        span.end + args.upstream_utr_length,
                        span.strand)
            # copy metadata attributes from old ORF
            logging.debug(my_orf.attr)
            new_transcript = Transcript(new_region, **my_orf.attr)
            logging.debug(new_transcript)
            output_file.write(new_transcript.as_gff3())


def parse_arguments(input_args=None):
    parser = argparse.ArgumentParser(
        description="Adds simluated UTR region to ORFS in GFF3 file",
        epilog="As an alternative to the commandline, params can be placed in "
        "a file, one per line, and specified on the commandline like "
        "'%(prog)s @params.conf'.",
        fromfile_prefix_chars='@')
    parser.add_argument("gff_file",
                        help="Source GFF3 file",
                        metavar="INPUT_FILE")
    parser.add_argument("output_file",
                        help="Output GFF3 file",
                        metavar="OUTPUT_FILE")
    parser.add_argument("--upstream_utr_length", type=int, default=50,
                        help="Length of simulated upstream UTR "
                        "(default: %(default)s)")
    parser.add_argument("--downstream_utr_length", type=int, default=50,
                        help="Length of simulated downstream UTR "
                        "(default: %(default)s)")
    parser.add_argument("-v", "--verbose",
                        help="increase output verbosity",
                        action="store_true")
    parser.add_argument("--version", action="version",
                        version="%(prog)s " + __version__)
    args = parser.parse_args(input_args)
    return(args)


try:
    snakemake
except NameError:
    snakemake = None

if snakemake is not None:
    args = parse_arguments(
            [snakemake.input[0], snakemake.output[0]])
    main(args, logging.INFO)
elif __name__ == '__main__':
    # Standard boilerplate to call the main() function to begin
    # the program.
    args = parse_arguments()
    # Setup logging
    if args.verbose:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.INFO

    main(args, loglevel)
