#!/usr/bin/env python


""" Determine codon counts for each transcript
"""

from __future__ import division

import argparse
import logging
import pprint
import numpy
from plastid import GFF3_TranscriptAssembler, GTF2_TranscriptAssembler
from plastid import BAMGenomeArray
from plastid import ThreePrimeMapFactory, FivePrimeMapFactory
from plastid import SizeFilterFactory
from Bio import SeqIO

__version__ = "0.1"
__author__ = "Lance Parsons"
__author_email__ = "lparsons@princeton.edu"
__copyright__ = "Copyright 2016, Lance Parsons, Princeton University"
__license__ = \
    "BSD 2-Clause License http://www.opensource.org/licenses/BSD-2-Clause"


def main(args, loglevel):

    logging.basicConfig(format="%(levelname)s: %(message)s", level=loglevel)

    logging.debug("Building sequence dictionary")
    seq_dict = SeqIO.index(args.fasta, "fasta")
    logging.debug("Reading Annotations")
    
    if args.gff:
        transcripts = list(
                GFF3_TranscriptAssembler(open(args.gff),
                                         add_three_for_stop=args.add_three))
    elif args.gtf:
        transcripts = list(
                GTF2_TranscriptAssembler(open(args.gtf),
                                         add_three_for_stop=args.add_three))
        
    logging.debug("Reading Alignments")
    alignments = BAMGenomeArray([args.bam])
    
    if sum([args.threeprime,args.fiveprime]) !=1:
        logging.error("Must specify only one and at least one mapping type (--fiveprime or --threeprime)")
        exit(1)
        
        
    if args.threeprime:
        alignments.set_mapping(ThreePrimeMapFactory(offset=args.offset))
    elif args.fiveprime:
        alignments.set_mapping(FivePrimeMapFactory(offset=args.offset))

    alignments.add_filter("size", SizeFilterFactory(min=args.min_length,
                                                    max=args.max_length))
    outfh = open(args.outfile, 'w')
    outfh.write("%s\n" % "\t".join(
        ("gene_id", "gene_name", "codon_seq", "codon_index", "codon_count_sum",
         "position_1_count", "position_2_count", "position_3_count")))
    for (i, transcript) in enumerate(transcripts):
        if(i == 0 or (i + 1) % 100 == 0):
            logging.info("Evaluated %s genes" % (i + 1))
        logging.debug(transcript.get_name())
        logging.debug(pprint.pformat(transcript.attr))
        if len(transcript) <= 0:
            logging.warn("Transcript %s is length zero (0), skipping!",
                         transcript.get_name())
            continue
        if transcript.attr.get("pseudo", None) == "true":
            logging.warn("Transcript %s is a pseudogene, skipping!",
                         transcript.get_name())
            continue
        transcript_seq = transcript.get_cds().get_sequence(seq_dict)
        transcript_counts = transcript.get_cds().get_counts(alignments)
        
        if len(transcript_seq)%3  != 0:
            logging.warn("Transcript %s length (%i) is not a multiple of "
                         "three, skipping!" % (transcript.get_name(),
                                               len(transcript_counts)))
            continue
        num_codons = len(transcript_seq)/3
        logging.debug("Trancript length %i basepairs, %f codons" %
                      (len(transcript_counts), num_codons))
        for codon_index in range(1, int(numpy.floor(num_codons))):
            codon_start = (codon_index - 1) * 3
            codon_stop = codon_start + 3
            codon_seq = transcript_seq[codon_start:codon_stop]
            codon_counts = transcript_counts[codon_start:codon_stop]
            codon_count_sum = sum(codon_counts)
            outfh.write("%s\t%s\t%s\t%i\t%i\t%i\t%i\t%i\n" %
                        (transcript.get_name(),
                         transcript.attr.get("gene", ""),
                         codon_seq, codon_index, codon_count_sum,
                         codon_counts[0], codon_counts[1],
                         codon_counts[2]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Calculates codon usage for Ribo-seq data",
        epilog="As an alternative to the commandline, params can be placed in "
        "a file, one per line, and specified on the commandline like "
        "'%(prog)s @params.conf'.",
        fromfile_prefix_chars='@')
    parser.add_argument("--bam", help="BAM file of mapped reads",
                        metavar="FILE", required=True)

    parser.add_argument("--gff", help="GFF file of gene annotations",
                        metavar="FILE")
    parser.add_argument("--gtf", help="GTF file of gene annotations",
                        metavar="FILE")
    
    parser.add_argument("--fasta", help="Fasta file of genomic sequence",
                        metavar="FILE", required=True)
    parser.add_argument("-o", "--outfile", help="Output filename",
                        metavar="FILE", required=True)
    parser.add_argument("--add_three", default=False, action="store_true",
                        help="If supplied, coding regions will be extended by "
                        "3 nucleotides at their 3\' ends (except for GTF2 "
                        "files that explicitly include `stop_codon` "
                        "features). Use if your annotation file excludes stop "
                        "codons from CDS.")
                        
    parser.add_argument("--threeprime",default=False,action = "store_true",help="Use three prime for the map")
    parser.add_argument("--fiveprime",default=False,action = "store_true",help="Use five prime for the map")

    parser.add_argument("--offset", help="Integer representing the offset "
                        "into the read, starting from the 3' end, at which "
                        "data should be plotted (default: %(default)s)",
                        type=int, default=0, metavar="OFFSET")
    parser.add_argument("--codon_buffer", help="Number of codons to ignore "
                        "after the start and before the stop (default: "
                        "%(default)s)"
                        "", type=int, default=4, metavar="N")
    parser.add_argument("--min_length", help="Minimum length of reads to "
                        "count (default: %(default)s)", type=int, default=29,
                        metavar="N")
    parser.add_argument("--max_length", help="Minimum length of reads to "
                        "count (default: %(default)s)", type=int, default=35,
                        metavar="N")
    parser.add_argument("-v", "--verbose",
                        help="increase output verbosity",
                        action="store_true")
    parser.add_argument("--version", action="version",
                        version="%(prog)s " + __version__)
    args = parser.parse_args()

    # Setup logging
    if args.verbose:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.INFO

    main(args, loglevel)
