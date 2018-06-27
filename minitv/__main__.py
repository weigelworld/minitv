#!/usr/bin/env python3
"""
minitv

Alignment viewer using AliTV.
"""

__author__ = "Travis Wrightsman"

import argparse
from datetime import datetime
import dendropy
from io import StringIO
import itertools
import logging
import minitv.align
import minitv.alitv
import minitv.annotation
import minitv.phylo
from minitv.seq import FASTAFile
import os
from pathlib import Path
from pprint import pformat
from subprocess import CalledProcessError
import sys
import tempfile


def run_minitv(minitv_args):
    """ Main entry point """
    logging.info('Starting minitv.')

    logging.debug("Arguments:\n{}".format(pformat(vars(minitv_args), indent=1)))

    fasta_file_paths = tuple([Path(minitv_args.referenceFASTA)] + [Path(query) for query in minitv_args.queryFASTA])
    logging.info('Parsing sequences')
    fasta_files = [FASTAFile(path) for path in fasta_file_paths]

    # start alignment subprocesses in background
    aligner_pool, aligner_results = minitv.align.async_align_all_pairwise(
        fasta_files, minitv_args.work_dir, minitv_args.aligner, minitv_args.aligner_args, minitv_args.concurrent_aln)

    # Begin loading AliTV data
    ali_tv = minitv.alitv.AliTV(fasta_files)

    # process annotation
    if minitv_args.annotation:
        annotations = []
        annotation_types = minitv_args.annotation_types.split(',')
        for i, annotation_file_path in enumerate(minitv_args.annotation):
            logging.info("Parsing annotation for {}".format(fasta_files[i].name))
            annotations.append(minitv.annotation.GFF3Annotation(annotation_file_path, annotation_types))
        ali_tv.load_annotations(annotations, minitv_args.annotation_types.split(','))

    # wait for aligner subprocesses
    aligner_pool.join()

    # parse the alignments
    fasta_file_pairs = itertools.combinations(fasta_files, r=2)
    alignment_groups = {}
    aligner_err_file_path = minitv_args.work_dir / "{aligner}/{aligner}.err".format(aligner=minitv_args.aligner)
    with open(aligner_err_file_path, 'ab') as aligner_err_file:
        all_alignments_successful = True
        for fasta_file_pair in fasta_file_pairs:
            target = fasta_file_pair[0]
            query = fasta_file_pair[1]
            aligner_result = aligner_results[(target.name, query.name)]
            try:
                aligner_out, aligner_err = aligner_result.get()
                alignment_groups[target.name, query.name] = getattr(minitv.align,
                                                                    'parse_output_{}'.format(minitv_args.aligner))(
                                                                    StringIO(aligner_out.decode('utf-8')))
            except CalledProcessError as e:
                aligner_err = e.stderr
                all_alignments_successful = False
            finally:
                aligner_err_file.write("--- {} alignment between {} and {} at {} ---\n".format(
                    minitv_args.aligner, target.name, query.name, datetime.now()).encode('utf-8'))
                aligner_err_file.write(aligner_err)
        if not all_alignments_successful:
            raise RuntimeError("One or more alignment processes failed. See the stderr output in {}.".format(
                aligner_err_file_path))

    # process the tree
    if minitv_args.tree:
        if not minitv_args.no_tree:
            tree = dendropy.Tree.get(path=minitv_args.tree, schema='newick')
        else:
            logging.warning("Tree file provided on command line but --no_tree was also specified, ignoring tree file.")
            tree = None
    elif not minitv_args.no_tree:
        tree = minitv.phylo.tree_from_alignment_groups(alignment_groups)
        logging.debug("Calculated UPGMA tree: {}".format(tree.as_string('newick', suppress_rooting=True)).rstrip())
    else:
        tree = None

    logging.info('Parsing tree')
    ali_tv.load_tree(tree)

    ali_tv.set_soft_filters(minitv_args.min_link_identity, minitv_args.min_link_length)
    # Hard filter the alignments
    if minitv_args.hard_filter:
        logging.info('Hard filtering alignments')
        filtered_alignment_groups = {}
        for target, query in alignment_groups:
            filtered_alignment_groups[target, query] = []
            for alignment in alignment_groups[target, query]:
                if (alignment.identity >= minitv_args.min_link_identity) and (
                        abs(alignment.query_end - alignment.query_start) >= minitv_args.min_link_length) and (
                        abs(alignment.target_end - alignment.target_start) >= minitv_args.min_link_length):
                    filtered_alignment_groups[target, query].append(alignment)
        alignment_groups = filtered_alignment_groups

    logging.info('Loading alignments')
    ali_tv.load_links(alignment_groups)
    ali_tv.order_and_orient_sequences(alignment_groups)

    if minitv_args.min_aln_cov:
        logging.info('Hiding chromosomes by alignment coverage')
        ali_tv.hide_chromosomes_by_alignment_coverage(minitv_args.min_aln_cov)

    if minitv_args.min_ref_cov:
        logging.info('Hiding chromosomes by reference alignment coverage')
        ali_tv.hide_chromosomes_by_reference_coverage(minitv_args.min_ref_cov)

    ali_tv.optimize_configuration()

    logging.info('Generating JSON output')

    return ali_tv.get_json(indent=1)


def main():

    parser = argparse.ArgumentParser(
        prog="python3 minitv.py"
    )

    parser.add_argument(
        "referenceFASTA",
        metavar="reference.fa",
        help="FASTA file of reference sequence"
    )

    parser.add_argument(
        "queryFASTA",
        metavar="query.fa",
        nargs='+',
        help="FASTA file of query sequence(s). This argument can instead take a single file of file names containing "
             "paths relative to the working directory."
    )

    parser.add_argument(
        "-r",
        "--region",
        help="Bounds of region to view, in standard samtools region format (refSeqID:start-end)"
    )

    parser.add_argument(
        '-t', '--tree',
        help="Path to file containing tree in Newick format"
    )

    parser.add_argument(
        '-a', '--annotation',
        help="Path to a GFF3-format annotation file. One argument per GFF file, in the same order of the genome FASTA "
             "files. This argument can also be supplied once with a file of filenames, including "
             "the path to the reference annotation at the top.",
        action="append"
    )

    parser.add_argument(
        '--annotation-types',
        help="Comma-separated list of feature types to display from annotation. Default: gene",
        default="gene"
    )

    parser.add_argument(
        '-d', '--debug',
        help="Output detailed debugging messages",
        action="store_const",
        dest="logLevel",
        const=logging.DEBUG,
        default=logging.WARNING
    )

    parser.add_argument(
        "-v", "--verbose",
        help="Output progress and other informative messages",
        action="store_const",
        dest="logLevel",
        const=logging.INFO
    )

    parser.add_argument(
        "--work_dir",
        help="Directory to save intermediate files to, such as alignment output and logs.",
        type=Path,
        default=None
    )

    parser.add_argument(
        "--aligner",
        help="Name of a minitv-supported aligner to use",
        type=str,
        choices=['minimap2'],
        default='minimap2'
    )

    parser.add_argument(
        "--aligner_args",
        help="A string of additional command line arguments to be passed to the aligner",
        type=str,
        default=''
    )

    parser.add_argument(
        "--concurrent_aln",
        help="Maximum number of concurrent aligner subprocesses. (This is NOT the number of threads per alignment "
        "subprocess)",
        type=int,
        default=4
    )

    parser.add_argument(
        "--min_link_identity",
        help="Minimum link identity for AliTV to initially display",
        type=float,
        default=0.0
    )

    parser.add_argument(
        "--min_link_length",
        help="Minimum link length for AliTV to initially display",
        type=int,
        default=0
    )

    parser.add_argument(
        "--min_aln_cov",
        help="Minimum percent coverage, from 0 to 100, by alignments that a sequence needs to be initially visible in "
             "AliTV.",
        type=float,
        default=0.0
    )

    parser.add_argument(
        "--min_ref_cov",
        help="Minimum percent coverage, from 0 to 100, of a query sequence by alignments to the reference needed for"
             "the query sequence to be initially visible in AliTV.",
        type=float,
        default=0.0
    )

    parser.add_argument(
        "--hard_filter",
        help="In addition to setting AliTV link filters, remove filtered links from JSON output",
        action="store_true"
    )

    parser.add_argument(
        "--no_tree",
        help="Disable automatic tree prediction and order genomes by their order on the command line",
        action="store_true"
    )

    args = parser.parse_args()

    # set up logging to stderr
    root = logging.getLogger()
    root.setLevel(args.logLevel)
    stderr_log_handler = logging.StreamHandler(sys.stderr)
    stderr_log_handler.setLevel(args.logLevel)
    stderr_log_formatter = logging.Formatter("{asctime} [{module}:{levelname}] {message}", style='{')
    stderr_log_handler.setFormatter(stderr_log_formatter)
    root.addHandler(stderr_log_handler)

    # Argument checking/processing
    #  queryFASTA
    if len(args.queryFASTA) == 1:
        # check for special case where queryFASTA is a file of file names
        possible_fofn_path = Path(args.queryFASTA[0])
        query_paths = []
        with open(possible_fofn_path) as possible_fofn:
            first_line = possible_fofn.readline()
            if not first_line.startswith('>'):
                # most likely not a FASTA file, but rather a file of file names
                query_paths.append(first_line.rstrip())
                for remaining_line in possible_fofn:
                    query_paths.append(remaining_line.rstrip())
                args.queryFASTA = query_paths

    #  annotation
    if args.annotation and (len(args.annotation) == 1):
        # check for special case where annotation in a file of file names
        possible_fofn_path = Path(args.annotation[0])
        annotation_paths = []
        with open(possible_fofn_path) as possible_fofn:
            first_line = possible_fofn.readline()
            if not (first_line.startswith('#') or ("\t" in first_line)):
                # most likely not a GFF file, but rather a file of file names
                annotation_paths.append(first_line.rstrip())
                for remaining_line in possible_fofn:
                    annotation_paths.append(remaining_line.rstrip())
                args.annotation = annotation_paths

    if args.annotation and (len(args.annotation) != len(args.queryFASTA) + 1):
        # Is the number of annotation files provided equal to the number of genome files provided?
        sys.exit('The number of genomes provided does not match the number of annotation files provided.')

    #  work_dir
    if not args.work_dir:
        # create a temporary working directory if none provided
        temp_work_dir = tempfile.TemporaryDirectory()
        args.work_dir = Path(temp_work_dir.name)
        logging.debug("No working directory specified. Created temporary working directory at '{}'".format(
            temp_work_dir.name
        ))
    else:
        args.work_dir = Path(os.path.expanduser(args.work_dir))
        if not args.work_dir.exists():
            sys.exit("Specified working directory doesn't exist: '{}'".format(args.work_dir))

    #  min_link_identity
    if 0 > args.min_link_identity > 100:
        sys.exit("min_link_identity must be between 0 and 100")

    #  min_link_length
    if 0 > args.min_link_length:
        sys.exit("min_link_length must be at least 0")

    #  min_aln_cov
    if 0 > args.min_aln_cov > 100:
        sys.exit("min_aln_cov must be between 0 and 100")

    #  min_ref_cov
    if 0 > args.min_ref_cov > 100:
        sys.exit("min_ref_cov must be between 0 and 100")

    # send the AliTV JSON to stdout
    print(run_minitv(args))


if __name__ == '__main__':
    main()
