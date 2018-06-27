"""
Handles alignment software and their outputs
"""

from enum import Enum
from itertools import combinations
import logging
from multiprocessing.pool import ThreadPool
from pathlib import Path
import subprocess


class Strand(Enum):
    """ File schema independent method to encode alignment strand. """
    FORWARD = 1
    UNKNOWN = 0
    REVERSE = -1


def async_align_all_pairwise(fasta_files, work_dir, aligner, extra_args, max_concurrent):
    """Calls the appropriate asynchronous alignment pipeline for the given aligner.

    Returns a tuple of two objects:

    1. a multiprocessing.pool.ThreadPool object containing all of the threads watching the alignment subprocesses
    2. a dictionary with keys of tuples of (target path stem, query path stem) and values of
    multiprocessing.pool.AsyncResult, which can be queried for the return values of the alignment subprocess

    **DO NOT** read the stdout and stderr of these AsyncResult objects until join() is called on the parent ThreadPool
    object to ensure that all subprocesses are finished.

    :param fasta_files: fasta files to be aligned
    :type fasta_files: collections.Iterable[minisyn.seq.FASTAFile]
    :param work_dir: working directory to output alignment files
    :type work_dir: pathlib.Path
    :param aligner: name of the supported aligner to use
    :type aligner: str
    :param extra_args: extra arguments to be passed to the aligner subprocess
    :type extra_args: str
    :param max_concurrent: maximum number of concurrent aligner subprocesses
    :type max_concurrent: int
    :return: alignment subprocess-watcher threads and their results
    :rtype: ThreadPool, dict
    """

    return globals()["async_align_all_pairwise_{}".format(aligner)](fasta_files, work_dir, extra_args, max_concurrent)


def async_align_all_pairwise_minimap2(fasta_files, work_dir, extra_args, max_concurrent):
    """Starts the asynchronous indexing and alignment minimap2 subprocesses and returns objects connected to their
    results.

    :param fasta_files: fasta files to be aligned
    :type fasta_files: collections.Iterable[minisyn.seq.FASTAFile]
    :param work_dir: working directory to output alignment files
    :type work_dir: pathlib.Path
    :param extra_args: extra arguments to be passed to the minimap2 subprocess
    :type extra_args: str
    :param max_concurrent: maximum number of concurrent minimap2 subprocesses
    :type max_concurrent: int
    :return: minimap2 subprocess-watcher threads and their results
    :rtype: ThreadPool, dict
    """

    fasta_file_pairs = combinations(fasta_files, r=2)
    indexing_results = {}
    mapping_results = {}

    minimap2_work_dir = work_dir / Path('minimap2')

    minimap2_work_dir.mkdir(exist_ok=True)

    index_tp = ThreadPool(max_concurrent)
    for fasta_file in fasta_files:
        index_path = minimap2_work_dir / Path("{}.mmi".format(fasta_file.path.stem))
        indexing_results[fasta_file.name] = index_tp.apply_async(
            minimap2_index,
            args=(fasta_file.path, index_path, extra_args))
    index_tp.close()

    map_tp = ThreadPool(max_concurrent)
    for fasta_file_pair in fasta_file_pairs:
        target = fasta_file_pair[0]
        query = fasta_file_pair[1]
        target_index_path = minimap2_work_dir / Path("{}.mmi".format(target.path.stem))
        mapping_results[tuple([target.name, query.name])] = map_tp.apply_async(
            minimap2_align,
            args=(target_index_path, query.path, extra_args, index_tp))
    map_tp.close()

    return map_tp, mapping_results


def minimap2_index(sequence_file_path, index_path, extra_args):
    """Indexes a given sequence file using minimap2.

    :param sequence_file_path: path to the sequence file to be indexed
    :type sequence_file_path: pathlib.Path
    :param index_path: path to the output index file
    :type sequence_file_path: pathlib.Path
    :param extra_args: extra arguments passed to minimap2
    :type extra_args: str
    :return: stdout and stderr of minimap2 index process
    :rtype: bytes, bytes
    """

    minimap2_index_command_line = "minimap2 {} -d {} {}".format(extra_args, index_path, sequence_file_path)

    logging.debug("Starting minimap2 indexing ({})".format(minimap2_index_command_line))

    minimap2_index_process = subprocess.Popen(
        minimap2_index_command_line,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    out, err = minimap2_index_process.communicate()

    if minimap2_index_process.returncode != 0:
        raise subprocess.CalledProcessError(
            returncode=minimap2_index_process.returncode,
            cmd=minimap2_index_command_line,
            output=out,
            stderr=err)

    logging.debug("minimap2 indexing of {} has finished.".format(sequence_file_path.stem))

    return out, err


def minimap2_align(target_path, query_path, extra_args, indexing_pool=None):
    """Runs minimap2 on given target and query.

    :param target_path: path to target genome FASTA file
    :type target_path: pathlib.Path
    :param query_path: path to query genome FASTA file
    :type query_path: pathlib.Path
    :param extra_args: extra arguments passed to minimap2
    :type extra_args: str
    :param indexing_pool: indexing ThreadPool object to wait on before starting alignment process
    :type indexing_pool: multiprocessing.pool.ThreadPool
    :return: stdout and stderr of minimap2 alignment process
    :rtype: bytes, bytes
    """

    minimap2_align_command_line = "minimap2 {} {} {}".format(extra_args, target_path, query_path)

    if indexing_pool:
        # wait for a ThreadPool of minimap2 indexing subprocesses
        indexing_pool.join()

    logging.debug("Starting minimap2 alignment ({}).".format(minimap2_align_command_line))

    minimap2_align_process = subprocess.Popen(
        minimap2_align_command_line,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    out, err = minimap2_align_process.communicate()

    if minimap2_align_process.returncode != 0:
        raise subprocess.CalledProcessError(
            returncode=minimap2_align_process.returncode,
            cmd=minimap2_align_command_line,
            output=out,
            stderr=err)

    logging.debug("minimap2 alignment between {} and {} has finished.".format(target_path.stem, query_path.stem))

    return out, err


class Alignment:
    """
    A single alignment region between two genomes.

    This class implements values required to be converted to AliTV Link objects.
    """

    def __init__(self, query_name, query_start, query_end,
                 target_name, target_start, target_end,
                 identity, strand, length):
        """

        :param query_name: Query sequence name
        :type query_name: str
        :param query_start: Query start (0-based)
        :type query_start: int
        :param query_end: Query end (0-based)
        :type query_end: int
        :param target_name: Target sequence name
        :type target_name: str
        :param target_start: Target start (0-based)
        :type target_start: int
        :param target_end: Target end (0-based)
        :type target_end: int
        :param identity: Percent identity of a sequence, represented as a float between 0 and 100.
        :type identity: float
        :param strand: Relative strand
        :type strand: Strand
        :param length: Alignment length
        :type length: int
        """

        self.query_name = query_name
        self.query_start = query_start
        self.query_end = query_end

        self.target_name = target_name
        self.target_start = target_start
        self.target_end = target_end

        self.identity = identity
        if isinstance(strand, Strand):
            self.strand = strand
        else:
            raise ValueError(
                "{} is not a valid Alignment.strand value. Alignment objects need to be instantiated with a member of"
                " minitv.align.Strand.".format(type(strand)))
        self.length = length

    def __len__(self):
        return self.length


class PAFAlignment(Alignment):
    """
    A single alignment between two genomes in PAF format.
    """

    def __init__(self, query_name, query_len, query_start, query_end, strand, target_name, target_len, target_start,
                 target_end, matches, length, mapping_quality, tags):
        """

        :param query_name: Query sequence name
        :type query_name: str
        :param query_len: Query sequence length
        :type query_len: int
        :param query_start: Query start (0-based)
        :type query_start: int
        :param query_end: Query end (0-based)
        :type query_end: int
        :param strand: Relative strand: "+" or "-"
        :type strand: str
        :param target_name: Target sequence name
        :type target_name: str
        :param target_len: Target sequence length
        :type target_len: int
        :param target_start: Target start on original strand (0-based)
        :type target_start: int
        :param target_end: Target end on original strand (0-based)
        :type target_end: int
        :param matches: Number of residue matches
        :type matches: int
        :param length: Alignment block length
        :type length: int
        :param mapping_quality: Mapping quality (0-255; 255 for missing)
        :type mapping_quality: int
        :param tags: A dictionary of SAM-like typed key-value pairs
        :type tags: dict
        """

        if strand == '+':
            strand = Strand.FORWARD
        elif strand == '-':
            strand = Strand.REVERSE
        else:
            raise ValueError("{} is an invalid PAF alignment strand value.".format(strand))

        identity = (matches / length) * 100
        super().__init__(query_name, query_start, query_end, target_name, target_start, target_end, identity,
                         strand, length)

        self.query_len = query_len
        self.target_len = target_len
        self.matches = matches

        if mapping_quality < 0 or mapping_quality > 255:
            raise ValueError("'{}' is an invalid mapping quality for the PAF format. Mapping quality must be between 0 "
                             "and 255.".format(mapping_quality))
        else:
            self.mapping_quality = mapping_quality

        self._tags = tags

    def get_tag(self, key):
        return self._tags[key]

    def set_tag(self, key, value):
        self._tags[key] = value


def parse_output_minimap2(output):
    """

    :param output: minimap2 stdout
    :type output: io.TextIOBase or io.StringIO
    :rtype: list[Alignment]
    """

    alignments = []

    for line in output:
        split_line = line.rstrip().split("\t")

        query_name = split_line[0]
        query_len = int(split_line[1])
        query_start = int(split_line[2])
        query_end = int(split_line[3])
        strand = split_line[4]
        target_name = split_line[5]
        target_len = int(split_line[6])
        target_start = int(split_line[7])
        target_end = int(split_line[8])
        matches = int(split_line[9])
        length = int(split_line[10])
        mapping_quality = int(split_line[11])
        tags = {}

        if len(split_line) > 12:
            # parse SAM tags
            key_value_pairs = split_line[12:]
            for key_value_pair in key_value_pairs:
                key, type_char, value = key_value_pair.split(':')
                if type_char == "i":
                    tags[key] = int(value)
                elif type_char == 'f':
                    tags[key] = float(value)
                else:
                    tags[key] = value

        new_alignment = PAFAlignment(query_name, query_len, query_start, query_end, strand, target_name,
                                     target_len, target_start, target_end, matches, length, mapping_quality,
                                     tags)
        alignments.append(new_alignment)

    return alignments
