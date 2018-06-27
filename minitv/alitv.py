"""
Handles AliTV integration
"""

from collections import defaultdict
import dendropy
import itertools
import json
import logging
from minitv.align import Strand
from minitv.phylo import tree_to_json_structure


def rec_dd():
    """default_factory for infinitely-recursive defaultdicts.

    Credit
    ======

    https://stackoverflow.com/questions/19189274/defaultdict-of-defaultdict-nested

    :return: infinitely-recursing defaultdict
    :rtype: collections.defaultdict
    """
    return defaultdict(rec_dd)


def merge_intervals(intervals):
    """Generator that yields overlapping intervals from a list.

    Credit
    ======

    https://codereview.stackexchange.com/questions/69242/merging-overlapping-intervals

    :param intervals: a list of tuples (x, y) representing intervals of [x, y]
    :type intervals: list[tuple(numbers.Real, numbers.Real)]
    :return: merged intervals
    :rtype: tuple(numbers.Real, numbers.Real)
    """

    sorted_by_lower_bound = sorted(intervals, key=lambda interval: interval[0])

    if not sorted_by_lower_bound:
        return

    # low and high represent the bounds of the current run of merges
    low, high = sorted_by_lower_bound.pop(0)

    for next_low, next_high in sorted_by_lower_bound:
        if next_low <= high:  # new interval overlaps current run
            high = max(high, next_high)  # merge with the current run
        else:  # current run is over
            yield low, high  # yield accumulated interval
            low, high = next_low, next_high  # start new run

    yield low, high


class Link:
    """
    Represents a connection between two link Features.
    """

    def __init__(self, identifier, target, source, identity):
        """

        :param identifier: unique identifier (within AliTV object) of the Link
        :type identifier: str
        :param target: the connected Feature on the target Chromosome
        :type target: Feature
        :param source: the connected Feature on the source Chromosome
        :type source: Feature
        :param identity: the percent identity of the sequence between the Features linked, multiplied by 100
        :type identity: float
        """

        self.id = identifier
        self.target = target
        self.source = source

        if (identity > 100.0) or (identity < 0.0):
            raise ValueError(
                "'{}' is not a valid identity value for a Link object (0 <= identity <= 100)".format(identity))
        self.identity = identity


class Feature:
    """
    Represents a region on a Chromosome.
    """

    def __init__(self, identifier, name, feature_type, start, end, chromosome):
        """

        :param identifier: unique identifier (within AliTV object) of the Feature
        :type identifier: str
        :param name: the name of the Feature
        :type name: str
        :param feature_type: the type/classification of the feature
        :type feature_type: str
        :param start: the starting base of the feature
        :type start: int
        :param end: the ending base of the feature
        :type end: int
        :param chromosome: the Chromosome containing this Feature
        :type chromosome: Chromosome
        """
        self.id = identifier
        self.name = name
        self.type = feature_type
        self.start = start
        self.end = end
        self.chromosome = chromosome

    def __len__(self):
        return (self.end - self.start) + 1


class Chromosome:
    """
    Represents a contiguous sequence drawn as one box in AliTV.
    """

    def __init__(self, genome, record):
        """

        :param genome: the Genome that this Chromosome belongs to
        :type genome: Genome
        :param record: the Bio.SeqRecord object containing the sequence information of the chromosome
        :type record: Bio.SeqRecord.SeqRecord
        """
        self.genome = genome
        self.seq = record.seq
        self.name = record.id
        self.id = "{}_{}".format(genome.id, self.name)

    def __len__(self):
        return len(self.seq)


class Genome:
    """
    Represents a group of Chromosome objects in a specific order.
    """

    def __init__(self, fasta_file, is_reference=False):
        """

        :param fasta_file: the FASTAFile object containing the Chromosome sequences
        :type fasta_file: minisyn.seq.FASTAFile
        """
        self.fasta = fasta_file
        self.id = fasta_file.name
        self.chromosomes = []
        self._chromosomes_dict = {}
        self.is_reference = is_reference

        for record in fasta_file.records:
            new_chromosome = Chromosome(self, record)
            self.chromosomes.append(new_chromosome)
            self._chromosomes_dict[new_chromosome.name] = new_chromosome

    def __iter__(self):
        return iter(self.chromosomes)

    def get_chromosome_names(self):
        return list(self._chromosomes_dict.keys())

    def get_chromosome_by_name(self, chromosome_name):
        return self._chromosomes_dict[chromosome_name]


class AliTV:
    """
    Represents the AliTV application.
    """

    def __init__(self, fasta_files):
        """Loads, processes, and converts various input data into the AliTV JSON format.

        :param fasta_files: list of FASTA files
        :type fasta_files: list[minitv.seq.FASTAFile]
        """

        self.links = []
        self.features = []

        self.fasta_files = fasta_files

        self._jsonDict = {
            'data': rec_dd(),
            'conf': rec_dd(),
            'filters': rec_dd()
        }

        # populate karyo/chromosomes
        self.genomes = []
        self._jsonDict['filters']['karyo']['genome_order'] = []
        self._jsonDict['filters']['karyo']['order'] = []

        self.genomes.append(Genome(fasta_files[0], is_reference=True))
        for fasta_file in fasta_files[1:]:
            genome = Genome(fasta_file)
            self.genomes.append(genome)

        for genome in self.genomes:
            for chromosome in genome:
                self._jsonDict['data']['karyo']['chromosomes'][chromosome.id] = {
                    'genome_id': genome.id,
                    'name': chromosome.id,
                    'length': len(chromosome),
                    'seq': None
                }
                self._jsonDict['filters']['karyo']['chromosomes'][chromosome.id] = {
                    'visible': True,
                    'reverse': False
                }

                self._jsonDict['filters']['karyo']['order'].append(chromosome.id)

        self._jsonDict['data']['tree'] = None
        self._jsonDict['conf']['tree']['drawTree'] = False
        self._jsonDict['filters']['karyo']['genome_order'] = [genome.id for genome in self.genomes]

        # default configuration and filters
        self._jsonDict['conf']['linear']['drawAllLinks'] = False
        self._jsonDict['conf']['linear']['startLineColor'] = '#1d91c0'
        self._jsonDict['conf']['linear']['endLineColor'] = '#1d91c0'

        self._jsonDict['conf']['graphicalParameters']['canvasWidth'] = 900
        self._jsonDict['conf']['graphicalParameters']['canvasHeight'] = 900
        self._jsonDict['conf']['graphicalParameters']['karyoHeight'] = 30
        self._jsonDict['conf']['graphicalParameters']['karyoDistance'] = 2000
        self._jsonDict['conf']['graphicalParameters']['linkKaryoDistance'] = 20
        self._jsonDict['conf']['graphicalParameters']['tickLabelFrequency'] = 10
        self._jsonDict['conf']['graphicalParameters']['tickDistance'] = 10000
        self._jsonDict['conf']['graphicalParameters']['treeWidth'] = 200
        self._jsonDict['conf']['graphicalParameters']['genomeLabelWidth'] = 200

        self._jsonDict['conf']['minLinkIdentity'] = 0
        self._jsonDict['conf']['midLinkIdentity'] = 50
        self._jsonDict['conf']['maxLinkIdentity'] = 100
        self._jsonDict['conf']['minLinkIdentityColor'] = '#D21414'
        self._jsonDict['conf']['midLinkIdentityColor'] = '#FFEE05'
        self._jsonDict['conf']['maxLinkIdentityColor'] = '#1DAD0A'
        self._jsonDict['conf']['minLinkLength'] = 100
        self._jsonDict['conf']['maxLinkLength'] = 5000
        self._jsonDict['conf']['layout'] = 'linear'
        self._jsonDict['conf']['features']['showAllFeatures'] = False
        self._jsonDict['conf']['features']['fallbackStyle'] = {
            'form': 'rect',
            'color': '#787878',
            'height': 30,
            'visible': False
        }

        self._jsonDict['conf']['labels']['ticks']['showTicks'] = True
        self._jsonDict['conf']['labels']['ticks']['showTickLabels'] = True
        self._jsonDict['conf']['labels']['chromosome']['showChromosomeLabels'] = False
        self._jsonDict['conf']['labels']['genome']['showGenomeLabels'] = True
        self._jsonDict['conf']['labels']['features']['showFeatureLabels'] = False

        self._jsonDict['filters']['links']['minLinkLength'] = 0
        self._jsonDict['filters']['links']['minLinkIdentity'] = 0
        self._jsonDict['filters']['links']['maxLinkIdentity'] = 100
        self._jsonDict['filters']['onlyShowAdjacentLinks'] = True
        self._jsonDict['filters']['showAllChromosomes'] = False

    def optimize_configuration(self):
        """Attempts to optimize AliTV drawing configuration to best display the data"""

        return None

    def set_soft_filters(self, minimum_link_identity, minimum_link_length):
        """Sets soft link filters for AliTV.

        :param minimum_link_identity: minimum link identity for initial AliTV display [0, 100]
        :type minimum_link_identity: float
        :param minimum_link_length: minimum link length for initial AliTV display [0, +inf]
        :type minimum_link_length: float
        """

        if 0 <= minimum_link_identity <= 100:
            self._jsonDict['filters']['links']['minLinkIdentity'] = minimum_link_identity
        else:
            raise ValueError("AliTV minimum link identity must be between 0 and 100")
        if 0 <= minimum_link_length:
            self._jsonDict['filters']['links']['minLinkLength'] = minimum_link_length
        else:
            raise ValueError("AliTV minimum link length must be at least 0")

    def hide_chromosomes_by_alignment_coverage(self, minimum_alignment_coverage):
        """Sets chromosomes invisible based on minimum percent coverage by links.

        :param minimum_alignment_coverage: minimum percent coverage by alignments/links to be visible [0, 100]
        :type minimum_alignment_coverage: float
        """

        alignment_intervals = {}
        for link in self.links:
            if link.target.chromosome not in alignment_intervals:
                alignment_intervals[link.target.chromosome] = []
            alignment_intervals[link.target.chromosome].append((link.target.start, link.target.end))
            if link.source.chromosome not in alignment_intervals:
                alignment_intervals[link.source.chromosome] = []
            alignment_intervals[link.source.chromosome].append((link.source.start, link.source.end))

        self.hide_chromosomes_by_region_coverage(alignment_intervals, minimum_alignment_coverage)

    def hide_chromosomes_by_reference_coverage(self, minimum_reference_coverage):
        """Sets chromosomes invisible based on minimum percent coverage by links to the reference.

        :param minimum_reference_coverage: minimum percent coverage by links to the reference to be visible [0, 100]
        :type minimum_reference_coverage: float
        """

        reference_link_intervals = {}
        for link in self.links:
            if link.target.chromosome.genome.is_reference:
                feature = link.source
                chromosome = feature.chromosome
            elif link.source.chromosome.genome.is_reference:
                feature = link.target
                chromosome = feature.chromosome
            else:
                continue
            if chromosome not in reference_link_intervals:
                reference_link_intervals[chromosome] = []
            reference_link_intervals[chromosome].append((feature.start, feature.end))

        self.hide_chromosomes_by_region_coverage(reference_link_intervals, minimum_reference_coverage)

    def hide_chromosomes_by_region_coverage(self, regions, minimum_region_coverage):
        """Sets chromosomes invisible based on minimum percent coverage by arbitrary regions.

        :param regions: dictionary with Chromosome objects as keys and values of lists of intervals that
            represent regions (overlaps are fine) on the respective Chromosome objects
        :type regions: dict[Chromosome, list[tuple(numbers.Real, numbers.Real)]
        :param minimum_region_coverage: minimum percent coverage by regions to be visible [0, 100]
        """

        for chromosome in regions:
            merged_intervals = merge_intervals(regions[chromosome])
            total_reference_link_coverage = 0
            for low, high in merged_intervals:
                total_reference_link_coverage += (high - low + 1)
            percent_region_coverage = (total_reference_link_coverage / len(chromosome.seq)) * 100
            if percent_region_coverage < minimum_region_coverage:
                logging.debug("Hiding chromosome '{}' because its reference coverage is too low ({} < {})".format(
                    chromosome.id, round(percent_region_coverage, 2), minimum_region_coverage
                ))
                self._jsonDict['filters']['karyo']['chromosomes'][chromosome.id]['visible'] = False

    def load_links(self, alignment_groups):
        """

        :param alignment_groups: a dictionary with keys of tuples (target, query) where target and query are names of
          the FASTA files of the target and query of the alignment, respectively
        :type alignment_groups: dict
        """

        # populate links
        genome_pairs = tuple(itertools.combinations(self.genomes, r=2))
        for genome_pair in genome_pairs:
            target_genome = genome_pair[0]
            query_genome = genome_pair[1]
            alignment_group = alignment_groups[target_genome.fasta.name, query_genome.fasta.name]
            for alignment in alignment_group:
                target_chromosome = target_genome.get_chromosome_by_name(alignment.target_name)
                query_chromosome = query_genome.get_chromosome_by_name(alignment.query_name)

                query_block_id = "b{}".format(str(len(self.features)))
                query_block = Feature(query_block_id, '', 'link', alignment.query_start, alignment.query_end,
                                      query_chromosome)
                self.features.append(query_block)

                target_block_id = "b{}".format(str(len(self.features)))
                target_block = Feature(target_block_id, '', 'link', alignment.target_start, alignment.target_end,
                                       target_chromosome)
                self.features.append(target_block)

                # orient query link block according to strand
                if alignment.strand == Strand.REVERSE:
                    self._jsonDict['data']['features']['link'][query_block.id] = {
                        'karyo': query_block.chromosome.id,
                        'start': query_block.end,
                        'end': query_block.start
                    }
                else:
                    self._jsonDict['data']['features']['link'][query_block.id] = {
                        'karyo': query_block.chromosome.id,
                        'start': query_block.start,
                        'end': query_block.end
                    }

                self._jsonDict['data']['features']['link'][target_block.id] = {
                    'karyo': target_block.chromosome.id,
                    'start': target_block.start,
                    'end': target_block.end
                }

                # AliTV only needs one-way links
                link_id = "link{}".format(len(self.links))
                link = Link(link_id, target_block, query_block, round(alignment.identity, 2))
                self.links.append(link)

                self._jsonDict['data']['links'][query_genome.id][target_genome.id][link_id] = {
                    'target': target_block.id,
                    'source': query_block.id,
                    'identity': link.identity
                }

    def order_and_orient_sequences(self, alignment_groups):
        """Orders and orients each Chromosome in non-reference Genomes based on the first longest alignment between each
        Chromosome and a reference Chromosome.

        The "reference" genome is currently the first genome in the list loaded into AliTV.

        :param alignment_groups: a dictionary with keys of tuples (target, query) where target and query are names of
          the FASTA files of the target and query of the alignment, respectively
        :type alignment_groups: dict
        """

        self._jsonDict['filters']['karyo']['order'] = []

        reference_genome = self.genomes[0]
        for chromosome in reference_genome:
            self._jsonDict['filters']['karyo']['order'].append(chromosome.id)

        for genome in self.genomes[1:]:
            try:
                alignment_group = alignment_groups[reference_genome.fasta.name, genome.fasta.name]
                reference_attr_prefix = 'target'
                non_reference_attr_prefix = 'query'
            except KeyError:
                alignment_group = alignment_groups[genome.fasta.name, reference_genome.fasta.name]
                reference_attr_prefix = 'query'
                non_reference_attr_prefix = 'target'

            longest_alignments = defaultdict(list)
            for alignment in alignment_group:

                # find the longest alignment for each non-reference chromosome
                non_reference_chromosome_name = alignment.__getattribute__(non_reference_attr_prefix + "_name")
                if len(longest_alignments[non_reference_chromosome_name]) < len(alignment):
                    longest_alignments[non_reference_chromosome_name] = alignment

            aligned_chromosomes = set(longest_alignments.keys())
            all_chromosome_names = set([chromosome.name for chromosome in genome])
            reversed_chromosomes = set()
            reference_query_chromosome_order = {}
            for chromosome in reference_genome:
                reference_query_chromosome_order[chromosome.name] = []

            for aligned_chromosome in aligned_chromosomes:
                longest_alignment = longest_alignments[aligned_chromosome]

                if longest_alignment.strand == '-':
                    reversed_chromosomes.add(aligned_chromosome)

                reference_chromosome_name = longest_alignment.__getattribute__(reference_attr_prefix + "_name")
                reference_chromosome_start = longest_alignment.__getattribute__(reference_attr_prefix + "_start")
                non_reference_chromosome_name = longest_alignment.__getattribute__(non_reference_attr_prefix + "_name")
                reference_query_chromosome_order[reference_chromosome_name].append(
                    (genome.get_chromosome_by_name(non_reference_chromosome_name).id,
                     reference_chromosome_start)
                )

            for reversed_chromosome in reversed_chromosomes:
                alitv_karyo_id = genome.get_chromosome_by_name(reversed_chromosome).id
                self._jsonDict['filters']['karyo']['chromosomes'][alitv_karyo_id]['reverse'] = True

            for reference_chromosome_name in reference_query_chromosome_order:
                sorted_query_order = sorted(reference_query_chromosome_order[reference_chromosome_name],
                                            key=lambda x: x[1])
                sorted_query_order_names = [pair[0] for pair in sorted_query_order]
                self._jsonDict['filters']['karyo']['order'] += sorted_query_order_names
            unaligned_chromosome_names = list(all_chromosome_names - aligned_chromosomes)
            unaligned_chromosome_ids = [genome.get_chromosome_by_name(chromosome_name).id for chromosome_name in
                                        unaligned_chromosome_names]
            self._jsonDict['filters']['karyo']['order'] += list(unaligned_chromosome_ids)

    def load_annotations(self, annotations, annotation_types):
        """

        :param annotations: iterable of minitv.annotation.GFF3Annotation objects from genome GFF3 files
        :type annotations: collections.Iterable[minitv.annotation.GFF3Annotation]
        :param annotation_types: iterable of strings of feature types to filter the annotation
        :type annotation_types: collections.Iterable[str]
        """

        for i, annotation in enumerate(annotations):
            genome = self.genomes[i]

            for feature in annotation.features:
                if feature.type not in self._jsonDict['data']['features']:
                    self._jsonDict['data']['features'][feature.type] = []

                try:
                    chromosome = genome.get_chromosome_by_name(feature.seqid)
                    new_feature = {
                        'karyo': chromosome.id,
                        'start': feature.start,
                        'end': feature.end
                    }
                    if hasattr(feature, 'name'):
                        new_feature['name'] = feature.name
                    self._jsonDict['data']['features'][feature.type].append(new_feature)
                except KeyError:
                    # Feature is on a chromosome not loaded into AliTV, skip
                    pass

        for annotation_type in annotation_types:
            self._jsonDict['conf']['features']['supportedFeatures'][annotation_type] = {
                'form': 'rect',
                'color': "#E2EDFF",
                'height': 30,
                'visible': True
            }

    def load_tree(self, tree):
        """Convert the tree into AliTV JSON and load it

        :param tree: tree for the AliTV genomes
        :type tree: dendropy.Tree or None
        """

        # process tree to AliTV JSON
        if tree:
            self._jsonDict['data']['tree'] = tree_to_json_structure(tree)
            self._jsonDict['conf']['tree']['drawTree'] = True
            self._jsonDict['filters']['karyo']['genome_order'] = [str(leaf.taxon).strip('"\'') for leaf in
                                                                  tree.leaf_node_iter()]
        else:
            self._jsonDict['filters']['karyo']['genome_order'] = [genome.id for genome in self.genomes]

    def get_json(self, indent):
        """

        :return: AliTV JSON data
        :rtype: str
        """

        return json.dumps(self._jsonDict, indent=indent)
