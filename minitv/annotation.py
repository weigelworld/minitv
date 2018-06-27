"""
Objects for parsing and representing annotation file formats
"""


class GFF3Feature:
    """
    Represents a single GFF3 feature.
    """

    def __init__(self, seqid, source, feature_type, start, end, score, strand, phase, attributes):
        """Converts split GFF3 line into python native types and stores them in the object

        :param seqid: "ID of the landmark used to establish the coordinate system for the current feature"
        :type seqid: str
        :param source: "free text qualifier [describing] the algorithm or [procedure] that generated this feature"
        :type source: str
        :param feature_type: type of feature, restricted to Sequence Ontology or SO accession number
        :type feature_type: str
        :param start: positive 1-based integer coordinate, relative to seqid; always < END
        :type start: str
        :param end: positive 1-based integer coordinate, relative to seqid, always > START
        :type end: str
        :param score: floating point number, E-value (sequences) or P-value (gene predictors)
        :type score: str
        :param strand: + for positive strand (relative to landmark), - for minus, ? for unknown, . for none
        :type strand: str
        :param phase: where feature begins with reference to reading phrase: 0, 1, or 2
        :type phase: str
        :param attributes: dictionary of tag, value pairs
        :type attributes: dict
        """

        self.seqid = seqid
        self.source = source
        self.type = feature_type
        self.start = int(start)
        self.end = int(end)
        self.score = float(score) if score != '.' else None
        self.strand = strand
        self.phase = int(phase) if phase != '.' else None
        self.attributes = attributes

        if 'ID' in attributes:
            self.ID = attributes['ID']
        if 'Name' in attributes:
            self.name = attributes['Name']


class GFF3Annotation:
    """
    Represents a GFF3 annotation file.
    """

    def __init__(self, path, types=None):
        """

        :param path: path to the GFF3 annotation file
        :type path: os.PathLike
        :param types: Optional feature types to filter the annotation by
        :type types: collections.Iterable[str]
        """

        self.features = []

        with open(path) as input_GFF3:
            for line in input_GFF3:
                if not line.startswith('#'):
                    seqid, source, feature_type, start, end, score, strand, phase, attributes = line.rstrip().split(
                        "\t")

                    #attributes = dict([key_value_pair.split('=') for key_value_pair in attributes.split(';')])
                    # FIXME: some annotation software (*ahem* CAT) doesn't output strict GFF3 and allows unescaped
                    # semi-colons in attribute values, this is a sloppy fix for this, discarding anything after the
                    # first semi-colon in a value
                    attributes = [key_value_pair.split('=') for key_value_pair in attributes.split(';')]
                    attributes = dict([ele for ele in attributes if len(ele) == 2])
                    # END FIXME

                    if (types and feature_type in types) or (not types):
                        feature = GFF3Feature(seqid, source, feature_type, start, end, score, strand, phase, attributes)
                        self.features.append(feature)
