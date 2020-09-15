"""Extract sequences from session files

TODO: Download hit cluster genomic regions from NCBI
"""


import logging
import re


from cblaster.classes import Session
from cblaster.helpers import efetch_sequences


LOG = logging.getLogger(__name__)


def parse_organisms(organisms):
    """Parses specified organisms and creates RegEx patterns."""
    return [
        re.compile(organism)
        for organism in organisms
    ]


def organism_matches(organism, patterns):
    """Tests organism filter RegEx patterns against a given organism name."""
    for pattern in patterns:
        if pattern.match(organism):
            return True
    return False


def parse_scaffolds(scaffolds):
    """Parses scaffold names and ranges

    e.g.
        scaf_123 --> {"name": "scaf_123"}
        scaf_123:520-62000 --> {"name": "scaf_123", "start": 520, "end": 62000}
    """
    records = {}
    for scaffold in scaffolds:
        name, *parts = scaffold.split(":")
        start = None
        end = None
        if parts:
            try:
                start, end = [int(p) for p in parts.split("-")]
            except ValueError:
                LOG.exception("Expected range in format start-end")
        records[name] = dict(start=start, end=name)
    return records


def out_of_bounds(subject, start, end):
    """Tests if a subject overlaps with or is outside of a given range."""
    return (
        subject.end < start
        or subject.start < start < subject.end
        or subject.start < end < subject.end
        or subject.start > end
    )


def flatten(array):
    """Flattens a list of lists.
    e.g. [[1, 2, 3], [4, 5, 6]] --> [1, 2, 3, 4, 5, 6]
    """
    flat = []
    for element in array:
        flat.extend(element)
    return flat


def record_to_fasta(record, delimiter=None, name_only=False):
    """Formats a given record as FASTA."""
    return ">{}\n{}".format(
        record_to_header(record, delimiter=delimiter, name_only=name_only),
        record.get("sequence")
    )


def record_to_header(record, delimiter=None, name_only=False):
    """Builds a header for a given record."""
    if name_only:
        return record["name"]
    fields = ["name", "organism", "scaffold", "start", "end"]
    values = [str(record[field]) for field in fields]
    if delimiter:
        return delimiter.join(values)
    return "{} [organism={}] [scaffold={}:{}-{}]".format(*values)


def format_records(records, delimiter=None, to_fasta=False, name_only=False):
    """Formats records """
    func = record_to_fasta if to_fasta else record_to_header
    return "\n".join(
        func(record, delimiter=delimiter, name_only=name_only)
        for record in records
    )


def extract_records(
    session,
    in_cluster=True,
    queries=None,
    organisms=None,
    scaffolds=None,
):
    """Extracts subject sequence names from a session file."""
    if organisms:
        organisms = parse_organisms(organisms)
    if scaffolds:
        scaffolds = parse_scaffolds(scaffolds)
    records = []
    for organism in session.organisms:
        if organisms and not organism_matches(organism.name, organisms):
            continue
        for accession, scaffold in organism.scaffolds.items():
            if scaffolds:
                if accession not in scaffolds:
                    continue
                record = scaffolds[accession]
                start = record.get("start")
                end = record.get("end")
            else:
                start = None
                end = None
            if in_cluster:
                subjects = flatten(cluster for cluster in scaffold.clusters)
            else:
                subjects = scaffold.subjects
            for subject in subjects:
                if (start and end) and out_of_bounds(subject, start, end):
                    continue
                if queries and not any(h.query in queries for h in subject.hits):
                    continue
                record = dict(
                    name=subject.name,
                    organism=organism.name,
                    scaffold=scaffold.accession,
                    start=subject.start,
                    end=subject.end,
                )
                records.append(record)
    return records


def extract(
    session,
    in_cluster=True,
    delimiter=None,
    name_only=False,
    download=False,
    output=None,
    queries=None,
    organisms=None,
    scaffolds=None,
):
    """Extract subject sequences from a cblaster session.

    Parameters:
        session (Session): cblaster Session object
        download (bool): Download hit sequences from NCBI
        output (str): Output file name
        queries (list): Query sequence names
        organisms (list): Organism filtering regular expressions
        scaffolds (list): Scaffold names and ranges
        delimiter (str): Sequence description delimiter character
        name_only (bool): Do not save sequence descriptions
    """
    LOG.info("Starting cblaster extraction")
    LOG.info("Loading session from: %s", session)
    with open(session) as fp:
        session = Session.from_json(fp)

    LOG.info("Extracting subject sequences matching filters")
    records = extract_records(
        session,
        in_cluster=in_cluster,
        queries=queries,
        organisms=organisms,
        scaffolds=scaffolds,
    )

    if download:
        LOG.info("Fetching %i sequences from NCBI", len(records))
        headers = [record.get("name") for record in records]
        sequences = efetch_sequences(headers)
        for record in records:
            record["sequence"] = sequences.get(record["name"])

    # FASTA format if downloading from NCBI, otherwise newline separated IDs
    text = format_records(
        records,
        delimiter=delimiter,
        to_fasta=download,
        name_only=name_only,
    )

    if output:
        with open(output, "w") as fp:
            LOG.info("Writing output to %s", fp.name)
            fp.write(text)
    else:
        print(text)

    LOG.info("Done!")
    return records
