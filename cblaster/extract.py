"""Extract sequences from session files

TODO: Download hit cluster genomic regions from NCBI
"""


import logging
import re


from cblaster.classes import Session
from cblaster.helpers import efetch_sequences
from cblaster.database import query_sequences


LOG = logging.getLogger(__name__)


def parse_organisms(organisms):
    """Parses specified organisms and creates RegEx patterns."""
    return [re.compile(organism) for organism in organisms]


def organism_matches(organism, patterns):
    """Tests organism filter RegEx patterns against a given organism name."""
    for pattern in patterns:
        if pattern.match(organism):
            return True
    return False


def parse_scaffolds(scaffolds):
    """Parses scaffold names and ranges

    e.g.
        scaf_123 --> {"scaf_123": {"start": None, "end": None}}
        scaf_123:520-62000 --> {"scaf_123": {"start": 520, "end": 62000}}

    Args:
        scaffolds (list): a list of scaffold names with ranges

    Returns:
        A dictionary keyed on scaffolds names containing a dictionary with start and end keys
    """
    records = {}
    for scaffold in scaffolds:
        name, *parts = scaffold.split(":")
        start = None
        end = None
        if parts:
            try:
                start, end = [int(p) for p in parts[0].split("-")]
            except ValueError:
                LOG.exception("Expected range in format start-end")
        records[name] = dict(start=start, end=end)
    return records


def out_of_bounds(subject, start, end):
    """Check if the subject is completely outside of a given range"""
    return subject.start < start or subject.end > end


def record_to_fasta(record, delimiter=None, name_only=False):
    """Formats a given record as FASTA."""
    return ">{}\n{}".format(
        record_to_header(record, delimiter=delimiter, name_only=name_only),
        record.get("sequence"),
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
        func(record, delimiter=delimiter, name_only=name_only) for record in records
    )


def extract_records(
    session, queries=None, organisms=None, scaffolds=None,
):
    """Extracts subject sequence names from a session file given a list of filters

    Args:
        session (Session): cblaster session object
        queries (List):
            a list of query sequences that a subject should match with in order to be included
        organisms (List):
            a list of organism names that a subject should be part of in order to be included
        scaffolds (List):
            a list of scaffold names and ranges that a subject should be part of
            in order to be included
    Returns:
        a List of dictionaries with relevant information from subject objects.
    """
    if organisms:
        organisms = parse_organisms(organisms)
    if scaffolds:
        scaffolds = parse_scaffolds(scaffolds)
    if queries:
        queries = set(queries)
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
            subjects = scaffold.subjects
            for subject in subjects:
                if (start and end) and out_of_bounds(subject, start, end):
                    continue
                if queries and not any(h.query in queries for h in subject.hits):
                    continue
                record = dict(
                    id=subject.id,
                    name=subject.name,
                    organism=organism.name,
                    scaffold=scaffold.accession,
                    start=subject.start,
                    end=subject.end,
                )
                records.append(record)
    return records


def extract_sequences(session, records):
    """Collects sequences for the given list of records.

    Query the offline database for a local mode session or query ncbi for a remote session.
    Sequences are collected into the records value

    Args:
        session (Session): cblaster Session object
        records (List): list of dictionaries with records requested by the user to be extracted
    """
    mode = session.params["mode"]

    if mode == "remote":
        LOG.info("Fetching %i sequences from NCBI", len(records))
        headers = [record.get("name") for record in records]
        sequences = efetch_sequences(headers)
        attr_name = "name"
    elif mode == "local":
        LOG.info("Fetching %i sequences from database", len(records))
        database = session.params["sqlite_db"]
        rowids = [record.get("id") for record in records]
        sequences = {
            rowid: sequence for rowid, sequence in query_sequences(rowids, database)
        }
        attr_name = "id"
    else:
        raise NotImplementedError(
            f"Sessions generated with mode {mode} are not supported yet."
        )
    for record in records:
        record["sequence"] = sequences.get(record[attr_name])


def extract(
    session,
    delimiter=None,
    name_only=False,
    extract_seqs=False,
    output=None,
    queries=None,
    organisms=None,
    scaffolds=None,
):
    """Extract subject sequences from a cblaster session.

    Args:
        session (str): path to json file encoding a cblaster Session object
        extract_seqs (bool): Put the sequences of the extracted proteins into a fasta file
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
        queries=queries,
        organisms=organisms,
        scaffolds=scaffolds,
    )

    if extract_seqs:
        extract_sequences(session, records)

    text = format_records(
        records,
        delimiter=delimiter,
        to_fasta=extract_seqs,
        name_only=name_only,
    )

    if output:
        with open(output, "w") as fp:
            LOG.info("Writing output to %s", fp.name)
            fp.write(text)
    else:
        print(text)

    LOG.info("Done!")
