#!/usr/bin/env python3

"""
This module handles all interaction with NCBI's BLAST API, including launching new
remote searches, polling for completion status, and retrieval of results.
"""

import re
import time
import logging
import requests

from cblaster import helpers
from cblaster.classes import Hit


LOG = logging.getLogger(__name__)

BLAST_API_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"


def _prepare_input(query_file=None, query_ids=None):
    """Prepare query input for BLAST POST request.
    Query file is read into str; ID list is converted to newline delimited str.
    """
    if query_file and not query_ids:
        with open(query_file) as handle:
            return handle.read()
    if query_ids:
        return "\n".join(query_ids)
    raise ValueError("Expected 'query_file' or 'query_ids'")


def start(
    query_file=None,
    query_ids=None,
    database="nr",
    program="blastp",
    megablast=False,
    filtering="F",
    evalue=10,
    nucl_reward=None,
    nucl_penalty=None,
    gap_costs="11 1",
    matrix="BLOSUM62",
    hitlist_size=500,
    threshold=11,
    word_size=6,
    comp_based_stats=2,
    entrez_query=None,
):
    """Launch a remote BLAST search using NCBI BLAST API.

    Note that the HITLIST_SIZE, ALIGNMENTS and DESCRIPTIONS parameters must all be set
    together in order to mimic 'max_target_seqs' behaviour.

    Usage guidelines:

        1. Don't contact server more than once every 10 seconds
        2. Don't poll for a single RID more than once a minute
        3. Use URL parameter email/tool
        4. Run scripts weekends or 9pm-5am Eastern time on weekdays if >50 searches

    For a full description of the parameters, see:

        1. `BLAST API documentation<https://ncbi.github.io/blast-cloud/dev/api.html>`
        2. `BLAST documentation
        <https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp>`

    Parameters
    ----------
    query_file: str
        Path to a query FASTA file
    query_ids: list, tuple
        Collection of NCBI sequence identifiers
    database: str
        Target NCBI BLAST database
    program: str
        BLAST variant to run
    megablast: bool
        Enable megaBLAST option (only with BLASTn)
    filtering: str
        Low complexity filtering
    evalue: float
        E-value cutoff
    nucl_reward: int
        Reward for matching bases (only with BLASTN/megaBLAST)
    nucl_penalty: int
        Penalty for mismatched bases (only with BLASTN/megaBLAST)
    gap_costs: str
        Gap existence and extension costs
    matrix: str
        Scoring matrix name
    hitlist_size: int
        Number of database sequences to keep
    threshold: int
        Neighbouring score for initial words
    word_size: int
        Size of word for initial matches
    comp_based_stats: int
        Composition based statistics algorithm
    entrez_query: str
        NCBI Entrez search term for pre-filtering the BLAST database

    Returns
    -------
    rid: str
        Request Identifier (RID) assigned to the search
    rtoe: int
        Request Time Of Execution (RTOE), estimated run time (in seconds) of the search
    """
    query = _prepare_input(query_file, query_ids)

    parameters = {
        "CMD": "PUT",
        "DATABASE": database,
        "PROGRAM": program,
        "FILTER": filtering,
        "EXPECT": evalue,
        "GAPCOSTS": gap_costs,
        "MATRIX": matrix,
        "HITLIST_SIZE": hitlist_size,
        "ALIGNMENTS": hitlist_size,
        "DESCRIPTIONS": hitlist_size,
        "WORD_SIZE": word_size,
        "COMPOSITION_BASED_STATISTICS": comp_based_stats,
    }

    if entrez_query:
        parameters["ENTREZ_QUERY"] = entrez_query

    if program == "blastn":
        if megablast:
            parameters["MEGABLAST"] = "on"
        if nucl_reward:
            parameters["NUCL_REWARD"] = nucl_reward
        if nucl_penalty:
            parameters["NUCL_PENALTY"] = nucl_penalty
    else:
        # Does not apply to blastn
        parameters["THRESHOLD"] = threshold

    response = requests.post(BLAST_API_URL, files={"QUERY": query}, params=parameters)

    LOG.debug("Search parameters: %s", parameters)
    LOG.debug("Search URL: %s", response.url)

    rid, rtoe = re.findall(r"(?:RID|RTOE) = (.+?)[\n\s]", response.text)
    return rid, int(rtoe)


def check(rid):
    """Check completion status of a BLAST search given a Request Identifier (RID).

    Returns
    -------
    True:
        If search completed successfully and hits were reported
    False:
        Search is still being run

    Raises
    ------
    ValueError
        If the search has failed. This can be caused either by program error (whereby
        NCBI requests you submit an error report with the RID) or expiration of the RID
        (only stored for 24 hours).
    ValueError
        Search has completed successfully, but no hits were reported.
    """
    parameters = {"CMD": "Get", "RID": rid, "FORMAT_OBJECT": "SearchInfo"}

    response = requests.get(BLAST_API_URL, params=parameters)

    LOG.debug(response.url)

    search = re.findall(r"(?:Status|ThereAreHits)=(.+?)[\n\s]", response.text)

    if len(search) == 1:
        status = search[0]
        if status in ("UNKNOWN", "FAILED"):
            raise ValueError(f"Search {rid} failed (status={status})")
        if status == "WAITING":
            return False

    if search == ["READY", "yes"]:
        return True

    raise ValueError("Search completed, but found no hits")


def retrieve(rid, hitlist_size=500):
    """Retrieve BLAST results corresponding to a given Request Identifier (RID).

    Returns
    -------
    list
        List containing BLAST search results, with non-TSV elements (HTML, #info lines)
        removed. Each element in the returned list corresponds to one row in the BLAST
        table (still need to be split by tab).
    """

    parameters = {
        "CMD": "Get",
        "RID": rid,
        "FORMAT_TYPE": "Tabular",
        "FORMAT_OBJECT": "Alignment",
        "HITLIST_SIZE": hitlist_size,
        "ALIGNMENTS": hitlist_size,
        "DESCRIPTIONS": hitlist_size,
        "NCBI_GI": "F",
    }

    LOG.debug(parameters)

    response = requests.get(BLAST_API_URL, params=parameters)

    LOG.debug(response.url)

    # Remove HTML junk and info lines
    # BLAST results are stored inside <PRE></PRE> tags
    return [
        line
        for line in re.search("<PRE>(.+?)</PRE>", response.text, re.DOTALL)
        .group(1)
        .split("\n")
        if line and not line.startswith("#")
    ]


def poll(rid, delay=60, max_retries=-1):
    """Poll BLAST API with given Request Identifier (RID) until results are returned.

    As per NCBI usage guidelines, this function will only poll once per minute; this is
    calculated each time such that wait is constant (i.e. accounts for differing
    response time on the status check).

    Returns
    -------
    list:
        Output of retrieve()
    """
    if delay < 60:
        raise ValueError("Delay must be at least 60s")

    retries, previous = 0, 0
    while True:
        current = time.time()
        wait = previous - current + delay
        if wait > 0:
            time.sleep(wait)
            previous = current + wait
        else:
            previous = current

        LOG.info("Checking search status...")

        if check(rid):
            LOG.info("Search has completed successfully!")
            return

        if max_retries > 0 and retries == max_retries:
            raise ValueError(f"Reached maximum retry limit {max_retries}")

        retries += 1


def parse(
    handle,
    query_file=None,
    query_ids=None,
    max_evalue=0.01,
    min_identity=30,
    min_coverage=50,
):
    """Parse Tabular results from remote BLAST search performed via API.

    Since the API provides no option for returning query coverage, which is a metric we
    want to use for filtering hits, query sequences must be passed to this function so
    that their lengths can be compared to the alignment length.

    Parameters
    ----------
    handle: open file handle
        File handle (or file handle-like) object corresponding to BLAST results. Note
        that this function expects an iterable of tab-delimited lines and performs no
        validation/error checking
    query_file: str
        Path to query file
    query_ids: list, tuple
        Collection of NCBI sequence identifiers
    max_evalue: float
        Maximum e-value
    min_identity: float
        Minimum percent identity
    min_coverage: float
        Minimum percent query coverage

    Returns
    -------
    list
        List of Hit instances corresponding to criteria passing BLAST hits
    """
    sequences = helpers.get_sequences(query_file, query_ids)

    hits = []
    for line in handle:
        qid, sid, pident, *_, qstart, qend, _, _, evalue, score, _ = line.split("\t")

        # Manually calculate query coverage
        coverage = (int(qend) - int(qstart) + 1) / len(sequences[qid]) * 100

        hit = Hit(
            query=qid,
            subject=sid,
            identity=pident,
            coverage=coverage,
            evalue=evalue,
            bitscore=score,
        )

        if (
            hit.identity > min_identity
            and hit.coverage > min_coverage
            and hit.evalue < max_evalue
        ):
            hits.append(hit)

    if len(hits) == 0:
        raise ValueError("No results found")

    return hits


def search(
    rid=None,
    query_file=None,
    query_ids=None,
    min_identity=0.3,
    min_coverage=0.5,
    max_evalue=0.01,
    **kwargs,
):
    """Perform a remote BLAST search via the NCBI's BLAST API.

    This function launches a new search given a query FASTA file or list of valid NCBI
    identifiers, polls the API to check the completion status of the search, then
    retrieves and parses the results.

    It is also possible to call other BLAST variants using the `program` argument.

    Parameters
    ----------
    rid: str
        Request Identifier (RID) of a web BLAST search
    query_file: str
        Path to query FASTA file
    query_ids: list, tuple
        Collection of NCBI sequence identifiers
    min_identity: float
        Minimum percent identity
    min_coverage: float
        Minimum percent query coverage
    max_evalue: float
        Maximum e-value

    Returns
    -------
    list
        List of Organism instances, containing Scaffold instances that store Hits and
        clusters of Hits.
    """
    if not rid:
        LOG.info("Launching new search")

        # Start search, get request identifier (RID) and completion ETA (RTOE)
        rid, rtoe = start(query_file=query_file, query_ids=query_ids, **kwargs)

        LOG.info("Request Identifier (RID): %s", rid)
        LOG.info("Request Time Of Execution (RTOE): %ss", rtoe)

        # Wait the RTOE (sec) before bothering to poll
        time.sleep(rtoe)

        LOG.info("Polling NCBI for completion status")
        poll(rid)

    LOG.info("Retrieving results for search %s", rid)
    results = retrieve(rid)

    # Parse results for hits
    LOG.info("Parsing results...")
    results = parse(
        results,
        query_file=query_file,
        query_ids=query_ids,
        max_evalue=max_evalue,
        min_identity=min_identity,
        min_coverage=min_coverage,
    )

    return rid, results
