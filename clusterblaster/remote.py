#!/usr/bin/env python3

"""
This module handles all interaction with NCBI's BLAST API, including launching new
remote searches, polling for completion status, and retrieval of results.
"""

import re
import time
import logging
import requests

from clusterblaster import helpers
from clusterblaster.classes import Hit


logging.basicConfig(level=logging.INFO)
LOG = logging.getLogger(__name__)

BLAST_API_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"


def start(
    query_file=None,
    query_ids=None,
    database="nr",
    program="blastp",
    megablast=False,
    filtering="F",
    evalue=0.01,
    nucl_reward=None,
    nucl_penalty=None,
    gap_costs="11 1",
    matrix="BLOSUM62",
    hitlist_size=0,
    threshold=11,
    word_size=6,
    comp_based_stats=2,
    entrez_query=None,
):
    """Launch a remote BLAST search using NCBI BLAST API.

    Usage guidelines:
        1. Don't contact server more than once every 10 seconds
        2. Don't poll for a single RID more than once a minute
        3. Use URL parameter email/tool
        4. Run scripts weekends or 9pm-5am Eastern time on weekdays if >50 searches

    Common URL API: https://ncbi.github.io/blast-cloud/dev/api.html

    Further BLAST documentation:
    https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp
    https://www.ncbi.nlm.nih.gov/books/NBK279684/
    """
    parameters = {
        "CMD": "PUT",
        "DATABASE": database,
        "PROGRAM": program,
        "FILTER": filtering,
        "EXPECT": evalue,
        "GAPCOSTS": gap_costs,
        "MATRIX": matrix,
        "HITLIST_SIZE": hitlist_size,
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

    LOG.debug(parameters)

    if query_file and not query_ids:
        with open(query_file) as handle:
            query = handle.read()
    elif query_ids:
        query = "\n".join(query_ids)

    response = requests.post(BLAST_API_URL, files={"QUERY": query}, params=parameters)

    LOG.debug(response.url)

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
    else:
        if search == ["READY", "yes"]:
            return True
        raise ValueError("Search completed, but found no hits")


def retrieve(rid):
    """Retrieve BLAST results corresponding to a given Request Identifier (RID)."""
    parameters = {
        "CMD": "Get",
        "RID": rid,
        "FORMAT_TYPE": "Tabular",
        "FORMAT_OBJECT": "Alignment",
        "HITLIST_SIZE": 0,
        "NCBI_GI": "F",
    }

    LOG.debug(parameters)

    response = requests.get(BLAST_API_URL, params=parameters)

    LOG.debug(response.url)

    # Remove non-TSV junk; 15:-3 removes HTML, then parse out info lines (#)
    return [
        line
        for line in response.text.split("\n")[15:-3]
        if line and not line.startswith("#")
    ]


def poll(rid):
    """Poll BLAST API with given Request Identifier (RID) until results are returned."""
    previous = 0
    while True:
        current = time.time()
        wait = previous - current + 60
        if wait > 0:
            time.sleep(wait)
            previous = current + wait
        else:
            previous = current

        LOG.info("Checking search status...")
        if check(rid):
            LOG.info("Search has completed successfully")
            break

    LOG.info("Retrieving results from NCBI")
    return retrieve(rid)


def parse(
    handle,
    query_file=None,
    query_ids=None,
    max_evalue=0.01,
    min_identity=0.3,
    min_coverage=0.5,
):
    """Parse Tabular results from remote BLAST search performed via API.

    Since the API provides no option for returning query coverage, which is a metric we
    want to use for filtering hits, query sequences must be passed to this function so
    that their lengths can be compared to the alignment length.
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
        raise SystemExit("No results found")

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
    """
    if not rid:
        LOG.info("Launching new remote BLAST search")

        # Start search, get request identifier (RID) and completion ETA (RTOE)
        rid, rtoe = start(query_file=query_file, query_ids=query_ids, **kwargs)
        LOG.debug("RID: %s, RTOE: %s", rid, rtoe)

        # Wait the RTOE (sec) before bothering to poll
        time.sleep(rtoe)

        LOG.info("Polling NCBI for completion status of search %s", rid)
        results = poll(rid)
    else:
        LOG.info("Retrieving results for search %s", rid)
        results = retrieve(rid)

    # Parse results for hits
    LOG.info("Parsing results")
    results = parse(
        results,
        query_file=query_file,
        query_ids=query_ids,
        max_evalue=max_evalue,
        min_identity=min_identity,
        min_coverage=min_coverage,
    )

    return results
