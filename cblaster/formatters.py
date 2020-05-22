"""
cblaster result formatters.
"""


import builtins


def get_maximum_row_lengths(rows):
    """Finds the longest lengths of fields per column in a collection of rows."""
    lengths, total = [], len(rows[0])
    for index in range(total):
        largest = max(len(row[index]) for row in rows)
        lengths.append(largest)
    return lengths


def add_field_whitespace(rows, lengths):
    """Fills table fields with whitespace to specified lengths."""
    result = []
    for row in rows:
        fmt = [f"{row[index]:{length}}" for index, length in enumerate(lengths)]
        result.append(fmt)
    return result


def humanise(rows):
    """Formats a collection of fields as human-readable."""
    lengths = get_maximum_row_lengths(rows)
    table = add_field_whitespace(rows, lengths)
    return table


def generate_header_string(text, symbol="-"):
    """Generates a 2-line header string with underlined text.

    >>> header = generate_header_string("header string", symbol="*")
    >>> print(header)
    header string
    *************
    """
    return f"{text}\n{symbol * len(text)}"


def count_query_hits(queries, hits):
    """Counts total hits per query in a colllection of Hit objects.

    >>> queries = ["query1", "query2", "query3"]
    >>> hits = [Hit(query="query1"), Hit(query="query2"), Hit(query="other")]
    >>> count_query_hits(queries, hits)
    [1, 1, 0]

    Args:
        hits (list): Hit objects
    Returns:
        List of per-query counts corresponding to input.
    """
    return [sum(query == hit.query for hit in hits) for query in queries]


def get_hit_identities(queries, hits, key=max):
    """Get the maximum hit identity per query in a collection of `Hit` objects.

    Argument passed to key must be a callable, e.g. max or sum.

    >>> queries = ["query1", "query2", "query3"]
    >>> hits = [
    ...     Hit(query="query1", identity=0.9),
    ...     Hit(query="query2", identity=0.4),
    ...     Hit(query="query2", identity=0.7),
    ... ]
    >>> get_hit_identities(queries, hits, key=max)
    [0.9, 0.7, 0]
    >>> get_hit_identities(queries, hits, key=sum)
    [0.9, 1.1, 0]
    """
    return [
        key([hit.identity if query == hit.query else 0 for hit in hits])
        for query in queries
    ]


def binary(session, headers=False, human=False, identity=None, delimiter=","):
    """Generates a binary summary table from a Session object."""
    if identity:
        value_fn = lambda q, c: get_hit_identities(q, c, key=identity)
    else:
        value_fn = count_query_hits
    rows = [
        [
            organism.full_name,
            accession,
            str(cluster[0].start),
            str(cluster[-1].end),
            *[str(value) for value in value_fn(session.queries, cluster)]
        ]
        for organism in session.organisms
        for accession, scaffold in organism.scaffolds.items()
        for cluster in scaffold.clusters
    ]
    if headers:
        rows.insert(0, ["Organism", "Scaffold", "Start", "End", *session.queries])
    if human:
        delimiter = "  "
        rows = humanise(rows)
    return "\n".join(delimiter.join(row) for row in rows)


def _summarise(
    iterable,
    block_fn,
    header_string,
    header_symbol="-",
    condition_fn=None,
    **kwargs,
):
    blocks = []
    for item in iterable:
        if condition_fn and not condition_fn(item):
            continue
        block = block_fn(item, **kwargs)
        blocks.append(block)
    report = "\n\n".join(blocks)
    if header_string:
        hdr = generate_header_string(header_string, header_symbol)
        return f"{hdr}\n{report}"
    return report


def summarise_organism(organism, headers=True, human=True, decimals=4, delimiter=","):
    return _summarise(
        organism.scaffolds.values(),
        summarise_scaffold,
        organism.full_name,
        condition_fn=lambda s: len(s.clusters) > 0,
        header_symbol="=",
        headers=headers,
        human=human,
        decimals=decimals,
        delimiter=delimiter,
    )


def summarise_scaffold(scaffold, headers=True, human=True, decimals=4, delimiter=","):
    return _summarise(
        scaffold.clusters,
        summarise_cluster,
        scaffold.accession,
        headers=headers,
        human=human,
        decimals=decimals,
        delimiter=delimiter,
    )


def summarise_cluster(hits, decimals=4, headers=True, human=True, delimiter=","):
    """Generates a summary table for a hit cluster.

    Args:
        hits (list): collection of Hit objects
        decimals (int): number of decimal points to show
        show_headers (bool): show column headers in output
        human (bool): use human-readable format
    Returns:
        summary table
    """
    rows = [h.values(decimals) for h in hits]
    if headers:
        hdrs = [
            "Query",
            "Subject",
            "Identity",
            "Coverage",
            "E-value",
            "Bitscore",
            "Start",
            "End",
            "Strand",
        ]
        rows.insert(0, hdrs)
    if human:
        delimiter = "  "
        rows = humanise(rows)
    return "\n".join(delimiter.join(hit) for hit in rows)


def summary(session, headers=True, human=True, decimals=4, delimiter=","):
    return _summarise(
        session.organisms,
        summarise_organism,
        "cblaster search",
        condition_fn=lambda o: o.total_hit_clusters > 0,
        headers=headers,
        human=human,
        decimals=decimals,
        delimiter=delimiter,
        header_symbol="=",
    )
