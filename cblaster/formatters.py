"""cblaster result formatters."""


from operator import attrgetter


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


def get_cell_values(queries, subjects, key=len, attr=None):
    """Generates the values of cells in the binary matrix.

    This function calls some specified key function (def. max) against all
    values of a specified attribute (def. None) from Hits inside Subjects which
    match each query. By default, this function will just count all matching Hits
    (i.e. len() is called on all Hits whose query attribute matches). To find
    maximum identities, for example, provide key=max and attr='identity' to this
    function.

    Parameters:
        queries (list): Names of query sequences.
        subjects (list): Subject objects to generate vlaues for.
        key (callable): Some callable that takes a list and produces a value.
        attr (str): A Hit attribute to calculate values with in key function.
    """
    result = [0] * len(queries)
    for index, query in enumerate(queries):
        values = [
            getattr(hit, attr) if attr else hit
            for subject in subjects
            for hit in subject.hits
            if hit.query == query
        ]
        result[index] = key(values)
    return result


def set_decimals(value, decimals=4):
    if isinstance(value, int):
        return str(value)
    return f"{value:.{decimals}f}"


def binary(
    session,
    hide_headers=False,
    delimiter=None,
    key=len,
    attr="identity",
    decimals=4
):
    """Generates a binary summary table from a Session object."""
    value_fn = lambda q, c: get_cell_values(q, c, key=identity, attr=attr)
    rows = [
        [
            organism.full_name,
            accession,
            str(cluster[0].start),
            str(cluster[-1].end),
            *[
                set_decimals(value)
                for value in get_cell_values(
                    session.queries,
                    cluster,
                    key=key,
                    attr=attr
                )
            ]
        ]
        for organism in session.organisms
        for accession, scaffold in organism.scaffolds.items()
        for cluster in scaffold.clusters
    ]
    if not hide_headers:
        rows.insert(0, ["Organism", "Scaffold", "Start", "End", *session.queries])
    if not delimiter:
        delimiter = "  "
        rows = humanise(rows)
    return "\n".join(delimiter.join(row) for row in rows)


def _summarise(
    iterable,
    block_fn,
    header_string,
    header_symbol="-",
    condition_fn=None,
    separator="\n\n",
    **kwargs,
):
    """Skeleton function for generating a summary of a given iterable."""
    blocks = []
    for item in iterable:
        if condition_fn and not condition_fn(item):
            continue
        block = block_fn(item, **kwargs)
        blocks.append(block)
    report = separator.join(blocks)
    if header_string:
        hdr = generate_header_string(header_string, header_symbol)
        return f"{hdr}\n{report}"
    return report


def summarise_organism(organism, hide_headers=True, delimiter=None, decimals=4):
    return _summarise(
        organism.scaffolds.values(),
        summarise_scaffold,
        organism.full_name,
        condition_fn=lambda s: len(s.clusters) > 0,
        header_symbol="=",
        hide_headers=hide_headers,
        delimiter=delimiter,
        decimals=decimals,
        separator="\n\n",
    )


def summarise_scaffold(scaffold, hide_headers=True, delimiter=None, decimals=4):
    return _summarise(
        scaffold.clusters,
        summarise_subjects,
        scaffold.accession,
        hide_headers=hide_headers,
        delimiter=delimiter,
        decimals=decimals,
    )


def summarise_subjects(subjects, decimals=4, hide_headers=True, delimiter=None):
    """Generates a summary table for a hit cluster.

    Args:
        hits (list): collection of Hit objects
        decimals (int): number of decimal points to show
        show_headers (bool): show column headers in output
        human (bool): use human-readable format
    Returns:
        summary table
    """
    subjects.sort(key=attrgetter("start"))
    rows = []
    for subject in subjects:
        values = subject.values(decimals)
        rows.extend(values)
    if not hide_headers:
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
    if not delimiter:
        delimiter = "  "
        rows = humanise(rows)
    return "\n".join(delimiter.join(hit) for hit in rows)


def summary(session, hide_headers=False, delimiter=None, decimals=4):
    return _summarise(
        session.organisms,
        summarise_organism,
        "cblaster search",
        condition_fn=lambda o: o.total_hit_clusters > 0,
        hide_headers=hide_headers,
        delimiter=delimiter,
        decimals=decimals,
        header_symbol="=",
        separator="\n\n\n",
    )


def summarise_gne(data, hide_headers=False, delimiter=None, decimals=4):
    rows = []
    hdrs = ["Gap", "Means", "Medians", "Clusters"]
    if not hide_headers:
        rows.append(hdrs)
    for row in data:
        values = [set_decimals(row.get(key.lower()), decimals) for key in hdrs]
        rows.append(values)
    if not delimiter:
        delimiter = "  "
        rows = humanise(rows)
    return "\n".join(delimiter.join(row) for row in rows)


def gne_summary(data, hide_headers=False, delimiter=None, decimals=4):
    return _summarise(
        data,
        summarise_gne,
        "cblaster gne",
        hide_headers=hide_headers,
        delimiter=delimiter,
        decimals=decimals,
    )
