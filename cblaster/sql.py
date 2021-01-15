SCHEMA = """\
CREATE TABLE gene (
    id              INTEGER PRIMARY KEY,
    name            TEXT,
    start_pos       INTEGER,
    end_pos         INTEGER,
    strand          INTEGER,
    translation     TEXT,
    scaffold        TEXT,
    organism        TEXT
);\
"""

ID_QUERY = """\
SELECT
    id,
    name,
    start_pos,
    end_pos,
    strand,
    scaffold,
    organism
FROM
    gene
WHERE
    id IN ({})\
"""

INCLUSIVE_NAME_QUERY = """\
SELECT
    name,
    translation
FROM
    gene
WHERE
    name IN ({}) \
"""

INTERMEDIATE_GENES_QUERY = """\
SELECT
    start_pos,
    end_pos,
    name,
    strand
FROM
    gene
WHERE
    name NOT IN ({}) AND start_pos >= ? AND end_pos <= ?\
"""

FASTA = 'SELECT ">"||gene.id||"\n"||gene.translation||"\n" FROM gene'

INSERT = """\
INSERT INTO gene (
    name,
    start_pos,
    end_pos,
    strand,
    translation,
    scaffold,
    organism
)
VALUES
    (?, ?, ?, ?, ?, ?, ?)\
"""
