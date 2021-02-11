# Creates the feature table containing genes and scaffolds
SCHEMA = """\
CREATE TABLE feature (
    id              INTEGER PRIMARY KEY,
    feature_type    TEXT,
    name            TEXT,
    start_pos       INTEGER,
    end_pos         INTEGER,
    strand          INTEGER,
    sequence        TEXT,
    scaffold        TEXT,
    organism        TEXT
);\
"""

INSERT = """\
INSERT INTO feature (
    name,
    feature_type,
    start_pos,
    end_pos,
    strand,
    sequence,
    scaffold,
    organism
)
VALUES
    (?, ?, ?, ?, ?, ?, ?, ?)\
"""

# Base for querying scaffolds in feature table
SCAFFOLD_QUERY = """\
SELECT
    substr(sequence, ?, ?)
FROM
    feature
WHERE
    feature_type = \"scaffold\"
    AND name = ?
    AND organism = ?\
"""

GENE_QUERY = """\
SELECT
    id,
    name,
    start_pos,
    end_pos,
    strand,
    scaffold,
    organism
FROM
    feature
WHERE
    feature_type = \"gene\"
    AND id IN ({})\
"""

# Query genes between some range not in a list
INTERMEDIATE_GENES_QUERY = """\
SELECT
    start_pos,
    end_pos,
    name,
    strand
FROM
    feature
WHERE
    feature_type = \"gene\"
    AND name NOT IN ({})
    AND scaffold = ?
    AND organism = ?
    AND start_pos >= ? AND end_pos <= ?\
"""

# Write database sequences in FASTA format
FASTA = 'SELECT ">"||gene.id||"\n"||gene.translation||"\n" FROM gene'
