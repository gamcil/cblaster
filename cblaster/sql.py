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

# Insert into the feature table
INSERT = """\
INSERT INTO feature (
    feature_type,
    name,
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

SCAFFOLD_QUERY = """\
SELECT
    substr(sequence, {}, {})
FROM
    feature
WHERE
    feature_type = 'scaffold'
    AND name = ?
    AND organism = ?\
"""


# Retrieve all information about genes
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

# Retrieve sequence of features in a list
SEQUENCE_QUERY = """\
SELECT
    id,
    sequence
FROM
    feature
WHERE
    id IN ({})\
"""

# Query genes between some range not in a list
INTERMEDIATE_GENES_QUERY = """\
SELECT
    start_pos,
    end_pos,
    id,
    name,
    strand
FROM
    feature
WHERE
    feature_type = \"gene\"
    AND id NOT IN ({})
    AND scaffold = ?
    AND organism = ?
    AND start_pos >= ? AND end_pos <= ?\
"""

# Write database sequences in FASTA format
FASTA = 'SELECT ">"||feature.id||"\n"||feature.sequence||"\n" FROM feature WHERE feature_type = \"gene\"'
