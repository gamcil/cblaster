"""
Parser for text format output for searches against ClusteredNR (nr_cluster_seq).
"""


import re
from itertools import pairwise


QUERY_RE = re.compile(r"^Query\s+#(?P<num>\d+):\s*(?P<name>.*?)\s+Query\s+ID:\s*(?P<qid>\S+)\s+Length:\s*(?P<qlen>\d+)", re.M)
CLUSTER_HEADER_RE = re.compile(r"^Cluster:\s*(?P<rep_accession>\S+)\s+(?P<rep_description>.+)$", re.M)
CLUSTER_BODY_RE = re.compile(
	r"""
    ^Cluster:.*\n
	Num\ Members:\s*(?P<num_members>\d+)\n
	Num\ Taxa:\s*(?P<num_taxa>\d+)\n
	Scientific\ Name:\s*(?P<scientific>.+)\n
	Common\ Name\s*:\s*(?P<common>.+)\n
	Taxid:\s*(?P<taxid>\d+)\n
	Highest\ Bit\ Score:\s*(?P<highest_bit>[\d.]+)\n
	Total\ Bit\ Score:\s*(?P<total_bit>[\d.]+)\n
	Percent\ Coverage:\s*(?P<coverage>[\d.]+)%?\n
	Evalue:\s*(?P<evalue>\S+)\n
	Percent\ Identity:\s*(?P<identity>[\d.]+)%?\n
	Accession\ Length:\s*(?P<acc_len>\d+)\n
	(?:\s*\n
	   (?P<mcount>\d+)\s+cluster\s+member\(s\):.*
	   (?P<mheader>Accession[^\n]*)\n
	   (?P<members>(?:.+\n)*?)
	)?
	\Z
    """,
    re.MULTILINE | re.DOTALL | re.VERBOSE
)
MEMBER_HEADER_RE = re.compile(r"^(Accession)\s+(Scientific)\s+(Common)\s+(Taxid)", re.MULTILINE)
ALIGNMENTS_RE = re.compile(r"^Alignments:", re.MULTILINE)


def iter_blocks(pattern, text, meta = lambda m: m.groups()):
    matches = list(pattern.finditer(text))
    if not matches:
        return
    starts = [m.start() for m in matches] + [len(text)]
    for m, (s, e) in zip(matches, pairwise(starts)):
        yield meta(m), text[s:e]


def parse_nr_cluster_text(text):
    for (num, name, qid, qlen), query_block in iter_blocks(
        QUERY_RE,
        text,
        meta=lambda m: (m['num'], m['name'], m['qid'], m['qlen'])
    ):
        # Trim pairwise alignments
        align = ALIGNMENTS_RE.search(query_block)
        if align:
            query_block = query_block[:align.start()]

        for cluster_meta, cluster_block in iter_blocks(CLUSTER_HEADER_RE, query_block, meta=lambda m: (m['rep_accession'], m['rep_description'])):
            # rep_accession, rep_description = cluster_meta

            body = CLUSTER_BODY_RE.search(cluster_block)
            if not body:
                continue

            member_header = body.group("mheader") or ""
            member_table = (body.group("members") or "").strip()
            if not member_header or not member_table:
                continue
            
            # Get field character spacing from headers
            header_match = MEMBER_HEADER_RE.search(member_header)
            if not header_match:
                continue

            col_starts = [header_match.span(num)[0] for num in range(1, 5)] + [len(member_header)]
            identity = body.group("identity").strip()
            coverage = body.group("coverage").strip()
            evalue   = body.group("evalue").strip()
            bitscore = body.group("highest_bit").strip()

            for row in member_table.splitlines():
                cols = [row[s:e] for (s, e) in pairwise(col_starts)]
                if len(cols) < 4:
                    continue
                member_accession, *_ = cols
                values = [name, member_accession.strip(), identity, coverage, evalue, bitscore]
                yield "\t".join(str(v) for v in values)