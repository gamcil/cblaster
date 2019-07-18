#!/usr/bin/env python3

"""
Find species containing query proteins in same genomic region.

Process:
    1. Read in query FASTA file
    2. Search against BLASTp/diamond database
    3. Parse hits for each query that meet thresholds
    4. Look up hits on NCBI for genomic context
    5. Iterate scaffolds, report hits if enough genes are conserved (-c),
       printing spaces between hits if gaps is above a threshold (-g)

Similar to MultiGeneBlast, but:
    1. Uses kblastp/diamond which are way quicker
    2. Performs no synteny analysis, instead querying either
       funcoDB/NCBI for genomic context and doing basic text
       based analysis to determine if genes are close
       (again, way quicker)

Cameron Gilchrist
2019-03-27
"""

import re
import shutil
import subprocess

import click


class Organism:
    """ Store hits on scaffolds, and functionality for finding clusters.
    """
    def __init__(self, name, strain):
        self.name = name
        self.strain = strain
        self.scaffolds = {}

    def __repr__(self):
        separator = '-' * len(self.full_name)
        scaffolds = '\n'.join([
            str(scaffold) + '\n' for scaffold in self.scaffolds.values()
            if scaffold.clusters
        ])
        if scaffolds:
            return f'{self.full_name}\n{separator}\n{scaffolds}\n'
        return f'{self.full_name} -- No hits\n'

    def sort_scaffolds(self):
        for scaffold in self.scaffolds.values():
            scaffold.sort_scaffold()

    def find_clusters(self, conserve, gap):
        for scaffold in self.scaffolds.values():
            scaffold.find_clusters(conserve, gap)

    @property
    def full_name(self):
        if self.strain in self.name or not self.strain:
            return f'{self.name}'
        return f'{self.name} {self.strain}'


class Scaffold:
    """ Store hits on Scaffolds and method for finding clusters.
    """
    def __init__(self, accession, hits=None):
        self.accession = accession
        self.hits = hits if hits else []
        self.clusters = []

    def __repr__(self):
        if self.clusters:
            clusters = '\n'.join(
                str(hit)
                for start, end in self.clusters
                for hit in self.hits[start:end]
            )
            seperator = '-' * len(self.accession)
            return f'{self.accession}\n{seperator}\n{clusters}'
        # TODO: refactor to reporting function so can
        #       toggle showing single hits?
        # if self.hits:
        #     hits = '\n'.join(str(hit) for hit in self.hits)
        #     return (f'{self.accession} -- Only single hits found:\n'
        #             f'{hits}\n')
        return ''

    @property
    def protein_uids(self):
        return ','.join(hit.target for hit in self.hits)

    def sort_scaffold(self):
        self.hits.sort(key=lambda hit: hit.start, reverse=False)

    def find_clusters(self, conserve, gap):
        """ Find clustered hit groups on this Scaffold.
        """
        total_hits = len(self.hits)
        if total_hits < conserve:
            return None
        i = 0
        while i < total_hits:
            reset, last = False, i
            for j in range(i + 1, total_hits):
                if self.hits[j].end - self.hits[j - 1].start > gap:
                    reset, last = True, j - 1  # Gap is too big
                elif j == total_hits - 1:
                    reset, last = True, j  # At the last element
                if reset:
                    if last + 1 - i >= conserve:
                        self.clusters.append((i, last + 1))  # Save indices
                    break
            i = last + 1  # Move index ahead of last block


class Hit:
    """ Handle BLAST hits.
    """
    def __init__(self, query, subject, identity, coverage, evalue):
        self.query = query

        if 'gb' in subject or 'ref' in subject:
            subject = re.search(r'\|([A-Za-z0-9\._]+)\|',
                                subject).group(1)
        self.subject = subject
        self.identity = float(identity)
        self.coverage = float(coverage)
        self.evalue = float(evalue)
        self.start = None
        self.end = None

    def __repr__(self):
        return (f'{self.query}\t'
                f'{self.subject}\t'
                f'{self.identity}\t'
                f'{self.coverage}\t'
                f'{self.start}-{self.end}\t'
                f'{self.evalue}')


@click.command()
@click.argument('query')  # help='Path to FASTA file with query sequences')
@click.argument('database')  # help='Path to BLAST database being searched')
@click.option('-g', '--gap',
              default=20000,
              help='Maximum gap length (bp) between hits')
@click.option('-c', '--conserve',
              default=3,
              help='Number of query sequences that must be conserved in hits')
@click.option('-p', '--program',
              default='diamond',
              type=click.Choice(['diamond', 'blastp']),
              help='Program to use in search (def. "diamond")')
@click.option('-e', '--evalue',
              default=0.001,
              help='E-value cutoff (def. <= 0.01)')
@click.option('-d', '--pident',
              default=30.00,
              help='Percent identity cutoff (def. >= 30.00)')
@click.option('-v', '--qcovhsp',
              default=50.00,
              help='HSP coverage cutoff (def. >= 50.00)')
@click.option('-c', '--cpus',
              default=1,
              help='Number of CPUs to use')
def clusterblaster(
    query, database, gap, conserve, program,
    evalue, pident, qcovhsp, cpus
):
    """ Run clusterblaster.
    """
    # Run the specified search program
    click.echo('Running BLAST against your database')
    if program == 'diamond':
        results = diamond(query, database, evalue, pident, qcovhsp,
                          cpus)
    else:
        results = blastp(query, database, evalue, qcovhsp, cpus)

    # Parse the results and save
    hits = parse_blast(results, pident, qcovhsp, evalue)

    if not hits:
        raise SystemExit('No hits found')

    organisms = ncbi_genomic_context(hits)
    for organism in organisms.values():
        organism.find_clusters(conserve, gap)
        print(organism)


def blastp(query, database, evalue, qcovhsp, cpus):
    """ Search database using BLASTp.
    """
    cmd = ' '.join([
        'blastp',
        '-query', query,
        '-db', database,
        '-evalue', str(evalue),
        '-qcov_hsp_perc', str(qcovhsp),
        '-num_threads', str(cpus),
        '-subject_besthit',
        '-max_hsps', '1',
        '-outfmt', '6 qseqid sseqid pident qcovhsp evalue',
    ])

    print(f'Running BLASTp with:\n\n\t{cmd}')
    return subprocess.check_output(cmd, shell=True).decode('utf-8')


def diamond(query, database, evalue, pident, qcovhsp, cpus):
    """ Search database using DIAMOND.
    """
    executable = shutil.which('diamond')
    if not executable:  # apt installed
        executable = shutil.which('diamond-aligner')

    if not executable:
        raise ValueError('Could not find diamond on your system PATH')

    cmd = ' '.join([
        executable,
        'blastp',
        '--query', query,
        '--db', database,
        '--id', str(pident),
        '--evalue', str(evalue),
        '--query-cover', str(qcovhsp),
        '--outfmt', '6 qseqid sseqid pident qcovhsp evalue',
        '--threads', str(cpus)
    ])

    click.echo(f'Running diamond with:\n\n{cmd}\n')
    return subprocess.check_output(cmd, shell=True).decode('utf-8')


def parse_blast(results, pident, qcovhsp, evalue):
    """Parse results from run_blast(), return hit sets for each query
    """
    hits = []
    for row in results.split('\n')[:-1]:
        hit = Hit(*row.split('\t'))  # handles type conversions
        if hit.identity > pident and \
                hit.coverage > qcovhsp and \
                hit.evalue < evalue:
            hits.append(hit)

    # If no hits for anything, quit
    if len(hits) == 0:
        raise SystemExit('No results found')
    return hits


def parse_fasta(fasta):
    """ Read in protein sequences from a FASTA file, return dict."""
    name, records = '', {}
    with open(fasta, 'r') as handle:
        for line in handle:
            line = line.strip()
            if line.startswith('>'):
                name = line.replace('>', '')
                records[name] = ''
            else:
                records[name] += line
    return records


def ncbi_genomic_context(hits):
    """ Query NCBI for genomic context of a list of proteins.

        Uses the identical proteins format to find nucleotide record
        corresponding to the matched protein.
    """
    hitmap = {hit.subject: hit for hit in hits}

    efetch = shutil.which('efetch')
    if not efetch:
        raise ValueError('Could not find efetch on your system PATH')

    cmd = [
        efetch,
        '-db', 'protein',
        '-id', ','.join(hit.subject for hit in hits),
        '-format', 'ipg'
    ]

    click.echo(f'Finding genomic context with:\n\n{cmd}\n')
    result = subprocess.run(cmd, encoding='utf-8',
                            check=True,
                            stdout=subprocess.PIPE)

    # Skip header and last (empty) line
    organisms = {}
    previous = ''
    for line in result.stdout.strip().split('\n')[1:-1]:
        try:
            ipg, _, accession, start, end, \
                _, protein, _, organism, \
                strain, _ = line.split('\t')
        except ValueError:
            # Catch <Error> lines
            continue

        if ipg == previous:
            continue
        previous = ipg

        if organism not in organisms:
            organisms[organism] = Organism(name=organism, strain=strain)

        if accession not in organisms[organism].scaffolds:
            organisms[organism].scaffolds[accession] = Scaffold(accession)

        try:
            hit = hitmap.pop(protein)
        except KeyError:
            continue

        hit.start = int(start)
        hit.end = int(end)
        organisms[organism].scaffolds[accession].hits.append(hit)

    # Sort scaffolds in organisms by hit locations
    for organism in organisms.values():
        organism.sort_scaffolds()
    return organisms


if __name__ == '__main__':
    clusterblaster()
