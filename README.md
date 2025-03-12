# cblaster
[![Python package](https://github.com/gamcil/cblaster/actions/workflows/pythonapp.yml/badge.svg)](https://github.com/gamcil/cblaster/actions/workflows/pythonapp.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyPI version](https://badge.fury.io/py/cblaster.svg)](https://badge.fury.io/py/cblaster)
[![Documentation Status](https://readthedocs.org/projects/cblaster/badge/?version=latest)](https://cblaster.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3660769.svg)](https://doi.org/10.5281/zenodo.3660769)

**>>> Both cblaster and clinker can now be used without installation on the [CAGECAT webserver](http://cagecat.bioinformatics.nl/). <<<**

**>>> The locally installable version of cblaster now also integrates [`cagecleaner`](https://github.com/LucoDevro/cagecleaner), a tool to reduce genomic redundancy among cblaster hits. <<<**

## Outline

`cblaster` is a tool for finding clusters of co-located homologous sequences
in BLAST searches.

<img src="docs/source/_static/workflow.png" alt="cblaster search workflow" width=600>

Given a collection of protein sequences, `cblaster` can search sequence databases
remotely (via NCBI BLAST API) or locally (via `DIAMOND`). Search results are parsed
and filtered based on user thresholds for identity, coverage and e-value. The genomic
coordinates of remaining hits are obtained from the NCBI's Identical Protein
Group (IPG) database (or a local database in local searches). Finally,
`cblaster` scans for instances of collocation and generates visualisations:

<img src="docs/source/_static/results.png" alt="cblaster search results" width=700>

## Installation
It is highly recommended to install `cblaster` in a fresh conda environment using the yaml file in this repo:

```bash
$ conda env create -f env.yml
```

Make sure you have installed the NCBI Entrez-Direct Utilities as well. This can be done as described [here](https://www.ncbi.nlm.nih.gov/books/NBK179288/).

Once installed, make sure you configure cblaster with your email address:

```bash
$ cblaster config --email name@domain.com
```

You can find example search files, along with generated output, in the [examples folder
of the repository](https://github.com/gamcil/cblaster/tree/master/example).
For example output of the `cagecleaner` dereplication tool, head over to the [examples folder
in its own repo](https://github.com/LucoDevro/cagecleaner/tree/main/examples).

## Dependencies
All dependencies should have been managed by conda if you chose to use a conda environment,
except for the NCBI Entrez-Direct utilities.

`cblaster` itself is tested on Python 3.6, and its only external Python dependency are
the `requests` module (used for interaction with NCBI APIs), and `cagecleaner` for 
dereplicating hits.
If you want to perform local searches, you should have `diamond` installed and available
on your system $PATH.
`cblaster` will throw an error if a local search is started but it cannot find
`diamond` or `diamond-aligner` (alias when installed via apt) on the system.

`cagecleaner` is tested on Python 3.10, and comes with its own dependencies. These should
have been installed by conda. Head over to the [cagecleaner wiki](https://github.com/LucoDevro/cagecleaner/wiki/Software-requirements) to find a complete list if you need to install these manually.

## Usage
`cblaster` accepts FASTA files and collections of valid NCBI sequence identifiers
(GIs, accession numbers) as input.
A remote search can be performed as simply as:

```bash
$ cblaster search --query_file query.fasta
```

For example, to remotely search the
[burnettramic acids gene cluster, *bua*](https://pubs.acs.org/doi/10.1021/acs.orglett.8b04042)
, against the NCBI's nr database:

```bash
$ cblaster search -qf bua.fasta

[12:14:17] INFO - Starting cblaster in remote mode
[12:14:17] INFO - Launching new search
[12:14:19] INFO - Request Identifier (RID): WHS0UGYJ015
[12:14:19] INFO - Request Time Of Execution (RTOE): 25s
[12:14:44] INFO - Polling NCBI for completion status
[12:14:44] INFO - Checking search status...
[12:15:44] INFO - Checking search status...
[12:16:44] INFO - Checking search status...
[12:16:46] INFO - Search has completed successfully!
[12:16:46] INFO - Retrieving results for search WHS0UGYJ015
[12:16:51] INFO - Parsing results...
[12:16:51] INFO - Found 3944 hits meeting score thresholds
[12:16:51] INFO - Fetching genomic context of hits
[12:17:14] INFO - Searching for clustered hits across 705 organisms
[12:17:14] INFO - Writing summary to <stdout>

Aspergillus mulundensis DSM 5745
================================
NW_020797889.1
--------------
Query       Subject         Identity  Coverage  E-value    Bitscore  Start    End      Strand
QBE85641.1  XP_026607259.1  75.56     99.5918   0          742       1717881  1719409  -
QBE85642.1  XP_026607260.1  89.916    100       0          667       1719650  1720797  +
QBE85643.1  XP_026607261.1  89.532    83.1169   0          832       1721494  1722934  +
QBE85644.1  XP_026607262.1  64.829    98.9218   6.51e-157  455       1723252  1724467  -
QBE85645.1  XP_026607263.1  69.97     100       6.93e-157  449       1725113  1726277  -
QBE85646.1  XP_026607264.1  82.759    96.8447   0          670       1726892  1728302  +
QBE85647.1  XP_026607265.1  72.674    99.2048   0          764       1729735  1731338  +
QBE85648.1  XP_026607266.1  56.098    98.324    4.24e-64   205       1731701  1732402  -
QBE85649.1  XP_026607267.1  79.623    99.8746   0          6573      1732820  1745289  +

...
```

A query sequence absence/presence matrix can be generated using the `--binary` argument:

```
Organism                                   Scaffold        Start    End      QBE85641.1  QBE85642.1  QBE85643.1  QBE85644.1  QBE85645.1  QBE85646.1  QBE85647.1  QBE85648.1  QBE85649.1
Aspergillus mulundensis DSM 5745           NW_020797889.1  1717881  1745289  1           1           1           1           1           1           1           1           1         
Aspergillus versicolor CBS 583.65          KV878126.1      3162095  3187090  1           1           1           0           1           1           1           1           1         
Pseudomassariella vexata CBS 129021        MCFJ01000004.1  1606356  1628483  1           1           1           0           0           1           0           1           1         
Hypoxylon sp. CO27-5                       KZ112517.1      92119    112957   1           1           1           0           0           0           1           0           1         
Hypoxylon sp. EC38                         KZ111255.1      514739   535366   1           1           1           0           0           0           1           0           1         
Epicoccum nigrum ICMP 19927                KZ107839.1      2116719  2142558  1           1           0           0           0           1           1           0           1         
Aureobasidium subglaciale EXF-2481         NW_013566983.1  700476   718693   1           1           0           0           0           1           1           0           0         
Aureobasidium pullulans EXF-6514           QZBF01000009.1  18721    34295    1           1           0           0           0           1           1           0           0         
Aureobasidium pullulans EXF-5628           QZBI01000512.1  329      13401    1           0           0           0           0           1           1           0           0         
```

`cblaster` can also generate fully interactive visualisations of the binary
table. To view an example, click [here](https://cblaster.readthedocs.io/en/latest/_static/example.html).

For further usage examples and API documentation, please refer to the
[documentation](https://cblaster.readthedocs.io/en/latest/).

## Citation
If you found this tool useful, please cite:

```text
Cameron L M Gilchrist, Thomas J Booth, Bram van Wersch, Liana van Grieken, Marnix H Medema, Yit-Heng Chooi, cblaster: a remote search tool for rapid identification and visualisation of homologous gene clusters, Bioinformatics Advances, 2021;, vbab016, https://doi.org/10.1093/bioadv/vbab016
```

as well as

```text
De Vrieze, L., Biltjes, M., Lukashevich, S., Tsurumi, K., Masschelein, J. (2025) CAGEcleaner: reducing genomic redundancy in gene cluster mining. bioRxiv https://doi.org/10.1101/2025.02.19.639057
```

`cblaster` makes use of the following tools:
```
Buchfink, B., Xie, C. & Huson, D. H. Fast and sensitive protein alignment using DIAMOND. Nat. Methods 12, 59–60 (2015).

Acland, A. et al. Database resources of the National Center for Biotechnology Information. Nucleic Acids Res. 42, 7–17 (2014).
```

`cagecleaner` makes use of the following:
```
Salamzade, R., & Kalan, L. R. (2023). skDER: microbial genome dereplication approaches for comparative and metagenomic applications. https://doi.org/10.1101/2023.09.27.559801`
Shaw, J., & Yu, Y. W. (2023). Fast and robust metagenomic sequence comparison through sparse chaining with skani. Nature Methods, 20(11), 1661–1665. https://doi.org/10.1038/s41592-023-02018-3
```
