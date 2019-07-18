# clusterblaster
Find clustered hits from a BLAST search

## Process
This tool leverages the identical proteins functionality of NCBI entrez to
fetch the genomic context of BLAST hit accessions, then reports groups of
hits that are co-located in their respective genomes.

## Dependencies
`clusterblaster` requires `diamond` (alternatively, `diamond-aligner` if installed via
apt) to be on your system `$PATH`.

The only external Python dependency is `click`, which is used for the command-line
interface.

## Usage
`$ python3 clusterblaster.py --help`

```
Usage: clusterblaster.py [OPTIONS] QUERY DATABASE

  Run clusterblaster.

Options:
  -g, --gap INTEGER               Maximum gap length (bp) between hits
  -c, --conserve INTEGER          Number of query sequences that must be
                                  conserved in hits
  -p, --program [diamond|blastp]  Program to use in search (def. "diamond")
  -e, --evalue FLOAT              E-value cutoff (def. <= 0.01)
  -d, --pident FLOAT              Percent identity cutoff (def. >= 30.00)
  -v, --qcovhsp FLOAT             HSP coverage cutoff (def. >= 50.00)
  -c, --cpus INTEGER              Number of CPUs to use
  --help                          Show this message and exit.
```

## Example
1. Build your diamond database
```
diamond makedb database.fasta
```

2. Search your query sequences against the database using `clusterblaster`.
```
python3 clusterblaster.py query.fasta database.dmnd \
            -g 20000 -c 3 -p diamond -e 0.05 -v 70.00 -c 3
```

3. Enjoy your output
```
...
Aspergillus arachidicola CBS 117610
-----------------------------------
NEXV01000673.1
--------------
AAS90113.1	PIG79658.1	95.4	100.0	190063-191663	4.8e-288
AAS90109.1	PIG79657.1	93.6	100.0	192013-192813	1.7e-143
AAS90108.1	PIG79656.1	86.9	97.6	192914-194891	8.4e-308
AAS90107.1	PIG79655.1	98.2	100.0	195444-197092	1.9e-300
AAS90106.1	PIG79682.1	96.9	100.0	197660-199644	0.0
AAS90105.1	PIG79681.1	96.0	100.0	200160-202105	1.3e-310
AAS90104.1	PIG79684.1	94.1	100.0	203424-204896	8.6e-230
AAS90103.1	PIG79683.1	95.3	98.7	206202-208342	8e-214
AAS90102.1	PIG79678.1	94.6	85.4	208741-213083	1.4e-231
AAS90100.1	PIG79677.1	95.0	93.9	213307-220990	1.8e-263
AAS90097.1	PIG79676.1	96.2	95.7	221773-222747	2.8e-144
AAS90096.1	PIG79675.1	98.4	100.0	223327-224790	4.3e-240
AAS90091.1	PIG79621.1	96.8	100.0	225531-226865	2.3e-225
AAS90095.1	PIG79622.1	96.1	100.0	228311-234156	0.0
AAS90111.1	PIG79619.1	96.8	100.0	234524-239981	0.0
AAS90110.1	PIG79620.1	97.4	98.9	241288-243323	9.4e-145
AAS90093.1	PIG79625.1	97.5	100.0	243964-250591	0.0
AAS90092.1	PIG79626.1	96.5	100.0	252282-254262	8.6e-283

Aspergillus parasiticus SU-1
----------------------------
JZEE01000728.1
--------------
AAS90092.1	KJK60792.1	97.4	100.0	4401-8760	1.5e-282
AAS90093.1	KJK60793.1	99.1	100.0	10427-17052	0.0
AAS90110.1	KJK60791.1	97.0	100.0	17695-19731	5.5e-145
AAS90111.1	KJK60794.1	97.0	100.0	21009-26450	0.0
AAS90095.1	KJK60796.1	97.3	100.0	26829-32674	0.0
JZEE01000729.1
--------------
AAS90096.1	KJK60754.1	97.9	100.0	1414-2857	2.9e-236
AAS90098.1	KJK60750.1	80.8	74.2	3425-6023	1.1e-104
AAS90094.1	KJK60773.1	93.5	95.9	6468-7674	4.9e-203
AAS90099.1	KJK60755.1	99.2	100.0	8568-9467	4.1e-145
AAS90100.1	KJK60752.1	94.7	99.8	10327-12030	1.3e-280
AAS90102.1	KJK60771.1	97.6	98.6	13082-16606	1.7e-285
AAS90104.1	KJK60770.1	96.6	99.3	18633-22769	2.4e-232
AAS90105.1	KJK60756.1	94.2	100.0	24090-26035	6.7e-304
AAS90106.1	KJK60751.1	97.5	100.0	26547-28531	0.0
AAS90107.1	KJK60772.1	98.0	100.0	29098-30746	2.5e-300
AAS90108.1	KJK60753.1	92.9	100.0	31301-33114	0.0
AAS90109.1	KJK60763.1	92.9	100.0	33376-34176	3.3e-142
AAS90113.1	KJK60766.1	93.3	99.0	34536-36128	3.5e-278
...
```
The above is an excerpt from a search of aflatoxin cluster protein sequences
against other Aspergillus genomes.

Output is formatted as follows (tab-separated):
```
Organism name
-------------
Scaffold
--------
query   subject     identity (%)   coverage (%)    start - end   e-value
...
```
