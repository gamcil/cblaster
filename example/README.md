## Example cblaster search

In this folder you will find cblaster output generated from searching
the burnettramic acids gene cluster from *Aspergillus burnettii*, *bua*,
remotely against the NCBI's NR database.

| File | Description |
| ---- | ----------- |
| ``bua.gbk`` | *bua* cluster GenBank file |
| ``bua.fasta`` | FASTA file containing proteins in *bua* |
| ``bua_session.json`` | cblaster session file |
| ``bua_binary.txt`` | Binary absence/presence table output |
| ``bua_summary.txt`` | Default results summary |
| ``bua_results.html`` | Interactive output visualisation |

These files were generated with the following command:

	cblaster search \
		--query_file bua.gbk \
		--session bua_session.json \
		--binary bua_binary.txt \
		--output bua_summary.txt \
		--plot bua_results.html

Searching with the GenBank file allows for synteny scoring of detected clusters.
However, you could also search with the FASTA file by specifying:

	--query_file bua.fasta
