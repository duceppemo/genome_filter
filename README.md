# genome_filter
Filter downloaded genomes using BUSCO, QUAST and CheckM

## Usage
```
usage: python genome_filter.py [-h] -i /input/folder -o /output/folder -l rhizobiales_odb10 [-t 64] -r
                        {life,domain,phylum,class,order,family,genus,species} -n "Bacillus cereus"

Asses genome assemblies for completeness and contamination

optional arguments:
  -h, --help            show this help message and exit
  -i /input/folder, --input /input/folder
                        Input folder containing the genomes in fasta format.
  -o /output/folder, --output /output/folder
                        Folder to hold the result files
  -l rhizobiales_odb10, --lineage rhizobiales_odb10
                        See "https://busco.ezlab.org/list_of_lineages.html" to help find the correct lineage.
  -t 64, --threads 64   Number of CPU. Default is maximum CPU available(64)
  -r {life,domain,phylum,class,order,family,genus,species}, --rank {life,domain,phylum,class,order,family,genus,species}
                        Taxonomic level to use wiht CheckM.
  -n "Bacillus cereus", --name "Bacillus cereus"
                        Name of the taxon. E.g. "Bacillus cereus" if "-r" was "species". You have to used the double quotes
                        for "species" becasue there is a space between the genus and species.
```
