#!/usr/bin/env python3
import os.path
from multiprocessing import cpu_count
from argparse import ArgumentParser
from genome_filter_methods import Methods

__author__ = 'duceppemo'
__version__ = 'v0.2'


# TODO
# Make CheckM work with gzipped files

class GenomeFilter(object):
    def __init__(self, args):
        # Paths
        self.input = os.path.abspath(args.input)
        self.output_folder = os.path.abspath(args.output)
        self.lineage = args.lineage
        self.checkm2_db = os.path.abspath(args.database)

        # Performance
        self.cpu = args.threads
        if self.cpu > cpu_count():
            self.cpu = cpu_count()

        # Run
        self.run()

    def run(self):
        # Output folder structure
        busco_folder = self.output_folder + '/busco/'
        quast_folder = self.output_folder + '/quast/'
        checkm2_folder = self.output_folder + '/checkm2/'

        Methods.make_folder(busco_folder)
        Methods.make_folder(quast_folder)
        Methods.make_folder(checkm2_folder)

        # Let's get to work
        if os.path.isdir(self.input):
            print('Getting fasta files...')
            fasta_list = Methods.get_fasta_files(self.input)
            n_fasta = len(fasta_list)
            print('\tFound {} genomes.'.format(n_fasta))
        elif os.path.isfile(self.input):
            fasta_list = [self.input]
        else:
            raise.Exception('Error with input folder / file.')

        # BUSCO
        print('Running BUSCO...')
        print('\tDownloading lineage database')
        lineage = Methods.get_busco_lineage(self.lineage, busco_folder)

        print('\tParallel processing genomes')
        Methods.run_busco_parallel(fasta_list, busco_folder, self.cpu, lineage)

        print('\tPreparing report')
        Methods.process_busco_reports(busco_folder)

        # QUAST
        print('Running QUAST...')
        Methods.run_quast(fasta_list, quast_folder, self.cpu)

        # CheckM2
        print('Running CheckM2...')
        Methods.run_checkm2(fasta_list, checkm2_folder, self.cpu, self.checkm2_db)

        print('Done!')


if __name__ == "__main__":
    max_cpu = cpu_count()

    parser = ArgumentParser(description='Assess genome assemblies for completeness and contamination.')
    parser.add_argument('-i', '--input', metavar='/input/folder or /input/file.fasta',
                        required=True, type=str,
                        help='Input fasta file or input folder containing the genomes in fasta format. Mandatory.')
    parser.add_argument('-o', '--output', metavar='/output/folder',
                        required=True,
                        help='Folder to hold the result files. Mandatory.')
    parser.add_argument('-l', '--lineage', metavar='rhizobiales_odb12',
                        required=True, type=str,
                        help='Run "busco --list-datasets" in the terminal for a complete list of lineages.')
    parser.add_argument('-d', '--database', metavar='/path/to/checkm2_db/uniref100.KO.1.dmnd',
                        required=False,
                        type=str,
                        help='Path to the "uniref100.KO.1.dmnd" CheckM2 database, in case you did not '
                             'previously run "checkm2 database --setdblocation". Optional.')
    parser.add_argument('-t', '--threads', metavar=str(max_cpu),
                        required=False, type=int, default=max_cpu,
                        help='Number of CPU. Optional.\nDefault is maximum CPU available ({}).'.format(max_cpu))

    # Get the arguments into an object
    arguments = parser.parse_args()

    GenomeFilter(arguments)
