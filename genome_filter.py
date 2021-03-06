#!/usr/local/bin python3

import os
import sys
import subprocess
from glob import glob
import pandas as pd
from multiprocessing import cpu_count
from argparse import ArgumentParser
from concurrent import futures
import urllib.request
import tarfile
import pathlib
from shutil import rmtree


# TODO -> add check for valid lineages?


class GenomeFilter(object):
    def __init__(self, args):
        # Paths
        self.input_folder = args.input
        self.output_folder = args.output
        self.lineage = args.lineage

        # Performance
        self.cpu = args.threads
        if self.cpu > cpu_count():
            self.cpu = cpu_count()

        # Variables
        self.fasta_list = list()
        self.report_list = list()
        self.busco_env_path = ''
        self.quast_env_path = ''

        # Run
        self.run()

    def run(self):
        if not GenomeFilter.is_conda_installed():
            raise Exception('Could not find a valid installation of conda. Please install conda.\n'
                            'https://docs.conda.io/projects/conda/en/latest/user-guide/install/')

        self.fasta_list = GenomeFilter.get_fasta_files(self.input_folder)

        # BUSCO
        busco_path = ''
        if not GenomeFilter.is_conda_env_installed('busco'):
            GenomeFilter.install_busco_env()
            self.busco_env_path = GenomeFilter.get_conda_env_path('busco')
        # if not GenomeFilter.is_env_activated('busco'):
        #     GenomeFilter.conda_activate('busco')

        lineage = GenomeFilter.get_lineage(self.lineage, self.output_folder)
        self.run_busco_parallel(self.fasta_list, self.output_folder, self.cpu, lineage, 'busco')
        self.process_busco_reports()
        # GenomeFilter.conda_deactivate()

        # QUAST
        quast_path = ''
        if not GenomeFilter.is_conda_env_installed('quast'):
            GenomeFilter.install_quast_env()
            self.quast_env_path = GenomeFilter.get_conda_env_path('quast')
        # if not GenomeFilter.is_env_activated('quast'):
        #     GenomeFilter.conda_activate('quast')
        self.run_quast(self.fasta_list, self.output_folder, self.cpu, 'quast')
        # GenomeFilter.conda_deactivate()

    # conda
    @staticmethod
    def is_conda_installed():
        return os.path.exists(os.path.join(sys.prefix, 'conda-meta'))

    @staticmethod
    def is_conda_env_installed(env):
        cmd = ['conda', 'info', '--envs']
        return env in subprocess.check_output(cmd).decode()

    @staticmethod
    def is_env_activated(env):
        return env in os.environ['CONDA_DEFAULT_ENV']

    @staticmethod
    def get_conda_env_path(env):
        cmd = ['conda', 'info', '--envs']
        env_dict = dict()
        env_out = subprocess.check_output(cmd).decode().split('\n')
        for line in env_out:
            if line.startswith('#') or line == '':
                continue
            e, p = line.split(' ', maxsplit=1)
            e = e.strip()
            p = p.strip()
            env_dict[e] = p

        return env_dict[env]

    # BUSCO
    @staticmethod
    def install_busco_env():
        cmd = ['conda', 'create', '-y', '-n', 'busco', '-c', 'bioconda', '-c' 'anaconda', 'busco', 'pandas']
        subprocess.run(cmd)

    @staticmethod
    def install_quast_env():
        cmd = ['conda', 'create', '-y', '-n', 'quast', '-c', 'bioconda', 'quast']
        subprocess.run(cmd)

    @staticmethod
    def conda_activate(env):
        cmd = ['conda', 'activate', env]
        subprocess.run(cmd)

    @staticmethod
    def conda_deactivate():
        cmd = ['conda', 'deactivate']
        subprocess.run(cmd)

    @staticmethod
    def get_lineage(lineage, output_folder):
        pathlib.Path(output_folder + '/busco_downloads/lineages').mkdir(parents=True, exist_ok=True)

        if not os.path.exists(output_folder + '/busco_downloads/file_versions.tsv'):
            # Get the list of all lineage files
            url = 'https://busco-data.ezlab.org/v4/data/file_versions.tsv'
            urllib.request.urlretrieve(url, output_folder + '/busco_downloads/file_versions.tsv')

        lineage_file = ''
        lin = ''
        version = ''
        with open(output_folder + '/busco_downloads/file_versions.tsv', 'r') as f:
            for line in f:
                if lineage in line:
                    lin, version = line.split('\t')[:2]
                    lineage_file = '.'.join(tuple((lin, version)))
                    break
        if not lineage_file:
            raise Exception('Lineage not found in valid BUSCO lineages.')

        tar_file = output_folder + '/busco_downloads/lineages/' + lineage_file + '.tar.gz'
        if not os.path.exists(tar_file):
            url = 'https://busco-data.ezlab.org/v4/data/lineages/' + lineage_file + '.tar.gz'
            urllib.request.urlretrieve(url, tar_file)

        tar = tarfile.open(tar_file)
        tar.extractall(path=output_folder + '/busco_downloads/lineages/')
        tar.close()
        return output_folder + '/busco_downloads/lineages/' + lin

    @staticmethod
    def run_busco(fasta_file, output_folder, cpu, lineage, env):
        sample_out = os.path.basename(fasta_file).split('.')[0]
        conda_run = ['conda', 'run', '-n', env]
        cmd = ['busco',
               '--force',
               '-i', fasta_file,
               '-m', 'genome',
               '--out_path', output_folder,
               '-o', sample_out,
               '-l', lineage,
               '--offline',
               '--quiet',
               '-c', str(cpu)]
        subprocess.run(conda_run + cmd)

    @staticmethod
    def run_busco_parallel(fasta_list, output_folder, cpu, lineage, env):
        n_cpu = int(cpu/4)  # run 4 parallel processes
        with futures.ThreadPoolExecutor(max_workers=n_cpu) as executor:
            args = ((fasta, output_folder, n_cpu, lineage, env) for fasta in fasta_list)
            for results in executor.map(lambda x: GenomeFilter.run_busco(*x), args):
                pass
        # for fasta_file in fasta_list:
        #     GenomeFilter.run_busco(fasta_file, output_folder, cpu, lineage)

    def process_busco_reports(self):
        # Put path of all reports in a list
        self.report_list = GenomeFilter.get_busco_reports(self.output_folder)
        # Parse individual busco reports, and write a single report tsv file
        GenomeFilter.merge_busco_reports(self.report_list, self.output_folder + '/busco_all.tsv')
        # Delete folders
        root, dir_list, file_list = next(os.walk(self.output_folder))
        for d in dir_list:
            if d != 'busco_downloads':
                rmtree(root + '/' + d)

    @staticmethod
    def get_fasta_files(fasta_folder):
        fasta_list = list()
        ext_tuple = ('.fasta', '.fna', '.fa')
        for root, dirs, files in os.walk(fasta_folder):
            for file in files:
                if file.endswith(ext_tuple):
                    fasta_list.append(root + os.sep + file)
        return fasta_list

    @staticmethod
    def get_busco_reports(report_folder):
        return [y for x in os.walk(report_folder)
                for y in glob(os.path.join(x[0], '*short_summary.specific*'))]

    # Parse report info into dictionary
    """
    # Create consolidated report
        # C:99.2%[S:99.2%,D:0.0%],F:0.3%,M:0.5%,n:639    
        # 634 Complete BUSCOs (C)            
        # 634 Complete and single-copy BUSCOs (S)    
        # 0   Complete and duplicated BUSCOs (D)     
        # 2   Fragmented BUSCOs (F)              
        # 3   Missing BUSCOs (M)             
        # 639 Total BUSCO groups searched
    """
    @staticmethod
    def merge_busco_reports(report_list, report_file):
        busco_dict = dict()
        for report in report_list:
            with open(report, 'r') as f:
                sample = os.path.basename(report).split('.')[-2]
                for line in f:
                    if 'C:' in line:
                        line = line.strip()
                        info_list = line.split('%')
                        completed = info_list[0].split(':')[1]
                        comp_single = info_list[1].split(':')[1]
                        com_dup = info_list[2].split(':')[1]
                        frag = info_list[3].split(':')[1]
                        miss = info_list[4].split(':')[1]
                        total = info_list[5].split(':')[1]
                        busco_dict[sample] = [completed, comp_single, com_dup, frag, miss, total]

        # Convert dictionary to pandas data frame
        df = pd.DataFrame.from_dict(busco_dict, orient='index',
                                    columns=['Completed', 'Completed_single_copy', 'Completed_duplicated',
                                             'Fragmented', 'Missing', 'Total'])
        # Write report
        df.to_csv(report_file, sep='\t')

    @staticmethod
    def run_quast(fasta_list, output_folder, cpu, env):
        quast_out = output_folder + '/quast'
        pathlib.Path(quast_out).mkdir(parents=True, exist_ok=True)
        conda_run = ['conda', 'run', '-n', env, ]
        cmd = ['quast',
               '--threads', str(cpu),
               '--min-contig', str(0),
               '--output-dir', quast_out,
               '--fast'] + fasta_list
        subprocess.run(conda_run + cmd)


if __name__ == "__main__":
    max_cpu = cpu_count()

    parser = ArgumentParser(description='Create dendrogram from a square distance matrix')
    parser.add_argument('-i', '--input', metavar='/input/folder',
                        required=True,
                        help='Input folder containing the genomes in fasta format')
    parser.add_argument('-o', '--output', metavar='/output/folder',
                        required=True,
                        help='Folder to hold the result files')
    parser.add_argument('-l', '--lineage', metavar='rhizobiales_odb10',
                        required=True,
                        help='See "https://busco.ezlab.org/list_of_lineages.html" to help find the correct lineage.')
    parser.add_argument('-t', '--threads', metavar=str(max_cpu),
                        required=False,
                        type=int, default=max_cpu,
                        help='Number of CPU. Default is maximum CPU available({})'.format(max_cpu))

    # Get the arguments into an object
    arguments = parser.parse_args()

    GenomeFilter(arguments)
