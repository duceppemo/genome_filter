import os
import sys
import urllib
import subprocess
from glob import glob
from concurrent import futures
import pandas as pd
import urllib.request
import tarfile
import pathlib
from shutil import rmtree
from tqdm import tqdm


class Methods(object):
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
    def make_folder(my_folder):
        pathlib.Path(my_folder).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def get_busco_lineage(lineage, output_folder):
        download_folder = output_folder + '/busco_downloads/'
        lineage_folder = download_folder + 'lineages/'
        Methods.make_folder(lineage_folder)

        if not os.path.exists(download_folder + 'file_versions.tsv'):
            # Get the list of all lineage files
            url = 'https://busco-data.ezlab.org/v5/data/file_versions.tsv'
            try:
                urllib.request.urlretrieve(url, download_folder + 'file_versions.tsv')
            except urllib.error.URLError as e:
                sys.exit('There was a problem downloading the lineage database file: {}'.format(e))

        lineage_file = ''
        lin = ''
        with open(download_folder + 'file_versions.tsv', 'r') as f:
            for line in f:
                if lineage in line:
                    lin, version = line.split('\t')[:2]
                    lineage_file = '.'.join(tuple((lin, version)))
                    break
        if not lineage_file:
            raise Exception('Lineage not found in valid BUSCO lineages.')

        tar_file = lineage_folder + lineage_file + '.tar.gz'
        if not os.path.exists(tar_file):
            url = 'https://busco-data.ezlab.org/v5/data/lineages/' + lineage_file + '.tar.gz'
            try:
                urllib.request.urlretrieve(url, tar_file)
            except urllib.error.URLError as e:
                sys.exit('There was a problem downloading the {} lineage file: {}'.format(lineage, e))

        tar = tarfile.open(tar_file)
        tar.extractall(path=lineage_folder)
        tar.close()
        return lineage_folder + lin

    @staticmethod
    def run_busco(fasta_file, output_folder, cpu, lineage):
        sample_out = os.path.basename(fasta_file).split('.')[0]

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

        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def run_busco_parallel(fasta_list, output_folder, cpu, lineage):
        n_cpu = int(cpu/4)  # run 4 parallel processes
        n_tasks = len(fasta_list)

        with tqdm(total=n_tasks, desc="\t\tBUSCO progress") as pbar:
            with futures.ThreadPoolExecutor(max_workers=n_cpu) as executor:
                args = ((fasta, output_folder, n_cpu, lineage) for fasta in fasta_list)
                for results in executor.map(lambda x: Methods.run_busco(*x), args):
                    pbar.update(1)

        # Cleanup


    @staticmethod
    def process_busco_reports(output_folder):
        # Put path of all reports in a list
        report_list = Methods.get_busco_reports(output_folder)
        # Parse individual busco reports, and write a single report tsv file
        Methods.merge_busco_reports(report_list, output_folder + '/busco_all.tsv')

        # Delete folders
        root, dir_list, file_list = next(os.walk(output_folder))
        for d in dir_list:
            # rmtree(root + '/' + d)
            if d != 'busco_downloads':
                rmtree(root + '/' + d)

    @staticmethod
    def get_busco_reports(report_folder):
        return [y for x in os.walk(report_folder)
                for y in glob(os.path.join(x[0], '*short_summary.specific*.txt'))]

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

        # Rename first column in dataframe
        df.index.name = 'Sample'

        # Write report
        df.to_csv(report_file, sep='\t')

    @staticmethod
    def run_quast(fasta_list, output_folder, cpu):
        cmd = ['quast',
               '--threads', str(cpu),
               '--min-contig', str(0),
               '--output-dir', output_folder,
               '--fast'] + fasta_list

        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def run_checkm2(fasta_list, output_folder, cpu, checkm2_db):
        # Put all input genome in a single folder (symbolic links) with same file extension
        checkm2_fasta_folder = output_folder + 'fasta/'
        checkm2_report_folder = output_folder + 'report/'

        Methods.make_folder(checkm2_fasta_folder)
        Methods.make_folder(checkm2_report_folder)

        for fasta_file in fasta_list:
            try:
                os.symlink(fasta_file, checkm2_fasta_folder + '.'.join(os.path.basename(fasta_file).split('.')[:-1])
                           + '.fasta')
            except FileExistsError:
                pass

        #  The 'general' gradient boost model is able to generalize well and is intended to be used on organisms
        #  not well represented in GenBank or RefSeq (roughly, when an organism is novel at the level of order,
        #  class or phylum).

        # The 'specific' neural network model is more accurate when predicting completeness of organisms more
        # closely related to the reference training set (roughly, when an organism belongs to a known species,
        # genus or family).

        cmd = ['checkm2', 'predict',
               '--force',
               '--remove_intermediates',
               '--threads', str(cpu),
               '--input', checkm2_fasta_folder,
               '--output_directory', checkm2_report_folder,
               '--extension', 'fasta',
               '--specific']
        if checkm2_db:
            cmd += ['--database_path', checkm2_db]

        subprocess.run(cmd)  #, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        # Cleanup
        os.remove(checkm2_report_folder + '/checkm2.log')
        rmtree(checkm2_fasta_folder)
