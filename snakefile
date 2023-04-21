import os
import pandas as pd
import re


class df_prebuild(object):
    def __init__(self, df):
        df = df.copy()
        self.__check(df)
        df = self.__interleave_paired_samples(df)
        self.df = df

    def __check(self, df):
        for idx, v in df.iterrows():
            if bool(v['R1']) ^ bool(v['R2']):
                raise ValueError("One of the paired-end reads is missing (R1 or R2).")
            if not bool(v['R1']) ^ bool(v['interleaved']):
                raise ValueError("R1/R2 and interleaved should be mutually exclusive.")

    def __interleave_paired_samples(self, df):
        paired_idx = df.query('interleaved == ""').index
        df.loc[paired_idx, 'interleaved'] = df.loc[paired_idx, 'sample'] + '_interleaved.fastq.gz'
        return df

    def samples(self):
        return self.df['sample'].drop_duplicates().tolist()

    def filepath_generator(self, sample, path="", column='interleaved'):
        return [os.path.join(path, x) for x in self.df.query('sample == @sample')[column].tolist()]


configfile: "config.yaml"
data = pd.read_csv(config['Metadata'], keep_default_na=False, na_values=['_'], comment="#")
data = df_prebuild(data)

############### Input settings #############
input_list = list()
include: "rules/preprocess.smk"
input_list.extend(["{dir}/preprocess_done/{sample}_R1.fastq.gz".format(dir=config['raw_data']['fastq'], sample=sample) for sample in data.samples()])
input_list.extend(["{dir}/preprocess_done/{sample}_R2.fastq.gz".format(dir=config['raw_data']['fastq'], sample=sample) for sample in data.samples()])

### FastQC
if config['fastqc']['run']:
    input_list.extend(["{dir}/pre_trim/{sample}_fastqc.html".format(dir=config['fastqc']['output'], sample=sample) for sample in data.samples()])
    input_list.extend(["{dir}/pre_trim/{sample}_fastqc.zip".format(dir=config['fastqc']['output'], sample=sample) for sample in data.samples()])

### Fastp and post-trim fastqc
if config['trimming']['run']:
    input_list.extend(["{dir}/trimmed/{sample}.fastq.gz".format(dir=config['raw_data']['fastq'], sample=sample) for sample in data.samples()])
    input_list.extend(["{dir}/trimmed/{sample}.fastq.gz".format(dir=config['raw_data']['fastq'], sample=sample) for sample in data.samples()])
    if config['fastqc']['run']:
        input_list.extend(["{dir}/post_trim/{sample}_fastqc.html".format(dir=config['fastqc']['output'], sample=sample) for sample in data.samples()])
        input_list.extend(["{dir}/post_trim/{sample}_fastqc.zip".format(dir=config['fastqc']['output'], sample=sample) for sample in data.samples()])

### Assembly
if config['assembly']['run']:
    include: "rules/assembly.smk"
    input_list.extend(["{dir}/{sample}_contigs.fasta".format(dir=config['assembly']['filtered_contigs'], sample=sample) for sample in data.samples()])
    input_list.extend(["{dir}/{sample}_filtered.fasta".format(dir=config['assembly']['filtered_contigs'], sample=sample) for sample in data.samples()])  # filtered_fasta

rule all:
    input:
        input_list
