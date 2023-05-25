localrules: contig_link

rule spades:
    input:
        reads = expand("{dir}/preprocess_done/{{sample}}_{R}.fastq.gz", dir=config['raw_data']['fastq'], R=["R1", "R2"])
    output:
        contigs = "{dir}/{{sample}}/contigs.fasta".format(dir=config['assembly']['output']),
        scaffolds = "{dir}/{{sample}}/scaffolds.fasta".format(dir=config['assembly']['output'])
    log:
        "logs/spades_{sample}.log"
    params:
        extra="--only-assembler"
    threads: 12
    resources:
        partition=config['partition']['highmem'],
        time="30-00:00:00",
        mem_mb=lambda w, input, attempt: min(max((input.size // 1000000) * 10 * (2 + attempt), 100000), 600000)
        # Set the mem as input_size(mb) * 10 * (3 for first try, 4 for second try and 5 for third try) or at least 100G
        # and the maximun usage would not excess 600000 (600G)
    conda:
        "envs/assembler.yaml"
    wrapper:
        "file:workflow/wrappers/metaspades"

rule contig_link:
    input:
        rules.spades.output.contigs
    output:
        "{dir}/{{sample}}_contigs.fasta".format(dir=config['assembly']['filtered_contigs'])
    shell:
        """
        ln -sr {input} {output}
        """

rule filter_contig_length:
    input:
        rules.contig_link.output
    output:
        "{dir}/{{sample}}_filtered.fasta".format(dir=config['assembly']['filtered_contigs'])
    params:
        min_contig_length = config['assembly']['min_contig_length']
    conda:
        "envs/preprocess.yaml"
    shell:
        """
        reformat.sh in={input} out={output} minlength={params.min_contig_length}
        """
