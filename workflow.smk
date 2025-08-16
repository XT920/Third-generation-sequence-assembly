# ""
# name: Third generation sequence assembly
# # date: 2023-04-24
# description: Snakemake pipeline to do Third generation sequence assembly
# author: Tong Xu
# dependencies:
#     - fastp = 0.23.2
#     - flye = 2.9.2
#     - minimap2 = 2.24
#     - samtools = 1.17
#     - pilon = 1.24

# rules:
#     - fastp
#     - flye
#     - minimap2
#     - pilon


# ""

import os
from os import path
import pandas as pd
from Bio import SeqIO


base_dir = "/home/map/snakemake/TG"
pilon_jar = "/home/xjm/software/pilon-1.24.jar"
SG_reads_dir = path.join(base_dir, "SG_reads")
TG_reads_dir = path.join(base_dir, "TG_reads")
SAMPLE = os.listdir(SG_reads_dir)

FASTP_OUT_PATH = "fastp_output"
flye_OUT_PATH = "flye_output"
minimap2_OUT_PATH = "minimap2_flye_output"
pilon_OUT_PATH = "pilon_flye_output"
final_genome_PATH = "TG_genome"
unfixed_genome_PATH = "TG_unfixed_genome"
prodigal_OUT_PATH = "prodigal_output"

shell.executable('/bin/zsh') # set a shell to run commands


rule all:
    input:
        expand(FASTP_OUT_PATH + "/{sample}/{sample}_1.fq", sample=SAMPLE),
        expand(FASTP_OUT_PATH + "/{sample}/{sample}_2.fq", sample=SAMPLE),
        expand(FASTP_OUT_PATH + "/{sample}/{sample}_fastp.json", sample=SAMPLE),
        expand(FASTP_OUT_PATH + "/{sample}/{sample}_fastp.html", sample=SAMPLE),
        expand(flye_OUT_PATH + "/{sample}/assembly.fasta", sample=SAMPLE),
        expand(minimap2_OUT_PATH + "/{sample}/{sample}_aln.bam", sample=SAMPLE),
        expand(pilon_OUT_PATH + "/{sample}/{sample}_pilon.fasta", sample=SAMPLE),
        expand(minimap2_OUT_PATH + "/{sample}/{sample}_aln_2.bam", sample=SAMPLE),
        expand(pilon_OUT_PATH + "/{sample}/{sample}_pilon_2.fasta",sample=SAMPLE),
        expand(minimap2_OUT_PATH + "/{sample}/{sample}_aln_3.bam", sample=SAMPLE),
        expand(pilon_OUT_PATH + "/{sample}/{sample}_pilon_3.fasta",sample=SAMPLE),
        expand(unfixed_genome_PATH + "/{sample}_unfixed.fasta",sample=SAMPLE),
        expand(prodigal_OUT_PATH + "/{sample}_unfixed.gff",sample=SAMPLE),
        expand(final_genome_PATH + "/{sample}.fasta",sample=SAMPLE),




rule fastp:
    input:
        r1 = SG_reads_dir + "/{sample}/{sample}_1.fq",
        r2 = SG_reads_dir + "/{sample}/{sample}_2.fq",
    output:
        r1 = FASTP_OUT_PATH + "/{sample}/{sample}_1.fq",
        r2 = FASTP_OUT_PATH + "/{sample}/{sample}_2.fq",
        json = ensure(FASTP_OUT_PATH + "/{sample}/{sample}_fastp.json", non_empty=True), # If the command somecommand happens to generate an empty output, the job will fail with an error listing the unexpected empty file.
        html = ensure(FASTP_OUT_PATH + "/{sample}/{sample}_fastp.html", non_empty=True), # Often, it is a good idea to combine ensure annotations with retry definitions, e.g. for retrying upon invalid checksums or empty files.
    log:
        FASTP_OUT_PATH + "/{sample}/{sample}_fastp.log"
    threads: 10
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} "
        "-j {output.json} -h {output.html} -w {threads} "
        "-q 20 " # the quality value that a base is qualified
        "-u 10 " # how many percents of bases are allowed to be unqualified (0~100) Default 40 means 40%.
        "-n 3 " # number of N allowed
        "-l 40 " # reads shorter than length_required will be discarded, particularly those with adapter.
        "-5 " # enable trimming in 5' ends.
        "-3 " # enable trimming in 3' ends.
        "&> {log}"


# mapping
rule flye:
    params:
        genomeSize = "4.5m", #不同样本可能具体匹配, flye可以不指定这个 对于Phaeobacter gallaeciensis 4.5m的组装效果比4.7m好
        outdir = flye_OUT_PATH + "/{sample}",
        asm_coverage = 300 #限制测序深度，测序深度过高会组装失败！Using longest 300x reads for contig assembly
    input:
        TG_reads_dir + "/{sample}/{sample}.filtered_reads.fq"
    output:
        flye_OUT_PATH + "/{sample}/assembly.fasta",
    threads: 96
    shell:
        "flye --nano-raw {input} -g {params.genomeSize} --asm-coverage {params.asm_coverage} --out-dir {params.outdir} --threads {threads}"



# 基于三代数据的自我校正 racon，大部分文章没有做这一步

# Correct the assembly

rule minimap2:
    input:
        genome=rules.flye.output,
        r1 = rules.fastp.output.r1,
        r2 = rules.fastp.output.r2,
    output:
        minimap2_OUT_PATH + "/{sample}/{sample}_aln.bam"
    threads: 10
    shell:
        "minimap2 -ax sr {input.genome} {input.r1} {input.r2} | samtools sort -@ {threads} -O bam -o {output};"
        "samtools index -@ {threads} {output};"

rule pilon:
    input:
        genome = rules.flye.output,
        bam = rules.minimap2.output
    output:
        pilon_OUT_PATH + "/{sample}/{sample}_pilon.fasta"
    params:
        prefix = "{sample}_pilon",
        outdir= pilon_OUT_PATH + "/{sample}",
        pilonJar = pilon_jar
    shell:
        "java -Xmx16G -jar {params.pilonJar} --genome {input.genome} --frags {input.bam} --output {params.prefix} --fix all --mindepth 10 --outdir {params.outdir} --changes --verbose"

# 只有大于这个测序深度才会被校正，可以指定具体深度，也可以指定一个比例小数，下面的解释
# --mindepth depth
#           Variants (snps and indels) will only be called if there is coverage of good pairs
#           at this depth or more; if this value is >= 1, it is an absolute depth, if it is a
#           fraction < 1, then minimum depth is computed by multiplying this value by the mean
#           coverage for the region, with a minumum value of 5 (default 0.1: min depth to call
#           is 10% of mean coverage or 5, whichever is greater).



# 第二轮校正，这轮基本就剩校正为数不多的snp, 一次跑的例子是11个
rule minimap2_2:
    input:
        genome=rules.pilon.output,
        r1 = rules.fastp.output.r1,
        r2 = rules.fastp.output.r2,
    output:
        minimap2_OUT_PATH + "/{sample}/{sample}_aln_2.bam"
    threads: 10
    shell:
        "minimap2 -ax sr {input.genome} {input.r1} {input.r2} | samtools sort -@ {threads} -O bam -o {output};"
        "samtools index -@ {threads} {output};"

rule pilon_2:
    input:
        genome = rules.pilon.output,
        bam = rules.minimap2_2.output
    output:
        pilon_OUT_PATH + "/{sample}/{sample}_pilon_2.fasta"
    params:
        prefix = "{sample}_pilon_2",
        outdir= pilon_OUT_PATH + "/{sample}",
        pilonJar= pilon_jar
    shell:
        "java -Xmx16G -jar {params.pilonJar} --genome {input.genome} --frags {input.bam} --output {params.prefix} --fix all --mindepth 10 --outdir {params.outdir} --changes --verbose"
# 第三轮校正
rule minimap2_3:
    input:
        genome=rules.pilon_2.output,
        r1 = rules.fastp.output.r1,
        r2 = rules.fastp.output.r2,
    output:
        minimap2_OUT_PATH + "/{sample}/{sample}_aln_3.bam"
    threads: 10
    shell:
        "minimap2 -ax sr {input.genome} {input.r1} {input.r2} | samtools sort -@ {threads} -O bam -o {output};"
        "samtools index -@ {threads} {output};"

rule pilon_3:
    input:
        genome = rules.pilon_2.output,
        bam = rules.minimap2_3.output
    output:
        pilon_OUT_PATH + "/{sample}/{sample}_pilon_3.fasta"
    params:
        prefix = "{sample}_pilon_3",
        outdir= pilon_OUT_PATH + "/{sample}",
        pilonJar= pilon_jar
    shell:
        "java -Xmx16G -jar {params.pilonJar} --genome {input.genome} --frags {input.bam} --output {params.prefix} --fix all --mindepth 10 --outdir {params.outdir} --changes --verbose"

# contig重命名
rule rename:
    input:
        rules.pilon_3.output,
        flye_OUT_PATH + "/{sample}/assembly_info.txt"
    output:
        unfixed_genome_PATH + "/{sample}_unfixed.fasta"
    params:
        "{sample}"
    run:
        # print(input.genome)
        df = pd.read_table(input[1])
        records = SeqIO.index(input[0], "fasta")
        plasmid_index = 1
        contig_index = 1
        records_list = []
        for i,row in df.iterrows():
            print(i)
            rec = records.get(f"{row['#seq_name']}_pilon_pilon_pilon")
            rec.description = f"{row['#seq_name']}_length_{len(rec.seq)}_cov_{row['cov.']}_circ_{row['circ.']}"
            # assembly_info.txt是有按照contig长短排序的，排第一个的是最长的染色体
            if i == 0:
                if row["circ."] == "Y":

                    rec.id = f"{params[0]}_chromosome_1"

                else:
                    print("##### chromosome is not circular!")
            else:
                if row["circ."] == "Y":
                    rec.id = f"{params[0]}_plasmid_{plasmid_index}"
                    plasmid_index += 1
                else:
                    rec.id = f"{params[0]}_contig_{contig_index}"
                    contig_index += 1


            records_list.append(rec)
        with open(output[0],"w") as f:
            SeqIO.write(records_list,f,"fasta")


# 修复成环contig两端可能存在的断裂基因
# Prodigal预测ORF，找到第一个没有ORF重叠的间隔区，把间隔区前面的序列挪到contig的尾部
rule prodigal:
    input:
        rules.rename.output
    output:
        prodigal_OUT_PATH + "/{sample}_unfixed.gff"
    shell:
        "prodigal -i {input} -f gff -o {output};"
        "sed -i '/^#/d' {output}" # 删除注释行 方便pandas读取



rule move_seq:
    input:
        rules.rename.output,
        rules.prodigal.output
    output:
        final_genome_PATH + "/{sample}.fasta"
    run:
        df_gff = pd.read_table(input[1], header=None)
        df_gff[3] = df_gff[3].astype(int)
        df_gff[4] = df_gff[4].astype(int)
        df_gff = df_gff.groupby(by=0)

        records = SeqIO.parse(input[0], "fasta")
        records_list = []
        for rec in records:

            if rec.description.endswith("circ_Y"):

                contig_id = rec.id
                group = df_gff.get_group(contig_id)
                row_index = 0
                while True:
                    if group.iloc[row_index,:][4] < group.iloc[row_index+1,:][3]:
                        start = group.iloc[row_index,:][3]
                        end = group.iloc[row_index,:][4]
                        rec.seq = rec.seq[end:] + rec.seq[0:end]
                        break
                    row_index += 1

            records_list.append(rec)
        with open(output[0],"w") as f:
            SeqIO.write(records_list,f,"fasta")



