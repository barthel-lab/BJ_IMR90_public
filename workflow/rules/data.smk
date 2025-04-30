rule trim_hic:
  input:
    r1=get_hic_r1,
    r2=get_hic_r2
  output:
    r1="{sample}/hic/{sample}.hic.{i,[0-9]+}.R1.fastq.gz",
    r2="{sample}/hic/{sample}.hic.{i,[0-9]+}.R2.fastq.gz"
  conda: "cutadaptenv"
  threads: 42
  resources:
    mem_mb=250000
  shell:
    "cutadapt --cores {threads} -u 5 -U 5 -o {output.r1} -p {output.r2} {input.r1} {input.r2}"

rule merge_hic:
  input:
    get_hic
  output:
    "{sample}/hic/{sample}.hic.all.R{j}.fastq.gz"
  threads: 1
  resources:
    mem_mb=50000
  run:
    if len(input) != 1:
      shell("cat {input} > {output}")
    else:
      shell("ln -sr {input} {output}")

rule copy_tags:
  output:
    "tags.txt"
  threads: 1
  resources:
    mem_mb=50000
  params:
    tags=Path(workflow.basedir)/"../tags.txt"
  shell:
    "cp {params.tags} {output}"

rule samtools_fastq_ONT:
  input:
    tags="tags.txt",
    reads=get_ubam_ONT
  output:
    temp("{sample}/ONT/{sample}.{type}.{i, [0-9]+}.fastq")
  threads: 42
  resources:
    mem_mb=250000
  shell:
    "samtools view -b -h -D dx:tags.txt {input.reads} | samtools fastq -@ {threads} - > {output}"

rule combine_fastq:
  input:
    assembly_ONT_fastq
  output:
    "{sample}/ONT/{sample}.all.ONT.fastq"
  threads: 1
  resources:
    mem_mb=150000
  run:
    if len(input) > 1:
      shell("cat {input} > {output}")
    else:
      shell("ln -sr {input} {output}")

rule index_fastq:
  input:
    "{sample}/ONT/{sample}.all.ONT.fastq"
  output:
    "{sample}/ONT/{sample}.all.ONT.fastq.fai"
  threads: 1
  resources:
    mem_mb=150000
  shell:
    "samtools fqidx {input}"

rule herro_all_cpu:
  input:
    fq="{sample}/ONT/{sample}.all.ONT.fastq",
    fai="{sample}/ONT/{sample}.all.ONT.fastq.fai"
  output:
    "{sample}/ONT/{sample}.all.ONT.overlaps.paf"
  threads: 42
  resources:
    mem_mb=384000,
    #slurm_partition="bigmem",
    runtime=1440
  shell:
    "dorado correct {input.fq} --to-paf > {output}"

rule herro_all_gpu:
  input:
    fq="{sample}/ONT/{sample}.all.ONT.fastq",
    paf="{sample}/ONT/{sample}.all.ONT.overlaps.paf"
  output:
    "{sample}/ONT/{sample}.all.ONT.corrected.fasta"
  threads: 32
  resources:
    mem_mb=384000,
    gpu=2,
    gpu_model="a100",
    slurm_partition="gpu-a100",
    runtime=4320,
    cpus_per_gpu=16
  shell:
    "dorado correct {input.fq} --from-paf {input.paf} > {output}"

#rule herro:
#  input:
#    fq="{sample}/ONT/{sample}.all.ONT.fastq",
#    fai="{sample}/ONT/{sample}.all.ONT.fastq.fai"
#  output:
#    "{sample}/ONT/{sample}.all.ONT.corrected.fasta"
#  threads: 32
#  resources:
#    mem_mb=384000,
#    slurm_partition="gpu-a100",
#    runtime=5760,
#    nodes=1,
#    tasks=1,
#    cpus_per_task=1,
#    slurm_extra="'--gres=gpu:A100:2 --cpus-per-gpu 16'"
#  shell:
#    "dorado correct {input.fq} > {output}"

rule samtools_fasta_hifi:
  input:
    get_ubam_hifi
  output:
    "{sample}/hifi/{sample}.{type}.{i, [0-9]+}.fasta"
  threads: 42
  resources:
    mem_mb=250000
  shell:
    "samtools fasta -@ {threads} {input} > {output}"

rule combine_fasta:
  input:
    hifi_fasta
  output:
    "{sample}/hifi/{sample}.all.hifi.fasta"
  threads: 1
  resources:
    mem_mb=50000
  run:
    if len(input) > 1:
      shell("cat {input} > {output}")
    else:
      shell("ln -sr {input} {output}")

rule index_fasta:
  input:
    "{sample}/hifi/{sample}.all.hifi.fasta"
  output:
    "{sample}/hifi/{sample}.all.hifi.fasta.fai"
  threads: 1
  resources:
    mem_mb=50000
  shell:
    "samtools faidx {input}"
