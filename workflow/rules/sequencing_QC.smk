#genomescope plots
rule kmc_ONT:
  input:
    "{sample}/ONT/{sample}.all.ONT.fastq"
  output:
    pre=temp("figures/sequencing_qc/genomescope/ONT/{sample}/KMC/{sample}.ONT.kmc_pre"),
    suf=temp("figures/sequencing_qc/genomescope/ONT/{sample}/KMC/{sample}.ONT.kmc_suf"),
  params:
    kmcdb="figures/sequencing_qc/genomescope/ONT/{sample}/KMC/{sample}.ONT",
    tmp="figures/sequencing_qc/genomescope/ONT/{sample}/KMC/tmp_{sample}.ONT"
  threads: 42
  resources:
    mem_mb=150000
  shell:
    """
    mkdir -p {params.tmp}
    kmc -k31 -t{threads} -m100 -sm -ci1 -cs100000000 -fq {input} {params.kmcdb} {params.tmp}
    """

rule kmc_hifi:
  input:
    "{sample}/hifi/{sample}.all.hifi.fasta"
  output:
    pre=temp("figures/sequencing_qc/genomescope/hifi/{sample}/KMC/{sample}.hifi.kmc_pre"),
    suf=temp("figures/sequencing_qc/genomescope/hifi/{sample}/KMC/{sample}.hifi.kmc_suf"),
  params:
    kmcdb="figures/sequencing_qc/genomescope/hifi/{sample}/KMC/{sample}.hifi",
    tmp="figures/sequencing_qc/genomescope/hifi/{sample}/KMC/tmp_{sample}.hifi"
  threads: 42
  resources:
    mem_mb=150000
  shell:
    """
    mkdir -p {params.tmp}
    kmc -k31 -t{threads} -m100 -sm -ci1 -cs100000000 -fm {input} {params.kmcdb} {params.tmp}
    """

rule nonzero_hist:
  input:
    pre="figures/sequencing_qc/genomescope/{type}/{sample}/KMC/{sample}.{type}.kmc_pre",
    suf="figures/sequencing_qc/genomescope/{type}/{sample}/KMC/{sample}.{type}.kmc_suf",
  output:
    "figures/sequencing_qc/genomescope/{type}/{sample}/KMC/{sample}.{type}.nonzero.hist"
  params:
    kmcdb="figures/sequencing_qc/genomescope/{type}/{sample}/KMC/{sample}.{type}",
    tmphist="figures/sequencing_qc/genomescope/{type}/{sample}/KMC/{sample}.{type}.zero.hist"
  threads: 1
  resources:
    mem_mb=150000
  shell:
    """
    kmc_tools transform {params.kmcdb} histogram {params.tmphist} -cx100000000
    grep -vP "\\t0" {params.tmphist} > {output}
    rm {params.tmphist}
    """

rule genomescope2:
  input:
    "figures/sequencing_qc/genomescope/{type}/{sample}/KMC/{sample}.{type}.nonzero.hist"
  output:
    "figures/sequencing_qc/genomescope/{type}/{sample}/linear_plot.png"
  params:
    outdir="figures/sequencing_qc/genomescope/{type}/{sample}"
  threads: 1
  resources:
    mem_mb=50000
  conda: "R_v4.4.1"
  shell:
    "genomescope.R -i {input} -o {params.outdir} -k31 --num_rounds 1"

#nanoplot plots
rule NanoPlot_combined_ONT:
  input:
    "{sample}/ONT/{sample}.all.ONT.fastq"
  output:
    "figures/sequencing_qc/NanoPlot/ONT/{sample}/{sample}.combined.ONT.NanoPlot-report.html",
    "figures/sequencing_qc/NanoPlot/ONT/{sample}/{sample}.combined.ONT.NanoStats.txt"
  params:
    outdir="figures/sequencing_qc/NanoPlot/ONT/{sample}/",
    prefix="{sample}.combined.ONT."
  threads: 42
  resources:
    mem_mb=150000
  shell:
    "NanoPlot --fastq {input} -t {threads} -o {params.outdir} -p {params.prefix} -f pdf --N50"

rule NanoPlot_combined_hifi:
  input:
    assembly_hifi_ubam
  output:
    "figures/sequencing_qc/NanoPlot/hifi/{sample}/{sample}.combined.hifi.NanoPlot-report.html",
    "figures/sequencing_qc/NanoPlot/hifi/{sample}/{sample}.combined.hifi.NanoStats.txt"
  params:
    outdir="figures/sequencing_qc/NanoPlot/hifi/{sample}/",
    prefix="{sample}.combined.hifi."
  threads: 42
  resources:
    mem_mb=150000
  shell:
    "NanoPlot --ubam {input} -t {threads} -o {params.outdir} -p {params.prefix} -f pdf --N50"

rule NanoPlot_combined_hic:
  input:
    hic_fastq
  output:
    "figures/sequencing_qc/NanoPlot/hic/{sample}/{sample}.combined.hic.NanoPlot-report.html",
    "figures/sequencing_qc/NanoPlot/hic/{sample}/{sample}.combined.hic.NanoStats.txt"
  params:
    outdir="figures/sequencing_qc/NanoPlot/hic/{sample}/",
    prefix="{sample}.combined.hic."
  threads: 42
  resources:
    mem_mb=150000
  shell:
    "NanoPlot --fastq {input} -t {threads} -o {params.outdir} -p {params.prefix} -f pdf --N50"

rule NanoPlot_combined_rna:
  input:
    rna_fastq
  output:
    "figures/sequencing_qc/NanoPlot/RNA/{sample}/{sample}.combined.RNA.NanoPlot-report.html",
    "figures/sequencing_qc/NanoPlot/RNA/{sample}/{sample}.combined.RNA.NanoStats.txt"
  params:
    outdir="figures/sequencing_qc/NanoPlot/RNA/{sample}/",
    prefix="{sample}.combined.RNA."
  threads: 42
  resources:
    mem_mb=150000
  shell:
    "NanoPlot --fastq {input} -t {threads} -o {params.outdir} -p {params.prefix} -f pdf --N50"

#podplot plots
rule make_fai_fofn_combined_ONT:
  input:
    expand("{sample}/ONT/{sample}.all.ONT.fastq.fai", sample=samples)
  output:
    "figures/sequencing_qc/podplot/ONT/all.combined.ONT.fastq.fai.fofn"
  params:
    samples=expand("{sample} ONT", sample=samples)
  threads: 1
  resources:
    mem_mb=50000
  script:
    "../scripts/make_fai_fofn.py"

rule make_fai_fofn_combined_hifi:
  input:
    expand("{sample}/hifi/{sample}.all.hifi.fasta.fai", sample=samples)
  output:
    "figures/sequencing_qc/podplot/hifi/all.combined.hifi.fasta.fai.fofn"
  params:
    samples=expand("{sample} HiFi", sample=samples)
  threads: 1
  resources:
    mem_mb=50000
  script:
    "../scripts/make_fai_fofn.py"

rule podplot_combined_ONT:
  input:
    "figures/sequencing_qc/podplot/ONT/all.combined.ONT.fastq.fai.fofn"
  output:
    "figures/sequencing_qc/podplot/ONT/all.combined.ONT.podplot.png"
  params:
    size=3117292070,
    script_path=Path(workflow.basedir)/"scripts/podplot.py"
  threads: 1
  resources:
    mem_mb=50000
  shell:
    "{params.script_path} -g {params.size} {input} {output}"

rule podplot_combined_hifi:
  input:
    "figures/sequencing_qc/podplot/hifi/all.combined.hifi.fasta.fai.fofn"
  output:
    "figures/sequencing_qc/podplot/hifi/all.combined.hifi.podplot.png"
  params:
    size=3117292070,
    script_path=Path(workflow.basedir)/"scripts/podplot.py"
  threads: 1
  resources:
    mem_mb=50000
  shell:
    "{params.script_path} -g {params.size} {input} {output}"

#hic_qc plots
#rule hic_qc:
#  input:
#    ""
#  output:
#    "figures/sequencing_qc/hic_qc/hic/{sample}/{sample}.combined.hic.range.png"
#  threads: 
#  resources:
#    mem_mb=
#  shell:
#    "python hic_qc.py -b input.bam -n num_reads_to_use -r"

#summary table
rule sequencing_qc_summary_table:
  input:
    expand("figures/sequencing_qc/NanoPlot/{kind}/{sample}/{sample}.combined.{kind}.NanoStats.txt", sample=samples, kind=types)
  output:
    "figures/sequencing_qc/summary/sequencing_qc_summary_table.tsv"
  threads: 1
  resources:
    mem_mb=50000
  params:
    samples=samples,
    kinds=types
  script:
    "../scripts/sequencing_qc_summary_table.py"
