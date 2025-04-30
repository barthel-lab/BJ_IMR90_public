rule meryl_fofn_parental:
  input:
    get_meryl_fofn_parental_inputs
  output:
    "results/{sample}/meryl/{type}/{type}.fofn"
  threads: 1
  resources:
    mem_mb=50000
  shell:
    "ls {input} > {output}"

rule meryl_fofn_proband:
  input:
    get_meryl_fofn_proband_inputs
  output:
    "results/{sample}/meryl/proband.fofn"
  threads: 1
  resources:
    mem_mb=50000
  shell:
    "ls {input} > {output}"

rule meryl_submit_parental:
  input:
    "results/{sample}/meryl/{type}/{type}.fofn"
  output:
    "results/{sample}/meryl/{type}/{type}_compress.k30.meryl/merylIndex"
  params:
    outdir="results/{sample}/meryl/{type}",
    inputmeryl="{type}.fofn",
    outprefix="{type}_compress",
    meryl_db="{type}/{type}_compress.k30.meryl"
  shell:
    """
    cd {params.outdir}
    _submit_build.sh -c 30 {params.inputmeryl} {params.outprefix}
    cd ..
    sleep 300
    ln -s {params.meryl_db}
    """

rule meryl_submit_proband:
  input:
    "results/{sample}/meryl/proband.fofn"
  output:
    "results/{sample}/meryl/proband_compress.k30.meryl/merylIndex"
  params:
    outdir="results/{sample}/meryl",
  shell:
    """
    cd {params.outdir}
    _submit_build.sh -c 30 proband.fofn proband_compress
    sleep 300
    """

rule hapmers_proband:
  input:
    "results/{sample}/meryl/maternal/maternal_compress.k30.meryl/merylIndex",
    "results/{sample}/meryl/paternal/paternal_compress.k30.meryl/merylIndex",
    "results/{sample}/meryl/proband_compress.k30.meryl/merylIndex"
  output:
    "results/{sample}/meryl/paternal_compress.k30.hapmer.meryl/merylIndex",
    "results/{sample}/meryl/maternal_compress.k30.hapmer.meryl/merylIndex"
  params:
    mathapdir="results/{sample}/meryl/maternal_compress.k30.hapmer.meryl",
    pathapdir="results/{sample}/meryl/paternal_compress.k30.hapmer.meryl",
    outdir="results/{sample}/meryl",
    mat="maternal_compress.k30.meryl",
    pat="paternal_compress.k30.meryl",
    proband="proband_compress.k30.meryl"
  shell:
    """
    rm -r {params.mathapdir}
    rm -r {params.pathapdir}
    cd {params.outdir}
    hapmers.sh {params.mat} {params.pat} {params.proband}
    """

rule hapmers_noproband:
  input:
    "results/{sample}/meryl/maternal/maternal_compress.k30.meryl/merylIndex",
    "results/{sample}/meryl/paternal/paternal_compress.k30.meryl/merylIndex"
  output:
    "results/{sample}/meryl/paternal_compress.k30.only.meryl/merylIndex",
    "results/{sample}/meryl/maternal_compress.k30.only.meryl/merylIndex"
  params:
    mathapdir="results/{sample}/meryl/maternal_compress.k30.only.meryl",
    pathapdir="results/{sample}/meryl/paternal_compress.k30.only.meryl",
    outdir="results/{sample}/meryl",
    mat="maternal_compress.k30.meryl",
    pat="paternal_compress.k30.meryl"
  shell:
    """
    rm -r {params.mathapdir}
    rm -r {params.pathapdir}
    cd {params.outdir}
    hapmers.sh {params.mat} {params.pat}
    """

rule verkko_hifi_error_corrected_ont_hic_phased:
  input:
    r1=assembly_hic_r1,
    r2=assembly_hic_r2,
    ont=assembly_ONT,
    hifi=assembly_hifi_error_corrected_ont
  output:
    protected("{sample}/assembly/verkko_error_corrected_ont/assembly.fasta"),
    protected("{sample}/assembly/verkko_error_corrected_ont/assembly.haplotype1.fasta"),
    protected("{sample}/assembly/verkko_error_corrected_ont/assembly.haplotype2.fasta"),
    protected("{sample}/assembly/verkko_error_corrected_ont/assembly.unassigned.fasta"),
    protected("{sample}/assembly/verkko_error_corrected_ont/assembly.disconnected.fasta")
  params:
    outdir="{sample}/assembly/verkko_error_corrected_ont",
  conda: "verkko"
  threads: 1
  resources:
    mem_mb=50000,
    runtime=4320
  shell:
    "verkko -d {params.outdir} --hifi {input.hifi} --nano {input.ont} --hic1 {input.r1} --hic2 {input.r2} --grid"

rule verkko_hifi_error_corrected_ont_trio_phased:
  input:
    pat=assembly_pat_hapmer_index,
    mat=assembly_mat_hapmer_index,
    ont=assembly_ONT,
    hifi=assembly_hifi_error_corrected_ont
  output:
    protected("{sample}/assembly/verkko_error_corrected_ont/assembly.fasta"),
    protected("{sample}/assembly/verkko_error_corrected_ont/assembly.haplotype1.fasta"),
    protected("{sample}/assembly/verkko_error_corrected_ont/assembly.haplotype2.fasta"),
    protected("{sample}/assembly/verkko_error_corrected_ont/assembly.unassigned.fasta"),
    protected("{sample}/assembly/verkko_error_corrected_ont/assembly.disconnected.fasta")
  params:
    outdir="{sample}/assembly/verkko_error_corrected_ont",
    pat_hapmer_indir=assembly_pat_hapmer_indir,
    mat_hapmer_indir=assembly_mat_hapmer_indir
  conda: "verkko"
  threads: 1
  resources:
    mem_mb=50000,
    runtime=4320
  shell:
    "verkko -d {params.outdir} --hifi {input.hifi} --nano {input.ont} --hap-kmers {params.pat_hapmer_indir} {params.mat_hapmer_indir} trio --grid"

rule verkko_hifi_ont_hic_phased:
  input:
    r1=assembly_hic_r1,
    r2=assembly_hic_r2,
    ont=assembly_ONT,
    hifi=assembly_hifi
  output:
    protected("{sample}/assembly/verkko/assembly.fasta"),
    protected("{sample}/assembly/verkko/assembly.haplotype1.fasta"),
    protected("{sample}/assembly/verkko/assembly.haplotype2.fasta"),
    protected("{sample}/assembly/verkko/assembly.unassigned.fasta"),
    protected("{sample}/assembly/verkko/assembly.disconnected.fasta")
  params:
    outdir="{sample}/assembly/verkko",
  conda: "verkko"
  threads: 1
  resources:
    mem_mb=50000,
    runtime=4320
  shell:
    "verkko -d {params.outdir} --hifi {input.hifi} --nano {input.ont} --hic1 {input.r1} --hic2 {input.r2} --grid"

rule verkko_hifi_ont_trio_phased:
  input:
    pat=assembly_pat_hapmer_index,
    mat=assembly_mat_hapmer_index,
    ont=assembly_ONT,
    hifi=assembly_hifi
  output:
    protected("{sample}/assembly/verkko/assembly.fasta"),
    protected("{sample}/assembly/verkko/assembly.haplotype1.fasta"),
    protected("{sample}/assembly/verkko/assembly.haplotype2.fasta"),
    protected("{sample}/assembly/verkko/assembly.unassigned.fasta"),
    protected("{sample}/assembly/verkko/assembly.disconnected.fasta")
  params:
    outdir="{sample}/assembly/verkko",
    pat_hapmer_indir=assembly_pat_hapmer_indir,
    mat_hapmer_indir=assembly_mat_hapmer_indir
  conda: "verkko"
  threads: 1
  resources:
    mem_mb=50000,
    runtime=4320
  shell:
    "verkko -d {params.outdir} --hifi {input.hifi} --nano {input.ont} --hap-kmers {params.pat_hapmer_indir} {params.mat_hapmer_indir} trio --grid"
