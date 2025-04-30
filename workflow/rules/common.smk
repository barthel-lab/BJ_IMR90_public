samples = list(config["data"].keys())
types = ["ONT", "hifi", "hic", "RNA"]
longtypes = ["ONT", "hifi"]
haps = ["hap1", "hap2", "dip"]

assembly_samples = []
for sample in samples:
  datatypes = config["data"][sample]
  if (("hifi_notfiber" in datatypes) or ("hifi_fiber" in datatypes)):
    has_hifi = True
  else:
    has_hifi = False
  #if ("ONT_UL" in datatypes):
  #  has_ONT_UL = True
  #else:
  #  has_ONT_UL = False
  if ("ONT_UL" in datatypes) or ("ONT_nonfrag" in datatypes) or ("ONT_frag" in datatypes):
    has_ONT = True
  else:
    has_ONT = False
  if ("hic" in datatypes):
    has_hic = True
  else:
    has_hic = False
  if ("trio" in datatypes):
    has_trio = True
  else:
    has_trio = False
  #this requirement will be relaxed when implement the other modes
  if (has_hifi and has_ONT) and (has_hic or has_trio):
    assembly_samples.append(sample)
  #if ("illumina" in datatypes):
  #  illumina_samples.append(sample)
  #  has_illumina = True
  #else:
  #  has_illumina = False
  #if (has_hifi) and (not has_ONT) and (has_hic):
  #  assembly_samples_no_ONT.append(sample)
  #if has_hifi and has_ONT_UL and has_hic:
  #  assembly_samples.append(sample)
  #  if has_illumina:
  #    meryl_samples.append(sample)

#update this dictionary with the sex of your samples
sample_to_sex = {"BJ": "male", "IMR90": "female"}

#Functions to define the input for verkko
def get_hic_r1(wildcards):
  return config["data"][wildcards.sample]["hic"][int(wildcards.i)]["R1"]

def get_hic_r2(wildcards):
  return config["data"][wildcards.sample]["hic"][int(wildcards.i)]["R2"]

def get_hic(wildcards):
  files = []
  for i in config["data"][wildcards.sample]["hic"]:
    files.append(f"{wildcards.sample}/hic/{wildcards.sample}.hic.{i}.R{wildcards.j}.fastq.gz")
  return files

def hic_fastq(wildcards):
  files = []
  for i in config["data"][wildcards.sample]["hic"]:
    files.append(f"{wildcards.sample}/hic/{wildcards.sample}.hic.{i}.R1.fastq.gz")
    files.append(f"{wildcards.sample}/hic/{wildcards.sample}.hic.{i}.R2.fastq.gz")
  return files

def assembly_hic_r1(wildcards):
  files = []
  for i in config["data"][wildcards.sample]["hic"]:
    files.append(f"{wildcards.sample}/hic/{wildcards.sample}.hic.{i}.R1.fastq.gz")
  return files

def assembly_hic_r2(wildcards):
  files = []
  for i in config["data"][wildcards.sample]["hic"]:
    files.append(f"{wildcards.sample}/hic/{wildcards.sample}.hic.{i}.R2.fastq.gz")
  return files

def get_meryl_fofn_parental_inputs(wildcards):
  files = []
  if not wildcards.type in config["data"][wildcards.sample]["trio"]:
    return files
  for i in config["data"][wildcards.sample]["trio"][wildcards.type]:
    for R in ["R1", "R2"]:
      filename = config["data"][wildcards.sample]["trio"][wildcards.type][i][R]
      files.append(filename)
  return files

def get_meryl_fofn_proband_inputs(wildcards):
  files = []
  if not "proband" in config["data"][wildcards.sample]["trio"]:
    return files
  for i in config["data"][wildcards.sample]["trio"]["proband"]:
    for R in ["R1", "R2"]:
      filename = config["data"][wildcards.sample]["trio"]["proband"][i][R]
      files.append(filename)
  return files

def assembly_pat_hapmer_index(wildcards):
  return f"results/{wildcards.sample}/meryl/paternal_compress.k30.hapmer.meryl/merylIndex"

def assembly_mat_hapmer_index(wildcards):
  return f"results/{wildcards.sample}/meryl/maternal_compress.k30.hapmer.meryl/merylIndex"

def assembly_pat_hapmer_indir(wildcards):
  return f"results/{wildcards.sample}/meryl/paternal_compress.k30.hapmer.meryl"

def assembly_mat_hapmer_indir(wildcards):
  return f"results/{wildcards.sample}/meryl/maternal_compress.k30.hapmer.meryl"

def rna_fastq(wildcards):
  files = []
  for i in config["data"][wildcards.sample]["rna_short"]:
    files.append(config["data"][wildcards.sample]["rna_short"][int(i)]["R1"])
    files.append(config["data"][wildcards.sample]["rna_short"][int(i)]["R2"])
  return files

def get_all_rna_r1(wildcards):
  files = []
  for i in config["data"][wildcards.sample]["rna_short"]:
    files.append(f"figures/external_validation/fastp/{wildcards.sample}/RNA/{wildcards.sample}.{i}.rna.paired_r1.fastp.fastq.gz")
  return files

def get_all_rna_r2(wildcards):
  files = []
  for i in config["data"][wildcards.sample]["rna_short"]:
    files.append(f"figures/external_validation/fastp/{wildcards.sample}/RNA/{wildcards.sample}.{i}.rna.paired_r2.fastp.fastq.gz")
  return files

def get_all_rna_unpaired(wildcards):
  files = []
  for i in config["data"][wildcards.sample]["rna_short"]:
    files.append(f"figures/external_validation/fastp/{wildcards.sample}/RNA/{wildcards.sample}.{i}.rna.unpaired_r1.fastp.fastq.gz")
  for i in config["data"][wildcards.sample]["rna_short"]:
    files.append(f"figures/external_validation/fastp/{wildcards.sample}/RNA/{wildcards.sample}.{i}.rna.unpaired_r2.fastp.fastq.gz")
  return files

def get_ubam_ONT(wildcards):
  return config["data"][wildcards.sample][wildcards.type][int(wildcards.i)]

def assembly_ONT_fastq(wildcards):
  files = []
  if "ONT_UL" in config["data"][wildcards.sample]:
    for i in config["data"][wildcards.sample]["ONT_UL"]:
      files.append(f"{wildcards.sample}/ONT/{wildcards.sample}.ONT_UL.{i}.fastq")
  return files

def assembly_ONT(wildcards):
  files = []
  files.append(f"{wildcards.sample}/ONT/{wildcards.sample}.all.ONT.fastq")
  return files

def get_ubam_hifi(wildcards):
  return config["data"][wildcards.sample][wildcards.type][int(wildcards.i)]

def hifi_fasta(wildcards):
  files = []
  if "hifi_notfiber" in config["data"][wildcards.sample]:
    for i in config["data"][wildcards.sample]["hifi_notfiber"]:
      files.append(f"{wildcards.sample}/hifi/{wildcards.sample}.hifi_notfiber.{i}.fasta")
  if "hifi_fiber" in config["data"][wildcards.sample]:
    for i in config["data"][wildcards.sample]["hifi_fiber"]:
      files.append(f"{wildcards.sample}/hifi/{wildcards.sample}.hifi_fiber.{i}.fasta")
  return files

def assembly_hifi_error_corrected_ont(wildcards):
  #We include all the hifi data plus the herro corrected ONT data as "hifi" input to verkko
  files = []
  if "hifi_notfiber" in config["data"][wildcards.sample]:
    for i in config["data"][wildcards.sample]["hifi_notfiber"]:
      files.append(f"{wildcards.sample}/hifi/{wildcards.sample}.hifi_notfiber.{i}.fasta")
  if "hifi_fiber" in config["data"][wildcards.sample]:
    for i in config["data"][wildcards.sample]["hifi_fiber"]:
      files.append(f"{wildcards.sample}/hifi/{wildcards.sample}.hifi_fiber.{i}.fasta")
  files.append(f"{wildcards.sample}/ONT/{wildcards.sample}.all.ONT.corrected.fasta")
  return files

def assembly_hifi(wildcards):
  files = []
  if "hifi_notfiber" in config["data"][wildcards.sample]:
    for i in config["data"][wildcards.sample]["hifi_notfiber"]:
      files.append(f"{wildcards.sample}/hifi/{wildcards.sample}.hifi_notfiber.{i}.fasta")
  if "hifi_fiber" in config["data"][wildcards.sample]:
    for i in config["data"][wildcards.sample]["hifi_fiber"]:
      files.append(f"{wildcards.sample}/hifi/{wildcards.sample}.hifi_fiber.{i}.fasta")
  return files

def assembly_hifi_ubam(wildcards):
  files = []
  if "hifi_notfiber" in config["data"][wildcards.sample]:
    for i in config["data"][wildcards.sample]["hifi_notfiber"]:
      files.append(config["data"][wildcards.sample]["hifi_notfiber"][i])
  if "hifi_fiber" in config["data"][wildcards.sample]:
    for i in config["data"][wildcards.sample]["hifi_fiber"]:
      files.append(config["data"][wildcards.sample]["hifi_fiber"][i])
  return files

def get_illumina(wildcards):
  return config["data"][wildcards.sample]["illumina"][int(wildcards.i)][f"R{wildcards.j}"]

def get_illumina_r1(wildcards):
  return config["data"][wildcards.sample]["illumina"][int(wildcards.i)]["R1"]

def get_illumina_r2(wildcards):
  return config["data"][wildcards.sample]["illumina"][int(wildcards.i)]["R2"]

def get_meryl_inputs(wildcards):
  inputs = []
  for i in config["data"][wildcards.sample]["illumina"]:
    fileinput = f"{wildcards.sample}/illumina/meryl/{wildcards.sample}.{i}.R1.meryl/merylIndex"
    inputs.append(fileinput)
    fileinput = f"{wildcards.sample}/illumina/meryl/{wildcards.sample}.{i}.R2.meryl/merylIndex"
    inputs.append(fileinput)
  return inputs

def get_meryl_indirs(wildcards):
  indirs = []
  for i in config["data"][wildcards.sample]["illumina"]:
    filename = f"{wildcards.sample}/illumina/meryl/{wildcards.sample}.{i}.R1.meryl"
    indirs.append(filename)
    filename = f"{wildcards.sample}/illumina/meryl/{wildcards.sample}.{i}.R2.meryl"
    indirs.append(filename)
  return indirs

merqury_summary_inputs = []
for sample in samples:
  for assembly in ["hg38", "CHM13", "precuration"]:
    merqury_summary_inputs.append(f"figures/precuration_qc/merqury/{sample}/{sample}.{assembly}.merqury.qv")
    merqury_summary_inputs.append(f"figures/precuration_qc/merqury/{sample}/{sample}.{assembly}.merqury.completeness.stats")
  merqury_summary_inputs.append(f"figures/postcuration_qc/merqury/{sample}/{sample}.postcuration.merqury.qv")
  merqury_summary_inputs.append(f"figures/postcuration_qc/merqury/{sample}/{sample}.postcuration.merqury.completeness.stats")

hisat2_parse_inputs = []
for sample in samples:
  hisat2_parse_inputs.append(f"figures/precuration_qc/hisat/{sample}/hg38/{sample}.rna.hg38.AlignMetrics.txt")
  hisat2_parse_inputs.append(f"figures/precuration_qc/hisat/{sample}/CHM13/{sample}.rna.CHM13.AlignMetrics.txt")
  hisat2_parse_inputs.append(f"figures/postcuration_qc/hisat/{sample}/dip/{sample}.rna.{sample}.dip.postcuration.AlignMetrics.txt")

chipseq_parse_inputs = []
for sample in samples:
  chipseq_parse_inputs.append(f"figures/precuration_qc/bwa/{sample}/hg38/{sample}.all.chipseq.hg38.AlignMetrics.txt")
  chipseq_parse_inputs.append(f"figures/precuration_qc/bwa/{sample}/CHM13/{sample}.all.chipseq.CHM13.AlignMetrics.txt")
  chipseq_parse_inputs.append(f"figures/postcuration_qc/bwa/{sample}/dip/{sample}.all.chipseq.{sample}.dip.postcuration.AlignMetrics.txt")

def get_rna_r1(wildcards):
  return config["data"][wildcards.sample]["rna_short"][int(wildcards.i)]["R1"]

def get_rna_r2(wildcards):
  return config["data"][wildcards.sample]["rna_short"][int(wildcards.i)]["R2"]

def get_chipseq(wildcards):
  return config["data"][wildcards.sample]["chipseq"][int(wildcards.i)]

def get_chipseq_fastp(wildcards):
  return f"figures/external_validation/fastp/{wildcards.sample}/chipseq/{wildcards.sample}.{wildcards.i}.chipseq.fastp.fastq.gz"

def get_chroms(wildcards):
  assert wildcards.sample in sample_to_sex, "Please update the sample_to_sex dictionary in common.smk with the sex of your samples."
  sex = sample_to_sex[wildcards.sample]
  autosomes_hap1 = [f"chr{i}_hap1" for i in range(1, 23)]
  autosomes_hap2 = [f"chr{i}_hap2" for i in range(1, 23)]
  if sex == "female":
    chroms_hap1 = autosomes_hap1 + ["chrX_hap1"]
  if sex == "male":
    chroms_hap1 = autosomes_hap1 + ["chrY_hap1"]
  chroms_hap2 = autosomes_hap2 + ["chrX_hap2"]
  chroms = chroms_hap1 + chroms_hap2
  if wildcards.hap == "hap1":
    return chroms_hap1
  elif wildcards.hap == "hap2":
    return chroms_hap2
  elif wildcards.hap == "dip":
    return chroms

def get_CHM13_hap(wildcards):
  sex = sample_to_sex[wildcards.sample]
  if sex == "male" and int(wildcards.i) == 1:
    return "/tgen_labs/barthel/references/CHM13v2/chm13v2.0.autosomes_plus_chrY.fasta"
  else:
    return "/tgen_labs/barthel/references/CHM13v2/chm13v2.0.autosomes_plus_chrX.fasta"
