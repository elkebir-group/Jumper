configfile: "config.yaml"
#configfile: "test.yaml"

def parse_inputFile(file):
    with open(file) as ifile:
        entries = [entry for entry in ifile.read().split('\n') if len(entry) > 0]
    return entries


rule all:
    input:
        expand("jumper/nedges_{nedges}/{sample}_decomposition.out", sample=parse_inputFile(config["sampleList"]), nedges=config["nedges"]),
        expand("jumper/nedges_{nedges}/{sample}_matching.out", sample=parse_inputFile(config["sampleList"]), nedges=config["nedges"]),
        expand("jumper/nedges_{nedges}/{sample}_phasing.out", sample=parse_inputFile(config["sampleList"]), nedges=config["nedges"]),
        expand("jumper/nedges_{nedges}/{sample}_breakpoints.out", sample=parse_inputFile(config["sampleList"]), nedges=config["nedges"]),
        expand("jumper/nedges_{nedges}/{sample}_edges.out", sample=parse_inputFile(config["sampleList"]), nedges=config["nedges"]),
        expand("jumper/nedges_{nedges}/{sample}_graph.dot", sample=parse_inputFile(config["sampleList"]), nedges=config["nedges"]),
        expand("jumper/nedges_{nedges}/{sample}_sj_reads.tsv", sample=parse_inputFile(config["sampleList"]), nedges=config["nedges"]),
        expand("jumper/nedges_{nedges}/{sample}_canonical_transcripts.out", sample=parse_inputFile(config["sampleList"]), nedges=config["nedges"]),


rule jumper:
    input:
        bam="bam/{sample}.Aligned.sortedByCoord.out.bam",
    output:
        decomposition="jumper/nedges_{nedges}/{sample}_decomposition.out",
        matching="jumper/nedges_{nedges}/{sample}_matching.out",
        phasing="jumper/nedges_{nedges}/{sample}_phasing.out",
        breakpoints="jumper/nedges_{nedges}/{sample}_breakpoints.out",
        edges="jumper/nedges_{nedges}/{sample}_edges.out",
        sj_reads="jumper/nedges_{nedges}/{sample}_sj_reads.tsv",
        graph="jumper/nedges_{nedges}/{sample}_graph.dot",
        transcripts="jumper/nedges_{nedges}/{sample}_canonical_transcripts.out"
    benchmark: "jumper/nedges_{nedges}/{sample}_benchmark.log"
    log:
        std="jumper/nedges_{nedges}/{sample}.log",
        err="jumper/nedges_{nedges}/{sample}.err.log"
    params:
        refFile=lambda wildcards: config["reference"][wildcards.sample],
        npaths=config["ntranscripts"],
        samplingFrequency=config["samplingFrequency"],
        timelimit=config["timelimit"],
        nthreads=config["nthreads"],
        is_paired=lambda wildcards: config["is_paired"][wildcards.sample]
    shell:
        #"jumper --fasta {params.refFile} --bam {input.bam} "
        "python ../jumper/jumper_main.py --fasta {params.refFile} --bam {input.bam} "
        "--outputBreakpoints {output.breakpoints} "
        "--outputEdges {output.edges} "
        "--inputCSV jumper_graph/{wildcards.sample}_sj_reads.tsv "
        "--outputCSV {output.sj_reads} "
        "--outputDOT {output.graph} "
        "--outputDecomposition {output.decomposition} "
        "--outputMatching {output.matching} "
        "--outputPhasing {output.phasing} "
        "--outputTranscripts {output.transcripts} "
        "--greedy True --paired {params.is_paired} --width 0 "
        "--samplingFrequency {params.samplingFrequency} "
        "--sj_threshold 0 -n {wildcards.nedges} -k {params.npaths} --threads {params.nthreads} "
        "--timelimit {params.timelimit} "
        "> {log.std} 2> {log.err}"
