configfile: "config.yaml"

rule all:
    input:
        expand("jumper/{sense}_{seed}/sample_{rep}_decomposition.out", sense=config["senses"], seed=config["seeds"], rep=config["rep"]),
        expand("jumper/{sense}_{seed}/sample_{rep}_matching.out", sense=config["senses"], seed=config["seeds"], rep=config["rep"]),
        expand("jumper/{sense}_{seed}/sample_{rep}_phasing.out", sense=config["senses"], seed=config["seeds"], rep=config["rep"]),
        expand("jumper/{sense}_{seed}/sample_{rep}_breakpoints.out", sense=config["senses"], seed=config["seeds"], rep=config["rep"]),
        expand("jumper/{sense}_{seed}/sample_{rep}_edges.out", sense=config["senses"], seed=config["seeds"], rep=config["rep"]),
        expand("jumper/{sense}_{seed}/sample_{rep}_graph.dot", sense=config["senses"], seed=config["seeds"], rep=config["rep"]),
        expand("jumper/{sense}_{seed}/sample_{rep}_sj_reads.tsv", sense=config["senses"], seed=config["seeds"], rep=config["rep"]),
        expand("jumper/{sense}_{seed}/sample_{rep}_canonical_transcripts.out", sense=config["senses"], seed=config["seeds"], rep=config["rep"]),


rule jumper:
    input:
        bam="bam/star/{sense}_{seed}/sample_{rep}.Aligned.sortedByCoord.out.bam",
    output:
        decomposition="jumper/{sense}_{seed}/sample_{rep}_decomposition.out",
        matching="jumper/{sense}_{seed}/sample_{rep}_matching.out",
        phasing="jumper/{sense}_{seed}/sample_{rep}_phasing.out",
        breakpoints="jumper/{sense}_{seed}/sample_{rep}_breakpoints.out",
        edges="jumper/{sense}_{seed}/sample_{rep}_edges.out",
        sj_reads="jumper/{sense}_{seed}/sample_{rep}_sj_reads.tsv",
        graph="jumper/{sense}_{seed}/sample_{rep}_graph.dot",
        transcripts="jumper/{sense}_{seed}/sample_{rep}_canonical_transcripts.out"
    log:
        std="jumper/{sense}_{seed}/sample_{rep}.log",
        err="jumper/{sense}_{seed}/sample_{rep}.err.log"
    params:
        refFile=config["ref"],
        npaths=config["ntranscripts"],
        samplingFrequency=config["samplingFrequency"],
        timelimit=config["timelimit"],
        nthreads=config["nthreads"],
        resolution=config["resolution"],
        sj_threshold=config["sj_threshold"]
    shell:
        "python jumper.py --fasta {params.refFile} --bam {input.bam} "
        "--outputBreakpoints {output.breakpoints} "
        "--outputEdges {output.edges} "
        "--outputCSV {output.sj_reads} "
        "--outputDOT {output.graph} "
        "--outputDecomposition {output.decomposition} "
        "--outputPhasing {output.phasing} "
        "--outputMatching {output.matching} "
        "--outputTranscripts {output.transcripts} "
        "--greedy True "
        "--paired True --sj_threshold {params.sj_threshold} -k {params.npaths} --threads {params.nthreads} "
        "--timelimit {params.timelimit} "
        "--samplingFrequency {params.samplingFrequency} "
        "> {log.std} 2> {log.err}"
