configfile: "config.yaml"

rule all:
    input:
        expand("abundance/stringtie/{sense}_{seed}/sample_{rep}.assem.gtf", sense=config["senses"], seed=config["seeds"], rep=config["rep"]),
        expand("abundance/scallop/{sense}_{seed}/sample_{rep}.gtf", sense=config["senses"], seed=config["seeds"], rep=config["rep"]),


rule stringtie_abund:
    input:
        bam="bam/star_XS/{sense}_{seed}/sample_{rep}.Aligned.sortedByCoord.out.bam",
    threads: 1,
    benchmark: "abundance/benchmark/stringtie_{sense}_{seed}_sample_{rep}.benchmark"
    output:
        assembled="abundance/stringtie/{sense}_{seed}/sample_{rep}.assem.gtf",
        abundance="abundance/stringtie/{sense}_{seed}/sample_{rep}.abund.tab",
    shell:
        "stringtie -p {threads} "
        "-o {output.assembled} "
        "-A {output.abundance} "
        "{input.bam} "

rule scallop_abund:
    input:
        bam="bam/star/{sense}_{seed}/sample_{rep}.Aligned.sortedByCoord.out.bam",
    log: "log/scallop_{sense}_{seed}_sample_{rep}.log"
    benchmark: "abundance/benchmark/scallop_{sense}_{seed}_sample_{rep}.benchmark"
    output:
        assembled="abundance/scallop/{sense}_{seed}/sample_{rep}.gtf",
    shell:
        "scallop -i {input.bam} "
        "-o {output.assembled} &> {log}"
