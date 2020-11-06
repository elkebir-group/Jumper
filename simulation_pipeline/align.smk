configfile: "config.yaml"

rule all:
    input:
        config["ref_dir"] + "/SAindex",
        expand("bam/star/{sense}_{seed}/sample_{rep}.Aligned.out.bam", sense=config["senses"], seed=config["seeds"], rep=config["rep"]),
        expand("bam/star/{sense}_{seed}/sample_{rep}.Aligned.sortedByCoord.out.bam.bai", sense=config["senses"], seed=config["seeds"], rep=config["rep"]),


rule STAR_index:
    input:
        ref=config["ref"],
        gtf=config["gtf_annot"],
    output: config["ref_dir"] + "/SAindex"
    params:
        length=config["read_length"]-1,
        path=config["ref_dir"],
    log: "log/star_index.log"
    threads: 24
    shell:
        "STAR --runMode genomeGenerate --genomeDir {params.path} "
        "--genomeFastaFiles {input.ref} --sjdbGTFfile {input.gtf} "
        "--sjdbOverhang {params.length} --runThreadN {threads} &> {log} "


rule STAR_align:
    input:
        fasta1="reads/{sense}_{seed}/sample_{rep}_1.fasta",
        fasta2="reads/{sense}_{seed}/sample_{rep}_2.fasta",
    threads: 8
    params:
        ref_dir=config["ref_dir"],
        sample_name="{sense}_{seed}_{rep}",
        outprefix="bam/star/{sense}_{seed}/sample_{rep}.", # cannot be an output, snakemake will wait for this file
    log:
        std="log/star_{sense}_{seed}_sample_{rep}.log",
        err="log/star_{sense}_{seed}_sample_{rep}.err.log",
    output:
        bam="bam/star/{sense}_{seed}/sample_{rep}.Aligned.out.bam",
        splice="bam/star/{sense}_{seed}/sample_{rep}.SJ.out.tab",
    shell:
        "STAR --runThreadN {threads} --genomeDir {params.ref_dir} "
        "--readFilesIn {input.fasta1} {input.fasta2} "
        "--twopassMode Basic "
        "--outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 "
        "--outSJfilterOverhangMin 12 12 12 12--outSJfilterCountUniqueMin 1 1 1 1 "
        "--outSJfilterCountTotalMin 1 1 1 1 --outSJfilterDistToOtherSJmin 0 0 0 0 "
        "--outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --scoreGapNoncan -4 "
        "--scoreGapATAC -4 --chimOutType WithinBAM HardClip --chimScoreJunctionNonGTAG 0 "
        "--alignSJstitchMismatchNmax -1 -1 -1 -1 --alignIntronMin 20 --alignIntronMax 1000000 "
        "--alignMatesGapMax 1000000 "
        "--outFileNamePrefix {params.outprefix} "
        "--outSAMtype BAM Unsorted "
        "--outSAMattrRGline ID:{params.sample_name} LB:{params.sample_name} PL:None SM:{params.sample_name} > {log.std} 2> {log.err}"

rule sort_star:
    input: "bam/star/{sense}_{seed}/sample_{rep}.Aligned.out.bam",
    output: "bam/star/{sense}_{seed}/sample_{rep}.Aligned.sortedByCoord.out.bam",
    threads: 4
    shell: "samtools sort -@ {threads} {input} -o {output} && samtools index {output}"

rule star_index:
    input:
        bam="bam/star/{sense}_{seed}/sample_{rep}.Aligned.sortedByCoord.out.bam",
    output:
        bai="bam/star/{sense}_{seed}/sample_{rep}.Aligned.sortedByCoord.out.bam.bai",
    shell:
        "samtools index {input.bam}"
