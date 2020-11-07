# Simulation Pipeline

The simulation pipeline generates short-read paired end bam files and runs three transcript assembly methods -- [Scallop](https://github.com/Kingsford-Group/scallop), [StringTie](https://github.com/gpertea/stringtie) and [Jumper](https://github.com/elkebir-group/Jumper) on the simulated bam files.
The transcripts and their abundances are simulated using a discontinuous transcription model.
The RNA-seq experiments are simulated using [polyester](https://github.com/alyssafrazee/polyester) and reads are aligned using [STAR](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html).
The resulting bam files then serve as input for the transcript assembly methods.
The results of the assembly methods can be compared using the `evaluation.py` code provided in the `jumper` folder.

Simulated bam files used to generate the figures in the paper 'JUMPER: Discontinuous Transcript Assembly in SARS-CoV-2' can be found at `https://databank.illinois.edu/datasets/IDB-6667667`.

The simulation parameters can be set in `config.yaml`.
Simulation parameters are
PARAMETER           | DESCRIPTION
--------------------|-------------
`sampleBreakpoints` | Junction locations of the input segment graph
`sampleEdges`       | Edges of the input segment graph
`npaths`            | Total number of transcript molecules (can be identical) generated in the simulation
`nreads`            | Total number of reads generated from the simulated transcripts
`nreplicates`       | Number of sequencing experiments for each simulated instance of transcripts and their abundances

## Steps of the pipeline

  1. [Simulation of transcripts and sequencing reads](#simulate)
  2. [Alignment](#align)
  3. [Assembly](#assembly)
     * [Scallop and StringTie](#existing)
     * [Jumper](#jumper)

<a name="simulate"></a>
## Simulation of transcripts and sequencing reads

    $ snakemake -j1

<a name="align"></a>
## Alignment of reads

    $ snakemake -j1 -s align.smk

<a name="assembly"></a>
## Assembly

<a name="existing"></a>
### Scallop and StringTie

    $ snakemake -j1 -s abundance.smk

<a name="jumper"></a>
### Jumper

    $ snakemake -j1 -s jumper.smk

