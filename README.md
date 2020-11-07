# Jumper

![Overview of Jumper](overview.png)
(a) Viruses in the order Nidovirales generate a set <img src="https://latex.codecogs.com/gif.latex?\mathcal{T}" /> of discontinuous transcripts with varying abundances <img src="https://latex.codecogs.com/gif.latex?\mathbf{c}" /> during infection.
(b) Next generation sequencing will produce an alignment <img src="https://latex.codecogs.com/gif.latex?\mathcal{R}" /> with two types of aligned reads: unphased reads that map to a contiguous genomic region (black) and phased reads that map to distinct genomic regions (red).
(c) From <img src="https://latex.codecogs.com/gif.latex?\mathcal{R}" /> we obtain the segment graph <img src="https://latex.codecogs.com/gif.latex?G" />, a directed acyclic graph with a unique Hamiltonian path. Jumper solves the Discontinuous Transciption Assembly problem to infer <img src="https://latex.codecogs.com/gif.latex?\mathcal{T}" /> and <img src="https://latex.codecogs.com/gif.latex?\mathbf{c}" /> with maximum likelihood.

## Contents

  1. [Pre-requisites](#pre-requisites)
  2. [Usage instcructions](#usage)
     * [I/O formats](#io)
     * [Jumper](#jumper)
     * [simulation](#simulation)

<a name="pre-requisites"></a>
## Pre-requisites
+ python3 (>=3.6)
+ [numpy](https://numpy.org/doc/)
+ [pysam](https://pysam.readthedocs.io/en/latest/)
+ [pandas](https://pandas.pydata.org/pandas-docs/stable/index.html)
+ [gurobipy](https://www.gurobi.com/documentation/9.0/quickstart_mac/py_python_interface.html)

<a name="usage"></a>
## Usage instructions

<a name="io"></a>
### I/O formats
The input for Jumper is bam file containing the sequencing data and a fasta file containing the reference genome.
The output is similar to a fasta file format, where each transcript name is followed by the edges in the corresponding path in the segment graph (see `data/sample_transcripts.out` for an example).

### Arguments
    usage: jumper.py [-h] [-b BAM] -f FASTA [-k NUMPATHS]
                     [--min-base-qual MIN_BASE_QUAL]
                     [--min-mapping-qual MIN_MAPPING_QUAL] [-w WIDTH]
                     [--samplingFrequency SAMPLINGFREQUENCY]
                     [--sj_threshold SJ_THRESHOLD]
                     [--phasing_threshold PHASING_THRESHOLD] [--greedy GREEDY]
                     [--outputCSV OUTPUTCSV] [--outputPhasing OUTPUTPHASING]
                     [--inputCSV INPUTCSV] [--inputPhasing INPUTPHASING]
                     [--inputBreakpoints INPUTBREAKPOINTS]
                     [--inputEdges INPUTEDGES] [--outputGraph OUTPUTGRAPH]
                     [--outputDOT OUTPUTDOT]
                     [--outputTranscripts OUTPUTTRANSCRIPTS]
                     [--outputBreakpoints OUTPUTBREAKPOINTS]
                     [--outputEdges OUTPUTEDGES]
                     [--outputDecomposition OUTPUTDECOMPOSITION]
                     [--outputMatching OUTPUTMATCHING] [--outputGTF OUTPUTGTF]
                     [--report REPORT] [--noverbose] [--threads THREADS]
                     [--timelimit TIMELIMIT] [--maxIter MAXITER]

    optional arguments:
      -h, --help            show this help message and exit
      -b BAM, --bam BAM     aligned bam file
      -f FASTA, --fasta FASTA
                            fasta file
      -k NUMPATHS           number of paths for the flow decomposition
      --min-base-qual MIN_BASE_QUAL
                            minimum base quality [20]
      --min-mapping-qual MIN_MAPPING_QUAL
                            minimum mapping quality [20]
      -w WIDTH, --width WIDTH
                            spliced junction width parameter [6]
      --samplingFrequency SAMPLINGFREQUENCY
                            number of sampling points for the likelihood function
      --sj_threshold SJ_THRESHOLD
                            minimum support for splicing junction [20]
      --phasing_threshold PHASING_THRESHOLD
                            coverage threshold for transcripts [0]
      --greedy GREEDY       set greedy flag to TRUE
      --outputCSV OUTPUTCSV
                            output csv file for sj reads
      --outputPhasing OUTPUTPHASING
                            output file containing phasing reads
      --inputCSV INPUTCSV   input csv file with sj reads
      --inputPhasing INPUTPHASING
                            input phasing file
      --inputBreakpoints INPUTBREAKPOINTS
                            input file containing breakpoints
      --inputEdges INPUTEDGES
                            input file containing graph edges
      --outputGraph OUTPUTGRAPH
                            output graph file
      --outputDOT OUTPUTDOT
                            output DOT file for splice graph
      --outputTranscripts OUTPUTTRANSCRIPTS
                            output file for transcripts
      --outputBreakpoints OUTPUTBREAKPOINTS
                            output file containing breakpoints
      --outputEdges OUTPUTEDGES
                            output file containing graph edges
      --outputDecomposition OUTPUTDECOMPOSITION
                            output file for the decomposed non-canonical
                            transcripts
      --outputMatching OUTPUTMATCHING
                            output file for the matching of phasing reads to
                            inferred transcripts
      --outputGTF OUTPUTGTF
                            output file in GTF format
      --report REPORT       output file for report on the splice graph
      --noverbose           do not output statements from internal solvers
                            [default is false]
      --threads THREADS     number of threads allowed to be used [1]
      --timelimit TIMELIMIT
                            time limt for the gurobi solvers in seconds [None]
      --maxIter MAXITER     maximum iterations for the greedy algorithm [100]

### Example
One way to see how to use Jumper is through `simulation_pipeline` for example cases of using Jumper on simulated bam files.
The Jumper usage is shown in the snakemake `simulation_pipeline/jumper.smk`.

Here we will run Jumper to reconstruct the transcripts on simulated phasing reads.

#### Simulate transcripts and phasing reads
    
    $ python simulate_segment_graph.py --seed 0 --sense neg --npaths 2 --inputBreakpoints ../data/sampleBreakpoints.out --inputEdges ../data/sampleEdges.out --outputPaths ../data/sample_transcripts.out --outputFasta ../data/sample_transcripts.fasta -f ../data/reference.fasta --outputReadCounts ../data/sample_readcounts.out --outputPhasing ../data/sample_phasing.out --nreads 1000

This command generates 1000 phasing reads in `../data/sample_phasing`.
The ground transcripts are written in `../data/sample_transcripts.out`.

#### Reconstruct the transcripts

    $ python jumper.py --inputBreakpoints ../data/sampleBreakpoints.out --inputEdges ../data/sampleEdges.out --inputPhasing ../data/sample_phasing.out --outputDecomposition ../data/sample_decomposition.out -k 50 -f ../data/reference.fasta --greedy True --outputMatching ../data/sample_matching.out > ../data/sample.log
  
The reconstructed transcripts are written to `../data/sample_decomposition.out`.

