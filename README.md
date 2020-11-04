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
+ [pysam](https://pysam.readthedocs.io/en/latest/)
+ [pandas](https://pandas.pydata.org/pandas-docs/stable/index.html)
+ [gurobipy](https://www.gurobi.com/documentation/9.0/quickstart_mac/py_python_interface.html)

<a name="usage"></a>
## Usage instructions

<a name="io"></a>
### I/O formats
The input for Jumper is bam file containing the sequencing data and a fasta file containing the reference genome.
