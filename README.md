# GenCore

## Overview

**GenCore** is a cutting-edge tool designed for the comparison of genomic sequences. 
Unlike traditional genomic comparison methods that rely on hashing techniques, **GenCore** utilizes Locally Consistent Parsing (LCP), a novel approach that processes strings recursively to identify cores adhering to the LCP rules. 
These cores enable highly accurate comparison of genomes to evaluate distances and ultimately construct phylogenetic trees.

## Features

- **Genome Comparison**: Efficiently compare genomes using distance metrics.

- **Distance Matrix Calculation**: Compute various similarity metrics including Jaccard, Dice, and normalized vector similarity.

- **Multi-threading Support**: Leverage multiple threads for faster processing.

- **Flexible Input Formats**: Support for FASTA (`.fa`), FASTQ (gzipped) (`.fq.gz`), and BAM (`.bam`) file formats.

- **Core Management**: Read and write core files for streamlined genomic analysis.

- **Customizable Parameters**: Options for adjusting the number of threads, LCP levels, and distance calculation methods.

## Getting Started

### Prerequisites

- A modern C++ compiler (e.g., g++ or clang++) capable of C++11 or later.

- GNU Make for building the program.

- Access to a Unix-like environment (Linux, macOS).

- `htslib`: Required for reading BAM and FASTQ files.

- `lcptools`: Required for locally consistent parsing.

- `biopython`, `matplotlib`, `PyQt5`, and `ete3`: Required for constructing phylogenetic tree from distance matrices.

### Installation

1) Download the source code to your local machine.
  
2) Extract the files (if compressed) and navigate to the project directory.

3) Navigate to `lcptools/program`.
    
4) Run `make` in the terminal within the project directory. This compiles the source code and generates an executable named `gencore`.

### Installation on Unix (Linux & macOS)

First, open your terminal. 
Then execute the following commands to shallow clone the repository and compile the source code:

- **Clone the Repository**:

```cpp
git clone --depth 1 --recursive https://github.com/BilkentCompGen/gencore.git
cd gencore
```

- **Build the Project**: Run the following command to compile the project and its dependencies:
```
make install
```

This will compile the `gencore` executable and install `htslib` and `lcptools` as required.

- **Clean Up**: To remove object files and the target binary, run:

```
make clean
```

These instructions assume that you have `git`, a C++ compiler, and `make` installed on your system. 

### Reinstalling Dependencies

If you need to reinstall `htslib` or `lcptools`, use:

```
make reinstall-htslib
make reinstall-lcptools
```

## Usage

The **GenCore** tool can be executed with various command-line options. Below are the primary usage patterns:

### Basic Command Structure


```
./gencore [OPTIONS]
```

### Options

- Read Cores:

```
-r              Read cores from specified files.
                Usage: ./gencore -r file1.cores,file2.cores
                or: ./gencore -r -f filenames.txt
```

- **File Formats**: Specify the format of the files you are processing:

```
[fa|fq|bam]     Execute program with specified files in the given format.
                Supported formats: [ fa | fq.gz | bam ]
                Usage: ./gencore fa ref1.fa,ref2.fa
                       ./gencore fq reads1.fq.gz,reads2.fq.gz
                       ./gencore bam aln1.bam,aln2.bam
```

- **Input File List**:

```
-f [filename]   Execute program with a file containing file names.
                Usage: ./gencore fa -f filenames.txt
```

- **LCP Level**:

```
-l [level]      Set LCP level. [Default: 4]
                Usage: ./gencore fa ref1.fa,ref2.fa -l 4
```

- **Number of Threads**:

```
-t [number]     Set number of threads. [Default: 8]
                Usage: ./gencore fa ref1.fa,ref2.fa -t 2
```

- **Set Type**:

```
[--set|--vec]   Choose calculation method based on set or vector of cores. [Default: vector]
                Usage: ./gencore fa ref1.fa,ref2.fa --set

```

- **Write Cores**:

```
-w [filenames]  Store cores processed from input files.
                Usage: ./gencore fa ref1.fa,ref2.fa -w ref1.cores,ref2.cores
                or: ./gencore fa ref1.fa,ref2.fa -w -f filenames.txt
```

## Input Files

The **GenCore** tool requires specific input files to process genomic data and compute distance matrices. 
Below are the details regarding the expected input file formats:

1) FASTA Files: Each file must be in FASTA format, where each sequence represents a genomic region or chromosome.where genome_short_name is a unique 

2) FASTQ Files: Each file must be in FASTQ format and should be compressed using gzip (.gz), where each sequence represents a genomic region.

3) BAM Files: Each file must be in BAM format, which is a binary version of the SAM format.

### File Input Options

You can provide input files using one of the following methods:

1) Command Line Arguments: Concatenate file names with commas, e.g., file1.fasta,file2.fastq,file3.bam.

2) File Option: Use the -f option to specify a file containing the names of input files, where each file name should be on a separate line. For example:

    ```
    file1.fasta
    file2.fasta
    file3.fasta
    ```

## Program Outputs

The **GenCore** tool generates several output files containing distance matrices calculated from the genomic cores processed. 
The outputs are saved in the specified prefix format, and the following files will be created based on the provided command-line arguments:

### Output Files

1) **Dice Distance Matrix**:

  - Filename: `<prefix>.dice.phy`

  - Format: The first line contains the number of genomes followed by a matrix of Dice distances. Each subsequent line starts with the short name of the genome, followed by the distances from that genome to all other genomes. The distance values are represented as floating-point numbers.

2) **Jaccard Similarity Matrix**:

  - Filename: `<prefix>.jaccard.phy`

  - Format: Similar to the Dice distance matrix, this file contains the number of genomes on the first line, followed by the Jaccard similarity values. Each genome’s line starts with its short name and is followed by its similarities to all other genomes. Similar to the Dice matrix, these values are also represented as floating-point numbers.

2) **Normalized Vector Similarity Matrix**:

  - Filename: `<prefix>.ns.phy`

  - Format: This file contains the number of genomes in the first line, followed by the normalized vector similarity distances. Each genome’s line starts with its short name and is followed by the similarity distances to all other genomes.

### Note 

The distances and similarities computed are useful for phylogenetic analysis, allowing researchers to understand the relationships and evolutionary distances between different genomes. 
Ensure to specify a valid prefix to easily manage output files and facilitate subsequent analyses.

Output files contain distances which is calculated by subtracting the similarity score from 1.

### Similarity Metrics

* `Jaccard Similarity`: This metric measures the similarity between the two genomes based on the intersection over the union of their features. A value closer to 1 indicates higher similarity. 

* `Dice Similarity`: Similar to the Jaccard Similarity, the Dice Similarity measures the overlap between two genomes but considers the size of the two sets in its calculation, leading to a more sensitive measure in certain contexts.

* `Distance Based Similarity`: This metric evaluates the similarity based on a distance function, with a value closer to 1 indicating minimal distance or high similarity between the genomes. It also include normalization of the given genomes having different read depths.

#### Distance Based Similarity

The Distance Based Similarity between two genomic sequences, `G1` and `G2`, is determined through a formula that considers the occurrences of each LCP core within the genomes and their sequencing depths. 
This approach aims to offer a detailed perspective on similarity, factoring in not only the presence of cores but also their abundance in relation to sequencing depth.

#### Formula Breakdown

Given two genomic sequences, `G1` and `G2`, the similarity is calculated using the following formula:

$$
1 - \frac{\sum_{i=1}^{n} | d_2 \cdot c_1^i - d_1 \cdot c_2^i |}{\sum_{i=1}^{n} d_2 \cdot c_1^i + d_1 \cdot c_2^i}
$$

Where:
- \( $n$ \): The upper limit of the summation (number of cores).
- \( $d_1, d_2$ \): represent the sequencing depths of genomes `G1` and `G2`, respectively.
- \( $c_1^i, c_2^i$ \): Corresponds to the number of occurrences of core `i` in genomes, `G1` and `G2`, respectively.


#### Interpretation

This formula calculates the similarity by evaluating the absolute difference in the weighted occurrences of each core between the two genomes, normalized by the sum of these occurrences and subtracting from *1*. 
The weights are assigned based on the sequencing depths of the genomes, allowing the similarity score to reflect variations in sequencing efforts accurately.

* A Similarity score close to 1 indicates a high degree of similarity, implying minimal differences in the adjusted occurrences of cores between the genomes.
* A score nearing 0 suggests greater disparity between the genomic sequences, as indicated by significant differences in the weighted core occurrences.

## License

**GenCore** is released under BSD 3-Clause License.

For more details, please visit the [license file](https://github.com/BilkentCompGen/gencore/blob/main/LICENSE).
