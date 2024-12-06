# GenCore

## Overview

**GenCore** is a cutting-edge tool designed for the comparison of genomic sequences. 
Unlike traditional genomic comparison methods that rely on hashing techniques, **GenCore** utilizes [*Locally Consistent Parsing (LCP)*](https://github.com/BilkentCompGen/lcptools), a novel approach that processes strings recursively. 
This technique allows for the identification of genomic cores that are more representative in sequence composition and structure, leading to highly accurate comparisons when evaluating genetic distances and constructing phylogenetic trees.

## Features

- **Genome Comparison**: Efficiently compare genomes using distance metrics.

- **Distance Matrix Calculation**: Compute various similarity metrics including Jaccard and Dice.

- **Multi-threading Support**: Leverage multiple threads for faster processing.

- **Flexible Input Formats**: Supports FASTA (`.fa`), gzipped FASTQ (`.fq.gz`), and BAM (`.bam`) file formats.

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
    
4) Run `make install` in the terminal within the project directory. This compiles the source code and generates an executable named `gencore`.

### Installation on Unix (Linux & macOS)

Open your terminal and execute the following commands to clone the repository and compile the source code:

- **Clone the Repository**:

```bash
# Clone the repository
git clone --depth 1 --recursive https://github.com/BilkentCompGen/gencore.git
cd gencore

# Install dependencies and build the project
make install
# This installs required dependencies (htslib, lcptools) and compiles the GenCore executable.
```

- **Clean Up**: To remove object files and the target binary, run: *`make clean`*

These instructions assume that you have `git`, a C++ compiler, and `make` installed on your system. 

### Reinstalling Dependencies

If you need to reinstall `htslib`, use *`make reinstall-htslib`*; tp reinstall `lcptools`, use *`make reinstall-lcptools`*

## Usage

The **GenCore** tool can be executed with various command-line options. Below are the primary usage patterns:

### Basic Command Structure


```
./gencore [OPTIONS]
```

### Options

- **Read Cores**: Read mode which is reading pre-processed cores from the file.

```
-r [filename]   Read cores from specified files.
                Usage: ./gencore -r file1.lcpt,file2.lcpt
                       ./gencore -r -f filenames.txt
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
-l [num]        Set LCP level. [Default: 4]
                Usage: ./gencore fa ref1.fa,ref2.fa -l 4
```

- **Number of Threads**:

```
-t [num]        Set number of threads. [Default: 8]
                Usage: ./gencore fa ref1.fa,ref2.fa -t 2
```

- **Minimum Core Count Value**:

```
--min-cc [num]  Set the minimum frequency (core count) for a core to be included in the similarity analysis. 
                Only cores with a frequency greater than or equal to the specified value will be considered.
                You can give multiple threashold in txt file that will be set to each genomes argument. 
                [Default: 1]
                Usage: ./gencore fa ref1.fa,ref2.fa --min-cc 2
                   or  ./gencore fa ref1.fa,ref2.fa --min-cc -f min_threasholds.txt
```

- **Maximum Core Count Value**:

```
--max-cc [num]  Set the maximum frequency (core count) for a core to be included in the similarity analysis. 
                Only cores with a frequency less than or equal to the specified value will be considered. 
                You can give multiple threashold in txt file that will be set to each genomes argument. 
                [Default: 4294967295]
                Usage: ./gencore fa ref1.fa,ref2.fa --max-cc 100
                   or  ./gencore fa ref1.fa,ref2.fa --max-cc -f max_threashols.txt
```

- **Set Type**:

```
[--set|--vec]   Specify the mode for treating LCP cores. [--set] focuses only on presence or absence.
                [--vec] includes frequencies includes in the similarity calculation. [Default: set]
```

- **Write Cores**:

```
-w [filename]   Store cores processed from input files.
                Usage: ./gencore fa ref1.fa,ref2.fa -w ref1.lcpt,ref2.lcpt
                       ./gencore fa ref1.fa,ref2.fa -w -f filenames.txt
```

- **Prefix**:

```
-p [prefix]     Set the prefix for the output of the similarity matrices results. [Default: gc] 
                Usage: ./gencore fa ref1.fa,ref2.fa -p simu
```

- **Short Names**:

```
-s [filename]   Set short names of input files. Default is first 10 characters of input file names.
                Usage: ./gencore fa ref1.fa,ref2.fa -s ref1,ref2
                       ./gencore fa ref1.fa,ref2.fa -s -f shortnames.txt
```

- **Verbose**:

```
-v              Verbose. [Default: false]
                Usage: ./gencore fa ref1.fa,ref2.fa -v
```

### Sample Command

```bash
./gencore fa ref1.fa,ref2.fa -t 4 -l 5 -p output_prefix -v
```

#### Explanation of the Command:

- `fa ref1.fa,ref2.fa`: This specifies the file format for the input files. In this case, fa stands for FASTA files. You could also use fq for FASTQ files or bam for BAM files, depending on the input format. Input files containing the genomic sequences in FASTA format that you want to compare. You can list multiple files, separated by commas, as shown in this example.
- `-t 4`: This option sets the number of threads to use. In this example, it specifies that GenCore should use 4 threads for parallel processing. This can help speed up the computation, especially when working with large datasets.
- `-l 5`: This option sets the LCP (Locally Consistent Parsing) level. The LCP level controls how deeply the tool looks into the sequences. A higher LCP level means more detailed analysis, it may also decrease processing time as less number of cores are compared while this will put more distance in between genomes. The default value is 4, but here it’s set to 5 for a more thorough comparison.
- `-p output_prefix`: This option sets a prefix for the output files. All the generated distance matrices and other output files will start with this prefix. In this case, the prefix is output_prefix, so the generated files might look like,`output_prefix.set.dice.lvl5.phy` and `output_prefix.set.jaccard.lvl5.phy`.
- `-v`: This option enables verbose output, which will display additional information about the processing steps. It’s useful for debugging or for gaining insights into the tool's progress and performance during execution.

## Input Files

The **GenCore** tool requires specific input files to process genomic data and compute distance matrices. 
Below are the details regarding the expected input file formats:

1) **FASTA Files**: Each file must be in FASTA format, where each sequence represents a genomic region or chromosome.

2) **FASTQ Files**: Each file must be in FASTQ format and should be compressed using gzip (`.gz`), where each sequence represents a genomic region.

3) **BAM Files**: Each file must be in BAM format, which is a binary version of the SAM format.

### File Input Options

You can provide input files using one of the following methods:

1) ``Command Line Arguments``: Provide the file names separated by commas, e.g., `ref1.fa,ref2.fa,ref3.fa`.

2) ``File Option``: Use the `-f` option to specify a file containing the names of input files, where each file name should be on a separate line. For example:

    ```
    ref1.fa
    ref2.fa
    ref3.fa
    ```

## Program Outputs

The **GenCore** tool generates several output files containing distance matrices calculated from the genomic cores processed. 
The outputs are saved in the specified prefix format, and the following files will be created based on the provided command-line arguments:

### Output Files

1) **Dice Distance Matrix**:

  - Filename: `gc.set.dice.lvl4.phy`

  - Formula: $DICE(A,B) = \frac{2\times|A\cap B|}{|A|+|B|}$

  - Format: The first line contains the number of genomes, followed by a matrix of Dice distances. Each subsequent line starts with the short name of the genome, followed by the distances from that genome to all other genomes. The distance values are represented as floating-point numbers.

2) **Jaccard Similarity Matrix**:

  - Filename: `gc.set.jaccard.lvl4.phy`

  - Formula: $JACCARD(A,B) = \frac{|A \cap B|}{|A \cup B|}$
  
  - Format: Similar to the Dice distance matrix, this file contains the number of genomes on the first line, followed by the Jaccard similarity values. Each genome’s line starts with its short name and is followed by its similarities to all other genomes. Similar to the Dice matrix, these values are also represented as floating-point numbers. The names are as the default values: prefix, mode, and lcp level.

These matrices are saved in `.phy` format, which is commonly used for storing distance matrices in phylogenetic analysis.

### Note 

The distances and similarities computed are useful for phylogenetic analysis, allowing researchers to understand the relationships and evolutionary distances between different genomes. 
Ensure to specify a valid prefix to easily manage output files and facilitate subsequent analyses.

Output files contain distances, which are calculated by subtracting the similarity score from 1.

### Similarity Metrics

* `Jaccard Similarity`: This metric measures the similarity between the two genomes based on the intersection over the union of their features. A value closer to 1 indicates higher similarity. 

* `Dice Similarity`: Similar to the Jaccard Similarity, the Dice Similarity measures the overlap between two genomes but considers the size of the two sets in its calculation, leading to a more sensitive measure in certain contexts.

## Additional Command for Phylogenetic Tree Construction:

```bash
python3 phylowizard.py output_prefix.set.dice.lvl5.phy
```

### Explanation of the Command:
  
  - ``python3``: This invokes Python 3 to run the script. Make sure Python 3 is installed on your system. If you're working within a virtual environment, activate it beforehand.
  
  - ``phylowizard.py``: This is the Python script that will generate the phylogenetic tree. It takes a distance matrix in `.phy` format as input and produces both a Newick format tree and a PNG image of the tree.
  
  - ``output_prefix.set.dice.lvl5.phy``: This file is the distance matrix, typically generated from sequence comparisons using a tool like **GenCore**. It is in `.phy` format, which contains the distances between the genomes or sequences being analyzed. The `output_prefix` should be replaced with the actual prefix used during the generation of the matrix.

### What This Command Does:

1. The `phylowizard.py` script reads the `.phy` distance matrix file.
2. It generates the phylogenetic tree based on this distance matrix.
3. The script outputs two files:
    - A Newick format tree file, which contains the tree structure in a textual format that can be used for further analysis or visualization.
    - A PNG image of the phylogenetic tree, which provides a visual representation of the evolutionary relationships.

## License

**GenCore** is released under the BSD 3-Clause License, which allows for redistribution and use in source and binary forms, with or without modification, under certain conditions. For more detailed terms, please refer to the [license file](https://github.com/BilkentCompGen/gencore/blob/main/LICENSE).
