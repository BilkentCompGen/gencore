# GENCORE

## Overview

**gencore** is a bioinformatics tool designed for the comparison of genomic sequences. 
Unlike traditional genomic comparison methods that rely on hashing techniques, **gencore** utilizes [*Locally Consistent Parsing (LCP)*](https://github.com/BilkentCompGen/lcptools), a novel approach that processes strings iteratively. 
This technique allows for the paritioning genomic data into small modules, called cores, leading to highly accurate comparisons when evaluating genetic distances and constructing phylogenetic trees.

## Features

- **Genome Comparison**: Efficiently compare genomes using distance metrics.

- **Distance Matrix Calculation**: Compute similarity metrics including Jukes-Cantor model.

- **Multi-threading Support**: Leverage multiple threads for faster processing.

- **Flexible Input Formats**: Supports FASTA (`.fa`), FASTQ (`fq`/`.fq.gz`).

## Getting Started

### Prerequisites

- A modern C compiler (e.g., gcc).

- GNU Make for building the program.

- Access to a Unix-like environment (Linux, macOS).

- `biopython`: Required for constructing phylogenetic tree from distance matrices.

---

### Installation on Unix (Linux & macOS)

Open your terminal and execute the following commands to clone the repository and compile the source code:

- **Clone the Repository**:

```bash
# clone the repository
git clone --depth 1 --recursive https://github.com/BilkentCompGen/gencore.git
cd gencore

# install dependencies and build the project
make install
```

## Usage

The general usage format is:

```bash
./gencore [PROGRAM] [OPTIONS]
```

### [PROGRAM]

1. **fa**: Processing assembled genomes.
2. **fq**: Processing genome reads.
3. **ld**: Processing precomputed cores.

For detailed options for each program, see the sections below.

---

### `fa`: Processing Assembled Genomes

```bash
./gencore fa [OPTIONS]
```

#### Options:

- **`-i [filename]`**: The file containing filenames of genome files (one per line).

- **`-l [num]`**: LCP-level (default: 4).

- **`-t [num]`**: Number of threads (default: 8).

- **`--min-cc [num]`**: Minimum frequency (core count) for a core (default: 1).

- **`--min-cc-file [filename]`**: File containing minimum core count thresholds for each input file (one per line).

- **`--max-cc [num]`**: Maximum frequency (core count) for a core (default: UINT64\_MAX).

- **`--max-cc-file [filename]`**: File containing maximum core count thresholds for each input file (one per line).

- **`[--set|--vec]`**: Compute distances based on set or vector of cores (default: set).

- **`-o [filename]`**: Output file to store cores.

- **`-p [prefix]`**: Prefix for the result files (default: gc).

- **`-s [filename]`**: File containing short names for input files (one name per line). Default is the first 10 characters of input filenames.

- **`-v`**: Enable verbose output (default: false).

---

### `fq`: Processing Genome Reads

```bash
./gencore fq [OPTIONS]
```

#### Options:

- **`-i [filename]`**: The file containing filenames of genome read files (one per line).

- **`-l [num]`**: LCP-level (default: 4).

- **`-t [num]`**: Number of threads (default: 8).

- **`--min-cc [num]`**: Minimum frequency (core count) for a core (default: 15).

- **`--min-cc-file [filename]`**: File containing minimum core count thresholds for each input file (one per line).

- **`--max-cc [num]`**: Maximum frequency (core count) for a core (default: 256).

- **`--max-cc-file [filename]`**: File containing maximum core count thresholds for each input file (one per line).

- **`[--set|--vec]`**: Compute distances based on set or vector of cores (default: set).

- **`-o [filename]`**: Output file to store cores.

- **`-p [prefix]`**: Prefix for the result files (default: gc).

- **`-s [filename]`**: File containing short names for input files (one name per line). Default is the first 10 characters of input filenames.

- **`-v`**: Enable verbose output (default: false).

---

### Example Command

To process assembled genomes listed in `genome_files.txt` with default settings:

```bash
./gencore fa -i genome_files.txt -o output.cores -p results_prefix
```

- `genome_files.txt`: Contains the genome filenames, one per line.
- `output.cores`: File to store the generated cores.
- `results_prefix`: Prefix for the result files (e.g., distance matrices).

If you want to use short names for the genomes from `short_names.txt`:

```bash
./gencore fa -i genome_files.txt -o output.cores -p results_prefix -s short_names.txt
```

For processing genome reads with specific core count limits:

```bash
./gencore fq -i reads_files.txt --min-cc 5 --max-cc 100 -o reads.cores -p reads_prefix
```

---

## Input Data Format

- **Input Files**: Text files where each line corresponds to a file or argument for the respective genome or read.
  - Example for `genome_files.txt`:
    ```
    genome1.fasta
    genome2.fasta
    genome3.fasta
    ```
  - The same line index in all input files corresponds to the same genome or read argument. For example:
    - Line 1 of `genome_files.txt` corresponds to Line 1 of `short_names.txt`.

The example files are provided under [example](https://github.com/BilkentCompGen/lcptools/tree/main/example) folder.

---

## Program Outputs

The **gencore** tool generates several output files containing distance matrices calculated from the genomic cores processed. 
The outputs are saved in the specified prefix format, and the following files will be created based on the provided command-line arguments:

### Output Files

1) **Dice Distance Matrix**:

  - Filename: `gc.set.dice.lvl4.phy`

  - Formula: $\text{Dice}(A,B) = 1 - \frac{2\times|A\cap B|}{|A|+|B|}$

  - Format: The first line contains the number of genomes, followed by a matrix of Dice distances. Each subsequent line starts with the short name of the genome, followed by the distances from that genome to all other genomes. The distance values are represented as floating-point numbers.

2) **Jaccard Distance Matrix**:

  - Filename: `gc.set.jaccard.lvl4.phy`

  - Formula: $\text{Jaccard}(A,B) = 1 - \frac{|A \cap B|}{|A \cup B|}$
  
  - Format: Similar to the Dice distance matrix, this file contains the number of genomes on the first line, followed by the Jaccard distance values. Each genome’s line starts with its short name and is followed by its distances to all other genomes. Similar to the Dice matrix, these values are also represented as floating-point numbers. The names are as the default values: prefix, mode, and lcp level.

3) **Jukes-Cantor Correction Matrix**:

  - Filename: `gc.set.jc.lvl4.phy`

  - Formula: This matrix is based on the Dice similarity values and includes corrections using the Jukes-Cantor model. The steps are as follows:

    1. **Calculate Hamming Distance**:  
      $\text{Hamming Distance} = 1 - (\text{Dice Similarity})^{\frac{1}{\text{Average Length}}}$

      (Where "Average Length" is the average size of cores being compared.)

    2. **Apply Jukes-Cantor Correction**:  
      $\text{Jukes-Cantor}(A,B) = -\frac{3}{4} \log\left(1 - \text{Hamming Distance} \times \frac{3}{4}\right)$

  - Format: The first line contains the number of genomes, followed by a matrix of Jukes-Cantor corrected distances. Each subsequent line starts with the short name of the genome, followed by the corrected distances from that genome to all other genomes. These values are represented as floating-point numbers.

---

## Additional Command for Phylogenetic Tree Construction:

```bash
python3 phylowizard.py gc.set.dice.lvl4.phy
```

The script outputs a Newick format tree files, `gc.set.dice.lvl4.upgma.newick` and `gc.set.dice.lvl4.nj.newick`, which contains the tree structures, constructed with upgma and neighbor joining algorithm, in a textual format that can be used for further analysis or visualization.


## License

**gencore** is released under the BSD 3-Clause License, which allows for redistribution and use in source and binary forms, with or without modification, under certain conditions. For more detailed terms, please refer to the [license file](https://github.com/BilkentCompGen/gencore/blob/main/LICENSE).
