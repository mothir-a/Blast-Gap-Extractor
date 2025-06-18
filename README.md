# Blast-Gap-Extractor: Uncovered Region Extractor from All-vs-All BLASTn Results

## Overview

This script identifies *uncovered (non-overlapping)* regions in a set of nucleotide genomes, based on the results of an all-vs-all BLASTn search. These regions represent stretches of each genome that do not show significant similarity to any of the others in the dataset.

This can be useful for:

- Identifying strain- or species-specific regions
- Designing unique primers or probes
- Exploring viral or microbial genomic diversity

## Input

### 1. BLASTn Output (`outfmt 6`)

An all-vs-all BLASTn comparison should be performed in advance using `makeblastdb` and `blastn` as shown below:

```zsh
cat *.fasta > all_genomes.mfa
```

# Make BLAST database
makeblastdb -in all_genomes.mfa -dbtype nucl

# Execute BLASTn
blastn \
  -query all_genomes.mfa \
  -db all_genomes.mfa \
  -out all_blastn_out.tsv \
  -outfmt 6

** Note:** Replace all_genomes.mfa and all_blastn_out.tsv with appropriate filenames for your data.

### 2. Genome Length Table
Prepare a tab-separated file listing each genome's accession ID and its length. The format is:

acc.ver    length    gc%    gcsi
PanV_1     184321    27.1   35.7
PanV_2     176543    28.2   34.8
...

No header is required. Only the first two columns are used.

acc.ver: Genome ID (must match query and subject fields in the BLAST output)
length: Genome length in base pairs

## How to Use

  1. Place the Python script in the same directory as your input files.
  2. Adjust the input file names inside the script if needed.
  3. Run the script:

```zsh
python uncovered_region_extractor.py
```

The script performs the following:

- Reads the BLASTn result file (all_blastn_out.tsv)
- Reads the genome length file (all_genomes_seqinfo.tsv)
- Removes self-comparisons from the BLAST results
- Calculates uncovered regions (regions not covered by any alignment with other genomes)
- Outputs these regions in a tab-separated format

## Output

The script will create a file named uncovered_regions.tsv with the following format:

genome    start    end
PanV_1    1        523
PanV_1    15423    15788
PanV_2    9873     11542
...

Each row indicates a region of the genome that is not covered by any BLAST hit to other genomes.

## Example: Extracting Uncovered Sequences >100 bp

You can extract specific uncovered regions using extractseq (part of the EMBOSS suite). For example:

```zsh
cat uncovered_regions.tsv | \
awk '/PanV_1/{print}' | \
awk '$3-$2>100{print}' | \
awk 'BEGIN{print "extractseq PanV_1.fasta uncovered_regions_100_PanV_1_extract.mfa -regions \""}\
{print $2"-"$3" \\"}\
END{print "\" -stdout -separate"}' > uncovered_regions_100_PanV_1_extract.sh
```

This generates a shell script to extract regions longer than 100 bp from PanV_1.fasta.

You can then run the generated script:

```zsh
sh uncovered_regions_100_PanV_1_extract.sh
```

## Dependencies

- Python 3.x
- pandas library

Install the required Python package (if not already installed)

```zsh
pip install pandas
```

## License
- GNU GENERAL PUBLIC LICENSE Version 3.

## Contact Information
- motohiro-akashi[at]st.seikei.ac.jp (M.A.)

**Note:** You can reach out for support or questions related to the script.

## Citation
Motohiro Akashi, Masaharu Takemura and Seiichi Suzuki, Genome-Wide Characterization of Non-Shared Sequences among Amphora-Shaped Giant Viruses, 2025, XXXXXXXX
