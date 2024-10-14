# seqdemu

Under development


## About

No nonse demultiplexer. 

## Features

- Demultiplex sequences with custom barcodes (both ends).
- Allow having barcodes in the middle of the string.
- Allow a mismatch.

## Installation

Dependencies etc.

```bash

# Clone the repository
git clone https://github.com/hsgweon/seqdemu.git
cd seqdemu
mamba create -n seqdemu_env -y -c conda-forge -c bioconda conda-forge::biopython progressbar2

## Activate seqdemu_ev
mamba activate seqdemu_env

## Add the directory to PATH or use the absolute path.

## To run
seqdemu.py -i gzipped_dorado_fastq_file -o output_filename -b barcode_file -m number_of_mismatch

## Run example data provided
cd test
../seqdemu.py -i seq.fastq.gz -o test_m0 -b ../barcodes.csv -m 0

