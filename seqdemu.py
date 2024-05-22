#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser("Reads and writes each entry as a single file.")
parser.add_argument("-i",
                    action = "store", 
                    dest = "infile", 
                    metavar = "infile",
                    help = "[REQUIRED]", 
                    required = True)
parser.add_argument("-o",
                    action = "store",
                    dest = "outfile",
                    metavar = "outfile",
                    help = "[REQUIRED]",
                    required = True)
parser.add_argument("-b",
                    action = "store",
                    dest = "barcodes",
                    metavar = "barcodes",
                    help = "[REQUIRED]",
                    required = True)
parser.add_argument("-m",
                    action = "store",
                    dest = "mismatch",
                    metavar = "mismatch",
                    help = "[REQUIRED]",
                    required = True)
options = parser.parse_args()


import gzip
import multiprocessing
import os, shutil
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
import progressbar
from multiprocessing import Manager


def boyer_moore_search_with_exact_mismatches(text, pattern, exact_mismatches):
    n = len(text)
    m = len(pattern)
    results = []

    # Preprocess pattern for Boyer-Moore
    bad_character = {}
    for i in range(m - 1):
        bad_character[ord(pattern[i])] = m - i - 1

    # Search with Boyer-Moore
    i = 0
    while i <= n - m:
        j = m - 1
        mismatches = 0
        while j >= 0:
            if pattern[j] != text[i + j]:
                mismatches += 1
                if mismatches > exact_mismatches:
                    break
            j -= 1

        if j == -1 and mismatches == exact_mismatches:
            results.append((i, i + m - 1))
            i += 1
        else:
            i += max(bad_character.get(ord(text[i + m - 1]), m), m - j - 1)

    return results

def process_chunk(args):
    
    chunk, progress_counter, lock, temp_file_prefix_cpu = args

    temp_files_full = [open(os.path.join(temp_dir, f"{temp_file_prefix_cpu}_full_{i}.fasta"), "w") for i in range(len(list_samplenames))]
    temp_files_barc = [open(os.path.join(temp_dir, f"{temp_file_prefix_cpu}_barc_{i}.fasta"), "w") for i in range(len(list_samplenames))]
    temp_files_noba = [open(os.path.join(temp_dir, f"{temp_file_prefix_cpu}_noba_{i}.fasta"), "w") for i in range(len(list_samplenames))]
    temp_files_mult = [open(os.path.join(temp_dir, f"{temp_file_prefix_cpu}_mult_{i}.fasta"), "w") for i in range(len(list_samplenames))]

    exact_mismatches = int(options.mismatch)
    
    for records in chunk:

        for record in records:
            
            for i in range(len(list_samplenames)):
                
                for pair in [0, 1]:
                    barcode_F = list_barcodes[i][pair][0]
                    barcode_RC = list_barcodes[i][pair][1]
                    sequence = str(record.seq).upper()

                    result_F = boyer_moore_search_with_exact_mismatches(sequence, barcode_F, exact_mismatches)
                    result_RC = boyer_moore_search_with_exact_mismatches(sequence, barcode_RC, exact_mismatches)

                    if len(result_F) == 1 and len(result_RC) == 1:
                        # output original sequence
                        temp_files_full[i].write(f">{record.description}\n{sequence}\n")
                        # output sequences with barcode intact
                        temp_files_barc[i].write(f">{record.description}\n{sequence[result_F[0][0]:result_RC[0][1] + 1]}\n")
                        # output sequence with barcodes trimmed
                        temp_files_noba[i].write(f">{record.description}\n{sequence[result_F[0][1]+1:result_RC[0][0]]}\n")

                    elif (len(result_F) == 1 and len(result_RC) > 1) or (len(result_F) > 1 and len(result_RC) == 1):
                        # print(result_F, result_RC)
                        temp_files_mult[i].write(f">{record.description}\n{sequence}\n")

            with lock:
                
                progress_counter.value += 1
                pbar.update(progress_counter.value)


def divide_file_chunks(infile, chunk_size):
    
    chunks = []
    with gzip.open(infile, "rt") as f:
        records = list(SeqIO.parse(f, "fastq"))
        # print(records[2500])
        for i in range(0, len(records), chunk_size):
            chunks.append(records[i:i + chunk_size])
    return chunks


if __name__ == "__main__":
    
    list_barcodes = []
    list_samplenames = []

    # Count number of sequences in FASTQ
    # infile = gzip.open(options.infile, "rt")
    with gzip.open(options.infile, 'rb') as f:
        for i, l in enumerate(f):
            pass
    number_of_sequences = int((i+1)/4)

    manager = Manager()
    progress_counter = manager.Value('i', 0)
    lock = manager.Lock()

    # Barcodes
    with open(options.barcodes, "r") as barcodes:
        for line in barcodes:
            parts = line.rstrip().split("\t")
            list_barcodes.append([(parts[0], str(Seq(parts[1]).reverse_complement())), (parts[1], str(Seq(parts[0]).reverse_complement()))])
            list_samplenames.append(parts[2])

    print("Calculating how to divide FASTQ file into chunks and assigning to multiple CPUs...")

    # Number of CPUs to use
    num_cpus = multiprocessing.cpu_count()
    num_cpus = 100
    chunk_size = 1000 # Number of sequences in a chunk
    chunks = divide_file_chunks(options.infile, chunk_size)

    # Specify the name of the temporary directory
    temp_dir = "temp"
    shutil.rmtree(temp_dir, ignore_errors=True)  # Remove directory if it exists
    os.makedirs(temp_dir)

    # Start progressbar
    pbar = progressbar.ProgressBar(max_value=number_of_sequences).start()

    # Pool
    pool_args = [(chunks[i::num_cpus], progress_counter, lock, f"cpu{i}") for i in range(num_cpus)]

    with multiprocessing.Pool(processes=num_cpus) as pool:
        pool.map(process_chunk, pool_args)

    outfiles = {}
    for name in list_samplenames:
        outfiles[name] = {
            "full": open(options.outfile + "_" + name + "_full.fasta", "w"),
            "barc": open(options.outfile + "_" + name + "_barc.fasta", "w"),
            "noba": open(options.outfile + "_" + name + "_noba.fasta", "w"),
            "mult": open(options.outfile + "_" + name + "_mult.fasta", "w")
        }

    for trim in ["full", "barc", "noba", "mult"]:
        for i in range(len(list_samplenames)):
            for cpu in range(num_cpus):
                with open(os.path.join(temp_dir, f"cpu{cpu}_{trim}_{i}.fasta"), "r") as temp_file:
                    for line in temp_file:
                        if line.strip():
                            outfiles[list_samplenames[i]][trim].write(line)

    # Close all files
    for name in list_samplenames:
        for trim in ["full", "barc", "noba"]:
            outfiles[name][trim].close()

    shutil.rmtree(temp_dir, ignore_errors=True)
     