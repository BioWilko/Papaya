#!/usr/bin/env python

import glob
import argparse
import pathlib
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument(
    "-b",
    "--barcode",
    help="Barcode as it appears within the 'fastq_pass' directory (e.g. 'barcode01')",
    type=str,
)
parser.add_argument(
    "-r",
    "--runpath",
    help="Path to run directory (the folder containing the sequencing summary)",
    type=pathlib.Path,
)
parser.add_argument(
    "-a",
    "--all",
    help="Return all barcodes to stdout separated by newlines",
    action="store_true",
)
parser.add_argument("-x", "--expunge", help="Parse a PAF file for reads with ref identity >= 90%", type=str)
args = parser.parse_args()


def fastq_list(runpath, barcode):
    fastq_path = os.path.join(runpath, f"fastq_pass/{barcode}/*.fastq")
    fastq_files = glob.glob(fastq_path)
    return fastq_files


def symlink_fastq(fastq_list):
    cwd = os.getcwd()
    for fastq in fastq_list:
        fastq_filename = fastq.split("/")[-1]
        symlink_path = os.path.join(cwd, fastq_filename)
        os.symlink(fastq, symlink_path)


def filter_barcodes(runpath):
    barcode_path = os.path.join(runpath, "fastq_pass/barcode*/")
    barcode_list = glob.glob(barcode_path)
    filtered_barcodes = []
    for barcode_dir in barcode_list:
        fastq_pattern = os.path.join(barcode_dir, "*.fastq")
        fastq_files = glob.glob(fastq_pattern)
        if len(fastq_files) > 1:
            filtered_barcodes.append(barcode_dir.split("/")[-2])
    return filtered_barcodes

def reads_to_expunge (paf_file):
    to_expunge = []
    with open(paf_file, "r") as f:
        reads = f.readlines()
    for read in reads:
        cols = read.split("\t")
        if int(cols[9]) / int(cols[10]) >= 0.8:
            to_expunge.append(cols[0])
    return set(to_expunge)
    
if args.expunge:
    for read in reads_to_expunge(args.expunge):
        sys.stdout.write(read + "\n")
    

if args.all:
    barcodes = filter_barcodes(args.runpath)
    for barcode in barcodes:
        sys.stdout.write(barcode + "\n")

if args.barcode:
    fastq_files = fastq_list(args.runpath, args.barcode)
    symlink_fastq(fastq_files)
