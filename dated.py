#!/usr/bin/env python

import argparse
import sys
import os
import warnings
import subprocess
import asyncio
import tempfile
import shutil
from multiprocessing import Process, Pool
from functools import partial
from Bio import SeqIO
from Bio.Phylo.PAML import codeml
from Bio.Align.Applications import ClustalwCommandline

"""DATED : DATing Evolutionary events and Divergence times."""


def read_sequence(filename):
    return SeqIO.to_dict(SeqIO.parse(filename, "fasta"))

def read_pairs(filename):
    ortho_pairs = []
    with open(filename, 'r') as pairs:
        for pair in pairs:
            ortho_pairs.append(pair.rstrip().split(","))
    return ortho_pairs

def calculate_ds(pairs, pep_seqdb, cds_seqdb):
    tmp_output_folder = tempfile.mkdtemp(prefix="dated_",dir="./")
    ds_func = partial(get_ds, pep_seqdb, cds_seqdb, tmp_output_folder)
    with Pool(processes=8) as pool:
        pool.map(ds_func, pairs)
    shutil.rmtree(tmp_output_folder)  # delete tmp output folder

def is_list(pairs):
    return isinstance(pairs,list) and len(pairs) == 2

def get_ds(pep_seqdb, cds_seqdb, tmp_output_folder, pairs):
    if(not is_list(pairs)):
        raise ValueError(f"Error: {pairs} is not a proper pair")
    [seqA, seqB] = pairs

    try:
        cwd = os.getcwd()
        tmp_folder = tempfile.mkdtemp(dir=tmp_output_folder)
        os.chdir(tmp_folder)
        pep_seq_path = "pep_pair.fasta"
        cds_seq_path = "cds_pair.fasta"
        aln_path = "pep_pair.aln"
        pal_path = "pep_pair.pal2nal"

        if(check_input_sequences(seqA, seqB, pep_seqdb, cds_seqdb)):
            run_clustalw(seqA, seqB, pep_seqdb, pep_seq_path)
            run_pal2nal(seqA, seqB, cds_seqdb, cds_seq_path,
                    aln_path, pal_path)
            ds = run_codeml(tmp_folder, pal_path)
            print(",".join(map(str,(seqA, seqB, ds))))
        else:
            print("WARNING: missing sequence in FASTA ", pairs)
    except Exception as e:
        print(e)
    finally:
        try:
            os.chdir(cwd)
            shutil.rmtree(tmp_folder)  # delete directory
        except OSError as exc:
            raise exc

def check_input_sequences(seqA, seqB, pep_seqdb, cds_seqdb):
    return seq_exists(seqA, pep_seqdb) and seq_exists(seqB, pep_seqdb) and \
           seq_exists(seqA, cds_seqdb) and seq_exists(seqB, cds_seqdb)


def get_exec_dir():
    # get the executable directory
    return os.path.dirname(os.path.abspath(__file__))

def seq_exists(seqA, seqdb):
    # check if a sequence exists in database
    return seqA in seqdb


def run_clustalw(seqA, seqB, seqdb, seq_path):
    with open(seq_path, "w") as seqfile:
        SeqIO.write(seqdb[seqA], seqfile, "fasta")
        SeqIO.write(seqdb[seqB], seqfile, "fasta")
    clustalw = ClustalwCommandline('clustalw', infile=seq_path)
    stdout, stderr = clustalw()

def run_pal2nal(seqA, seqB, seqdb, cds_seq_path, aln_path, pal_path):
    with open(cds_seq_path, "w") as seqfile:
        SeqIO.write(seqdb[seqA], seqfile, "fasta")
        SeqIO.write(seqdb[seqB], seqfile, "fasta")
    f = open(pal_path, "w")
    subprocess.run(["pal2nal.pl",aln_path, cds_seq_path,
                    "-nogap", "-output", "paml"], stdout=f)

def run_codeml(tmp_folder, pal_path):
    cml = codeml.Codeml(alignment=pal_path,
                        out_file="pair.ks",
                        tree="pep_pair.dnd")
    exec_dir = get_exec_dir()
    cml.read_ctl_file(exec_dir+"/config/codeml.ctl")
    results = cml.run().get("pairwise")
    prot1 = next(iter(results.values()))
    for prot2, attributes in prot1.items():
        ds = (attributes.get("dS"))
    return ds


def main():
    parser = argparse.ArgumentParser(
        description='Given a list of ')
    parser.add_argument('pep_filename', type=str)
    parser.add_argument('cds_filename', type=str)
    parser.add_argument('pair_list', type=str)
    args = parser.parse_args()
    pep_seqdb = read_sequence(args.pep_filename)
    cds_seqdb = read_sequence(args.cds_filename)
    pairs = read_pairs(args.pair_list)
    calculate_ds(pairs, pep_seqdb, cds_seqdb)




if __name__ == '__main__':
    main()