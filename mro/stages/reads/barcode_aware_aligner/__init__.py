#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Sort BAM file by sorting buckets, then concatenating bucket files
#

import subprocess
import multiprocessing
import glob
import os.path
import json
import tenkit.reference as tk_ref

__MRO__ = """
stage BARCODE_AWARE_ALIGNER(
    in  fastq.gz[] reads,
    in  string     reference_path,
    in  string     read_group,
    in  string     sample_id,
    in  bool       exclude_non_bc_reads,
    in  string[]   read_groups,
    out bam        bc_sorted_bam,
    out json       position_chunks,
    src py         "stages/reads/barcode_aware_aligner",
) split using (
    in  fastq.gz   barcode_reads,
)
"""
def split(args):
    mem_in_gb = 16
    threads = 4
    chunk_defs = [{'barcode_reads': read, '__mem_gb': mem_in_gb,'__threads':threads} for read in args.reads]
    return {'chunks': chunk_defs}

def main(args, outs):
    # run lariat
    out_base = os.path.dirname(outs.bc_sorted_bam)
    read_groups_arg = ','.join(args.read_groups) # TODO this string could get pretty big if there are too many read groups
    centromere_file = tk_ref.get_centromere_regions(args.reference_path)
    subprocess.check_call(['lariat',
                           '-reads='+args.barcode_reads,
                           '-read_groups='+read_groups_arg,
                           '-genome='+args.reference_path+"/fasta/genome.fa",
                           '-sample_id='+str(args.sample_id),
                           '-threads='+str(args.__threads),
                           '-centromeres='+centromere_file,
                           '-trim_length='+str(args.trim_length),
                           '-output='+out_base])

    # sort output
    outs.position_chunks = []
    cmds = []
    for pos_bucketed_bam in sorted(glob.glob(out_base+"/*pos_bucketed*")):
        outs.position_chunks.append(pos_bucketed_bam+"_pos_sorted.bam")
        cmds.append(['samtools','sort',pos_bucketed_bam,pos_bucketed_bam+"_pos_sorted"])

    pool = multiprocessing.Pool(args.__threads)
    try:
        pool.map(subprocess.check_call, cmds)
    finally:
        pool.close()
        pool.join()

def join(args, outs, chunk_defs, chunk_outs):
    outs.bc_sorted_bam = [chunk_out.bc_sorted_bam for chunk_out in chunk_outs]
    position_chunks = {}
    for chunk_out in chunk_outs:
        for index, filename in enumerate(chunk_out.position_chunks):
            position_chunks.setdefault(index,[])
            position_chunks[index].append(filename)
    with open(outs.position_chunks,'w') as out_file:
        json.dump(position_chunks, out_file)
