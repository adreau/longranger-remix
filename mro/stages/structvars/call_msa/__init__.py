#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import os
import numpy as np
import pandas as pd
import subprocess
import longranger.sv.io as sv_io
import tenkit.reference

__MRO__ = """
stage CALL_MSA(
    in  bam     possorted_bam,
    in  string  reference_path,
    in  bedpe   candidates,
    in  int     padding_abs,
    in  float   padding_fract,
    in  int     min_coverage,
    in  int     mismatch_penalty,
    in  int     min_overlap_score,
    in  float   min_flanking_cov,
    out bedpe   sv_calls,
    src py     "stages/structvars/call_msa",
) split using (
    in  bed chunk_bed,
)
"""


def split(args):
    sv_df = sv_io.read_sv_bedpe_to_df(args.candidates)

    if len(sv_df) == 0:
        chunk_defs = [{'chunk_bed': None}]
        return {'chunks': chunk_defs}

    types = [sv_io.get_sv_type(row.info) for (_, row) in sv_df.iterrows()]
    sv_df["type"] = types

    sv_df["size"] = sv_df["stop2"] - sv_df["start1"]
    sv_df["padding"] = np.round(sv_df["size"] * args.padding_fract + args.padding_abs).astype(np.int32)
    sv_df["new_start"] = np.maximum(0, sv_df["start1"] - sv_df["padding"])
    sv_df["new_stop"] = sv_df["stop2"] + sv_df["padding"]
    sv_df["new_size"] = sv_df["new_stop"] - sv_df["new_start"]
    sv_df = sv_df[(sv_df["new_size"] < 4000) | (sv_df["type"] == "DEL")]

    sv_df.sort(['chrom1', 'start1'], inplace=True)

    nsvs = sv_df.shape[0]
    nsvs_per_chunk = max(100, int(np.ceil(nsvs / 200.0)))
    nchunks = int(np.ceil(nsvs / float(nsvs_per_chunk)))
    chunk_defs = []

    for i in range(nchunks):
        chunk_start = i * nsvs_per_chunk
        chunk_end =  min(nsvs, (i + 1) * nsvs_per_chunk)
        subset = sv_df[chunk_start:chunk_end][["chrom1", "new_start", "new_stop"]]

        # Figure out correct padding
        fn = os.path.join(os.getcwd(), "chunk_%d.bed" % i)
        subset.to_csv(fn, header=False, sep="\t", index=False)
        chunk_defs.append({'chunk_bed': fn, "__mem_gb": 3.0})


    return {'chunks': chunk_defs}


def main(args, outs):

    rust_env = os.environ.copy()
    rust_env["RUST_BACKTRACE"] = "1"

    if args.chunk_bed is None:
        sv_io.write_sv_df_to_bedpe(None, outs.sv_calls)
        return

    # Run PVC
    fasta = tenkit.reference.get_fasta(args.reference_path)
    pvc_args = ['pvc', '--min-overlap', str(args.min_overlap_score), '--min-coverage', str(args.min_coverage), '--mismatch', str(args.mismatch_penalty), '--flanking-cov', str(args.min_flanking_cov)]
    pvc_args.extend(["call-bed", "-o", outs.sv_calls, fasta, args.possorted_bam, args.chunk_bed])
    subprocess.check_call(pvc_args, env=rust_env)


def join(args, outs, chunk_defs, chunk_outs):
    join_df = None
    for chunk in chunk_outs:
        df = sv_io.read_sv_bedpe_to_df(chunk.sv_calls)
        join_df = pd.concat([join_df, df], ignore_index = True)

    sv_io.write_sv_df_to_bedpe(join_df, outs.sv_calls)