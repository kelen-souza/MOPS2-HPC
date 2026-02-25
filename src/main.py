from pycompss.api.api import compss_wait_on, compss_barrier
from apps import *
import argparse
import os
from pathlib import Path


def main(input_alignment_file, base_work_dir, max_sequences, pastar_threads, similar):
    align_out = "alignment.00.txt"
    sequence_dir = compss_wait_on(create_dir(
        os.path.join(base_work_dir, "sequences")))
    split_files = compss_wait_on(split_sequences(
        input_alignment_file, sequence_dir))
    pairs_rootdir = compss_wait_on(create_dir(
        os.path.join(base_work_dir, "pair_sequences")))
    metrics = list()
    folders_to_create = list()
    for i in range(len(split_files)):
        for j in range(i + 1, len(split_files)):
            seq1 = split_files[i]
            seq2 = split_files[j]
            pair_id = seq1.split(".")[-2] + "_" + seq2.split(".")[-2]
            folders_to_create.append(create_dir(
                os.path.join(pairs_rootdir, pair_id)))
    # sync after creating all the folders to avoid error on pycompss
    compss_barrier()
    for i in range(len(split_files)):
        for j in range(i + 1, len(split_files)):
            seq1 = split_files[i]
            seq2 = split_files[j]
            pair_id = seq1.split(".")[-2] + "_" + seq2.split(".")[-2]
            pair_dir = os.path.join(pairs_rootdir, pair_id)
            seq1f = os.path.join(sequence_dir, seq1)
            seq2f = os.path.join(sequence_dir, seq2)
            alignf = Path(os.path.join(pair_dir, align_out))
            # print(pair_id)
            # print(seq1f)
            # print(seq2f)
            # print(alignf)
            masa(pair_dir, seq1f, seq2f, alignf)
            metrics.append(get_metrics(pair_dir, alignf, seq1, seq2))
    joined_sequences = os.path.join(base_work_dir, "selected_sequences.fasta")
    # compss_wait_on(metrics)
    filter_sequences(metrics, sequence_dir, max_sequences,
                     similar, joined_sequences)
    msa_alignment = os.path.join(base_work_dir, "msa_alignment.fasta")
    result = pastar(msa_alignment, pastar_threads, joined_sequences)
    compss_wait_on(result)
    if os.path.exists(msa_alignment):
        print("Workflow finished with success")
    else:
        print("Workflow failed!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str,
                        required=True, help="Multifasta file used as input")
    parser.add_argument("-w", "--workdir", type=str, required=False, default=os.getcwd(),
                        help="Working directory, where all the outputs will be stored")
    parser.add_argument("-m", "--max_seqs", type=int, required=False, default=5,
                        help="Maximum number of sequences that will undergo the multi-sequence alignment with pastar")
    parser.add_argument("-p", "--pastar_threads", type=int, required=False,
                        default=1, help="Number of threads to be used in pastar")
    parser.add_argument("-s", "--similar", help="By default the workflow will select the most divergent sequences",
                        type=bool, action=argparse.BooleanOptionalAction, default=False)
    args = parser.parse_args()
    main(args.input, args.workdir, args.max_seqs,
         args.pastar_threads, args.similar)
