from pycompss.api.api import compss_wait_on
from apps import *
import argparse
import os
from pathlib import Path
import csv

def main(input_alignment_file, base_work_dir, max_sequences, pastar_threads, mode):
    align_out = "alignment.00.txt"
    sequence_dir = compss_wait_on(
        create_dir(os.path.join(base_work_dir, "sequences")))
    split_files = compss_wait_on(
        split_sequences(input_alignment_file, sequence_dir))
    pairs_rootdir = compss_wait_on(
        create_dir(os.path.join(base_work_dir, "pair_sequences")))
    pair_ids = []
    folders = []
    for i in range(len(split_files)):
        for j in range(i + 1, len(split_files)):
            pair_id = (
                split_files[i].split(".")[-2]
                + "_"
                + split_files[j].split(".")[-2])
            folders.append(
                create_dir(os.path.join(pairs_rootdir, pair_id)))
            pair_ids.append((split_files[i], split_files[j]))
    compss_wait_on(folders)
    raw_metrics = []
    for seq1, seq2 in pair_ids:
        pair_id = seq1.split(".")[-2] + "_" + seq2.split(".")[-2]
        pair_dir = os.path.join(pairs_rootdir, pair_id)
        seq1f = os.path.join(sequence_dir, seq1)
        seq2f = os.path.join(sequence_dir, seq2)
        alignf = Path(os.path.join(pair_dir, align_out))
        masa(pair_dir, seq1f, seq2f, alignf)
        metrics = get_metrics(pair_dir, alignf, seq1, seq2)
        raw_metrics.append(metrics)
    identities = []
    for m in raw_metrics:
        identities.append(compute_identity(m))
    identities = compss_wait_on(identities)
    pairs = []
    for (seq1, seq2), identity in zip(pair_ids, identities):
        pairs.append((seq1, seq2, identity))
    similarity = {s: {} for s in split_files}
    for s1, s2, identity in pairs:
        similarity[s1][s2] = identity
        similarity[s2][s1] = identity

    if mode == "divergent":
        distance = {
            s: {
                t: 100.0 - similarity[s][t]
                for t in similarity[s]
            }
            for s in similarity
        }
        metric_name = "average pairwise distance (100 - identity)"
    else:
        distance = similarity
        metric_name = "average pairwise identity"

    avg_metric = {}
    for s in split_files:
        if distance[s]:
            avg_metric[s] = sum(distance[s].values()) / len(distance[s])
        else:
            avg_metric[s] = 0.0

    first_selected = max(avg_metric, key=avg_metric.get)

    explanation_txt = os.path.join(
        base_work_dir,
        "first_minmax_selection.txt"
    )

    with open(explanation_txt, "w") as f:
        f.write(f"First sequence selected: {first_selected}\n")
        f.write(f"Mode: {mode}\n")
        f.write(f"Metric used: {metric_name}\n")
        f.write(f"Average metric value: {avg_metric[first_selected]:.6f}\n\n")
        f.write("Pairwise values:\n")
        for other in distance[first_selected]:
            f.write(
                f"{first_selected} vs {other}: "
                f"{distance[first_selected][other]:.6f}\n"
            )

    csv_out = os.path.join(base_work_dir, "pairwise_identity.csv")
    write_pairwise_csv(pairs, csv_out)

    selected = maxmin_selection(
        pairs,
        split_files,
        max_sequences,
        mode=mode)

    joined_sequences = os.path.join(
        base_work_dir,
        "selected_sequences.fasta")

    write_selected_sequences(
        selected,
        sequence_dir,
        joined_sequences)

    msa_alignment = os.path.join(
        base_work_dir,
        "msa_alignment.fasta")

    compss_wait_on(
        pastar(msa_alignment, pastar_threads, joined_sequences))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True, help="Multifasta file used as input")
    parser.add_argument("-w", "--workdir", type=str, required=False, default=os.getcwd(), help="Working directory, where all the outputs will be stored")
    parser.add_argument("-m", "--max_seqs", type=int, required=False, default=5, help="Maximum number of sequences that will undergo the multi-sequence alignment with pastar")
    parser.add_argument("-p", "--pastar_threads", type=int, required=False, default=1, help="Number of threads to be used in pastar")
    parser.add_argument(
        "--mode",
        choices=["divergent", "similar"],
        default="divergent")
    args = parser.parse_args()
    main(
        args.input,
        args.workdir,
        args.max_seqs,
        args.pastar_threads,
        args.mode)
