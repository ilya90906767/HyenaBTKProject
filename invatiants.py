import argparse
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import numpy as np
import os

def get_arguments():
    parser = argparse.ArgumentParser(description="Filter invariant sites and delete windows with >50% gaps.")
    parser.add_argument("-i", "--input", type=str, required=True, help="Input folder containing multi-sequence alignment files in FASTA format")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output folder for filtered alignment files in PHYLIP format")
    parser.add_argument("-w", "--window_size", type=int, default=100000, help="Window size for gap checking")
    parser.add_argument("-g", "--gap_threshold", type=float, default=0.5, help="Threshold for gap percentage to delete a window")
    return parser.parse_args()

def filter_invariants(alignment):
    if alignment.get_alignment_length() == 0:
        return alignment

    # Convert alignment to NumPy array with 'S1' type
    align_array = np.array([list(rec) for rec in alignment], dtype='S1')
    
    # Find invariant columns
    invariant_mask = np.all(align_array == align_array[:, 0, None], axis=0)
    
    # Select variant columns
    variant_columns = align_array[:, ~invariant_mask]
    
    # Create filtered alignment
    filtered_records = []
    for idx, record in enumerate(alignment):
        seq = ''.join(variant_columns[idx, :].astype(str))
        filtered_records.append(SeqRecord(seq=Seq(seq), id=record.id, description=record.description))
    
    return MultipleSeqAlignment(filtered_records)

def filter_gaps(alignment, window_size, gap_threshold):
    align_length = alignment.get_alignment_length()
    records = list(alignment)
    
    filtered_sequences = {record.id: [] for record in records}
    
    for start in range(0, align_length, window_size):
        end = min(start + window_size, align_length)
        window = alignment[:, start:end]
        
        gap_counts = np.array([rec.seq.count('-') for rec in window])
        gap_percentage = gap_counts / (end - start)
        
        if np.any(gap_percentage < gap_threshold):
            for rec in window:
                filtered_sequences[rec.id].append(str(rec.seq))
    
    final_records = []
    for record in records:
        filtered_seq = ''.join(filtered_sequences[record.id])
        final_records.append(SeqRecord(Seq(filtered_seq), id=record.id, description=record.description))
    
    return MultipleSeqAlignment(final_records)

def write_alignment(alignment, file):
    with open(file, "w") as output_handle:
        AlignIO.write(alignment, output_handle, "phylip-relaxed")

def read_alignment(file):
    try:
        alignment = SeqIO.parse(file, "fasta")
        return MultipleSeqAlignment(alignment)
    except ValueError as e:
        print(f"Error reading file {file}: {e}")
        return None

def convert_fa_to_phy(fa_file, phy_file):
    with open(fa_file, "r") as fa_handle, open(phy_file, "w") as phy_handle:
        alignment = SeqIO.parse(fa_handle, "fasta")
        AlignIO.write(alignment, phy_handle, "phylip-relaxed")

def record_to_binary(record):
    seq = record.seq
    binary_seq = ''.join(['0' if (seq[i] == 'A' and seq[i+1] == 'G') or (seq[i] == 'G' and seq[i+1] == 'A') or
                          (seq[i] == 'C' and seq[i+1] == 'T') or (seq[i] == 'T' and seq[i+1] == 'C')
                          else '1' for i in range(len(seq)-1)])
    return SeqRecord(seq=Seq(binary_seq), id=record.id, description=record.description)

def make_names_unique(alignment):
    names = [rec.id for rec in alignment]
    unique_names = []
    for name in names:
        original_name = name
        if name in unique_names:
            i = 1
            while True:
                new_name = f"{original_name}_{i}"
                if new_name not in unique_names:
                    name = new_name
                    break
                i += 1
        unique_names.append(name)
    for i, rec in enumerate(alignment):
        rec.id = unique_names[i]
    return alignment

def main():
    args = get_arguments()
    input_folder = args.input
    output_folder = args.output

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for file in os.listdir(input_folder):
        if file.endswith(".fa"):
            print(file)
            input_file = os.path.join(input_folder, file)
            output_file = os.path.join(output_folder, file.replace('.fa', '.phy'))
            name_map_file = os.path.join(output_folder, file.replace('.fa', '_name_map.txt'))

            alignment = read_alignment(input_file)

            if alignment is None or alignment.get_alignment_length() == 0:
                print(f"Skipping empty file: {file}")
                continue

            if len(alignment) < 8:
                print(f"Skipping file with less than 8 records: {file}")
                continue

            print(f"Processing file {file}...")
            print(f"Original alignment length: {alignment.get_alignment_length()}")

            # Filter invariant sites
            filtered_alignment = filter_invariants(alignment)

            print(f"Alignment length after filtering invariants: {filtered_alignment.get_alignment_length()}")

            # Filter windows with >50% gaps
            final_alignment = filter_gaps(filtered_alignment, args.window_size, args.gap_threshold)

            print(f"Alignment length after filtering gaps: {final_alignment.get_alignment_length()}")

            # Rename sequence names to numerical IDs
            name_map = {}
            for i, rec in enumerate(final_alignment):
                original_name = rec.id
                new_name = str(i + 1)
                rec.id = new_name
                name_map[new_name] = original_name

            # Write the name map to a text file
            with open(name_map_file, "w") as f:
                for new_name, original_name in name_map.items():
                    f.write(f"{new_name}\t{original_name}\n")

            # Convert alignment to binary
            binary_records = [record_to_binary(rec) for rec in final_alignment]

            # Write the filtered binary alignment to the output file
            with open(output_file, "w") as f:
                SeqIO.write(binary_records, f, "phylip")

            print(f"Output written to {output_file}\n")

if __name__ == "__main__":
    main()