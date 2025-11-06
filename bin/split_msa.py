#!/usr/bin/env python3
# Modified from ProteinFold project [https://github.com/nf-core/proteinfold/blob/dev/bin/msa_manager.py](https://github.com/nf-core/proteinfold/blob/dev/bin/msa_manager.py)


import os
import string
import argparse
from pathlib import Path



MAX_MSA_SEQS = 16384
MAX_PAIRED_SEQS = 8192
ID_CHARS = list(string.ascii_uppercase) + list(string.ascii_lowercase) + [str(x) for x in range(10)]



def get_sub_sequences(seq_lengths, whole_seq):
    out_seqs = []
    curr_seq = ""
    curr_seq_itr = 0
    total_letters = 0
    for letter in whole_seq:
        curr_seq += letter
        if letter.isupper() or letter == "-":
            total_letters += 1
        if total_letters == seq_lengths[curr_seq_itr]:
            out_seqs.append(curr_seq)
            curr_seq = ""
            curr_seq_itr += 1
            total_letters = 0


    if len(out_seqs) != len(seq_lengths):
        print("Something wrong in the input file, could not generate the required number of sequences")
        exit(1)


    return out_seqs



def parse_msa(msa_path, output_dir, base_filename):
    os.makedirs(output_dir, exist_ok=True)
    homolog = ""
    section_index = 0


    with open(msa_path, "r") as file:
        first_line = file.readline()
        if not first_line.startswith("#"):
            homologs_lengths = [len(file.readline().strip('\n'))]
            sequence_groups = [[[],[]]]
            is_multimer = False
        else:
            homologs_lengths = [int(x.strip()) for x in first_line.replace('#',"").split()[0].split(",")]
            sequence_groups = [[[], []] for _ in range(len(homologs_lengths))]
            is_multimer = True


    with open(msa_path, "r") as file:
        if is_multimer:
            file.readline()
        header_line = file.readline().strip()[1:]
        expected_section_headers = [x.strip() for x in header_line.split()]
        current_header = header_line
        first_seq = False
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if homolog:
                    if first_seq and section_index > 0:
                        first_seq = False
                    else:
                        sub_sequences = get_sub_sequences(homologs_lengths, homolog)
                        for seq_index in range(len(homologs_lengths)):
                            if section_index == 0:
                                if len(sequence_groups[seq_index][0]) < MAX_PAIRED_SEQS:
                                    sequence_groups[seq_index][0].append(sub_sequences[seq_index])
                            else:
                                if seq_index == section_index - 1:
                                    if len(sequence_groups[seq_index][1]) + len(sequence_groups[seq_index][0]) < MAX_MSA_SEQS:
                                        sequence_groups[seq_index][1].append(sub_sequences[seq_index])


                homolog = ""
                current_header = line[1:].strip()


                if section_index < len(homologs_lengths) and current_header == expected_section_headers[section_index]:
                    section_index += 1
                    first_seq = True
            else:
                homolog += line
        if homolog:
            if first_seq and section_index > 0:
                first_seq = False
            else:
                sub_sequences = get_sub_sequences(homologs_lengths, homolog)
                for seq_index in range(len(homologs_lengths)):
                    if section_index == 0:
                        if len(sequence_groups[seq_index][0]) < MAX_PAIRED_SEQS:
                            sequence_groups[seq_index][0].append(sub_sequences[seq_index])
                    else:
                        if seq_index == section_index - 1:
                            if len(sequence_groups[seq_index][1]) + len(sequence_groups[seq_index][0]) < MAX_MSA_SEQS:
                                sequence_groups[seq_index][1].append(sub_sequences[seq_index])


    for seq_index in range(len(homologs_lengths)):
        filename = os.path.join(output_dir, f"{base_filename}_msa_{seq_index}.csv")
        with open(filename, "w") as out_file:
            out_file.write("key,sequence\n")
            if len(homologs_lengths)==1: #Homo-oligomer: all sequences are paired
                paired_sequences = sequence_groups[seq_index][0]+sequence_groups[seq_index][1]
                for i, seq in enumerate(paired_sequences):
                    out_file.write(f"{i},{seq}\n")
            else:
                paired_sequences = sequence_groups[seq_index][0]
                for i, seq in enumerate(paired_sequences, start=1):
                    out_file.write(f"{i},{seq}\n")


                unpaired_sequences = sequence_groups[seq_index][1]
                for seq in unpaired_sequences:
                    out_file.write(f"-1,{seq}\n")



def main():
    parser = argparse.ArgumentParser(description="Split multi-A3M files in a directory into CSV sequences per section.")
    parser.add_argument("-i","--input_dir", help="Directory containing input .a3m files")
    parser.add_argument("-o", "--output_dir", default="output_msa", help="Directory to write output CSVs")


    args = parser.parse_args()
    
    # Get all .a3m files in the input directory
    input_path = Path(args.input_dir)
    if not input_path.is_dir():
        print(f"Error: {args.input_dir} is not a valid directory")
        exit(1)
    
    a3m_files = list(input_path.glob("*.a3m"))
    
    if not a3m_files:
        print(f"No .a3m files found in {args.input_dir}")
        exit(1)
    
    print(f"Found {len(a3m_files)} .a3m file(s) to process")
    
    # Process each .a3m file
    for a3m_file in a3m_files:
        base_filename = a3m_file.stem  # Gets filename without extension
        print(f"Processing {a3m_file.name}...")
        parse_msa(str(a3m_file), args.output_dir, base_filename)
    
    print(f"Done! Output written to {args.output_dir}")



if __name__ == "__main__":
    main()
