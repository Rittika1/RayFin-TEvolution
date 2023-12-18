import re
import sys

def read_fasta(fasta_file):
    sequences = []
    current_sequence = None
    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_sequence:
                    sequences.append(current_sequence)
                current_sequence = {'header': line, 'sequence': ''}
            else:
                current_sequence['sequence'] += line
    if current_sequence:
        sequences.append(current_sequence)
    return sequences

def write_fasta(sequences, output_file):
    with open(output_file, 'w') as file:
        for sequence in sequences:
            file.write(sequence['header'] + '\n')
            sequence_lines = [sequence['sequence'][i:i+70] for i in range(0, len(sequence['sequence']), 70)]
            file.write('\n'.join(sequence_lines) + '\n')

def exclude_sequences(contamination_file, fasta_file):
    exclude_seqs = set()
    found_exclude = False
    with open(contamination_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line == "Exclude:":
                found_exclude = True
            elif line == "Trim:":
                found_exclude = False
            elif found_exclude:
                parts = line.split('\t')
                sequence_name = parts[0]
                exclude_seqs.add(sequence_name)
    return exclude_seqs

def trim_sequences(contamination_file, fasta_file):
    trim_info = {}
    found_trim = False
    with open(contamination_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line == "Trim:":
                found_trim = True
            elif line == "Duplicated:":
                found_trim = False
            elif found_trim:
                parts = line.split('\t')
                if len(parts) >= 3:
                    sequence_name = parts[0]
                    spans = parts[2]
                    start, end = map(int, spans.split(',')[0].split('..'))
                    trim_info[sequence_name] = (start, end)
    return trim_info

def keep_duplicated_sequences(contamination_file, fasta_file):
    keep_seqs = set()
    found_duplicate = False
    with open(contamination_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line == "Duplicated:":
                found_duplicate = True
            elif found_duplicate:
                parts = line.split()
                if parts:
                    sequence_name = parts[0]
                    keep_seqs.add(sequence_name)
    return keep_seqs

def filter_sequences_by_length(sequences):
    filtered_sequences = []
    for sequence in sequences:
        if len(sequence['sequence']) >= 200:
            filtered_sequences.append(sequence)
    return filtered_sequences

def modify_sequence_headers(sequences, sample_id):
    modified_sequences = []
    for sequence in sequences:
        header = sequence['header'][1:]  # Remove ">" from the header
        if sample_id not in header:
            if "TRINITY" in header:
                new_header = ">" + header.replace("TRINITY", sample_id)
            else:
                new_header = ">" + sample_id + "_" + header
            sequence['header'] = new_header
        modified_sequences.append(sequence)
    return modified_sequences

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py contamination_file.fasta input_fasta_file.fasta sample_id")
        sys.exit(1)
    
    contamination_file = sys.argv[1]
    fasta_file = sys.argv[2]
    sample_id = sys.argv[3]

    exclude_seqs = exclude_sequences(contamination_file, fasta_file)
    trim_info = trim_sequences(contamination_file, fasta_file)
    keep_seqs = keep_duplicated_sequences(contamination_file, fasta_file)
    
    sequences = read_fasta(fasta_file)

    filtered_sequences = []
    
    for sequence in sequences:
        header = sequence['header'][1:]  # Remove ">" from the header
        if header not in exclude_seqs:
            if header in trim_info:
                start, end = trim_info[header]
                sequence['sequence'] = sequence['sequence'][:start - 1] + sequence['sequence'][end:]
            if header not in keep_seqs:
                filtered_sequences.append(sequence)
    
    # Filter sequences by length
    filtered_sequences = filter_sequences_by_length(filtered_sequences)
    
    # Modify sequence headers
    modified_sequences = modify_sequence_headers(filtered_sequences, sample_id)
    
    # Write the final output to a single file
    output_file = f"{sample_id}_filtered_and_modified_absolutefinal.fasta"
    write_fasta(modified_sequences, output_file)
    
    print(f"Filtered and modified sequences saved to {output_file}")

if __name__ == "__main__":
    main()
