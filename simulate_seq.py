limport random
from Bio import SeqIO

def simulate_rna_seq_reads(input_fasta, output_fasta, read_length, num_reads):
    # Load the cDNA sequence from the input FASTA file
    with open(input_fasta, "r") as input_handle:
        cDNA_record = next(SeqIO.parse(input_handle, "fasta"))
        cDNA_sequence = str(cDNA_record.seq)

    sequence_length = len(cDNA_sequence)
    reads = []

    # Generate reads
    for i in range(num_reads):
        # Randomly select a start position for the read
        start_pos = random.randint(0, sequence_length - read_length)
        read_seq = cDNA_sequence[start_pos:start_pos + read_length]
        read_id = f"read_{i+1}"
        reads.append((read_id, read_seq))

    # Write the reads to an output FASTA file
    with open(output_fasta, "w") as output_handle:
        for read_id, read_seq in reads:
            output_handle.write(f">{read_id}\n{read_seq}\n")

    print(f"Generated {num_reads} reads of length {read_length} from the input cDNA sequence.")
    print(f"Output saved to {output_fasta}")

# Example usage
input_fasta = "path/to/fasta/input_cdna.fasta"   # Path to the input cDNA FASTA file
output_fasta = "path/to/fasta/simulated_reads_1000.fasta"  # Path to the output FASTA file for simulated reads
read_length = 100  # Length of each read
num_reads = 10000   # Number of reads to generate

simulate_rna_seq_reads(input_fasta, output_fasta, read_length, num_reads)
