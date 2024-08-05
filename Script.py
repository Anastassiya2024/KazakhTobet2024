import pandas as pd
from Bio import SeqIO

# Load BLAST results
blast_results = pd.read_csv('blast_results.txt', sep='\t', header=None)
blast_results.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

# Filter for mitochondrial sequences (assuming identity less than 95% with nuclear genome indicates true mitochondrial origin)
mito_seqs = set(blast_results[blast_results['pident'] < 95]['qseqid'])

# Load original FASTA and filter sequences
with open('filtered_mitochondrial.fasta', 'w') as output_handle:
    for record in SeqIO.parse("mitochondrial.fasta", "fasta"):
        if record.id in mito_seqs:
            SeqIO.write(record, output_handle, "fasta")

# Print summary of filtering
print(f"Total sequences in original FASTA: {len(list(SeqIO.parse('mitochondrial.fasta', 'fasta')))}")
print(f"Total sequences retained after filtering: {len(mito_seqs)}")
