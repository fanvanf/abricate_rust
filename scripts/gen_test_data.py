#!/usr/bin/env python3
"""Generate test AMR database sequences and a test genome."""
import random
import os

random.seed(42)

def random_dna(length):
    return ''.join(random.choice('ATGC') for _ in range(length))

# Generate 3 AMR gene sequences (~500bp each)
gene1 = random_dna(500)  # blaTEM-1
gene2 = random_dna(450)  # aac3-Ia
gene3 = random_dna(400)  # tetA

# Write database sequences file
db_dir = os.path.join(os.path.dirname(__file__), '..', 'db', 'test_amr')
os.makedirs(db_dir, exist_ok=True)

with open(os.path.join(db_dir, 'sequences'), 'w') as f:
    f.write(f">test_amr~~~blaTEM-1~~~TEST001~~~BETA-LACTAM class A beta-lactamase TEM-1\n")
    for i in range(0, len(gene1), 70):
        f.write(gene1[i:i+70] + '\n')
    
    f.write(f">test_amr~~~aac3-Ia~~~TEST002~~~GENTAMICIN aminoglycoside N-acetyltransferase AAC(3)-Ia\n")
    for i in range(0, len(gene2), 70):
        f.write(gene2[i:i+70] + '\n')
    
    f.write(f">test_amr~~~tetA~~~TEST003~~~TETRACYCLINE tetracycline efflux protein TetA\n")
    for i in range(0, len(gene3), 70):
        f.write(gene3[i:i+70] + '\n')

# Create test genome: contains gene1 and gene2 exactly, plus some random flanking
genome_parts = [
    random_dna(200),  # flanking
    gene1,             # blaTEM-1
    random_dna(300),  # flanking
    gene2,             # aac3-Ia
    random_dna(150),  # flanking
]

test_dir = os.path.join(os.path.dirname(__file__), '..', 'test')
os.makedirs(test_dir, exist_ok=True)

with open(os.path.join(test_dir, 'test_genome.fasta'), 'w') as f:
    f.write(">contig_1\n")
    genome_seq = ''.join(genome_parts)
    for i in range(0, len(genome_seq), 70):
        f.write(genome_seq[i:i+70] + '\n')

print(f"Database sequences written to {db_dir}/sequences")
print(f"Test genome written to {test_dir}/test_genome.fasta")
print(f"Gene1 (blaTEM-1): {len(gene1)}bp")
print(f"Gene2 (aac3-Ia): {len(gene2)}bp")
print(f"Gene3 (tetA): {len(gene3)}bp")
print(f"Genome contig_1: {len(genome_seq)}bp (contains blaTEM-1 and aac3-Ia)")
