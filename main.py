from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np

start_codon = "ATG"
stop_codons = ["TAG", "TAA", "TGA"]
complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
amino_acids = "ACDEFGHIKLMNPQRSTVWY"

fasta_files = [
    "./data/bacterial1.fasta",
    "./data/bacterial2.fasta",
    "./data/bacterial3.fasta",
    "./data/bacterial4.fasta",
    "./data/mamalian1.fasta",
    "./data/mamalian2.fasta",
    "./data/mamalian3.fasta",
    "./data/mamalian4.fasta"
]


def reverse_complement(seq):
  reverse_compliment_seq = ''.join(complement[nucleotide] for nucleotide in seq.upper())
  return reverse_compliment_seq[::-1]


def find_dna_pairs(seq):
  orf_positions = []  #(start_index, end_index)
  seq_length = len(seq)

  # Loop over all 3 possible reading frames
  for start_pos in range(3):
    current_pos = start_pos

    # Walk through the sequence in steps of 3 (codon length)
    while current_pos + 3 <= seq_length:
      codon = seq[current_pos:current_pos+3]
      if codon == start_codon:
        start_index = current_pos
        current_pos += 3

        # Continue reading codons until a stop codon is found
        while current_pos + 3 <= seq_length:
          codon = seq[current_pos:current_pos+3]
          current_pos += 3

          if codon in stop_codons:
            end_index = current_pos
            orf_positions.append((start_index, end_index))
            break
      else:
        current_pos += 3

  return orf_positions

def print_genome_info(genome):
  name, pairs, sequence, rev_pairs, rev_seq = genome
  count = 1;

  print(f"Genome: {name}")
  for (start, finish) in pairs:
    orf_sequence = sequence[start:finish]
    print(f"{count}) {orf_sequence}")
    count += 1

  for (start, finish) in rev_pairs:
    orf_sequence = rev_seq[start:finish]
    print(f"{count}) {orf_sequence}")
    count += 1

def get_all_genome_dna_info():
  genomes = []
  for file_name in fasta_files:
    seq_record = SeqIO.read(file_name, "fasta")

    sequence = str(seq_record.seq)
    reverse_complement_sequence = reverse_complement(sequence)
    genome_name = seq_record.id

    forward_pairs = find_dna_pairs(sequence)
    filtered_forward_pairs = [(start, end) for start, end in forward_pairs if (end - start) >= 100]

    reverse_pairs = find_dna_pairs(reverse_complement_sequence)
    filtered_reverse_pairs = [(start, end) for start, end in reverse_pairs if (end - start) >= 100]

    genomes.append((genome_name, filtered_forward_pairs, sequence, filtered_reverse_pairs, reverse_complement_sequence))

  return genomes

def dna_to_protein_genomes_info(genomes_info):

  translated_info = []
  for genome_info in genomes_info:
    name, pairs, sequence, rev_pairs, rev_seq = genome_info
    translated_sequence = [Seq(sequence[start:stop]).translate(table="Bacterial", to_stop=True) for (start, stop) in pairs] + [Seq(rev_seq[start:stop]).translate(table="Bacterial", to_stop=True) for (start, stop) in rev_pairs]
    translated_info.append((name, translated_sequence))

  return translated_info

def count_codon_frequencies(protein_seqs):
  codon_hash = dict.fromkeys(amino_acids, 0)
  total_codons = 0

  for seq in protein_seqs:
    for codon in seq:
      if codon in codon_hash:
        codon_hash[codon] += 1
        total_codons += 1

  if total_codons > 0:
        codon_hash = {codon: count / total_codons for codon, count in codon_hash.items()}

  return codon_hash

def count_dicodon_frequencies(protein_seqs):
  dicodon_hash = dict.fromkeys([f"{a}{b}" for a in amino_acids for b in amino_acids], 0)
  total_dicodons = 0

  for seq in protein_seqs:
    for i in range(len(seq) - 1):
      dicodon = seq[i:i + 2]
      if dicodon in dicodon_hash:
        dicodon_hash[dicodon] += 1
        total_dicodons += 1

  if total_dicodons > 0:
        dicodon_hash = {dicodon: count / total_dicodons for dicodon, count in dicodon_hash.items()}

  return dicodon_hash

def get_all_genome_frequencies(all_genome_protein_info):
  genomes_protein_frequency = []
  for name, translated_sequence in all_genome_protein_info:
    genomes_protein_frequency.append((name, count_codon_frequencies(translated_sequence), count_dicodon_frequencies(translated_sequence)))

  return genomes_protein_frequency

#for name, codons_hash, dicodons_hash in genome_frequencies:
  #print(name)
  #for codon, frequency in codons_hash.items():
    #print(f"{codon} {frequency}")

def compute_distance_matrix(frequencies_list):
  n = len(frequencies_list)
  distance_matrix = np.zeros((n, n))

  for i in range(n):
    for j in range(i, n):
      freq1 = list(frequencies_list[i].values())
      freq2 = list(frequencies_list[j].values())
      dist = np.sqrt(sum((freq1[k] - freq2[k])**2 for k in range(len(freq1))))
      distance_matrix[i][j] = dist
      distance_matrix[j][i] = dist

  return distance_matrix

def print_distance_matrix(matrix, names):
  n = len(names)
  print(n)

  for i in range(n):
    print(names[i], end="")
    for j in range(n):
      print(f"{matrix[i][j]:7.3f}", end="")
    print()

def write_distance_matrix_to_phylip(matrix, names, file_name):
  n = len(names)

  with open(file_name, 'w') as f:
    f.write(f"{n}\n")

    for i in range(n):
      f.write(f"{names[i]}")
      f.write(" ".join(f"{matrix[i][j]:7.3f}" for j in range(n)))
      f.write("\n")

all_genome_info = get_all_genome_dna_info()
translated_genome_info = dna_to_protein_genomes_info(all_genome_info)
genome_frequencies = get_all_genome_frequencies(translated_genome_info)

codon_freqency_hash = [codon_hash for _, codon_hash, _ in genome_frequencies]
dicodon_frequency_hash = [dicodon_hash for _, _, dicodon_hash in genome_frequencies]
genome_names = [names for names, _, _ in genome_frequencies]

codon_distance_matrix = compute_distance_matrix(codon_freqency_hash)
dicodon_distance_matrix = compute_distance_matrix(dicodon_frequency_hash)

# Print Codon Distance Matrix
print("Codon Distance Matrix:")
print_distance_matrix(codon_distance_matrix, genome_names)

# Print Dicodon Distance Matrix
print("\nDicodon Distance Matrix:")
print_distance_matrix(dicodon_distance_matrix, genome_names)

write_distance_matrix_to_phylip(codon_distance_matrix, genome_names, "outputs/codon_distance_matrix.phy")
write_distance_matrix_to_phylip(dicodon_distance_matrix, genome_names, "outputs/dicodon_distance_matrix.phy")