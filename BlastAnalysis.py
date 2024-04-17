from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

UNIPROT_ID = 'P68871'

def blast_sequence_and_get_top_hits(fasta_path, output_path, top_hits=10):
    # Read the fasta sequence
    print("0\n")

    record = SeqIO.read(fasta_path, format="fasta")
    print("11\n")
    
    # Perform BLAST search
    result_handle = NCBIWWW.qblast("blastp", "nr", record.seq)
    print("2\n")
    
    # Parse BLAST results
    blast_records = NCBIXML.parse(result_handle)
    print("3\n")
    # Extract top hits for different species
    hits = {}
    print(blast_records)
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if len(hits) >= top_hits:
                    break
                species_name = alignment.hit_def.split(">")[0].split("[")[-1].replace("]", "")
                #
                # if species_name not in hits and species_name != "Homo sapiens":
                if species_name != "Homo sapiens":
                    hits[species_name] = hsp.sbjct
    print("4\n")

    # Write top hits to a file
    with open(output_path, "w") as output_file:
        for species, sequence in hits.items():
            output_file.write(f">{species}\n{sequence}\n")
    print("3\n")

    return output_path

# The fasta file should already contain the sequence you want to blast
fasta_path = f"Sequences/sequence_{UNIPROT_ID}.fasta"
output_path = f"Phase2/sequences_{UNIPROT_ID}.fasta"

# Run the function
output_file = blast_sequence_and_get_top_hits(fasta_path, output_path)
print(f"Top hits written to {output_file}")
