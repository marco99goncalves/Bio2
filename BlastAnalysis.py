from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
import Logger

logger = Logger.setup_logger()

def blast_sequence_and_get_top_hits(fasta_path, output_path, top_hits=11):
    # Read the fasta sequence

    record = SeqIO.read(fasta_path, format="fasta")
    
    # Perform BLAST search

    #result_handle = NCBIWWW.qblast("blastp", "nr", record.seq, hitlist_size=1000)

    result_handle = NCBIWWW.qblast( program="blastp", database="nr", sequence=record.seq,
                                    expect=0.05, word_size=5,
                                    matrix_name="BLOSUM62", gapcosts="11 1",
                                    hitlist_size=1000)

    logger.info("BLAST search completed. Parsing results...\n")

    blast_records = NCBIXML.parse(result_handle)
    logger.info("BLAST search parsed.\n")

    # Extract top hits for different species
    hits = {}
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                species_name = alignment.hit_def.split(">")[0].split("[")[-1].replace("]", "").strip().lower()
                if len(hits) >= top_hits:
                    break

                if species_name not in hits and species_name.find("synthetic") == -1 and species_name.find("homo") == -1:
                    hits[species_name.replace(" ", "_")] = (hsp.sbjct, hsp.score)

    # Write top hits to a file
    with open(output_path, "w") as output_file:
        for species, (sequence, score) in hits.items():
            output_file.write(f">{species}\n{sequence}\n")

    logger.info("Sequences written to {output_path}\n")

    return output_path