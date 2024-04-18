import DataCollection
import BlastAnalysis
import Logger
import sys
import MSA

logger = Logger.setup_logger()

#get the arugment from the command line
UNIPROT_ID = sys.argv[1]
FASTA_PATH = "sequence.fasta"
PHASE2_OUTPUT_PATH = "sequences_to_analyse.fasta"
PHASE3_OUTPUT_PATH = "alignment.txt"

logger.info(f"Collecting the sequence for UniProt ID {UNIPROT_ID}... \n")
DataCollection.fetch_uniprot_sequence(uniprot_id=UNIPROT_ID, fasta_path=FASTA_PATH)

logger.info(f"Performing BLAST analysis for UniProt ID {UNIPROT_ID}... (This may take some time)\n")
#output_file = BlastAnalysis.blast_sequence_and_get_top_hits(fasta_path=FASTA_PATH, output_path=PHASE2_OUTPUT_PATH)
output_file = PHASE2_OUTPUT_PATH
logger.info(f"BLAST analysis completed. Results written to {output_file}\n")

# Run the MSA
MSA.align_sequences(PHASE2_OUTPUT_PATH, PHASE3_OUTPUT_PATH)