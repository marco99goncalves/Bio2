import DataCollection
import BlastAnalysis
import Logger
import sys
import MSA
import CreatePhylogeneticTree

logger = Logger.setup_logger()

#get the arugment from the command line
UNIPROT_ID = sys.argv[1]
FASTA_PATH = "sequence.fasta"
PHASE2_OUTPUT_PATH = "sequences_to_analyse.fasta"
PHASE3_OUTPUT_PATH = "alignment.txt"
PHYLO_FILE_NAME = "PhylogeneticTree"


RUN_BLASTA = True
RUN_MSA = True

MATRIX="BLOSUM62"

logger.info(f"Collecting the sequence for UniProt ID {UNIPROT_ID}... \n")
DataCollection.fetch_uniprot_sequence(uniprot_id=UNIPROT_ID, fasta_path=FASTA_PATH)

if RUN_BLASTA:
    logger.info(f"Performing BLAST analysis for UniProt ID {UNIPROT_ID}... (This may take some time)\n")
    BlastAnalysis.blast_sequence_and_get_top_hits(fasta_path=FASTA_PATH, output_path=PHASE2_OUTPUT_PATH, matrix_name=MATRIX)
    logger.info(f"BLAST analysis completed. Results written to {PHASE2_OUTPUT_PATH}\n")


    # Append the results from phase 1 to the file from blasta analysis
    logger.info("Appending the results from phase 1 to the BLAST results...\n")
    with open(PHASE2_OUTPUT_PATH, 'a') as file:
        with open(FASTA_PATH, 'r') as phase1_file:
            file.write(phase1_file.read())
else:
    logger.info("Skipping BLAST analysis step.\n")

# Run the MSA
if RUN_MSA:
    logger.info("Running the multiple sequence alignment (MSA)...\n")
    MSA.align_sequences(PHASE2_OUTPUT_PATH, PHASE3_OUTPUT_PATH)
    logger.info(f"MSA completed. Results written to {PHASE3_OUTPUT_PATH}\n")
else:
    logger.info("Skipping MSA step.\n")

# Build the phylogenetic tree
logger.info("Building the phylogenetic tree...\n")
<<<<<<< HEAD
CreatePhylogeneticTree.build(PHASE3_OUTPUT_PATH, matrix_name=MATRIX)
=======
CreatePhylogeneticTree.build(PHASE3_OUTPUT_PATH, PHYLO_FILE_NAME, UNIPROT_ID, export_formats=["pdf"])
>>>>>>> refs/remotes/origin/main
logger.info("Phylogenetic tree built successfully.\n")