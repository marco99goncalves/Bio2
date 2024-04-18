import requests

import Logger
logger = Logger.setup_logger()

def fetch_uniprot_sequence(uniprot_id, fasta_path):
    """
    Fetches the sequence and additional information for a given UniProt ID.
    Args:
    - uniprot_id: The identifier of the UniProt entry.
    
    Returns:
    - Writes the sequence to a fasta file and prints additional information.
    """
    base_url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = requests.get(base_url)
    
    if response.ok:
        fasta_data = response.text

        # Write the sequence to a file
        with open(fasta_path, 'w') as file:
            file.write(fasta_data)
        logger.info(f"Sequence written to {fasta_path}\n")
    else:
        logger.error(f"Failed to retrieve data for UniProt ID {uniprot_id}\n")