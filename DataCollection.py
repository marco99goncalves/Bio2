import requests

PATH = "Sequences"
UNIPROT_ID = 'P68871'

def fetch_uniprot_sequence(uniprot_id):
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
        print(fasta_data)  # For demonstration, printing the data
        
        # Write the sequence to a file
        with open(f"{PATH}/sequence_{uniprot_id}.fasta", 'w') as file:
            file.write(fasta_data)
        print(f"Sequence written to {PATH}/{uniprot_id}.fasta")
    else:
        print(f"Failed to retrieve data for UniProt ID {uniprot_id}")

# Example usage:
fetch_uniprot_sequence(UNIPROT_ID)  # Hemoglobin subunit beta (This ID might change)
