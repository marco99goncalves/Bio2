from Bio import ExPASy
from Bio import SeqIO
from Bio import AlignIO
import requests
import time
import Logger

logger = Logger.setup_logger()

def align_sequences(input_file, output_file):
    # The URL to the EBI Clustal Omega REST service
    url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"
    headers = {'Content-Type': 'application/x-www-form-urlencoded'}

    with open(input_file, 'r') as file:
        content = file.read().strip()


    # Prepare the payload. Join sequences into a single string for submission.
    # The sequences parameter should be a list of sequence strings in FASTA format.
    payload = {
        'email': 'marco99goncalves@gmail.com',  # Use your actual email here
        'sequence': content
    }

    # Submit the job
    response = requests.post(url + "/run", headers=headers, data=payload)
    if not response.ok:
        logger.error("Error submitting job")
        logger.error(response.text)
        return
    


    job_id = response.text
    logger.info(f"Job submitted successfully. Job ID: {job_id}")

    waitTime = 1
    timesWaited = 0
    # Check job status
    while True:
        status_response = requests.get(url + f"/status/{job_id}")
        status = status_response.text
        logger.info(f"Job status: {status}")
        if status == 'FINISHED':
            break
        timesWaited += 1
        if(timesWaited == 5): #Serves way to prevent too much logging
            waitTime *= 2
            timesWaited = 0
        time.sleep(waitTime)

    # Retrieve the results
    result_url = url + f"/result/{job_id}/fa"
    result_response = requests.get(result_url)
    if result_response.ok:
        with open(output_file, 'w') as file:
            file.write(result_response.text)
        logger.info(f"Alignment saved to {output_file}")
    else:
        logger.error("Error retrieving results")
        logger.error(result_response.text)
