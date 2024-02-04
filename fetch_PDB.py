import os
import datetime
import requests

# Create a directory with today's date for storing downloaded PDB files
directory = "pdbDataset_" + datetime.date.today().strftime('%d%m%Y')

# Check if the directory already exists, and create it if not
if not os.path.exists(directory):
    os.mkdir(directory)

def fetch(id):
    """
    Fetch PDB file content for a given PDB ID.

    Args:
        id (str): PDB ID to fetch.

    Returns:
        bytes: Content of the PDB file as bytes.
    """
    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % id
    try:
        response = requests.get(url, verify=True)
        response.raise_for_status()  # Raise an error for bad responses
        return response.content
    except requests.exceptions.RequestException as e:
        print(f"Error fetching PDB file for ID {id}: {e}")
        return None

# Specify the filename containing a list of PDB IDs
pdb_list_filename = 'pdbList(11April2017).txt'

# Open the file containing PDB IDs and iterate through each ID
with open(pdb_list_filename) as pdb_list:
    for pdb_id in pdb_list:
        pdb_id = pdb_id.strip()
        page = fetch(pdb_id)

        # If PDB file content is successfully fetched, save it to a file
        if page is not None:
            output_file_path = os.path.join(directory, f"{pdb_id}.pdb")
            with open(output_file_path, "wb") as output_file:
                output_file.write(page)
            print(f"Downloaded {pdb_id}.pdb successfully.")
