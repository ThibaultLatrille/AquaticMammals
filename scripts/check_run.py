import os
import tarfile

# Define the target folder
base_folder = 'data_processed/OrthoMam'
output_archive = 'filtered_folders.tar.gz'

# List to store folders that meet the criteria
folders_to_archive = []

# Iterate through all subfolders
for root, dirs, files in os.walk(base_folder):
    # Check if the required file exists and the excluded file does not exist
    if 'placnr.rootree' in files and 'nodeomega_1.cov' not in files:
        folders_to_archive.append(root)

print(f"Folders to archive: {folders_to_archive}")

# Create a tar.gz archive for the selected folders
with tarfile.open(output_archive, 'w:gz') as archive:
    for folder in folders_to_archive:
        archive.add(folder, arcname=os.path.relpath(folder, base_folder))

print(f"Archive created: {output_archive}")