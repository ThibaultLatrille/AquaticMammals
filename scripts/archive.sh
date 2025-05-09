#!/bin/bash
# Find and archive cov and nhx files from Bayescode in the cluster
find . -type f \( -name "*.cov" -o -name "*.nhx" \) -print0 | tar --null -czf cov_archive.tar.gz --files-from=-

find . -type f \( -name "*.RELAX.json" \) -print0 | tar --null -czf relax_archive.tar.gz --files-from=-


# Change the working directory for all .param files
old_string="/work/FAC/FBM/DBC/nsalamin/software/nsalamin/default/tlatrill"
new_string="/Users/tlatrille/Documents"
# Find all .param files and replace the string
find "." -type f -name "*.param" | while read -r file; do
    sed -i.bak "s|$old_string|$new_string|g" "$file"
done

# Compress hyphy preprocessed files local
tar -cvzf "data_processed.tar.gz" ./data_processed/

# Uncompress the tar file in remote
tar -xvzf "data_processed.tar.gz"
rm -rf ./data_processed/*/._*
rm -rf ./data_processed/*/*/._*
# touch ./data_processed/*/*/*
