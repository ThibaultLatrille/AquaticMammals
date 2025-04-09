#!/bin/bash
# Find and archive the files
tar -czf "cov.tar.gz" ./data_processed/OrthoMam/*/*.cov ./data_processed/OrthoMam/*/*.nhx

find . -type f \( -name "*.cov" -o -name "*.nhx" \) -print0 | tar --null -czf cov_archive.tar.gz --files-from=-

# Define the directory to search and the strings to replace
base_dir="."
old_string="/work/FAC/FBM/DBC/nsalamin/software/nsalamin/default/tlatrill"
new_string="/Users/tlatrille/Documents"

# Find all .param files and replace the string
find "$base_dir" -type f -name "*.param" | while read -r file; do
    sed -i.bak "s|$old_string|$new_string|g" "$file"
done

echo "Replacement completed."