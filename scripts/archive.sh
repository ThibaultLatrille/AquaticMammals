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
tar -cvzf "hyphy.tar.gz" ./data_processed/Hyphy/

# Uncompress the tar file in remote
tar -xvzf "hyphy.tar.gz"
rm -rf ./data_processed/Hyphy/._*
rm -rf ./data_processed/Hyphy/*/._*
touch ./data_processed/Hyphy/*/*

/Users/tlatrille/Documents/AquaticMammals/utils/BayesCode/bin/nodeomega /Users/tlatrille/Documents/AquaticMammals/data_processed/OrthoMam/23345/nodeomega_1
/Users/tlatrille/Documents/AquaticMammals/utils/BayesCode/bin/readnodeomega --every 1 --until 2000 --burnin 1000 --cov /Users/tlatrille/Documents/AquaticMammals/data_processed/OrthoMam/23345/nodeomega_1
/Users/tlatrille/Documents/AquaticMammals/utils/BayesCode/bin/readnodeomega --every 1 --until 2000 --burnin 1000 --newick /Users/tlatrille/Documents/AquaticMammals/data_processed/OrthoMam/23345/nodeomega_1