
#!/bin/bash

echo "Enterotool conda Environment Setup:"
echo "Hi, this is the installer, If you proceed, it will create a conda enviorment with all the required packages!!"
echo "If you have all the depencies requiered, you dont need to run this script."
echo "Make sure you have the following tools installed: seqkit, trimmomatic, vsearch, blastn, matplotlib, pandas"
echo "if not, just run this script and it will install all the required packages."
echo "Do you want to proceed? (y/n)"

read -r response 
if [[ ! "$response" =~ ^([yY][eE][sS]|[yY])$ ]]; then
    echo "Exiting..."
    exit 0
fi

#in case you don like the name a came up wt.
ENV_NAME="enterotool"

if ! command -v conda &> /dev/null; then
    echo "Conda is not installed. Please install Miniconda or Anaconda and try again." >&2
    exit 1
fi

echo "Creating conda environment: $ENV_NAME"
conda create -y -n $ENV_NAME python=3.9

echo "Activating environment: $ENV_NAME"
source $(conda info --base)/etc/profile.d/conda.sh
conda activate $ENV_NAME

echo "Installing required packages..."
conda install -y -c bioconda seqkit trimmomatic vsearch blast
conda install -y matplotlib pandas

echo "Verifying installations..."
if ! command -v seqkit &> /dev/null || ! command -v trimmomatic &> /dev/null || \
   ! command -v vsearch &> /dev/null || ! command -v blastn &> /dev/null; then
    echo "Some tools were not installed correctly. Please check the logs and try again." >&2
    exit 1
fi

echo "Make sure you have the database files in the database folder,"
echo " if not, download them from the Enterotool GitHub repository."
echo "or made your own database files with blast -makeblastdb command."
echo "#####################################################"
echo "Environment setup completed successfully!"
echo "Activate the environment with: conda activate $ENV_NAME"
echo "You can now run the Enterotool-pipeline script."
