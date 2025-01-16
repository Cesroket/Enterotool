
#!/bin/bash

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
echo "Environment setup completed successfully!"
echo "Activate the environment with: conda activate $ENV_NAME"
echo "You can now run the Enterotool-pipeline script."