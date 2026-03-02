# MEIsensor (v1.0)
MEIsensor is a deep-learning based detection and calssification of mobile element insertions from long-read sequencing data.
## Overview
MEIsensor is a deep-learning based framework for the detection and classification of mobile element insertions (MEIs) from long-read sequencing data. It is designed to accurately identify and subtype Alu, LINE1, and SVA insertions.
![Overview of MEIsensor](figure/figure1.png)
## License

## Installation
We recommend using the conda virtual environment to install MEIsensor (Platform: Linux).
```bash
git clone https://github.com/Yudong-XJTU/MEIsensor.git
cd MEIsensor
```
If your CUDA version is higher than 12.8, you can directly install the environment using:
```bash
conda env create -f environment.yml -n your_env_name
```
Alternatively, you can follow the steps below to install the environment manually. This is especially recommended for users with lower CUDA versions, as you may need to manually adjust the PyTorch version and installation source.
```bash
# Create a conda environment for MEIsensor
conda create -n MEIsensor python=3.10
conda activate MEIsensor
pip install torch==2.7.0 torchvision==0.22.0 torchaudio==2.7.0 --index-url https://download.pytorch.org/whl/cu128
# Install other packages
pip install -r requirements.txt
```
Check if CUDA is available:
```bash
python -c "import torch; print(torch.__version__, torch.version.cuda, torch.cuda.is_available())"
```
## Usage

```bash
python src/main.py -i $BAM_PATH -v $OUTPUT_PATH --reference $REF_PATH -m $MODEL_PATH
```
### required parameters
```bash
-i BAM_PATH                       Absolute path to output
-v OUTPUT_PATH                    Absolute path to bam file
--reference REF_PATH(Optional)    Absolute path to your reference genome
-m MODEL_PATH                     Absolute path to trained model
-t THREAD                         Number of parallel threads to use (speed-up for multi-core CPUs)
```