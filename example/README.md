## Setup

1. clone this repository with `git clone https://github.com/vanheeringen-lab/anansnake.git`
2. `cd` into the anansnake directory
3. create a conda environment with anansnake (`mamba create -n anansnake anansnake` or `mamba env create -n anansnake -f requirements.yaml`)
4. activate the conda environment with `mamba activate anansnake`

## Running anansnake on the example data
To check if everything is set up right, we can do a dry run:
```bash
anansnake --configfile example/config.yaml --dry-run
```
If you get an error, be sure to check the red text!
I've added human-readable feedback where I could.

With the example data, you still need to provide a genome and gene annotation.
You can download the GRCz11 genome and gene annotation with 
```bash
genomepy install GRCz11 --annotation
```
If you do another dry run you should have no more errors, and see what is going to happen.

To do the real run, you need to specify how many cores you want to use, and how much RAM you have:
```bash
anansnake --configfile example/config.yaml --resources mem_mb=48_000 --cores 12
```
