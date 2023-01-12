# Anansnake

Link seq2science output to ANANSE with 2 sample tables and one config file.

![](docs/img/anansnake.PNG)

## Installation
```bash
mamba create -n anansnake -c bioconda anansnake
```

Don't forget to activate the conda environment with `mamba activate anansnake`.

## Running anansnake on the example data
The anansnake github repository contains an `example` folder which can be downloaded to try the workflow.
Here we assume you've downloaded the folder in your current working directory.
Check [it's README](https://github.com/vanheeringen-lab/anansnake/blob/master/example/README.md) for additional details!

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

## Running anansnake
Anansnake works with seq2science in- & output: The RNA- and ATAC-seq `samples.tsv` files are the same you've used for seq2science, with one addition (see below).
The counts tables are output files without any changes.

The RNA- and ATAC-seq samples are combined via a shared column in the samples.tsv files.
In the example data, this is the `anansnake` column.
Which conditions from the `anansnake` column are compared is set in the `config.yaml` file, under `contrasts`. 

For files and settings you can check out the example folder.

## Troubleshooting
ANANSE can take tonnes of memory. If your machine freezes, reduce the number of threads or mem_mb.
