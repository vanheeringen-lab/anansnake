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
Check [it's README](https://github.com/vanheeringen-lab/anansnake/blob/master/example/README.md) for details!

## Running anansnake

Anansnake works with seq2science in- & output: The RNA- and ATAC-seq `samples.tsv` files are the same you've used for seq2science, with one addition (see below).
The counts tables are output files without any changes.

The RNA- and ATAC-seq samples are combined via a shared column in the samples.tsv files.
In the example data, this is the `anansnake` column.
Which conditions from the `anansnake` column are compared is set in the `config.yaml` file, under `contrasts`. 

For files and settings & command line examples you can check out the [example folder](https://github.com/vanheeringen-lab/anansnake/blob/master/example).

## Troubleshooting

ANANSE can take tonnes of memory. If your machine freezes, reduce the number of threads or mem_mb.
