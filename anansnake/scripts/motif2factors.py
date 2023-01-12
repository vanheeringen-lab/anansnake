from contextlib import redirect_stdout, redirect_stderr
from os.path import dirname, join
from shutil import copyfile


# log errors
with open(str(snakemake.log), "w") as f:
    with redirect_stdout(f), redirect_stderr(f):

        if snakemake.params.get_orthologs is False:
            # copy the default m2f
            print(f"Copying the motif database ({snakemake.params.database})")
            from gimmemotifs.utils import pfmfile_location

            in_pfmfile = pfmfile_location(snakemake.params.database)
            out_pfmfile = snakemake.output[0]

            in_m2ffile = in_pfmfile.replace(".pfm", ".motif2factors.txt")
            out_m2ffile = out_pfmfile.replace(".pfm", ".motif2factors.txt")

            copyfile(in_pfmfile, out_pfmfile)
            copyfile(in_m2ffile, out_m2ffile)

        else:
            # create an ortholog m2f
            from gimmemotifs.orthologs import motif2factor_from_orthologs

            tmpdir = snakemake.params.tmpdir
            if tmpdir is not None:
                tmpdir = join(tmpdir, "motif2factors")
            
            motif2factor_from_orthologs(
                database=snakemake.params.database,
                new_reference=[snakemake.input.genome],
                genomes_dir=snakemake.params.genomes_dir,
                outdir=dirname(snakemake.output[0]),
                tmpdir=tmpdir,
                keep_intermediate=snakemake.params.keeptmp,
                threads=snakemake.threads,
            )
