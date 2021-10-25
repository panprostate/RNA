# PPCG RNA working group

This repository contains the pipeline used by the Pan Prostate Cancer Group.

## **Requirements**
All tools used in this pipeline are installed at runtime through a singularity container. The only requirements are the ones to run the Snakemake workflow.

- mamba >= 1.2 (optional, for faster installation)
- Snakemake >= 6.2
- Singularity >=3.7 (Should work on previous versions although not tested)
- Cookiecutter >= 1.7.2


With exception of Singularity, these can be easily installed using [conda](https://docs.conda.io/en/latest/miniconda.html).
**Please do not install Singularity through conda**. 

Singularity must be owned by `root` in order to run this workflow. If you do not have access to a singularity installation owned by `root` [try requesting to your system administrator](https://sylabs.io/guides/3.7/user-guide/quick_start.html#installation-request). 

## **Before you begin**
When executing snakemake, try to avoid canceling it with `ctrl + c`. If an error occurs, snakemake will clean the files from the failed step before shutting down. Forcing it to shutdown will cause corrupted files to be left behind, potentially causing troubles. 
On the same note, avoid moving/creating/removing files under the `results/` folder. Snakemake relies on the time stamp of the files to coordinate the order of execution of each step. If you must remove files, only remove files resulting from the last step ran, otherwise this will cause snakemake to reprocess all steps that relies on these files.

## **How to run**
1. Download the snakemake workflow.
```
wget https://github.com/panprostate/RNA/releases/download/v1.2.2/PPCG_RNA_v1.2.2.tar.gz
tar -xzvf PPCG_RNA_v1.2.2.tar.gz
```

2. This workflow requires a sample metadata table (tab separated) containing the following fields:
- sample_name = The sample name, by default anything preceding the "_L00" pattern.
- unit_name = The unit/lane name, by default it matches the following pattern ".+\_(L00.)\.+".
- fq1 = Absolute path to the R1 file of the fastq pair.
- fq2 = Absolute path to the R2 file of the fastq pair. For single-ended libraries this field should contain `NA`.

Create this table under the `config/` folder.
As reference, a helper script `createTable.sh` is included under `config/` to create this table. Make sure the table format is correct before proceding.

3. If your computing enviroment is managed by Slurm, a cluster profile is included at `slurm_profile/`. There are two files that can be adjusted:
- `slurm_profile/config.yaml`

Most likely you will want to adjust only two settings in this file:
```
# Number of cores used by Snakemake to manage the workflow - if increasing, be careful to not have your job killed if not running in an interactive session.

local-cores: 1

# Number of maximum simultaneous jobs that can be run.

jobs: 5
```
- `slurm_profile/settings.json`

This file can be used to include center-specific parameters for your submissions. For example, to specify a partition `panda` on all jobs submissions.
```
"SBATCH_DEFAULTS": "partition=panda job-name={rule}_{wildcards} output=job_logs/slurm-%j.out"
```

4. Edit the amount of resources requested and other general settings in `config/config.yaml`. Some of the options are:
**NEW** - library_type: Process single-ended libraries. Please note that a mixed sample table is not supported at the moment, make sure all samples are SE or PE in your `samples` variable.
- samples: path to the metadata table described in step 2.
- buildIndex: True or False wether to build indexes from files or download built indexes.
- adapters: adapter sequences used by cutadapt.
- mem_*: Amount of RAM memory requested for the job in MB.
- ncpus_*: Number of cores requested for the job.
- rt_*: Maximum time a job is allowed to run before it is killed.
- compression: Compression level from 1 to 9 for intermediated steps. Final outputs are always compressed at the default level.

5. From the base directory execute Snakemake:
```
snakemake --use-conda --use-singularity --profile slurm_profile
```

The workflow will automatically pull and create a Singularity container from `docker://condaforge/mambaforge:4.10.1-0` and create environments with the necessary packages for each module.

6. If debugging is necessary, the folder `logs/` contains the error log generated by the software run in each step. In some cases the errors are not directly caused by a software, but by the execution of a script or the cluster resource manager - in those cases the error will be output in `job_logs/` with the name `slurm-<job-number>.out`.

## Troubleshooting
Depending on the number of samples and read length, you might end up with an excessive amount of SJs detect. In this case you will encounter the error bellow:
`Fatal LIMIT error: the number of junctions to be inserted on the fly =XXXXXXX is larger than the limitSjdbInsertNsj=1000000`
While you can increase this limit, this might cause you to encounter another error:
`EXITING because of FATAL ERROR: cannot insert junctions on the fly because of strand GstrandBit problem`
The solution for this is to decrease the number of SJs being inserted on the fly. You can do this by turning on the SJ filter in the `config.yaml` file, which by default will remove SJs supported by only 1 read (uniquely mapped).
`SJ_filter: True`
`SJ_minCounts: 1`

# **Workflow overview**
Assuming a sample named `S1` which has been split in two files by lane (`L001` and `L002`).

![overview](https://github.com/panprostate/RNA/blob/master/workflow-schematics.jpg?raw=true)
