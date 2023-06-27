# WGS fusion validation pipeline


## Requirements

A linux OS is required to run Singularity. The pipeline can be run without Singularity but requires an OS capable of running samtools, bedtools and novoalign. It is possible to set up a conda environment identical to the one in the singularity containers by copying the conda installs section of the singularity .def file in a new environment:

```
conda install -c conda-forge mamba

mamba install -c bioconda -c conda-forge -c anaconda -y samtools=1.13 bedtools=2.29.2 r-tidyverse=1.3.1 bioconductor-iranges=2.28.0 novoalign=3.09.00
```

## Input files

The user needs to supply two input files. The first is an indexed fasta file of the reference genome that the WGS data to be queried is aligned to. These should be located under input/ref_genome/  in .fa and .fai formats.

The second file to be supplied is a table of fusions that are to be validated. This should be a tab-delimited table with the following columns (case sensitive!):

|Input column|Description|Example|
|------------|-----------|-------|
|fusion_id | Unique identifier for each fusion transcript | 1|
|sample_id | Unique identifier / barcode for each sample | sample_01|
|path |Absolute or bound path to the raw WGS files | /disk/data/sample1/sample1.bam|
|fiveprime_chr |5' chromosome| 8|
|fiveprime_strand |5' strand| -|
|fiveprime_search_start | 5' fusion search region start coordinate| 54956927|
|fiveprime_search_end | 5' fusion search region end coordinate |55013968|
|threeprime_chr | 3' chromosome |8|
|threeprime_strand | 3' strand | + |
|threeprime_search_start | 3' fusion search region start coordinate | 61427495 |
|threeprime_search_end | 3' fusion search region end coordinate| 61531639|

Other columns are accepted but will be ignored. The order of the columns does not matter. An example of an input table can be found in the input/ directory.

In order for breakpoints to be identified, the aligned reads in the WGS data must not hard-clipped. Data preprocessing is a crucial step for the pipeline to run correctly. Coordinates for fusion transcripts and gene annotation must be in the same assembly of the human genome as the aligned WGS data, and tools such as LiftOver from UCSC can be used to convert gene and fusion junction coordinates to match WGS data.

## Running the pipeline

If running locally, the recommended way of installing singularity is via conda. A simple environment can be create using the following command:

`conda create -n singularity singularity=3.8.6`

The singularity container can subsequently be built using the following:

`singularity build --fakeroot singularity/fusion_pipeline.sif singularity/fusion_pipeline.def`

The full fusion validation pipeline can then be run using:

`singularity exec singularity/fusion_pipeline.sif Rscript scripts/fusion_validation_pipeline.R`

Singularity containers, by default, provide an isolated environment with limited access to the host system's files and resources for security reasons. If the raw sequencing data files are stored outside of this environment they will have to be bound to be accessed by singularity. This can be done by modifying the .def file before building the container, or by adding the optional  `--bind` option to the container.

As an example, assume the WGS data is stored in the directory "data" located under "/disk/projects/data". To bind the data to the singularity, the following command would be run:

`singularity exec --bind /disk/projects/data:/data singularity/fusion_pipeline.sif Rscript scripts/fusion_validation_pipeline.R`

In this example, the "--bind" option is used to bind the "/disk/projects/data" directory on the host system to the "/data" directory inside the container. By doing so, the container will have access to the files within that directory. The "path" column in the input file fusion_search_table.txt should reflect the bound data structure.
