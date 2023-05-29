# WGS fusion validation pipeline

## Running the pipeline

The recommended way of installing singularity is via conda. A simple environment can be create using the following command:
`conda create -n singularity singularity=3.8.6`

The singularity container can subsequently be build using the following:
`singularity build –fakeroot singularity/fusion_pipeline.sif singularity/fusion_pipeline.def`

The full fusion validation pipeline can then be run using
`singularity exec singularity/fusion_pipeline.sif Rscript scripts/run_pipeline.R` 

Singularity containers, by default, provide an isolated environment with limited access to the host system's files and resources for security reasons. If the raw sequencing data files are stored outside of this environment they will have to be bound to be accessed by singularity. This can be done by modifying the .def file before building the container, or by adding the optional  `--bind` option to the container. 

As an example, assume the WGS data is stored in the directory “data” located under “/disk/projects/data”. To bind the data to the singularity, the following command would be run:

`singularity exec --bind /disk/projects/data:/data fusion_pipeline.sif Rscript scripts/run_pipeline.R`

In this example, the "--bind" option is used to bind the "/disk/projects/data" directory on the host system to the "/data" directory inside the container. By doing so, the container will have access to the files within that directory. The “path” column in the input file fusion_search_table.txt should reflect the bound data structure. 



