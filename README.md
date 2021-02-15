# Functions for V(D)J sequencing data for SingleCellExperiment objects: `SCEdjvdj`

The folks at the [RNA Bioscience Initiative](https://github.com/rnabioco) have created the [`djvdj`](https://github.com/rnabioco/djvdj) package that contains many useful functions for exploring V(D)J Sequencing data, usually obtained via 10X Genomics and CellRanger.
Here, I've done some very light editing to make the functions work with with `SingleCellExperiment` objects. Note that they will *only* work with SCE objects; if you isnist on Seurat objects, use the original package.


