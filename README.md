# `ld12` version 1.0.1

Automated analysis and visualization generation of the LD12 genomes and the freshwater genomes comparison paper.

* All code written by [Lucas Sinclair](http://envonautics.com/#lucas).

### Publication
The published paper for which this pipeline was made can be found here:

http://www.nature.com/ismej/journal/vaop/ncurrent/full/ismej2015260a.html

### Flowchart
The code is documented with many docstrings, in addition here is overview of what happens in the analysis:

![Flowchart](/../master/documentation/flowchart.png?raw=true "Flowchart")

### Optional arguments
The command line tool supports a few optional arguments:

* `e_value` : Minimum e-value in similarity search. Defaults to 0.0001.
* `mcl_factor` : The MCL clustering factor. Defaults to 1.5.
* `seq_type` : Either `nucl` or `prot`. Defaults to `prot`.
* `num_threads` : Number of threads to use. Default to the number of cores on the current machine.
* `min_identity` : Minimum identity in similarity search. Defaults to 0.97.
* `min_coverage` : Minimum query coverage in similarity search. Defaults to 0.97.