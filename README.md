<!--
Copyright (C) 2024 Roberto Rossini <roberros@uio.no>

SPDX-License-Identifier: MIT
-->

# Nextflow workflow to run NCHG

[![CI](https://github.com/paulsengroup/nchg-nf/actions/workflows/ci.yml/badge.svg)](https://github.com/paulsengroup/nchg-nf/actions/workflows/ci.yml)

This repository hosts a Nextflow workflow to identify statistically significant interactions from Hi-C data using [NCHG](https://github.com/paulsengroup/NCHG).


## Requirements

### Software requirements

- Nextflow (at least version: v22.10.8. Pipeline was developed using v24.04.2)
- Docker or Singularity/Apptainer

### Required input files

The workflow can be run in two ways:
1. Using a sample sheet (recommended, supports processing multiple samples at once)
2. By specifying options directly on the CLI or using a config

#### Using a samplesheet

The samplesheet should be a TSV file with the following columns:

| sample       | hic_file                               | resolution | domains                  | mask              |
|--------------|----------------------------------------|------------|--------------------------|-------------------|
| sample_name  | myfile.hic                             | 50000      | tads.bed                 | mask.bed          |
| 4DNFI74YHN5W | 4DNFI74YHN5W.mcool::/resolutions/50000 | 50000      | 4DNFI74YHN5W_domains.bed | assembly_gaps.bed |


- __sample__: Sample names/ids. This field will be used as prefix to in the output file names (see [below](#running-the-workflow)).
- __hic_file__: Path to a file in .hic or Cooler format.
- __resolution__: Resolution to be used for the data analysis (50-100kbp are good starting points).
- __domains__ (optional) : path to a BED3+ file with a list of pre-computed domains (e.g. TADs). When provided, NCHG will aggregate interactions between all possible pairs of domains, and compute their statistical significance.
- __mask__ (optional): path to a BED3+ file with the list of regions to be masked out.

URI syntax for multi-resolution Cooler files is supported (e.g. `myfile.mcool::/resolutions/bin_size`).

Furthermore, all contact matrices (as well as domain and mask files when provided) should use the same reference genome assembly.

<details>
<summary> <b>Without using a samplesheet</b> </summary>

To run the workflow without a samplesheet is not available, the following parameters are required:

- __sample__
- __hic_file__
- __resolution__

Parameters have the same meaning as the header fields outlined in the [previous section](#using-a-samplesheet).

The above parameters can be passed directly through the CLI when calling `nextflow run`:

```bash
nextflow run --sample='4DNFI74YHN5W' \
             --hic_file='data/4DNFI74YHN5W.mcool' \
             --resolution=50000
             ...
```

Alternatively, parameters can be written to a `config` file:
```console
user@dev:/tmp$ cat myconfig.txt

sample       = '4DNFI74YHN5W'
hic_file     = 'data/4DNFI74YHN5W.mcool'
resolution   = 50000
```

and the `config` file is then passed to `nextflow run`:
``` bash
nextflow run -c myconfig.txt ...
```

</details>

### Optional files and parameters

In addition to the mandatory parameters, the pipeline accepts the following parameters:

- __cytoband__: path to a [cytoband](https://software.broadinstitute.org/software/igv/cytoband) file. Used to mask centromeric regions.
- __assembly_gaps__: path to a BED file with the list of assembly gaps/unmappable regions.

Note that NCHG by default uses the `MAD-max` filter to remove bins with suspiciously high or low marginals, so providing the above files is usually not requirerd.

- __mad_max__: cutoff used by NCHG when performing the `MAD-max` filtering.
- __bad_bin_fraction__: bad bin fraction used by NCHG to discard domains overlapping with a high fraction of bad bins.
- __fdr_cis__: adjusted pvalue used by NCHG to filter significant cis interactions.
- __log_ratio_cis__: log ratio used by NCHG to filter significant cis interactions.
- __fdr_trans__: adjusted pvalue used by NCHG to filter significant trans interactions.
- __log_ratio_trans__: log ratio used by NCHG to filter significant trans interactions.
- __use_cis_interactions__: use interactions from the cis portion of the Hi-C matrix.
- __use_trans_interactions__: use interactions from the trans portion of the Hi-C matrix.

By default, the workflow results are published under `result/`. The output folder can be customized through the __outdir__ parameter.

For a complete list of parameters supported by the workflow refer to the workflow main [config](nextflow.config) file.

## Running the workflow

First, download the example datasets using script `utils/download_example_datasets.sh`.

```bash
# This will download files inside folder data/
utils/download_example_datasets.sh data/
```

Next, create a `samplesheet.tsv` file like the follwing (make sure you are using tabs, not spaces!)

```tsv
sample   hic_file      resolution    domains    mask
example  data/4DNFI74YHN5W.mcool   50000    
```

Finally, run the workflow with:
```console
user@dev:/tmp$ nextflow run --max_cpus=8 \
                            --max_memory=16.GB \
                            --max_time=2.h \
                            --sample_sheet=samplesheet.tsv \
                            --outdir=data/results/ \
                            https://github.com/paulsengroup/nchg-nf \
                            -with-singularity  # Replace this with -with-docker to use Docker instead

 N E X T F L O W   ~  version 24.04.2

Launching `./main.nf` [golden_blackwell] DSL2 - revision: e923a03f8f

-- PARAMETERS
-- sample_sheet: samplesheet.tsv
-- outdir: results/
-- publish_dir_mode: copy
-- cytoband: null
-- assembly_gaps: null
-- mad_max: 5
-- bad_bin_fraction: 0.1
-- fdr_cis: 0.01
-- log_ratio_cis: 1.5
-- fdr_trans: 0.01
-- log_ratio_trans: 1.5
-- plot_format: png
-- hic_tgt_resolution_plots: 500000
-- plot_sig_interactions_cmap_lb: null
-- plot_sig_interactions_cmap_ub: 2.0
-- skip_expected_plots: false
-- skip_sign_interaction_plots: false
executor >  local (2)
[03/3b02c4] SAMPLESHEET:CHECK_SYNTAX                 [100%] 1 of 1, cached: 1 ✔
[da/c8ed7a] SAMPLESHEET:CHECK_FILES                  [100%] 1 of 1, cached: 1 ✔
[97/828730] NCHG:GENERATE_MASK (example)             [100%] 1 of 1, cached: 1 ✔
[3e/0ee8aa] NCHG:EXPECTED (example)                  [100%] 1 of 1, cached: 1 ✔
[f2/084a7b] NCHG:GENERATE_CHROMOSOME_PAIRS (example) [100%] 2 of 2, cached: 2 ✔
[ed/2b3e57] NCHG:DUMP_CHROM_SIZES (example)          [100%] 1 of 1, cached: 1 ✔
[60/84727e] NCHG:COMPUTE (example (chr18:chrY))      [100%] 231 of 231, cached: 231 ✔
[37/0fd95d] NCHG:MERGE (example (trans))             [100%] 2 of 2, cached: 2 ✔
[af/f0ba9a] NCHG:FILTER (example (cis))              [100%] 2 of 2, cached: 2 ✔
[6f/600f8c] NCHG:VIEW (example (cis))                [100%] 2 of 2, cached: 2 ✔
[de/a0069a] NCHG:CONCAT (example)                    [100%] 1 of 1, cached: 1 ✔
[66/880021] NCHG:PLOT_EXPECTED (example)             [100%] 1 of 1 ✔
[9a/0075b8] NCHG:GET_HIC_PLOT_RESOLUTION (example)   [100%] 1 of 1, cached: 1 ✔
[a3/779d15] NCHG:PLOT_SIGNIFICANT (example)          [100%] 1 of 1 ✔
Completed at: 03-Jul-2024 18:48:25
Duration    : 1m 42s
CPU hours   : 0.6 (23.7% cached)
Succeeded   : 2
Cached      : 246
```

This will create a `data/results/` folder with the following files:
- `example.filtered.tsv.gz` - TSV with the statistically significant interactions detected by NCHG.
- `expected_values_example.cis.h5` - HDF5 file with the expected values computed by NCHG.
- `plots/example/example.*.*.png` - Plots showing the log ratio computed by NCHG for each chromosome pair analyzed.
- `plots/example/example_cis.png` - Plot showing the expected value profile computed by NCHG.

<details>
<summary>Troubleshooting</summary>

If you get permission errors when using `-with-docker`:
- Pass option `-process.containerOptions="--user root"` to `nextflow run`

If you get an error similar to:
```
Cannot find revision `v0.4.0` -- Make sure that it exists in the remote repository `https://github.com/paulsengroup/nchg-nf`
```

try to remove folder `~/.nextflow/assets/paulsengroup/nchg-nf` before running the workflow

</details>

## Getting help

If you are having trouble running the workflow feel free to reach out by starting a new discussion [here](https://github.com/paulsengroup/nchg-nf/discussions).

Bug reports and feature requests can be submitted by opening an [issue](https://github.com/paulsengroup/nchg-nf/issues).
