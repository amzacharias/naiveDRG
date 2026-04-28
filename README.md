-   <a href="#preface" id="toc-preface">Preface</a>
-   <a href="#setup" id="toc-setup">Setup</a>
-   <a href="#figure-5a-workflow" id="toc-figure-5a-workflow">Figure 5a
    workflow</a>
    -   <a href="#helpers" id="toc-helpers">Helpers</a>
    -   <a href="#quantify-transcript-counts"
        id="toc-quantify-transcript-counts">Quantify transcript counts</a>
    -   <a href="#data-preparation" id="toc-data-preparation">Data
        preparation</a>
    -   <a href="#from-candidates-identify-differentially-expressed-genes"
        id="toc-from-candidates-identify-differentially-expressed-genes">From
        candidates, identify differentially expressed genes</a>
-   <a href="#peer-review-supplementary-figure-2-workflow"
    id="toc-peer-review-supplementary-figure-2-workflow">Peer Review
    Supplementary Figure 2 workflow</a>
    -   <a href="#helpers-1" id="toc-helpers-1">Helpers</a>
-   <a href="#scripts" id="toc-scripts">Scripts</a>
-   <a href="#done" id="toc-done">Done!</a>

------------------------------------------------------------------------

## Preface

Data analysis for the manuscript, [“Nociceptor clock genes control
excitability and pain perception in a sex and time-dependent
manner”](https://doi.org/10.1101/2025.04.07.646998)

------------------------------------------------------------------------

## Setup

**Important**:

-   Consider reading the `README.html` file which has a floating table
    of contents.
-   This project assumes you are using resources from the [The Centre
    for Advanced Computing](https://cac.queensu.ca/), which uses a SLURM
    job scheduler.
    -   It is highly recommended that you use a cloud computing system.
        You may need to edit scripts to load dependencies in a manner
        compatible with your system.
-   Ensure all scripts and data are stored in an R project folder.
-   Script names are numbered so the order of execution is more obvious.
-   Set the R current working directory to the project working
    directory. Most scripts assume that the project directory is the
    current working directory.
-   <mark>Caution! Some scripts use absolute paths (especially bash
    scripts) </mark>
    -   Run the following commands in the terminal to replace the
        `absolutePath` spaceholder found in scripts with your absolute
        path to the project directory.

    <!-- -->

        find . -type f -name "*.sh" -exec sed -i'' -e 's#absolutePath#/my/custom/path#g' {} +
        find . -type f -name "*.R" -exec sed -i'' -e 's#absolutePath#/my/custom/path#g' {} +

**Figure 5a session info**:

-   R version 3.6.0 (2019-04-26)
-   Platform: x86_64-redhat-linux-gnu (64-bit)
-   Running under: CentOS Linux 7 (Core)
-   Matrix products: default
-   BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so

**Peer Review Supplementary Figure 2 session info**:

-   R version r/4.5.0 (2019-04-26)
-   Platform: x86_64-pc-linux-gnu
-   Running under: CRocky Linux 8.10 (Green Obsidian)
-   Matrix products: default
-   BLAS/LAPACK: FlexiBLAS IMKL; LAPACK version 3.11.0

------------------------------------------------------------------------

## Figure 5a workflow

> Figure5a directory.

### Helpers

Notice the `../0_helpers` folder. This directory contains many R
functions that minimize repetition of code and are generally helpful.

### Quantify transcript counts

1.  Navigate to `./1_stringtie`

2.  Run `1_writePass1Scripts.R` to write individual scripts for pass 1.
    Execute scripts in the `pass1IndivScripts` directory. Use the
    `2_checkSuccess.R` and `jobsToRun.sh` scripts to monitor progress.

        # 1 cpu, 5 GB memory 
        # REF_GTF is the full GTF file from Gencode
        module load StdEnv/2020 stringtie/2.1.5
        stringtie $INPUT -p 5 -G $REF_GTF -o $OUT_GTF

3.  Run `3_writeGtfLists.R` to prepare the merging of individual GTFS
    from pass 1.

4.  Run `*.sh*` files in the `3_merge` folder to execute the merging of
    GTF files.

        # 5 cpu, 3 GB memory
        module load StdEnv/2020 stringtie/2.1.5
        stringtie --merge -p 20 -o $OUTPUT -G $REF_GTF $GTFS_LIST

5.  Evaluate StringTie performance with `4.1_writeGffCompareScripts.R`
    and `*.sh` scripts in the `4_gffCompare` folder.

6.  Run `5_writePass2Scripts.R` to write individual scripts for pass 2.
    Execute scripts in the `pass2IndivScripts` directory. Use the
    `2_checkSuccess.R` and `jobsToRun.sh` scripts to monitor progress.

        # 1 cpu, 5 GB memory
        # REF_GTF is the merged gtf that corresponds to this sample's tissue
        module load StdEnv/2020 stringtie/2.1.5
        stringtie $INPUT -b $BALL -e -p 5 -G $REF_GTF -o $OUT_GTF

7.  To get transcript id to gene name mapping, run
    `6_isoformAnalyzeR/isoformAnalyzeR.R`.

### Data preparation

1.  Navigate to `2_dataPrep`
2.  Clean the count matrix

<!-- -->

1.  Run `0_id2name.R` to get a dataframe with ensembl ID to gene
    name/symbol conversion information.
2.  Run `1_outlierRemoval.R` to … a. Perform outlier detection with
    *arrayQualityMetrics*. A sample is considered an outlier if
    -   it is marked as an outlier before and after normalization by the
        same outlier detection metrics, and/or,
    -   it is marked as an outlier by multiple outlier detection metrics
        after normalization
    -   Note: No samples were considered outliers and removed.
        2.  Normalize counts with the [median of ratios
            method](https://doi.org/10.1186/gb-2010-11-10-r106)
3.  Run `2_filtering.R` to perform non-specific filtering to remove
    lowly expressed features. - This and the last step are performed
    more-so to optimize the filtering threshold that will be used for
    DeSeq2, and prepare counts for unknown future analyses.

### From candidates, identify differentially expressed genes

1.  Navigate to `./3_deseqCandidate`
2.  Prepare candidate gene lists. Download gene lists from KEGG with
    `1.0_getKeggGenes.R`. Run `1.1_prepareCanddiates.R` to clean
    candidate lists.
3.  Run `2_writeScripts.R` to write a bash script for each analysis.
    -   Candidate genes are removed after the lowly expressed genes are
        removed.
4.  Run `*.sh` scripts in the `bash` folder to execute analyses.

## Peer Review Supplementary Figure 2 workflow

> PR_SFigure2 directory.

### Helpers

Notice the `workflow/scripts/helpers` folder. This directory contains
many R functions that minimize repetition of code and are generally
helpful.

## Scripts

Execute computational demanding scripts with bash scripts (`.sh`).

Rscripts being run: `dea_Markers.R` and `dea_nav18.R`

## Done!
