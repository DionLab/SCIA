Single Clone based Instability Assay (SCIA) GUI

SCIA GUI is a cross-platform graphical application for analyzing repeat instability from long read sequencing data. 
It integrates the deterministic Repeat Detector (RD) algorithm with downstream instability quantification, delta distribution analysis, statistical comparison and visualization.

This cross-platform software provides complete workflow from raw FASTA files to instability indices and comparitive statistical outputs without the need of command line input(CLI).

The current repeat detector profile is configured for Huntington disease(CAG repeat disorder).

**Core Features**

**Repeat Detector Integration**

-Supports restrictive and permissive profiles

-Accepts unaligned FASTA input

-Uses deterministic profile weighting (pfsearch-based)

-Reverse complement analysis supported

**Outputs:**

Histogram of repeat size distribution

Threshold-filtered histograms

Histogram plots

**Instability Index Analysis**

Quantifies repeat instability using weighted distribution shifts relative to a reference mode.

Three instability modes supported:

Fixed Mode from Control (D0):Uses control mode to analyse changes in the treatment samples

Per-Sample Mode: Uses detected mode independently for each sample

Manual Mode: User-specified repeat size

Outputs:

