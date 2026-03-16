Single Clone based Instability Assay (SCIA) GUI

SCIA GUI is a cross-platform graphical application for analyzing repeat instability from long read sequencing data. 
It integrates the deterministic Repeat Detector (RD) algorithm with downstream instability quantification, delta distribution analysis, statistical comparison and visualization.

This cross-platform software provides complete workflow from raw FASTA files to instability indices and comparitive statistical outputs without the need of command line input(CLI).

The current repeat detector profile is configured for Huntington disease(CAG repeat disorder).

**Core Features**

**Repeat Detector tab**

Supports restrictive and permissive profiles

Accepts unaligned FASTA input

Uses deterministic profile weighting (pfsearch-based)

Reverse complement analysis supported

**Instability Index Analysis**

Quantifies repeat instability using weighted distribution shifts relative to a reference mode

Three instability modes supported:

Fixed Mode from Control (D0):Uses control mode to analyse changes in the treatment samples

Per-Sample Mode: Uses detected mode independently for each sample

Manual Mode: User-specified repeat size

**Outputs**

1. Histogram folder with histograms of repeat size distribution and threshold-filtered histograms
2. Within histograms folder, automatic folder organisation based on the sample name (eg. D0, day0, con, ctrl, control and rest as treated)
3. Raw histogram plots
4. instability_results
5. d0_mode_info

**Delta plot tab**

Supports plotting of one or overlaying two delta plots for comparison

Delta plot 1: Accepts folders containing histograms of control and treatment samples

Delta plot 2:when enabled, accepts folders containing histograms of control and treatment samples

User can modify the legend label, plot settings and output folder paths.

**Outputs**

For each dataset:
1. Binned by offset file containing bin offset from the chosen mode, bin center, delta sum value aling with area-under-curve(AUC)
2. Normalized with offset file with freq of repeats normalised by total number of reads per sample, average, sd, sem, delta and bin offset values
3. Stats file with expansion and contraction mean, AUC, overall mean and bias ratio
4. Kolmogorov–Smirnov test results
5. delta plot as png


**GUI for Windows and macOS**

**Windows**

1. Download [Docker Desktop for Windows](https://docs.docker.com/desktop/install/windows-install/), install it, and ensure it is running.
2. Download the Repeat Detector Docker image from [Zenodo]((https://zenodo.org/records/18863035/files/repeat-detector.tar?download=1))
3. Download GUI  [Zenodo](https://zenodo.org/records/18863035/files/GUI_publish_v2_Mar2026_Win.zip?download=1)
4. Copy the GUI files and run them locally from your computer rather than from networked folders.
5. Locate, unzip 'GUI_publish_v2_Mar2026_Win.exe' and right-click and **Run as administrator**
6. If Windows shows a security warning, click **More info**, then click **Run anyway**.
7. When the app opens, browse and locate the downloaded `repeat-detector.tar` file
8. Next, browse and locate the FASTA file. If you have FASTQ/FASTQ.gz files, use `seqkit fq2fa` to convert them to FASTA format.
9. Choose output base folder.
10. Next, select D0 or control fasta files, make sure the name has the string (D0, Day0, Ctrl, control, ctl, this is needed for automatic folder organisation of folders)
11. Next choose appropriate profile- restrictive or permissive
12. Choose appropriate Instability mode and max repeat size and thershold for filtering.
13. Click **Run RD Program**
14. It may take a minute or two for the Docker image to load.
15. In delta plot tab, choose appropriate control and treated folders for Dataset 1 and Dataset 2.
16. Adjust plot parameters as per the experiment.
17. Choose the result folder, this is where all the results will be stored.
18. Click **Run Delta + KS Test** button to process. (based on screen size, this button might be hidden, to fix this, adjust the scale: System>Display> Scale to 100%)



**macOS**

1. Download Docker Desktop for macOS from [Docker Desktop for Mac](https://docs.docker.com/desktop/install/mac-install/), install it, and ensure it is running.
2. Download the Repeat Detector Docker image from [Zenodo]((https://zenodo.org/records/18863035/files/repeat-detector.tar?download=1))
3. Dowload GUI [Zenodo](https://zenodo.org/records/18863035/files/GUI_publish_v1_Feb2026_mac.zip?download=1)
4. Copy the GUI files and run them locally from your computer rather than from networked folders.
5. Locate, unzip and double-click 'GUI_publish_v1_Feb2026_mac'.
6. If macOS shows a security warning, go to System Preferences > Security & Privacy, and allow the app to run by clicking **Open anyway**
7. It might take a few minutes to open, if it does not open with first double clicking, repeat the clicking and wait
8. When the app opens, browse and locate the downloaded `repeat-detector.tar` file
9. Next, browse and locate the FASTA file. If you have FASTQ/FASTQ.gz files, use `seqkit fq2fa` to convert them to FASTA format.
10. Choose output base folder.
11. Next, select D0 or control fasta files, make sure the name has the string (D0, Day0, Ctrl, control, ctl, this is needed for automatic folder organisation of folders)
12. Next choose appropriate profile- restrictive or permissive
13. Choose appropriate Instability mode and max repeat size and thershold for filtering.
14. Click **Run RD Program**
15. It may take a minute or two for the Docker image to load.
16. In delta plot tab, choose appropriate control and treated folders for Dataset 1 and Dataset 2.
17. Adjust plot parameters as per the experiment.
18. Choose the result folder, this is where all the results will be stored.
19. Click **Run Delta + KS Test** button to process.



