#!/usr/bin/env python3
"""
Enhanced RD + Delta Plot GUI (MODIFIED VERSION)

Tab 1 NEW ADDITIONS (2026-02-13):
- Instability Index can now run in 3 modes:
  1) Fixed mode from D0 (existing behaviour)
  2) Normal (mode per sample)
  3) Manual forced mode (user-entered integer)

2026-02-25 PATCH:
- Fully supports BOTH .rest.histogram (restrictive) and .per.histogram (permissive)
  across:
    * Tab 1 histogram generation + copy + cleanup
    * Enhanced plots
    * Thresholded histogram generation
    * Instability runner (all 3 mode strategies)
    * Organize histograms into Control/Treated for delta tab
- Delta tab already reads *.histogram so it works for both, but label cleaning
  now strips ".rest" / ".per" tags from sample stems for cleaner legends.
"""

import importlib.util
import subprocess
import sys
import time


def install_and_import(package, import_name=None):
    if import_name is None:
        import_name = package

    # Check if running as a PyInstaller frozen executable
    if getattr(sys, 'frozen', False):
        # Running as compiled .exe - just import, don't install
        try:
            return __import__(import_name)
        except ImportError:
            print(f"ERROR: Required package '{package}' not found in executable!")
            return None
    else:
        # Running as normal Python script - can install if needed
        if importlib.util.find_spec(import_name) is None:
            print(f"Package '{package}' not found. Installing...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])
        return __import__(import_name)


print("Loading GUI components...")
wx = install_and_import("wxpython", "wx")
pubsub = install_and_import("PyPubSub", "pubsub")

import os
import re
import glob
import platform
import threading
import warnings
import shutil

print("Loading data processing libraries...")
import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.backends.backend_pdf import PdfPages
from pubsub import pub

scipy = install_and_import("scipy")
from scipy.stats import ks_2samp
from scipy.signal import find_peaks

statsmodels = install_and_import("statsmodels")
from statsmodels.stats.multitest import multipletests

# Install and import scienceplots for publication-ready styling
try:
    import scienceplots  # registers the styles
    available = set(plt.style.available)
    if "science" in available:
        plt.style.use(["science", "no-latex"])
        print("Using scienceplots style for publication-ready plots")
    else:
        print("scienceplots imported but 'science' style not available in this environment. Falling back.")
        plt.style.use("default")
except Exception as e:
    print(f"scienceplots unavailable ({e}). Falling back to default style.")
    plt.style.use("default")

# Set Arial font for publication-ready plots
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 10
plt.rcParams["axes.linewidth"] = 1.0

print("GUI ready!")


# -----------------------------------
# Histogram suffix helper
# -----------------------------------
def hist_suffix_for_profile(profile: str) -> str:
    """Return histogram suffix based on RD profile selection."""
    return ".rest.histogram" if profile == "restrictive" else ".per.histogram"


# -----------------------------------
# Screen utilities
# -----------------------------------
def get_screen_dimensions():
    try:
        display = wx.Display(0)
        screen_rect = display.GetGeometry()
        width = int(screen_rect.width * 0.85)
        height = int(screen_rect.height * 0.85)
        width = max(1000, min(width, 1600))
        height = max(700, min(height, 1200))
        return width, height, screen_rect
    except Exception as e:
        print(f"Error getting screen dimensions: {e}")
        return 1200, 800, None


def center_window_on_screen(window):
    try:
        display = wx.Display(0)
        screen_rect = display.GetGeometry()
        w = window.GetSize()
        x = screen_rect.x + (screen_rect.width - w.width) // 2
        y = screen_rect.y + (screen_rect.height - w.height) // 2
        window.SetPosition((x, y))
    except Exception:
        window.Center()


def created_by_banner(parent, text="Created by Ruban Rex and Angus Dixon"):
    banner_panel = wx.Panel(parent)
    banner_panel.SetBackgroundColour("#D3D3D3")
    sizer = wx.BoxSizer(wx.HORIZONTAL)
    credit_text = wx.StaticText(banner_panel, label=text)
    credit_text.SetForegroundColour("#333333")
    credit_text.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
    sizer.AddStretchSpacer()
    sizer.Add(credit_text, flag=wx.ALIGN_CENTER_VERTICAL | wx.ALL, border=5)
    sizer.AddStretchSpacer()
    banner_panel.SetSizer(sizer)
    return banner_panel


# -----------------------------------
# Clean up histograms
# -----------------------------------
def cleanup_rest_histograms_from_fasta(fasta_dir):
    try:
        for file in os.listdir(fasta_dir):
            # supports both outputs
            if file.endswith((".rest.histogram", ".per.histogram")):
                file_path = os.path.join(fasta_dir, file)
                try:
                    os.remove(file_path)
                    print(f"✓ Removed from FASTA directory: {file_path}")
                except Exception as e:
                    print(f"✗ Failed to remove {file_path}: {e}")
    except Exception as e:
        print(f"Error during cleanup from FASTA directory: {e}")


def cleanup_rest_histograms_from_output(output_folder):
    try:
        for root, dirs, files in os.walk(output_folder):
            for file in files:
                # supports both outputs
                if file.endswith((".rest.histogram", ".per.histogram")):
                    file_path = os.path.join(root, file)
                    if not any(folder in root for folder in ["histograms", "raw_histogram_plots"]):
                        try:
                            os.remove(file_path)
                            print(f"✓ Removed from output folder: {file_path}")
                        except Exception as e:
                            print(f"✗ Failed to remove {file_path}: {e}")
    except Exception as e:
        print(f"Error during cleanup from output folder: {e}")


# -----------------------------------
# Instability utilities
# -----------------------------------
def load_histogram(path):
    df = pd.read_csv(path, sep="\t", skiprows=6, header=None, names=["repeat_size", "reads"])
    df = df.dropna()
    df["repeat_size"] = df["repeat_size"].astype(int)
    df["reads"] = df["reads"].astype(int)
    return df


def detect_mode_robust(df, fallback_to_max=True):
    try:
        if df.empty:
            return (0, 0) if fallback_to_max else (_ for _ in ()).throw(ValueError("Empty histogram"))
        if df["reads"].sum() == 0:
            return (0, 0) if fallback_to_max else (_ for _ in ()).throw(ValueError("No reads in histogram"))

        data_range = df["reads"].max() - df["reads"].min()
        prominence_threshold = 1.0 if data_range < 10 else max(0.01 * data_range, 1.0)

        peaks, properties = find_peaks(df["reads"], prominence=prominence_threshold, distance=5)

        if len(peaks) == 0:
            if fallback_to_max:
                max_idx = df["reads"].idxmax()
                mode = df.loc[max_idx, "repeat_size"]
                peak_reads = df.loc[max_idx, "reads"]
                print(f"⚠️ No peaks detected, using maximum at repeat_size={mode} reads={peak_reads}")
                return mode, peak_reads
            raise ValueError("No peaks detected")

        peaks_sorted = sorted(
            [(df.iloc[p]["repeat_size"], df.iloc[p]["reads"]) for p in peaks],
            key=lambda x: x[1],
            reverse=True,
        )
        return peaks_sorted[0][0], peaks_sorted[0][1]

    except Exception as e:
        if fallback_to_max and (not df.empty) and (df["reads"].sum() > 0):
            max_idx = df["reads"].idxmax()
            mode = df.loc[max_idx, "repeat_size"]
            peak_reads = df.loc[max_idx, "reads"]
            print(f"⚠️ Peak detection failed, using maximum at repeat_size={mode} reads={peak_reads}")
            return mode, peak_reads
        print("⚠️ Peak detection failed, returning 0, 0")
        return 0, 0


def compute_instability_fixed_mode(df, fixed_mode, threshold_percent=5):
    if df.empty or fixed_mode == 0:
        return 0.0, 0.0, 0.0, 0, 0, 0.0, 0.0

    if fixed_mode not in df["repeat_size"].values:
        closest_idx = (df["repeat_size"] - fixed_mode).abs().idxmin()
        fixed_mode = df.loc[closest_idx, "repeat_size"]
        print(f"⚠️ Fixed mode not found, using closest: {fixed_mode}")

    peak_reads = df.loc[df["repeat_size"] == fixed_mode, "reads"].iloc[0]
    threshold_value = (threshold_percent / 100.0) * peak_reads

    df_sorted = df.sort_values("repeat_size").reset_index(drop=True)
    mode_idx = df_sorted[df_sorted["repeat_size"] == fixed_mode].index[0]

    left = None
    for i in range(mode_idx - 1, -1, -1):
        if df_sorted.loc[i, "reads"] <= threshold_value:
            left = int(df_sorted.loc[i, "repeat_size"])
            break
    if left is None:
        left = max(df_sorted["repeat_size"].min(), fixed_mode - 5)

    right = None
    for i in range(mode_idx + 1, len(df_sorted)):
        if df_sorted.loc[i, "reads"] <= threshold_value:
            right = int(df_sorted.loc[i, "repeat_size"])
            break
    if right is None:
        right = min(df_sorted["repeat_size"].max(), fixed_mode + 5)

    filtered = df[(df["repeat_size"] >= left) & (df["repeat_size"] <= right)].copy()
    if filtered["reads"].sum() == 0:
        return 0.0, 0.0, 0.0, left, right, threshold_value, peak_reads

    filtered["norm_height"] = filtered["reads"] / filtered["reads"].sum()
    filtered["delta"] = filtered["repeat_size"] - fixed_mode
    filtered["weighted"] = filtered["norm_height"] * filtered["delta"]

    instability = filtered["weighted"].sum()
    contraction = filtered[filtered["delta"] < 0]["weighted"].sum()
    expansion = filtered[filtered["delta"] > 0]["weighted"].sum()

    return instability, contraction, expansion, left, right, threshold_value, peak_reads


def get_average_mode_from_multiple_d0_files(d0_files):
    modes = []
    for d0_file in d0_files:
        try:
            df = load_histogram(d0_file)
            mode, _ = detect_mode_robust(df)
            modes.append(mode)
            print(f"✓ File {os.path.basename(d0_file)}: mode = {mode}")
        except Exception as e:
            print(f"⚠️ Error processing {os.path.basename(d0_file)}: {e}")
            continue

    if not modes:
        raise ValueError("❌ No valid D0 files provided")

    avg_mode = int(round(np.mean(modes)))
    print(f"✓ Average mode from {len(modes)} D0 files: {avg_mode}")
    print(f"✓ Individual modes: {modes}")
    return avg_mode, modes


def organize_histograms_by_group(histogram_folder, control_fasta_files):
    """
    Copies histogram files in histogram_folder into histogram_folder/Control and histogram_folder/Treated.
    Supports .rest.histogram and .per.histogram plus _thresholded.histogram.
    """
    control_folder = os.path.join(histogram_folder, "Control")
    treated_folder = os.path.join(histogram_folder, "Treated")
    os.makedirs(control_folder, exist_ok=True)
    os.makedirs(treated_folder, exist_ok=True)

    all_hist_files = [
        f
        for f in os.listdir(histogram_folder)
        if (
            f.endswith(".rest.histogram")
            or f.endswith(".per.histogram")
            or f.endswith("_thresholded.histogram")
        )
        and not f.startswith(".")
    ]

    control_base_names = [os.path.splitext(os.path.basename(f))[0] for f in control_fasta_files]

    control_hist_files, treated_hist_files = [], []
    for hist_file in all_hist_files:
        hist_path = os.path.join(histogram_folder, hist_file)

        hist_base = hist_file
        if hist_base.endswith("_thresholded.histogram"):
            hist_base = hist_base.replace("_thresholded.histogram", "")
        elif hist_base.endswith(".rest.histogram"):
            hist_base = hist_base.replace(".rest.histogram", "")
        elif hist_base.endswith(".per.histogram"):
            hist_base = hist_base.replace(".per.histogram", "")

        is_control = any(cb in hist_base for cb in control_base_names)

        if is_control:
            dest_path = os.path.join(control_folder, hist_file)
            shutil.copy2(hist_path, dest_path)
            control_hist_files.append(dest_path)
            print(f"✓ Organized {hist_file} -> Control")
        else:
            dest_path = os.path.join(treated_folder, hist_file)
            shutil.copy2(hist_path, dest_path)
            treated_hist_files.append(dest_path)
            print(f"✓ Organized {hist_file} -> Treated")

    return control_hist_files, treated_hist_files, control_folder, treated_folder


# -----------------------------
# Instability runner with 3 strategies (supports rest/per)
# -----------------------------
def run_instability_analysis(
    histogram_folder,
    output_dir,
    threshold_percent=5.0,
    d0_hist_files=None,
    mode_strategy="fixed_from_d0",  # "fixed_from_d0" | "per_sample" | "manual"
    manual_mode=None,
    control_fasta_files=None,
    hist_ext=(".rest.histogram", ".per.histogram"),
):
    """
    Runs instability index with different mode strategies:
      - fixed_from_d0: compute average mode from provided D0 histogram files, use for all samples
      - per_sample: detect mode per histogram and use that mode for that sample
      - manual: use manual_mode for all samples

    Supports BOTH .rest.histogram and .per.histogram (hist_ext).

    Saves:
      - output_dir/instability_results.csv
      - output_dir/d0_mode_info.csv (only if fixed_from_d0)
    """
    results = []
    os.makedirs(output_dir, exist_ok=True)

    if isinstance(hist_ext, str):
        hist_ext = (hist_ext,)

    hist_files = [
        os.path.join(histogram_folder, f)
        for f in os.listdir(histogram_folder)
        if f.endswith(tuple(hist_ext)) and not f.startswith(".")
    ]
    if not hist_files:
        raise ValueError(f"❌ No histogram files found in histogram folder matching: {hist_ext}")

    d0_hist_files = d0_hist_files or []
    d0_hist_set = set(d0_hist_files)

    forced_mode = None
    d0_mode_info = None

    if mode_strategy == "fixed_from_d0":
        if not d0_hist_files:
            raise ValueError("❌ mode_strategy='fixed_from_d0' requires D0 histogram files.")
        valid_d0 = [f for f in d0_hist_files if os.path.exists(f)]
        if not valid_d0:
            raise ValueError("❌ No valid D0 histogram files found.")
        forced_mode, individual_modes = get_average_mode_from_multiple_d0_files(valid_d0)
        d0_mode_info = pd.DataFrame(
            {
                "D0_Histogram_File": [os.path.basename(f) for f in valid_d0],
                "Individual_Mode": individual_modes,
                "Average_Mode_Used": [forced_mode] * len(valid_d0),
            }
        )
    elif mode_strategy == "manual":
        if manual_mode is None:
            raise ValueError("❌ mode_strategy='manual' requires manual_mode.")
        forced_mode = int(manual_mode)
    elif mode_strategy == "per_sample":
        forced_mode = None
    else:
        raise ValueError(f"Unknown mode_strategy: {mode_strategy}")

    for hist_path in hist_files:
        hist_file = os.path.basename(hist_path)
        is_d0 = hist_path in d0_hist_set
        group = "Control" if is_d0 else "Treated"

        df = load_histogram(hist_path)

        if mode_strategy == "per_sample":
            mode, _ = detect_mode_robust(df)
            used_mode = int(mode)
        else:
            used_mode = int(forced_mode)

        instab, con, exp, left, right, thr, peak_reads = compute_instability_fixed_mode(
            df, used_mode, threshold_percent
        )

        # Derive a nice sample_name regardless of rest/per suffix
        sample_stem = hist_file
        if sample_stem.endswith(".rest.histogram"):
            sample_stem = sample_stem.replace(".rest.histogram", "")
        elif sample_stem.endswith(".per.histogram"):
            sample_stem = sample_stem.replace(".per.histogram", "")
        elif sample_stem.endswith("_thresholded.histogram"):
            sample_stem = sample_stem.replace("_thresholded.histogram", "")

        results.append(
            {
                "Sample_name": sample_stem,
                "Dataset": os.path.basename(histogram_folder),
                "Group": group,
                "Mode_Strategy": mode_strategy,
                "Forced_Mode": used_mode,
                "Instability_Index": float(f"{instab:.4f}"),
                "Contraction_Index": float(f"{con:.4f}"),
                "Expansion_Index": float(f"{exp:.4f}"),
                "Left_Boundary": int(left),
                "Right_Boundary": int(right),
                "ThresholdReads": float(f"{thr:.4f}"),
                "PeakReads_at_Mode": float(f"{peak_reads:.4f}"),
            }
        )

    out_df = pd.DataFrame(results)
    output_csv = os.path.join(output_dir, "instability_results.csv")
    out_df.to_csv(output_csv, index=False)

    if d0_mode_info is not None:
        d0_info_csv = os.path.join(output_dir, "d0_mode_info.csv")
        d0_mode_info.to_csv(d0_info_csv, index=False)

    print(f"✓ Instability results saved to: {output_csv}")

    # Optional: organize files into Control/Treated for delta tab usage
    if control_fasta_files:
        print("Organizing histogram files into Control and Treated folders...")
        organize_histograms_by_group(histogram_folder, control_fasta_files)

    return out_df


# -----------------------------------
# Enhanced raw histogram plotting (supports rest/per)
# -----------------------------------
def create_enhanced_raw_histogram_plot(histogram_path, threshold_percent=5, output_dir=None):
    try:
        if not os.path.exists(histogram_path):
            return None

        df = pd.read_table(histogram_path, skiprows=range(0, 6), names=["repeat_size", "reads"])
        sample_name = (
            os.path.basename(histogram_path)
            .replace(".rest.histogram", "")
            .replace(".per.histogram", "")
        )

        data_range = df["reads"].max() - df["reads"].min()
        prominence_threshold = 0.01 * data_range
        peaks, properties = find_peaks(df["reads"], prominence=prominence_threshold, distance=10)

        if len(peaks) == 0:
            return None

        peak_info = []
        for i, peak_idx in enumerate(peaks):
            repeat_size = df.iloc[peak_idx]["repeat_size"]
            reads = df.iloc[peak_idx]["reads"]
            prominence = properties["prominences"][i]
            peak_info.append({"repeat_size": int(repeat_size), "reads": int(reads), "prominence": prominence})

        peak_info = sorted(peak_info, key=lambda x: x["reads"], reverse=True)
        selected_peak = peak_info[0]
        selected_mode = selected_peak["repeat_size"]
        selected_peak_reads = selected_peak["reads"]

        threshold_value = (threshold_percent / 100.0) * selected_peak_reads

        df_sorted = df.sort_values("repeat_size").reset_index(drop=True)
        mode_idx = df_sorted[df_sorted["repeat_size"] == selected_mode].index[0]

        left_boundary = None
        for i in range(mode_idx - 1, -1, -1):
            if df_sorted.loc[i, "reads"] <= threshold_value:
                left_boundary = df_sorted.loc[i, "repeat_size"]
                break
        if left_boundary is None:
            left_boundary = df_sorted["repeat_size"].min()

        right_boundary = None
        for i in range(mode_idx + 1, len(df_sorted)):
            if df_sorted.loc[i, "reads"] <= threshold_value:
                right_boundary = df_sorted.loc[i, "repeat_size"]
                break
        if right_boundary is None:
            right_boundary = df_sorted["repeat_size"].max()

        filtered_df = df[(df["repeat_size"] >= left_boundary) & (df["repeat_size"] <= right_boundary)].copy()

        fig = Figure(figsize=(12, 8))
        ax = fig.add_subplot(111)

        excluded_df = df[~df["repeat_size"].isin(filtered_df["repeat_size"])]
        if not excluded_df.empty:
            ax.bar(
                excluded_df["repeat_size"],
                excluded_df["reads"],
                edgecolor="lightgray",
                facecolor="lightgray",
                alpha=0.5,
                label="Excluded Data",
            )

        ax.bar(
            filtered_df["repeat_size"],
            filtered_df["reads"],
            edgecolor="red",
            facecolor="red",
            alpha=0.8,
            label="Included Data (Used for Analysis)",
        )

        for peak in peak_info:
            if peak["repeat_size"] != selected_mode:
                ax.axvline(x=peak["repeat_size"], color="purple", linestyle=":", linewidth=1, alpha=0.5)
                ax.plot(peak["repeat_size"], peak["reads"], "v", color="purple", markersize=8, alpha=0.7)

        ax.axhline(
            y=threshold_value,
            color="orange",
            linestyle="--",
            linewidth=2,
            label=f"{threshold_percent}% Threshold: {threshold_value:.1f}",
        )
        ax.axvline(x=left_boundary, color="blue", linestyle=":", linewidth=2, label=f"Left Boundary: {left_boundary}")
        ax.axvline(x=right_boundary, color="blue", linestyle=":", linewidth=2, label=f"Right Boundary: {right_boundary}")
        ax.axvline(x=selected_mode, color="green", linestyle="-", linewidth=3, label=f"Selected Peak: {selected_mode}")
        ax.plot(selected_mode, selected_peak_reads, "v", color="green", markersize=8)

        ax.set_title(f"{sample_name}: Enhanced Instability Index Analysis\n")
        ax.set_xlabel("Repeat Size")
        ax.set_ylabel("Reads")
        ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
        ax.grid(True, alpha=0.3)

        if output_dir is None:
            output_dir = os.path.dirname(histogram_path)

        plots_dir = os.path.join(output_dir, "raw_histogram_plots")
        os.makedirs(plots_dir, exist_ok=True)

        plot_filename = f"{sample_name}_enhanced_plot.png"
        plot_path = os.path.join(plots_dir, plot_filename)
        canvas = FigureCanvasAgg(fig)
        fig.savefig(plot_path, bbox_inches="tight", dpi=300)
        plt.close(fig)

        filtered_csv_path = os.path.join(plots_dir, f"{sample_name}_filtered_data.csv")
        filtered_df.to_csv(filtered_csv_path, index=False)

        peaks_csv_path = os.path.join(plots_dir, f"{sample_name}_detected_peaks.csv")
        pd.DataFrame(peak_info).to_csv(peaks_csv_path, index=False)

        return plot_path

    except Exception as e:
        print(f"Error creating enhanced raw histogram plot: {e}")
        return None


# -----------------------------------
# Thresholded histogram generation (supports rest/per)
# -----------------------------------
def apply_threshold_filtering(df, selected_mode, selected_peak_reads, threshold_percent=5):
    threshold_value = threshold_percent / 100.0 * selected_peak_reads
    df_sorted = df.sort_values("repeat_size").reset_index(drop=True)
    mode_position = df_sorted[df_sorted["repeat_size"] == selected_mode].index

    if len(mode_position) == 0:
        left_boundary = selected_mode - 5
        right_boundary = selected_mode + 5
        filtered_df = df[(df["repeat_size"] >= left_boundary) & (df["repeat_size"] <= right_boundary)].copy()
        return filtered_df, left_boundary, right_boundary, threshold_value

    mode_idx = mode_position[0]

    left_intersection = None
    for i in range(mode_idx - 1, -1, -1):
        if df_sorted.iloc[i]["reads"] <= threshold_value:
            left_intersection = df_sorted.iloc[i]["repeat_size"]
            break

    right_intersection = None
    for i in range(mode_idx + 1, len(df_sorted)):
        if df_sorted.iloc[i]["reads"] <= threshold_value:
            right_intersection = df_sorted.iloc[i]["repeat_size"]
            break

    if left_intersection is None:
        for i in range(len(df_sorted)):
            if df_sorted.iloc[i]["reads"] > threshold_value:
                left_intersection = df_sorted.iloc[i]["repeat_size"]
                break
        if left_intersection is None:
            left_intersection = selected_mode - 5

    if right_intersection is None:
        for i in range(len(df_sorted) - 1, -1, -1):
            if df_sorted.iloc[i]["reads"] > threshold_value:
                right_intersection = df_sorted.iloc[i]["repeat_size"]
                break
        if right_intersection is None:
            right_intersection = selected_mode + 5

    left_boundary = left_intersection
    right_boundary = right_intersection

    filtered_df = df[(df["repeat_size"] >= left_boundary) & (df["repeat_size"] <= right_boundary)].copy()
    return filtered_df, left_boundary, right_boundary, threshold_value


def generate_thresholded_histogram(histogram_path, output_dir, threshold_percent=5):
    try:
        if not os.path.exists(histogram_path):
            return None

        df = pd.read_table(histogram_path, skiprows=range(0, 6), names=["repeat_size", "reads"])
        sample_name = (
            os.path.basename(histogram_path)
            .replace(".rest.histogram", "")
            .replace(".per.histogram", "")
        )

        data_range = df["reads"].max() - df["reads"].min()
        prominence_threshold = 0.01 * data_range
        peaks, properties = find_peaks(df["reads"], prominence=prominence_threshold, distance=10)

        if len(peaks) == 0:
            return None

        peak_info = []
        for i, peak_idx in enumerate(peaks):
            repeat_size = df.iloc[peak_idx]["repeat_size"]
            reads = df.iloc[peak_idx]["reads"]
            prominence = properties["prominences"][i]
            peak_info.append({"repeat_size": int(repeat_size), "reads": int(reads), "prominence": prominence})

        peak_info = sorted(peak_info, key=lambda x: x["reads"], reverse=True)
        selected_peak = peak_info[0]
        selected_mode = selected_peak["repeat_size"]
        selected_peak_reads = selected_peak["reads"]

        filtered_df, left_boundary, right_boundary, threshold_value = apply_threshold_filtering(
            df, selected_mode, selected_peak_reads, threshold_percent=threshold_percent
        )

        thresholded_path = os.path.join(output_dir, f"{sample_name}_thresholded.histogram")
        with open(thresholded_path, "w") as f:
            f.write(f"# Thresholded Histogram - {threshold_percent}% threshold\n")
            f.write(f"# Sample: {sample_name}\n")
            f.write(f"# Selected Peak: {selected_mode}\n")
            f.write(f"# Selected Peak Reads: {selected_peak_reads}\n")
            f.write(f"# {threshold_percent}% Threshold: {threshold_value:.1f}\n")
            f.write(f"# Boundaries: [{left_boundary}, {right_boundary}]\n")
            f.write("#\n")
            for _, row in filtered_df.iterrows():
                f.write(f"{int(row['repeat_size'])}\t{int(row['reads'])}\n")

        return thresholded_path

    except Exception as e:
        print(f"Error generating thresholded histogram: {e}")
        return None


# -----------------------------------
# TAB 1: RD + Organization + Instability Index
# -----------------------------------
class GenerateHistogramsPanel(wx.Panel):
    def __init__(self, parent):
        super().__init__(parent)
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(created_by_banner(self), flag=wx.EXPAND | wx.TOP | wx.BOTTOM, border=5)

        # Docker tar
        row_tar = wx.BoxSizer(wx.HORIZONTAL)
        row_tar.Add(wx.StaticText(self, label="Docker .tar file:"), 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 8)
        self.tarFileCtrl = wx.TextCtrl(self)
        row_tar.Add(self.tarFileCtrl, 1, wx.EXPAND)
        tarBrowseBtn = wx.Button(self, label="Browse", size=(80, 30))
        tarBrowseBtn.Bind(wx.EVT_BUTTON, self.on_browse_tar_file)
        row_tar.Add(tarBrowseBtn, 0, wx.LEFT, 8)
        vbox.Add(row_tar, 0, wx.EXPAND | wx.ALL, 8)

        # FASTA dir
        row_fa = wx.BoxSizer(wx.HORIZONTAL)
        row_fa.Add(wx.StaticText(self, label="FASTA Files Directory:"), 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 8)
        self.fastaDirCtrl = wx.TextCtrl(self)
        row_fa.Add(self.fastaDirCtrl, 1, wx.EXPAND)
        faBrowse = wx.Button(self, label="Browse", size=(80, 30))
        faBrowse.Bind(wx.EVT_BUTTON, self.on_browse_fasta)
        row_fa.Add(faBrowse, 0, wx.LEFT, 8)
        vbox.Add(row_fa, 0, wx.EXPAND | wx.ALL, 8)

        # Output base
        row_out = wx.BoxSizer(wx.HORIZONTAL)
        row_out.Add(wx.StaticText(self, label="Output Base Folder:"), 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 8)
        self.outDirCtrl = wx.TextCtrl(self)
        row_out.Add(self.outDirCtrl, 1, wx.EXPAND)
        outBrowse = wx.Button(self, label="Browse", size=(80, 30))
        outBrowse.Bind(wx.EVT_BUTTON, self.on_browse_out)
        row_out.Add(outBrowse, 0, wx.LEFT, 8)
        vbox.Add(row_out, 0, wx.EXPAND | wx.ALL, 8)

        # D0 FASTA selection (for fixed-from-D0 mode)
        row_d0 = wx.BoxSizer(wx.HORIZONTAL)
        row_d0.Add(
            wx.StaticText(self, label="D0/Control FASTA Files (multiple allowed):"),
            0,
            wx.RIGHT | wx.ALIGN_CENTER_VERTICAL,
            8,
        )
        self.d0FilesCtrl = wx.TextCtrl(self)
        row_d0.Add(self.d0FilesCtrl, 1, wx.EXPAND)
        d0BrowseBtn = wx.Button(self, label="Browse", size=(80, 30))
        d0BrowseBtn.Bind(wx.EVT_BUTTON, self.on_browse_multiple_d0_files)
        row_d0.Add(d0BrowseBtn, 0, wx.LEFT, 8)
        vbox.Add(row_d0, 0, wx.EXPAND | wx.ALL, 8)

        # Options row
        options_row = wx.BoxSizer(wx.HORIZONTAL)

        profile_box = wx.StaticBox(self, label="Profile Mode")
        profile_sizer = wx.StaticBoxSizer(profile_box, wx.HORIZONTAL)
        self.restrictive_rb = wx.RadioButton(self, label="Restrictive", style=wx.RB_GROUP)
        self.permissive_rb = wx.RadioButton(self, label="Permissive")
        self.restrictive_rb.SetValue(True)
        profile_sizer.Add(self.restrictive_rb, 0, wx.ALL, 5)
        profile_sizer.Add(self.permissive_rb, 0, wx.ALL, 5)
        options_row.Add(profile_sizer, 0, wx.RIGHT, 15)

        # Instability enabled
        self.instability_cb = wx.CheckBox(self, label="Run Instability Index")
        self.instability_cb.SetValue(True)
        options_row.Add(self.instability_cb, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 10)

        # Instability mode selector
        inst_box = wx.StaticBox(self, label="Instability mode")
        inst_sizer = wx.StaticBoxSizer(inst_box, wx.VERTICAL)

        self.inst_mode_fixed_rb = wx.RadioButton(self, label="Fixed mode from D0", style=wx.RB_GROUP)
        self.inst_mode_normal_rb = wx.RadioButton(self, label="Normal (mode per sample)")
        self.inst_mode_manual_rb = wx.RadioButton(self, label="Manual forced mode")
        self.inst_mode_fixed_rb.SetValue(True)

        inst_sizer.Add(self.inst_mode_fixed_rb, 0, wx.ALL, 2)
        inst_sizer.Add(self.inst_mode_normal_rb, 0, wx.ALL, 2)

        manual_row = wx.BoxSizer(wx.HORIZONTAL)
        manual_row.Add(self.inst_mode_manual_rb, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
        manual_row.Add(wx.StaticText(self, label="Mode:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        self.manual_mode_spin = wx.SpinCtrl(self, value="0", min=0, max=5000, initial=0, size=(90, -1))
        manual_row.Add(self.manual_mode_spin, 0)
        inst_sizer.Add(manual_row, 0, wx.ALL, 2)

        options_row.Add(inst_sizer, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 10)

        self.generate_thresholded_cb = wx.CheckBox(self, label="Generate thresholded histograms")
        self.generate_thresholded_cb.SetValue(True)
        options_row.Add(self.generate_thresholded_cb, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 10)

        self.enhanced_plots_cb = wx.CheckBox(self, label="Generate enhanced raw histogram plots")
        self.enhanced_plots_cb.SetValue(True)
        options_row.Add(self.enhanced_plots_cb, 0, wx.ALIGN_CENTER_VERTICAL)

        vbox.Add(options_row, 0, wx.ALL, 8)

        settings_row = wx.BoxSizer(wx.HORIZONTAL)
        settings_row.Add(wx.StaticText(self, label="Repeat Size for Histogram:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        self.repeat_size_spin = wx.SpinCtrl(self, value="250", min=50, max=1000, initial=250, size=(80, -1))
        settings_row.Add(self.repeat_size_spin, 0, wx.RIGHT, 15)

        settings_row.Add(wx.StaticText(self, label="Threshold % for filtering:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        self.threshold_spin = wx.SpinCtrlDouble(self, value="5.0", min=0.1, max=50.0, inc=0.1, size=(80, -1))
        settings_row.Add(self.threshold_spin, 0, wx.RIGHT, 10)
        vbox.Add(settings_row, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 8)

        self.run_btn = wx.Button(self, label="Run RD Program", size=(160, 38))
        self.run_btn.Bind(wx.EVT_BUTTON, self.on_run)
        vbox.Add(self.run_btn, 0, wx.ALL, 8)

        self.gauge = wx.Gauge(self, range=100)
        vbox.Add(self.gauge, 0, wx.EXPAND | wx.ALL, 8)

        vbox.Add(wx.StaticText(self, label="Run log:"), 0, wx.LEFT | wx.TOP, 8)
        self.outputCtrl = wx.TextCtrl(self, style=wx.TE_MULTILINE | wx.TE_READONLY)
        vbox.Add(self.outputCtrl, 1, wx.EXPAND | wx.ALL, 8)

        self.SetSizer(vbox)

    def on_browse_tar_file(self, _):
        with wx.FileDialog(
            self,
            "Choose Docker .tar file:",
            wildcard="Tar files (*.tar)|*.tar",
            style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST,
        ) as d:
            if d.ShowModal() == wx.ID_OK:
                self.tarFileCtrl.SetValue(d.GetPath())

    def on_browse_fasta(self, _):
        with wx.DirDialog(self, "Choose the FASTA files directory:") as d:
            if d.ShowModal() == wx.ID_OK:
                self.fastaDirCtrl.SetValue(d.GetPath())

    def on_browse_out(self, _):
        with wx.DirDialog(self, "Choose the output base directory:") as d:
            if d.ShowModal() == wx.ID_OK:
                self.outDirCtrl.SetValue(d.GetPath())

    def on_browse_multiple_d0_files(self, _):
        with wx.FileDialog(
            self,
            "Choose D0/Control FASTA files (multiple allowed):",
            wildcard="FASTA files (*.fasta)|*.fasta|All files (*.*)|*.*",
            style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST | wx.FD_MULTIPLE,
        ) as d:
            if d.ShowModal() == wx.ID_OK:
                selected_files = d.GetPaths()
                self.d0FilesCtrl.SetValue(";".join(selected_files))

    def log(self, msg):
        wx.CallAfter(self.outputCtrl.AppendText, msg + "\n")

    def on_run(self, _):
        tar_file = self.tarFileCtrl.GetValue().strip()
        fasta_dir = self.fastaDirCtrl.GetValue().strip()
        output_dir = self.outDirCtrl.GetValue().strip()
        d0_files_str = self.d0FilesCtrl.GetValue().strip()

        profile = "restrictive" if self.restrictive_rb.GetValue() else "permissive"
        generate_thresholded = self.generate_thresholded_cb.GetValue()
        run_instability = self.instability_cb.GetValue()
        generate_enhanced_plots = self.enhanced_plots_cb.GetValue()
        threshold_percent = self.threshold_spin.GetValue()
        repeat_size = self.repeat_size_spin.GetValue()

        # choose mode strategy
        mode_strategy = "fixed_from_d0"
        manual_mode = None
        if self.inst_mode_normal_rb.GetValue():
            mode_strategy = "per_sample"
        elif self.inst_mode_manual_rb.GetValue():
            mode_strategy = "manual"
            manual_mode = int(self.manual_mode_spin.GetValue())

        if not fasta_dir or not os.path.isdir(fasta_dir):
            wx.MessageBox("Please select a valid FASTA files directory.", "Error", wx.OK | wx.ICON_ERROR)
            return

        if not output_dir:
            output_dir = fasta_dir

        d0_files = []
        if d0_files_str:
            d0_files = [f.strip() for f in d0_files_str.split(";") if f.strip()]

        if run_instability and mode_strategy == "fixed_from_d0" and not d0_files:
            wx.MessageBox(
                "You selected 'Fixed mode from D0' for instability. Please select at least one D0/Control FASTA file.",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            return

        if run_instability and mode_strategy == "manual" and (manual_mode is None or manual_mode <= 0):
            wx.MessageBox(
                "You selected 'Manual forced mode'. Please enter a positive mode value.",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            return

        self.run_btn.Enable(False)
        self.gauge.SetValue(0)

        thread = threading.Thread(
            target=self._run_all,
            args=(
                tar_file,
                fasta_dir,
                output_dir,
                d0_files,
                profile,
                generate_thresholded,
                run_instability,
                generate_enhanced_plots,
                threshold_percent,
                repeat_size,
                mode_strategy,
                manual_mode,
            ),
            daemon=True,
        )
        thread.start()

    def _run_all(
        self,
        tar_file,
        fasta_dir,
        output_dir,
        d0_files,
        profile,
        generate_thresholded,
        run_instability,
        generate_enhanced_plots,
        threshold_percent,
        repeat_size,
        mode_strategy,
        manual_mode,
    ):
        try:
            self.log(f"Starting RD with profile={profile}, repeat_range=[0:{repeat_size}]")
            wx.CallAfter(self.gauge.SetValue, 5)

            hist_dir = os.path.join(output_dir, "histograms")
            os.makedirs(hist_dir, exist_ok=True)

            # Load Docker image from .tar if provided
            if tar_file and os.path.isfile(tar_file) and tar_file.endswith(".tar"):
                self.log(f"Loading Docker image from {tar_file}...")

                # Use cross-platform Docker command
                docker_cmd = "docker" if platform.system() == "Windows" else "/usr/local/bin/docker"
                load_cmd = [docker_cmd, "load", "-i", tar_file]

                result = subprocess.run(load_cmd, capture_output=True, text=True)
                if result.returncode != 0:
                    self.log(f"❌ Error loading Docker image: {result.stderr}")
                    wx.CallAfter(self._processing_complete, False, None)
                    return
                self.log("Docker image loaded successfully")
            elif tar_file:
                self.log(f"Warning: Invalid tar file path: {tar_file}")
                self.log("Continuing without Docker image...")

            wx.CallAfter(self.gauge.SetValue, 15)

            all_fasta_files = [os.path.join(fasta_dir, f) for f in os.listdir(fasta_dir) if f.endswith(".fasta")]
            if not all_fasta_files:
                self.log("No FASTA files found.")
                wx.CallAfter(self._processing_complete, False, None)
                return

            prf_file = (
                "/app/RepeatDetector/Profiles/CAG/Annex10_cag_correctedFreq_notlog_AND_Complete.prf"
                if profile == "restrictive"
                else "/app/RepeatDetector/Profiles/CAG/Annex2_cag.prf"
            )
            output_suffix = ".rest" if profile == "restrictive" else ".per"
            hist_ext = hist_suffix_for_profile(profile)  # .rest.histogram or .per.histogram

            total = len(all_fasta_files)
            processed_histograms = []

            for i, fasta_file in enumerate(all_fasta_files, 1):
                base = os.path.splitext(os.path.basename(fasta_file))[0]
                out_name = f"{base}{output_suffix}"
                self.log(f"[{i}/{total}] {os.path.basename(fasta_file)}")

                docker_exec = "docker" if platform.system() == "Windows" else "/usr/local/bin/docker"
                docker_cmd = (
                    f"""{docker_exec} run --rm -v "{fasta_dir}:/mnt/fasta" """
                    f"""-e LD_LIBRARY_PATH=/app/RepeatDetector/build/external/htslib/src/htslib """
                    f"""repeat-detector /usr/local/bin/RepeatDetecter --prf {prf_file} /mnt/fasta/{os.path.basename(fasta_file)} """
                    f"""--output-name "/mnt/fasta/{out_name}" -o histogram --with-revcomp """
                    f"""--cycle-range [0:{repeat_size}] --verbose"""
                )
                result = subprocess.run(docker_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                if result.returncode != 0:
                    self.log(f"❌ Error: {result.stderr.decode('utf-8')}")
                else:
                    # sample.rest.histogram or sample.per.histogram
                    hist_path = os.path.join(fasta_dir, f"{out_name}.histogram")
                    if os.path.exists(hist_path):
                        processed_histograms.append(hist_path)

                wx.CallAfter(self.gauge.SetValue, 15 + int(55 * (i / total)))

            self.log("✔ RD histogram generation complete.")
            wx.CallAfter(self.gauge.SetValue, 75)

            self.log("Copying histograms to output folder...")
            for hist_path in processed_histograms:
                if os.path.exists(hist_path):
                    shutil.copy2(hist_path, hist_dir)

            self.log("Cleaning up temporary histogram files from FASTA directory...")
            cleanup_rest_histograms_from_fasta(fasta_dir)

            # Enhanced plots
            if generate_enhanced_plots:
                self.log("Generating enhanced raw histogram plots...")
                wx.CallAfter(self.gauge.SetValue, 80)

                raw_plots_dir = os.path.join(output_dir, "raw_histogram_plots")
                os.makedirs(raw_plots_dir, exist_ok=True)

                # Use hist_ext (rest OR per)
                hist_files = [os.path.join(hist_dir, f) for f in os.listdir(hist_dir) if f.endswith(hist_ext)]
                self.log(f"Found {len(hist_files)} {hist_ext} files for enhanced plotting")

                for i, histogram_path in enumerate(hist_files):
                    try:
                        plot_path = create_enhanced_raw_histogram_plot(
                            histogram_path, threshold_percent=threshold_percent, output_dir=raw_plots_dir
                        )
                        if plot_path:
                            sample_name = (
                                os.path.basename(histogram_path)
                                .replace(".rest.histogram", "")
                                .replace(".per.histogram", "")
                            )
                            self.log(f"✓ Generated enhanced plot for {sample_name}")
                    except Exception as e:
                        self.log(f"❌ Failed to generate enhanced plot: {e}")

                    if hist_files:
                        wx.CallAfter(self.gauge.SetValue, 80 + int(5 * (i / len(hist_files))))

            # Thresholded
            if generate_thresholded:
                self.log("Generating thresholded histograms...")
                wx.CallAfter(self.gauge.SetValue, 85)

                # Use hist_ext (rest OR per)
                hist_files = [os.path.join(hist_dir, f) for f in os.listdir(hist_dir) if f.endswith(hist_ext)]
                self.log(f"Found {len(hist_files)} {hist_ext} files for thresholding")

                for i, histogram_path in enumerate(hist_files):
                    try:
                        generate_thresholded_histogram(histogram_path, hist_dir, threshold_percent)
                        sample_name = (
                            os.path.basename(histogram_path)
                            .replace(".rest.histogram", "")
                            .replace(".per.histogram", "")
                        )
                        self.log(f"✓ Generated thresholded histogram for {sample_name}")
                    except Exception as e:
                        self.log(f"❌ Failed to generate thresholded histogram: {e}")

                    if hist_files:
                        wx.CallAfter(self.gauge.SetValue, 85 + int(5 * (i / len(hist_files))))

            # Instability (3 modes)
            if run_instability:
                self.log(f"Running Instability Index (mode_strategy={mode_strategy})...")
                wx.CallAfter(self.gauge.SetValue, 90)

                # Find D0 histogram files in hist_dir (matching current profile ext)
                d0_hist_files = []
                if d0_files:
                    for d0_fasta in d0_files:
                        d0_base = os.path.splitext(os.path.basename(d0_fasta))[0]
                        for hist_file in os.listdir(hist_dir):
                            # Match rest or per depending on chosen profile
                            if hist_file.endswith(hist_ext) and d0_base in hist_file:
                                d0_hist_files.append(os.path.join(hist_dir, hist_file))
                                self.log(f"✓ Found D0 histogram: {hist_file}")

                if mode_strategy == "fixed_from_d0" and not d0_hist_files:
                    self.log("❌ No D0 histogram files found. Cannot run fixed-from-D0 instability.")
                else:
                    run_instability_analysis(
                        histogram_folder=hist_dir,
                        output_dir=output_dir,
                        threshold_percent=threshold_percent,
                        d0_hist_files=d0_hist_files,
                        mode_strategy=mode_strategy,
                        manual_mode=manual_mode,
                        control_fasta_files=d0_files if d0_files else None,
                        hist_ext=hist_ext,
                    )
                    self.log(f"✓ Instability analysis complete (strategy={mode_strategy}). Results in: {output_dir}")

            self.log("Final cleanup of output folder...")
            cleanup_rest_histograms_from_output(output_dir)

            wx.CallAfter(self.gauge.SetValue, 100)
            self.log("All processing complete.")
            wx.CallAfter(self._processing_complete, True, None)

        except Exception as e:
            self.log(f"❌ Error: {e}")
            import traceback

            self.log(f"Traceback: {traceback.format_exc()}")
            wx.CallAfter(self._processing_complete, False, None)

    def _processing_complete(self, success, ds_map):
        self.run_btn.Enable(True)
        self.gauge.SetValue(100 if success else 0)


# -------------------------------------------------------------------------
# TAB 2: Delta plot panel
# -------------------------------------------------------------------------
def normalize_by_column_sum(df, start_col=1):
    df_norm = df.copy()
    for col in df_norm.columns[start_col:]:
        col_sum = df_norm[col].sum()
        df_norm[col] = df_norm[col] / col_sum if col_sum > 0 else 0
    return df_norm


def compute_control_treatment_stats(df, control_cols, treatment_cols):
    df_stats = df.copy()
    df_stats["control_avg"] = df_stats[control_cols].mean(axis=1)
    if len(control_cols) > 1:
        df_stats["control_sd"] = df_stats[control_cols].std(axis=1)
        df_stats["control_sem"] = df_stats["control_sd"] / np.sqrt(len(control_cols))
    else:
        df_stats["control_sd"] = 0
        df_stats["control_sem"] = 0

    df_stats["treatment_avg"] = df_stats[treatment_cols].mean(axis=1)
    if len(treatment_cols) > 1:
        df_stats["treatment_sd"] = df_stats[treatment_cols].std(axis=1)
        df_stats["treatment_sem"] = df_stats["treatment_sd"] / np.sqrt(len(treatment_cols))
    else:
        df_stats["treatment_sd"] = 0
        df_stats["treatment_sem"] = 0

    df_stats["Delta(Treated-Control)"] = df_stats["treatment_avg"] - df_stats["control_avg"]
    return df_stats


def align_to_global_mode(df):
    bin_col = df.columns[0]
    sample_cols = df.columns[1:]

    control_cols = [col for col in sample_cols if any(x in col.upper() for x in ["D0", "DAY0", "CONTROL", "CTRL"])]
    if not control_cols:
        control_cols = sample_cols[:1]

    control_avg = df[control_cols].mean(axis=1)
    global_mode_idx = control_avg.idxmax()
    global_mode = df.iloc[global_mode_idx][bin_col]
    df["bin_offset"] = df[bin_col] - global_mode
    return df, global_mode


def bin_delta_by_offset_with_auc(df, bin_size=1):
    if "bin_offset" not in df.columns or "Delta(Treated-Control)" not in df.columns:
        raise ValueError("DataFrame must have 'bin_offset' and 'Delta(Treated-Control)' columns")

    df_binned = df.copy()
    df_binned["bin"] = (df_binned["bin_offset"] // bin_size) * bin_size
    grouped = (
        df_binned.groupby("bin")
        .agg({"Delta(Treated-Control)": ["sum", "std", "count"], "bin_offset": "first"})
        .reset_index()
    )
    grouped.columns = ["bin", "Delta_sum", "Delta_std", "n", "Bin_Offset"]
    grouped["AUC"] = grouped["Delta_sum"].apply(lambda x: max(x, 0))
    grouped["Delta_sum_se"] = np.sqrt(grouped["n"]) * grouped["Delta_std"]
    grouped["Bin_Center"] = grouped["Bin_Offset"] + bin_size / 2
    return grouped[["Bin_Offset", "Bin_Center", "Delta_sum", "Delta_std", "Delta_sum_se", "n", "AUC"]]


def resolve_uploaded_paths(paths, include_thresholded=False):
    resolved = []
    for path in paths:
        if os.path.isdir(path):
            files = (
                glob.glob(os.path.join(path, "*.histogram"))
                + glob.glob(os.path.join(path, "*.csv"))
                + glob.glob(os.path.join(path, "*.tsv"))
                + glob.glob(os.path.join(path, "*.txt"))
            )
            if include_thresholded:
                thresholded_files = [f for f in files if f.endswith("_thresholded.histogram")]
                resolved.extend(thresholded_files if thresholded_files else files)
            else:
                regular_files = [f for f in files if not f.endswith("_thresholded.histogram")]
                resolved.extend(regular_files)
        else:
            resolved.append(path)
    return resolved


def read_histogram_files_with_folder_prefix(file_paths, folder_path, group_type, use_thresholded=False):
    merged_df = None
    folder_name = os.path.basename(os.path.normpath(folder_path))

    for path in file_paths:
        ext = os.path.splitext(path)[1].lower()
        is_hist = (ext == ".histogram")
        is_thresholded = path.endswith("_thresholded.histogram")
        delimiter = "\t" if ext in (".tsv", ".histogram") else ","

        if is_thresholded:
            skip = 7
        elif is_hist:
            skip = 6
        else:
            skip = 0

        try:
            df = pd.read_csv(path, header=None, usecols=[0, 1], sep=delimiter, skiprows=skip, engine="python")
        except Exception:
            try:
                df = pd.read_csv(path, header=None, usecols=[0, 1], sep=None, skiprows=skip, engine="python")
            except Exception as e2:
                print(f"Could not read {path}: {e2}")
                continue

        df = df.dropna()
        df.iloc[:, 1] = pd.to_numeric(df.iloc[:, 1], errors="coerce")
        df = df.dropna()
        if df.shape[1] < 2 or df.empty:
            continue

        base_name = os.path.basename(path)

        # Strip ".rest" / ".per" tags from the stem for nicer labels
        if is_thresholded:
            stem = base_name.replace("_thresholded.histogram", "")
        else:
            stem = os.path.splitext(base_name)[0]  # includes ".rest" or ".per" in the stem if present

        stem = stem.replace(".rest", "").replace(".per", "")
        clean_name = stem

        if is_thresholded:
            label = f"{folder_name}_{clean_name}_thresholded"
        else:
            label = f"{folder_name}_{clean_name}"

        series = pd.Series(df.iloc[:, 1].values, index=pd.to_numeric(df.iloc[:, 0].values, errors="coerce"))
        series = series.dropna()
        series.index.name = "repeat length"
        series.name = label

        if merged_df is None:
            merged_df = series.to_frame()
        else:
            merged_df = merged_df.join(series, how="outer")

    if merged_df is None or merged_df.empty:
        raise ValueError(f"No valid files found in {folder_path}")

    merged_df = merged_df.fillna(0).reset_index()
    return merged_df


def weighted_mean_and_sd(x, w):
    if np.sum(w) == 0:
        return 0, 0
    mean = np.sum(w * x) / np.sum(w)
    variance = np.sum(w * (x - mean) ** 2) / np.sum(w)
    return mean, np.sqrt(variance)


def analyze_auc_stats(binned_data):
    if binned_data.empty or "AUC" not in binned_data.columns:
        return {}
    df = binned_data.copy()
    df = df[(df["AUC"] > 0) & df["Bin_Offset"].notna()]
    results = {}

    exp = df[df["Bin_Offset"] > 0]
    if len(exp) > 0:
        mean_exp, sd_exp = weighted_mean_and_sd(exp["Bin_Offset"].values, exp["AUC"].values)
        results["expansion_mean"] = float(f"{mean_exp:.4f}")
        results["expansion_sd"] = float(f"{sd_exp:.4f}")
        results["expansion_AUC"] = float(f"{exp['AUC'].sum():.4f}")

    con = df[df["Bin_Offset"] < 0]
    if len(con) > 0:
        mean_con, sd_con = weighted_mean_and_sd(con["Bin_Offset"].values, con["AUC"].values)
        results["contraction_mean"] = float(f"{mean_con:.4f}")
        results["contraction_sd"] = float(f"{sd_con:.4f}")
        results["contraction_AUC"] = float(f"{con['AUC'].sum():.4f}")

    overall = df[df["Bin_Offset"] != 0]
    if len(overall) > 0:
        mean_all, sd_all = weighted_mean_and_sd(overall["Bin_Offset"].values, overall["AUC"].values)
        results["overall_mean"] = float(f"{mean_all:.4f}")
        results["overall_sd"] = float(f"{sd_all:.4f}")

    if "expansion_AUC" in results and "contraction_AUC" in results:
        if results["contraction_AUC"] > 0:
            results["bias_ratio"] = float(f"{results['expansion_AUC'] / results['contraction_AUC']:.4f}")
        else:
            results["bias_ratio"] = float("inf") if results["expansion_AUC"] > 0 else 0.0

    return results


def calculate_delta_stats(binned_data):
    if binned_data.empty:
        return {}
    return {
        "mean_delta": float(f"{binned_data['Delta_sum'].mean():.4f}"),
        "std_delta": float(f"{binned_data['Delta_sum'].std():.4f}"),
        "max_delta": float(f"{binned_data['Delta_sum'].max():.4f}"),
        "min_delta": float(f"{binned_data['Delta_sum'].min():.4f}"),
        "total_bins": len(binned_data),
        "positive_bins": len(binned_data[binned_data["Delta_sum"] > 0]),
        "negative_bins": len(binned_data[binned_data["Delta_sum"] < 0]),
    }


def generate_delta_plot_png(
    binned_data1,
    binned_data2=None,
    output_path="delta_plot.png",
    title="Delta Sum Plot",
    color1="red",
    color2="blue",
    label1="Dataset 1",
    label2="Dataset 2",
    xlim_left=-100,
    xlim_right=100,
    ylim=None,
    bin_size=1,
    show_d0_line=True,
    ylabel="Frequency(D42-D0)",
):
    plt.figure(figsize=(8, 5))
    plt.plot(binned_data1["Bin_Offset"], binned_data1["Delta_sum"], color=color1, linewidth=2, label=label1)

    if binned_data2 is not None:
        plt.plot(binned_data2["Bin_Offset"], binned_data2["Delta_sum"], color=color2, linewidth=2, label=label2)

    if show_d0_line:
        plt.axvline(x=0, color="black", linestyle="--", linewidth=1.5, alpha=0.7)

    plt.axhline(y=0.0, color="black", linewidth=1.5, alpha=0.7)

    plt.xlim(xlim_left, xlim_right)
    if ylim:
        plt.ylim(ylim[0], ylim[1])

    plt.xlabel("Change in repeat size", fontsize=11)
    plt.ylabel(ylabel, fontsize=11)
    plt.title(title, fontsize=12, pad=10)
    plt.legend(loc="best", fontsize=10, frameon=True)
    plt.grid(False)
    plt.gcf().set_facecolor("white")
    plt.gca().set_facecolor("white")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close()


def clean_label(s):
    s = (s or "").strip()
    if not s:
        return "Dataset"
    s = s.replace(" ", "_")
    return "".join(ch if ch.isalnum() or ch in ["_", "-"] else "_" for ch in s)


def align_tailtrim_and_ks_on_offsets(
    ds1_binned_path,
    ds2_binned_path,
    out_dir,
    name1="Dataset1",
    name2="Dataset2",
    tol=0.0,
    out_aligned_csv="realigned_zero_filled_deltas.csv",
    out_ks_csv="ks_results.csv",
):
    os.makedirs(out_dir, exist_ok=True)

    col1 = f"Delta_sum_{clean_label(name1)}"
    col2 = f"Delta_sum_{clean_label(name2)}"

    ds1 = pd.read_csv(ds1_binned_path)
    ds2 = pd.read_csv(ds2_binned_path)

    off = "Bin_Offset"
    dcol = "Delta_sum"

    p1 = ds1[[off, dcol]].copy()
    p2 = ds2[[off, dcol]].copy()
    p1.columns = ["Bin_Offset", col1]
    p2.columns = ["Bin_Offset", col2]

    p1["Bin_Offset"] = pd.to_numeric(p1["Bin_Offset"], errors="coerce")
    p2["Bin_Offset"] = pd.to_numeric(p2["Bin_Offset"], errors="coerce")
    p1[col1] = pd.to_numeric(p1[col1], errors="coerce")
    p2[col2] = pd.to_numeric(p2[col2], errors="coerce")

    aligned = pd.merge(p1, p2, on="Bin_Offset", how="outer").sort_values("Bin_Offset").reset_index(drop=True)
    aligned[col1] = aligned[col1].fillna(0.0)
    aligned[col2] = aligned[col2].fillna(0.0)

    nonzero_mask = (aligned[col1].abs() > tol) | (aligned[col2].abs() > tol)
    if nonzero_mask.any():
        last_idx = nonzero_mask[nonzero_mask].index[-1]
        aligned = aligned.loc[:last_idx].reset_index(drop=True)
    else:
        aligned = aligned.iloc[0:0].copy()

    aligned_path = os.path.join(out_dir, out_aligned_csv)
    aligned.to_csv(aligned_path, index=False)

    delta1 = aligned[col1].to_numpy()
    delta2 = aligned[col2].to_numpy()

    if len(delta1) == 0 or len(delta2) == 0:
        ks_stat, ks_p, adj_p = np.nan, np.nan, np.nan
        sig = False
        max_off = None
        n_bins = 0
    else:
        ks_stat, ks_p = ks_2samp(delta1, delta2)
        adj_p = multipletests([ks_p], method="fdr_bh")[1][0]
        sig = bool(adj_p < 0.05)
        max_off = int(aligned["Bin_Offset"].max()) if len(aligned) else None
        n_bins = int(len(aligned))

    ks_results = pd.DataFrame(
        [
            {
                "dataset1_name": name1,
                "dataset2_name": name2,
                "delta_col_1": col1,
                "delta_col_2": col2,
                "n_bins_after_tail_trim": n_bins,
                "max_bin_offset_after_tail_trim": max_off,
                "ks_statistic": (float(ks_stat) if np.isfinite(ks_stat) else np.nan),
                "raw_pvalue": (float(ks_p) if np.isfinite(ks_p) else np.nan),
                "adjusted_pvalue_fdr_bh": (float(adj_p) if np.isfinite(adj_p) else np.nan),
                "significant_at_0p05": sig,
                "tol_used_for_nonzero": float(tol),
                "aligned_csv": os.path.basename(aligned_path),
            }
        ]
    )

    ks_path = os.path.join(out_dir, out_ks_csv)
    ks_results.to_csv(ks_path, index=False)
    return aligned_path, ks_path, ks_results.iloc[0].to_dict()


# -------------------------------------------------------------------------
# TAB 2: Delta GUI panel
# -------------------------------------------------------------------------
class RepeatDetectorPanel(wx.Panel):
    def __init__(self, parent):
        super().__init__(parent)
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(created_by_banner(self), flag=wx.EXPAND | wx.TOP | wx.BOTTOM, border=5)

        ds1_box = wx.StaticBox(self, label="Dataset 1 (User organizes D0 and D42/different days)")
        ds1_sizer = wx.StaticBoxSizer(ds1_box, wx.VERTICAL)

        r1 = wx.BoxSizer(wx.HORIZONTAL)
        r1.Add(wx.StaticText(self, label="Control Folder (D0):"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        self.control_folder_ctrl1 = wx.TextCtrl(self)
        r1.Add(self.control_folder_ctrl1, 1, wx.RIGHT, 5)
        b1 = wx.Button(self, label="Browse")
        b1.Bind(wx.EVT_BUTTON, lambda evt: self._browse_dir_into(self.control_folder_ctrl1))
        r1.Add(b1, 0)
        ds1_sizer.Add(r1, 0, wx.EXPAND | wx.ALL, 5)

        r2 = wx.BoxSizer(wx.HORIZONTAL)
        r2.Add(wx.StaticText(self, label="Treated Folder (D42 or other treatments):"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        self.treatment_folder_ctrl1 = wx.TextCtrl(self)
        r2.Add(self.treatment_folder_ctrl1, 1, wx.RIGHT, 5)
        b2 = wx.Button(self, label="Browse")
        b2.Bind(wx.EVT_BUTTON, lambda evt: self._browse_dir_into(self.treatment_folder_ctrl1))
        r2.Add(b2, 0)
        ds1_sizer.Add(r2, 0, wx.EXPAND | wx.ALL, 5)

        r1_label = wx.BoxSizer(wx.HORIZONTAL)
        r1_label.Add(wx.StaticText(self, label="Legend Label for Dataset 1:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        self.legend_label1_ctrl = wx.TextCtrl(self, value="Dataset 1", size=(200, -1))
        r1_label.Add(self.legend_label1_ctrl, 0, wx.RIGHT, 5)
        ds1_sizer.Add(r1_label, 0, wx.EXPAND | wx.ALL, 5)

        vbox.Add(ds1_sizer, 0, wx.EXPAND | wx.ALL, 10)

        ds2_box = wx.StaticBox(self, label="Dataset 2 (optional overlay)")
        ds2_sizer = wx.StaticBoxSizer(ds2_box, wx.VERTICAL)

        self.enable_dataset2_cb = wx.CheckBox(self, label="Enable Dataset 2 overlay")
        ds2_sizer.Add(self.enable_dataset2_cb, 0, wx.ALL, 5)

        r3 = wx.BoxSizer(wx.HORIZONTAL)
        r3.Add(wx.StaticText(self, label="Control Folder (D0):"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        self.control_folder_ctrl2 = wx.TextCtrl(self)
        r3.Add(self.control_folder_ctrl2, 1, wx.RIGHT, 5)
        b3 = wx.Button(self, label="Browse")
        b3.Bind(wx.EVT_BUTTON, lambda evt: self._browse_dir_into(self.control_folder_ctrl2))
        r3.Add(b3, 0)
        ds2_sizer.Add(r3, 0, wx.EXPAND | wx.ALL, 5)

        r4 = wx.BoxSizer(wx.HORIZONTAL)
        r4.Add(wx.StaticText(self, label="Treated Folder (D42 or other treatments):"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        self.treatment_folder_ctrl2 = wx.TextCtrl(self)
        r4.Add(self.treatment_folder_ctrl2, 1, wx.RIGHT, 5)
        b4 = wx.Button(self, label="Browse")
        b4.Bind(wx.EVT_BUTTON, lambda evt: self._browse_dir_into(self.treatment_folder_ctrl2))
        r4.Add(b4, 0)
        ds2_sizer.Add(r4, 0, wx.EXPAND | wx.ALL, 5)

        r2_label = wx.BoxSizer(wx.HORIZONTAL)
        r2_label.Add(wx.StaticText(self, label="Legend Label for Dataset 2:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        self.legend_label2_ctrl = wx.TextCtrl(self, value="Dataset 2", size=(200, -1))
        r2_label.Add(self.legend_label2_ctrl, 0, wx.RIGHT, 5)
        ds2_sizer.Add(r2_label, 0, wx.EXPAND | wx.ALL, 5)

        vbox.Add(ds2_sizer, 0, wx.EXPAND | wx.ALL, 10)

        plot_box = wx.StaticBox(self, label="Plot Settings (Bin Offset)")
        plot_sizer = wx.StaticBoxSizer(plot_box, wx.VERTICAL)

        settings_row1 = wx.BoxSizer(wx.HORIZONTAL)
        settings_row1.Add(wx.StaticText(self, label="X-axis left limit:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        self.xlim_left_ctrl = wx.SpinCtrl(self, value="-100", min=-200, max=0, initial=-100, size=(80, -1))
        settings_row1.Add(self.xlim_left_ctrl, 0, wx.RIGHT, 15)

        settings_row1.Add(wx.StaticText(self, label="X-axis right limit:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        self.xlim_right_ctrl = wx.SpinCtrl(self, value="100", min=0, max=200, initial=100, size=(80, -1))
        settings_row1.Add(self.xlim_right_ctrl, 0, wx.RIGHT, 15)

        settings_row1.Add(wx.StaticText(self, label="Bin size:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        self.bin_size_ctrl = wx.SpinCtrl(self, value="1", min=1, max=20, initial=1, size=(60, -1))
        settings_row1.Add(self.bin_size_ctrl, 0)
        plot_sizer.Add(settings_row1, 0, wx.EXPAND | wx.ALL, 5)

        colors_row = wx.BoxSizer(wx.HORIZONTAL)
        colors_row.Add(wx.StaticText(self, label="Dataset 1 color:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        self.color1_choice = wx.Choice(self, choices=["red", "blue", "green", "purple", "orange"])
        self.color1_choice.SetSelection(0)
        colors_row.Add(self.color1_choice, 0, wx.RIGHT, 15)

        colors_row.Add(wx.StaticText(self, label="Dataset 2 color:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        self.color2_choice = wx.Choice(self, choices=["blue", "red", "green", "purple", "orange"])
        self.color2_choice.SetSelection(0)
        colors_row.Add(self.color2_choice, 0)
        plot_sizer.Add(colors_row, 0, wx.EXPAND | wx.ALL, 5)

        settings_row2 = wx.BoxSizer(wx.HORIZONTAL)
        settings_row2.Add(wx.StaticText(self, label="Y-axis min:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        self.ymin_ctrl = wx.TextCtrl(self, value="", size=(60, -1))
        settings_row2.Add(self.ymin_ctrl, 0, wx.RIGHT, 5)

        settings_row2.Add(wx.StaticText(self, label="Y-axis max:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        self.ymax_ctrl = wx.TextCtrl(self, value="", size=(60, -1))
        settings_row2.Add(self.ymax_ctrl, 0, wx.RIGHT, 15)

        settings_row2.Add(wx.StaticText(self, label="Y-axis label:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        self.ylabel_ctrl = wx.TextCtrl(self, value="Frequency(D42-D0)", size=(160, -1))
        settings_row2.Add(self.ylabel_ctrl, 0)

        plot_sizer.Add(settings_row2, 0, wx.EXPAND | wx.ALL, 5)

        self.show_d0_cb = wx.CheckBox(self, label="Show D0 vertical line (bin_offset = 0)")
        self.show_d0_cb.SetValue(True)
        plot_sizer.Add(self.show_d0_cb, 0, wx.ALL, 5)

        vbox.Add(plot_sizer, 0, wx.EXPAND | wx.ALL, 10)

        out_box = wx.StaticBox(self, label="Output")
        out_sizer = wx.StaticBoxSizer(out_box, wx.VERTICAL)

        o1 = wx.BoxSizer(wx.HORIZONTAL)
        o1.Add(wx.StaticText(self, label="Results folder:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        self.results_dir_ctrl = wx.TextCtrl(self, value=os.path.join(os.path.expanduser("~"), "rd_delta_results"))
        o1.Add(self.results_dir_ctrl, 1, wx.RIGHT, 5)
        b5 = wx.Button(self, label="Browse")
        b5.Bind(wx.EVT_BUTTON, lambda evt: self._browse_dir_into(self.results_dir_ctrl))
        o1.Add(b5, 0)
        out_sizer.Add(o1, 0, wx.EXPAND | wx.ALL, 5)

        self.png_ctrl = wx.TextCtrl(self, value="delta_plot.png")
        out_sizer.Add(wx.StaticText(self, label="PNG filename (inside results folder):"), 0, wx.LEFT | wx.TOP, 5)
        out_sizer.Add(self.png_ctrl, 0, wx.EXPAND | wx.ALL, 5)

        vbox.Add(out_sizer, 0, wx.EXPAND | wx.ALL, 10)

        btn_row = wx.BoxSizer(wx.HORIZONTAL)
        self.run_btn = wx.Button(self, label="Run Delta + KS Test", size=(200, 35))
        self.run_btn.Bind(wx.EVT_BUTTON, self.on_run)
        btn_row.Add(self.run_btn, 0, wx.RIGHT, 10)

        self.clear_btn = wx.Button(self, label="Clear all", size=(120, 35))
        self.clear_btn.Bind(wx.EVT_BUTTON, self.on_clear_all)
        btn_row.Add(self.clear_btn, 0, wx.RIGHT, 10)

        self.status = wx.StaticText(self, label="Browse and choose Control (D0) and Treated (D42) folders")
        btn_row.Add(self.status, 1, wx.ALIGN_CENTER_VERTICAL)
        vbox.Add(btn_row, 0, wx.EXPAND | wx.ALL, 10)

        vbox.Add(wx.StaticLine(self), 0, wx.EXPAND | wx.ALL, 5)
        vbox.Add(wx.StaticText(self, label="Plot Preview:"), 0, wx.LEFT, 10)

        self.plot_panel = wx.Panel(self)
        self.plot_panel.SetMinSize((600, 400))
        self.plot_canvas = None
        vbox.Add(self.plot_panel, 1, wx.EXPAND | wx.ALL, 10)

        self.SetSizer(vbox)

    def _browse_dir_into(self, ctrl):
        with wx.DirDialog(self, "Select folder") as d:
            if d.ShowModal() == wx.ID_OK:
                ctrl.SetValue(d.GetPath())

    def on_clear_all(self, _):
        self.control_folder_ctrl1.SetValue("")
        self.treatment_folder_ctrl1.SetValue("")
        self.legend_label1_ctrl.SetValue("Dataset 1")

        self.enable_dataset2_cb.SetValue(False)
        self.control_folder_ctrl2.SetValue("")
        self.treatment_folder_ctrl2.SetValue("")
        self.legend_label2_ctrl.SetValue("Dataset 2")

        self.xlim_left_ctrl.SetValue(-100)
        self.xlim_right_ctrl.SetValue(100)
        self.bin_size_ctrl.SetValue(1)
        self.color1_choice.SetSelection(0)
        self.color2_choice.SetSelection(0)
        self.ymin_ctrl.SetValue("")
        self.ymax_ctrl.SetValue("")
        self.ylabel_ctrl.SetValue("Frequency(D42-D0)")
        self.show_d0_cb.SetValue(True)

        self.png_ctrl.SetValue("delta_plot.png")

        if self.plot_canvas:
            self.plot_canvas.Destroy()
            self.plot_canvas = None
            self.plot_panel.Layout()

        self.status.SetLabel("Cleared. Browse and choose Control (D0) and Treated (D42) folders")

    def on_run(self, _):
        try:
            out_dir = self.results_dir_ctrl.GetValue().strip()
            os.makedirs(out_dir, exist_ok=True)
            png_path = os.path.join(out_dir, self.png_ctrl.GetValue().strip() or "delta_plot.png")

            c1 = self.control_folder_ctrl1.GetValue().strip()
            t1 = self.treatment_folder_ctrl1.GetValue().strip()
            if not (c1 and t1 and os.path.isdir(c1) and os.path.isdir(t1)):
                wx.MessageBox("Dataset1 Control/Treated folders are not valid.", "Error", wx.OK | wx.ICON_ERROR)
                return

            control_files1 = resolve_uploaded_paths([c1], include_thresholded=False)
            treatment_files1 = resolve_uploaded_paths([t1], include_thresholded=False)

            df_control1 = read_histogram_files_with_folder_prefix(control_files1, c1, "control", use_thresholded=False)
            df_treatment1 = read_histogram_files_with_folder_prefix(treatment_files1, t1, "treatment", use_thresholded=False)

            merged1 = pd.merge(df_control1, df_treatment1, on="repeat length", how="inner")
            merged1 = normalize_by_column_sum(merged1, start_col=1)

            control_folder_name1 = os.path.basename(os.path.normpath(c1))
            treatment_folder_name1 = os.path.basename(os.path.normpath(t1))

            control_cols1 = [c for c in merged1.columns if c.startswith(f"{control_folder_name1}_")]
            treatment_cols1 = [c for c in merged1.columns if c.startswith(f"{treatment_folder_name1}_")]

            stats_df1 = compute_control_treatment_stats(merged1, control_cols1, treatment_cols1)
            aligned_df1, global_mode1 = align_to_global_mode(stats_df1)

            norm_csv1 = os.path.join(out_dir, "Dataset1_normalized_with_offset.csv")
            aligned_df1.to_csv(norm_csv1, index=False)

            bin_size = self.bin_size_ctrl.GetValue()
            binned_data1 = bin_delta_by_offset_with_auc(aligned_df1, bin_size=bin_size)

            binned_csv1 = os.path.join(out_dir, "Dataset1_binned_by_offset.csv")
            binned_data1.to_csv(binned_csv1, index=False)

            basic_stats1 = calculate_delta_stats(binned_data1)
            basic_stats1["global_mode"] = global_mode1
            auc_stats1 = analyze_auc_stats(binned_data1)
            all_stats1 = {**basic_stats1, **auc_stats1}
            pd.DataFrame([all_stats1]).to_csv(os.path.join(out_dir, "Dataset1_stats.csv"), index=False)

            binned_data2 = None
            all_stats2 = {}

            if self.enable_dataset2_cb.GetValue():
                c2 = self.control_folder_ctrl2.GetValue().strip()
                t2 = self.treatment_folder_ctrl2.GetValue().strip()
                if c2 and t2 and os.path.isdir(c2) and os.path.isdir(t2):
                    control_files2 = resolve_uploaded_paths([c2], include_thresholded=False)
                    treatment_files2 = resolve_uploaded_paths([t2], include_thresholded=False)

                    df_control2 = read_histogram_files_with_folder_prefix(control_files2, c2, "control", use_thresholded=False)
                    df_treatment2 = read_histogram_files_with_folder_prefix(treatment_files2, t2, "treatment", use_thresholded=False)

                    merged2 = pd.merge(df_control2, df_treatment2, on="repeat length", how="inner")
                    merged2 = normalize_by_column_sum(merged2, start_col=1)

                    control_folder_name2 = os.path.basename(os.path.normpath(c2))
                    treatment_folder_name2 = os.path.basename(os.path.normpath(t2))

                    control_cols2 = [c for c in merged2.columns if c.startswith(f"{control_folder_name2}_")]
                    treatment_cols2 = [c for c in merged2.columns if c.startswith(f"{treatment_folder_name2}_")]

                    stats_df2 = compute_control_treatment_stats(merged2, control_cols2, treatment_cols2)
                    aligned_df2, global_mode2 = align_to_global_mode(stats_df2)

                    aligned_df2.to_csv(os.path.join(out_dir, "Dataset2_normalized_with_offset.csv"), index=False)

                    binned_data2 = bin_delta_by_offset_with_auc(aligned_df2, bin_size=bin_size)
                    binned_csv2 = os.path.join(out_dir, "Dataset2_binned_by_offset.csv")
                    binned_data2.to_csv(binned_csv2, index=False)

                    basic_stats2 = calculate_delta_stats(binned_data2)
                    basic_stats2["global_mode"] = global_mode2
                    auc_stats2 = analyze_auc_stats(binned_data2)
                    all_stats2 = {**basic_stats2, **auc_stats2}
                    pd.DataFrame([all_stats2]).to_csv(os.path.join(out_dir, "Dataset2_stats.csv"), index=False)

            # KS(offset-aligned) on binned CSVs if dataset2 enabled
            ks_offset_result = None
            aligned_offsets_csv = None
            ks_offsets_csv = None
            if self.enable_dataset2_cb.GetValue() and binned_data2 is not None:
                name1_for_cols = self.legend_label1_ctrl.GetValue().strip() or "Dataset1"
                name2_for_cols = self.legend_label2_ctrl.GetValue().strip() or "Dataset2"
                aligned_offsets_csv, ks_offsets_csv, ks_offset_result = align_tailtrim_and_ks_on_offsets(
                    ds1_binned_path=binned_csv1,
                    ds2_binned_path=binned_csv2,
                    out_dir=out_dir,
                    name1=name1_for_cols,
                    name2=name2_for_cols,
                    tol=0.0,
                )

            xlim_left = self.xlim_left_ctrl.GetValue()
            xlim_right = self.xlim_right_ctrl.GetValue()
            color1 = self.color1_choice.GetString(self.color1_choice.GetSelection())
            color2 = self.color2_choice.GetString(self.color2_choice.GetSelection())
            ylabel = self.ylabel_ctrl.GetValue().strip() or "Frequency(D42-D0)"

            legend_label1 = self.legend_label1_ctrl.GetValue().strip() or "Dataset 1"
            legend_label2 = self.legend_label2_ctrl.GetValue().strip() or "Dataset 2"

            ylim = None
            try:
                ymin = float(self.ymin_ctrl.GetValue()) if self.ymin_ctrl.GetValue().strip() else None
                ymax = float(self.ymax_ctrl.GetValue()) if self.ymax_ctrl.GetValue().strip() else None
                if ymin is not None and ymax is not None:
                    ylim = (ymin, ymax)
            except Exception:
                pass

            generate_delta_plot_png(
                binned_data1,
                binned_data2,
                output_path=png_path,
                title="Delta plot",
                color1=color1,
                color2=color2,
                label1=legend_label1,
                label2=legend_label2,
                xlim_left=xlim_left,
                xlim_right=xlim_right,
                ylim=ylim,
                bin_size=bin_size,
                show_d0_line=self.show_d0_cb.GetValue(),
                ylabel=ylabel,
            )

            self._show_preview(binned_data1, binned_data2, xlim_left, xlim_right, ylim, color1, color2, legend_label1, legend_label2, ylabel)

            status_msg = f"Done. Saved: {png_path} | Normalized CSVs | Binned CSVs w/ AUC | Stats CSVs"
            if ks_offset_result:
                status_msg += (
                    f" | KS(offset-aligned) p={ks_offset_result['raw_pvalue']:.4g}"
                    f" (adj={ks_offset_result['adjusted_pvalue_fdr_bh']:.4g})"
                    f" n_bins={ks_offset_result['n_bins_after_tail_trim']}"
                )
            if "bias_ratio" in all_stats1:
                status_msg += f" | Dataset1 Bias Ratio={all_stats1['bias_ratio']:.4f}"
            if "bias_ratio" in all_stats2:
                status_msg += f" | Dataset2 Bias Ratio={all_stats2['bias_ratio']:.4f}"

            self.status.SetLabel(status_msg)

        except Exception as e:
            wx.MessageBox(str(e), "Error", wx.OK | wx.ICON_ERROR)
            self.status.SetLabel("Failed.")

    def _show_preview(
        self,
        binned_data1,
        binned_data2=None,
        xlim_left=-100,
        xlim_right=100,
        ylim=None,
        color1="red",
        color2="blue",
        legend_label1="Dataset 1",
        legend_label2="Dataset 2",
        ylabel="Frequency(D42-D0)",
    ):
        if self.plot_canvas:
            self.plot_canvas.Destroy()
            self.plot_canvas = None

        fig = Figure(figsize=(10, 5), dpi=100)
        ax = fig.add_subplot(111)

        ax.plot(binned_data1["Bin_Offset"], binned_data1["Delta_sum"], color=color1, linewidth=2, label=legend_label1)
        if binned_data2 is not None:
            ax.plot(binned_data2["Bin_Offset"], binned_data2["Delta_sum"], color=color2, linewidth=2, label=legend_label2)

        if self.show_d0_cb.GetValue():
            ax.axvline(x=0, color="black", linestyle="--", linewidth=1, alpha=0.7)
        ax.axhline(y=0.0, color="black", linewidth=1, alpha=0.7)

        ax.set_xlim(xlim_left, xlim_right)
        if ylim:
            ax.set_ylim(ylim[0], ylim[1])

        ax.set_title("Delta plot", fontweight="bold")
        ax.set_xlabel("change in repeat size")
        ax.set_ylabel(ylabel)
        ax.grid(False)
        fig.patch.set_facecolor("white")
        ax.set_facecolor("white")

        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(handles, labels, loc="best")

        self.plot_canvas = FigureCanvas(self.plot_panel, -1, fig)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.plot_canvas, 1, wx.EXPAND)
        self.plot_panel.SetSizer(sizer)
        self.plot_panel.Layout()


# -------------------------
# Main frame
# -------------------------
class MainFrame(wx.Frame):
    def __init__(self):
        window_width, window_height, _ = get_screen_dimensions()
        super().__init__(None, title="RD and Delta Plot GUI", size=(window_width, window_height))
        self.SetMinSize((1000, 700))

        nb = wx.Notebook(self)
        self.gen_panel = GenerateHistogramsPanel(nb)
        self.analyze_panel = RepeatDetectorPanel(nb)

        nb.AddPage(self.gen_panel, "Repeat Detector")
        nb.AddPage(self.analyze_panel, "Delta Plot")

        menubar = wx.MenuBar()
        file_menu = wx.Menu()
        file_menu.Append(wx.ID_EXIT, "Exit\tCtrl+Q")
        self.Bind(wx.EVT_MENU, lambda evt: self.Close(), id=wx.ID_EXIT)
        menubar.Append(file_menu, "File")

        help_menu = wx.Menu()
        menubar.Append(help_menu, "Help")

        self.SetMenuBar(menubar)
        center_window_on_screen(self)
        self.Show()


class App(wx.App):
    def OnInit(self):
        print("Starting app...")
        try:
            print(f"Platform: {platform.system()} {platform.release()}")
            print(f"Python: {platform.python_version()}")
        except Exception:
            pass
        MainFrame()
        return True


if __name__ == "__main__":
    App(False).MainLoop()