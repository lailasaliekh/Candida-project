#!/usr/bin/env python3

import os
import re
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import networkx as nx
from scipy.optimize import curve_fit
from RodShapedBacteria import RodShapedBacterium
from findCriticalLengthProb import getCells
import utilities as ut
import fastCellPlotting as fcp
import math
from collections import defaultdict
ut.setMPL()

DEFAULT_OUTPUT_DIR = "/Users/s2507701/OneDrive - University of Edinburgh/leonardo/RUNS/data_production/PaperFigures_April/"

def plot_counts_over_time(data_dir, save_path=None, output_dir=DEFAULT_OUTPUT_DIR):
    file_pattern = os.path.join(data_dir, "biofilm_*.dat")
    files = sorted(glob.glob(file_pattern), key=lambda x: int(re.search(r'(\d+)', x).group(1)))

    time_steps = []
    ca_hyphal_counts = []
    ca_mutant_counts = []
    candida_counts = []
    pseudomonas_counts = []

    is_cam = 'CAm' in data_dir  # check if path indicates mutant mode
    is_pa= 'PA' in data_dir  # check if path indicates PA mode
    is_ca= 'CA' in data_dir  # check if path indicates CA mode
    for file_path in files:
        match = re.search(r'(\d+)', os.path.basename(file_path))
        if not match:
            continue
        time_step = int(match.group(1))
        df = pd.read_csv(file_path, sep="\t")

        pseudomonas_count = df[df["radius"] == 0.5].shape[0]
        pseudomonas_counts.append(pseudomonas_count)

        if is_cam:
            cells = getCells(file_path)
            ca_hyphal = 0
            ca_mutant = 0
            for cell in cells:
                if cell.radius >= 1.0:
                    if (cell.upper_link is None or np.isnan(cell.upper_link) or cell.upper_link == "None") and (cell.lower_link is None or np.isnan(cell.lower_link) or cell.lower_link == "None"):
                        ca_mutant += 1
                    else:
                        ca_hyphal += 1
            ca_hyphal_counts.append(ca_hyphal)
            ca_mutant_counts.append(ca_mutant)
        else:
            candida_count = df[df["radius"] >= 1.00855].shape[0]
            candida_counts.append(candida_count)

        time_steps.append(time_step * 0.1)
    plt.figure(figsize=(5, 3.5))

    if is_cam and is_ca:
        if len(ca_hyphal_counts) == len(time_steps):
            plt.plot(time_steps, ca_hyphal_counts, 'o', color='red', label="CA Hyphal")
        if len(ca_mutant_counts) == len(time_steps):
            plt.plot(time_steps, ca_mutant_counts, 'o', color='#9e003a', label="CA Mutant")

    elif is_cam and is_pa:
        if len(ca_mutant_counts) == len(time_steps):
            plt.plot(time_steps, ca_mutant_counts, 'o', color='#9e003a', label="CA Mutant")
        else:
            print("⚠️ Length mismatch: ca_mutant_counts")
        if len(pseudomonas_counts) == len(time_steps):
            plt.plot(time_steps, pseudomonas_counts, 's', color='cyan', label="PA-SA")
        else:
            print("⚠️ Length mismatch: pseudomonas_counts")

    elif is_ca and is_pa:
        if len(candida_counts) == len(time_steps):
            plt.plot(time_steps, candida_counts, 'o', color='red', label="CA Hyphal")
        if len(pseudomonas_counts) == len(time_steps):
            plt.plot(time_steps, pseudomonas_counts, 's', color='cyan', label="PA-SA")
    elif is_ca and not is_cam and not is_pa:
        if len(candida_counts) == len(time_steps):
            plt.plot(time_steps, candida_counts, 'o', color='red', label="CA Hyphal")
        else:
            print("⚠️ Length mismatch: candida_counts")
    elif is_pa and not is_cam and not is_ca:
        if len(pseudomonas_counts) == len(time_steps):
            plt.plot(time_steps, pseudomonas_counts, 's', color='cyan', label="PA-SA")
        else:
            print("⚠️ Length mismatch: pseudomonas_counts")
    elif is_cam and not is_ca and not is_pa:
        if len(ca_mutant_counts) == len(time_steps):
            plt.plot(time_steps, ca_mutant_counts, 'o', color='#9e003a', label="CA Mutant")
        else:
            print("⚠️ Length mismatch: ca_mutant_counts")
    else:
        print("⚠️ Plot skipped: unknown or unsupported combination of conditions.")

    plt.xlabel("Time (h)")
    plt.ylabel("Cell/Segment Count")
    plt.yscale("log")
    plt.legend()
    plt.tight_layout()

    # plt.figure(figsize=(5, 3.5))
    # if is_cam and is_ca:
    #     plt.plot(time_steps, ca_hyphal_counts, 'o', color='red', label="CA Hyphal")
    #     plt.plot(time_steps, ca_mutant_counts, 'o', color='#9e003a', label="CA Mutant")
    # if is_cam and is_pa:
    #     plt.plot(time_steps, ca_mutant_counts, 'o', color='#9e003a', label="CA Mutant")
    #     plt.plot(time_steps, pseudomonas_counts, 's', color='cyan', label="PA-SA")
    # if is_ca and is_pa:
    #     plt.plot(time_steps, candida_counts, 'o', color='red', label="CA Hyphal")
    #     plt.plot(time_steps, pseudomonas_counts, 's', color='cyan', label="PA-SA")
    # plt.xlabel("Time (h)")
    # plt.ylabel("Cell/Segment Count")
    # plt.yscale("log")
    # plt.legend()
    # plt.tight_layout()

    if save_path is None:
        parts = os.path.normpath(data_dir).split(os.sep)
        tag = "_".join(parts[-2:])
        os.makedirs(output_dir, exist_ok=True)
        save_path = os.path.join(output_dir, f"counts_{tag}.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved cell count plot to: {save_path}")
    plt.show()

def plot_average_counts_over_repeats(parent_dir, output_dir=DEFAULT_OUTPUT_DIR):
    repeat_dirs = sorted(glob.glob(os.path.join(parent_dir, "repeat*")))
    all_ca_hyphal = {}
    all_ca_mutant = {}
    all_candida = {}
    all_pseudomonas = {}
    all_times = set()

    is_cam = 'CAm' in parent_dir  # check if path indicates mutant mode
    is_pa= 'PA' in parent_dir  # check if path indicates PA mode
    is_ca= 'CA' in parent_dir  # check if path indicates CA mode


    for repeat_dir in repeat_dirs:
        file_pattern = os.path.join(repeat_dir, "biofilm_*.dat")
        files = sorted(glob.glob(file_pattern), key=lambda x: int(re.search(r'(\d+)', x).group(1)))
        for file_path in files:
            match = re.search(r'(\d+)', os.path.basename(file_path))
            if not match:
                continue
            time_step = int(match.group(1))
            time_h = time_step * 0.1
            df = pd.read_csv(file_path, sep="\t")
            pseudomonas = df[df["radius"] == 0.5].shape[0]
            all_pseudomonas.setdefault(time_h, []).append(pseudomonas)

            if is_cam:
                cells = getCells(file_path)
                hyphal = 0
                mutant = 0
                for cell in cells:
                    if cell.radius >= 1.0:
                        if (cell.upper_link is None or np.isnan(cell.upper_link) or cell.upper_link == "None") and (cell.lower_link is None or np.isnan(cell.lower_link) or cell.lower_link == "None"):
                            mutant += 1
                        else:
                            hyphal += 1
                all_ca_hyphal.setdefault(time_h, []).append(hyphal)
                all_ca_mutant.setdefault(time_h, []).append(mutant)
            else:
                candida = df[df["radius"] >= 1.00855].shape[0]
                all_candida.setdefault(time_h, []).append(candida)

            all_times.add(time_h)

    sorted_times = sorted(all_times)
    avg_ca_hyphal, sem_ca_hyphal = [], []
    avg_ca_mutant, sem_ca_mutant = [], []
    avg_candida, sem_candida = [], []
    avg_pseudomonas, sem_pseudomonas = [], []

    for t in sorted_times:
        p_vals = all_pseudomonas.get(t, [])
        while len(p_vals) < len(repeat_dirs):
            p_vals.append(0)
        avg_pseudomonas.append(np.mean(p_vals))
        sem_pseudomonas.append(np.std(p_vals, ddof=1) / np.sqrt(len(p_vals)))

        if is_cam:
            h_vals = all_ca_hyphal.get(t, [])
            m_vals = all_ca_mutant.get(t, [])
            while len(h_vals) < len(repeat_dirs):
                h_vals.append(0)
            while len(m_vals) < len(repeat_dirs):
                m_vals.append(0)
            avg_ca_hyphal.append(np.mean(h_vals))
            sem_ca_hyphal.append(np.std(h_vals, ddof=1) / np.sqrt(len(h_vals)))
            avg_ca_mutant.append(np.mean(m_vals))
            sem_ca_mutant.append(np.std(m_vals, ddof=1) / np.sqrt(len(m_vals)))
        else:
            c_vals = all_candida.get(t, [])
            while len(c_vals) < len(repeat_dirs):
                c_vals.append(0)
            avg_candida.append(np.mean(c_vals))
            sem_candida.append(np.std(c_vals, ddof=1) / np.sqrt(len(c_vals)))

    plt.figure(figsize=(5.5, 4))
    # Only plot what is present and valid
    if is_cam and is_ca:
        if len(avg_ca_hyphal) == len(sorted_times):
            plt.errorbar(sorted_times, avg_ca_hyphal, yerr=sem_ca_hyphal, fmt='o', color='red', label="CA Hyphal")
        if len(avg_ca_mutant) == len(sorted_times):
            plt.errorbar(sorted_times, avg_ca_mutant, yerr=sem_ca_mutant, fmt='o', color='#9e003a', label="CA Mutant")

    elif is_cam and is_pa:
        if len(avg_ca_mutant) == len(sorted_times):
            plt.errorbar(sorted_times, avg_ca_mutant, yerr=sem_ca_mutant, fmt='o', color='#9e003a', label="CA Mutant")
        if len(avg_pseudomonas) == len(sorted_times):
            plt.errorbar(sorted_times, avg_pseudomonas, yerr=sem_pseudomonas, fmt='s', color='cyan', label="PA-SA")

    elif is_ca and is_pa:
        if len(avg_candida) == len(sorted_times):
            plt.errorbar(sorted_times, avg_candida, yerr=sem_candida, fmt='o', color='red', label="CA Hyphal")
        if len(avg_pseudomonas) == len(sorted_times):
            plt.errorbar(sorted_times, avg_pseudomonas, yerr=sem_pseudomonas, fmt='s', color='cyan', label="PA-SA")

    else:
        print("⚠️ Plot skipped: unknown or unsupported combination of conditions.")

    # if is_cam and is_ca:
    #     plt.errorbar(sorted_times, avg_ca_hyphal, yerr=sem_ca_hyphal, fmt='o', color='red', label="CA Hyphal")
    #     plt.errorbar(sorted_times, avg_ca_mutant, yerr=sem_ca_mutant, fmt='o', color='#9e003a', label="CA Mutant")
    # if is_cam and is_pa:
    #     plt.errorbar(sorted_times, avg_ca_mutant, yerr=sem_ca_mutant, fmt='o', color='#9e003a', label="CA Mutant")
    #     plt.errorbar(sorted_times, avg_pseudomonas, yerr=sem_pseudomonas, fmt='s', color='cyan', label="PA-SA")
    # if is_pa and is_ca:
    #     plt.errorbar(sorted_times, avg_candida, yerr=sem_candida, fmt='o', color='red', label="CA Hyphal")
    #     plt.errorbar(sorted_times, avg_pseudomonas, yerr=sem_pseudomonas, fmt='s', color='cyan', label="PA-SA")
    plt.xlabel("Time (h)")
    plt.ylabel("Average Cell/Segment Count")
    plt.legend()
    plt.tight_layout()

    parts = os.path.normpath(parent_dir).split(os.sep)
    tag = "_".join(parts[-2:])
    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"avg_counts_{tag}.pdf")
    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved average cell count plot to: {save_path}")
    plt.show()


def plot_cells_grid(parent_dir, output_path=None, num_snapshots=7, output_dir=DEFAULT_OUTPUT_DIR):
    repeat_dirs = sorted(glob.glob(os.path.join(parent_dir, "repeat*")))
    num_repeats = len(repeat_dirs)

    fig, axes = plt.subplots(num_repeats, num_snapshots, figsize=(15, 3 * num_repeats),
                             constrained_layout=True, facecolor='w')

    def plotCells(ax, file):
        RodShapedBacterium.sus_vis_radius_factor = 0.7
        dat = pd.read_csv(file, sep='\t')
        cells = getCells(file)
        maxx, minx = dat['pos_x'].max() + 3, dat['pos_x'].min() - 3
        maxy, miny = dat['pos_y'].max() + 3, dat['pos_y'].min() - 3
        X, Y = maxx - minx, maxy - miny
        X_c, Y_c = 0.5 * (maxx + minx), 0.5 * (maxy + miny)
        if X >= Y:
            Y = X
            miny, maxy = Y_c - 0.5 * Y, Y_c + 0.5 * Y
        else:
            X = Y
            minx, maxx = X_c - 0.5 * X, X_c + 0.5 * X

        fcp.addAllCellsToPlot(cells, ax, ax_rng=maxx - minx, show_id=False, ec='w')
        fcp.addChainLinks(cells, ax, ax_rng=maxx - minx, show_id=False, colour='r')

        scale = 1.17
        wall_color = 'k'
        ax.plot([-75.1/scale, 75.1/scale], [-0.1/scale, -0.1/scale], color='k', linestyle='dashed', alpha=0.5)
        ax.plot([-75.1/scale, 75.1/scale], [150/scale, 150/scale], color=wall_color, alpha=0.6)
        ax.plot([-75.1/scale, -75.1/scale], [-0.1/scale, 150/scale], color=wall_color, alpha=0.6)
        ax.plot([75.1/scale, 75.1/scale], [-0.1/scale, 150/scale], color=wall_color, alpha=0.6)

        ax.set_xlim([-80, 80])
        ax.set_ylim([-5, 160])
        ax.axis('scaled')
        ax.axis('off')

    for r, repeat_dir in enumerate(repeat_dirs):
        file_pattern = os.path.join(repeat_dir, "biofilm_*.dat")
        files = sorted(glob.glob(file_pattern))
        if len(files) < num_snapshots:
            selected_files = files
        else:
            selected_indices = np.linspace(5, len(files) - 1, num_snapshots, dtype=int)
            selected_files = [files[i] for i in selected_indices]


            # # Let's take 4 snapshots spaced by len(files) // 4
            # selected_indices = [5]  # start from the beginning (or you can use 1 if you want to skip t=0)

            # # Add 3 more points evenly spaced
            # quarter = len(files) // 4
            # selected_indices += [quarter, 2 * quarter, 3 * quarter]

            # Get corresponding files
            selected_files = [files[i] for i in selected_indices]
        for c, (ax, file) in enumerate(zip(axes[r], selected_files)):
            plotCells(ax, file)
            match = re.search(r'biofilm_(\d+)\.dat$', os.path.basename(file))
            if match:
                frame_number = int(match.group(1))
                time_hours = frame_number * 0.1
                ax.set_title(f"{time_hours:.1f} h", color='k', fontsize=12)
        axes[r, 0].set_ylabel(f"Repeat {r+1}", fontsize=12, color='k')

    if output_path is None:
        parts = os.path.normpath(parent_dir).split(os.sep)
        tag = "_".join(parts[-2:])
        os.makedirs(output_dir, exist_ok=True)
        output_path_pdf = os.path.join(output_dir, f"snapshots_{tag}.pdf")
        output_path_png = os.path.join(output_dir, f"snapshots_{tag}.png")

    fig.savefig(output_path_pdf, format="pdf", dpi=600, bbox_inches="tight")
    fig.savefig(output_path_png, dpi=600, bbox_inches="tight")
    print(f"Saved snapshot grid to: {output_path_pdf}")
    plt.show()

def plot_cells_grid1(parent_dir, output_path=None, num_snapshots=4, output_dir=DEFAULT_OUTPUT_DIR):
    repeat_dirs = sorted(glob.glob(os.path.join(parent_dir, "repeat*")))
    num_repeats = len(repeat_dirs)

    fig, axes = plt.subplots(num_repeats, num_snapshots, figsize=(15, 3 * num_repeats),
                             constrained_layout=True, facecolor='w')

    # Ensure axes is always a 2D array
    if num_repeats == 1 and num_snapshots == 1:
        axes = np.array([[axes]])
    elif num_repeats == 1:
        axes = np.expand_dims(axes, axis=0)  # shape (1, num_snapshots)
    elif num_snapshots == 1:
        axes = np.expand_dims(axes, axis=1)  # shape (num_repeats, 1)

    def plotCells(ax, file):
        RodShapedBacterium.sus_vis_radius_factor = 0.7
        dat = pd.read_csv(file, sep='\t')
        cells = getCells(file)
        maxx, minx = dat['pos_x'].max() + 3, dat['pos_x'].min() - 3
        maxy, miny = dat['pos_y'].max() + 3, dat['pos_y'].min() - 3
        X, Y = maxx - minx, maxy - miny
        X_c, Y_c = 0.5 * (maxx + minx), 0.5 * (maxy + miny)
        if X >= Y:
            Y = X
            miny, maxy = Y_c - 0.5 * Y, Y_c + 0.5 * Y
        else:
            X = Y
            minx, maxx = X_c - 0.5 * X, X_c + 0.5 * X

        fcp.addAllCellsToPlot(cells, ax, ax_rng=maxx - minx, show_id=False, ec='w')
        fcp.addChainLinks(cells, ax, ax_rng=maxx - minx, show_id=False, colour='r')

        scale = 1.17
        wall_color = 'k'
        ax.plot([-75.1/scale, 75.1/scale], [-0.1/scale, -0.1/scale], color='k', linestyle='dashed', alpha=0.5)
        ax.plot([-75.1/scale, 75.1/scale], [150/scale, 150/scale], color=wall_color, alpha=0.6)
        ax.plot([-75.1/scale, -75.1/scale], [-0.1/scale, 150/scale], color=wall_color, alpha=0.6)
        ax.plot([75.1/scale, 75.1/scale], [-0.1/scale, 150/scale], color=wall_color, alpha=0.6)

        ax.set_xlim([-80, 80])
        ax.set_ylim([-5, 160])
        ax.axis('scaled')
        ax.axis('off')

    for r, repeat_dir in enumerate(repeat_dirs):
        file_pattern = os.path.join(repeat_dir, "biofilm_*.dat")
        files = sorted(glob.glob(file_pattern))
        if len(files) < num_snapshots:
            selected_files = files
        else:
            selected_indices = [5]  # start from a consistent early point
            quarter = len(files) // 4
            selected_indices += [quarter, 2 * quarter, -1]
            selected_files = [files[i] for i in selected_indices]
        for c, (ax, file) in enumerate(zip(axes[r], selected_files)):
            plotCells(ax, file)
            match = re.search(r'biofilm_(\d+)\.dat$', os.path.basename(file))
            if match:
                frame_number = int(match.group(1))
                time_hours = frame_number * 0.1
                ax.set_title(f"{time_hours:.1f} h", color='k', fontsize=12)
        axes[r, 0].set_ylabel(f"Repeat {r+1}", fontsize=12, color='k')
        # Turn off unused subplots in this row
        for c in range(len(selected_files), num_snapshots):
            axes[r, c].axis('off')
    if output_path is None:
        parts = os.path.normpath(parent_dir).split(os.sep)
        tag = "_".join(parts[-2:])
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"snapshots_{tag}.png")

    fig.savefig(output_path, format="png", dpi=600, bbox_inches="tight")
    print(f"Saved snapshot grid to: {output_path}")
    plt.show()

def plot_chain_growth(parent_dir, output_dir=DEFAULT_OUTPUT_DIR):
    import networkx as nx
    from scipy.optimize import curve_fit

    def growth_model(t, B0, lambda_):
        return B0 * np.exp(lambda_ * t)

    is_cam = 'CAm' in parent_dir
    is_ca = 'CA' in parent_dir
    is_pa = 'PA' in parent_dir

    repeat_dirs = sorted(glob.glob(os.path.join(parent_dir, "repeat*")))
    if not repeat_dirs:
        print("No repeat directories found.")
        return

    time_to_chains = {}
    max_time = 0
    if is_ca:
        for repeat_dir in repeat_dirs:
            files = sorted(glob.glob(os.path.join(repeat_dir, "biofilm_*.dat")), key=lambda x: int(re.search(r'biofilm_(\d+)', os.path.basename(x)).group(1)))
            for file in files:
                match = re.search(r'biofilm_(\d+)', os.path.basename(file))
                if not match:
                    continue
                time = int(match.group(1)) * 0.1
                max_time = max(max_time, time)

                data = pd.read_csv(file, sep='\t')
                if 'upper_link' not in data.columns or 'lower_link' not in data.columns:
                    continue
                def is_valid_link(val):
                    return not (val is None or pd.isna(val) or val == "None" or val == -1)

                G = nx.Graph()
                for idx, row in data.iterrows():
                    upper = row['upper_link']
                    lower = row['lower_link']
                    if is_valid_link(upper):
                        G.add_edge(idx, int(upper))
                    if is_valid_link(lower):
                        G.add_edge(idx, int(lower))
                # valid_data = data[(data['upper_link'] != -1) | (data['lower_link'] != -1)]

                # G = nx.Graph()
                # for idx, row in valid_data.iterrows():
                #     if row['upper_link'] != -1:
                #         G.add_edge(idx, row['upper_link'])
                #     if row['lower_link'] != -1:
                #         G.add_edge(idx, row['lower_link'])

                connected_components = list(nx.connected_components(G))
                num_chains = sum(1 for comp in connected_components if len(comp) >= 2)

                time_to_chains.setdefault(time, []).append(num_chains)

        sorted_times = sorted(time_to_chains.keys())
        avg_chains = []
        sem_chains = []
        for t in sorted_times:
            vals = time_to_chains[t]
            avg_chains.append(np.mean(vals))
            sem_chains.append(np.std(vals, ddof=1) / np.sqrt(len(vals)))

        time_points = np.array(sorted_times)
        avg_chains = np.array(avg_chains)

        if len(time_points) == 0:
            print("No valid chain data found.")
            return

        params, _ = curve_fit(growth_model, time_points, avg_chains, p0=[avg_chains[0], 0.5])
        B0_fit, lambda_fit = params 
        time_doubling= np.log(2)/lambda_fit
        plt.figure(figsize=(5, 3))
        plt.errorbar(time_points, avg_chains, yerr=sem_chains, fmt='o', label="Mean ± SEM", color='red')
        plt.plot(time_points, growth_model(time_points, *params), label=f"Fit: Doubling Time = {time_doubling:.4f} hr", color='k')
        plt.xlabel("Time (h)")
        plt.ylabel("Average Number of Chains")
        plt.yscale('log')
        plt.legend()
        plt.tight_layout()

        tag = "_".join(os.path.normpath(parent_dir).split(os.sep)[-2:])
        os.makedirs(output_dir, exist_ok=True)
        save_path = os.path.join(output_dir, f"chain_growth_{tag}.pdf")
        plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
        print(f"Saved chain growth plot to: {save_path}")
        print(f"Fitted Growth Rate (lambda): {lambda_fit:.4f} per hour")
        print(f"Initial Number of Chains (B0): {B0_fit:.2f}")
        plt.show()

def plot_chain_growth1(parent_dir, output_dir=DEFAULT_OUTPUT_DIR):


    def growth_model(t, B0, lambda_):
        return B0 * np.exp(lambda_ * t)

    is_cam = 'CAm' in parent_dir
    is_ca = 'CA' in parent_dir
    is_pa = 'PA' in parent_dir

    repeat_dirs = sorted(glob.glob(os.path.join(parent_dir, "repeat*")))
    all_files = []
    for repeat_dir in repeat_dirs:
        files = glob.glob(os.path.join(repeat_dir, "biofilm_*.dat"))
        all_files.extend(files)
    if is_ca:
        all_files = [f for f in all_files if "biofilm" in os.path.basename(f)]
        files = sorted(all_files, key=lambda x: int(re.search(r'biofilm_(\d+)', x).group(1)))

        time_points = []
        chain_counts = []

        for file in files:
            match = re.search(r'biofilm_(\d+)', os.path.basename(file))
            if not match:
                continue
            time = int(match.group(1)) * 0.1

            data = pd.read_csv(file, sep='\t')
            if 'upper_link' not in data.columns or 'lower_link' not in data.columns:
                continue

            valid_data = data[(data['upper_link'] != -1) | (data['lower_link'] != -1)]

            G = nx.Graph()
            for idx, row in valid_data.iterrows():
                if row['upper_link'] != -1:
                    G.add_edge(idx, row['upper_link'])
                if row['lower_link'] != -1:
                    G.add_edge(idx, row['lower_link'])

            connected_components = list(nx.connected_components(G))
            num_chains = sum(1 for comp in connected_components if len(comp) >= 2)

            time_points.append(time)
            chain_counts.append(num_chains)

        if len(time_points) == 0:
            print("No valid chain data found.")
            return

        time_points = np.array(time_points)
        chain_counts = np.array(chain_counts)

        params, _ = curve_fit(growth_model, time_points, chain_counts, p0=[chain_counts[0], 0.5])
        B0_fit, lambda_fit = params

        plt.figure(figsize=(5, 3))
        plt.scatter(time_points, chain_counts, label="Data", color='red')
        plt.plot(time_points, growth_model(time_points, *params), label=f"Fit: growth = {lambda_fit:.4f}", color='k')
        plt.xlabel("Time (hours)")
        plt.ylabel("Number of Chains")
        plt.yscale('log')
        plt.legend()
        plt.tight_layout()
        
        tag = "_".join(os.path.normpath(parent_dir).split(os.sep)[-2:])
        os.makedirs(output_dir, exist_ok=True)
        save_path = os.path.join(output_dir, f"chain_growth_{tag}.pdf")
        plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
        print(f"Saved chain growth plot to: {save_path}")
        print(f"Fitted Growth Rate (lambda): {lambda_fit:.4f} per hour")
        print(f"Initial Number of Chains (B0): {B0_fit:.2f}")
        plt.show()
    # else:
    #     print("Chain growth analysis is not available for this dataset.")   
        

def plot_kymograph_average(parent_dir, output_dir=DEFAULT_OUTPUT_DIR):
    grid_size = 2
    repeat_dirs = sorted(glob.glob(os.path.join(parent_dir, "repeat*")))
    x_min, x_max, y_min, y_max = None, None, None, None
    repeat_time_steps = []

    is_cam = 'CAm' in parent_dir
    is_ca = 'CA' in parent_dir
    is_pa = 'PA' in parent_dir

    for repeat_dir in repeat_dirs:
        file_pattern = os.path.join(repeat_dir, "biofilm_*.dat")
        files = sorted(glob.glob(file_pattern), key=lambda x: int(re.search(r'biofilm_(\d+)\.dat$', x).group(1)))
        repeat_time_steps.append(len(files))

        for file_path in files:
            df = pd.read_csv(file_path, sep="\t")
            x_min = min(x_min, df["pos_x"].min()) if x_min is not None else df["pos_x"].min()
            x_max = max(x_max, df["pos_x"].max()) if x_max is not None else df["pos_x"].max()
            y_min = min(y_min, df["pos_y"].min()) if y_min is not None else df["pos_y"].min()
            y_max = max(y_max, df["pos_y"].max()) if y_max is not None else df["pos_y"].max()

    x_bins = np.arange(x_min, x_max + grid_size, grid_size)
    y_bins = np.arange(y_min, y_max + grid_size, grid_size)
    grid_shape = (len(y_bins) - 1, len(x_bins) - 1)

    max_time_steps = max(repeat_time_steps)
    min_time_steps = min(repeat_time_steps)
    truncate_data = False
    num_time_steps = min_time_steps if truncate_data else max_time_steps
    total_cell_counts = np.zeros((len(y_bins) - 1, num_time_steps, 3))  # (y_bins, time_steps, [CA, CAm, PA])
    num_repeats = len(repeat_dirs)

    for repeat_dir in repeat_dirs:
        file_pattern = os.path.join(repeat_dir, "biofilm_*.dat")
        files = sorted(glob.glob(file_pattern), key=lambda x: int(re.search(r'biofilm_(\d+)\.dat$', x).group(1)))
        padded_files = files[:min_time_steps] if truncate_data else files + [None] * (max_time_steps - len(files))

        for t_index, file_path in enumerate(padded_files):
            if file_path is None:
                continue
            df = pd.read_csv(file_path, sep="\t")
            for _, row in df.iterrows():
                x_index = np.digitize(row["pos_x"], x_bins) - 1
                y_index = np.digitize(row["pos_y"], y_bins) - 1
                if 0 <= x_index < grid_shape[1] and 0 <= y_index < grid_shape[0]:
                    if row["radius"] >= 1.0:
                        if is_cam and is_ca:
                            has_links = (row.upper_link is None or np.isnan(row.upper_link) or row.upper_link == "None") and (row.lower_link is None or np.isnan(row.lower_link) or row.lower_link == "None")
                            if has_links:
                                total_cell_counts[y_index, t_index, 0] += 1  # CA Hyphal (red)
                            else:
                                total_cell_counts[y_index, t_index, 1] += 1  # CA Mutant (cranberry)
                        else:
                            total_cell_counts[y_index, t_index, 0] += 1  # CA or CAm alone
                    elif row["radius"] == 0.5:
                        total_cell_counts[y_index, t_index, 2] += 1  # PA (blue)

    avg_ca = total_cell_counts[:, :, 0] / num_repeats
    avg_cam = total_cell_counts[:, :, 1] / num_repeats
    avg_pa = total_cell_counts[:, :, 2] / num_repeats

    bin_factor = 2
    summed_ca = np.add.reduceat(avg_ca, np.arange(0, avg_ca.shape[0], bin_factor), axis=0)
    summed_cam = np.add.reduceat(avg_cam, np.arange(0, avg_cam.shape[0], bin_factor), axis=0)
    summed_pa = np.add.reduceat(avg_pa, np.arange(0, avg_pa.shape[0], bin_factor), axis=0)

    num_boxes_per_row = grid_shape[1]
    red_channel = 255 * (summed_ca / num_boxes_per_row)
    cranberry_channel = 100 * (summed_cam / num_boxes_per_row)  # less intense red for mutant
    blue_channel = 255 * (summed_pa / num_boxes_per_row)
    red_channel = np.nan_to_num(np.clip(red_channel, 0, 255))
    cranberry_channel = np.nan_to_num(np.clip(cranberry_channel, 0, 255))
    blue_channel = np.nan_to_num(np.clip(blue_channel, 0, 255))

    rgb_image = np.zeros((summed_ca.shape[0], num_time_steps, 3), dtype=np.uint8)
    rgb_image[:, :, 0] = red_channel
    rgb_image[:, :, 1] = cranberry_channel  # show mutant Candida in cranberry tone
    rgb_image[:, :, 2] = blue_channel

    fig, ax = plt.subplots(figsize=(2.5, 2.2))
    im = ax.imshow(rgb_image, aspect="auto", extent=[0, num_time_steps * 0.1, y_min * 1.17, y_max * 1.17], origin="lower")
    ax.set_xlabel("Time (h)")
    ax.set_ylabel("Y Position ($\\mu$m)")

    tag = "_".join(os.path.normpath(parent_dir).split(os.sep)[-2:])
    title = tag.replace("_", " ")
    ax.set_title(title)
    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"kymograph_avg_{tag}.pdf")
    fig.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved averaged kymograph to: {save_path}")
    plt.show()



### Averaging over existing frames and not filling with zeroes
    

def plot_average_counts_over_repeats_without_zerosfilling(parent_dir, output_dir=DEFAULT_OUTPUT_DIR):
    repeat_dirs = sorted(glob.glob(os.path.join(parent_dir, "repeat*")))

    counts_data = {}

    is_cam = 'CAm' in parent_dir
    is_pa = 'PA' in parent_dir
    is_ca = 'CA' in parent_dir

    for repeat_dir in repeat_dirs:
        file_pattern = os.path.join(repeat_dir, "biofilm_*.dat")
        files = sorted(glob.glob(file_pattern), key=lambda x: int(re.search(r'(\d+)', x).group(1)))
        for file_path in files:
            match = re.search(r'(\d+)', os.path.basename(file_path))
            if not match:
                continue
            time_step = int(match.group(1))
            time_h = time_step * 0.1
            df = pd.read_csv(file_path, sep="\t")

            counts = counts_data.setdefault(time_h, {"pseudomonas": [], "ca_hyphal": [], "ca_mutant": [], "candida": []})

            pseudomonas = df[df["radius"] == 0.5].shape[0]
            counts["pseudomonas"].append(pseudomonas)

            if is_cam:
                cells = getCells(file_path)
                hyphal, mutant = 0, 0
                for cell in cells:
                    if cell.radius >= 1.0:
                        if (cell.upper_link is None or np.isnan(cell.upper_link) or cell.upper_link == "None") and \
                           (cell.lower_link is None or np.isnan(cell.lower_link) or cell.lower_link == "None"):
                            mutant += 1
                        else:
                            hyphal += 1
                counts["ca_hyphal"].append(hyphal)
                counts["ca_mutant"].append(mutant)
            else:
                candida = df[df["radius"] >= 1.00855].shape[0]
                counts["candida"].append(candida)

    sorted_times = sorted(counts_data.keys())

    plt.figure(figsize=(5.5, 4))

    for label, color, key in [("CA Hyphal", 'red', "ca_hyphal"), ("CA Mutant", '#9e003a', "ca_mutant"),
                              ("Candida", 'red', "candida"), ("PA-SA", 'cyan', "pseudomonas")]:
        avg, sem, times = [], [], []
        for t in sorted_times:
            vals = counts_data[t][key]
            if vals:  # Only consider if data is present
                avg.append(np.mean(vals))
                sem.append(np.std(vals, ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0)
                times.append(t)

        if len(avg) == 0:
            continue  # Skip if no data

        if ((key == "ca_hyphal" and is_cam and is_ca) or
            (key == "ca_mutant" and is_cam and (is_ca or is_pa)) or
            (key == "candida" and not is_cam and is_ca) or
            (key == "pseudomonas" and is_pa)):
            plt.errorbar(times, avg, yerr=sem, fmt='o', color=color, label=label)

    plt.xlabel("Time (h)")
    plt.ylabel("Average Cell/Segment Count")
    plt.yscale("log")
    plt.legend()
    plt.tight_layout()

    parts = os.path.normpath(parent_dir).split(os.sep)
    tag = "_".join(parts[-2:])
    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"avg_counts_{tag}_without_filling.pdf")
    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved average cell count plot to: {save_path}")
    plt.show()

def compute_kymograph_average_without_zerosfilling(parent_dir, output_dir=DEFAULT_OUTPUT_DIR):
    repeat_dirs = sorted(glob.glob(os.path.join(parent_dir, "repeat*")))
    is_cam = 'CAm' in parent_dir
    is_ca = 'CA' in parent_dir
    is_pa = 'PA' in parent_dir

    all_data = {}
    x_min, x_max, y_min, y_max = None, None, None, None

    for repeat_dir in repeat_dirs:
        files = sorted(glob.glob(os.path.join(repeat_dir, "biofilm_*.dat")), key=lambda x: int(re.search(r'biofilm_(\d+)', os.path.basename(x)).group(1)))
        for file in files:
            time = int(re.search(r'biofilm_(\d+)', os.path.basename(file)).group(1)) * 0.1
            df = pd.read_csv(file, sep='\t')
            x_min = min(x_min, df["pos_x"].min()) if x_min is not None else df["pos_x"].min()
            x_max = max(x_max, df["pos_x"].max()) if x_max is not None else df["pos_x"].max()
            y_min = min(y_min, df["pos_y"].min()) if y_min is not None else df["pos_y"].min()
            y_max = max(y_max, df["pos_y"].max()) if y_max is not None else df["pos_y"].max()
            all_data.setdefault(time, []).append(df)

    x_bins = np.arange(x_min, x_max + 2, 2)
    y_bins = np.arange(y_min, y_max , 2)
    grid_shape = (len(y_bins) - 1, len(x_bins) - 1)

    def build_channel_matrix(data_list, key_func):
        time_steps = sorted(data_list.keys())
        kymo = []
        for t in time_steps:
            stacked = np.zeros((grid_shape[0],))
            counts = 0
            for df in data_list[t]:
                frame = np.zeros(grid_shape)
                for _, row in df.iterrows():
                    if key_func(row):
                        x_idx = np.digitize(row["pos_x"], x_bins) - 1
                        y_idx = np.digitize(row["pos_y"], y_bins) - 1
                        if 0 <= x_idx < grid_shape[1] and 0 <= y_idx < grid_shape[0]:
                            frame[y_idx, x_idx] += 1
                stacked += frame.sum(axis=1)
                counts += 1
            kymo.append(stacked / counts if counts > 0 else stacked)
        return np.array(kymo).T, time_steps

    ca_mask = lambda row: row["radius"] >= 1.0 and not (is_cam and ((row.upper_link is None or row.upper_link == "None" or np.isnan(row.upper_link)) and (row.lower_link is None or row.lower_link == "None" or np.isnan(row.lower_link))))
    cam_mask = lambda row: row["radius"] >= 1.0 and (row.upper_link is None or row.upper_link == "None" or np.isnan(row.upper_link)) and (row.lower_link is None or row.lower_link == "None" or np.isnan(row.lower_link))
    pa_mask = lambda row: row["radius"] == 0.5

    ca_img, time_steps = build_channel_matrix(all_data, ca_mask)
    cam_img, _ = build_channel_matrix(all_data, cam_mask)
    pa_img, _ = build_channel_matrix(all_data, pa_mask)
## revised colour scheme
    # rgb_img = np.zeros((ca_img.shape[0], ca_img.shape[1], 3), dtype=np.uint8)
    #     # Red = CA hyphal
    # rgb_img[:, :, 0] = np.clip(ca_img * 10 + cam_img * 5, 0, 255).astype(np.uint8)  # red + cranberry contribution
    # # Green = muted / none
    # rgb_img[:, :, 1] = np.zeros_like(ca_img, dtype=np.uint8)  # No green
    # # Blue = PA + part of CAm to give cranberry tint
    # rgb_img[:, :, 2] = np.clip(pa_img * 10 + cam_img * 3, 0, 255).astype(np.uint8)

    # Revised color scheme for white CAm
    rgb_img = np.zeros((ca_img.shape[0], ca_img.shape[1], 3), dtype=np.uint8)

    # Red channel: CA hyphal + CAm
    rgb_img[:, :, 0] = np.clip(ca_img * 10 + cam_img * 10, 0, 255).astype(np.uint8)  

    # Green channel: only CAm
    rgb_img[:, :, 1] = np.clip(cam_img * 10, 0, 255).astype(np.uint8)  

    # Blue channel: PA + CAm
    rgb_img[:, :, 2] = np.clip(pa_img * 10 + cam_img * 10, 0, 255).astype(np.uint8)


    plt.figure(figsize=(2.5, 2.2))
    plt.imshow(rgb_img, aspect='auto', extent=[0, max(time_steps), y_min * 1.17, y_max * 1.17], origin='lower')
    plt.xlabel("Time (h)")
    plt.ylabel("Y Position ($\mu$m)")
    plt.xlim(0,11.1)
    plt.ylim(0,150)
    tag = "_".join(os.path.normpath(parent_dir).split(os.sep)[-2:])
    # plt.title(tag.replace("_", " "))
    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"kymograph_avg_{tag}_without_filling_until_11hours_bin_tuesday1.pdf")
    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved averaged kymograph to: {save_path}")
    plt.show()


def plot_chain_growth_without_zero_padding(parent_dir, output_dir=DEFAULT_OUTPUT_DIR):
    def growth_model(t, B0, lambda_):
        return B0 * np.exp(lambda_ * t)

    is_ca = 'CA' in parent_dir
    repeat_dirs = sorted(glob.glob(os.path.join(parent_dir, "repeat*")))
    if not repeat_dirs:
        print("No repeat directories found.")
        return

    time_to_chains = {}
    for repeat_dir in repeat_dirs:
        files = sorted(glob.glob(os.path.join(repeat_dir, "biofilm_*.dat")),
                      key=lambda x: int(re.search(r'biofilm_(\d+)', os.path.basename(x)).group(1)))
        for file in files:
            match = re.search(r'biofilm_(\d+)', os.path.basename(file))
            if not match:
                continue
            time = int(match.group(1)) * 0.1
            data = pd.read_csv(file, sep='\t')

            if 'upper_link' not in data.columns or 'lower_link' not in data.columns:
                continue
            def is_valid_link(val):
                return not (val is None or pd.isna(val) or val == "None" or val == -1)

            G = nx.Graph()
            for idx, row in data.iterrows():
                upper = row['upper_link']
                lower = row['lower_link']
                if is_valid_link(upper):
                    G.add_edge(idx, int(upper))
                if is_valid_link(lower):
                    G.add_edge(idx, int(lower))
            # valid_data = data[(data['upper_link'] != -1) | (data['lower_link'] != -1)]

            # G = nx.Graph()
            # for idx, row in valid_data.iterrows():
            #     if row['upper_link'] != -1:
            #         G.add_edge(idx, int(row['upper_link']))
            #     if row['lower_link'] != -1:
            #         G.add_edge(idx, int(row['lower_link']))

            connected_components = list(nx.connected_components(G))
            num_chains = sum(1 for comp in connected_components if len(comp) >= 2)

            time_to_chains.setdefault(time, []).append(num_chains)

    sorted_times = sorted(time_to_chains.keys())
    avg_chains = []
    sem_chains = []
    time_points = []

    for t in sorted_times:
        vals = time_to_chains[t]
        if vals:
            avg_chains.append(np.mean(vals))
            sem_chains.append(np.std(vals, ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0)
            time_points.append(t)

    if not time_points:
        print("No valid chain data found.")
        return

    time_points = np.array(time_points)
    avg_chains = np.array(avg_chains)

    params, _ = curve_fit(growth_model, time_points, avg_chains, p0=[avg_chains[0], 0.5])
    B0_fit, lambda_fit = params
    time_doubling = np.log(2) / lambda_fit

    plt.figure(figsize=(5, 3))
    plt.errorbar(time_points, avg_chains, yerr=sem_chains, fmt='o', label="Mean ± SEM", color='red')
    plt.plot(time_points, growth_model(time_points, *params), label=f"Fit: Doubling Time = {time_doubling:.4f} hr", color='k')
    plt.xlabel("Time (h)")
    plt.ylabel("Average Number of Chains")
    plt.yscale('log')
    plt.legend()
    plt.tight_layout()

    tag = "_".join(os.path.normpath(parent_dir).split(os.sep)[-2:])
    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"chain_growth_{tag}_no_zeropadding.pdf")
    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved chain growth plot to: {save_path}")
    print(f"Fitted Growth Rate (lambda): {lambda_fit:.4f} per hour")
    print(f"Initial Number of Chains (B0): {B0_fit:.2f}")
    plt.show()

def plot_total_and_deleted_cells(data_dir, log_file_path=None, save_path=None, output_dir=DEFAULT_OUTPUT_DIR):
    import matplotlib.pyplot as plt
    from collections import defaultdict

    # Auto-detect log file if not provided
    if log_file_path is None:
        repeat_dir = os.path.dirname(data_dir)
        log_candidates = [f for f in os.listdir(repeat_dir) if f.endswith('.txt')]
        if log_candidates:
            log_file_path = os.path.join(repeat_dir, log_candidates[0])
            print(f"Detected log file: {log_file_path}")
        else:
            print(" No log file found in repeat directory.")
            return

    # Parse deleted cells from log file
    deleted_cells_by_time = defaultdict(list)
    current_time = None

    with open(log_file_path, "r") as file:
        for line in file:
            if "saving at time:" in line:
                parts = line.strip().split()
                try:
                    current_time = float(parts[3])
                except (IndexError, ValueError):
                    continue
            if "Deleting cell with ID" in line:
                match = re.search(r"Deleting cell with ID: (\d+)", line)
                if match and current_time is not None:
                    cell_id = int(match.group(1))
                    deleted_cells_by_time[current_time].append(cell_id)

    file_pattern = os.path.join(data_dir, "biofilm_*.dat")
    frame_files = sorted(glob.glob(file_pattern), key=lambda x: int(re.search(r"(\d+)", os.path.basename(x)).group(1)))

    times = []
    total_cell_counts = []
    deleted_cell_counts = []

    for frame_path in frame_files:
        match = re.search(r"(\d+)", os.path.basename(frame_path))
        if not match:
            continue
        time = int(match.group(1)) * 0.1
        df = pd.read_csv(frame_path, sep='\t')

        total_cells = len(df)
        deleted_ids = deleted_cells_by_time.get(time, [])

        times.append(time)
        total_cell_counts.append(total_cells)
        deleted_cell_counts.append(len(deleted_ids))

    # Plotting
    fig, ax1 = plt.subplots(figsize=(6, 3.5))

    color = 'tab:blue'
    ax1.set_xlabel('Time (hours)')
    ax1.set_ylabel('Total Cell Count', color=color)
    ax1.plot(times, total_cell_counts, marker='o', color=color, label="Total cells")
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()
    color = 'k'
    ax2.set_ylabel('Deleted Cell Count', color=color)
    ax2.plot(times, deleted_cell_counts, marker='x', linestyle='--', color=color, label="Deleted cells")
    ax2.tick_params(axis='y', labelcolor=color)
    fig.tight_layout()

    if save_path is None or os.path.isdir(save_path):
        parts = os.path.normpath(data_dir).split(os.sep)
        tag = "_".join(parts[-2:])
        os.makedirs(output_dir, exist_ok=True)
        save_path = os.path.join(output_dir, f"cell_vs_deleted_{tag}.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved total vs deleted cell count plot to: {save_path}")
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Automate Candida/Pseudomonas counts and snapshot plotting.")
    subparsers = parser.add_subparsers(dest="command")

    count_parser = subparsers.add_parser("count", help="Plot cell counts over time")
    count_parser.add_argument("--data_dir", required=True)
    count_parser.add_argument("--save_path", required=False)
    count_parser.add_argument("--output_dir", required=False)

    snap_parser = subparsers.add_parser("snapshots", help="Plot snapshot grid")
    snap_parser.add_argument("--parent_dir", required=True)
    snap_parser.add_argument("--output_path", required=False)
    snap_parser.add_argument("--num_snapshots", type=int, default=7)
    snap_parser.add_argument("--output_dir", required=False)

    both_parser = subparsers.add_parser("both", help="Plot counts and snapshots")
    both_parser.add_argument("--data_dir", required=True)
    both_parser.add_argument("--parent_dir", required=True)
    both_parser.add_argument("--num_snapshots", type=int, default=7)
    both_parser.add_argument("--output_dir", required=False)

    avg_parser = subparsers.add_parser("avg_counts", help="Plot averaged cell counts with SEM")
    avg_parser.add_argument("--parent_dir", required=True)
    avg_parser.add_argument("--output_dir", required=False)

    all_parser = subparsers.add_parser("all", help="Plot all: counts,growth, snapshots, averages")
    all_parser.add_argument("--data_dir", required=True)
    all_parser.add_argument("--parent_dir", required=True)
    all_parser.add_argument("--num_snapshots", type=int, default=7)
    all_parser.add_argument("--output_dir", required=False)

    args = parser.parse_args()
    output_dir = args.output_dir or DEFAULT_OUTPUT_DIR

    if args.command == "count":
        plot_counts_over_time(args.data_dir, args.save_path, output_dir)
    elif args.command == "snapshots":
        plot_cells_grid(args.parent_dir, args.output_path, args.num_snapshots, output_dir)
    elif args.command == "both":
        plot_counts_over_time(args.data_dir, None, output_dir)
        plot_cells_grid(args.parent_dir, None, args.num_snapshots, output_dir)
    elif args.command == "avg_counts":
        plot_average_counts_over_repeats(args.parent_dir, output_dir)

        
    elif args.command == "all":
        # plot_counts_over_time(args.data_dir, None, output_dir)
        # plot_total_and_deleted_cells(args.data_dir, None, output_dir)
        # plot_average_counts_over_repeats(args.parent_dir, output_dir)
        # plot_chain_growth(args.parent_dir, output_dir)
        # plot_kymograph_average(args.parent_dir, output_dir)
        compute_kymograph_average_without_zerosfilling(args.parent_dir, output_dir)
        # plot_average_counts_over_repeats_without_zerosfilling(args.parent_dir, output_dir)
        # plot_chain_growth_without_zero_padding(args.parent_dir, output_dir)
        # plot_cells_grid(args.parent_dir, None, args.num_snapshots, output_dir) 
        
    else:
        parser.print_help()

if __name__ == "__main__":
    main()

