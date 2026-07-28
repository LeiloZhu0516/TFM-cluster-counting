# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 18:02:25 2026

@author: Pcx
"""
import os
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
import linecache
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LogNorm
import scienceplots

# =============================================================================
# Matplotlib Configuration
# =============================================================================
plt.style.use(['science', 'ieee'])
plt.rcParams['text.usetex'] = False 
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'Times', 'DejaVu Serif']

# =============================================================================
# Helper Functions
# =============================================================================
def analyze_clusters(df, L=40.0, d=1.3, tolerance=0.7):
    #df contains the positions of all particles of the system at a specific time
    #d is the maximum distance between two particles to consider them aggregated
    #tolerance is to ensure the aggregates are chains aligned with the field, 0.7 is approx sin(45)
    
    # 1. Extract coordinates and force them into the [0, L] range
    # This handles negative values and values > L
    coords = df[['x', 'y', 'z']].astype(float).values % L 
    
    # 2. The KDTree will accept the data
    tree = cKDTree(coords, boxsize=L)   # Reads the whole system and enforces PBC
    pairs = tree.query_pairs(r=d)       # Find all pairs within distance d
    
    # 3. Disjoint Set Union (DSU) to group clusters
    parent = list(range(len(df)))   # Assign a cluster ID to each particle, all particles within the same cluster share the same cluster ID (root)
    size = [1] * len(df)            # Tracks the size of each cluster
    degree = [0] * len(df)          # Tracks the number of bonds of each particle
    
    def find(i): # Looks for the cluster ID (root) of the particle i
        root = i
        while parent[root] != root: 
            root = parent[root] 
        # Path compression: wire everything directly to the root
        while parent[i] != root: 
            next_node = parent[i]
            parent[i] = root
            i = next_node
        return root

    def union(root_i, root_j): # Union by Size logic (avoids linked lists)
        if root_i != root_j:
            if size[root_i] < size[root_j]:
                parent[root_i] = root_j
                size[root_j] += size[root_i]
            else:
                parent[root_j] = root_i
                size[root_i] += size[root_j]

    for i, j in pairs: # Puts all particles whithin the same cluster under teh same cluster ID
        # First we ensure that the visualized pairs can form horizontal chains
        dy = coords[i,1] - coords[j,1]
        dy = dy - L * np.round(dy / L)
        if abs(dy) > tolerance: 
            continue
        
        root_i = find(i)
        root_j = find(j)
        # Connect particles if they are not already connected and they have degree<2
        if root_i != root_j:
            if degree[i] < 2 and degree[j] < 2:
                union(root_i, root_j)
                degree[i] += 1
                degree[j] += 1
                
    # 4. Extract Clusters
    clusters = {}
    for i in range(len(df)):
        root = find(i)
        if root not in clusters: 
            clusters[root] = []
        clusters[root].append(i)
    
    chain_list = list(clusters.values())
    
    # 5. Statistics
    sizes = [len(c) for c in chain_list]
    counts = pd.Series(sizes).value_counts().sort_index()
    avg_size = np.mean(sizes)
    
    return chain_list, avg_size, counts

def get_raddist(df, L, N, dr=0.1, max_d=6.5):
    # df contains the positions of all particles of the system at a specific time
    # L is the size of the simulation box
    # N is the total number of particles
    # dr is the width of each bar of the histogram
    # max_d is the maximum distance considered
    coords = df[['x', 'y']].astype(float).values % L
    rho = N / (L**2)
    
    # 1. Build the tree with periodic boundaries
    tree = cKDTree(coords, boxsize=[L, L])
    # 2. Find pairs within max_d
    pairs = tree.query_pairs(r=max_d, output_type='ndarray')
    
    if len(pairs) == 0:
        return np.zeros(int(max_d/dr)), np.arange(dr/2, max_d, dr)
    
    # 3. Calculate distances
    p1 = coords[pairs[:, 0]]
    p2 = coords[pairs[:, 1]]
    dx = p1[:, 0] - p2[:, 0]
    dy = p1[:, 1] - p2[:, 1]
    dx = dx - L * np.round(dx / L)
    dy = dy - L * np.round(dy / L)
    distances = np.sqrt(dx**2 + dy**2)
    
    # 4. Bin the results
    bins = np.arange(0, max_d + dr, dr)
    hist, bin_edges = np.histogram(distances, bins=bins)
    
    # 5. Correct Normalization
    # Each distance is counted once by query_pairs, so multiply by 2 (symmetry)
    shell_counts = hist * 2
    # Calculate area of each ring: pi * (r_outer^2 - r_inner^2)
    r_inner = bin_edges[:-1]
    r_outer = bin_edges[1:]
    shell_areas = np.pi * (r_outer**2 - r_inner**2)
    # g(r) = (Number of pairs in shell) / (N * density * area of shell)
    g_r = shell_counts / (N * rho * shell_areas)
    radii = (r_inner + r_outer) / 2
    
    return g_r, radii

# =============================================================================
# Main Execution
# =============================================================================
if __name__ == "__main__":
    # --- 1. Define Input Parameters ---
    Gammalist = [15.68]
    Pelist = [2, 3, 4, 5, 6, 40]
    
    # Toggle for dynamic tauR vs explicit tauR list
    use_dynamic_tauR = True  
    
    k_multiplier_list = np.arange(0.1, 3.01, 0.1) # Used if use_dynamic_tauR = True
    tauR_explicit_list = [0.1, 0.5, 1.0, 2.0]     # Used if use_dynamic_tauR = False

    # --- 2. Toggles for Specific Case Graphs ---
    plot_flags = {
        'size_evol': False,
        'size_dist': False,
        'rad_dist': False,
        'energy_evol': False
    }

    # --- 3. Toggles for Aggregate Graphs ---
    aggregate_flags = {
        'N_vs_tauR_tauNR': True,
        'N_vs_Gamma': True,  
        'N_vs_Pe': True,
        'std_vs_tauR_tauNR': True,
        'heatmap_Pe_Gamma': False
    }

    window_size = 21 #size of the time window around which the average is made to smoothen the data
    steadystate_start = 300 #a rough estimate of the time step from which the system has entered steady-state

    # --- Directory Setup ---
    directories = ['size evolution graph', 'size distribution graph', 'radial distribution', 'energy graph', 'aggregate graphs']
    for d in directories:
        os.makedirs(d, exist_ok=True)

    # --- Data Storage ---
    results_data = []

    # Iterate over Pe and Gamma
    for Pe in Pelist:
        for gamma in Gammalist:
            tau_NR = (gamma**(1/3)) / Pe
            
            # Select the third parameter list based on the configuration mode
            third_params = k_multiplier_list if use_dynamic_tauR else tauR_explicit_list
            
            for val in third_params:
                if use_dynamic_tauR:
                    k = val
                    tauR = k * tau_NR
                else:
                    tauR = val
                    k = tauR / tau_NR 
                
                myPestr = f"{Pe:.2f}".replace('.', 'p')
                mygammastr = f"{gamma:.2f}".replace('.', 'p')
                mytauRstr = f'{tauR:.2f}'.replace('.', 'p')
                
                filename = f"dump_gamma_{mygammastr}_Pe_{myPestr}_tauR_{mytauRstr}.lammpstrj"
                
                try:
                    N = int(linecache.getline(filename, 4))
                    line = linecache.getline(filename, 6)
                    a, b = map(float, line.split())
                    L = b - a
                    sizelist = []
                    xtotal = 0 #counts the total time step
                    
                    total_counts = pd.Series(dtype=float)
                    num_steady_frames = 0
                    # Load Data
                    my_headers = ['type', 'x', 'y', 'z']
                    reader = pd.read_csv(filename, sep=r'\s+', chunksize=N+9, names=my_headers, usecols=range(len(my_headers)))
                except FileNotFoundError:
                    print(f"Skipping {filename} (Not found)")
                    continue
                
                for chunk in reader:
                    df = chunk.iloc[9:].copy()
                    clusters, avg_size, counts = analyze_clusters(df, L)
                    sizelist.append(avg_size)
                    
                    if xtotal >= steadystate_start:
                        total_counts = total_counts.add(counts, fill_value=0)
                        num_steady_frames += 1
                    xtotal += 1
                #averages the counts of clusters over the whole steady state
                avg_counts = (total_counts / num_steady_frames) if num_steady_frames > 0 else total_counts
                
                steady_state_sizes = sizelist[steadystate_start:]
                if len(steady_state_sizes) > 0:
                    big_avg_size = np.mean(steady_state_sizes)
                    std_size = np.std(steady_state_sizes)
                else:
                    big_avg_size, std_size = 0, 0
                    
                results_data.append({
                    'Pe': Pe,
                    'Gamma': gamma,
                    'tauR': tauR,
                    'tauNR': tau_NR,
                    'k_multiplier': k,
                    'tauR_ratio': k,
                    'N_mean': big_avg_size,
                    'N_std': std_size
                })
                
                # --- 4. Plot Specific Case Graphs ---
                if plot_flags['size_evol']:
                    x = np.arange(xtotal)
                    weights = np.ones(window_size) / window_size
                    y_smooth = np.convolve(sizelist, weights, mode='valid')
                    plt.figure()
                    plt.plot(x[window_size-1:], y_smooth)
                    plt.axhline(big_avg_size, color='r', linestyle='--', label=f"N={big_avg_size:.2f}")
                    plt.ylabel('Average size')
                    plt.xlabel('time')
                    plt.title(rf"Average size ($Pe$={Pe:.2f}, $\Gamma$={gamma:.2f}, $\tau_R$={tauR:.2f})")
                    plt.legend()
                    plt.savefig(f'size evolution graph/N_gamma_{mygammastr}_Pe_{myPestr}_tauR_{mytauRstr}.pdf', dpi=300, bbox_inches='tight')
                    plt.close()

                if plot_flags['size_dist']:
                    avg_counts_norm = avg_counts / N
                    plt.figure()
                    avg_counts_norm.plot(kind='bar', color='skyblue', width=1.00, edgecolor='black', linewidth=0.5)
                    plt.title(rf'Chain size distribution ($Pe$={Pe:.2f}, $\Gamma$={gamma:.2f})')
                    plt.xlabel('Size of the chain')
                    plt.ylabel('Counts')
                    plt.xticks(rotation=45)
                    if avg_counts_norm.index.max() >= 30:
                        step_size = 5  
                        positions = range(0, len(avg_counts_norm), step_size)
                        labels = avg_counts_norm.index[positions]
                        plt.xticks(positions, labels, rotation=45)
                    plt.tight_layout()
                    plt.savefig(f'size distribution graph/N_gamma_{mygammastr}_Pe_{myPestr}_tauR_{mytauRstr}.pdf', dpi=300, bbox_inches='tight')
                    plt.close()

                if plot_flags['rad_dist']:
                    raddist_avg, radii = get_raddist(df, L, N)
                    plt.figure()
                    plt.plot(radii, raddist_avg)
                    plt.title(rf"Radial distribution ($Pe$={Pe:.2f}, $\Gamma$={gamma:.2f}, $\tau_R$={tauR:.2f})")
                    plt.xlabel('Distance')
                    plt.ylabel('Counts')
                    plt.savefig(f'radial distribution/raddist_gamma_{mygammastr}_Pe_{myPestr}_tauR_{mytauRstr}.pdf', dpi=300, bbox_inches='tight')
                    plt.close()
                    
                if plot_flags['energy_evol']:
                    energyfile = f"energies_gamma_{mygammastr}_Pe_{myPestr}.data"
                    try:
                        energy_headers = ['time', 'Etotal', 'PE', 'KE']
                        energy_reader = pd.read_csv(energyfile, sep=r'\s+', skiprows=2, names=energy_headers, usecols=range(len(energy_headers)))
                        avg_E = energy_reader['Etotal'].iloc[300:].mean()
                        weights = np.ones(window_size) / window_size
                        E_smooth = np.convolve(energy_reader['Etotal'], weights, mode='valid')
                        
                        x = np.arange(len(energy_reader['Etotal']))
                        plt.figure()
                        plt.plot(x[window_size-1:], E_smooth)
                        plt.axhline(avg_E, color='r', linestyle='--', label=f"E={avg_E:.2f}")
                        plt.ylabel('Total Energy')
                        plt.xlabel('time')
                        plt.title(rf"Total energy with $Pe$={Pe:.2f}")
                        plt.legend()
                        plt.savefig(f'energy graph/energies_gamma_{mygammastr}_Pe_{myPestr}.pdf', dpi=300, bbox_inches='tight')
                        plt.close()
                    except FileNotFoundError:
                        print(f"Warning: {energyfile} not found. Skipping energy plot.")
                
                linecache.clearcache()

    # =============================================================================
    # 5. Aggregate Graphs 
    # =============================================================================
    df_res = pd.DataFrame(results_data)
    
    if df_res.empty:
        print("No data processed. Ensure your files exist and match the naming format.")
        exit()

    cmap = plt.get_cmap('tab10')
    
    # Determine what to group by for secondary axes
    group_var = 'k_multiplier' if use_dynamic_tauR else 'tauR'
    group_label_str = r'\tau_R/\tau_{esc}' if use_dynamic_tauR else r'\tau_R'

    if aggregate_flags['N_vs_tauR_tauNR']:
        plt.figure()
        groups = df_res.groupby(['Pe', 'Gamma'])
        for i, ((pe_val, gamma_val), group) in enumerate(groups):
            group = group.sort_values('tauR_ratio')
            plt.errorbar(group['tauR_ratio'], group['N_mean'], yerr=group['N_std'], 
                         fmt='-o', color=cmap(i % 10), capsize=3, markersize=3, 
                         label=rf'$Pe$={pe_val}, $\Gamma$={gamma_val}')
            
        plt.ylabel(r'Steady state mean chain size ($\langle N \rangle$)')
        plt.xlabel(r'Characteristic time ($\tau_R/\tau_{esc}$)')
        plt.legend(ncol=2, bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.savefig('aggregate graphs/N_vs_tauNR.pdf', dpi=300, bbox_inches='tight')
        plt.close()
        
    if aggregate_flags['std_vs_tauR_tauNR']:
        plt.figure()
        groups = df_res.groupby(['Pe', 'Gamma'])
        for i, ((pe_val, gamma_val), group) in enumerate(groups):
            group = group.sort_values('tauR_ratio')
            plt.plot(group['tauR_ratio'], group['N_std'], '-o', color=cmap(i % 10),
                     markersize=3, label=rf'$Pe$={pe_val}, $\Gamma$={gamma_val}')
        plt.ylabel('std')
        plt.xlabel(r'$\tau_R/\tau_{esc}$')
        plt.legend(title="Values", ncol=2, bbox_to_anchor=(1.05, 1), loc='upper left', frameon=True)
        plt.savefig('aggregate graphs/std_vs_tauNR.pdf', dpi=300, bbox_inches='tight')
        plt.close()

    if aggregate_flags['N_vs_Gamma']:
        plt.figure()
        groups = df_res.groupby(['Pe', group_var])
        for i, ((pe_val, g_val), group) in enumerate(groups):
            group = group.sort_values('Gamma')
            plt.errorbar(group['Gamma'], group['N_mean'], yerr=group['N_std'], 
                         fmt='-o', color=cmap(i % 10), capsize=3, markersize=3, 
                         label=rf'$Pe$={pe_val}, ${group_label_str}$={g_val:.2f}')
        plt.xlabel(r'$\Gamma$')
        plt.ylabel(r'$\langle N \rangle$')
        plt.legend(ncol=2, bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.savefig('aggregate graphs/N_vs_Gamma.pdf', dpi=300, bbox_inches='tight')
        plt.close()

    if aggregate_flags['N_vs_Pe']:
        plt.figure()
        groups = df_res.groupby(['Gamma', group_var])
        for i, ((gamma_val, g_val), group) in enumerate(groups):
            group = group.sort_values('Pe')
            plt.errorbar(group['Pe'], group['N_mean'], yerr=group['N_std'], 
                         fmt='-o', color=cmap(i % 10), capsize=3, markersize=3, 
                         label=rf'$\Gamma$={gamma_val}, ${group_label_str}$={g_val:.2f}')
        plt.xlabel(r'$Pe$')
        plt.ylabel(r'$\langle N \rangle$')
        plt.legend(ncol=2, bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.savefig('aggregate graphs/N_vs_Pe.pdf', dpi=300, bbox_inches='tight')
        plt.close()

    if aggregate_flags['heatmap_Pe_Gamma']:
        # Isolating the first element of the active third parameter list for the heatmap
        fixed_val = third_params[0]
        subset = df_res[df_res[group_var] == fixed_val]
        
        if not subset.empty:
            pivot_table = subset.pivot(index='Pe', columns='Gamma', values='N_mean')
            
            plt.figure()
            ax = sns.heatmap(pivot_table, cbar=True, square=True, cmap="inferno", norm=LogNorm(vmin=1, vmax=20))
            ax.invert_yaxis()
            plt.ylabel(r'$Pe$')
            plt.xlabel(r'$\Gamma$')
            plt.title(rf'Mean Size (Log Scale) at ${group_label_str}$={fixed_val:.2f}')
            plt.savefig('aggregate graphs/Heatmap.pdf', dpi=300, bbox_inches='tight')
            plt.close()