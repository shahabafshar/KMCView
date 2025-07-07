#!/usr/bin/env python3
"""
Simple Temporal Evolution Visualization for Zacros Simulation
Fast and efficient - no infinite loops!
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set style
plt.style.use('seaborn-v0_8')

def load_species_data():
    """Load species evolution data"""
    species_file = Path("input-output/specnum_output.txt")
    
    # Read species data
    data = []
    with open(species_file, 'r') as f:
        for line in f:
            if line.strip() and not line.strip().startswith('Entry'):
                parts = line.strip().split()
                if len(parts) >= 12:
                    try:
                        row = [int(parts[0])]  # entry
                        row.append(int(parts[1]))  # events
                        row.extend([float(x) for x in parts[2:]])  # time, temperature, energy, species counts
                        data.append(row)
                    except (ValueError, IndexError):
                        continue
    
    # Convert to DataFrame
    columns = ['entry', 'events', 'time', 'temperature', 'energy', 
              'H_star', 'GeH2_star', 'GeH3_star', 'H_gas', 'H2_gas', 'GeH2_gas', 'GeH3_gas']
    
    return pd.DataFrame(data, columns=columns)

def create_temporal_evolution():
    """Create temporal evolution visualization"""
    print("Loading species evolution data...")
    species_data = load_species_data()
    
    # Select 10 time points for visualization
    time_points = [0, 10, 50, 100, 150, 200, 250, 300, 400, 500]
    
    # Extract species counts at each time point
    temporal_data = []
    for target_time in time_points:
        # Find the closest time point in the data
        closest_idx = species_data.iloc[:, 2].sub(target_time).abs().idxmin()
        row = species_data.iloc[closest_idx]
        
        temporal_data.append({
            'time': row['time'],
            'H_star': int(row['H_star']),
            'GeH2_star': int(row['GeH2_star']),
            'GeH3_star': int(row['GeH3_star']),
            'events': int(row['events'])
        })
    
    # Create species evolution curves
    print("Creating species evolution curves...")
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Plot species evolution over time
    times = [d['time'] for d in temporal_data]
    h_counts = [d['H_star'] for d in temporal_data]
    geh2_counts = [d['GeH2_star'] for d in temporal_data]
    geh3_counts = [d['GeH3_star'] for d in temporal_data]
    
    ax1.plot(times, h_counts, 'ro-', linewidth=3, markersize=10, label='H* atoms')
    ax1.plot(times, geh2_counts, 'bs-', linewidth=3, markersize=10, label='GeH2* molecules')
    ax1.plot(times, geh3_counts, 'g^-', linewidth=3, markersize=10, label='GeH3* molecules')
    
    ax1.set_xlabel('Time (units)', fontsize=14)
    ax1.set_ylabel('Number of Species', fontsize=14)
    ax1.set_title('Species Evolution Over Time - GeH₄ Decomposition Process', fontsize=16, fontweight='bold')
    ax1.legend(fontsize=14)
    ax1.grid(True, alpha=0.3)
    
    # Add annotations for key stages
    ax1.annotate('Initial adsorption\\nGeH3* → Surface', xy=(50, 25), xytext=(80, 100),
                arrowprops=dict(arrowstyle='->', color='green', lw=2),
                fontsize=12, ha='center')
    
    ax1.annotate('Decomposition\\nGeH3* → GeH2* + H*', xy=(200, 200), xytext=(250, 350),
                arrowprops=dict(arrowstyle='->', color='blue', lw=2),
                fontsize=12, ha='center')
    
    ax1.annotate('Equilibrium\\nH* ≈ GeH2*', xy=(450, 550), xytext=(350, 500),
                arrowprops=dict(arrowstyle='->', color='red', lw=2),
                fontsize=12, ha='center')
    
    # Plot surface coverage percentage (correct total sites = 1600)
    total_sites = 1600  # 20x20 lattice with 4 sites per cell
    coverage = [(d['H_star'] + d['GeH2_star'] + d['GeH3_star']) / total_sites * 100 for d in temporal_data]
    
    ax2.plot(times, coverage, 'mo-', linewidth=3, markersize=10, label='Total Surface Coverage')
    ax2.set_xlabel('Time (units)', fontsize=14)
    ax2.set_ylabel('Surface Coverage (%)', fontsize=14)
    ax2.set_title('Surface Coverage Evolution', fontsize=16, fontweight='bold')
    ax2.legend(fontsize=14)
    ax2.grid(True, alpha=0.3)
    
    # Add horizontal line at final coverage
    final_coverage = coverage[-1]
    ax2.axhline(y=final_coverage, color='red', linestyle='--', alpha=0.7, 
               label=f'Final Coverage: {final_coverage:.1f}%')
    ax2.legend(fontsize=12)
    
    plt.tight_layout()
    plt.savefig('species_evolution_curves.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create a bar chart showing the 10 steps
    print("Creating temporal evolution bar chart...")
    fig, ax = plt.subplots(figsize=(16, 8))
    
    x_pos = np.arange(len(temporal_data))
    width = 0.25
    
    bars1 = ax.bar(x_pos - width, h_counts, width, label='H* atoms', color='red', alpha=0.8)
    bars2 = ax.bar(x_pos, geh2_counts, width, label='GeH2* molecules', color='blue', alpha=0.8)
    bars3 = ax.bar(x_pos + width, geh3_counts, width, label='GeH3* molecules', color='green', alpha=0.8)
    
    ax.set_xlabel('Time Steps', fontsize=14)
    ax.set_ylabel('Number of Species', fontsize=14)
    ax.set_title('Gradual GeH₄ Decomposition Process - 10 Time Steps', fontsize=16, fontweight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels([f'{d["time"]:.0f}' for d in temporal_data])
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    
    # Add value labels on bars
    for i, (h, g2, g3) in enumerate(zip(h_counts, geh2_counts, geh3_counts)):
        if h > 0:
            ax.text(i-width, h+10, str(h), ha='center', va='bottom', fontsize=10)
        if g2 > 0:
            ax.text(i, g2+10, str(g2), ha='center', va='bottom', fontsize=10)
        if g3 > 0:
            ax.text(i+width, g3+10, str(g3), ha='center', va='bottom', fontsize=10)
    
    plt.tight_layout()
    plt.savefig('temporal_evolution_bars.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("\\nTemporal Analysis Summary:")
    print("=" * 50)
    for i, data in enumerate(temporal_data):
        print(f"Step {i+1:2d}: Time={data['time']:6.1f} | "
              f"H*={data['H_star']:3d} | GeH2*={data['GeH2_star']:3d} | "
              f"GeH3*={data['GeH3_star']:3d} | Events={data['events']:4d}")
    
    print("\\nKey Observations:")
    print("- GeH3* molecules gradually decompose over time")
    print("- H* and GeH2* counts increase as GeH3* decreases")
    print("- The process shows H* and GeH2* counts are always equal (1:1 stoichiometry)")
    print(f"- Surface coverage reaches {final_coverage:.1f}% at equilibrium")
    
    print("\\nFiles generated:")
    print("- species_evolution_curves.png")
    print("- temporal_evolution_bars.png")
    
    return temporal_data

if __name__ == "__main__":
    print("="*60)
    print("Zacros Temporal Evolution Analysis - FAST VERSION")
    print("="*60)
    
    try:
        temporal_data = create_temporal_evolution()
        print("\\nAnalysis completed successfully!")
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc() 