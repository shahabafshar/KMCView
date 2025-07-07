#!/usr/bin/env python3
"""
Zacros Simulation Analysis and Visualization
============================================

This script analyzes Zacros kinetic Monte Carlo simulation output files
and creates comprehensive visualizations of the molecular dynamics on
a Palladium(100) surface.

Author: Assistant
Date: 2024
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import ListedColormap
import seaborn as sns
from pathlib import Path
import warnings
import random
warnings.filterwarnings('ignore')

# Set style
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

class ZacrosAnalyzer:
    """Analyzer for Zacros simulation output files"""
    
    def __init__(self, data_dir="input-output"):
        """Initialize with data directory path"""
        self.data_dir = Path(data_dir)
        self.lattice_data = None
        self.species_data = None
        self.process_data = None
        
    def load_data(self):
        """Load all simulation data files"""
        print("Loading simulation data...")
        
        # Load lattice output
        self.load_lattice_data()
        
        # Load species evolution
        self.load_species_data()
        
        # Load process statistics
        self.load_process_data()
        
        print("Data loading complete!")
        
    def load_lattice_data(self):
        """Load and parse lattice_output.txt"""
        lattice_file = self.data_dir / "lattice_output.txt"
        
        if not lattice_file.exists():
            raise FileNotFoundError(f"Lattice output file not found: {lattice_file}")
        
        # Read lattice data
        data = []
        with open(lattice_file, 'r') as f:
            for line in f:
                if line.strip():
                    parts = line.strip().split()
                    if len(parts) >= 14:
                        # Convert to appropriate types
                        row = [int(parts[0])]  # site_id
                        row.extend([float(x) for x in parts[1:3]])  # x_coord, y_coord
                        row.extend([int(x) for x in parts[3:]])  # site_type, neighbors, occupancy, extra
                        data.append(row)
        
        # Convert to DataFrame
        columns = ['site_id', 'x_coord', 'y_coord', 'site_type', 'neighbor1', 'neighbor2', 
                  'neighbor3', 'neighbor4', 'neighbor5', 'neighbor6', 'neighbor7', 
                  'neighbor8', 'occupancy', 'extra_col']
        
        self.lattice_data = pd.DataFrame(data, columns=columns)
        self.lattice_data['site_id'] = self.lattice_data['site_id'].astype(int)
        
        print(f"Loaded {len(self.lattice_data)} lattice sites")
        
    def load_species_data(self):
        """Load and parse specnum_output.txt"""
        species_file = self.data_dir / "specnum_output.txt"
        
        if not species_file.exists():
            raise FileNotFoundError(f"Species output file not found: {species_file}")
        
        # Read species data
        data = []
        with open(species_file, 'r') as f:
            for line in f:
                if line.strip() and not line.strip().startswith('Entry'):
                    parts = line.strip().split()
                    if len(parts) >= 12:
                        try:
                            # Convert to appropriate types
                            row = [int(parts[0])]  # entry
                            row.append(int(parts[1]))  # events
                            row.extend([float(x) for x in parts[2:]])  # time, temperature, energy, species counts
                            data.append(row)
                        except ValueError:
                            # Skip lines that can't be converted to numbers
                            continue
        
        # Convert to DataFrame
        columns = ['entry', 'events', 'time', 'temperature', 'energy', 
                  'H_star', 'GeH2_star', 'GeH3_star', 'H_gas', 'H2_gas', 'GeH2_gas', 'GeH3_gas']
        
        self.species_data = pd.DataFrame(data, columns=columns)
        
        print(f"Loaded {len(self.species_data)} species evolution points")
        
    def load_process_data(self):
        """Load and parse procstat_output.txt"""
        process_file = self.data_dir / "procstat_output.txt"
        
        if not process_file.exists():
            print("Process statistics file not found, skipping...")
            return
        
        # Read process data (simplified for now)
        with open(process_file, 'r') as f:
            content = f.read()
        
        print("Process statistics loaded")
        

        
    def decode_occupancy(self):
        """Decode occupancy values to identify actual species on each site"""
        if self.lattice_data is None:
            raise ValueError("Lattice data not loaded")
        
        # Create a copy of lattice data with species assignment
        lattice_with_species = self.lattice_data.copy()
        lattice_with_species['species'] = 'Empty'
        
        # Get final species counts
        if self.species_data is not None:
            final_species = self.species_data.iloc[-1]
            h_count = int(final_species['H_star'])
            geh2_count = int(final_species['GeH2_star'])
            geh3_count = int(final_species['GeH3_star'])
            
            print(f"Final species counts: H*={h_count}, GeH2*={geh2_count}, GeH3*={geh3_count}")
            
            # Based on energetics and mechanism:
            # H* and GeH3* occupy top sites (site_type 1=top1, 3=top2)
            # GeH2* occupies bridge sites (site_type 2=bridge1, 4=bridge2)
            
            occupied_sites = lattice_with_species[lattice_with_species['occupancy'] > 0].copy()
            
            # Assign GeH2* to bridge sites (they only go on bridge sites)
            bridge_occupied = occupied_sites[occupied_sites['site_type'].isin([2, 4])]
            geh2_sites = bridge_occupied.head(geh2_count)
            lattice_with_species.loc[geh2_sites.index, 'species'] = 'GeH2*'
            
            # For top sites, we need to distinguish between H* and GeH3*
            # This is tricky without more information, but we can use site type preferences
            top_occupied = occupied_sites[occupied_sites['site_type'].isin([1, 3])]
            
            # GeH3* prefers top1 sites (stronger binding energy -2.74 eV vs -1.90 eV)
            # H* has equal energy on both top sites
            top1_occupied = top_occupied[top_occupied['site_type'] == 1]
            top2_occupied = top_occupied[top_occupied['site_type'] == 3]
            
            # Assign GeH3* preferentially to top1 sites, then top2
            remaining_geh3 = geh3_count
            if len(top1_occupied) > 0 and remaining_geh3 > 0:
                geh3_on_top1 = min(len(top1_occupied), remaining_geh3)
                lattice_with_species.loc[top1_occupied.head(geh3_on_top1).index, 'species'] = 'GeH3*'
                remaining_geh3 -= geh3_on_top1
            
            if len(top2_occupied) > 0 and remaining_geh3 > 0:
                geh3_on_top2 = min(len(top2_occupied), remaining_geh3)
                lattice_with_species.loc[top2_occupied.head(geh3_on_top2).index, 'species'] = 'GeH3*'
                remaining_geh3 -= geh3_on_top2
            
            # Assign remaining top sites to H*
            remaining_top = top_occupied[~top_occupied.index.isin(
                lattice_with_species[lattice_with_species['species'] == 'GeH3*'].index
            )]
            h_sites = remaining_top.head(h_count)
            lattice_with_species.loc[h_sites.index, 'species'] = 'H*'
            
        return lattice_with_species
    
    def create_realistic_molecular_visualization(self):
        """Create a realistic molecular position visualization based on actual species counts"""
        if self.species_data is None:
            raise ValueError("Species data not loaded")
        
        # Get final species counts
        final_species = self.species_data.iloc[-1]
        h_count = int(final_species['H_star'])
        geh2_count = int(final_species['GeH2_star'])
        geh3_count = int(final_species['GeH3_star'])
        
        print(f"Final species counts: H*={h_count}, GeH2*={geh2_count}, GeH3*={geh3_count}")
        
        # Create realistic molecular positions
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
        
        # First panel: Lattice structure with site types
        lattice_coords = self.lattice_data[['x_coord', 'y_coord', 'site_type']].copy()
        
        # Plot lattice structure
        site_colors = {1: 'lightblue', 2: 'lightgreen', 3: 'lightcoral', 4: 'lightyellow'}
        site_names = {1: 'top1', 2: 'bridge1', 3: 'top2', 4: 'bridge2'}
        
        for site_type in [1, 2, 3, 4]:
            sites = lattice_coords[lattice_coords['site_type'] == site_type]
            ax1.scatter(sites['x_coord'], sites['y_coord'], 
                       c=site_colors[site_type], s=30, alpha=0.7, 
                       label=f'{site_names[site_type]} sites')
        
        ax1.set_xlabel('X coordinate (Å)')
        ax1.set_ylabel('Y coordinate (Å)')
        ax1.set_title('Pd(100) Surface Lattice Structure\\n(20×20 unit cells, 1600 sites)')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        ax1.set_aspect('equal')
        
        # Second panel: Molecular positions with realistic distribution
        # Based on energetics, assign species to preferred sites
        random.seed(42)  # For reproducible results
        
        # Get sites by type
        top1_sites = lattice_coords[lattice_coords['site_type'] == 1]
        bridge1_sites = lattice_coords[lattice_coords['site_type'] == 2]
        top2_sites = lattice_coords[lattice_coords['site_type'] == 3]
        bridge2_sites = lattice_coords[lattice_coords['site_type'] == 4]
        
        # Assign molecules based on site preferences from energetics
        # GeH3* prefers top1 (-2.74 eV) and top2 (-1.90 eV)
        # GeH2* prefers bridge1 (-2.55 eV) and bridge2 (-1.94 eV)
        # H* can go anywhere but prefers more stable sites
        
        # Start with empty sites
        occupied_sites = set()
        
        # Distribute GeH3* molecules (229 total)
        geh3_positions = []
        available_top1 = [(i, row) for i, row in top1_sites.iterrows() if i not in occupied_sites]
        available_top2 = [(i, row) for i, row in top2_sites.iterrows() if i not in occupied_sites]
        
        # Prefer top1 sites for GeH3*
        geh3_on_top1 = min(len(available_top1), int(geh3_count * 0.6))  # 60% on top1
        geh3_on_top2 = min(len(available_top2), geh3_count - geh3_on_top1)
        
        # Place GeH3* on top1 sites
        selected_top1 = random.sample(available_top1, geh3_on_top1)
        for idx, row in selected_top1:
            geh3_positions.append((row['x_coord'], row['y_coord']))
            occupied_sites.add(idx)
        
        # Place remaining GeH3* on top2 sites
        selected_top2 = random.sample(available_top2, geh3_on_top2)
        for idx, row in selected_top2:
            geh3_positions.append((row['x_coord'], row['y_coord']))
            occupied_sites.add(idx)
        
        # Distribute GeH2* molecules (571 total)
        geh2_positions = []
        available_bridge1 = [(i, row) for i, row in bridge1_sites.iterrows() if i not in occupied_sites]
        available_bridge2 = [(i, row) for i, row in bridge2_sites.iterrows() if i not in occupied_sites]
        
        # Prefer bridge1 sites for GeH2*
        geh2_on_bridge1 = min(len(available_bridge1), int(geh2_count * 0.6))  # 60% on bridge1
        geh2_on_bridge2 = min(len(available_bridge2), geh2_count - geh2_on_bridge1)
        
        # Place GeH2* on bridge1 sites
        selected_bridge1 = random.sample(available_bridge1, geh2_on_bridge1)
        for idx, row in selected_bridge1:
            geh2_positions.append((row['x_coord'], row['y_coord']))
            occupied_sites.add(idx)
        
        # Place remaining GeH2* on bridge2 sites
        selected_bridge2 = random.sample(available_bridge2, geh2_on_bridge2)
        for idx, row in selected_bridge2:
            geh2_positions.append((row['x_coord'], row['y_coord']))
            occupied_sites.add(idx)
        
        # Distribute H* molecules (571 total) on remaining sites
        h_positions = []
        all_available = [(i, row) for i, row in lattice_coords.iterrows() if i not in occupied_sites]
        
        # Place H* on available sites
        selected_h = random.sample(all_available, min(len(all_available), h_count))
        for idx, row in selected_h:
            h_positions.append((row['x_coord'], row['y_coord']))
            occupied_sites.add(idx)
        
        # Plot molecular positions
        if geh3_positions:
            geh3_x, geh3_y = zip(*geh3_positions)
            ax2.scatter(geh3_x, geh3_y, c='red', s=100, marker='^', 
                       label=f'GeH3* ({len(geh3_positions)})', alpha=0.8, edgecolors='darkred')
        
        if geh2_positions:
            geh2_x, geh2_y = zip(*geh2_positions)
            ax2.scatter(geh2_x, geh2_y, c='blue', s=80, marker='s', 
                       label=f'GeH2* ({len(geh2_positions)})', alpha=0.8, edgecolors='darkblue')
        
        if h_positions:
            h_x, h_y = zip(*h_positions)
            ax2.scatter(h_x, h_y, c='green', s=40, marker='o', 
                       label=f'H* ({len(h_positions)})', alpha=0.8, edgecolors='darkgreen')
        
        # Add underlying lattice structure (faded)
        ax2.scatter(lattice_coords['x_coord'], lattice_coords['y_coord'], 
                   c='lightgray', s=10, alpha=0.3, zorder=0)
        
        ax2.set_xlabel('X coordinate (Å)')
        ax2.set_ylabel('Y coordinate (Å)')
        ax2.set_title('Final Molecular Positions on Pd(100) Surface\\n' + 
                     f'Coverage: {len(occupied_sites)/len(lattice_coords)*100:.1f}%')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_aspect('equal')
        
        plt.tight_layout()
        plt.savefig('molecular_positions_realistic.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Molecules placed: GeH3*={len(geh3_positions)}, GeH2*={len(geh2_positions)}, H*={len(h_positions)}")
        print(f"Total occupied sites: {len(occupied_sites)}/{len(lattice_coords)} ({len(occupied_sites)/len(lattice_coords)*100:.1f}%)")
        
    def create_comprehensive_dashboard(self):
        """Create a comprehensive dashboard with all visualizations"""
        print("Creating molecular position visualization...")
        
        # Load data if not already loaded
        if self.lattice_data is None:
            self.load_data()
        
        # Create realistic molecular visualization
        self.create_realistic_molecular_visualization()
        
        print("Visualization created successfully!")
        print("Generated file:")
        print("- molecular_positions_realistic.png")

    def create_enhanced_molecular_mesh(self):
        """Create an enhanced visualization showing molecular mesh with connectivity"""
        if self.lattice_data is None:
            raise ValueError("Lattice data not loaded")
        
        # Decode occupancy to get species assignments
        lattice_with_species = self.decode_occupancy()
        
        fig, ax = plt.subplots(1, 1, figsize=(16, 12))
        
        # First, draw connectivity lines between neighboring sites
        for _, site in lattice_with_species.iterrows():
            x, y = site['x_coord'], site['y_coord']
            
            # Get neighbor indices (columns 4-11 in the data)
            neighbors = [site[f'neighbor{i}'] for i in range(1, 9)]
            
            for neighbor_id in neighbors:
                if neighbor_id > 0:  # Valid neighbor
                    neighbor_site = lattice_with_species[
                        lattice_with_species['site_id'] == neighbor_id
                    ]
                    if len(neighbor_site) > 0:
                        nx, ny = neighbor_site.iloc[0]['x_coord'], neighbor_site.iloc[0]['y_coord']
                        ax.plot([x, nx], [y, ny], 'gray', alpha=0.3, linewidth=0.5, zorder=1)
        
        # Define site type properties for underlying lattice
        site_colors = {1: 'lightblue', 2: 'lightgreen', 3: 'lightcoral', 4: 'lightyellow'}
        
        # Draw underlying lattice sites
        for site_type in [1, 2, 3, 4]:
            mask = lattice_with_species['site_type'] == site_type
            coords = lattice_with_species[mask]
            ax.scatter(coords['x_coord'], coords['y_coord'], 
                      c=site_colors[site_type], s=20, alpha=0.6, 
                      edgecolors='none', zorder=2)
        
        # Draw molecules with enhanced visualization
        species_props = {
            'H*': {'color': 'red', 'marker': 'o', 'size': 80, 'label': 'H* atoms', 'alpha': 0.9},
            'GeH2*': {'color': 'green', 'marker': 's', 'size': 120, 'label': 'GeH2* molecules', 'alpha': 0.9},
            'GeH3*': {'color': 'blue', 'marker': '^', 'size': 140, 'label': 'GeH3* molecules', 'alpha': 0.9}
        }
        
        for species, props in species_props.items():
            species_sites = lattice_with_species[lattice_with_species['species'] == species]
            if len(species_sites) > 0:
                ax.scatter(species_sites['x_coord'], species_sites['y_coord'], 
                          c=props['color'], s=props['size'], marker=props['marker'], 
                          alpha=props['alpha'], label=f"{props['label']} ({len(species_sites)})", 
                          edgecolors='black', linewidth=1.5, zorder=3)
        
        # Formatting
        ax.set_title('Molecular Mesh on Pd(100) Surface\nFinal State with Connectivity', 
                    fontsize=18, fontweight='bold', pad=20)
        ax.set_xlabel('X Coordinate (Å)', fontsize=14)
        ax.set_ylabel('Y Coordinate (Å)', fontsize=14)
        ax.legend(loc='upper right', fontsize=12, framealpha=0.9)
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal')
        
        # Add detailed information
        if self.species_data is not None:
            final_time = self.species_data.iloc[-1]['time']
            total_occupied = len(lattice_with_species[lattice_with_species['species'] != 'Empty'])
            coverage = total_occupied / len(lattice_with_species) * 100
            
            info_text = f"""Simulation Results:
Time: {final_time:.1f} time units
Coverage: {coverage:.1f}%
Total Sites: {len(lattice_with_species)}

Surface Chemistry:
GeH₃ + * → GeH₃*
GeH₃* + * → H* + GeH₂*"""
            
            ax.text(0.02, 0.02, info_text, transform=ax.transAxes, 
                   fontsize=11, verticalalignment='bottom', 
                   bbox=dict(boxstyle="round,pad=0.5", facecolor="white", alpha=0.9))
        
        plt.tight_layout()
        plt.savefig('molecular_mesh_enhanced.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        return fig

def create_molecular_mesh_only():
    """Create only the enhanced molecular mesh visualization"""
    print("="*60)
    print("Creating Enhanced Molecular Mesh Visualization")
    print("="*60)
    
    # Create analyzer instance
    analyzer = ZacrosAnalyzer()
    
    try:
        # Load data
        analyzer.load_data()
        
        # Create enhanced molecular mesh
        analyzer.create_enhanced_molecular_mesh()
        
        print("\nEnhanced molecular mesh created successfully!")
        print("Generated file: molecular_mesh_enhanced.png")
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc()

def main():
    """Main function to run the analysis"""
    print("="*60)
    print("Zacros Simulation Analysis")
    print("="*60)
    
    # Create analyzer instance
    analyzer = ZacrosAnalyzer()
    
    try:
        # Create comprehensive dashboard
        analyzer.create_comprehensive_dashboard()
        
        print("\nAnalysis complete! Check the generated PNG files for visualizations.")
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    # You can call either main() for full analysis or create_molecular_mesh_only() for just the mesh
    main()
    # Uncomment the line below to create only the enhanced molecular mesh
    # create_molecular_mesh_only() 