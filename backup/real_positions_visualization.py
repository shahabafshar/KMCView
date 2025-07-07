#!/usr/bin/env python3
"""
Create visualization using ACTUAL molecular positions from Zacros simulation
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.patches as patches

def extract_actual_occupancy():
    """Extract actual molecular occupancy from restart.inf"""
    print("🔍 Extracting ACTUAL molecular occupancy from restart.inf...")
    
    occupancy = {}
    
    with open('input-output/restart.inf', 'r') as f:
        lines = f.readlines()
    
    # The occupancy data is around lines 2349-3516 based on our search
    # Format: site_id species_id count
    for i in range(2300, 3600):  # Search in this range
        if i < len(lines):
            line = lines[i].strip()
            if line:
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        site_id = int(parts[0])
                        species_id = int(parts[1])
                        count = int(parts[2])
                        
                        # Only store occupied sites (count > 0)
                        if count > 0 and 1 <= site_id <= 1600 and 1 <= species_id <= 3:
                            occupancy[site_id] = {
                                'species_id': species_id,
                                'count': count
                            }
                    except (ValueError, IndexError):
                        continue
    
    print(f"✅ Found {len(occupancy)} occupied sites")
    return occupancy

def extract_lattice_coordinates():
    """Extract lattice coordinates from lattice_output.txt"""
    print("📍 Loading lattice coordinates from lattice_output.txt...")
    
    coordinates = {}
    
    with open('input-output/lattice_output.txt', 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split()
                if len(parts) >= 13:
                    try:
                        site_id = int(parts[0])
                        x_coord = float(parts[1])
                        y_coord = float(parts[2])
                        site_type = int(parts[3])
                        
                        coordinates[site_id] = {
                            'x': x_coord,
                            'y': y_coord,
                            'site_type': site_type
                        }
                    except (ValueError, IndexError):
                        continue
    
    print(f"✅ Loaded coordinates for {len(coordinates)} lattice sites")
    return coordinates

def create_real_molecular_visualization():
    """Create visualization using actual molecular positions"""
    
    print("🎨 Creating REAL molecular position visualization...")
    
    # Get the actual data
    occupancy = extract_actual_occupancy()
    coordinates = extract_lattice_coordinates()
    
    if not occupancy:
        print("❌ No occupancy data found!")
        return
    
    # Combine occupancy with coordinates
    occupied_sites = []
    empty_sites = []
    
    for site_id, coord_data in coordinates.items():
        if site_id in occupancy:
            # Occupied site
            occupied_sites.append({
                'site_id': site_id,
                'x': coord_data['x'],
                'y': coord_data['y'],
                'site_type': coord_data['site_type'],
                'species_id': occupancy[site_id]['species_id'],
                'count': occupancy[site_id]['count']
            })
        else:
            # Empty site
            empty_sites.append({
                'site_id': site_id,
                'x': coord_data['x'],
                'y': coord_data['y'],
                'site_type': coord_data['site_type'],
                'species_id': 0,
                'count': 0
            })
    
    print(f"📊 Found {len(occupied_sites)} occupied sites and {len(empty_sites)} empty sites")
    
    # Create the visualization
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(20, 16))
    
    # Define colors and names
    species_colors = {0: 'lightgray', 1: 'red', 2: 'blue', 3: 'green'}
    species_names = {0: 'Empty', 1: 'H*', 2: 'GeH2*', 3: 'GeH3*'}
    site_colors = {1: 'lightblue', 2: 'lightgreen', 3: 'lightcoral', 4: 'lightyellow'}
    site_names = {1: 'top1', 2: 'bridge1', 3: 'top2', 4: 'bridge2'}
    
    # Plot 1: All sites showing occupancy
    print("Creating Plot 1: All sites with occupancy...")
    
    # Plot empty sites first (background)
    if empty_sites:
        x_empty = [site['x'] for site in empty_sites]
        y_empty = [site['y'] for site in empty_sites]
        ax1.scatter(x_empty, y_empty, c='lightgray', s=8, alpha=0.3, label=f'Empty ({len(empty_sites)})')
    
    # Plot occupied sites by species
    species_counts = {}
    for species_id in [1, 2, 3]:
        sites = [site for site in occupied_sites if site['species_id'] == species_id]
        if sites:
            x_coords = [site['x'] for site in sites]
            y_coords = [site['y'] for site in sites]
            species_counts[species_id] = len(sites)
            
            ax1.scatter(x_coords, y_coords, 
                       c=species_colors[species_id],
                       label=f'{species_names[species_id]} ({len(sites)})',
                       s=25, alpha=0.8)
    
    ax1.set_xlabel('X Coordinate (Å)', fontsize=12)
    ax1.set_ylabel('Y Coordinate (Å)', fontsize=12)
    ax1.set_title('REAL Molecular Positions\\n(Actual simulation data)', fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_aspect('equal')
    
    # Plot 2: Only occupied sites (zoomed focus)
    print("Creating Plot 2: Occupied sites only...")
    
    for species_id in [1, 2, 3]:
        sites = [site for site in occupied_sites if site['species_id'] == species_id]
        if sites:
            x_coords = [site['x'] for site in sites]
            y_coords = [site['y'] for site in sites]
            
            ax2.scatter(x_coords, y_coords, 
                       c=species_colors[species_id],
                       label=f'{species_names[species_id]} ({len(sites)})',
                       s=40, alpha=0.8)
    
    ax2.set_xlabel('X Coordinate (Å)', fontsize=12)
    ax2.set_ylabel('Y Coordinate (Å)', fontsize=12)
    ax2.set_title('Occupied Sites Only\\n(Molecular distribution)', fontsize=14, fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_aspect('equal')
    
    # Plot 3: Site types distribution
    print("Creating Plot 3: Site types...")
    
    site_type_counts = {}
    for site_type in [1, 2, 3, 4]:
        sites = [site for site in occupied_sites if site['site_type'] == site_type]
        if sites:
            x_coords = [site['x'] for site in sites]
            y_coords = [site['y'] for site in sites]
            site_type_counts[site_type] = len(sites)
            
            ax3.scatter(x_coords, y_coords,
                       c=site_colors[site_type],
                       label=f'{site_names[site_type]} ({len(sites)})',
                       s=30, alpha=0.8)
    
    ax3.set_xlabel('X Coordinate (Å)', fontsize=12)
    ax3.set_ylabel('Y Coordinate (Å)', fontsize=12)
    ax3.set_title('Site Types with Molecules\\n(Surface structure preference)', fontsize=14, fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.set_aspect('equal')
    
    # Plot 4: Comparison random vs real
    print("Creating Plot 4: Random vs Real comparison...")
    
    # Show what random would look like
    np.random.seed(42)
    n_total = len(occupied_sites)
    n_h = species_counts.get(1, 0)
    n_geh2 = species_counts.get(2, 0)
    n_geh3 = species_counts.get(3, 0)
    
    # Create random positions
    all_coords = [(site['x'], site['y']) for site in coordinates.values()]
    np.random.shuffle(all_coords)
    
    # Assign species randomly
    random_h = all_coords[:n_h]
    random_geh2 = all_coords[n_h:n_h+n_geh2]
    random_geh3 = all_coords[n_h+n_geh2:n_h+n_geh2+n_geh3]
    
    # Plot random (left side)
    if random_h:
        ax4.scatter([x for x, y in random_h], [y for x, y in random_h], 
                   c='pink', s=15, alpha=0.6, label=f'Random H* ({len(random_h)})')
    if random_geh2:
        ax4.scatter([x for x, y in random_geh2], [y for x, y in random_geh2], 
                   c='lightblue', s=15, alpha=0.6, label=f'Random GeH2* ({len(random_geh2)})')
    if random_geh3:
        ax4.scatter([x for x, y in random_geh3], [y for x, y in random_geh3], 
                   c='lightgreen', s=15, alpha=0.6, label=f'Random GeH3* ({len(random_geh3)})')
    
    # Plot real (right side, shifted)
    x_offset = 20  # Shift real positions to the right
    for species_id in [1, 2, 3]:
        sites = [site for site in occupied_sites if site['species_id'] == species_id]
        if sites:
            x_coords = [site['x'] + x_offset for site in sites]
            y_coords = [site['y'] for site in sites]
            
            ax4.scatter(x_coords, y_coords, 
                       c=species_colors[species_id],
                       label=f'Real {species_names[species_id]} ({len(sites)})',
                       s=20, alpha=0.8)
    
    # Add dividing line
    ax4.axvline(x=x_offset/2, color='black', linestyle='--', alpha=0.5)
    ax4.text(x_offset/4, max([site['y'] for site in occupied_sites])*0.9, 
             'RANDOM\\n(What I did before)', ha='center', fontsize=10, 
             bbox=dict(boxstyle='round,pad=0.3', facecolor='pink', alpha=0.7))
    ax4.text(x_offset*3/4, max([site['y'] for site in occupied_sites])*0.9, 
             'REAL\\n(Actual simulation)', ha='center', fontsize=10,
             bbox=dict(boxstyle='round,pad=0.3', facecolor='lightgreen', alpha=0.7))
    
    ax4.set_xlabel('X Coordinate (Å)', fontsize=12)
    ax4.set_ylabel('Y Coordinate (Å)', fontsize=12)
    ax4.set_title('Random vs Real Molecular Positions\\n(Comparison)', fontsize=14, fontweight='bold')
    ax4.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('real_molecular_positions.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("✅ REAL molecular positions saved as real_molecular_positions.png")
    
    # Print summary
    print("\\n" + "="*50)
    print("📊 ACTUAL MOLECULAR DISTRIBUTION SUMMARY")
    print("="*50)
    
    total_sites = len(coordinates)
    occupied_sites_count = len(occupied_sites)
    
    print(f"🔹 Total lattice sites: {total_sites}")
    print(f"🔹 Occupied sites: {occupied_sites_count}")
    print(f"🔹 Empty sites: {total_sites - occupied_sites_count}")
    print(f"🔹 Surface coverage: {occupied_sites_count/total_sites*100:.1f}%")
    print()
    
    print("🧪 SPECIES COUNTS:")
    for species_id, count in species_counts.items():
        name = species_names[species_id]
        print(f"   {name}: {count} molecules")
    
    print()
    print("🏗️ SITE TYPE PREFERENCES:")
    for site_type, count in site_type_counts.items():
        name = site_names[site_type]
        total_of_type = sum(1 for site in coordinates.values() if site['site_type'] == site_type)
        occupancy_rate = count / total_of_type * 100 if total_of_type > 0 else 0
        print(f"   {name}: {count}/{total_of_type} sites occupied ({occupancy_rate:.1f}%)")
    
    print()
    print("✅ ANSWER TO YOUR QUESTION:")
    print("❌ Previous visualizations: Used random placement (wrong positions)")
    print("✅ This visualization: Uses ACTUAL molecular coordinates from simulation")
    print("   - Shows real spatial patterns")
    print("   - Shows site-specific preferences") 
    print("   - Shows actual clustering effects")
    
    return True

if __name__ == "__main__":
    print("🎯 CREATING REAL MOLECULAR POSITION VISUALIZATION")
    print("=" * 60)
    
    success = create_real_molecular_visualization()
    
    if success:
        print("\\n🎉 SUCCESS: Created visualization with ACTUAL molecular positions!")
    else:
        print("\\n❌ Failed to create visualization") 