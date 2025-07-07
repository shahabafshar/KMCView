#!/usr/bin/env python3
"""
Extract REAL molecular positions from Zacros simulation data
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def extract_real_positions_from_restart():
    """Extract actual molecular positions from restart.inf"""
    print("🔍 Extracting ACTUAL molecular positions from restart.inf...")
    
    real_positions = {}
    
    with open('input-output/restart.inf', 'r') as f:
        lines = f.readlines()
    
    # Parse the restart file - look for occupancy data
    parsing_sites = False
    
    for i, line in enumerate(lines):
        line = line.strip()
        if not line:
            continue
            
        # Look for the section with site occupancy data
        if 'site' in line.lower() and ('occupancy' in line.lower() or 'entity' in line.lower()):
            parsing_sites = True
            continue
            
        if parsing_sites and line:
            parts = line.split()
            if len(parts) >= 3:
                try:
                    # Format might be: site_id species_id count
                    site_id = int(parts[0])
                    species_id = int(parts[1])
                    count = int(parts[2]) if len(parts) > 2 else 1
                    
                    if species_id > 0:  # Occupied site
                        real_positions[site_id] = {
                            'species_id': species_id,
                            'count': count
                        }
                except (ValueError, IndexError):
                    continue
    
    print(f"Found {len(real_positions)} occupied sites from restart.inf")
    return real_positions

def extract_lattice_coordinates():
    """Extract site coordinates from lattice_output.txt"""
    print("📍 Loading lattice coordinates...")
    
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
    
    print(f"Loaded coordinates for {len(coordinates)} lattice sites")
    return coordinates

def create_real_position_visualization():
    """Create visualization using real molecular positions"""
    
    # Get occupancy data
    occupancy = extract_real_positions_from_restart()
    
    # Get coordinates
    coordinates = extract_lattice_coordinates()
    
    if not occupancy:
        print("❌ No occupancy data found - creating alternative visualization")
        return create_alternative_visualization()
    
    # Combine occupancy with coordinates
    real_molecules = []
    
    for site_id, occ_data in occupancy.items():
        if site_id in coordinates:
            coord_data = coordinates[site_id]
            real_molecules.append({
                'x': coord_data['x'],
                'y': coord_data['y'],
                'site_type': coord_data['site_type'],
                'species_id': occ_data['species_id'],
                'count': occ_data['count']
            })
    
    if not real_molecules:
        print("❌ Could not match occupancy with coordinates")
        return create_alternative_visualization()
    
    # Create the visualization
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Plot 1: All sites with species
    species_colors = {1: 'red', 2: 'blue', 3: 'green', 0: 'lightgray'}
    species_names = {1: 'H*', 2: 'GeH2*', 3: 'GeH3*'}
    
    for species_id in [0, 1, 2, 3]:
        molecules = [m for m in real_molecules if m['species_id'] == species_id]
        if molecules:
            x_coords = [m['x'] for m in molecules]
            y_coords = [m['y'] for m in molecules]
            
            label = species_names.get(species_id, 'Empty')
            if species_id > 0:
                label += f' ({len(molecules)})'
            
            ax1.scatter(x_coords, y_coords, 
                       c=species_colors[species_id],
                       label=label, s=20, alpha=0.7)
    
    ax1.set_xlabel('X Coordinate (Å)')
    ax1.set_ylabel('Y Coordinate (Å)')
    ax1.set_title('REAL Molecular Positions\\n(Actual simulation data)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_aspect('equal')
    
    # Plot 2: Site types
    site_colors = {1: 'lightblue', 2: 'lightgreen', 3: 'lightcoral', 4: 'lightyellow'}
    site_names = {1: 'top1', 2: 'bridge1', 3: 'top2', 4: 'bridge2'}
    
    for site_type in [1, 2, 3, 4]:
        molecules = [m for m in real_molecules if m['site_type'] == site_type]
        if molecules:
            x_coords = [m['x'] for m in molecules]
            y_coords = [m['y'] for m in molecules]
            
            ax2.scatter(x_coords, y_coords,
                       c=site_colors[site_type],
                       label=f'{site_names[site_type]} ({len(molecules)})',
                       s=20, alpha=0.7)
    
    ax2.set_xlabel('X Coordinate (Å)')
    ax2.set_ylabel('Y Coordinate (Å)')
    ax2.set_title('Site Types Distribution\\n(Surface structure)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig('real_molecular_positions.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("✅ REAL molecular positions saved as real_molecular_positions.png")
    
    # Print summary
    species_counts = {}
    for mol in real_molecules:
        species_id = mol['species_id']
        species_counts[species_id] = species_counts.get(species_id, 0) + mol['count']
    
    print("\n📊 ACTUAL molecular distribution:")
    for species_id, count in species_counts.items():
        name = species_names.get(species_id, f'Species_{species_id}')
        print(f"   {name}: {count} molecules")
    
    return True

def create_alternative_visualization():
    """Create visualization showing the comparison between random vs should-be-real"""
    print("Creating comparison: Random filling vs Should-be-real positions")
    
    # This is what I did before (random) vs what we should show
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Left: What I did (random filling)
    np.random.seed(42)
    lattice_size = 20
    grid = np.zeros((lattice_size, lattice_size))
    
    # Fill with random positions (representative)
    positions = [1]*200 + [2]*200 + [3]*100  # Simplified counts
    np.random.shuffle(positions)
    
    for i, pos in enumerate(positions):
        if i < lattice_size * lattice_size:
            x, y = i // lattice_size, i % lattice_size
            grid[x, y] = pos
    
    im1 = ax1.imshow(grid, cmap=plt.cm.Set3, origin='lower')
    ax1.set_title('What I Did:\\nRANDOM Filling\\n(Correct counts, wrong positions)')
    ax1.set_xticks([])
    ax1.set_yticks([])
    
    # Right: What should be shown
    ax2.text(0.5, 0.5, 'What SHOULD be shown:\\n\\n' +
             '✅ Real molecular coordinates\\n' +
             '✅ Actual spatial clustering\\n' +
             '✅ Site-specific preferences\\n' +
             '✅ True surface patterns\\n\\n' +
             'FROM:\\n' +
             '• lattice_output.txt (coordinates)\\n' +
             '• restart.inf (final occupancy)\\n\\n' +
             'This would show the ACTUAL\\n' +
             'molecular arrangement that\\n' +
             'resulted from the simulation.',
             transform=ax2.transAxes, fontsize=12,
             ha='center', va='center',
             bbox=dict(boxstyle='round,pad=0.5', facecolor='lightblue'))
    
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.set_title('What SHOULD be shown:\\nREAL Positions\\n(Actual simulation results)')
    
    plt.tight_layout()
    plt.savefig('comparison_random_vs_real.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("✅ Comparison saved as comparison_random_vs_real.png")
    return False

if __name__ == "__main__":
    print("🎯 EXTRACTING REAL MOLECULAR POSITIONS")
    print("=" * 50)
    
    success = create_real_position_visualization()
    
    if success:
        print("\\n✅ SUCCESS: Used actual simulation data!")
    else:
        print("\\n⚠️  Used alternative approach - may need manual data extraction")
    
    print("\\nTo answer your question:")
    print("❌ Previous visualizations: Random filling (correct counts, wrong positions)")
    print("✅ This attempt: Extracting real positions from simulation files") 