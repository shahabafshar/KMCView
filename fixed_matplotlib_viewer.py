#!/usr/bin/env python3
"""
🔬 Enhanced Matplotlib Molecular Evolution Viewer
===============================================

Features:
• Dynamic lattice structure parsing from lattice_input.dat
• Works with any set of Zacros input files
• Configurable file paths and parameters
• Proper 2D lattice mesh visualization
• Working play/pause animation controls
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, TextBox
from matplotlib.animation import FuncAnimation
import os
import re
import warnings
warnings.filterwarnings('ignore')

class LatticeParser:
    """Parse lattice_input.dat file to extract structure parameters"""
    
    def __init__(self, lattice_file='input-output/lattice_input.dat'):
        self.lattice_file = lattice_file
        self.cell_vectors = None
        self.repeat_cell = None
        self.site_type_names = []
        self.site_coordinates = []
        self.site_types = []
        self.parse_lattice_file()
    
    def parse_lattice_file(self):
        """Parse the lattice input file"""
        if not os.path.exists(self.lattice_file):
            print(f"❌ Lattice file not found: {self.lattice_file}")
            print("   Using default parameters...")
            self.set_default_parameters()
            return
        
        print(f"📖 Reading lattice structure from: {self.lattice_file}")
        
        try:
            with open(self.lattice_file, 'r') as f:
                lines = f.readlines()
            
            i = 0
            while i < len(lines):
                line = lines[i].strip()
                
                if line.startswith('cell_vectors'):
                    # Read cell vectors (2x2 matrix)
                    self.cell_vectors = []
                    i += 1
                    for _ in range(2):
                        if i < len(lines):
                            vector_line = lines[i].strip()
                            if vector_line and not vector_line.startswith('#'):
                                try:
                                    values = [float(x) for x in vector_line.split()]
                                    if len(values) >= 2:
                                        self.cell_vectors.append(values[:2])
                                except ValueError:
                                    pass
                            i += 1
                
                elif line.startswith('repeat_cell'):
                    # Read repeat cell dimensions
                    parts = line.split()
                    if len(parts) >= 3:
                        try:
                            self.repeat_cell = [int(parts[1]), int(parts[2])]
                        except ValueError:
                            pass
                
                elif line.startswith('site_type_names'):
                    # Read site type names
                    i += 1
                    while i < len(lines):
                        type_line = lines[i].strip()
                        if type_line and not type_line.startswith('#') and not type_line.startswith('n_'):
                            self.site_type_names.extend(type_line.split())
                            break
                        i += 1
                
                elif line.startswith('site_types'):
                    # Read site types for each position
                    i += 1
                    while i < len(lines):
                        type_line = lines[i].strip()
                        if type_line and not type_line.startswith('#') and not type_line.startswith('site_coordinates'):
                            self.site_types.extend(type_line.split())
                            break
                        i += 1
                
                elif line.startswith('site_coordinates'):
                    # Read fractional coordinates
                    i += 1
                    while i < len(lines):
                        coord_line = lines[i].strip()
                        if coord_line and not coord_line.startswith('#'):
                            if coord_line.startswith('neighboring_structure'):
                                break
                            try:
                                values = [float(x) for x in coord_line.split()]
                                if len(values) >= 2:
                                    self.site_coordinates.append((values[0], values[1]))
                            except ValueError:
                                pass
                        i += 1
                
                i += 1
                
        except Exception as e:
            print(f"❌ Error parsing lattice file: {e}")
            print("   Using default parameters...")
            self.set_default_parameters()
    
    def set_default_parameters(self):
        """Set default parameters if file parsing fails"""
        self.cell_vectors = [[6.133780000000000, 0.0], [0.0, 6.133780000000000]]
        self.repeat_cell = [20, 20]
        self.site_type_names = ['top1', 'bridge1', 'top2', 'bridge2']
        self.site_types = ['top1', 'bridge1', 'top2', 'bridge2']
        self.site_coordinates = [(0.0, 0.0), (0.0, 0.25), (0.0, 0.5), (0.0, 0.75)]
    
    def get_lattice_info(self):
        """Return parsed lattice information"""
        return {
            'cell_vectors': self.cell_vectors,
            'repeat_cell': self.repeat_cell,
            'site_type_names': self.site_type_names,
            'site_types': self.site_types,
            'site_coordinates': self.site_coordinates
        }

class FixedMolecularViewer:
    def __init__(self, data_dir='input-output'):
        """Initialize the viewer with configurable data directory"""
        self.data_dir = data_dir
        self.evolution_df = pd.DataFrame()
        self.coordinates = {}
        self.lattice_bounds = {}
        self.lattice_info = {}
        
        # Animation state
        self.current_step = 0
        self.is_playing = False
        self.animation_obj = None
        
        # Initialize
        self.load_data()
        self.load_positions()
        self.setup_figure()
        
        plt.show()
    
    def load_data(self):
        """Load species evolution data from specnum_output.txt"""
        data_file = os.path.join(self.data_dir, 'specnum_output.txt')
        
        if not os.path.exists(data_file):
            print(f"❌ Data file not found: {data_file}")
            self.evolution_df = pd.DataFrame({'time': [0.0], 'H*': [0.0], 'GeH2*': [0.0], 'GeH3*': [0.0]})
            return
        
        try:
            with open(data_file, 'r') as f:
                lines = f.readlines()
            
            data = []
            for line in lines:
                if line.strip() and not line.startswith('#'):
                    try:
                        parts = line.strip().split()
                        if len(parts) >= 5:
                            float(parts[0])  # Test if numeric
                            step = int(parts[0])
                            nevents = int(parts[1])
                            time = float(parts[2])
                            h_star = float(parts[5]) if len(parts) > 5 else 0.0
                            geh2_star = float(parts[6]) if len(parts) > 6 else 0.0
                            geh3_star = float(parts[7]) if len(parts) > 7 else 0.0
                            
                            data.append({
                                'step': step,
                                'nevents': nevents,
                                'time': time,
                                'H*': h_star,
                                'GeH2*': geh2_star,
                                'GeH3*': geh3_star
                            })
                    except (ValueError, IndexError):
                        continue
            
            self.evolution_df = pd.DataFrame(data)
            print(f"✅ Loaded {len(self.evolution_df)} time points from {data_file}")
            
        except Exception as e:
            print(f"❌ Error loading data: {e}")
            self.evolution_df = pd.DataFrame({'time': [0.0], 'H*': [0.0], 'GeH2*': [0.0], 'GeH3*': [0.0]})
    
    def load_positions(self):
        """Generate lattice positions from lattice_input.dat structure"""
        lattice_file = os.path.join(self.data_dir, 'lattice_input.dat')
        
        # Parse lattice structure
        parser = LatticeParser(lattice_file)
        self.lattice_info = parser.get_lattice_info()
        
        # Extract parameters
        cell_vectors = self.lattice_info['cell_vectors']
        repeat_cell = self.lattice_info['repeat_cell']
        site_type_names = self.lattice_info['site_type_names']
        site_coordinates = self.lattice_info['site_coordinates']
        
        # Calculate cell dimensions
        cell_vector_x = cell_vectors[0][0]  # x-component of first vector
        cell_vector_y = cell_vectors[1][1]  # y-component of second vector
        repeat_x, repeat_y = repeat_cell
        
        print(f"   Cell vectors: {cell_vector_x:.6f} × {cell_vector_y:.6f} Å")
        print(f"   Repeat pattern: {repeat_x} × {repeat_y}")
        print(f"   Site types: {site_type_names}")
        
        # Generate all lattice sites
        site_id = 1
        for i in range(repeat_x):
            for j in range(repeat_y):
                # Base coordinates for this unit cell
                base_x = i * cell_vector_x
                base_y = j * cell_vector_y
                
                # Add each site in the unit cell
                for k, (frac_x, frac_y) in enumerate(site_coordinates):
                    # Convert fractional to absolute coordinates
                    abs_x = base_x + frac_x * cell_vector_x
                    abs_y = base_y + frac_y * cell_vector_y
                    
                    site_type = site_type_names[k] if k < len(site_type_names) else f'site_{k}'
                    
                    self.coordinates[site_id] = {
                        'x': abs_x,
                        'y': abs_y,
                        'site_type': site_type,
                        'unit_cell': (i, j),
                        'local_site': k
                    }
                    site_id += 1
        
        # Calculate lattice boundaries
        self.lattice_bounds = {
            'x_min': 0.0,
            'x_max': repeat_x * cell_vector_x,
            'y_min': 0.0,
            'y_max': repeat_y * cell_vector_y,
            'cell_size_x': cell_vector_x,
            'cell_size_y': cell_vector_y
        }
        
        print(f"✅ Generated {len(self.coordinates)} lattice sites")
        print(f"   Lattice dimensions: {self.lattice_bounds['x_max']:.2f} × {self.lattice_bounds['y_max']:.2f} Å")
        print(f"   Sites per unit cell: {len(site_coordinates)}")

    def setup_figure(self):
        """Create the figure and controls"""
        self.fig, self.ax = plt.subplots(figsize=(14, 10))
        plt.subplots_adjust(bottom=0.25, top=0.9)
        
        # Info text
        self.info_text = self.fig.text(0.5, 0.95, '', ha='center', fontsize=12, 
                                      bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.8))
        
        # Initial plot
        self.update_plot(0)
        
        # Set up controls
        self.setup_controls()
        
        # Set up animation (but don't start it)
        self.setup_animation()
    
    def setup_controls(self):
        """Set up interactive controls"""
        # Slider
        ax_slider = plt.axes([0.2, 0.1, 0.5, 0.03])
        self.slider = Slider(ax_slider, 'Time Step', 0, len(self.evolution_df)-1, 
                           valinit=0, valfmt='%d', valstep=1)
        self.slider.on_changed(self.on_slider_change)
        
        # Play button
        ax_play = plt.axes([0.75, 0.1, 0.08, 0.04])
        self.play_button = Button(ax_play, 'Play')
        self.play_button.on_clicked(self.toggle_play)
        
        # Reset button
        ax_reset = plt.axes([0.85, 0.1, 0.08, 0.04])
        self.reset_button = Button(ax_reset, 'Reset')
        self.reset_button.on_clicked(self.reset_animation)
        
        # Text input
        ax_text = plt.axes([0.2, 0.05, 0.1, 0.03])
        self.text_box = TextBox(ax_text, 'Step:', initial='0')
        self.text_box.on_submit(self.on_text_submit)
    
    def setup_animation(self):
        """Set up the FuncAnimation object"""
        def animate_frame(frame):
            if self.is_playing:
                self.current_step = (self.current_step + 1) % len(self.evolution_df)
                self.slider.set_val(self.current_step)
                self.update_plot(self.current_step)
                return []
            return []
        
        # Create animation object but don't start it
        self.animation_obj = FuncAnimation(self.fig, animate_frame, 
                                         interval=300, blit=False, repeat=True)
        self.animation_obj.event_source.stop()  # Stop initially
    
    def update_plot(self, step):
        """Update the plot for a given time step"""
        self.ax.clear()
        
        if step >= len(self.evolution_df):
            step = len(self.evolution_df) - 1
        
        # Get current data
        current_data = self.evolution_df.iloc[step]
        current_time = current_data['time']
        h_count = int(current_data['H*'])
        geh2_count = int(current_data['GeH2*'])
        geh3_count = int(current_data['GeH3*'])
        
        # Generate positions
        positions = self.generate_positions_for_step(step)
        
        # Plot molecules
        colors = {'H*': 'red', 'GeH2*': 'blue', 'GeH3*': 'green'}
        sizes = {'H*': 30, 'GeH2*': 60, 'GeH3*': 90}
        
        for species, coords in positions.items():
            if coords:
                x_coords = [c['x'] for c in coords]
                y_coords = [c['y'] for c in coords]
                self.ax.scatter(x_coords, y_coords, 
                              c=colors[species], s=sizes[species], 
                              alpha=0.7, label=f'{species} ({len(coords)})')
        
        # Set axis limits using lattice bounds
        if hasattr(self, 'lattice_bounds'):
            x_margin = (self.lattice_bounds['x_max'] - self.lattice_bounds['x_min']) * 0.05
            y_margin = (self.lattice_bounds['y_max'] - self.lattice_bounds['y_min']) * 0.05
            
            self.ax.set_xlim(self.lattice_bounds['x_min'] - x_margin, self.lattice_bounds['x_max'] + x_margin)
            self.ax.set_ylim(self.lattice_bounds['y_min'] - y_margin, self.lattice_bounds['y_max'] + y_margin)
        
        # Add lattice structure visualization (light grid)
        if hasattr(self, 'lattice_bounds') and len(self.coordinates) > 0:
            # Draw unit cell boundaries
            cell_size_x = self.lattice_bounds['cell_size_x']
            cell_size_y = self.lattice_bounds['cell_size_y']
            repeat_x = int(self.lattice_bounds['x_max'] / cell_size_x)
            repeat_y = int(self.lattice_bounds['y_max'] / cell_size_y)
            
            for i in range(repeat_x + 1):
                x_line = i * cell_size_x
                if x_line <= self.lattice_bounds['x_max']:
                    self.ax.axvline(x=x_line, color='lightgray', alpha=0.3, linewidth=0.5)
            
            for j in range(repeat_y + 1):
                y_line = j * cell_size_y
                if y_line <= self.lattice_bounds['y_max']:
                    self.ax.axhline(y=y_line, color='lightgray', alpha=0.3, linewidth=0.5)
                    
            # Show empty lattice sites as small dots
            empty_x = []
            empty_y = []
            for coord in self.coordinates.values():
                empty_x.append(coord['x'])
                empty_y.append(coord['y'])
            if empty_x:
                self.ax.scatter(empty_x, empty_y, c='lightgray', s=1, alpha=0.4, zorder=1)
        
        # Labels and formatting
        self.ax.set_xlabel('X Coordinate (Å)', fontsize=12)
        self.ax.set_ylabel('Y Coordinate (Å)', fontsize=12)
        self.ax.set_aspect('equal')
        self.ax.grid(True, alpha=0.3)
        self.ax.legend(loc='upper left')
        
        # Update info text
        info_str = (f"Time: {current_time:.1f} | Step: {step+1}/{len(self.evolution_df)} | "
                   f"H*: {h_count} | GeH2*: {geh2_count} | GeH3*: {geh3_count}")
        self.info_text.set_text(info_str)
        
        # Only draw if not animating (to avoid conflicts)
        if not self.is_playing:
            self.fig.canvas.draw_idle()
    
    def generate_positions_for_step(self, step):
        """Generate molecular positions for a given time step"""
        if step >= len(self.evolution_df):
            return {}
        
        current_data = self.evolution_df.iloc[step]
        h_count = int(current_data['H*'])
        geh2_count = int(current_data['GeH2*'])
        geh3_count = int(current_data['GeH3*'])
        
        total_molecules = h_count + geh2_count + geh3_count
        
        if total_molecules == 0:
            return {}
        
        # Use available site coordinates
        all_sites = list(self.coordinates.keys())
        if not all_sites:
            return {}
        
        # Generate consistent positions
        np.random.seed(42 + step)
        occupied_sites = np.random.choice(all_sites, 
                                        size=min(total_molecules, len(all_sites)), 
                                        replace=False)
        
        positions = {}
        idx = 0
        
        # Assign positions to species
        for species, count in [('H*', h_count), ('GeH2*', geh2_count), ('GeH3*', geh3_count)]:
            positions[species] = []
            for _ in range(count):
                if idx < len(occupied_sites):
                    site_id = occupied_sites[idx]
                    positions[species].append(self.coordinates[site_id])
                    idx += 1
        
        return positions
    
    def on_slider_change(self, val):
        """Handle slider value changes"""
        step = int(val)
        self.current_step = step
        self.update_plot(step)
        self.text_box.set_val(str(step))
    
    def toggle_play(self, event):
        """Toggle play/pause animation"""
        if self.is_playing:
            print("⏸️ Animation paused")
            self.is_playing = False
            self.animation_obj.event_source.stop()
            self.play_button.label.set_text('Play')
        else:
            print("▶️ Animation playing")
            self.is_playing = True
            self.animation_obj.event_source.start()
            self.play_button.label.set_text('Pause')
        
        self.fig.canvas.draw_idle()
    
    def reset_animation(self, event):
        """Reset animation to beginning"""
        print("🔄 Animation reset")
        self.is_playing = False
        self.animation_obj.event_source.stop()
        self.current_step = 0
        self.slider.set_val(0)
        self.update_plot(0)
        self.text_box.set_val('0')
        self.play_button.label.set_text('Play')
        self.fig.canvas.draw_idle()
    
    def on_text_submit(self, text):
        """Handle text input submission"""
        try:
            step = int(text)
            step = max(0, min(step, len(self.evolution_df) - 1))
            self.current_step = step
            self.slider.set_val(step)
            self.update_plot(step)
        except ValueError:
            print(f"❌ Invalid step number: {text}")

if __name__ == '__main__':
    # Usage: python fixed_matplotlib_viewer.py [data_directory]
    import sys
    
    if len(sys.argv) > 1:
        data_dir = sys.argv[1]
    else:
        data_dir = 'input-output'
    
    print(f"🔬 Starting Enhanced Molecular Viewer with data from: {data_dir}")
    viewer = FixedMolecularViewer(data_dir=data_dir) 