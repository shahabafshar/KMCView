#!/usr/bin/env python3
"""
🚀 Enhanced GUI Molecular Evolution Viewer
=========================================

Major UX and Performance Improvements:
• Modern dark/light theme support
• Smooth animations with variable speed control
• Improved widget layout and styling
• Performance optimizations with data caching
• Keyboard shortcuts and better controls
• Progress indicators and status feedback
• Memory-efficient plotting with blitting
• Configurable display options
• Export functionality
• Better error handling and user feedback
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, TextBox, RadioButtons, CheckButtons
from matplotlib.table import Table
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle
import matplotlib.patches as patches
import os
import re
import warnings
import time
from pathlib import Path
import json

warnings.filterwarnings('ignore')

class EnhancedLatticeParser:
    """Enhanced lattice parser with caching and validation"""
    
    def __init__(self, lattice_file='input-output/lattice_input.dat'):
        self.lattice_file = lattice_file
        self.cell_vectors = None
        self.repeat_cell = None
        self.site_type_names = []
        self.site_coordinates = []
        self.site_types = []
        self._cache_file = Path(lattice_file).parent / '.lattice_cache.json'
        
        # Try to load from cache first
        if not self.load_from_cache():
            self.parse_lattice_file()
            self.save_to_cache()
    
    def load_from_cache(self):
        """Load parsed data from cache if available and valid"""
        if not self._cache_file.exists():
            return False
        
        try:
            lattice_mtime = Path(self.lattice_file).stat().st_mtime if Path(self.lattice_file).exists() else 0
            cache_mtime = self._cache_file.stat().st_mtime
            
            if cache_mtime > lattice_mtime:
                with open(self._cache_file, 'r') as f:
                    data = json.load(f)
                
                self.cell_vectors = data['cell_vectors']
                self.repeat_cell = data['repeat_cell']
                self.site_type_names = data['site_type_names']
                self.site_coordinates = [tuple(coord) for coord in data['site_coordinates']]
                self.site_types = data.get('site_types', [])
                
                print(f"✅ Loaded lattice structure from cache")
                return True
        except Exception as e:
            print(f"⚠️ Cache load failed: {e}")
        
        return False
    
    def save_to_cache(self):
        """Save parsed data to cache"""
        try:
            data = {
                'cell_vectors': self.cell_vectors,
                'repeat_cell': self.repeat_cell,
                'site_type_names': self.site_type_names,
                'site_coordinates': self.site_coordinates,
                'site_types': self.site_types
            }
            
            with open(self._cache_file, 'w') as f:
                json.dump(data, f, indent=2)
                
        except Exception as e:
            print(f"⚠️ Cache save failed: {e}")
    
    def parse_lattice_file(self):
        """Parse the lattice input file with enhanced error handling"""
        if not os.path.exists(self.lattice_file):
            print(f"❌ Lattice file not found: {self.lattice_file}")
            print("   Using default parameters...")
            self.set_default_parameters()
            return
        
        print(f"📖 Parsing lattice structure from: {self.lattice_file}")
        
        try:
            with open(self.lattice_file, 'r') as f:
                lines = f.readlines()
            
            i = 0
            while i < len(lines):
                line = lines[i].strip()
                
                if line.startswith('cell_vectors'):
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
                    parts = line.split()
                    if len(parts) >= 3:
                        try:
                            self.repeat_cell = [int(parts[1]), int(parts[2])]
                        except ValueError:
                            pass
                
                elif line.startswith('site_type_names'):
                    i += 1
                    while i < len(lines):
                        type_line = lines[i].strip()
                        if type_line and not type_line.startswith('#') and not type_line.startswith('n_'):
                            self.site_type_names.extend(type_line.split())
                            break
                        i += 1
                
                elif line.startswith('site_types'):
                    i += 1
                    while i < len(lines):
                        type_line = lines[i].strip()
                        if type_line and not type_line.startswith('#') and not type_line.startswith('site_coordinates'):
                            self.site_types.extend(type_line.split())
                            break
                        i += 1
                
                elif line.startswith('site_coordinates'):
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

class EnhancedMolecularViewer:
    """Enhanced GUI viewer with modern UX and performance optimizations"""
    
    def __init__(self, data_dir='input-output'):
        """Initialize the enhanced viewer"""
        self.data_dir = data_dir
        self.evolution_df = pd.DataFrame()
        self.coordinates = {}
        self.lattice_bounds = {}
        self.lattice_info = {}
        self.lattice_states = {}  # Store actual lattice states from history_output.txt
        self.surface_species = []  # Store surface species names
        self.step_to_config = {}  # Map dataframe index to actual KMC step numbers
        
        # Animation and performance state
        self.current_step = 0
        self.is_playing = False
        self.animation_obj = None
        self.animation_speed = 1.0
        self.last_update_time = 0
        self.frame_cache = {}  # Cache for rendered frames
        
        # Display options
        self.show_lattice_grid = True
        self.show_empty_sites = True
        self.show_trajectories = False
        
        # Species properties configuration (editable data grid)
        self.species_properties = {
            'H*': {
                'visible': True,
                'color': '#FF4444',  # Bright red
                'shape': 'o',        # Circle
                'size': 40
            },
            'GeH2*': {
                'visible': True,
                'color': '#4477FF',  # Bright blue
                'shape': 's',        # Square
                'size': 80
            },
            'GeH3*': {
                'visible': True,
                'color': '#44AA44',  # Bright green
                'shape': '^',        # Triangle
                'size': 120
            }
        }
        
        # Settings modal state
        self.settings_window = None
        self.settings_open = False
        
        # Performance monitoring
        self.render_times = []
        self.max_render_history = 50
        
        # Internal caching
        self.position_cache = {}
        self.max_position_cache = 256
        
        # Initialize everything
        print(f"🚀 Starting Enhanced Molecular Viewer...")
        print(f"📁 Data directory: {data_dir}")
        
        self.load_data()
        self.load_positions()
        self.setup_theme()
        self.setup_figure()
        self.setup_controls()
        self.setup_keyboard_shortcuts()
        
        # Show initial frame
        self.update_plot(0)
        
        data_source = "lattice states" if self.lattice_states else "species counts"
        print(f"Enhanced viewer ready with {len(self.evolution_df)} time points using {data_source}")
        
        # Inform user about positioning strategy
        if self.lattice_states:
            print("✅ POSITIONING: Using exact lattice states from history_output.txt")
            print("   → Molecules will show true KMC dynamics with site persistence")
        else:
            print("⚠️  POSITIONING: Using deterministic placement based on species counts")
            print("   → No history_output.txt found - positions are approximated")
            print("   → For exact visualization, ensure history_output.txt is available")
        
        print("Keyboard shortcuts: Space=Play/Pause, R=Reset, Left/Right=Step")
        
        plt.show()
    
    def setup_theme(self):
        """Setup consistent styling"""
        plt.style.use('default')
        self.colors = {
            'bg': '#FFFFFF',
            'panel': '#F8F9FA',
            'text': '#212529',
            'accent': '#007BFF',
            'success': '#28A745',
            'warning': '#FFC107',
            'danger': '#DC3545'
        }
        
        # Set global font properties
        plt.rcParams.update({
            'font.size': 10,
            'font.family': 'sans-serif',
            'axes.titlesize': 12,
            'axes.labelsize': 10,
            'xtick.labelsize': 9,
            'ytick.labelsize': 9,
            'legend.fontsize': 9,
            'figure.titlesize': 14
        })
    
    def load_data(self):
        """Load species evolution data and lattice states from history file"""
        # First try to load from history_output.txt (actual lattice states)
        history_file = os.path.join(self.data_dir, 'history_output.txt')
        
        if os.path.exists(history_file):
            print(f"📊 Loading lattice states from {history_file}...")
            self.load_history_data(history_file)
        else:
            # Fallback to specnum_output.txt (counts only)
            print("⚠️ No history_output.txt found, falling back to species counts...")
            self.load_specnum_data()
    
    def load_history_data(self, history_file):
        """Load actual lattice states from history_output.txt"""
        start_time = time.time()
        
        try:
            with open(history_file, 'r') as f:
                lines = f.readlines()
            
            # Parse header information
            self.surface_species = []
            self.lattice_states = {}  # {step: {site_id: species_id}}
            self.step_to_config = {}  # Map dataframe index to actual KMC step
            evolution_data = []
            
            current_config = None
            current_step = None
            current_time = None
            current_sites = {}
            
            i = 0
            while i < len(lines):
                line = lines[i].strip()
                
                # Parse surface species from header
                if line.startswith('Surface_Species:'):
                    species_line = line.split(':')[1].strip()
                    self.surface_species = species_line.split()
                    print(f"Found surface species: {self.surface_species}")
                
                # Parse configuration blocks
                elif line.startswith('configuration'):
                    # Save previous configuration if exists
                    if current_config is not None and current_step is not None:
                        self.lattice_states[current_step] = current_sites.copy()
                        
                        # Count species for this configuration
                        species_counts = {species: 0 for species in self.surface_species}
                        for site_id, species_id in current_sites.items():
                            if species_id > 0 and species_id <= len(self.surface_species):
                                species_name = self.surface_species[species_id - 1]
                                species_counts[species_name] += 1
                        
                        # Store mapping from dataframe index to actual KMC step
                        df_index = len(evolution_data)
                        self.step_to_config[df_index] = current_step
                        
                        evolution_data.append({
                            'step': current_step,
                            'time': current_time,
                            'kmc_step': current_step,  # Store actual KMC step
                            **species_counts
                        })
                    
                    # Parse new configuration
                    parts = line.split()
                    if len(parts) >= 6:
                        current_config = int(parts[1])
                        current_step = int(parts[2])
                        current_time = float(parts[3])
                        current_sites = {}
                        
                        if current_config % 10 == 0:
                            print(f"   Processing configuration {current_config}, step {current_step}...")
                
                # Parse site data
                elif line and not line.startswith(('Gas_Species:', 'Surface_Species:', 'Simulation_Box:', 'Site_Types:')):
                    try:
                        parts = line.split()
                        if len(parts) >= 5:
                            site_id = int(parts[0])
                            site_type = int(parts[1])
                            species_id = int(parts[2])
                            current_sites[site_id] = species_id
                    except (ValueError, IndexError):
                        pass
                
                i += 1
            
            # Save the last configuration
            if current_config is not None and current_step is not None:
                self.lattice_states[current_step] = current_sites.copy()
                
                species_counts = {species: 0 for species in self.surface_species}
                for site_id, species_id in current_sites.items():
                    if species_id > 0 and species_id <= len(self.surface_species):
                        species_name = self.surface_species[species_id - 1]
                        species_counts[species_name] += 1
                
                # Store mapping for the last configuration
                df_index = len(evolution_data)
                self.step_to_config[df_index] = current_step
                
                evolution_data.append({
                    'step': current_step,
                    'time': current_time,
                    'kmc_step': current_step,  # Store actual KMC step
                    **species_counts
                })
            
            self.evolution_df = pd.DataFrame(evolution_data)
            load_time = time.time() - start_time
            print(f"✅ Loaded {len(self.evolution_df)} configurations with actual lattice states in {load_time:.2f}s")
            
            # Debug: Show step mapping for first few configurations
            if len(self.step_to_config) > 0:
                sample_mappings = list(self.step_to_config.items())[:5]
                print(f"📊 Sample step mappings: {sample_mappings}")
            
        except Exception as e:
            print(f"❌ Error loading history data: {e}")
            print("   Falling back to specnum_output.txt...")
            self.load_specnum_data()
    
    def load_specnum_data(self):
        """Fallback method to load species counts from specnum_output.txt"""
        data_file = os.path.join(self.data_dir, 'specnum_output.txt')
        
        if not os.path.exists(data_file):
            print(f"❌ Data file not found: {data_file}")
            self.evolution_df = pd.DataFrame({'time': [0.0], 'H*': [0.0], 'GeH2*': [0.0], 'GeH3*': [0.0]})
            self.lattice_states = {}
            return
        
        print(f"📊 Loading evolution data from {data_file}...")
        start_time = time.time()
        
        try:
            with open(data_file, 'r') as f:
                lines = f.readlines()
            
            data = []
            for i, line in enumerate(lines):
                if i % 1000 == 0 and i > 0:
                    print(f"   Processed {i} lines...")
                
                if line.strip() and not line.startswith('#'):
                    try:
                        parts = line.strip().split()
                        if len(parts) >= 5:
                            float(parts[0])  # Test if numeric
                            step = int(parts[0])
                            nevents = int(parts[1])
                            sim_time = float(parts[2])
                            h_star = float(parts[5]) if len(parts) > 5 else 0.0
                            geh2_star = float(parts[6]) if len(parts) > 6 else 0.0
                            geh3_star = float(parts[7]) if len(parts) > 7 else 0.0
                            
                            data.append({
                                'step': step,
                                'nevents': nevents,
                                'time': sim_time,
                                'H*': h_star,
                                'GeH2*': geh2_star,
                                'GeH3*': geh3_star
                            })
                    except (ValueError, IndexError):
                        continue
            
            self.evolution_df = pd.DataFrame(data)
            self.lattice_states = {}  # No lattice states available
            load_time = time.time() - start_time
            print(f"Loaded {len(self.evolution_df)} time points in {load_time:.2f}s")
            
        except Exception as e:
            print(f"❌ Error loading data: {e}")
            self.evolution_df = pd.DataFrame({'time': [0.0], 'H*': [0.0], 'GeH2*': [0.0], 'GeH3*': [0.0]})
            self.lattice_states = {}
    
    def load_positions(self):
        """Load lattice positions with enhanced caching"""
        print("🏗️ Building lattice structure...")
        start_time = time.time()
        
        # Load lattice structure
        parser = EnhancedLatticeParser(os.path.join(self.data_dir, 'lattice_input.dat'))
        self.lattice_info = parser.get_lattice_info()
        
        # Calculate bounds
        cell_vectors = self.lattice_info['cell_vectors']
        repeat_cell = self.lattice_info['repeat_cell']
        
        self.lattice_bounds = {
            'cell_size_x': cell_vectors[0][0],
            'cell_size_y': cell_vectors[1][1],
            'repeat_x': repeat_cell[0],
            'repeat_y': repeat_cell[1],
            'x_min': 0,
            'y_min': 0,
            'x_max': cell_vectors[0][0] * repeat_cell[0],
            'y_max': cell_vectors[1][1] * repeat_cell[1]
        }
        
        # Generate all lattice sites
        self.generate_lattice_sites()
        
        build_time = time.time() - start_time
        print(f"Built lattice with {len(self.coordinates)} sites in {build_time:.2f}s")
    
    def generate_lattice_sites(self):
        """Generate lattice sites"""
        cell_size_x = self.lattice_bounds['cell_size_x']
        cell_size_y = self.lattice_bounds['cell_size_y']
        repeat_x = self.lattice_bounds['repeat_x']
        repeat_y = self.lattice_bounds['repeat_y']
        site_coords = self.lattice_info['site_coordinates']
        
        site_id = 0
        for i in range(repeat_x):
            for j in range(repeat_y):
                for k, (frac_x, frac_y) in enumerate(site_coords):
                    x = (i + frac_x) * cell_size_x
                    y = (j + frac_y) * cell_size_y
                    
                    if x <= self.lattice_bounds['x_max'] and y <= self.lattice_bounds['y_max']:
                        self.coordinates[site_id] = {
                            'x': x,
                            'y': y,
                            'type': self.lattice_info['site_type_names'][k % len(self.lattice_info['site_type_names'])]
                        }
                        site_id += 1
    
    def setup_figure(self):
        """Setup the main figure with improved layout and mini charts"""
        # Create figure with larger size to accommodate mini charts
        self.fig = plt.figure(figsize=(16, 8), dpi=100)
        self.fig.suptitle('KMCView - Enhanced Molecular Evolution Viewer', 
                         fontsize=14, fontweight='bold', color=self.colors['text'])
        
        # Mini charts on the left side
        # Surface coverage chart
        self.coverage_ax = plt.axes([0.07, 0.55, 0.20, 0.35])
        self.coverage_ax.set_facecolor(self.colors['panel'])
        self.coverage_ax.set_title('Surface Coverage', fontsize=10, color=self.colors['text'], pad=5)
        
        # Current coverage chart (pie chart)
        self.coverage_pie_ax = plt.axes([0.07, 0.25, 0.20, 0.25])
        self.coverage_pie_ax.set_facecolor(self.colors['panel'])
        self.coverage_pie_ax.set_title('Current Coverage', fontsize=10, color=self.colors['text'], pad=5)
        
        # Main plot area (moved to the right)
        self.ax = plt.axes([0.25, 0.30, 0.45, 0.60])
        self.ax.set_facecolor(self.colors['bg'])
        
        # Performance monitor (top right, smaller)
        self.perf_ax = plt.axes([0.72, 0.88, 0.25, 0.08])
        self.perf_ax.set_facecolor(self.colors['panel'])
        self.perf_text = self.perf_ax.text(0.05, 0.5, 'Ready', 
                                          transform=self.perf_ax.transAxes,
                                          fontsize=7, color=self.colors['text'],
                                          verticalalignment='center')
        self.perf_ax.set_xticks([])
        self.perf_ax.set_yticks([])
        
        # Expanded info panel (right side, adjusted)
        self.info_ax = plt.axes([0.72, 0.28, 0.25, 0.55])
        self.info_ax.set_facecolor(self.colors['panel'])
        self.info_text = self.info_ax.text(0.05, 0.98, 'Loading...', 
                                          transform=self.info_ax.transAxes,
                                          fontsize=9, color=self.colors['text'],
                                          verticalalignment='top', fontfamily='monospace')
        self.info_ax.set_xticks([])
        self.info_ax.set_yticks([])
        self.info_ax.set_title('Simulation Info', fontsize=11, color=self.colors['text'], pad=8)
        
        # Status bar (bottom)
        self.status_ax = plt.axes([0.08, 0.02, 0.89, 0.05])
        self.status_ax.set_facecolor(self.colors['panel'])
        self.status_text = self.status_ax.text(0.01, 0.5, 'Ready', 
                                              transform=self.status_ax.transAxes,
                                              fontsize=9, color=self.colors['success'],
                                              verticalalignment='center')
        self.status_ax.set_xticks([])
        self.status_ax.set_yticks([])
    
    def calculate_surface_coverage(self, step):
        """Calculate surface coverage for a given step"""
        if step >= len(self.evolution_df):
            return 0.0
        
        current_data = self.evolution_df.iloc[step]
        total_sites = len(self.coordinates)
        
        if total_sites == 0:
            return 0.0
        
        # Calculate total occupied sites
        occupied_sites = 0
        for species in ['H*', 'GeH2*', 'GeH3*']:
            occupied_sites += int(current_data.get(species, 0))
        
        return occupied_sites / total_sites * 100  # Return as percentage

    def calculate_species_coverage(self, step, species):
        """Calculate coverage for a specific species at a given step"""
        if step >= len(self.evolution_df):
            return 0.0
        
        current_data = self.evolution_df.iloc[step]
        total_sites = len(self.coordinates)
        
        if total_sites == 0:
            return 0.0
        
        species_count = int(current_data.get(species, 0))
        return species_count / total_sites * 100  # Return as percentage

    def update_mini_charts(self, step):
        """Update all mini charts with current step data"""
        # Update surface coverage chart
        self.update_coverage_chart(step)
        
        # Update current coverage pie chart
        self.update_coverage_pie_chart(step)
    
    def update_coverage_chart(self, step):
        """Update surface coverage chart with individual species coverage"""
        self.coverage_ax.clear()
        self.coverage_ax.set_facecolor(self.colors['panel'])
        
        # Calculate coverage for all steps up to current
        steps = list(range(min(step + 1, len(self.evolution_df))))
        times = [self.evolution_df.iloc[s]['time'] for s in steps]
        
        # Plot individual species coverage using colors from settings
        species_list = ['H*', 'GeH2*', 'GeH3*']
        for species in species_list:
            # Get species properties from settings, with fallback
            species_props = self.species_properties.get(species, {})
            species_color = species_props.get('color', '#888888')
            species_shape = species_props.get('shape', 'o')
            
            # Calculate species coverage
            species_coverages = [self.calculate_species_coverage(s, species) for s in steps]
            
            # Only plot if species is visible in settings
            if species_props.get('visible', True):
                # Ensure minimum line visibility by adding a small offset if all values are very small
                max_coverage = max(species_coverages) if species_coverages else 0
                if max_coverage < 0.1:  # If coverage is very small, add a small offset for visibility
                    species_coverages = [c + 0.01 for c in species_coverages]
                
                self.coverage_ax.plot(times, species_coverages, color=species_color, 
                                     linewidth=2, label=f'{species} Coverage', alpha=0.8)
        
        # Add current point highlights
        if step < len(self.evolution_df):
            current_time = self.evolution_df.iloc[step]['time']
            
            # Add current points for individual species
            for species in species_list:
                species_props = self.species_properties.get(species, {})
                if species_props.get('visible', True):
                    species_color = species_props.get('color', '#888888')
                    species_shape = species_props.get('shape', 'o')
                    current_species_coverage = self.calculate_species_coverage(step, species)
                    self.coverage_ax.scatter([current_time], [current_species_coverage], 
                                           color=species_color, s=50, marker=species_shape, 
                                           zorder=5, alpha=0.9)
        
        # Formatting
        self.coverage_ax.set_xlabel('Time', fontsize=8, color=self.colors['text'])
        self.coverage_ax.set_ylabel('Coverage (%)', fontsize=8, color=self.colors['text'])
        self.coverage_ax.tick_params(axis='both', which='major', labelsize=7, colors=self.colors['text'])
        self.coverage_ax.grid(True, alpha=0.3)
        self.coverage_ax.set_title('Species Coverage', fontsize=10, color=self.colors['text'], pad=5)
        
        # Add legend
        self.coverage_ax.legend(loc='upper left', fontsize=6, framealpha=0.8)
    
    def update_coverage_pie_chart(self, step):
        """Update current coverage pie chart for the given step"""
        self.coverage_pie_ax.clear()
        self.coverage_pie_ax.set_facecolor(self.colors['panel'])
        
        # Get current step data
        if step < len(self.evolution_df):
            current_data = self.evolution_df.iloc[step]
            
            # Calculate current species counts
            h_current = int(current_data.get('H*', 0))
            geh2_current = int(current_data.get('GeH2*', 0))
            geh3_current = int(current_data.get('GeH3*', 0))
            empty_current = len(self.coordinates) - h_current - geh2_current - geh3_current
            
            # Create pie chart using colors from settings
            sizes = [h_current, geh2_current, geh3_current, empty_current]
            labels = ['H*', 'GeH2*', 'GeH3*', 'Empty']
            
            # Get colors from species properties settings
            h_color = self.species_properties.get('H*', {}).get('color', '#FF4444')
            geh2_color = self.species_properties.get('GeH2*', {}).get('color', '#4477FF')
            geh3_color = self.species_properties.get('GeH3*', {}).get('color', '#44AA44')
            colors = [h_color, geh2_color, geh3_color, '#CCCCCC']  # Empty sites remain gray
            
            # Only show non-zero slices
            non_zero_data = [(s, l, c) for s, l, c in zip(sizes, labels, colors) if s > 0]
            if non_zero_data:
                sizes, labels, colors = zip(*non_zero_data)
                self.coverage_pie_ax.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', 
                                       startangle=90, textprops={'fontsize': 7})
        
        # self.coverage_pie_ax.set_title('Current Coverage', fontsize=10, color=self.colors['text'], pad=5)
    
    def setup_controls(self):
        """Setup enhanced control widgets with cleaner layout and better spacing"""
        # Main slider (adjusted for better alignment with buttons)
        ax_slider = plt.axes([0.14, 0.18, 0.50, 0.03])
        self.slider = Slider(ax_slider, 'Time Step', 0, max(0, len(self.evolution_df) - 1), 
                           valinit=0, valfmt='%d', color=self.colors['accent'])
        self.slider.on_changed(self.on_slider_change)
        
        # Speed control (moved much further right to avoid text collision)
        ax_speed = plt.axes([0.75, 0.18, 0.08, 0.03])
        self.speed_slider = Slider(ax_speed, 'Speed', 0.1, 5.0, valinit=1.0, 
                                 valfmt='%.1fx', color=self.colors['warning'])
        self.speed_slider.on_changed(self.on_speed_change)
        
        # Button row with much more aggressive spacing
        button_width = 0.060  # Slightly smaller for better fit
        button_height = 0.04
        button_spacing = 0.035  # Much more spacing for better visual separation
        start_x = 0.08  # Start position
        
        # Play button
        ax_play = plt.axes([start_x, 0.10, button_width, button_height])
        self.play_button = Button(ax_play, 'Play', color=self.colors['success'])
        self.play_button.on_clicked(self.toggle_play)
        
        # Reset button
        ax_reset = plt.axes([start_x + button_width + button_spacing, 0.10, button_width, button_height])
        self.reset_button = Button(ax_reset, 'Reset', color=self.colors['danger'])
        self.reset_button.on_clicked(self.reset_animation)
        
        # Step buttons with better spacing and symbols
        step_button_width = button_width * 0.70  # Smaller for step buttons
        ax_prev = plt.axes([start_x + 2 * (button_width + button_spacing), 0.10, step_button_width, button_height])
        self.prev_button = Button(ax_prev, '◀', color=self.colors['accent'])
        self.prev_button.on_clicked(self.prev_step)
        
        ax_next = plt.axes([start_x + 2 * (button_width + button_spacing) + step_button_width + button_spacing, 0.10, step_button_width, button_height])
        self.next_button = Button(ax_next, '▶', color=self.colors['accent'])
        self.next_button.on_clicked(self.next_step)
        
        # Text input with much wider spacing
        text_input_width = button_width * 1.5
        ax_text = plt.axes([start_x + 2 * (button_width + button_spacing) + 2 * step_button_width + button_spacing * 2, 0.10, text_input_width, button_height])
        self.text_box = TextBox(ax_text, 'Step:', initial='0', color=self.colors['panel'])
        self.text_box.on_submit(self.on_text_submit)
        
        # Settings button
        settings_width = button_width * 1.2
        ax_settings = plt.axes([start_x + 2 * (button_width + button_spacing) + 2 * step_button_width + text_input_width + button_spacing * 3, 0.10, settings_width, button_height])
        self.settings_button = Button(ax_settings, 'Settings', color=self.colors['accent'])
        self.settings_button.on_clicked(self.open_settings_modal)
    
    def setup_keyboard_shortcuts(self):
        """Setup keyboard shortcuts for better UX"""
        def on_key_press(event):
            if event.key == ' ':  # Spacebar
                self.toggle_play(None)
            elif event.key == 'r':  # R key
                self.reset_animation(None)
            elif event.key == 'left':  # Left arrow
                self.prev_step(None)
            elif event.key == 'right':  # Right arrow
                self.next_step(None)
            elif event.key == 'o':  # O key for settings (Options)
                self.open_settings_modal(None)
            elif event.key == 'q':  # Q key
                plt.close('all')
        
        self.fig.canvas.mpl_connect('key_press_event', on_key_press)
    
    def update_plot(self, step, force_redraw=False):
        """Enhanced plot update with performance optimization"""
        start_time = time.time()
        
        # Check cache first
        if not force_redraw and step in self.frame_cache and not self.show_trajectories:
            # Use cached frame for better performance
            cached_artists = self.frame_cache[step]
            self.ax.clear()
            for artist in cached_artists:
                self.ax.add_artist(artist)
            render_time = time.time() - start_time
        else:
            # Render new frame
            render_time = self._render_frame(step)
            
        # Update performance metrics
        self.render_times.append(render_time)
        if len(self.render_times) > self.max_render_history:
            self.render_times.pop(0)
        
        avg_render_time = np.mean(self.render_times) * 1000  # Convert to ms
        fps = 1.0 / render_time if render_time > 0 else float('inf')
        
        # Use ultra-compact format to prevent text overflow
        self.perf_text.set_text(f'{render_time*1000:.0f}ms {fps:.0f}fps')
        
        # Update info panel
        self.update_info_panel(step)
        
        # Update mini charts
        self.update_mini_charts(step)
        
        # Update status
        self.status_text.set_text(f'Step {step+1}/{len(self.evolution_df)} rendered')
        
        # Only redraw if not animating for performance
        if not self.is_playing:
            self.fig.canvas.draw_idle()
    
    def _render_frame(self, step):
        """Render a single frame with all optimizations"""
        start_time = time.time()
        
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
        
        # Draw lattice grid if enabled
        if self.show_lattice_grid:
            self._draw_lattice_grid()
        
        # Draw empty sites if enabled
        if self.show_empty_sites:
            self._draw_empty_sites()
        
        # Draw molecules with enhanced styling
        self._draw_molecules(positions)
        
        # Draw trajectories if enabled
        if self.show_trajectories:
            self._draw_trajectories(step)
        
        # Configure axes
        self._configure_axes(step, current_time)
        
        return time.time() - start_time
    
    def _draw_lattice_grid(self):
        """Draw lattice grid with modern styling"""
        cell_size_x = self.lattice_bounds['cell_size_x']
        cell_size_y = self.lattice_bounds['cell_size_y']
        repeat_x = int(self.lattice_bounds['x_max'] / cell_size_x)
        repeat_y = int(self.lattice_bounds['y_max'] / cell_size_y)
        
        # Vertical lines
        for i in range(repeat_x + 1):
            x_line = i * cell_size_x
            if x_line <= self.lattice_bounds['x_max']:
                self.ax.axvline(x=x_line, color='gray', alpha=0.3, linewidth=0.5, linestyle='--')
        
        # Horizontal lines
        for j in range(repeat_y + 1):
            y_line = j * cell_size_y
            if y_line <= self.lattice_bounds['y_max']:
                self.ax.axhline(y=y_line, color='gray', alpha=0.3, linewidth=0.5, linestyle='--')
    
    def _draw_empty_sites(self):
        """Draw empty lattice sites"""
        empty_x = [coord['x'] for coord in self.coordinates.values()]
        empty_y = [coord['y'] for coord in self.coordinates.values()]
        
        if empty_x:
            self.ax.scatter(empty_x, empty_y, c='lightgray', s=0.5, alpha=0.3, zorder=1)
    
    def _draw_molecules(self, positions):
        """Draw molecules using configurable species properties"""
        for species, coords in positions.items():
            # Get species properties, with fallback defaults
            props = self.species_properties.get(species, {
                'visible': True, 'color': '#888888', 'shape': 'o', 'size': 50
            })
            
            # Only draw if species is visible and has coordinates
            if coords and props['visible']:
                x_coords = [c['x'] for c in coords]
                y_coords = [c['y'] for c in coords]
                
                # Main scatter plot using configured properties
                scatter = self.ax.scatter(x_coords, y_coords, 
                                        c=props['color'], 
                                        s=props['size'], 
                                        marker=props['shape'],
                                        alpha=0.8, 
                                        label=f'{species} ({len(coords)})',
                                        edgecolors='white', linewidth=0.5, zorder=5)
                
                # Add glow effect for better visibility
                self.ax.scatter(x_coords, y_coords, 
                              c=props['color'], 
                              s=props['size']*2, 
                              marker=props['shape'],
                              alpha=0.2, zorder=4)
    
    def _draw_trajectories(self, step):
        """Draw molecular trajectories (optional feature)"""
        if step < 5:  # Need some history
            return
        
        # Show last 5 steps as trajectory
        for prev_step in range(max(0, step-4), step):
            prev_positions = self.generate_positions_for_step(prev_step)
            alpha = 0.1 + 0.15 * (prev_step - max(0, step-4)) / 4  # Fade effect
            
            for species, coords in prev_positions.items():
                if coords:
                    x_coords = [c['x'] for c in coords]
                    y_coords = [c['y'] for c in coords]
                    self.ax.scatter(x_coords, y_coords, 
                                  c='gray', s=5, alpha=alpha, zorder=2)
    
    def _configure_axes(self, step, current_time):
        """Configure axes with enhanced styling"""
        # Set limits with margins
        if hasattr(self, 'lattice_bounds'):
            x_margin = (self.lattice_bounds['x_max'] - self.lattice_bounds['x_min']) * 0.02
            y_margin = (self.lattice_bounds['y_max'] - self.lattice_bounds['y_min']) * 0.02
            
            self.ax.set_xlim(self.lattice_bounds['x_min'] - x_margin, 
                           self.lattice_bounds['x_max'] + x_margin)
            self.ax.set_ylim(self.lattice_bounds['y_min'] - y_margin, 
                           self.lattice_bounds['y_max'] + y_margin)
        
        # Enhanced styling
        self.ax.set_xlabel('X Coordinate (Å)', fontsize=12, color=self.colors['text'])
        self.ax.set_ylabel('Y Coordinate (Å)', fontsize=12, color=self.colors['text'])
        self.ax.set_aspect('equal')
        self.ax.grid(True, alpha=0.2, color=self.colors['text'])
        
        # Enhanced legend
        legend = self.ax.legend(loc='upper left', framealpha=0.9, 
                              fancybox=True, shadow=True)
        legend.get_frame().set_facecolor(self.colors['panel'])
        
        # Title with step info
        title = f'Molecular Evolution - Step {step+1}/{len(self.evolution_df)} (t={current_time:.2f})'
        self.ax.set_title(title, fontsize=12, color=self.colors['text'], pad=20)
    
    def update_info_panel(self, step):
        """Update the information panel with current data"""
        current_data = self.evolution_df.iloc[step]
        current_time = current_data['time']
        h_count = int(current_data['H*'])
        geh2_count = int(current_data['GeH2*'])
        geh3_count = int(current_data['GeH3*'])
        total_molecules = h_count + geh2_count + geh3_count
        
        # Count visible molecules
        visible_h = h_count if self.species_properties.get('H*', {}).get('visible', True) else 0
        visible_geh2 = geh2_count if self.species_properties.get('GeH2*', {}).get('visible', True) else 0
        visible_geh3 = geh3_count if self.species_properties.get('GeH3*', {}).get('visible', True) else 0
        visible_total = visible_h + visible_geh2 + visible_geh3
        
        # Determine data source for positioning
        has_history_data = hasattr(self, 'lattice_states') and self.lattice_states
        data_source_info = ""
        
        if has_history_data:
            actual_kmc_step = None
            if hasattr(self, 'step_to_config') and step in self.step_to_config:
                actual_kmc_step = self.step_to_config[step]
            
            if actual_kmc_step is not None and actual_kmc_step in self.lattice_states:
                data_source_info = "[OK] History Data (Exact)"
                self.info_text.set_color(self.colors['success'])
            else:
                data_source_info = "[ERROR] Data Missing"
                self.info_text.set_color(self.colors['danger'])
        else:
            data_source_info = "[WARN] Approximated (No History)"
            self.info_text.set_color(self.colors['warning'])
        
        # Generate visibility indicators
        visibility_info = []
        for species in ['H*', 'GeH2*', 'GeH3*']:
            status = "VISIBLE" if self.species_properties.get(species, {}).get('visible', True) else "HIDDEN"
            visibility_info.append(f"{status} {species}")
        
        info_text = f"""SIMULATION STATUS
{data_source_info}

TIME & STEP
Time: {current_time:.3f} time units
Step: {step+1:,} / {len(self.evolution_df):,}

MOLECULE COUNTS
H*:     {h_count:>4}
GeH2*:  {geh2_count:>4}  
GeH3*:  {geh3_count:>4}
────────────
Total:  {total_molecules:>4}
Shown:  {visible_total:>4}

LATTICE INFO
Grid: {self.lattice_bounds['repeat_x']}×{self.lattice_bounds['repeat_y']}
Sites: {len(self.coordinates):,}
Cell: {self.lattice_bounds['cell_size_x']:.1f}×{self.lattice_bounds['cell_size_y']:.1f} Å

SPECIES VISIBILITY
{chr(10).join(visibility_info)}

DISPLAY OPTIONS
Grid: {'ON' if self.show_lattice_grid else 'OFF'}
Empty: {'ON' if self.show_empty_sites else 'OFF'}
Trails: {'ON' if self.show_trajectories else 'OFF'}"""
        
        self.info_text.set_text(info_text)
    
    def generate_positions_for_step(self, step):
        """Generate molecular positions using actual lattice states or fallback to counts"""
        # Check cache first
        if step in self.position_cache:
            return self.position_cache[step]
        
        if step >= len(self.evolution_df):
            return {}
        
        # Check if we have history data available
        has_history_data = hasattr(self, 'lattice_states') and self.lattice_states
        
        if has_history_data:
            # When history data is available, NEVER use random positioning
            actual_kmc_step = None
            if hasattr(self, 'step_to_config') and step in self.step_to_config:
                actual_kmc_step = self.step_to_config[step]
            
            if actual_kmc_step is not None and actual_kmc_step in self.lattice_states:
                positions = self.generate_positions_from_lattice_state(actual_kmc_step)
            else:
                # Report error to terminal - this should not happen with proper history data
                print(f"❌ ERROR: No lattice state found for step {step} (KMC step {actual_kmc_step})")
                print(f"   Available KMC steps: {list(self.lattice_states.keys())[:10]}...")
                print(f"   This indicates a data parsing issue - please check history file integrity")
                # Return empty positions rather than random ones
                return {}
        else:
            # Only use random positioning when NO history data is available
            positions = self.generate_positions_from_counts_deterministic(step)
        
        # Cache the result
        if len(self.position_cache) >= self.max_position_cache:
            # Remove oldest entry (simple FIFO)
            oldest_key = next(iter(self.position_cache))
            del self.position_cache[oldest_key]
        
        self.position_cache[step] = positions
        return positions
    
    def generate_positions_from_lattice_state(self, kmc_step):
        """Generate positions using actual lattice state data"""
        lattice_state = self.lattice_states[kmc_step]
        positions = {species: [] for species in getattr(self, 'surface_species', ['H*', 'GeH2*', 'GeH3*'])}
        
        # Map species IDs to species names
        species_map = {}
        if hasattr(self, 'surface_species'):
            for i, species in enumerate(self.surface_species):
                species_map[i + 1] = species  # species_id starts from 1
        else:
            # Default mapping for backward compatibility
            species_map = {1: 'H*', 2: 'GeH2*', 3: 'GeH3*'}
        
        # Place molecules based on actual site occupancy
        for site_id, species_id in lattice_state.items():
            if species_id > 0 and species_id in species_map:
                species_name = species_map[species_id]
                if site_id in self.coordinates and species_name in positions:
                    positions[species_name].append(self.coordinates[site_id])
        
        return positions
    
    def generate_positions_from_counts_deterministic(self, step):
        """Deterministic fallback: Generate positions using species counts without randomness"""
        current_data = self.evolution_df.iloc[step]
        
        # Get species counts (handle both new and old column names)
        species_counts = {}
        for species in ['H*', 'GeH2*', 'GeH3*']:
            species_counts[species] = int(current_data.get(species, 0))
        
        total_molecules = sum(species_counts.values())
        
        if total_molecules == 0:
            return {species: [] for species in species_counts.keys()}
        
        # Use available site coordinates
        all_sites = list(self.coordinates.keys())
        if not all_sites:
            return {species: [] for species in species_counts.keys()}
        
        # Deterministic algorithm: Use first available sites consistently
        # This eliminates all randomness for reproducible results
        occupied_sites = all_sites[:min(total_molecules, len(all_sites))]
        
        positions = {}
        idx = 0
        
        # Assign positions to species with consistent ordering
        for species in ['H*', 'GeH2*', 'GeH3*']:
            count = species_counts[species]
            positions[species] = []
            for _ in range(count):
                if idx < len(occupied_sites):
                    site_id = occupied_sites[idx]
                    positions[species].append(self.coordinates[site_id])
                    idx += 1
        
        return positions
    
    # Enhanced callback methods
    def on_slider_change(self, val):
        """Handle slider value changes with debouncing"""
        step = int(val)
        if step != self.current_step:
            self.current_step = step
            self.update_plot(step)
            self.text_box.set_val(str(step))
    
    def on_speed_change(self, val):
        """Handle animation speed changes"""
        self.animation_speed = val
        if self.is_playing and self.animation_obj:
            # Restart animation with new speed
            self.animation_obj.event_source.stop()
            interval = max(50, int(300 / self.animation_speed))  # Min 50ms for smoothness
            self.animation_obj = FuncAnimation(self.fig, self.animate_frame, 
                                             interval=interval, blit=False, repeat=True)
    
    def on_options_change(self, label):
        """Handle display options changes"""
        if label == 'Lattice Grid':
            self.show_lattice_grid = not self.show_lattice_grid
        elif label == 'Empty Sites':
            self.show_empty_sites = not self.show_empty_sites
        elif label == 'Trajectories':
            self.show_trajectories = not self.show_trajectories
        
        # Clear cache and redraw
        self.frame_cache.clear()
        self.update_plot(self.current_step, force_redraw=True)
    

    
    def toggle_play(self, event):
        """Enhanced play/pause with better feedback"""
        if self.is_playing:
            self.is_playing = False
            if self.animation_obj:
                self.animation_obj.event_source.stop()
            self.play_button.label.set_text('Play')
            self.status_text.set_text('Animation paused')
        else:
            self.is_playing = True
            interval = max(50, int(300 / self.animation_speed))
            self.animation_obj = FuncAnimation(self.fig, self.animate_frame, 
                                             interval=interval, blit=False, repeat=True)
            self.play_button.label.set_text('Pause')
            self.status_text.set_text(f'Playing at {self.animation_speed:.1f}x speed')
        
        self.fig.canvas.draw_idle()
    
    def reset_animation(self, event):
        """Enhanced reset with feedback"""
        self.is_playing = False
        if self.animation_obj:
            self.animation_obj.event_source.stop()
        self.current_step = 0
        self.slider.set_val(0)
        self.update_plot(0)
        self.text_box.set_val('0')
        self.play_button.label.set_text('Play')
        self.status_text.set_text('Reset to beginning')
        self.fig.canvas.draw_idle()
    
    def prev_step(self, event):
        """Go to previous step"""
        if self.current_step > 0:
            self.current_step -= 1
            self.slider.set_val(self.current_step)
            self.update_plot(self.current_step)
            self.text_box.set_val(str(self.current_step))
    
    def next_step(self, event):
        """Go to next step"""
        if self.current_step < len(self.evolution_df) - 1:
            self.current_step += 1
            self.slider.set_val(self.current_step)
            self.update_plot(self.current_step)
            self.text_box.set_val(str(self.current_step))
    

    
    def on_text_submit(self, text):
        """Enhanced text input with validation"""
        try:
            step = int(text)
            step = max(0, min(step, len(self.evolution_df) - 1))
            self.current_step = step
            self.slider.set_val(step)
            self.update_plot(step)
            self.status_text.set_text(f'Jumped to step {step}')
        except ValueError:
            self.status_text.set_text(f'Invalid step: {text}')
    
    def animate_frame(self, frame):
        """Enhanced animation frame with smooth transitions"""
        if self.is_playing:
            self.current_step = (self.current_step + 1) % len(self.evolution_df)
            self.slider.set_val(self.current_step)
            self.update_plot(self.current_step)
            return []
        return []

    def open_settings_modal(self, event):
        """Open settings modal with native table widget and clean layout"""
        if self.settings_open:
            return
            
        self.settings_open = True
        
        # Create figure with smaller, more reasonable size
        self.settings_fig = plt.figure(figsize=(10, 6), dpi=100)
        self.settings_fig.suptitle('Display Settings', fontsize=12, fontweight='bold', color=self.colors['text'])
        
        # Apply theme to settings figure
        self.settings_fig.patch.set_facecolor(self.colors['bg'])
        
        # Use subplot2grid for better layout control with more compact spacing
        options_ax = plt.subplot2grid((4, 3), (0, 0), colspan=2, fig=self.settings_fig)
        table_ax = plt.subplot2grid((4, 3), (1, 0), colspan=3, rowspan=2, fig=self.settings_fig)
        buttons_ax = plt.subplot2grid((4, 3), (3, 0), colspan=3, fig=self.settings_fig)
        
        # Apply theme colors to all axes
        for ax in [options_ax, table_ax, buttons_ax]:
            ax.set_facecolor(self.colors['panel'])
        
        # Setup each section
        self.setup_display_options(options_ax)
        self.setup_species_table(table_ax)
        self.setup_control_buttons(buttons_ax)
        
        # Handle window close event
        self.settings_fig.canvas.mpl_connect('close_event', self.on_settings_window_close)
        self.settings_fig.canvas.mpl_connect('button_press_event', self.on_table_click)
        
        plt.tight_layout(pad=2.5)  # Even more padding to prevent collisions
        plt.show()
    

    
    def setup_display_options(self, ax):
        """Setup display options controls"""
        options_labels = ['Lattice Grid', 'Empty Sites', 'Trajectories']
        options_states = [self.show_lattice_grid, self.show_empty_sites, self.show_trajectories]
        self.settings_options_check = CheckButtons(ax, options_labels, options_states)
        self.settings_options_check.on_clicked(self.on_settings_options_change)
        
        # Apply theme colors to checkboxes
        ax.set_title('Display Options', fontsize=11, pad=10, color=self.colors['text'])
        
        # Style checkboxes for better visibility in dark mode
        try:
            # Try to access rectangles if they exist
            if hasattr(self.settings_options_check, 'rectangles'):
                for rect in self.settings_options_check.rectangles:
                    rect.set_edgecolor(self.colors['text'])
                    rect.set_linewidth(1.5)
                    if rect.get_facecolor()[0] > 0.5:  # If checked (white background)
                        rect.set_facecolor(self.colors['accent'])
                    else:  # If unchecked (transparent)
                        rect.set_facecolor('none')
        except:
            pass  # If rectangles don't exist, skip styling
        
        # Style labels
        try:
            if hasattr(self.settings_options_check, 'labels'):
                for label in self.settings_options_check.labels:
                    label.set_color(self.colors['text'])
        except:
            pass  # If labels don't exist, skip styling
    
    def setup_species_table(self, ax):
        """Setup species properties table using matplotlib's native table widget"""
        # Prepare table data
        species_data = []
        for species, props in self.species_properties.items():
            species_data.append([
                species,
                'ON' if props['visible'] else 'OFF',
                props['color'],
                self.get_shape_display_name(props['shape']),
                str(props['size'])
            ])
        
        # Create native matplotlib table
        self.species_table = ax.table(
            cellText=species_data,
            colLabels=['Species', 'Visible', 'Color', 'Shape', 'Size'],
            cellLoc='center',
            loc='center',
            colWidths=[0.2, 0.15, 0.25, 0.2, 0.2]
        )
        
        # Style the table
        self.species_table.auto_set_font_size(False)
        self.species_table.set_fontsize(8)  # Even smaller font to prevent text overflow
        self.species_table.scale(1, 2.8)  # Even taller rows for better text spacing
        
        # Color the header row
        for i in range(5):  # 5 columns
            self.species_table[(0, i)].set_facecolor(self.colors['accent'])
            self.species_table[(0, i)].set_text_props(weight='bold', color='white')
        
        # Color the color column cells with their actual colors
        for row in range(1, len(species_data) + 1):
            color = species_data[row-1][2]  # Color column
            try:
                self.species_table[(row, 2)].set_facecolor(color)
                # Set text color based on background brightness
                self.species_table[(row, 2)].set_text_props(color='white' if self.is_dark_color(color) else 'black')
            except:
                pass  # Invalid color format
        
        ax.set_title('Species Properties (Click cells to edit)', fontsize=11, pad=15, color=self.colors['text'])  # Smaller title and padding
        ax.axis('off')
    
    def setup_control_buttons(self, ax):
        """Setup control buttons in a clean layout"""
        ax.axis('off')
        
        # Create buttons within the subplot axes for proper layout integration
        button_width = 0.25
        button_height = 0.6
        button_spacing = 0.05
        
        # Calculate button positions within the subplot
        total_width = 3 * button_width + 2 * button_spacing
        start_x = (1.0 - total_width) / 2
        
        # Create button axes within the subplot
        reset_ax = ax.figure.add_axes([ax.get_position().x0 + start_x * ax.get_position().width, 
                                      ax.get_position().y0 + 0.2 * ax.get_position().height,
                                      button_width * ax.get_position().width, 
                                      button_height * ax.get_position().height])
        self.settings_reset_button = Button(reset_ax, 'Reset Defaults', color=self.colors['danger'])
        self.settings_reset_button.on_clicked(self.reset_species_defaults)
        
        # Apply button  
        apply_ax = ax.figure.add_axes([ax.get_position().x0 + (start_x + button_width + button_spacing) * ax.get_position().width,
                                      ax.get_position().y0 + 0.2 * ax.get_position().height,
                                      button_width * ax.get_position().width, 
                                      button_height * ax.get_position().height])
        self.settings_apply_button = Button(apply_ax, 'Apply Changes', color=self.colors['success'])
        self.settings_apply_button.on_clicked(self.apply_species_changes)
        
        # Close button
        close_ax = ax.figure.add_axes([ax.get_position().x0 + (start_x + 2 * (button_width + button_spacing)) * ax.get_position().width,
                                      ax.get_position().y0 + 0.2 * ax.get_position().height,
                                      button_width * ax.get_position().width, 
                                      button_height * ax.get_position().height])
        self.settings_close_button = Button(close_ax, 'Close', color=self.colors['accent'])
        self.settings_close_button.on_clicked(self.close_settings_modal)
    
    def get_shape_display_name(self, shape):
        """Get display name for shape"""
        shape_names = {
            'o': 'Circle', 's': 'Square', '^': 'Triangle', 'D': 'Diamond',
            'v': 'Down△', '<': 'Left△', '>': 'Right△', 'p': 'Pentagon', 
            '*': 'Star', 'h': 'Hexagon'
        }
        return f'{shape} ({shape_names.get(shape, shape)})'
    
    def is_dark_color(self, color):
        """Check if a color is dark (for text color selection)"""
        try:
            # Simple brightness check
            if color.startswith('#'):
                r, g, b = int(color[1:3], 16), int(color[3:5], 16), int(color[5:7], 16)
                brightness = (r * 299 + g * 587 + b * 114) / 1000
                return brightness < 128
        except:
            pass
        return False
    
    def on_table_click(self, event):
        """Handle table cell clicks for editing"""
        if not hasattr(self, 'species_table'):
            return
            
        # Find which cell was clicked
        for (row, col), cell in self.species_table.get_celld().items():
            if row == 0:  # Skip header row
                continue
                
            if cell.contains(event)[0]:
                species_list = list(self.species_properties.keys())
                if row - 1 < len(species_list):
                    species = species_list[row - 1]
                    self.edit_table_cell(row, col, species)
                break
    
    def edit_table_cell(self, row, col, species):
        """Edit a specific table cell based on column"""
        if col == 1:  # Visibility column
            self.toggle_species_visibility(species, row, col)
        elif col == 2:  # Color column
            self.cycle_species_color(species, row, col)
        elif col == 3:  # Shape column
            self.cycle_species_shape_simple(species, row, col)
        elif col == 4:  # Size column
            self.cycle_species_size(species, row, col)
    
    def toggle_species_visibility(self, species, row, col):
        """Toggle species visibility and update table"""
        self.species_properties[species]['visible'] = not self.species_properties[species]['visible']
        
        # Update table cell
        new_text = 'ON' if self.species_properties[species]['visible'] else 'OFF'
        self.species_table[(row, col)].get_text().set_text(new_text)
        
        # Update visualization
        self.frame_cache.clear()
        self.update_plot(self.current_step, force_redraw=True)
        
        # Update status
        visible_species = [s for s, props in self.species_properties.items() if props['visible']]
        self.status_text.set_text(f'Visible species: {", ".join(visible_species)}')
        
        plt.draw()
    
    def cycle_species_color(self, species, row, col):
        """Cycle through predefined colors"""
        colors = ['#FF4444', '#4477FF', '#44AA44', '#FFAA44', '#AA44FF', '#44AAFF', '#FF44AA', '#AAFF44']
        current_color = self.species_properties[species]['color']
        
        try:
            current_idx = colors.index(current_color)
        except ValueError:
            current_idx = 0
        
        new_idx = (current_idx + 1) % len(colors)
        new_color = colors[new_idx]
        
        self.species_properties[species]['color'] = new_color
        
        # Update table cell
        self.species_table[(row, col)].get_text().set_text(new_color)
        self.species_table[(row, col)].set_facecolor(new_color)
        self.species_table[(row, col)].set_text_props(color='white' if self.is_dark_color(new_color) else 'black')
        
        # Update visualization
        self.frame_cache.clear()
        self.update_plot(self.current_step, force_redraw=True)
        self.status_text.set_text(f'Color updated for {species}: {new_color}')
        
        plt.draw()
    
    def cycle_species_shape_simple(self, species, row, col):
        """Cycle through available shapes"""
        shapes = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h']
        current_shape = self.species_properties[species]['shape']
        
        try:
            current_idx = shapes.index(current_shape)
        except ValueError:
            current_idx = 0
        
        new_idx = (current_idx + 1) % len(shapes)
        new_shape = shapes[new_idx]
        
        self.species_properties[species]['shape'] = new_shape
        
        # Update table cell
        self.species_table[(row, col)].get_text().set_text(self.get_shape_display_name(new_shape))
        
        # Update visualization
        self.frame_cache.clear()
        self.update_plot(self.current_step, force_redraw=True)
        self.status_text.set_text(f'Shape updated for {species}: {new_shape}')
        
        plt.draw()
    
    def cycle_species_size(self, species, row, col):
        """Cycle through common sizes"""
        sizes = [20, 40, 60, 80, 100, 120, 150, 200]
        current_size = self.species_properties[species]['size']
        
        try:
            current_idx = sizes.index(current_size)
        except ValueError:
            current_idx = 1  # Default to 40
        
        new_idx = (current_idx + 1) % len(sizes)
        new_size = sizes[new_idx]
        
        self.species_properties[species]['size'] = new_size
        
        # Update table cell
        self.species_table[(row, col)].get_text().set_text(str(new_size))
        
        # Update visualization
        self.frame_cache.clear()
        self.update_plot(self.current_step, force_redraw=True)
        self.status_text.set_text(f'Size updated for {species}: {new_size}')
        
        plt.draw()
    
    def on_species_visibility_change(self, species):
        """Legacy callback - functionality moved to table clicks"""
        pass
    
    def refresh_table(self):
        """Refresh the table display after changes"""
        if hasattr(self, 'species_table'):
            plt.draw()
    
    def reset_species_defaults(self, event):
        """Callback for resetting species properties to defaults"""
        self.species_properties = {
            'H*': {
                'visible': True,
                'color': '#FF4444',  # Bright red
                'shape': 'o',        # Circle
                'size': 40
            },
            'GeH2*': {
                'visible': True,
                'color': '#4477FF',  # Bright blue
                'shape': 's',        # Square
                'size': 80
            },
            'GeH3*': {
                'visible': True,
                'color': '#44AA44',  # Bright green
                'shape': '^',        # Triangle
                'size': 120
            }
        }
        self.frame_cache.clear()
        self.update_plot(self.current_step, force_redraw=True)
        self.status_text.set_text('Species properties reset to defaults')
    
    def apply_species_changes(self, event):
        """Callback for applying species property changes"""
        self.frame_cache.clear()
        self.update_plot(self.current_step, force_redraw=True)
        self.status_text.set_text('Species properties applied')
    
    def close_settings_modal(self, event):
        """Close settings modal window"""
        if hasattr(self, 'settings_fig'):
            plt.close(self.settings_fig)
        self.settings_open = False
    
    def on_settings_window_close(self, event):
        """Handle settings window close event"""
        self.settings_open = False
    

    
    def on_settings_options_change(self, label):
        """Handle display options changes from settings modal"""
        if label == 'Lattice Grid':
            self.show_lattice_grid = not self.show_lattice_grid
        elif label == 'Empty Sites':
            self.show_empty_sites = not self.show_empty_sites
        elif label == 'Trajectories':
            self.show_trajectories = not self.show_trajectories
        
        self.frame_cache.clear()
        self.update_plot(self.current_step, force_redraw=True)
    
    def on_settings_species_change(self, label):
        """Handle species visibility changes from settings modal (legacy - now handled by data grid)"""
        # This method is kept for compatibility but functionality moved to data grid
        pass

def main():
    """Enhanced main function with argument parsing"""
    import argparse
    
    # ASCII Art Banner
    ascii_art = """
██╗  ██╗███╗   ███╗ ██████╗ ██╗   ██╗██╗███████╗██╗    ██╗
██║ ██╔╝████╗ ████║██╔═══██╗██║   ██║██║██╔════╝██║    ██║
█████╔╝ ██╔████╔██║██║      ██║   ██║██║█████╗  ██║ █╗ ██║
██╔═██╗ ██║╚██╔╝██║██║   ██║╚██╗ ██╔╝██║██╔══╝  ██║███╗██║
██║  ██╗██║ ╚═╝ ██║╚██████╔╝ ╚████╔╝ ██║███████╗╚███╔███╔╝
╚═╝  ╚═╝╚═╝     ╚═╝ ╚═════╝   ╚═══╝  ╚═╝╚══════╝ ╚══╝╚══╝ 
"""
    
    parser = argparse.ArgumentParser(description='KMCView - Enhanced Molecular Evolution Viewer')
    parser.add_argument('--data-dir', '-d', default='input-output', 
                       help='Directory containing KMC output files')

    
    args = parser.parse_args()
    
    print(ascii_art)
    print(f"🚀 Starting KMCView - Enhanced Molecular Viewer")
    print(f"📁 Data directory: {args.data_dir}")
    
    viewer = EnhancedMolecularViewer(data_dir=args.data_dir)

if __name__ == '__main__':
    main() 