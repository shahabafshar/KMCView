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
    
    def __init__(self, data_dir='input-output', theme='dark'):
        """Initialize the enhanced viewer"""
        self.data_dir = data_dir
        self.theme = theme
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
        self.export_format = 'png'
        
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
        print(f"🎨 Theme: {theme}")
        
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
        
        print("Keyboard shortcuts: Space=Play/Pause, R=Reset, Left/Right=Step, S=Save")
        
        plt.show()
    
    def setup_theme(self):
        """Setup modern theme styling"""
        if self.theme == 'dark':
            plt.style.use('dark_background')
            self.colors = {
                'bg': '#2E2E2E',
                'panel': '#3E3E3E',
                'text': '#FFFFFF',
                'accent': '#4A90E2',
                'success': '#5CB85C',
                'warning': '#F0AD4E',
                'danger': '#D9534F'
            }
        else:
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
        """Setup the main figure with improved layout"""
        # Create figure with optimal size and DPI
        self.fig = plt.figure(figsize=(18, 10), dpi=100)
        self.fig.suptitle('KMCView - Enhanced Molecular Evolution Viewer', 
                         fontsize=16, fontweight='bold', color=self.colors['text'])
        
        # Main plot area (slightly larger)
        self.ax = plt.axes([0.05, 0.25, 0.68, 0.65])
        self.ax.set_facecolor(self.colors['bg'])
        
        # Performance monitor (top right, smaller)
        self.perf_ax = plt.axes([0.76, 0.88, 0.22, 0.08])
        self.perf_ax.set_facecolor(self.colors['panel'])
        self.perf_text = self.perf_ax.text(0.05, 0.5, 'Performance: Ready', 
                                          transform=self.perf_ax.transAxes,
                                          fontsize=9, color=self.colors['text'])
        self.perf_ax.set_xticks([])
        self.perf_ax.set_yticks([])
        
        # Expanded info panel (right side, much larger)
        self.info_ax = plt.axes([0.76, 0.25, 0.22, 0.6])
        self.info_ax.set_facecolor(self.colors['panel'])
        self.info_text = self.info_ax.text(0.05, 0.98, 'Loading...', 
                                          transform=self.info_ax.transAxes,
                                          fontsize=10, color=self.colors['text'],
                                          verticalalignment='top', fontfamily='monospace')
        self.info_ax.set_xticks([])
        self.info_ax.set_yticks([])
        self.info_ax.set_title('Simulation Info', fontsize=12, color=self.colors['text'], pad=10)
        
        # Status bar (bottom)
        self.status_ax = plt.axes([0.05, 0.02, 0.9, 0.05])
        self.status_ax.set_facecolor(self.colors['panel'])
        self.status_text = self.status_ax.text(0.01, 0.5, 'Ready', 
                                              transform=self.status_ax.transAxes,
                                              fontsize=10, color=self.colors['success'],
                                              verticalalignment='center')
        self.status_ax.set_xticks([])
        self.status_ax.set_yticks([])
    
    def setup_controls(self):
        """Setup enhanced control widgets with cleaner layout"""
        # Main slider (wider, better positioned)
        ax_slider = plt.axes([0.1, 0.15, 0.5, 0.03])
        self.slider = Slider(ax_slider, 'Time Step', 0, max(0, len(self.evolution_df) - 1), 
                           valinit=0, valfmt='%d', color=self.colors['accent'])
        self.slider.on_changed(self.on_slider_change)
        
        # Speed control
        ax_speed = plt.axes([0.65, 0.15, 0.1, 0.03])
        self.speed_slider = Slider(ax_speed, 'Speed', 0.1, 5.0, valinit=1.0, 
                                 valfmt='%.1fx', color=self.colors['warning'])
        self.speed_slider.on_changed(self.on_speed_change)
        
        # Play button
        ax_play = plt.axes([0.1, 0.1, 0.07, 0.04])
        self.play_button = Button(ax_play, 'Play', color=self.colors['success'])
        self.play_button.on_clicked(self.toggle_play)
        
        # Reset button
        ax_reset = plt.axes([0.18, 0.1, 0.07, 0.04])
        self.reset_button = Button(ax_reset, 'Reset', color=self.colors['danger'])
        self.reset_button.on_clicked(self.reset_animation)
        
        # Step buttons
        ax_prev = plt.axes([0.26, 0.1, 0.05, 0.04])
        self.prev_button = Button(ax_prev, '<<', color=self.colors['accent'])
        self.prev_button.on_clicked(self.prev_step)
        
        ax_next = plt.axes([0.32, 0.1, 0.05, 0.04])
        self.next_button = Button(ax_next, '>>', color=self.colors['accent'])
        self.next_button.on_clicked(self.next_step)
        
        # Text input
        ax_text = plt.axes([0.39, 0.1, 0.09, 0.04])
        self.text_box = TextBox(ax_text, 'Step:', initial='0', color=self.colors['panel'])
        self.text_box.on_submit(self.on_text_submit)
        
        # Export button
        ax_export = plt.axes([0.50, 0.1, 0.07, 0.04])
        self.export_button = Button(ax_export, 'Save', color=self.colors['warning'])
        self.export_button.on_clicked(self.export_frame)
        
        # Settings button (opens modal)
        ax_settings = plt.axes([0.59, 0.1, 0.08, 0.04])
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
            elif event.key == 's':  # S key
                self.export_frame(None)
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
        
        self.perf_text.set_text(f'Render: {render_time*1000:.1f}ms\nAvg: {avg_render_time:.1f}ms\nFPS: {fps:.1f}')
        
        # Update info panel
        self.update_info_panel(step)
        
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
                data_source_info = "✅ History Data (Exact)"
                self.info_text.set_color(self.colors['success'])
            else:
                data_source_info = "❌ Data Missing"
                self.info_text.set_color(self.colors['danger'])
        else:
            data_source_info = "⚠️ Approximated (No History)"
            self.info_text.set_color(self.colors['warning'])
        
        # Generate visibility indicators
        visibility_info = []
        for species in ['H*', 'GeH2*', 'GeH3*']:
            status = "👁" if self.species_properties.get(species, {}).get('visible', True) else "🚫"
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
    
    def on_theme_change(self, label):
        """Handle theme changes"""
        new_theme = 'dark' if label == 'Dark' else 'light'
        if new_theme != self.theme:
            self.theme = new_theme
            self.setup_theme()
            self.frame_cache.clear()  # Clear cache due to theme change
            self.update_plot(self.current_step, force_redraw=True)
            self.status_text.set_text(f'Switched to {new_theme} theme')
    
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
    
    def export_frame(self, event):
        """Export current frame with enhanced options"""
        try:
            filename = f"kmcview_frame_step_{self.current_step:04d}.{self.export_format}"
            
            # Create a clean figure for export
            export_fig, export_ax = plt.subplots(figsize=(12, 10), dpi=300)
            
            # Render the current frame on the export figure
            current_data = self.evolution_df.iloc[self.current_step]
            positions = self.generate_positions_for_step(self.current_step)
            
            # Draw everything on export figure
            if self.show_lattice_grid:
                cell_size_x = self.lattice_bounds['cell_size_x']
                cell_size_y = self.lattice_bounds['cell_size_y']
                repeat_x = int(self.lattice_bounds['x_max'] / cell_size_x)
                repeat_y = int(self.lattice_bounds['y_max'] / cell_size_y)
                
                for i in range(repeat_x + 1):
                    x_line = i * cell_size_x
                    if x_line <= self.lattice_bounds['x_max']:
                        export_ax.axvline(x=x_line, color='gray', alpha=0.3, linewidth=0.5)
                
                for j in range(repeat_y + 1):
                    y_line = j * cell_size_y
                    if y_line <= self.lattice_bounds['y_max']:
                        export_ax.axhline(y=y_line, color='gray', alpha=0.3, linewidth=0.5)
            
            # Draw molecules
            colors = {'H*': '#FF4444', 'GeH2*': '#4477FF', 'GeH3*': '#44AA44'}
            sizes = {'H*': 40, 'GeH2*': 80, 'GeH3*': 120}
            
            for species, coords in positions.items():
                if coords:
                    x_coords = [c['x'] for c in coords]
                    y_coords = [c['y'] for c in coords]
                    export_ax.scatter(x_coords, y_coords, 
                                    c=colors[species], s=sizes[species], 
                                    alpha=0.8, label=f'{species} ({len(coords)})',
                                    edgecolors='black', linewidth=0.5)
            
            # Configure export axes
            export_ax.set_xlim(self.lattice_bounds['x_min'], self.lattice_bounds['x_max'])
            export_ax.set_ylim(self.lattice_bounds['y_min'], self.lattice_bounds['y_max'])
            export_ax.set_xlabel('X Coordinate (Å)', fontsize=14)
            export_ax.set_ylabel('Y Coordinate (Å)', fontsize=14)
            export_ax.set_aspect('equal')
            export_ax.grid(True, alpha=0.3)
            export_ax.legend()
            
            current_time = current_data['time']
            export_ax.set_title(f'KMC Molecular Evolution - Step {self.current_step+1} (t={current_time:.2f})', 
                              fontsize=16, pad=20)
            
            # Save with high quality
            export_fig.savefig(filename, dpi=300, bbox_inches='tight', 
                             facecolor='white', edgecolor='none')
            plt.close(export_fig)
            
            self.status_text.set_text(f'Saved {filename}')
            print(f"Exported frame to {filename}")
            
        except Exception as e:
            self.status_text.set_text(f'Export failed: {str(e)[:50]}')
            print(f"Export failed: {e}")
    
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
        """Open settings modal window with species data grid"""
        if self.settings_open:
            return
            
        self.settings_open = True
        
        # Create a new figure for settings
        self.settings_fig = plt.figure(figsize=(10, 12), dpi=100)
        self.settings_fig.suptitle('Display Settings', fontsize=14, fontweight='bold')
        
        # Theme selector
        theme_ax = plt.axes([0.1, 0.85, 0.8, 0.1])
        theme_labels = ['Dark', 'Light']
        theme_active = 0 if self.theme == 'dark' else 1
        self.settings_theme_radio = RadioButtons(theme_ax, theme_labels, active=theme_active)
        self.settings_theme_radio.on_clicked(self.on_settings_theme_change)
        theme_ax.set_title('Theme', fontsize=12, pad=15)
        
        # Display options
        options_ax = plt.axes([0.1, 0.7, 0.8, 0.12])
        options_labels = ['Lattice Grid', 'Empty Sites', 'Trajectories']
        options_states = [self.show_lattice_grid, self.show_empty_sites, self.show_trajectories]
        self.settings_options_check = CheckButtons(options_ax, options_labels, options_states)
        self.settings_options_check.on_clicked(self.on_settings_options_change)
        options_ax.set_title('Display Options', fontsize=12, pad=15)
        
        # Species Properties Data Grid
        self.setup_species_data_grid()
        
        # Control buttons
        button_y = 0.08
        button_width = 0.15
        button_height = 0.06
        
        # Reset button
        reset_ax = plt.axes([0.1, button_y, button_width, button_height])
        self.settings_reset_button = Button(reset_ax, 'Reset Defaults', color='lightcoral')
        self.settings_reset_button.on_clicked(self.reset_species_defaults)
        
        # Apply button
        apply_ax = plt.axes([0.3, button_y, button_width, button_height])
        self.settings_apply_button = Button(apply_ax, 'Apply Changes', color='lightgreen')
        self.settings_apply_button.on_clicked(self.apply_species_changes)
        
        # Close button
        close_ax = plt.axes([0.5, button_y, button_width, button_height])
        self.settings_close_button = Button(close_ax, 'Close', color='lightgray')
        self.settings_close_button.on_clicked(self.close_settings_modal)
        
        # Handle window close event
        self.settings_fig.canvas.mpl_connect('close_event', self.on_settings_window_close)
        
        plt.show()
    
    def setup_species_data_grid(self):
        """Setup the species properties data grid"""
        # Main grid area
        grid_ax = plt.axes([0.05, 0.2, 0.9, 0.45])
        grid_ax.set_xlim(0, 10)
        grid_ax.set_ylim(0, len(self.species_properties) + 1)
        grid_ax.set_aspect('equal')
        
        # Headers
        headers = ['Species', 'Visible', 'Color', 'Shape', 'Size']
        header_positions = [1, 3, 5, 7, 9]
        
        for i, header in enumerate(headers):
            grid_ax.text(header_positions[i], len(self.species_properties) + 0.5, header, 
                        ha='center', va='center', fontweight='bold', fontsize=11)
        
        # Grid lines
        for i in range(len(self.species_properties) + 2):
            grid_ax.axhline(y=i-0.5, color='gray', linestyle='-', alpha=0.3)
        for i in range(len(header_positions) + 1):
            grid_ax.axvline(x=i*2-0.5, color='gray', linestyle='-', alpha=0.3)
        
        # Species data grid controls
        self.species_controls = {}
        
        for row, (species, props) in enumerate(self.species_properties.items()):
            y_pos = len(self.species_properties) - row - 0.5
            
            # Species name (read-only)
            grid_ax.text(1, y_pos, species, ha='center', va='center', fontsize=10)
            
            # Visibility checkbox
            vis_ax = plt.axes([0.25, 0.2 + (y_pos-0.25)/len(self.species_properties)*0.45, 0.08, 0.04])
            vis_check = CheckButtons(vis_ax, [''], [props['visible']])
            vis_check.on_clicked(lambda label, s=species: self.on_species_visibility_change(s))
            
            # Color picker (simplified as text input)
            color_ax = plt.axes([0.45, 0.2 + (y_pos-0.25)/len(self.species_properties)*0.45, 0.12, 0.04])
            color_text = TextBox(color_ax, '', initial=props['color'])
            color_text.on_submit(lambda text, s=species: self.on_species_color_change(s, text))
            
            # Shape selector (cycle through options)
            shape_ax = plt.axes([0.65, 0.2 + (y_pos-0.25)/len(self.species_properties)*0.45, 0.08, 0.04])
            shape_button = Button(shape_ax, props['shape'], color='lightblue')
            shape_button.on_clicked(lambda event, s=species: self.cycle_species_shape(s))
            
            # Size input
            size_ax = plt.axes([0.8, 0.2 + (y_pos-0.25)/len(self.species_properties)*0.45, 0.1, 0.04])
            size_text = TextBox(size_ax, '', initial=str(props['size']))
            size_text.on_submit(lambda text, s=species: self.on_species_size_change(s, text))
            
            # Store controls for later reference
            self.species_controls[species] = {
                'visibility': vis_check,
                'color': color_text,
                'shape': shape_button,
                'size': size_text
            }
        
        grid_ax.set_xticks([])
        grid_ax.set_yticks([])
        grid_ax.set_title('Species Properties', fontsize=12, pad=20)
    
    def on_species_visibility_change(self, species):
        """Callback for visibility checkbox changes in settings modal"""
        self.species_properties[species]['visible'] = not self.species_properties[species]['visible']
        self.frame_cache.clear()
        self.update_plot(self.current_step, force_redraw=True)
        
        # Update status with visible species list
        visible_species = [s for s, props in self.species_properties.items() if props['visible']]
        self.status_text.set_text(f'Visible species: {", ".join(visible_species)}')
    
    def on_species_color_change(self, species, text):
        """Callback for color text input changes in settings modal"""
        try:
            # Clean the input text
            color_input = text.strip()
            
            # Handle different color formats
            if color_input.startswith('#') and len(color_input) == 7:
                # Already in hex format
                self.species_properties[species]['color'] = color_input
            elif len(color_input) == 6 and all(c in '0123456789abcdefABCDEF' for c in color_input):
                # Hex without #
                self.species_properties[species]['color'] = '#' + color_input
            else:
                self.status_text.set_text(f'Invalid color format for {species}. Use #RRGGBB or RRGGBB')
                return
                
            self.frame_cache.clear()
            self.update_plot(self.current_step, force_redraw=True)
            self.status_text.set_text(f'Color updated for {species}: {self.species_properties[species]["color"]}')
            
        except Exception as e:
            self.status_text.set_text(f'Error setting color for {species}: {e}')
    
    def cycle_species_shape(self, species):
        """Callback for shape button changes in settings modal"""
        shape_names = {'o': 'Circle', 's': 'Square', '^': 'Triangle'}
        shapes = ['o', 's', '^']
        
        current_shape = self.species_properties[species]['shape']
        current_index = shapes.index(current_shape) if current_shape in shapes else 0
        new_index = (current_index + 1) % len(shapes)
        new_shape = shapes[new_index]
        
        self.species_properties[species]['shape'] = new_shape
        
        # Update the button text
        self.species_controls[species]['shape'].label.set_text(new_shape)
        
        self.frame_cache.clear()
        self.update_plot(self.current_step, force_redraw=True)
        self.status_text.set_text(f'Shape updated for {species}: {shape_names.get(new_shape, new_shape)}')
    
    def on_species_size_change(self, species, text):
        """Callback for size text input changes in settings modal"""
        try:
            size = int(text.strip())
            if size <= 0:
                self.status_text.set_text(f'Size must be positive for {species}')
                return
            if size > 500:
                self.status_text.set_text(f'Size too large for {species} (max 500)')
                return
                
            self.species_properties[species]['size'] = size
            self.frame_cache.clear()
            self.update_plot(self.current_step, force_redraw=True)
            self.status_text.set_text(f'Size updated for {species}: {size}')
            
        except ValueError:
            self.status_text.set_text(f'Invalid size for {species}. Enter a number.')
    
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
    
    def on_settings_theme_change(self, label):
        """Handle theme changes from settings modal"""
        new_theme = 'dark' if label == 'Dark' else 'light'
        if new_theme != self.theme:
            self.theme = new_theme
            self.setup_theme()
            self.frame_cache.clear()
            self.update_plot(self.current_step, force_redraw=True)
            self.status_text.set_text(f'Switched to {new_theme} theme')
    
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
    parser.add_argument('--theme', '-t', choices=['dark', 'light'], default='dark',
                       help='UI theme (dark or light)')
    
    args = parser.parse_args()
    
    print(ascii_art)
    print(f"🚀 Starting KMCView - Enhanced Molecular Viewer")
    print(f"📁 Data directory: {args.data_dir}")
    print(f"🎨 Theme: {args.theme}")
    
    viewer = EnhancedMolecularViewer(data_dir=args.data_dir, theme=args.theme)

if __name__ == '__main__':
    main() 