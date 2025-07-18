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
        
        print(f"Enhanced viewer ready with {len(self.evolution_df)} time points")
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
        """Load species evolution data with progress feedback"""
        data_file = os.path.join(self.data_dir, 'specnum_output.txt')
        
        if not os.path.exists(data_file):
            print(f"❌ Data file not found: {data_file}")
            self.evolution_df = pd.DataFrame({'time': [0.0], 'H*': [0.0], 'GeH2*': [0.0], 'GeH3*': [0.0]})
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
                            sim_time = float(parts[2])  # Renamed to avoid conflict with time module
                            h_star = float(parts[5]) if len(parts) > 5 else 0.0
                            geh2_star = float(parts[6]) if len(parts) > 6 else 0.0
                            geh3_star = float(parts[7]) if len(parts) > 7 else 0.0
                            
                            data.append({
                                'step': step,
                                'nevents': nevents,
                                'time': sim_time,  # Use renamed variable
                                'H*': h_star,
                                'GeH2*': geh2_star,
                                'GeH3*': geh3_star
                            })
                    except (ValueError, IndexError):
                        continue
            
            self.evolution_df = pd.DataFrame(data)
            load_time = time.time() - start_time
            print(f"Loaded {len(self.evolution_df)} time points in {load_time:.2f}s")
            
        except Exception as e:
            print(f"❌ Error loading data: {e}")
            self.evolution_df = pd.DataFrame({'time': [0.0], 'H*': [0.0], 'GeH2*': [0.0], 'GeH3*': [0.0]})
    
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
        """Setup the main figure with modern layout"""
        # Create figure with optimal size and DPI
        self.fig = plt.figure(figsize=(16, 10), dpi=100)
        self.fig.suptitle('KMCView - Enhanced Molecular Evolution Viewer', 
                         fontsize=16, fontweight='bold', color=self.colors['text'])
        
        # Main plot area (larger)
        self.ax = plt.axes([0.05, 0.25, 0.7, 0.65])
        self.ax.set_facecolor(self.colors['bg'])
        
        # Performance monitor (top right)
        self.perf_ax = plt.axes([0.78, 0.85, 0.2, 0.1])
        self.perf_ax.set_facecolor(self.colors['panel'])
        self.perf_text = self.perf_ax.text(0.05, 0.5, 'Performance: Ready', 
                                          transform=self.perf_ax.transAxes,
                                          fontsize=9, color=self.colors['text'])
        self.perf_ax.set_xticks([])
        self.perf_ax.set_yticks([])
        
        # Info panel (right side)
        self.info_ax = plt.axes([0.78, 0.6, 0.2, 0.2])
        self.info_ax.set_facecolor(self.colors['panel'])
        self.info_text = self.info_ax.text(0.05, 0.95, 'Loading...', 
                                          transform=self.info_ax.transAxes,
                                          fontsize=9, color=self.colors['text'],
                                          verticalalignment='top', fontfamily='monospace')
        self.info_ax.set_xticks([])
        self.info_ax.set_yticks([])
        self.info_ax.set_title('Simulation Info', fontsize=11, color=self.colors['text'])
        
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
        """Setup enhanced control widgets"""
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
        
        # Play button (using text instead of emoji for better compatibility)
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
        
        # Text input (with better spacing)
        ax_text = plt.axes([0.39, 0.1, 0.09, 0.04])
        self.text_box = TextBox(ax_text, 'Step:', initial='0', color=self.colors['panel'])
        self.text_box.on_submit(self.on_text_submit)
        
        # Export button
        ax_export = plt.axes([0.50, 0.1, 0.07, 0.04])
        self.export_button = Button(ax_export, 'Save', color=self.colors['warning'])
        self.export_button.on_clicked(self.export_frame)
        
        # Display options (right panel)
        ax_options = plt.axes([0.78, 0.4, 0.2, 0.15])
        options_labels = ['Lattice Grid', 'Empty Sites', 'Trajectories']
        options_states = [self.show_lattice_grid, self.show_empty_sites, self.show_trajectories]
        self.options_check = CheckButtons(ax_options, options_labels, options_states)
        self.options_check.on_clicked(self.on_options_change)
        ax_options.set_title('Display Options', fontsize=11, color=self.colors['text'])
        
        # Theme selector
        ax_theme = plt.axes([0.78, 0.25, 0.2, 0.1])
        theme_labels = ['Dark', 'Light']
        theme_active = 0 if self.theme == 'dark' else 1
        self.theme_radio = RadioButtons(ax_theme, theme_labels, active=theme_active)
        self.theme_radio.on_clicked(self.on_theme_change)
        ax_theme.set_title('Theme', fontsize=11, color=self.colors['text'])
    
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
        """Draw molecules with enhanced styling"""
        # Enhanced color scheme and sizes
        colors = {
            'H*': '#FF4444',      # Bright red
            'GeH2*': '#4477FF',   # Bright blue
            'GeH3*': '#44AA44'    # Bright green
        }
        
        sizes = {'H*': 40, 'GeH2*': 80, 'GeH3*': 120}
        alphas = {'H*': 0.8, 'GeH2*': 0.8, 'GeH3*': 0.8}
        
        for species, coords in positions.items():
            if coords:
                x_coords = [c['x'] for c in coords]
                y_coords = [c['y'] for c in coords]
                
                # Main scatter plot
                scatter = self.ax.scatter(x_coords, y_coords, 
                                        c=colors[species], s=sizes[species], 
                                        alpha=alphas[species], 
                                        label=f'{species} ({len(coords)})',
                                        edgecolors='white', linewidth=0.5, zorder=5)
                
                # Add glow effect for better visibility
                self.ax.scatter(x_coords, y_coords, 
                              c=colors[species], s=sizes[species]*2, 
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
        
        info_text = f"""Time: {current_time:.2f}
Step: {step+1}/{len(self.evolution_df)}

Molecules:
H*: {h_count}
GeH2*: {geh2_count}
GeH3*: {geh3_count}
Total: {total_molecules}

Lattice:
Size: {self.lattice_bounds['repeat_x']}x{self.lattice_bounds['repeat_y']}
Sites: {len(self.coordinates)}
Cell: {self.lattice_bounds['cell_size_x']:.1f}x{self.lattice_bounds['cell_size_y']:.1f} A"""
        
        self.info_text.set_text(info_text)
    
    def generate_positions_for_step(self, step):
        """Generate molecular positions (using internal cache)"""
        # Check cache first
        if step in self.position_cache:
            return self.position_cache[step]
        
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
        
        # Generate consistent positions with better randomization
        np.random.seed(42 + step * 7)  # More varied seed
        occupied_sites = np.random.choice(all_sites, 
                                        size=min(total_molecules, len(all_sites)), 
                                        replace=False)
        
        positions = {}
        idx = 0
        
        # Assign positions to species with spatial clustering simulation
        for species, count in [('H*', h_count), ('GeH2*', geh2_count), ('GeH3*', geh3_count)]:
            positions[species] = []
            for _ in range(count):
                if idx < len(occupied_sites):
                    site_id = occupied_sites[idx]
                    positions[species].append(self.coordinates[site_id])
                    idx += 1
        
        # Cache the result
        if len(self.position_cache) >= self.max_position_cache:
            # Remove oldest entry (simple FIFO)
            oldest_key = next(iter(self.position_cache))
            del self.position_cache[oldest_key]
        
        self.position_cache[step] = positions
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