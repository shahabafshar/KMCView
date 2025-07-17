#!/usr/bin/env python3
"""
🌐 Zacros Colab-Compatible Molecular Evolution Viewer
====================================================

Google Colab version with:
• Web viewer using Colab proxy tunneling
• Static matplotlib viewer for notebook display
• Sample data generation for demonstration
• No localhost dependencies
"""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
from dash import Dash, dcc, html, Input, Output, State, callback_context
import matplotlib.pyplot as plt
import os
import warnings
import random
from google.colab import output
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
        self.site_coordinates = [(0.0, 0.0), (0.0, 0.25), (0.0, 0.5), (0.0, 0.75)]
    
    def get_lattice_info(self):
        """Return parsed lattice information"""
        return {
            'cell_vectors': self.cell_vectors,
            'repeat_cell': self.repeat_cell,
            'site_type_names': self.site_type_names,
            'site_coordinates': self.site_coordinates
        }

class SampleDataGenerator:
    """Generate sample Zacros data for Colab demonstration"""
    
    @staticmethod
    def generate_sample_evolution_data(n_steps=100):
        """Generate realistic sample evolution data"""
        data = []
        time_step = 0.1
        
        for step in range(n_steps):
            time = step * time_step
            
            # Generate realistic molecular evolution patterns
            h_star = max(0, int(50 + 20 * np.sin(time * 0.5) + random.gauss(0, 5)))
            geh2_star = max(0, int(30 + 15 * np.sin(time * 0.3 + 1) + random.gauss(0, 3)))
            geh3_star = max(0, int(20 + 10 * np.sin(time * 0.7 + 2) + random.gauss(0, 2)))
            
            data.append({
                'step': step,
                'nevents': random.randint(10, 50),
                'time': time,
                'H*': h_star,
                'GeH2*': geh2_star,
                'GeH3*': geh3_star
            })
        
        return pd.DataFrame(data)
    
    @staticmethod
    def create_sample_files():
        """Create sample input files for demonstration"""
        # Create input-output directory
        os.makedirs('input-output', exist_ok=True)
        
        # Generate sample evolution data
        df = SampleDataGenerator.generate_sample_evolution_data()
        
        # Save as specnum_output.txt
        with open('input-output/specnum_output.txt', 'w') as f:
            f.write("# Step Nevents Time Temperature H* GeH2* GeH3*\n")
            for _, row in df.iterrows():
                f.write(f"{row['step']} {row['nevents']} {row['time']:.6f} 300.0 {row['H*']} {row['GeH2*']} {row['GeH3*']}\n")
        
        # Create sample lattice file
        with open('input-output/lattice_input.dat', 'w') as f:
            f.write("""# Lattice structure for demonstration
cell_vectors
6.133780000000000 0.000000000000000
0.000000000000000 6.133780000000000
repeat_cell 20 20
site_type_names
top1 bridge1 top2 bridge2
site_coordinates
0.000000000000000 0.000000000000000
0.000000000000000 0.250000000000000
0.000000000000000 0.500000000000000
0.000000000000000 0.750000000000000
""")
        
        print("✅ Created sample data files in input-output/")
        return df

class ColabMatplotlibViewer:
    """Static matplotlib viewer for Colab notebooks"""
    
    def __init__(self, data_dir='input-output'):
        self.data_dir = data_dir
        self.evolution_df = pd.DataFrame()
        self.coordinates = {}
        self.lattice_bounds = {}
        self.lattice_info = {}
        
        # Load data
        self.load_evolution_data()
        self.load_lattice_structure()
        self.generate_lattice_sites()
    
    def load_evolution_data(self):
        """Load species evolution data"""
        data_file = os.path.join(self.data_dir, 'specnum_output.txt')
        
        if not os.path.exists(data_file):
            print(f"❌ Data file not found: {data_file}")
            self.evolution_df = SampleDataGenerator.generate_sample_evolution_data()
            return
        
        try:
            data = []
            with open(data_file, 'r') as f:
                lines = f.readlines()
            
            for line in lines:
                line = line.strip()
                if line and not line.startswith('#'):
                    try:
                        parts = line.split()
                        if len(parts) >= 8:
                            float(parts[0])
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
            print(f"✅ Loaded {len(self.evolution_df)} time points")
            
        except Exception as e:
            print(f"❌ Error loading data: {e}")
            self.evolution_df = SampleDataGenerator.generate_sample_evolution_data()
    
    def load_lattice_structure(self):
        """Load lattice structure from parser"""
        parser = LatticeParser(os.path.join(self.data_dir, 'lattice_input.dat'))
        self.lattice_info = parser.get_lattice_info()
        
        # Calculate lattice bounds
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
    
    def generate_lattice_sites(self):
        """Generate all lattice site coordinates"""
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
        
        print(f"✅ Generated {len(self.coordinates)} lattice sites")
    
    def generate_molecular_positions(self, step):
        """Generate molecular positions for a given step"""
        if step >= len(self.evolution_df):
            return {}
        
        current_data = self.evolution_df.iloc[step]
        h_count = int(current_data['H*'])
        geh2_count = int(current_data['GeH2*'])
        geh3_count = int(current_data['GeH3*'])
        
        positions = {'H*': [], 'GeH2*': [], 'GeH3*': []}
        
        # Get available sites
        available_sites = list(self.coordinates.values())
        random.shuffle(available_sites)
        
        # Assign molecules to sites
        site_index = 0
        
        for species, count in [('H*', h_count), ('GeH2*', geh2_count), ('GeH3*', geh3_count)]:
            for _ in range(count):
                if site_index < len(available_sites):
                    positions[species].append(available_sites[site_index])
                    site_index += 1
        
        return positions
    
    def plot_step(self, step):
        """Create a static plot for a specific step"""
        if step >= len(self.evolution_df):
            step = len(self.evolution_df) - 1
        
        current_data = self.evolution_df.iloc[step]
        current_time = current_data['time']
        h_count = int(current_data['H*'])
        geh2_count = int(current_data['GeH2*'])
        geh3_count = int(current_data['GeH3*'])
        
        positions = self.generate_molecular_positions(step)
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Plot lattice grid
        cell_size_x = self.lattice_bounds['cell_size_x']
        cell_size_y = self.lattice_bounds['cell_size_y']
        repeat_x = self.lattice_bounds['repeat_x']
        repeat_y = self.lattice_bounds['repeat_y']
        
        # Vertical lines
        for i in range(repeat_x + 1):
            x_line = i * cell_size_x
            if x_line <= self.lattice_bounds['x_max']:
                ax.axvline(x=x_line, color='lightgray', alpha=0.5, linewidth=0.5)
        
        # Horizontal lines
        for j in range(repeat_y + 1):
            y_line = j * cell_size_y
            if y_line <= self.lattice_bounds['y_max']:
                ax.axhline(y=y_line, color='lightgray', alpha=0.5, linewidth=0.5)
        
        # Plot empty sites
        all_x = [coord['x'] for coord in self.coordinates.values()]
        all_y = [coord['y'] for coord in self.coordinates.values()]
        ax.scatter(all_x, all_y, c='lightgray', s=1, alpha=0.4, label='Empty sites')
        
        # Plot molecules
        colors = {'H*': 'red', 'GeH2*': 'blue', 'GeH3*': 'green'}
        sizes = {'H*': 50, 'GeH2*': 100, 'GeH3*': 150}
        
        for species, coords in positions.items():
            if coords:
                x_coords = [c['x'] for c in coords]
                y_coords = [c['y'] for c in coords]
                ax.scatter(x_coords, y_coords, c=colors[species], s=sizes[species], 
                          alpha=0.7, label=f'{species} ({len(coords)})')
        
        # Configure plot
        ax.set_xlim(self.lattice_bounds['x_min'] - 1, self.lattice_bounds['x_max'] + 1)
        ax.set_ylim(self.lattice_bounds['y_min'] - 1, self.lattice_bounds['y_max'] + 1)
        ax.set_xlabel('X Coordinate (Å)')
        ax.set_ylabel('Y Coordinate (Å)')
        ax.set_title(f'Molecular Evolution - Step {step+1}/{len(self.evolution_df)} (Time: {current_time:.1f})')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal')
        
        plt.tight_layout()
        return fig
    
    def plot_evolution_timeline(self):
        """Plot the evolution timeline"""
        fig, ax = plt.subplots(figsize=(12, 6))
        
        ax.plot(self.evolution_df['time'], self.evolution_df['H*'], 'r-', label='H*', linewidth=2)
        ax.plot(self.evolution_df['time'], self.evolution_df['GeH2*'], 'b-', label='GeH2*', linewidth=2)
        ax.plot(self.evolution_df['time'], self.evolution_df['GeH3*'], 'g-', label='GeH3*', linewidth=2)
        
        ax.set_xlabel('Time')
        ax.set_ylabel('Number of Molecules')
        ax.set_title('Molecular Evolution Timeline')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        return fig

class ColabWebViewer:
    """Web viewer with Colab proxy support"""
    
    def __init__(self, data_dir='input-output', port=8050):
        self.data_dir = data_dir
        self.port = port
        
        # Load data
        self.evolution_df = pd.DataFrame()
        self.coordinates = {}
        self.lattice_bounds = {}
        self.lattice_info = {}
        
        self.load_evolution_data()
        self.load_lattice_structure()
        self.generate_lattice_sites()
        
        # Set up Dash app
        self.app = Dash(__name__)
        self.setup_layout()
        self.setup_callbacks()
    
    def load_evolution_data(self):
        """Load species evolution data"""
        data_file = os.path.join(self.data_dir, 'specnum_output.txt')
        
        if not os.path.exists(data_file):
            print(f"❌ Data file not found: {data_file}")
            self.evolution_df = SampleDataGenerator.generate_sample_evolution_data()
            return
        
        try:
            data = []
            with open(data_file, 'r') as f:
                lines = f.readlines()
            
            for line in lines:
                line = line.strip()
                if line and not line.startswith('#'):
                    try:
                        parts = line.split()
                        if len(parts) >= 8:
                            float(parts[0])
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
            print(f"✅ Loaded {len(self.evolution_df)} time points")
            
        except Exception as e:
            print(f"❌ Error loading data: {e}")
            self.evolution_df = SampleDataGenerator.generate_sample_evolution_data()
    
    def load_lattice_structure(self):
        """Load lattice structure from parser"""
        parser = LatticeParser(os.path.join(self.data_dir, 'lattice_input.dat'))
        self.lattice_info = parser.get_lattice_info()
        
        # Calculate lattice bounds
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
    
    def generate_lattice_sites(self):
        """Generate all lattice site coordinates"""
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
        
        print(f"✅ Generated {len(self.coordinates)} lattice sites")
    
    def generate_molecular_positions(self, step):
        """Generate molecular positions for a given step"""
        if step >= len(self.evolution_df):
            return {}
        
        current_data = self.evolution_df.iloc[step]
        h_count = int(current_data['H*'])
        geh2_count = int(current_data['GeH2*'])
        geh3_count = int(current_data['GeH3*'])
        
        positions = {'H*': [], 'GeH2*': [], 'GeH3*': []}
        
        # Get available sites
        available_sites = list(self.coordinates.values())
        random.shuffle(available_sites)
        
        # Assign molecules to sites
        site_index = 0
        
        for species, count in [('H*', h_count), ('GeH2*', geh2_count), ('GeH3*', geh3_count)]:
            for _ in range(count):
                if site_index < len(available_sites):
                    positions[species].append(available_sites[site_index])
                    site_index += 1
        
        return positions
    
    def setup_layout(self):
        """Set up the Dash layout"""
        self.app.layout = html.Div([
            html.H1("🌐 Zacros Molecular Evolution Viewer (Colab)", 
                   style={'textAlign': 'center', 'color': '#1976d2'}),
            
            # Control panel
            html.Div([
                html.Div([
                    html.Label("Time Step:", style={'fontWeight': 'bold'}),
                    dcc.Slider(
                        id='time-slider',
                        min=0,
                        max=max(0, len(self.evolution_df) - 1),
                        value=0,
                        marks={i: str(i) for i in range(0, len(self.evolution_df), max(1, len(self.evolution_df)//10))},
                        step=1
                    )
                ], style={'width': '60%', 'display': 'inline-block'}),
                
                html.Div([
                    html.Label("Step Input:", style={'fontWeight': 'bold'}),
                    dcc.Input(
                        id='step-input',
                        type='number',
                        value=0,
                        min=0,
                        max=max(0, len(self.evolution_df) - 1),
                        style={'width': '80px'}
                    )
                ], style={'width': '20%', 'display': 'inline-block', 'marginLeft': '20px'}),
                
                html.Div([
                    html.Button('Play', id='play-button', n_clicks=0, 
                               style={'backgroundColor': '#4CAF50', 'color': 'white', 'border': 'none', 'padding': '10px 20px'}),
                    html.Button('Reset', id='reset-button', n_clicks=0,
                               style={'backgroundColor': '#f44336', 'color': 'white', 'border': 'none', 'padding': '10px 20px', 'marginLeft': '10px'})
                ], style={'width': '20%', 'display': 'inline-block', 'textAlign': 'right'})
            ], style={'marginBottom': '20px', 'padding': '20px', 'backgroundColor': '#f5f5f5', 'borderRadius': '5px'}),
            
            # Animation interval
            dcc.Interval(
                id='animation-interval',
                interval=500,  # milliseconds
                n_intervals=0,
                disabled=True
            ),
            
            # Hidden stores
            dcc.Store(id='animation-state', data='stopped'),
            dcc.Store(id='current-step-store', data='0'),
            
            # Main plot
            dcc.Graph(id='molecular-plot', style={'height': '600px'}),
            
            # Info display
            html.Div(id='info-display', style={'textAlign': 'center', 'marginTop': '10px', 'fontSize': '16px'})
        ], style={'padding': '20px'})
    
    def setup_callbacks(self):
        """Set up Dash callbacks"""
        @self.app.callback(
            [Output("molecular-plot", "figure"),
             Output("info-display", "children"),
             Output("time-slider", "value"),
             Output("step-input", "value"),
             Output("animation-interval", "disabled"),
             Output("animation-state", "children"),
             Output("play-button", "children"),
             Output("play-button", "style"),
             Output("current-step-store", "children")],
            [Input("time-slider", "value"),
             Input("step-input", "value"),
             Input("play-button", "n_clicks"),
             Input("reset-button", "n_clicks"),
             Input("animation-interval", "n_intervals")],
            [State("animation-state", "children"),
             State("current-step-store", "children"),
             State("animation-interval", "disabled")]
        )
        def update_all_components(slider_value, input_value, play_clicks, reset_clicks, 
                                n_intervals, animation_state, current_step_store, interval_disabled):
            
            # Determine what triggered the callback
            ctx = callback_context
            trigger_id = ctx.triggered[0]['prop_id'].split('.')[0] if ctx.triggered else 'no-trigger'
            
            # Initialize styles
            play_style = {'backgroundColor': '#4CAF50', 'color': 'white', 'border': 'none', 'padding': '10px 20px'}
            pause_style = {'backgroundColor': '#FF9800', 'color': 'white', 'border': 'none', 'padding': '10px 20px'}
            
            # Handle different triggers
            if trigger_id == 'play-button':
                # Toggle play/pause
                if animation_state == 'playing':
                    is_playing = False
                    button_text = 'Play'
                    button_style = play_style
                else:
                    is_playing = True
                    button_text = 'Pause'
                    button_style = pause_style
                step = int(current_step_store or 0)
                
            elif trigger_id == 'reset-button':
                # Reset to beginning
                step = 0
                is_playing = False
                button_text = 'Play'
                button_style = play_style
                
            elif trigger_id == 'animation-interval':
                # Animation tick
                if animation_state == 'playing':
                    step = min(int(current_step_store or 0) + 1, len(self.evolution_df) - 1)
                    is_playing = step < len(self.evolution_df) - 1
                    button_text = 'Pause' if is_playing else 'Play'
                    button_style = pause_style if is_playing else play_style
                else:
                    step = int(current_step_store or 0)
                    is_playing = False
                    button_text = 'Play'
                    button_style = play_style
                    
            elif trigger_id == 'step-input':
                # Manual step input
                step = max(0, min(int(input_value or 0), len(self.evolution_df) - 1))
                is_playing = animation_state == 'playing'
                button_text = 'Pause' if is_playing else 'Play'
                button_style = pause_style if is_playing else play_style
                
            elif trigger_id == 'time-slider':
                # Slider change
                step = int(slider_value or 0)
                is_playing = animation_state == 'playing'
                button_text = 'Pause' if is_playing else 'Play'
                button_style = pause_style if is_playing else play_style
                
            else:
                # Default case
                step = int(slider_value or 0)
                is_playing = animation_state == 'playing'
                button_text = 'Pause' if is_playing else 'Play'
                button_style = pause_style if is_playing else play_style
            
            # Ensure step is within bounds
            step = max(0, min(step, len(self.evolution_df) - 1))
            
            # Get current simulation data
            current_data = self.evolution_df.iloc[step]
            current_time = current_data['time']
            h_count = int(current_data['H*'])
            geh2_count = int(current_data['GeH2*'])
            geh3_count = int(current_data['GeH3*'])
            
            # Generate molecular positions
            positions = self.generate_molecular_positions(step)
            
            # Create the figure
            fig = go.Figure()
            
            # Add lattice grid
            self.add_lattice_grid(fig)
            
            # Add empty lattice sites
            self.add_empty_sites(fig)
            
            # Add molecules
            self.add_molecules(fig, positions)
            
            # Configure plot layout
            self.configure_plot_layout(fig, step, current_time)
            
            # Create info display
            info_content = html.Div([
                html.Span(f"Time: {current_time:.1f} | Step: {step+1}/{len(self.evolution_df)} | ", 
                         style={'color': '#1976d2', 'fontWeight': 'bold'}),
                html.Span(f"H*: {h_count} | ", style={'color': 'red', 'fontWeight': 'bold'}),
                html.Span(f"GeH2*: {geh2_count} | ", style={'color': 'blue', 'fontWeight': 'bold'}),
                html.Span(f"GeH3*: {geh3_count}", style={'color': 'green', 'fontWeight': 'bold'})
            ])
            
            return (
                fig,                           # molecular-plot figure
                info_content,                  # info-display children
                step,                          # time-slider value
                step,                          # step-input value
                not is_playing,                # animation-interval disabled
                'playing' if is_playing else 'stopped',  # animation-state children
                button_text,                   # play-button children
                button_style,                  # play-button style
                str(step)                      # current-step-store children
            )
    
    def add_lattice_grid(self, fig):
        """Add lattice grid lines"""
        cell_size_x = self.lattice_bounds['cell_size_x']
        cell_size_y = self.lattice_bounds['cell_size_y']
        repeat_x = self.lattice_bounds['repeat_x']
        repeat_y = self.lattice_bounds['repeat_y']
        
        # Vertical lines
        for i in range(repeat_x + 1):
            x_line = i * cell_size_x
            if x_line <= self.lattice_bounds['x_max']:
                fig.add_shape(
                    type="line",
                    x0=x_line, y0=self.lattice_bounds['y_min'],
                    x1=x_line, y1=self.lattice_bounds['y_max'],
                    line=dict(color="lightgray", width=0.5)
                )
        
        # Horizontal lines
        for j in range(repeat_y + 1):
            y_line = j * cell_size_y
            if y_line <= self.lattice_bounds['y_max']:
                fig.add_shape(
                    type="line",
                    x0=self.lattice_bounds['x_min'], y0=y_line,
                    x1=self.lattice_bounds['x_max'], y1=y_line,
                    line=dict(color="lightgray", width=0.5)
                )
    
    def add_empty_sites(self, fig):
        """Add empty lattice sites"""
        all_x = [coord['x'] for coord in self.coordinates.values()]
        all_y = [coord['y'] for coord in self.coordinates.values()]
        
        if all_x and all_y:
            fig.add_trace(go.Scatter(
                x=all_x, y=all_y,
                mode='markers',
                marker=dict(color='lightgray', size=1, opacity=0.4),
                name='Empty sites',
                hoverinfo='skip',
                showlegend=False
            ))
    
    def add_molecules(self, fig, positions):
        """Add molecules"""
        colors = {'H*': 'red', 'GeH2*': 'blue', 'GeH3*': 'green'}
        sizes = {'H*': 8, 'GeH2*': 16, 'GeH3*': 24}
        
        for species, coords in positions.items():
            if coords:
                x_coords = [c['x'] for c in coords]
                y_coords = [c['y'] for c in coords]
                
                fig.add_trace(go.Scatter(
                    x=x_coords, y=y_coords,
                    mode='markers',
                    marker=dict(
                        color=colors[species], 
                        size=sizes[species],
                        opacity=0.7
                    ),
                    name=f'{species} ({len(coords)})',
                    hovertemplate=f'{species}<br>X: %{{x:.1f}} Å<br>Y: %{{y:.1f}} Å<extra></extra>'
                ))
    
    def configure_plot_layout(self, fig, step, current_time):
        """Configure plot layout"""
        x_margin = (self.lattice_bounds['x_max'] - self.lattice_bounds['x_min']) * 0.05
        y_margin = (self.lattice_bounds['y_max'] - self.lattice_bounds['y_min']) * 0.05
        
        fig.update_layout(
            title=f"Molecular Evolution - Step {step+1}/{len(self.evolution_df)} (Time: {current_time:.1f})",
            xaxis_title="X Coordinate (Å)",
            yaxis_title="Y Coordinate (Å)",
            xaxis=dict(
                range=[self.lattice_bounds['x_min'] - x_margin, self.lattice_bounds['x_max'] + x_margin],
                showgrid=True,
                gridwidth=1,
                gridcolor='lightgray'
            ),
            yaxis=dict(
                range=[self.lattice_bounds['y_min'] - y_margin, self.lattice_bounds['y_max'] + y_margin],
                showgrid=True,
                gridwidth=1,
                gridcolor='lightgray',
                scaleanchor="x",
                scaleratio=1
            ),
            showlegend=True,
            legend=dict(x=0.02, y=0.98),
            hovermode='closest',
            margin=dict(l=60, r=20, t=80, b=60),
            plot_bgcolor='white',
            paper_bgcolor='white'
        )
    
    def run(self):
        """Start the web server with Colab proxy"""
        print("🚀 Starting Colab-compatible web viewer...")
        print("📱 The viewer will be available via Colab's proxy tunnel")
        
        # Run with host='0.0.0.0' for Colab compatibility
        self.app.run(debug=False, port=self.port, host='0.0.0.0')

def create_sample_data():
    """Create sample data for demonstration"""
    return SampleDataGenerator.create_sample_files()

def run_matplotlib_viewer(data_dir='input-output'):
    """Run the matplotlib viewer in Colab"""
    viewer = ColabMatplotlibViewer(data_dir)
    
    # Show first step
    fig1 = viewer.plot_step(0)
    plt.show()
    
    # Show evolution timeline
    fig2 = viewer.plot_evolution_timeline()
    plt.show()
    
    return viewer

def run_web_viewer(data_dir='input-output', port=8050):
    """Run the web viewer with Colab proxy"""
    viewer = ColabWebViewer(data_dir, port)
    viewer.run()
    return viewer

# Main execution for Colab
if __name__ == '__main__':
    # Create sample data if needed
    create_sample_data()
    
    # Run matplotlib viewer
    print("🔬 Running Matplotlib Viewer...")
    viewer = run_matplotlib_viewer()
    
    # Run web viewer
    print("\n🌐 Running Web Viewer...")
    web_viewer = run_web_viewer() 