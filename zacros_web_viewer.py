#!/usr/bin/env python3
"""
🌐 Zacros Web Molecular Evolution Viewer
=======================================

Complete web-based recreation of the matplotlib GUI viewer with full feature parity:
• Dynamic lattice parsing from lattice_input.dat
• Real-time molecular evolution visualization
• Working animation controls (play/pause/reset)
• Proper 2D lattice mesh display
• Configurable data directories
• Interactive controls and navigation
"""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
from dash import Dash, dcc, html, Input, Output, State, callback_context
import os
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

class ZacrosWebViewer:
    """Complete web-based molecular evolution viewer with GUI feature parity"""
    
    def __init__(self, data_dir='input-output', port=8050):
        """Initialize the web viewer with configurable data directory"""
        self.data_dir = data_dir
        self.port = port
        
        print(f"🌐 Starting Zacros Web Viewer...")
        print(f"📁 Data directory: {data_dir}")
        print(f"📈 Loading molecular evolution data...")
        
        # Initialize data containers
        self.evolution_df = pd.DataFrame()
        self.coordinates = {}
        self.lattice_bounds = {}
        self.lattice_info = {}
        
        # Load all data
        self.load_evolution_data()
        self.load_lattice_structure()
        self.generate_lattice_sites()
        
        # Set up Dash app
        self.app = Dash(__name__)
        self.setup_layout()
        self.setup_callbacks()
        
        # Print status
        print(f"✅ Loaded {len(self.evolution_df)} time points")
        print(f"✅ Generated {len(self.coordinates)} lattice sites")
        print(f"   Lattice dimensions: {self.lattice_bounds['x_max']:.2f} × {self.lattice_bounds['y_max']:.2f} Å")
        
        # Debug first few data points
        if len(self.evolution_df) > 0:
            print(f"   First data point: {self.evolution_df.iloc[0].to_dict()}")
            positions = self.generate_molecular_positions(0)
            print(f"   Generated positions for step 0: {positions}")
        
        # Start server
        print("🚀 Starting web server...")
        print(f"📱 Open your browser to http://127.0.0.1:{port}")
        print("🎯 Web viewer has complete GUI feature parity!")
        
    def load_evolution_data(self):
        """Load species evolution data from specnum_output.txt"""
        data_file = os.path.join(self.data_dir, 'specnum_output.txt')
        
        if not os.path.exists(data_file):
            print(f"❌ Data file not found: {data_file}")
            self.evolution_df = pd.DataFrame({'time': [0.0], 'H*': [0.0], 'GeH2*': [0.0], 'GeH3*': [0.0]})
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
                        if len(parts) >= 8:  # Need at least 8 columns
                            # Test if first part is numeric
                            float(parts[0])
                            step = int(parts[0])
                            nevents = int(parts[1])
                            time = float(parts[2])
                            
                            # Extract species counts (columns 5, 6, 7)
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
            print(f"❌ Error loading evolution data: {e}")
            self.evolution_df = pd.DataFrame({'time': [0.0], 'H*': [0.0], 'GeH2*': [0.0], 'GeH3*': [0.0]})
    
    def load_lattice_structure(self):
        """Load and parse lattice structure from lattice_input.dat"""
        lattice_file = os.path.join(self.data_dir, 'lattice_input.dat')
        
        # Parse lattice structure
        parser = LatticeParser(lattice_file)
        self.lattice_info = parser.get_lattice_info()
        
        # Print lattice information
        cell_vectors = self.lattice_info['cell_vectors']
        repeat_cell = self.lattice_info['repeat_cell']
        site_type_names = self.lattice_info['site_type_names']
        site_coordinates = self.lattice_info['site_coordinates']
        
        cell_vector_x = cell_vectors[0][0]
        cell_vector_y = cell_vectors[1][1]
        repeat_x, repeat_y = repeat_cell
        
        print(f"   Cell vectors: {cell_vector_x:.6f} × {cell_vector_y:.6f} Å")
        print(f"   Repeat pattern: {repeat_x} × {repeat_y}")
        print(f"   Site types: {site_type_names}")
        print(f"   Sites per unit cell: {len(site_coordinates)}")
    
    def generate_lattice_sites(self):
        """Generate all lattice sites from structure parameters"""
        # Extract parameters
        cell_vectors = self.lattice_info['cell_vectors']
        repeat_cell = self.lattice_info['repeat_cell']
        site_type_names = self.lattice_info['site_type_names']
        site_coordinates = self.lattice_info['site_coordinates']
        
        # Calculate cell dimensions
        cell_vector_x = cell_vectors[0][0]
        cell_vector_y = cell_vectors[1][1]
        repeat_x, repeat_y = repeat_cell
        
        # Generate all lattice sites
        site_id = 1
        self.coordinates = {}
        
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
            'cell_size_y': cell_vector_y,
            'repeat_x': repeat_x,
            'repeat_y': repeat_y
        }
    
    def generate_molecular_positions(self, step):
        """Generate molecular positions for a given time step (same as GUI version)"""
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
        
        # Generate consistent positions (same seed as GUI version)
        np.random.seed(42 + step)
        occupied_sites = np.random.choice(all_sites, 
                                        size=min(total_molecules, len(all_sites)), 
                                        replace=False)
        
        positions = {}
        idx = 0
        
        # Assign positions to species (same order as GUI)
        for species, count in [('H*', h_count), ('GeH2*', geh2_count), ('GeH3*', geh3_count)]:
            positions[species] = []
            for _ in range(count):
                if idx < len(occupied_sites):
                    site_id = occupied_sites[idx]
                    positions[species].append(self.coordinates[site_id])
                    idx += 1
        
        return positions
    
    def setup_layout(self):
        """Set up the Dash app layout with complete GUI feature parity"""
        max_steps = len(self.evolution_df) - 1 if len(self.evolution_df) > 0 else 0
        
        self.app.layout = html.Div([
            # Title
            html.H1("🔬 Zacros Molecular Evolution Viewer", 
                   style={'textAlign': 'center', 'marginBottom': '20px', 'color': '#2c3e50'}),
            
            # Info display (matches GUI info text)
            html.Div(id="info-display", style={
                'textAlign': 'center',
                'padding': '15px',
                'backgroundColor': '#e3f2fd',
                'borderRadius': '10px',
                'margin': '0 20px 20px 20px',
                'fontSize': '16px',
                'fontWeight': 'bold',
                'border': '2px solid #1976d2'
            }),
            
            # Control panel (matches GUI controls)
            html.Div([
                # Button controls row
                html.Div([
                    html.Button("Play", id="play-button", n_clicks=0,
                              style={
                                  'margin': '5px', 
                                  'padding': '12px 24px', 
                                  'fontSize': '16px',
                                  'backgroundColor': '#4caf50',
                                  'color': 'white',
                                  'border': 'none',
                                  'borderRadius': '5px',
                                  'cursor': 'pointer'
                              }),
                    html.Button("Reset", id="reset-button", n_clicks=0,
                              style={
                                  'margin': '5px', 
                                  'padding': '12px 24px', 
                                  'fontSize': '16px',
                                  'backgroundColor': '#f44336',
                                  'color': 'white',
                                  'border': 'none',
                                  'borderRadius': '5px',
                                  'cursor': 'pointer'
                              }),
                    html.Div([
                        html.Span("Step: ", style={'margin': '5px', 'fontSize': '16px', 'fontWeight': 'bold'}),
                        dcc.Input(id="step-input", type="number", value=0, min=0, 
                                 max=max_steps, 
                                 style={
                                     'margin': '5px', 
                                     'padding': '8px', 
                                     'fontSize': '16px',
                                     'width': '100px',
                                     'borderRadius': '3px',
                                     'border': '1px solid #ccc'
                                 }),
                    ], style={'display': 'inline-block', 'verticalAlign': 'middle'})
                ], style={'textAlign': 'center', 'marginBottom': '20px'}),
                
                # Time slider (matches GUI slider)
                html.Div([
                    html.Label("Time Step Navigation:", style={'fontSize': '16px', 'fontWeight': 'bold', 'marginBottom': '10px'}),
                    dcc.Slider(
                        id='time-slider',
                        min=0,
                        max=max_steps,
                        value=0,
                        marks={i: f"{i}" for i in range(0, max_steps + 1, max(1, max_steps//10))},
                        step=1,
                        tooltip={"placement": "bottom", "always_visible": True},
                        updatemode='drag'
                    )
                ], style={'margin': '20px 40px'}),
            ]),
            
            # Main plot area
            dcc.Graph(id="molecular-plot", style={'height': '70vh', 'margin': '20px'}),
            
            # Animation interval component (hidden)
            dcc.Interval(id="animation-interval", interval=500, n_intervals=0, disabled=True),
            
            # Hidden storage for animation state
            html.Div(id="animation-state", style={'display': 'none'}, children='stopped'),
            html.Div(id="current-step-store", style={'display': 'none'}, children='0')
        ])
    
    def setup_callbacks(self):
        """Set up all Dash callbacks with proper synchronization to avoid double-click issues"""
        
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
            """Single callback to handle all UI updates and avoid synchronization issues"""
            
            ctx = callback_context
            
            # Default button styles
            play_style = {
                'margin': '5px', 'padding': '12px 24px', 'fontSize': '16px',
                'backgroundColor': '#4caf50', 'color': 'white', 'border': 'none',
                'borderRadius': '5px', 'cursor': 'pointer'
            }
            pause_style = {
                'margin': '5px', 'padding': '12px 24px', 'fontSize': '16px',
                'backgroundColor': '#ff9800', 'color': 'white', 'border': 'none',
                'borderRadius': '5px', 'cursor': 'pointer'
            }
            
            # Determine what triggered the callback and what step to display
            if not ctx.triggered:
                # Initial load
                step = 0
                is_playing = False
                button_text = 'Play'
                button_style = play_style
                
            else:
                trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
                
                if trigger_id == 'reset-button':
                    # Reset to beginning
                    step = 0
                    is_playing = False
                    button_text = 'Play'
                    button_style = play_style
                    print("🔄 Animation reset")
                    
                elif trigger_id == 'play-button':
                    # Toggle play/pause
                    current_step = int(current_step_store) if current_step_store else int(slider_value or 0)
                    
                    if interval_disabled:  # Currently stopped, start playing
                        step = current_step
                        is_playing = True
                        button_text = 'Pause'
                        button_style = pause_style
                        print("▶️ Animation playing")
                    else:  # Currently playing, pause
                        step = current_step
                        is_playing = False
                        button_text = 'Play'
                        button_style = play_style
                        print("⏸️ Animation paused")
                        
                elif trigger_id == 'animation-interval':
                    # Animation tick - advance to next step
                    if animation_state == 'playing':
                        current_step = int(current_step_store) if current_step_store else 0
                        step = (current_step + 1) % len(self.evolution_df)
                        is_playing = True
                        button_text = 'Pause'
                        button_style = pause_style
                    else:
                        # If not playing, keep current step
                        step = int(current_step_store) if current_step_store else 0
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
            
            # Return all outputs in a single callback to avoid synchronization issues
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
        # All callback functionality is now consolidated into the single update_all_components callback above
        # This eliminates synchronization issues and ensures single-click responsiveness
    
    def add_lattice_grid(self, fig):
        """Add lattice grid lines (same as GUI version)"""
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
        """Add empty lattice sites (same as GUI version)"""
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
        """Add molecules with same colors/sizes as GUI version"""
        # Same colors and sizes as GUI (maintaining exact size ratios)
        colors = {'H*': 'red', 'GeH2*': 'blue', 'GeH3*': 'green'}
        sizes = {'H*': 8, 'GeH2*': 16, 'GeH3*': 24}  # Ratios: 2.0, 3.0 (same as GUI)
        
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
        """Configure plot layout (same as GUI version)"""
        # Set axis limits with margins (same as GUI)
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
        """Start the web server"""
        self.app.run(debug=False, port=self.port, host='127.0.0.1')

def main():
    """Main entry point for command line usage"""
    import sys
    
    data_dir = 'input-output'
    port = 8050
    
    if len(sys.argv) > 1:
        data_dir = sys.argv[1]
    if len(sys.argv) > 2:
        try:
            port = int(sys.argv[2])
        except ValueError:
            port = 8050
    
    viewer = ZacrosWebViewer(data_dir=data_dir, port=port)
    viewer.run()

if __name__ == '__main__':
    main() 