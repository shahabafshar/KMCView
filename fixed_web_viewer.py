#!/usr/bin/env python3
"""
🌐 Enhanced Web-Based Molecular Evolution Viewer
==============================================

Features:
• Dynamic lattice structure parsing from lattice_input.dat
• Works with any set of Zacros input files
• Configurable file paths and parameters
• Proper 2D lattice mesh visualization
• Working play/pause animation controls
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

def load_positions(data_dir='input-output'):
    """Generate lattice positions from lattice_input.dat structure"""
    lattice_file = os.path.join(data_dir, 'lattice_input.dat')
    
    # Parse lattice structure
    parser = LatticeParser(lattice_file)
    lattice_info = parser.get_lattice_info()
    
    # Extract parameters
    cell_vectors = lattice_info['cell_vectors']
    repeat_cell = lattice_info['repeat_cell']
    site_type_names = lattice_info['site_type_names']
    site_coordinates = lattice_info['site_coordinates']
    
    # Calculate cell dimensions
    cell_vector_x = cell_vectors[0][0]  # x-component of first vector
    cell_vector_y = cell_vectors[1][1]  # y-component of second vector
    repeat_x, repeat_y = repeat_cell
    
    print(f"   Cell vectors: {cell_vector_x:.6f} × {cell_vector_y:.6f} Å")
    print(f"   Repeat pattern: {repeat_x} × {repeat_y}")
    print(f"   Site types: {site_type_names}")
    
    # Generate all lattice sites
    coordinates = {}
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
                
                coordinates[site_id] = {
                    'x': abs_x,
                    'y': abs_y,
                    'site_type': site_type,
                    'unit_cell': (i, j),
                    'local_site': k
                }
                site_id += 1
    
    # Calculate lattice boundaries
    lattice_bounds = {
        'x_min': 0.0,
        'x_max': repeat_x * cell_vector_x,
        'y_min': 0.0,
        'y_max': repeat_y * cell_vector_y,
        'cell_size_x': cell_vector_x,
        'cell_size_y': cell_vector_y
    }
    
    print(f"✅ Generated {len(coordinates)} lattice sites")
    print(f"   Lattice dimensions: {lattice_bounds['x_max']:.2f} × {lattice_bounds['y_max']:.2f} Å")
    print(f"   Sites per unit cell: {len(site_coordinates)}")
    
    return coordinates, lattice_bounds

class FixedWebViewer:
    def __init__(self, data_dir='input-output', port=8050):
        """Initialize the web viewer with configurable data directory"""
        self.data_dir = data_dir
        self.port = port
        
        print(f"🌐 Starting Enhanced Web Viewer with data from: {data_dir}")
        print("📈 Loading molecular evolution data...")
        
        # Load data
        self.evolution_df = self.load_data()
        self.coordinates, self.lattice_bounds = load_positions(data_dir)
        
        # Set up Dash app
        self.app = Dash(__name__)
        self.setup_layout()
        self.setup_callbacks()
        
        print(f"✅ Loaded {len(self.evolution_df)} time points")
        print(f"✅ Generated {len(self.coordinates)} lattice sites")
        
        # Start server
        print("🚀 Starting web server...")
        print(f"📱 Open your browser to http://127.0.0.1:{port}")
        print("🎯 Play button should work now!")
        self.app.run(debug=False, port=port, host='127.0.0.1')
    
    def load_data(self):
        """Load species evolution data from specnum_output.txt"""
        data_file = os.path.join(self.data_dir, 'specnum_output.txt')
        
        if not os.path.exists(data_file):
            print(f"❌ Data file not found: {data_file}")
            return pd.DataFrame({'time': [0.0], 'H*': [0.0], 'GeH2*': [0.0], 'GeH3*': [0.0]})
        
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
            
            df = pd.DataFrame(data)
            print(f"✅ Loaded {len(df)} time points from {data_file}")
            return df
            
        except Exception as e:
            print(f"❌ Error loading data: {e}")
            return pd.DataFrame({'time': [0.0], 'H*': [0.0], 'GeH2*': [0.0], 'GeH3*': [0.0]})

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
    
    def setup_layout(self):
        """Set up the Dash app layout"""
        self.app.layout = html.Div([
            html.H1("🎬 Enhanced Interactive Molecular Evolution Viewer", 
                   style={'textAlign': 'center', 'marginBottom': '20px'}),
            
            # Control panel
            html.Div([
                html.Div([
                    html.Button("Play", id="play-button", n_clicks=0,
                              style={'margin': '5px', 'padding': '10px 20px', 'fontSize': '16px'}),
                    html.Button("Reset", id="reset-button", n_clicks=0,
                              style={'margin': '5px', 'padding': '10px 20px', 'fontSize': '16px'}),
                    dcc.Input(id="step-input", type="number", value=0, min=0, 
                             max=len(self.evolution_df)-1, 
                             style={'margin': '5px', 'padding': '10px', 'fontSize': '16px'}),
                    html.Span("Step", style={'margin': '5px', 'fontSize': '16px'})
                ], style={'textAlign': 'center', 'marginBottom': '20px'}),
                
                # Slider
                html.Div([
                    dcc.Slider(
                        id='time-slider',
                        min=0,
                        max=len(self.evolution_df)-1,
                        value=0,
                        marks={i: str(i) for i in range(0, len(self.evolution_df), 50)},
                        step=1,
                        tooltip={"placement": "bottom", "always_visible": True}
                    )
                ], style={'margin': '20px'}),
                
                # Info panel
                html.Div(id="info-panel", style={
                    'textAlign': 'center',
                    'padding': '15px',
                    'backgroundColor': '#f0f0f0',
                    'borderRadius': '5px',
                    'margin': '20px',
                    'fontSize': '16px'
                })
            ]),
            
            # Plot
            dcc.Graph(id="molecular-plot", style={'height': '70vh'}),
            
            # Interval for animation - THIS IS THE KEY COMPONENT
            dcc.Interval(id="interval-component", interval=300, n_intervals=0, disabled=True),
            
            # Hidden div to store current step
            html.Div(id="current-step-store", style={'display': 'none'}, children='0')
        ])
    
    def setup_callbacks(self):
        """Set up all Dash callbacks"""
        
        @self.app.callback(
            [Output("molecular-plot", "figure"),
             Output("info-panel", "children")],
            [Input("time-slider", "value"),
             Input("step-input", "value")]
        )
        def update_plot(slider_value, input_value):
            """Update the plot based on slider or input"""
            # Use input value if it triggered the callback, otherwise use slider
            ctx = callback_context
            if ctx.triggered and ctx.triggered[0]['prop_id'] == 'step-input.value':
                step = int(input_value) if input_value is not None else 0
            else:
                step = int(slider_value) if slider_value is not None else 0
            
            # Ensure step is within bounds
            step = max(0, min(step, len(self.evolution_df) - 1))
            
            # Get current data
            current_data = self.evolution_df.iloc[step]
            current_time = current_data['time']
            h_count = int(current_data['H*'])
            geh2_count = int(current_data['GeH2*'])
            geh3_count = int(current_data['GeH3*'])
            
            # Generate positions for this step
            positions = self.generate_positions_for_step(step)
            
            # Create figure
            fig = go.Figure()
            
            # Add lattice grid
            cell_size_x = self.lattice_bounds['cell_size_x']
            cell_size_y = self.lattice_bounds['cell_size_y']
            repeat_x = int(self.lattice_bounds['x_max'] / cell_size_x)
            repeat_y = int(self.lattice_bounds['y_max'] / cell_size_y)
            
            for i in range(repeat_x + 1):
                x_line = i * cell_size_x
                if x_line <= self.lattice_bounds['x_max']:
                    fig.add_shape(
                        type="line",
                        x0=x_line, y0=self.lattice_bounds['y_min'],
                        x1=x_line, y1=self.lattice_bounds['y_max'],
                        line=dict(color="lightgray", width=0.5)
                    )
            
            for j in range(repeat_y + 1):
                y_line = j * cell_size_y
                if y_line <= self.lattice_bounds['y_max']:
                    fig.add_shape(
                        type="line",
                        x0=self.lattice_bounds['x_min'], y0=y_line,
                        x1=self.lattice_bounds['x_max'], y1=y_line,
                        line=dict(color="lightgray", width=0.5)
                    )
            
            # Add empty lattice sites
            all_x = [coord['x'] for coord in self.coordinates.values()]
            all_y = [coord['y'] for coord in self.coordinates.values()]
            fig.add_trace(go.Scatter(
                x=all_x, y=all_y,
                mode='markers',
                marker=dict(color='lightgray', size=2, opacity=0.3),
                name='Empty sites',
                hoverinfo='skip'
            ))
            
            # Add molecules
            colors = {'H*': 'red', 'GeH2*': 'blue', 'GeH3*': 'green'}
            sizes = {'H*': 8, 'GeH2*': 12, 'GeH3*': 16}
            
            for species, coords in positions.items():
                if coords:
                    x_coords = [c['x'] for c in coords]
                    y_coords = [c['y'] for c in coords]
                    fig.add_trace(go.Scatter(
                        x=x_coords, y=y_coords,
                        mode='markers',
                        marker=dict(color=colors[species], size=sizes[species]),
                        name=f'{species} ({len(coords)})',
                        hovertemplate=f'{species}<br>X: %{{x:.1f}} Å<br>Y: %{{y:.1f}} Å<extra></extra>'
                    ))
            
            # Update layout
            fig.update_layout(
                title=f"Molecular Evolution - Step {step+1}/{len(self.evolution_df)} (Time: {current_time:.1f})",
                xaxis_title="X Coordinate (Å)",
                yaxis_title="Y Coordinate (Å)",
                xaxis=dict(
                    range=[self.lattice_bounds['x_min'] - 5, self.lattice_bounds['x_max'] + 5],
                    showgrid=True,
                    gridwidth=1,
                    gridcolor='lightgray'
                ),
                yaxis=dict(
                    range=[self.lattice_bounds['y_min'] - 5, self.lattice_bounds['y_max'] + 5],
                    showgrid=True,
                    gridwidth=1,
                    gridcolor='lightgray',
                    scaleanchor="x",
                    scaleratio=1
                ),
                showlegend=True,
                hovermode='closest',
                margin=dict(l=50, r=50, t=80, b=50)
            )
            
            # Info panel content
            info_content = [
                html.Div([
                    html.Strong(f"Step {step+1}/{len(self.evolution_df)} | Time: {current_time:.1f} | "),
                    html.Span(f"H*: {h_count} | ", style={'color': 'red'}),
                    html.Span(f"GeH2*: {geh2_count} | ", style={'color': 'blue'}),
                    html.Span(f"GeH3*: {geh3_count}", style={'color': 'green'})
                ])
            ]
            
            return fig, info_content
        
        @self.app.callback(
            [Output("interval-component", "disabled"),
             Output("play-button", "children"),
             Output("current-step-store", "children")],
            [Input("play-button", "n_clicks"),
             Input("reset-button", "n_clicks")],
            [State("interval-component", "disabled"),
             State("time-slider", "value")]
        )
        def control_animation(play_clicks, reset_clicks, interval_disabled, current_slider_value):
            """Control play/pause and reset functionality"""
            ctx = callback_context
            if not ctx.triggered:
                return True, "Play", "0"
            
            button_id = ctx.triggered[0]['prop_id'].split('.')[0]
            
            if button_id == "reset-button":
                print("🔄 Reset clicked")
                return True, "Play", "0"
            elif button_id == "play-button":
                if interval_disabled:
                    print("▶️ Play started")
                    return False, "Pause", str(current_slider_value)
                else:
                    print("⏸️ Play paused")
                    return True, "Play", str(current_slider_value)
            
            return interval_disabled, "Play" if interval_disabled else "Pause", str(current_slider_value)
        
        @self.app.callback(
            Output("time-slider", "value"),
            [Input("interval-component", "n_intervals"),
             Input("reset-button", "n_clicks")],
            [State("current-step-store", "children"),
             State("interval-component", "disabled")]
        )
        def advance_animation(n_intervals, reset_clicks, current_step_str, interval_disabled):
            """Advance animation automatically"""
            ctx = callback_context
            if not ctx.triggered:
                return 0
            
            button_id = ctx.triggered[0]['prop_id'].split('.')[0]
            
            if button_id == "reset-button":
                return 0
            elif button_id == "interval-component" and not interval_disabled:
                current_step = int(current_step_str) if current_step_str else 0
                next_step = (current_step + 1) % len(self.evolution_df)
                return next_step
            
            return int(current_step_str) if current_step_str else 0
        
        @self.app.callback(
            Output("current-step-store", "children"),
            [Input("time-slider", "value")]
        )
        def update_current_step(slider_value):
            """Update stored current step"""
            return str(slider_value)

if __name__ == '__main__':
    # Usage: python fixed_web_viewer.py [data_directory] [port]
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
    
    viewer = FixedWebViewer(data_dir=data_dir, port=port) 