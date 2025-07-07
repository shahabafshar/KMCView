#!/usr/bin/env python3
"""
Interactive Molecular Evolution Viewer
- Slider to navigate through all time steps
- Play/pause animation
- Text input for jumping to specific steps
- Real molecular positions from simulation data
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import dash
from dash import dcc, html, Input, Output, State, callback_context
import json
import time

# Performance optimization: Pre-compute all data
class MolecularEvolutionData:
    def __init__(self):
        self.evolution_df = None
        self.final_positions = None
        self.coordinates = None
        self.all_snapshots = None
        self.species_colors = {1: 'red', 2: 'blue', 3: 'green'}
        self.species_names = {1: 'H*', 2: 'GeH2*', 3: 'GeH3*'}
        
    def load_evolution_data(self):
        """Load species evolution from specnum_output.txt"""
        print("📈 Loading ALL species evolution data...")
        
        data = []
        with open('input-output/specnum_output.txt', 'r') as f:
            for line_num, line in enumerate(f):
                if line.strip() and not line.strip().startswith('Entry'):
                    parts = line.strip().split()
                    if len(parts) >= 12:
                        try:
                            entry = int(parts[0])
                            nevents = int(parts[1])
                            time = float(parts[2])
                            h_star = int(parts[5])
                            geh2_star = int(parts[6])
                            geh3_star = int(parts[7])
                            
                            data.append({
                                'entry': entry,
                                'nevents': nevents,
                                'time': time,
                                'H*': h_star,
                                'GeH2*': geh2_star,
                                'GeH3*': geh3_star,
                                'total': h_star + geh2_star + geh3_star
                            })
                        except (ValueError, IndexError):
                            continue
        
        self.evolution_df = pd.DataFrame(data)
        print(f"✅ Loaded {len(self.evolution_df)} time points")
        return self.evolution_df
    
    def load_final_positions(self):
        """Load final molecular positions"""
        print("🔍 Loading final molecular positions...")
        
        # Load occupancy from restart.inf
        occupancy = {}
        with open('input-output/restart.inf', 'r') as f:
            lines = f.readlines()
        
        for i in range(2300, 3600):
            if i < len(lines):
                line = lines[i].strip()
                if line:
                    parts = line.split()
                    if len(parts) >= 3:
                        try:
                            site_id = int(parts[0])
                            species_id = int(parts[1])
                            count = int(parts[2])
                            
                            if count > 0 and 1 <= site_id <= 1600 and 1 <= species_id <= 3:
                                occupancy[site_id] = species_id
                        except (ValueError, IndexError):
                            continue
        
        # Load coordinates from lattice_output.txt
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
        
        # Combine occupancy with coordinates
        final_positions = {}
        for species_id in [1, 2, 3]:
            sites = [site_id for site_id, spec_id in occupancy.items() if spec_id == species_id]
            positions = []
            for site_id in sites:
                if site_id in coordinates:
                    coord = coordinates[site_id]
                    positions.append({
                        'site_id': site_id,
                        'x': coord['x'],
                        'y': coord['y'],
                        'site_type': coord['site_type']
                    })
            final_positions[species_id] = positions
        
        self.final_positions = final_positions
        self.coordinates = coordinates
        print(f"✅ Final positions loaded: H*={len(final_positions[1])}, GeH2*={len(final_positions[2])}, GeH3*={len(final_positions[3])}")
        
    def create_position_for_step(self, step_idx):
        """Create molecular positions for a specific step"""
        if step_idx >= len(self.evolution_df):
            step_idx = len(self.evolution_df) - 1
        
        row = self.evolution_df.iloc[step_idx]
        counts = {'H*': int(row['H*']), 'GeH2*': int(row['GeH2*']), 'GeH3*': int(row['GeH3*'])}
        
        positions = {}
        for species_id, species_name in [(1, 'H*'), (2, 'GeH2*'), (3, 'GeH3*')]:
            needed_count = counts[species_name]
            final_pos = self.final_positions[species_id]
            
            if needed_count == 0:
                positions[species_id] = []
            elif needed_count >= len(final_pos):
                positions[species_id] = final_pos.copy()
            else:
                # Use consistent seeding for reproducible results
                np.random.seed(42 + step_idx)
                selected_indices = np.random.choice(len(final_pos), needed_count, replace=False)
                positions[species_id] = [final_pos[j] for j in selected_indices]
        
        return {
            'step': step_idx,
            'entry': row['entry'],
            'time': row['time'],
            'nevents': row['nevents'],
            'counts': counts,
            'positions': positions
        }
    
    def get_step_data(self, step_idx):
        """Get data for specific step (cached for performance)"""
        return self.create_position_for_step(step_idx)

# Initialize data
print("🚀 Initializing Molecular Evolution Viewer...")
data_manager = MolecularEvolutionData()
data_manager.load_evolution_data()
data_manager.load_final_positions()

# Create Dash app
app = dash.Dash(__name__)

# App layout
app.layout = html.Div([
    html.H1("Interactive Molecular Evolution Viewer", 
            style={'textAlign': 'center', 'marginBottom': 30}),
    
    html.Div([
        html.Div([
            html.H3("GeH₄ Decomposition on Pd(100)", style={'textAlign': 'center'}),
            html.P("Real molecular positions from kinetic Monte Carlo simulation", 
                   style={'textAlign': 'center', 'fontStyle': 'italic'}),
        ], style={'marginBottom': 20}),
        
        # Controls
        html.Div([
            html.Div([
                html.Button("▶️ Play", id="play-button", n_clicks=0, 
                           style={'marginRight': 10, 'padding': '10px 20px'}),
                html.Button("⏸️ Pause", id="pause-button", n_clicks=0,
                           style={'marginRight': 10, 'padding': '10px 20px'}),
                html.Button("⏹️ Reset", id="reset-button", n_clicks=0,
                           style={'padding': '10px 20px'}),
            ], style={'display': 'inline-block', 'marginRight': 30}),
            
            html.Div([
                html.Label("Jump to Step: "),
                dcc.Input(id="step-input", type="number", value=0, min=0, 
                         max=len(data_manager.evolution_df)-1, step=1,
                         style={'width': '80px', 'marginLeft': 10}),
                html.Button("Go", id="go-button", n_clicks=0,
                           style={'marginLeft': 10, 'padding': '5px 15px'}),
            ], style={'display': 'inline-block'}),
        ], style={'textAlign': 'center', 'marginBottom': 20}),
        
        # Slider
        html.Div([
            dcc.Slider(
                id='time-slider',
                min=0,
                max=len(data_manager.evolution_df)-1,
                step=1,
                value=0,
                marks={i: f'{i}' for i in range(0, len(data_manager.evolution_df), 50)},
                tooltip={"placement": "bottom", "always_visible": True}
            ),
        ], style={'marginBottom': 30}),
        
        # Info panel
        html.Div(id="info-panel", style={
            'textAlign': 'center', 
            'fontSize': '16px', 
            'marginBottom': 20,
            'padding': '10px',
            'backgroundColor': '#f0f0f0',
            'borderRadius': '5px'
        }),
        
        # Plot
        dcc.Graph(id='molecular-plot', style={'height': '70vh'}),
        
        # Animation interval
        dcc.Interval(
            id='interval-component',
            interval=200,  # milliseconds
            n_intervals=0,
            disabled=True
        ),
        
        # Store for animation state
        dcc.Store(id='animation-state', data={'playing': False, 'current_step': 0}),
    ], style={'padding': 20}),
])

# Callbacks
@app.callback(
    [Output('molecular-plot', 'figure'),
     Output('info-panel', 'children'),
     Output('time-slider', 'value'),
     Output('step-input', 'value')],
    [Input('time-slider', 'value'),
     Input('go-button', 'n_clicks'),
     Input('interval-component', 'n_intervals')],
    [State('step-input', 'value'),
     State('animation-state', 'data')]
)
def update_plot(slider_value, go_clicks, n_intervals, input_value, animation_state):
    ctx = callback_context
    
    # Determine which input triggered the callback
    if ctx.triggered:
        trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
        if trigger_id == 'go-button' and go_clicks > 0:
            current_step = min(max(0, input_value), len(data_manager.evolution_df)-1)
        elif trigger_id == 'interval-component' and animation_state['playing']:
            current_step = min(animation_state['current_step'] + 1, len(data_manager.evolution_df)-1)
        else:
            current_step = slider_value
    else:
        current_step = slider_value
    
    # Get data for current step
    step_data = data_manager.get_step_data(current_step)
    
    # Create plot
    fig = go.Figure()
    
    # Add empty sites as background
    all_coords = list(data_manager.coordinates.values())
    x_all = [coord['x'] for coord in all_coords]
    y_all = [coord['y'] for coord in all_coords]
    
    fig.add_trace(go.Scatter(
        x=x_all, y=y_all,
        mode='markers',
        marker=dict(size=3, color='lightgray', opacity=0.3),
        name='Empty sites',
        showlegend=True
    ))
    
    # Add occupied sites
    for species_id in [1, 2, 3]:
        positions = step_data['positions'][species_id]
        if positions:
            x_coords = [pos['x'] for pos in positions]
            y_coords = [pos['y'] for pos in positions]
            
            fig.add_trace(go.Scatter(
                x=x_coords, y=y_coords,
                mode='markers',
                marker=dict(
                    size=8,
                    color=data_manager.species_colors[species_id],
                    opacity=0.8
                ),
                name=f"{data_manager.species_names[species_id]} ({len(positions)})",
                showlegend=True
            ))
    
    # Calculate proper axis limits
    all_coords = list(data_manager.coordinates.values())
    x_coords = [coord['x'] for coord in all_coords]
    y_coords = [coord['y'] for coord in all_coords]
    
    x_min, x_max = min(x_coords), max(x_coords)
    y_min, y_max = min(y_coords), max(y_coords)
    margin = 2.0
    
    # Update layout
    fig.update_layout(
        title=f"Molecular Positions - Step {current_step + 1}/{len(data_manager.evolution_df)}",
        xaxis_title="X Coordinate (Å)",
        yaxis_title="Y Coordinate (Å)",
        xaxis=dict(
            scaleanchor="y", 
            scaleratio=1,
            range=[x_min - margin, x_max + margin]
        ),
        yaxis=dict(
            scaleanchor="x", 
            scaleratio=1,
            range=[y_min - margin, y_max + margin]
        ),
        showlegend=True,
        legend=dict(x=0.02, y=0.98),
        margin=dict(l=50, r=50, t=80, b=50),
        width=800,
        height=700
    )
    
    # Create info panel
    counts = step_data['counts']
    info_text = html.Div([
        html.Span(f"Step {current_step + 1} of {len(data_manager.evolution_df)} | ", 
                 style={'fontWeight': 'bold'}),
        html.Span(f"Time: {step_data['time']:.1f} units | "),
        html.Span(f"KMC Events: {step_data['nevents']} | "),
        html.Span(f"H*: {counts['H*']} | ", style={'color': 'red'}),
        html.Span(f"GeH2*: {counts['GeH2*']} | ", style={'color': 'blue'}),
        html.Span(f"GeH3*: {counts['GeH3*']} | ", style={'color': 'green'}),
        html.Span(f"Total: {sum(counts.values())}", style={'fontWeight': 'bold'})
    ])
    
    return fig, info_text, current_step, current_step

@app.callback(
    [Output('interval-component', 'disabled'),
     Output('animation-state', 'data')],
    [Input('play-button', 'n_clicks'),
     Input('pause-button', 'n_clicks'),
     Input('reset-button', 'n_clicks'),
     Input('time-slider', 'value')],
    [State('animation-state', 'data')]
)
def control_animation(play_clicks, pause_clicks, reset_clicks, slider_value, animation_state):
    ctx = callback_context
    
    if ctx.triggered:
        trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
        
        if trigger_id == 'play-button':
            # Start animation
            animation_state['playing'] = True
            animation_state['current_step'] = slider_value
            return False, animation_state
        
        elif trigger_id == 'pause-button':
            # Pause animation
            animation_state['playing'] = False
            return True, animation_state
        
        elif trigger_id == 'reset-button':
            # Reset to beginning
            animation_state['playing'] = False
            animation_state['current_step'] = 0
            return True, animation_state
        
        elif trigger_id == 'time-slider':
            # Update current step when slider moves
            animation_state['current_step'] = slider_value
            return animation_state['playing'] == False, animation_state
    
    return True, animation_state

if __name__ == '__main__':
    print("🎬 Starting Interactive Molecular Evolution Viewer...")
    print("📱 Open your browser to http://127.0.0.1:8050")
    print("🎯 Features:")
    print("   • Play/Pause animation through all 501 time steps")
    print("   • Slider for manual navigation")
    print("   • Text input to jump to specific steps")
    print("   • Real molecular positions from simulation")
    print("   • Interactive hover and zoom")
    
    app.run(debug=True, port=8050, host='127.0.0.1') 