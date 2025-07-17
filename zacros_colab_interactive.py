#!/usr/bin/env python3
"""
🔬 Colab Interactive Molecular Evolution Viewer
==============================================

Uses IPython widgets for interactive controls in Google Colab:
• Interactive sliders and buttons that work in notebooks
• Real-time plot updates
• No GUI backend required
• Full Colab compatibility
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import warnings
import random
from IPython.display import display, clear_output
import ipywidgets as widgets
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

class ColabInteractiveViewer:
    """Interactive viewer using IPython widgets for Colab"""
    
    def __init__(self, data_dir='input-output'):
        self.data_dir = data_dir
        self.evolution_df = pd.DataFrame()
        self.coordinates = {}
        self.lattice_bounds = {}
        self.lattice_info = {}
        
        # Animation state
        self.current_step = 0
        self.is_playing = False
        self.animation_timer = None
        
        # Load data
        self.load_evolution_data()
        self.load_lattice_structure()
        self.generate_lattice_sites()
        
        # Create widgets
        self.create_widgets()
        self.setup_callbacks()
        
        # Display the interface
        self.display_interface()
    
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
    
    def create_widgets(self):
        """Create IPython widgets for the interface"""
        # Title
        self.title_widget = widgets.HTML(
            value="<h2>🔬 Interactive Zacros Molecular Evolution Viewer</h2>",
            layout=widgets.Layout(margin='10px 0px')
        )
        
        # Step slider
        self.step_slider = widgets.IntSlider(
            value=0,
            min=0,
            max=max(0, len(self.evolution_df) - 1),
            step=1,
            description='Step:',
            continuous_update=False,
            layout=widgets.Layout(width='60%')
        )
        
        # Step input
        self.step_input = widgets.IntText(
            value=0,
            min=0,
            max=max(0, len(self.evolution_df) - 1),
            description='Step:',
            layout=widgets.Layout(width='200px')
        )
        
        # Control buttons
        self.play_button = widgets.Button(
            description='▶️ Play',
            button_style='success',
            layout=widgets.Layout(width='100px')
        )
        
        self.pause_button = widgets.Button(
            description='⏸️ Pause',
            button_style='warning',
            layout=widgets.Layout(width='100px')
        )
        
        self.reset_button = widgets.Button(
            description='🔄 Reset',
            button_style='danger',
            layout=widgets.Layout(width='100px')
        )
        
        # Animation speed slider
        self.speed_slider = widgets.FloatSlider(
            value=0.5,
            min=0.1,
            max=2.0,
            step=0.1,
            description='Speed:',
            layout=widgets.Layout(width='200px')
        )
        
        # Info display
        self.info_display = widgets.HTML(
            value="<p>Ready to start...</p>",
            layout=widgets.Layout(margin='10px 0px')
        )
        
        # Output widget for plots
        self.plot_output = widgets.Output()
    
    def setup_callbacks(self):
        """Set up widget callbacks"""
        # Step slider callback
        self.step_slider.observe(self.on_step_change, names='value')
        
        # Step input callback
        self.step_input.observe(self.on_step_input_change, names='value')
        
        # Button callbacks
        self.play_button.on_click(self.on_play_click)
        self.pause_button.on_click(self.on_pause_click)
        self.reset_button.on_click(self.on_reset_click)
        
        # Speed slider callback
        self.speed_slider.observe(self.on_speed_change, names='value')
    
    def on_step_change(self, change):
        """Handle step slider change"""
        self.current_step = change['new']
        self.step_input.value = self.current_step
        self.update_plot()
        self.update_info()
    
    def on_step_input_change(self, change):
        """Handle step input change"""
        self.current_step = max(0, min(change['new'], len(self.evolution_df) - 1))
        self.step_slider.value = self.current_step
        self.update_plot()
        self.update_info()
    
    def on_play_click(self, b):
        """Handle play button click"""
        self.is_playing = True
        self.start_animation()
    
    def on_pause_click(self, b):
        """Handle pause button click"""
        self.is_playing = False
        self.stop_animation()
    
    def on_reset_click(self, b):
        """Handle reset button click"""
        self.is_playing = False
        self.stop_animation()
        self.current_step = 0
        self.step_slider.value = 0
        self.step_input.value = 0
        self.update_plot()
        self.update_info()
    
    def on_speed_change(self, change):
        """Handle speed slider change"""
        if self.is_playing:
            self.stop_animation()
            self.start_animation()
    
    def start_animation(self):
        """Start the animation"""
        if not self.is_playing:
            return
        
        # Calculate interval based on speed
        interval = int(1000 / self.speed_slider.value)  # Convert to milliseconds
        
        # Use JavaScript timer for animation
        js_code = f"""
        function animateStep() {{
            if (window.isPlaying) {{
                // Trigger step increment
                IPython.notebook.kernel.execute('viewer.animate_step()');
                setTimeout(animateStep, {interval});
            }}
        }}
        window.isPlaying = true;
        animateStep();
        """
        
        # Execute JavaScript
        output.eval_js(js_code)
    
    def stop_animation(self):
        """Stop the animation"""
        js_code = "window.isPlaying = false;"
        output.eval_js(js_code)
    
    def animate_step(self):
        """Increment step for animation"""
        if self.is_playing and self.current_step < len(self.evolution_df) - 1:
            self.current_step += 1
            self.step_slider.value = self.current_step
            self.step_input.value = self.current_step
            self.update_plot()
            self.update_info()
        elif self.current_step >= len(self.evolution_df) - 1:
            self.is_playing = False
            self.stop_animation()
    
    def update_plot(self):
        """Update the molecular plot"""
        with self.plot_output:
            clear_output(wait=True)
            
            # Create the plot
            fig, ax = plt.subplots(figsize=(12, 10))
            
            # Get current data
            current_data = self.evolution_df.iloc[self.current_step]
            current_time = current_data['time']
            h_count = int(current_data['H*'])
            geh2_count = int(current_data['GeH2*'])
            geh3_count = int(current_data['GeH3*'])
            
            # Generate positions
            positions = self.generate_molecular_positions(self.current_step)
            
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
            ax.set_title(f'Molecular Evolution - Step {self.current_step+1}/{len(self.evolution_df)} (Time: {current_time:.1f})')
            ax.legend()
            ax.grid(True, alpha=0.3)
            ax.set_aspect('equal')
            
            plt.tight_layout()
            plt.show()
    
    def update_info(self):
        """Update the information display"""
        if self.current_step < len(self.evolution_df):
            current_data = self.evolution_df.iloc[self.current_step]
            current_time = current_data['time']
            h_count = int(current_data['H*'])
            geh2_count = int(current_data['GeH2*'])
            geh3_count = int(current_data['GeH3*'])
            
            info_html = f"""
            <div style="background-color: #f0f0f0; padding: 10px; border-radius: 5px;">
                <strong>Time:</strong> {current_time:.1f} | 
                <strong>Step:</strong> {self.current_step+1}/{len(self.evolution_df)} | 
                <span style="color: red;"><strong>H*:</strong> {h_count}</span> | 
                <span style="color: blue;"><strong>GeH2*:</strong> {geh2_count}</span> | 
                <span style="color: green;"><strong>GeH3*:</strong> {geh3_count}</span>
            </div>
            """
            
            self.info_display.value = info_html
    
    def display_interface(self):
        """Display the complete interface"""
        # Create layout
        controls = widgets.HBox([
            self.step_slider,
            widgets.VBox([
                self.step_input,
                self.speed_slider
            ])
        ])
        
        buttons = widgets.HBox([
            self.play_button,
            self.pause_button,
            self.reset_button
        ])
        
        # Display everything
        display(self.title_widget)
        display(controls)
        display(buttons)
        display(self.info_display)
        display(self.plot_output)
        
        # Show initial plot
        self.update_plot()
        self.update_info()

def run_interactive_viewer(data_dir='input-output'):
    """Run the interactive viewer in Colab"""
    print("🔬 Starting Interactive Colab Viewer...")
    print("📱 This viewer uses IPython widgets for notebook compatibility")
    
    viewer = ColabInteractiveViewer(data_dir)
    return viewer

# Main execution for Colab
if __name__ == '__main__':
    # Create sample data if needed
    SampleDataGenerator.create_sample_files()
    
    # Run interactive viewer
    viewer = run_interactive_viewer() 