# 🌐 Zacros Molecular Evolution Viewer - Google Colab Guide

This guide explains how to run the Zacros molecular evolution viewer in Google Colab notebooks.

## 🚀 Quick Start

### Option 1: Use the Demo Notebook (Recommended)
1. Open the `colab_demo.ipynb` file in Google Colab
2. Follow the step-by-step instructions
3. Run each cell to see the viewer in action

### Option 2: Manual Setup
1. Upload `zacros_colab_viewer.py` to your Colab notebook
2. Run the installation and setup cells
3. Launch the viewer

## 📋 Requirements

The viewer requires these Python packages (automatically installed in Colab):
- `dash` - Web framework for interactive viewer
- `plotly` - Interactive plotting library
- `pandas` - Data manipulation
- `numpy` - Numerical computing
- `matplotlib` - Static plotting
- `scipy` - Scientific computing

## 🔧 Installation

```python
# Install required packages
!pip install dash plotly pandas numpy matplotlib scipy

# Upload the viewer code
from google.colab import files
uploaded = files.upload()

# Import the viewer
from zacros_colab_viewer import *
```

## 📊 Usage

### 1. Generate Sample Data
```python
# Create sample Zacros simulation data
sample_df = create_sample_data()
print(f"Created {len(sample_df)} time points")
```

### 2. Run Matplotlib Viewer
```python
# Static plots for notebook display
viewer = run_matplotlib_viewer()
```

### 3. Launch Web Viewer
```python
# Interactive web application
web_viewer = run_web_viewer(port=8050)
```

## 🌐 Web Viewer Access

### Automatic Proxy (Recommended)
When you run the web viewer, Colab will automatically provide a proxy link in the output. Click this link to access the viewer.

### Manual Proxy Setup
If automatic proxy doesn't work:

```python
# Install ngrok for manual tunneling
!pip install pyngrok
from pyngrok import ngrok

# Create tunnel
public_url = ngrok.connect(8050)
print(f"Web viewer available at: {public_url}")
```

## 📁 Using Your Own Data

### Upload Zacros Output Files
1. Upload your `specnum_output.txt` and `lattice_input.dat` files
2. Place them in the `input-output/` directory
3. Run the viewer with your data

```python
# Upload your files
from google.colab import files
uploaded = files.upload()

# Move to input-output directory
!mkdir -p input-output
!mv specnum_output.txt input-output/
!mv lattice_input.dat input-output/

# Run viewer with your data
viewer = run_matplotlib_viewer()
web_viewer = run_web_viewer()
```

### File Format Requirements

#### specnum_output.txt
```
# Step Nevents Time Temperature H* GeH2* GeH3*
0 25 0.000000 300.0 50 30 20
1 28 0.100000 300.0 52 28 22
2 31 0.200000 300.0 48 32 18
...
```

#### lattice_input.dat
```
# Lattice structure
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
```

## 🎮 Web Viewer Features

### Interactive Controls
- **Play/Pause**: Start/stop animation
- **Reset**: Return to beginning
- **Time Slider**: Navigate through time steps
- **Step Input**: Jump to specific step

### Visualization Elements
- **Lattice Grid**: Crystal structure background
- **Molecular Positions**: Real-time molecule display
- **Color Coding**: 
  - 🔴 Red: H* molecules
  - 🔵 Blue: GeH2* molecules
  - 🟢 Green: GeH3* molecules
- **Size Scaling**: Larger molecules appear bigger
- **Hover Information**: Click molecules for details

### Information Display
- Current simulation time
- Step number and total steps
- Real-time molecular counts

## 🔧 Troubleshooting

### Common Issues

#### Web Viewer Not Loading
**Problem**: The web viewer doesn't open or shows an error.

**Solutions**:
1. Check that the server started successfully
2. Look for the proxy link in the output
3. Try refreshing the page
4. Use manual ngrok proxy setup

```python
# Manual proxy setup
!pip install pyngrok
from pyngrok import ngrok
public_url = ngrok.connect(8050)
print(f"Access at: {public_url}")
```

#### No Data Displayed
**Problem**: The viewer shows empty plots or no molecules.

**Solutions**:
1. Verify `specnum_output.txt` exists in `input-output/` directory
2. Check file format matches expected structure
3. Ensure data contains non-zero molecular counts

```python
# Check data file
import pandas as pd
df = pd.read_csv('input-output/specnum_output.txt', sep=' ', comment='#')
print(df.head())
print(f"Total steps: {len(df)}")
```

#### Animation Not Working
**Problem**: Play button doesn't start animation.

**Solutions**:
1. Ensure you clicked the Play button
2. Check that animation interval is enabled
3. Verify data contains multiple time steps

#### Port Already in Use
**Problem**: Error about port 8050 being in use.

**Solutions**:
1. Use a different port number
2. Kill existing processes

```python
# Use different port
web_viewer = run_web_viewer(port=8051)

# Or kill existing processes
!pkill -f "python.*dash"
```

### Error Messages

#### "ModuleNotFoundError: No module named 'dash'"
```python
!pip install dash plotly
```

#### "FileNotFoundError: [Errno 2] No such file or directory"
```python
# Check if files exist
!ls -la input-output/
```

#### "Connection refused" or "Address already in use"
```python
# Use different port
web_viewer = run_web_viewer(port=8051)
```

## 📈 Performance Tips

### For Large Datasets
1. **Reduce animation speed**: Modify the interval in the web viewer
2. **Use static plots**: Use matplotlib viewer for large datasets
3. **Sample data**: Plot every nth step instead of all steps

```python
# Modify animation speed
# In the web viewer, change the interval parameter
dcc.Interval(id='animation-interval', interval=1000)  # 1 second instead of 500ms
```

### Memory Optimization
1. **Close unused plots**: Use `plt.close()` after displaying
2. **Clear variables**: Delete large dataframes when done
3. **Restart runtime**: If memory issues persist

```python
# Clear memory
import gc
plt.close('all')
del large_dataframe
gc.collect()
```

## 🔄 Colab Runtime Management

### Restart Runtime
If you encounter issues:
1. Go to Runtime → Restart runtime
2. Re-run installation cells
3. Re-upload your files

### Factory Reset
For persistent issues:
1. Go to Runtime → Factory reset runtime
2. Start fresh with the demo notebook

### GPU/TPU Acceleration
The viewer works on CPU, but you can enable GPU acceleration:
1. Go to Runtime → Change runtime type
2. Select GPU or TPU
3. Note: This won't significantly improve performance for this viewer

## 📚 Additional Resources

### Documentation
- [Dash Documentation](https://dash.plotly.com/)
- [Plotly Documentation](https://plotly.com/python/)
- [Google Colab Guide](https://colab.research.google.com/notebooks/basic_features_overview.ipynb)

### Example Workflows
1. **Data Analysis**: Use matplotlib viewer for static analysis
2. **Presentation**: Use web viewer for interactive demos
3. **Batch Processing**: Modify code for multiple datasets

### Customization
You can customize the viewer by modifying:
- Colors and sizes in `add_molecules()` method
- Animation speed in the interval parameter
- Plot layout in `configure_plot_layout()` method
- Data parsing in `load_evolution_data()` method

## 🆘 Getting Help

### Debug Information
Enable debug output by modifying the viewer code:

```python
# Add debug prints
print(f"Debug: Loaded {len(self.evolution_df)} data points")
print(f"Debug: Generated {len(self.coordinates)} lattice sites")
```

### Common Debug Commands
```python
# Check data structure
print(viewer.evolution_df.info())
print(viewer.lattice_info)

# Check file contents
!head -10 input-output/specnum_output.txt
!cat input-output/lattice_input.dat
```

### Support
For issues specific to this viewer:
1. Check the troubleshooting section above
2. Verify your data format matches the examples
3. Try the sample data first to confirm setup works

---

**Happy visualizing! 🚀**

The Zacros Molecular Evolution Viewer is now ready for your kinetic Monte Carlo simulation analysis in Google Colab. 