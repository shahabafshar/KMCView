# KMCView Project Structure (Minimalist)

This project is now streamlined to focus on the **Enhanced GUI Viewer** for kinetic Monte Carlo simulation data.

## 🎯 Current Active Files

```
KMCView/
├── 📄 enhanced_gui_viewer.py    # Main enhanced GUI application
├── 📄 test_enhanced_gui.py      # Tests for the GUI viewer
├── 📄 requirements.txt          # Python dependencies
├── 📄 README.md                 # Project documentation
├── 📄 LICENSE                   # MIT License
├── 📄 .gitignore               # Git ignore rules
├── 📁 input-output/            # KMC simulation data
├── 📁 backup/                  # Archive of old versions
├── 📁 .git/                    # Git repository
├── 📁 .venv/                   # Python virtual environment
└── 📁 .vscode/                 # VS Code settings
```

## 🚀 Usage

### Quick Start
```bash
# Install dependencies
pip install -r requirements.txt

# Run the enhanced GUI viewer
python enhanced_gui_viewer.py

# Or with options
python enhanced_gui_viewer.py --data-dir input-output --theme dark
```

### Features
- **Modern UI** with dark/light themes
- **Performance optimized** with caching
- **Keyboard shortcuts** (Space=Play/Pause, R=Reset, ←/→=Step, S=Save)
- **Interactive controls** (speed slider, step navigation)
- **Export functionality** for high-quality plots
- **Real-time performance monitoring**

## 📁 Backup Directory

All previous versions, Colab implementations, packaging files, and analysis scripts have been moved to the `backup/` directory to keep the main project clean and focused.

See `backup/README_BACKUP.md` for details on archived content.

## 🧪 Testing

```bash
# Run tests
python test_enhanced_gui.py
```

## 📋 Dependencies

- numpy >= 1.21.0
- pandas >= 1.3.0
- matplotlib >= 3.5.0
- scipy >= 1.7.0

## 🎨 Themes

The enhanced GUI viewer supports both dark and light themes:
- **Dark theme** (default): Modern dark interface
- **Light theme**: Clean light interface

Switch themes using the Theme selector in the GUI or command line:
```bash
python enhanced_gui_viewer.py --theme light
```

---

**Focus**: This streamlined structure allows developers to focus on the core enhanced GUI viewer while preserving development history in the backup directory. 