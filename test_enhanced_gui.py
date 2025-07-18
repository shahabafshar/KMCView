#!/usr/bin/env python3
"""
Test script for Enhanced GUI Viewer
"""

import os
import sys

def test_enhanced_gui():
    """Test the enhanced GUI viewer initialization"""
    print("🧪 Testing Enhanced GUI Viewer...")
    
    try:
        # Import the module
        print("📦 Importing enhanced_gui_viewer...")
        import enhanced_gui_viewer
        print("✅ Import successful")
        
        # Test data loading without GUI
        print("📊 Testing data loading...")
        viewer = enhanced_gui_viewer.EnhancedMolecularViewer.__new__(enhanced_gui_viewer.EnhancedMolecularViewer)
        
        # Initialize basic attributes
        viewer.data_dir = 'input-output'
        viewer.theme = 'dark'
        viewer.evolution_df = None
        viewer.coordinates = {}
        viewer.lattice_bounds = {}
        viewer.lattice_info = {}
        viewer.position_cache = {}
        viewer.max_position_cache = 256
        
        # Test data loading
        try:
            viewer.load_data()
            print(f"✅ Data loading test passed - {len(viewer.evolution_df)} time points loaded")
        except Exception as e:
            print(f"❌ Data loading failed: {e}")
            return False
        
        # Test lattice loading
        try:
            viewer.load_positions()
            print(f"✅ Lattice loading test passed - {len(viewer.coordinates)} sites generated")
        except Exception as e:
            print(f"❌ Lattice loading failed: {e}")
            return False
        
        print("🎉 All tests passed!")
        return True
        
    except Exception as e:
        print(f"❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == '__main__':
    success = test_enhanced_gui()
    sys.exit(0 if success else 1) 