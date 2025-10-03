#!/usr/bin/env python3
"""
DNA Analysis Application Launcher
==================================

This script launches the DNA Analysis GUI application.

Requirements:
- Python 3.6+
- tkinter (usually included with Python)
- pandas
- plotly
- sqlite3 (included with Python)

Usage:
    python run_dna_app.py

Features:
- Import DNA segment data from CSV files (FTDNA format)
- Store data in persistent SQLite database
- Interactive data management (view, edit, delete segments)
- Chromosome browser visualization
- Triangulation analysis for identifying ancestral lines
- Family line assignment and color coding
"""

import sys
import os

def check_dependencies():
    """Check if required packages are installed"""
    required_packages = ['pandas', 'plotly']
    missing_packages = []
    
    for package in required_packages:
        try:
            __import__(package)
        except ImportError:
            missing_packages.append(package)
    
    if missing_packages:
        print("Missing required packages:")
        for package in missing_packages:
            print(f"  - {package}")
        print("\nInstall missing packages with:")
        print(f"  pip install {' '.join(missing_packages)}")
        return False
    
    return True

def check_files():
    """Check if required files exist - skip check when running as exe"""
    # Check if running as a PyInstaller bundle
    if getattr(sys, 'frozen', False):
        # Running as compiled executable - skip file checks
        return True
    
    # Running as normal Python script - check for files
    required_files = ['dna_database.py', 'painter.py', 'dna_gui_main.py']
    missing_files = []
    
    for file in required_files:
        if not os.path.exists(file):
            missing_files.append(file)
    
    if missing_files:
        print("Missing required files:")
        for file in missing_files:
            print(f"  - {file}")
        print("\nMake sure all files are in the same directory as this launcher.")
        return False
    
    return True

def main():
    """Main launcher function"""
    print("DNA Analysis Application")
    print("=" * 40)
    
    # Check dependencies
    print("Checking dependencies...")
    if not check_dependencies():
        sys.exit(1)
    
    # Check files (only when not compiled)
    if not getattr(sys, 'frozen', False):
        print("Checking required files...")
        if not check_files():
            sys.exit(1)
    
    print("Starting DNA Analysis GUI...")
    
    try:
        # Import and run the GUI
        from dna_gui_main import main as gui_main
        gui_main()
    except ImportError as e:
        print(f"Error importing GUI module: {e}")
        print("Make sure all required files are in the same directory.")
        sys.exit(1)
    except Exception as e:
        print(f"Error starting application: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()  # Required for PyInstaller on Windows
    main()