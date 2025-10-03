# dna-bukkets
A comprehensive genealogical DNA analysis application for identifying ancestral lines through chromosome segment triangulation.

## Overview

This tool was created to overcome limitations in commercial DNA analysis platforms (particularly FTDNA) by providing:
- **Unlimited data import** - Import all DNA matches at once, not just 7 at a time
- **Persistent database** - All data saved locally in SQLite for repeated analysis
- **Advanced triangulation** - Identify which DNA segments come from specific ancestral lines
- **Flexible filtering** - Filter by cM range, chromosome, family name, and more
- **Interactive visualization** - Generate chromosome browser visualizations with color-coded family lines
- **Manual ancestor assignment** - Review triangulation groups and assign them to specific ancestors

## Project History

This tool was developed iteratively to address specific genealogical research needs:
1. Started with basic CSV import and chromosome visualization
2. Added persistent database storage for managing large datasets (14,000+ matches)
3. Implemented triangulation detection to identify shared ancestral segments
4. Optimized for performance using multiprocessing algorithms
5. Added filtering capabilities to focus on relevant relationship ranges
6. Created simple debugging tools to validate triangulation accuracy
7. Built GUI for ease of use and packaged as standalone Windows executable

## Features

### Data Management
- Import single or multiple CSV files from FTDNA chromosome browser exports
- Automatic duplicate detection and removal
- SQLite database for persistent storage
- Edit family assignments for segments
- Delete unwanted segments
- Search and filter by match name, family, or chromosome
- Track import history

### Visualization
- Interactive chromosome browser using Plotly
- Color-coded segments by family name
- Filter by chromosome or family
- Hover for detailed segment information
- Zoom and pan capabilities
- Export to HTML for sharing

### Triangulation Analysis
- Multiple algorithm options:
  - **Simple (Debug)** - Step-by-step overlap detection with detailed logging
  - **Optimized (Fast)** - Multiprocessed analysis for large datasets
  - **Original (Slow)** - Legacy algorithm for comparison
- Adjustable parameters:
  - Minimum overlap (Mb) - How much DNA must overlap to group segments
  - Minimum group size - How many people must share a segment
  - Minimum cM - Exclude distant matches (noise reduction)
  - Maximum cM - Exclude close relatives (parent/child/sibling)
- Results sorted by group size (largest first)
- Family groupings within each triangulation group
- Detailed segment positions for verification
- Export results to CSV

### Ancestor Assignment Workflow
1. Import your DNA match data from FTDNA
2. Set filters to focus on relevant relationships:
   - Min cM: 35-50 (removes noise from distant cousins)
   - Max cM: 150 (excludes immediate family)
   - Min Overlap: 10-25 Mb (requires substantial overlap)
3. Run triangulation analysis
4. Review groups sorted by size (largest = strongest ancestral lines)
5. Manually assign each group to known ancestors based on:
   - Family names present in the group
   - Chromosome location consistency
   - Relationship strength (cM values)
6. Export results for documentation

## Installation

### Running from Source (Python Required)

**Requirements:**
- Python 3.6 or higher
- Required packages: pandas, plotly, numpy, tkinter (usually included with Python)

**Setup:**
```bash
# Install dependencies
pip install pandas plotly numpy

# Run the application
python run_dna_app.py
```

### Standalone Executable (No Python Required)

Download `DNA_Analysis_Tool.exe` - no installation or dependencies needed. Just run the executable.

## Usage Guide

### Getting Your Data from FTDNA

1. Log into your FTDNA account
2. Go to Family Finder > Chromosome Browser
3. Select matches you want to analyze (you can select all)
4. Click "Download Segments" to get CSV file
5. Repeat if needed - the tool will merge multiple exports automatically

### Importing Data

1. Launch DNA Analysis Tool
2. Go to **Data Management** tab
3. Click **Import CSV File** (single file) or **Import Multiple Files** (batch import)
4. Select your FTDNA CSV file(s)
5. The tool will import segments and show statistics

### Running Triangulation Analysis

1. Go to **Analysis** tab
2. Set parameters based on your goals:
   - **For strong ancestral lines:** Min cM: 50, Max cM: 150, Min Overlap: 15-25 Mb
   - **For medium strength lines:** Min cM: 35, Max cM: 150, Min Overlap: 10 Mb
   - **For all relationships:** Min cM: 20, Max cM: 0 (no max), Min Overlap: 5 Mb
3. Select algorithm: **Simple (Debug)** for detailed analysis or **Optimized (Fast)** for large datasets
4. Choose report format: **Detailed (For Assignment)** shows full names and positions
5. Click **Run Triangulation Analysis**
6. Review results - groups are sorted by size (largest first)

### Understanding Results

**Triangulation Group Example:**
```
RANK #1: 247 PEOPLE SHARE DNA HERE
Location: Chromosome 2, 16,451,390 - 48,303,182 bp (31.9 Mb)
Strength: Average 15.2 cM per person, Total 3,754 cM

PEOPLE IN THIS GROUP (by family name):
  Smith family (156 people):
    • John Smith (45.2 cM, 179.6-216.8Mb)
    • Mary Smith (32.1 cM, 180.2-216.9Mb)
    ...
```

**What this means:**
- 247 people share overlapping DNA in this chromosome region
- This likely represents a common ancestor
- Most are from Smith family - suggests a Smith ancestral line
- You can manually assign this to "Great-great-grandfather James Smith" or similar

### Visualizing Results

1. Go to **Visualization** tab
2. Select chromosome(s) to view (or "Autosomal Only" for all)
3. Optionally filter by family name
4. Click **Generate Chromosome Browser**
5. Interactive plot opens in your web browser
6. Colors represent different family lines
7. Hover over segments for details

### Best Practices

**Filtering Strategy:**
- Use Max cM to exclude immediate family (they create "bridges" between unrelated lines)
- Start with higher Min cM (50+) to see strongest patterns first
- Gradually lower Min cM to discover more distant relationships
- Increase Min Overlap for more stringent grouping

**Data Quality:**
- Regularly review and correct family assignments in Data Management tab
- FTDNA surname extraction isn't always accurate
- Clean data = better triangulation results

**Interpreting Groups:**
- Large groups (100+ people) usually represent well-documented lines
- Pure family groups (all one surname) are easiest to assign
- Mixed family groups may represent marriages between families
- Very large region spans (100+ Mb) may indicate algorithm issues - increase Min Overlap

## Technical Architecture

### Components

**dna_database.py** - SQLite database management
- Segment storage and retrieval
- Family assignment tracking
- Import history
- Search and filtering

**painter.py** - Visualization engine
- Chromosome browser plotting
- Color assignment by family
- Plotly-based interactive charts

**optimized_triangulation.py** - Fast triangulation engine
- Multiprocessed chromosome analysis
- Interval overlap algorithms
- Connected components detection
- Optimized for large datasets (10,000+ segments)

**simple_triangulation.py** - Debug triangulation engine
- Step-by-step overlap validation
- Detailed logging for troubleshooting
- Clear, verifiable logic
- Used for algorithm verification

**dna_gui_main.py** - Main GUI application
- tkinter-based interface
- Tab-based navigation
- Threading for long-running operations
- Progress indicators

**run_dna_app.py** - Application launcher
- Dependency checking
- PyInstaller compatibility
- Multiprocessing support

### Database Schema

**dna_segments table:**
- Segment data from FTDNA exports
- Family assignments
- Import metadata
- Search indices

**family_assignments table:**
- Surname definitions
- Color codes
- Ancestral line notes

**triangulation_groups table:**
- Detected triangulation groups
- Group statistics
- Creation metadata

## Troubleshooting

**Issue: Groups include unrelated people**
- Solution: Increase Min Overlap parameter (try 15-25 Mb)
- Check if immediate family is creating "bridges" - use Max cM to exclude them

**Issue: No groups found**
- Solution: Lower Min Overlap (try 5 Mb) or Min cM threshold
- Ensure you have enough data imported (need multiple matches per chromosome)

**Issue: Analysis is very slow**
- Solution: Use "Optimized (Fast)" algorithm instead of "Original (Slow)"
- Filter data before analysis (use Min/Max cM to reduce dataset size)

**Issue: Exe won't start**
- Solution: Windows may flag PyInstaller executables - add antivirus exception
- Try running with admin privileges
- Check if antivirus quarantined the file

**Issue: Visualization doesn't open**
- Solution: Ensure default web browser is set
- Check if popup blocker is preventing browser window
- Look for HTML file in system temp folder

## Building from Source

To create the standalone executable:

```bash
# Install PyInstaller
pip install pyinstaller

# Remove incompatible package if using Anaconda
conda remove pathlib

# Build executable
pyinstaller dna_tool.spec

# Output will be in dist/DNA_Analysis_Tool.exe
```

## File Formats

### Input: FTDNA CSV Format
```csv
Match Name,Chromosome,Start Location,End Location,Centimorgans,Matching SNPs
John Smith,1,12345678,23456789,25.5,5000
```

### Output: Triangulation Results CSV
```csv
Group_ID,Group_Size,Match_Name,Surname,Chromosome,Start_Location,End_Location,Centimorgans,Matching_SNPs
1,247,John Smith,Smith,2,16451390,48303182,45.2,8798
```

## Data Privacy

- All data is stored locally on your computer
- No internet connection required after initial setup
- Database file (dna_analysis.db) stays on your machine
- Visualization files are temporary and can be deleted
- No data is sent to external servers

## Limitations

- Designed for FTDNA chromosome browser CSV format
- Other platforms (23andMe, AncestryDNA) require format conversion
- Triangulation assumes segments that overlap come from same ancestor (not always true)
- Cannot detect false positives from identical-by-state (IBS) segments
- Requires manual verification and ancestor assignment
- Large datasets (50,000+ segments) may require significant processing time

## Future Enhancements

Potential features for future versions:
- GEDCOM import for automatic ancestor assignment
- Tree visualization showing ancestral lines
- Segment painting on chromosome ideograms
- Integration with other DNA platforms
- Cloud storage option for collaboration
- Statistical confidence scores for triangulation groups
- Automatic ancestor suggestion using known relationships

## Credits

Developed through iterative collaboration focused on:
- Practical genealogical research needs
- Handling large-scale DNA match datasets
- Accurate triangulation algorithms
- User-friendly interface design
- Performance optimization for real-world data

## License

This tool is provided as-is for personal genealogical research use.

## Support

For questions or issues:
1. Check the Troubleshooting section above
2. Review the console output (if available) for error messages
3. Try the Simple (Debug) algorithm to see detailed analysis steps
4. Verify your CSV file format matches FTDNA chromosome browser exports

## Version History

**v1.0** - Initial release
- Basic CSV import and visualization
- Simple triangulation detection
- Manual family assignment

**v1.5** - Database and optimization
- SQLite persistent storage
- Multiprocessed triangulation
- cM filtering
- Export capabilities

**v2.0** - GUI and packaging
- Full tkinter interface
- Tab-based navigation
- Standalone executable
- Debug tools
- Advanced filtering (Min/Max cM)
- Multiple algorithm options

---

**Happy ancestor hunting!**
