import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import numpy as np
from collections import defaultdict
import colorsys
import os

class DNAChromosomeBrowser:
    def __init__(self):
        self.data = None
        self.family_colors = {}
        self.triangulation_groups = []
        self.chromosome_lengths = {
            1: 247249719, 2: 242193529, 3: 199501827, 4: 191273063, 5: 180857866,
            6: 170899992, 7: 158821424, 8: 146274826, 9: 140273252, 10: 135374737,
            11: 134452384, 12: 132349534, 13: 114142980, 14: 106368585, 15: 100338915,
            16: 88827254, 17: 78774742, 18: 76117153, 19: 63811651, 20: 62435964,
            21: 46944323, 22: 49691432, 'X': 154913754, 'Y': 57772954
        }
    
    def load_csv(self, file_path):
        """Load DNA segment data from CSV file"""
        try:
            self.data = pd.read_csv(file_path)
            print(f"Loaded {len(self.data)} segments from {file_path}")
            print(f"Columns: {list(self.data.columns)}")
            print(f"Unique matches: {self.data['Match Name'].nunique()}")
            print(f"Chromosomes: {sorted(self.data['Chromosome'].unique())}")
            return True
        except Exception as e:
            print(f"Error loading CSV: {e}")
            return False
    
    def merge_csv_files(self, file_paths):
        """Merge multiple CSV files from FTDNA exports"""
        all_data = []
        for file_path in file_paths:
            try:
                df = pd.read_csv(file_path)
                df['Source_File'] = os.path.basename(file_path)
                all_data.append(df)
                print(f"Loaded {len(df)} segments from {file_path}")
            except Exception as e:
                print(f"Error loading {file_path}: {e}")
        
        if all_data:
            self.data = pd.concat(all_data, ignore_index=True)
            # Remove duplicates based on all columns except Source_File
            cols_to_check = [col for col in self.data.columns if col != 'Source_File']
            self.data = self.data.drop_duplicates(subset=cols_to_check)
            print(f"Total segments after merging and deduplication: {len(self.data)}")
            return True
        return False
    
    def extract_surname(self, match_name):
        """Extract surname from match name"""
        # Simple surname extraction - take the last word
        return match_name.split()[-1] if match_name else "Unknown"
    
    def assign_family_colors(self):
        """Assign unique colors to each family/surname"""
        surnames = set(self.extract_surname(name) for name in self.data['Match Name'].unique())
        
        # Generate distinct colors
        n_colors = len(surnames)
        colors = []
        for i in range(n_colors):
            hue = i / n_colors
            rgb = colorsys.hsv_to_rgb(hue, 0.8, 0.9)
            colors.append(f'rgb({int(rgb[0]*255)}, {int(rgb[1]*255)}, {int(rgb[2]*255)})')
        
        self.family_colors = dict(zip(surnames, colors))
        return self.family_colors
    
    def detect_triangulation_groups(self, min_overlap_bp=1000000, min_group_size=2):
        """Detect groups of segments that likely triangulate (same ancestral line)"""
        triangulation_groups = []
        
        # Group by chromosome
        for chromosome in self.data['Chromosome'].unique():
            chr_data = self.data[self.data['Chromosome'] == chromosome].copy()
            chr_data = chr_data.sort_values('Start Location')
            
            # Find overlapping segments
            groups = []
            for i, row1 in chr_data.iterrows():
                current_group = [row1]
                
                for j, row2 in chr_data.iterrows():
                    if i >= j:
                        continue
                    
                    # Check for overlap
                    overlap_start = max(row1['Start Location'], row2['Start Location'])
                    overlap_end = min(row1['End Location'], row2['End Location'])
                    overlap_length = max(0, overlap_end - overlap_start)
                    
                    if overlap_length >= min_overlap_bp:
                        # Check if row2 is already in a group with row1
                        found_in_group = False
                        for member in current_group:
                            if member.name == row2.name:
                                found_in_group = True
                                break
                        
                        if not found_in_group:
                            current_group.append(row2)
                
                if len(current_group) >= min_group_size:
                    # Check if this group is already captured
                    group_names = set(member['Match Name'] for member in current_group)
                    is_duplicate = False
                    for existing_group in groups:
                        existing_names = set(member['Match Name'] for member in existing_group)
                        if group_names == existing_names:
                            is_duplicate = True
                            break
                    
                    if not is_duplicate:
                        groups.append(current_group)
            
            triangulation_groups.extend(groups)
        
        self.triangulation_groups = triangulation_groups
        print(f"Found {len(triangulation_groups)} potential triangulation groups")
        return triangulation_groups
    
    def create_chromosome_browser(self, chromosomes=None, show_triangulation=True):
        """Create interactive chromosome browser visualization"""
        if self.data is None:
            print("No data loaded. Please load CSV first.")
            return None
        
        if chromosomes is None:
            chromosomes = sorted([c for c in self.data['Chromosome'].unique() if c != 'Y'])
        
        self.assign_family_colors()
        
        # Create subplots for each chromosome
        n_chromosomes = len(chromosomes)
        fig = make_subplots(
            rows=n_chromosomes, cols=1,
            subplot_titles=[f'Chromosome {chr}' for chr in chromosomes],
            vertical_spacing=0.02,
            shared_xaxes=False
        )
        
        # Add segments for each chromosome
        for i, chromosome in enumerate(chromosomes, 1):
            chr_data = self.data[self.data['Chromosome'] == chromosome]
            chr_length = self.chromosome_lengths.get(chromosome, 250000000)
            
            # Add chromosome background
            fig.add_trace(
                go.Scatter(
                    x=[0, chr_length],
                    y=[i-0.4, i-0.4],
                    mode='lines',
                    line=dict(color='lightgray', width=20),
                    showlegend=False,
                    hoverinfo='skip'
                ),
                row=i, col=1
            )
            
            # Add DNA segments
            for idx, segment in chr_data.iterrows():
                surname = self.extract_surname(segment['Match Name'])
                color = self.family_colors.get(surname, 'blue')
                
                # Create hover text
                hover_text = (
                    f"Match: {segment['Match Name']}<br>"
                    f"Family: {surname}<br>"
                    f"Chr {chromosome}: {segment['Start Location']:,} - {segment['End Location']:,}<br>"
                    f"Length: {(segment['End Location'] - segment['Start Location'])/1000000:.1f} Mb<br>"
                    f"cM: {segment['Centimorgans']}<br>"
                    f"SNPs: {segment['Matching SNPs']:,}"
                )
                
                fig.add_trace(
                    go.Scatter(
                        x=[segment['Start Location'], segment['End Location']],
                        y=[i, i],
                        mode='lines',
                        line=dict(color=color, width=8),
                        name=surname,
                        showlegend=(chromosome == chromosomes[0]),  # Only show legend for first chromosome
                        hovertemplate=hover_text + '<extra></extra>',
                        legendgroup=surname
                    ),
                    row=i, col=1
                )
        
        # Update layout
        fig.update_layout(
            height=150 * n_chromosomes,
            title="DNA Chromosome Browser - Segments by Family Line",
            showlegend=True,
            legend=dict(
                orientation="v",
                yanchor="top",
                y=1,
                xanchor="left",
                x=1.02
            )
        )
        
        # Update x-axes
        for i in range(1, n_chromosomes + 1):
            fig.update_xaxes(
                title_text="Position (bp)" if i == n_chromosomes else "",
                row=i, col=1,
                tickformat=",d"
            )
            fig.update_yaxes(
                showticklabels=False,
                row=i, col=1
            )
        
        return fig
    
    def create_triangulation_report(self):
        """Create a report of detected triangulation groups"""
        if not self.triangulation_groups:
            print("No triangulation groups detected. Run detect_triangulation_groups() first.")
            return
        
        print(f"\n=== TRIANGULATION ANALYSIS REPORT ===")
        print(f"Found {len(self.triangulation_groups)} potential triangulation groups\n")
        
        for i, group in enumerate(self.triangulation_groups, 1):
            print(f"Group {i}:")
            print(f"  Chromosome: {group[0]['Chromosome']}")
            
            # Calculate group boundaries
            min_start = min(member['Start Location'] for member in group)
            max_end = max(member['End Location'] for member in group)
            
            print(f"  Region: {min_start:,} - {max_end:,} bp")
            print(f"  Length: {(max_end - min_start)/1000000:.1f} Mb")
            print(f"  Members ({len(group)}):")
            
            for member in group:
                surname = self.extract_surname(member['Match Name'])
                print(f"    - {member['Match Name']} ({surname})")
                print(f"      {member['Start Location']:,} - {member['End Location']:,} bp, {member['Centimorgans']} cM")
            
            # Check if all members are from same family
            surnames = set(self.extract_surname(member['Match Name']) for member in group)
            if len(surnames) == 1:
                print(f"  ✓ All members from {list(surnames)[0]} family")
            else:
                print(f"  ⚠ Mixed families: {', '.join(surnames)}")
            
            print()
    
    def export_triangulation_csv(self, filename="triangulation_groups.csv"):
        """Export triangulation groups to CSV"""
        if not self.triangulation_groups:
            print("No triangulation groups to export.")
            return
        
        export_data = []
        for group_id, group in enumerate(self.triangulation_groups, 1):
            for member in group:
                export_data.append({
                    'Group_ID': group_id,
                    'Match_Name': member['Match Name'],
                    'Surname': self.extract_surname(member['Match Name']),
                    'Chromosome': member['Chromosome'],
                    'Start_Location': member['Start Location'],
                    'End_Location': member['End Location'],
                    'Centimorgans': member['Centimorgans'],
                    'Matching_SNPs': member['Matching SNPs']
                })
        
        df = pd.DataFrame(export_data)
        df.to_csv(filename, index=False)
        print(f"Exported {len(export_data)} triangulation records to {filename}")

# Example usage
if __name__ == "__main__":
    # Initialize browser
    browser = DNAChromosomeBrowser()
    
    # Try to load the CSV file automatically
    csv_file = "201189_Chromosome_Browser_Results_20250928.csv"
    
    if browser.load_csv(csv_file):
        print("✓ CSV loaded successfully!")
        
        # Detect triangulation groups
        print("Detecting triangulation groups...")
        browser.detect_triangulation_groups(min_overlap_bp=1000000, min_group_size=2)
        
        # Create visualization
        print("Creating chromosome browser...")
        fig = browser.create_chromosome_browser()
        
        # Show the plot
        print("Opening visualization in browser...")
        fig.show()
        
        # Generate triangulation report
        print("\nGenerating triangulation report...")
        browser.create_triangulation_report()
        
        # Export results
        browser.export_triangulation_csv("triangulation_results.csv")
        
    else:
        print(f"❌ Could not find {csv_file}")
        print("Available options:")
        print("1. Put your CSV file in the same directory and rename it to '201189_Chromosome_Browser_Results_20250928.csv'")
        print("2. Or run interactively:")
        print("   python -i painter.py")
        print("   browser.load_csv('your_actual_filename.csv')")
        print("   browser.detect_triangulation_groups()")
        print("   fig = browser.create_chromosome_browser()")
        print("   fig.show()")