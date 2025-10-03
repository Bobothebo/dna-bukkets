import tkinter as tk
from tkinter import ttk, filedialog, messagebox, simpledialog
import pandas as pd
from dna_database import DNADatabase
from painter import DNAChromosomeBrowser
from optimized_triangulation import OptimizedTriangulationEngine
from simple_triangulation import SimpleTriangulationEngine
import webbrowser
import tempfile
import os
import threading

class DNAAnalysisGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("DNA Analysis Tool")
        self.root.geometry("1200x800")
        
        # Initialize database and analysis engines
        self.db = DNADatabase()
        self.browser = DNAChromosomeBrowser()
        self.optimized_engine = OptimizedTriangulationEngine()
        self.simple_engine = SimpleTriangulationEngine()
        
        # Create GUI
        self.create_widgets()
        self.refresh_data()
    
    def create_widgets(self):
        """Create the main GUI interface"""
        
        # Create notebook for tabs
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Tab 1: Data Management
        self.create_data_tab()
        
        # Tab 2: Visualization
        self.create_visualization_tab()
        
        # Tab 3: Analysis
        self.create_analysis_tab()
        
        # Tab 4: Database Info
        self.create_info_tab()
    
    def create_data_tab(self):
        """Create data management tab"""
        data_frame = ttk.Frame(self.notebook)
        self.notebook.add(data_frame, text="Data Management")
        
        # Import section
        import_frame = ttk.LabelFrame(data_frame, text="Import Data", padding=10)
        import_frame.pack(fill=tk.X, padx=10, pady=5)
        
        ttk.Button(import_frame, text="Import CSV File", 
                  command=self.import_csv_file).pack(side=tk.LEFT, padx=5)
        ttk.Button(import_frame, text="Import Multiple Files", 
                  command=self.import_multiple_files).pack(side=tk.LEFT, padx=5)
        
        # Data view section
        view_frame = ttk.LabelFrame(data_frame, text="Data View", padding=10)
        view_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        # Filter controls
        filter_frame = ttk.Frame(view_frame)
        filter_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(filter_frame, text="Filter by:").pack(side=tk.LEFT, padx=5)
        
        self.filter_type = ttk.Combobox(filter_frame, values=["All", "Family", "Match Name", "Chromosome"], 
                                       state="readonly", width=15)
        self.filter_type.set("All")
        self.filter_type.pack(side=tk.LEFT, padx=5)
        self.filter_type.bind('<<ComboboxSelected>>', self.on_filter_change)
        
        self.filter_value = ttk.Combobox(filter_frame, width=20)
        self.filter_value.pack(side=tk.LEFT, padx=5)
        self.filter_value.bind('<Return>', self.apply_filter)
        
        ttk.Button(filter_frame, text="Apply Filter", 
                  command=self.apply_filter).pack(side=tk.LEFT, padx=5)
        ttk.Button(filter_frame, text="Clear Filter", 
                  command=self.clear_filter).pack(side=tk.LEFT, padx=5)
        
        # Treeview for data display
        tree_frame = ttk.Frame(view_frame)
        tree_frame.pack(fill=tk.BOTH, expand=True, pady=5)
        
        # Create treeview with scrollbars
        self.tree = ttk.Treeview(tree_frame, columns=('ID', 'Match', 'Chr', 'Start', 'End', 'cM', 'SNPs', 'Family'), show='headings')
        
        # Define column headings and widths
        columns_config = [
            ('ID', 50), ('Match', 200), ('Chr', 60), ('Start', 100), 
            ('End', 100), ('cM', 80), ('SNPs', 80), ('Family', 120)
        ]
        
        for col, width in columns_config:
            self.tree.heading(col, text=col)
            self.tree.column(col, width=width)
        
        # Scrollbars
        v_scrollbar = ttk.Scrollbar(tree_frame, orient=tk.VERTICAL, command=self.tree.yview)
        h_scrollbar = ttk.Scrollbar(tree_frame, orient=tk.HORIZONTAL, command=self.tree.xview)
        self.tree.configure(yscrollcommand=v_scrollbar.set, xscrollcommand=h_scrollbar.set)
        
        # Pack treeview and scrollbars
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        v_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        h_scrollbar.pack(side=tk.BOTTOM, fill=tk.X)
        
        # Context menu for treeview
        self.tree.bind("<Button-3>", self.show_context_menu)  # Right click
        self.tree.bind("<Double-1>", self.edit_segment)  # Double click
        
        # Edit controls
        edit_frame = ttk.Frame(view_frame)
        edit_frame.pack(fill=tk.X, pady=5)
        
        ttk.Button(edit_frame, text="Edit Selected", 
                  command=self.edit_segment).pack(side=tk.LEFT, padx=5)
        ttk.Button(edit_frame, text="Delete Selected", 
                  command=self.delete_segment).pack(side=tk.LEFT, padx=5)
        ttk.Button(edit_frame, text="Refresh", 
                  command=self.refresh_data).pack(side=tk.LEFT, padx=5)
    
    def create_visualization_tab(self):
        """Create visualization tab"""
        viz_frame = ttk.Frame(self.notebook)
        self.notebook.add(viz_frame, text="Visualization")
        
        # Controls
        control_frame = ttk.LabelFrame(viz_frame, text="Visualization Controls", padding=10)
        control_frame.pack(fill=tk.X, padx=10, pady=5)
        
        # Chromosome selection
        ttk.Label(control_frame, text="Chromosomes:").pack(side=tk.LEFT, padx=5)
        self.chr_selection = ttk.Combobox(control_frame, values=["All", "Autosomal Only", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"], 
                                         state="readonly", width=15)
        self.chr_selection.set("Autosomal Only")
        self.chr_selection.pack(side=tk.LEFT, padx=5)
        
        # Family filter
        ttk.Label(control_frame, text="Family:").pack(side=tk.LEFT, padx=5)
        self.viz_family_filter = ttk.Combobox(control_frame, width=20)
        self.viz_family_filter.pack(side=tk.LEFT, padx=5)
        
        ttk.Button(control_frame, text="Generate Chromosome Browser", 
                  command=self.generate_chromosome_browser).pack(side=tk.LEFT, padx=10)
        
        # Info display
        info_frame = ttk.LabelFrame(viz_frame, text="Visualization Info", padding=10)
        info_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        self.viz_info = tk.Text(info_frame, wrap=tk.WORD, height=10)
        viz_scrollbar = ttk.Scrollbar(info_frame, orient=tk.VERTICAL, command=self.viz_info.yview)
        self.viz_info.configure(yscrollcommand=viz_scrollbar.set)
        
        self.viz_info.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        viz_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
    
    def create_analysis_tab(self):
        """Create analysis tab"""
        analysis_frame = ttk.Frame(self.notebook)
        self.notebook.add(analysis_frame, text="Analysis")
        
        # Triangulation controls
        triangulation_frame = ttk.LabelFrame(analysis_frame, text="Triangulation Analysis", padding=10)
        triangulation_frame.pack(fill=tk.X, padx=10, pady=5)
        
        # Parameters
        param_frame = ttk.Frame(triangulation_frame)
        param_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(param_frame, text="Min Overlap (Mb):").pack(side=tk.LEFT, padx=5)
        self.min_overlap = tk.StringVar(value="1.0")
        overlap_entry = ttk.Entry(param_frame, textvariable=self.min_overlap, width=10)
        overlap_entry.pack(side=tk.LEFT, padx=5)
        
        ttk.Label(param_frame, text="Min Group Size:").pack(side=tk.LEFT, padx=5)
        self.min_group_size = tk.StringVar(value="2")
        size_entry = ttk.Entry(param_frame, textvariable=self.min_group_size, width=10)
        size_entry.pack(side=tk.LEFT, padx=5)
        
        ttk.Label(param_frame, text="Min cM:").pack(side=tk.LEFT, padx=5)
        self.min_cm = tk.StringVar(value="0")
        cm_entry = ttk.Entry(param_frame, textvariable=self.min_cm, width=10)
        cm_entry.pack(side=tk.LEFT, padx=5)
        
        ttk.Label(param_frame, text="Max cM:").pack(side=tk.LEFT, padx=5)
        self.max_cm = tk.StringVar(value="0")
        max_cm_entry = ttk.Entry(param_frame, textvariable=self.max_cm, width=10)
        max_cm_entry.pack(side=tk.LEFT, padx=5)
        
        # Report format selection
        ttk.Label(param_frame, text="Report:").pack(side=tk.LEFT, padx=5)
        self.report_format = ttk.Combobox(param_frame, values=["Detailed (For Assignment)", "Quick Summary"], 
                                         state="readonly", width=20)
        self.report_format.set("Detailed (For Assignment)")
        self.report_format.pack(side=tk.LEFT, padx=5)
        
        ttk.Button(param_frame, text="Run Triangulation Analysis", 
                  command=self.run_triangulation).pack(side=tk.LEFT, padx=10)
        
        # Algorithm selection
        ttk.Label(param_frame, text="Algorithm:").pack(side=tk.LEFT, padx=5)
        self.algorithm_choice = ttk.Combobox(param_frame, values=["Simple (Debug)", "Optimized (Fast)", "Original (Slow)"], 
                                           state="readonly", width=15)
        self.algorithm_choice.set("Simple (Debug)")
        self.algorithm_choice.pack(side=tk.LEFT, padx=5)
        
        # Progress bar
        self.progress_var = tk.DoubleVar()
        self.progress_bar = ttk.Progressbar(triangulation_frame, variable=self.progress_var, 
                                          mode='indeterminate')
        self.progress_bar.pack(fill=tk.X, pady=5)
        
        # Export button
        export_frame = ttk.Frame(triangulation_frame)
        export_frame.pack(fill=tk.X, pady=5)
        
        ttk.Button(export_frame, text="Export Results to CSV", 
                  command=self.export_triangulation_results).pack(side=tk.LEFT, padx=5)
        
        # Results display
        results_frame = ttk.LabelFrame(analysis_frame, text="Analysis Results", padding=10)
        results_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        self.analysis_results = tk.Text(results_frame, wrap=tk.WORD)
        analysis_scrollbar = ttk.Scrollbar(results_frame, orient=tk.VERTICAL, command=self.analysis_results.yview)
        self.analysis_results.configure(yscrollcommand=analysis_scrollbar.set)
        
        self.analysis_results.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        analysis_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
    
    def create_info_tab(self):
        """Create database info tab"""
        info_frame = ttk.Frame(self.notebook)
        self.notebook.add(info_frame, text="Database Info")
        
        # Statistics
        stats_frame = ttk.LabelFrame(info_frame, text="Database Statistics", padding=10)
        stats_frame.pack(fill=tk.X, padx=10, pady=5)
        
        self.stats_text = tk.Text(stats_frame, height=8, wrap=tk.WORD)
        self.stats_text.pack(fill=tk.X)
        
        # Import history
        history_frame = ttk.LabelFrame(info_frame, text="Import History", padding=10)
        history_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        self.history_tree = ttk.Treeview(history_frame, columns=('File', 'Segments', 'Date'), show='headings')
        self.history_tree.heading('File', text='Filename')
        self.history_tree.heading('Segments', text='Segments Imported')
        self.history_tree.heading('Date', text='Import Date')
        
        self.history_tree.column('File', width=300)
        self.history_tree.column('Segments', width=150)
        self.history_tree.column('Date', width=200)
        
        self.history_tree.pack(fill=tk.BOTH, expand=True)
        
        # Refresh button
        ttk.Button(info_frame, text="Refresh Info", 
                  command=self.refresh_info).pack(pady=10)
    
    def import_csv_file(self):
        """Import single CSV file"""
        file_path = filedialog.askopenfilename(
            title="Select DNA CSV File",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        
        if file_path:
            try:
                imported, duplicates = self.db.import_csv_data(file_path)
                messagebox.showinfo("Import Complete", 
                                  f"Imported {imported} new segments\nSkipped {duplicates} duplicates")
                self.refresh_data()
                self.refresh_info()
            except Exception as e:
                messagebox.showerror("Import Error", f"Error importing file:\n{str(e)}")
    
    def import_multiple_files(self):
        """Import multiple CSV files"""
        file_paths = filedialog.askopenfilenames(
            title="Select DNA CSV Files",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        
        if file_paths:
            total_imported = 0
            total_duplicates = 0
            
            for file_path in file_paths:
                try:
                    imported, duplicates = self.db.import_csv_data(file_path)
                    total_imported += imported
                    total_duplicates += duplicates
                except Exception as e:
                    messagebox.showerror("Import Error", 
                                       f"Error importing {file_path}:\n{str(e)}")
            
            messagebox.showinfo("Import Complete", 
                              f"Total imported: {total_imported} new segments\nTotal skipped: {total_duplicates} duplicates")
            self.refresh_data()
            self.refresh_info()
    
    def refresh_data(self):
        """Refresh the data display"""
        # Clear existing data
        for item in self.tree.get_children():
            self.tree.delete(item)
        
        # Get data from database
        df = self.db.get_all_segments()
        
        # Populate treeview
        for _, row in df.iterrows():
            self.tree.insert('', 'end', values=(
                row['id'], row['match_name'], row['chromosome'],
                f"{row['start_location']:,}", f"{row['end_location']:,}",
                row['centimorgans'], f"{row['matching_snps']:,}", 
                row['assigned_family'] or ""
            ))
        
        # Update filter dropdown values
        families = self.db.get_all_families()
        self.filter_value['values'] = families
        self.viz_family_filter['values'] = ["All"] + families
        if not self.viz_family_filter.get():
            self.viz_family_filter.set("All")
    
    def refresh_info(self):
        """Refresh database info tab"""
        # Update statistics
        stats = self.db.get_database_stats()
        stats_text = "Database Statistics:\n\n"
        for key, value in stats.items():
            formatted_key = key.replace('_', ' ').title()
            if key == 'total_centimorgans':
                stats_text += f"{formatted_key}: {value:.2f} cM\n"
            else:
                stats_text += f"{formatted_key}: {value:,}\n"
        
        self.stats_text.delete('1.0', tk.END)
        self.stats_text.insert('1.0', stats_text)
        
        # Update import history
        for item in self.history_tree.get_children():
            self.history_tree.delete(item)
        
        history_df = self.db.get_import_history()
        for _, row in history_df.iterrows():
            self.history_tree.insert('', 'end', values=(
                row['filename'], row['segments_imported'], row['import_date']
            ))
    
    def on_filter_change(self, event=None):
        """Handle filter type change"""
        filter_type = self.filter_type.get()
        
        if filter_type == "Family":
            families = self.db.get_all_families()
            self.filter_value['values'] = families
        elif filter_type == "Chromosome":
            chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                          "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
            self.filter_value['values'] = chromosomes
        else:
            self.filter_value['values'] = []
        
        self.filter_value.set("")
    
    def apply_filter(self, event=None):
        """Apply filter to data view"""
        filter_type = self.filter_type.get()
        filter_value = self.filter_value.get()
        
        # Clear existing data
        for item in self.tree.get_children():
            self.tree.delete(item)
        
        if filter_type == "All" or not filter_value:
            df = self.db.get_all_segments()
        else:
            if filter_type == "Family":
                df = self.db.search_segments(filter_value, "family")
            elif filter_type == "Match Name":
                df = self.db.search_segments(filter_value, "match_name")
            elif filter_type == "Chromosome":
                df = self.db.search_segments(filter_value, "chromosome")
            else:
                df = self.db.get_all_segments()
        
        # Populate treeview with filtered data
        for _, row in df.iterrows():
            self.tree.insert('', 'end', values=(
                row['id'], row['match_name'], row['chromosome'],
                f"{row['start_location']:,}", f"{row['end_location']:,}",
                row['centimorgans'], f"{row['matching_snps']:,}", 
                row['assigned_family'] or ""
            ))
    
    def clear_filter(self):
        """Clear filter and show all data"""
        self.filter_type.set("All")
        self.filter_value.set("")
        self.refresh_data()
    
    def edit_segment(self, event=None):
        """Edit selected segment"""
        selection = self.tree.selection()
        if not selection:
            messagebox.showwarning("No Selection", "Please select a segment to edit")
            return
        
        item = self.tree.item(selection[0])
        segment_id = item['values'][0]
        current_family = item['values'][7]
        
        # Simple edit dialog
        new_family = tk.simpledialog.askstring("Edit Family Assignment", 
                                              f"Current family: {current_family}\nEnter new family:",
                                              initialvalue=current_family)
        
        if new_family is not None:
            self.db.update_segment_family(segment_id, new_family)
            self.refresh_data()
    
    def delete_segment(self):
        """Delete selected segment"""
        selection = self.tree.selection()
        if not selection:
            messagebox.showwarning("No Selection", "Please select a segment to delete")
            return
        
        item = self.tree.item(selection[0])
        segment_id = item['values'][0]
        match_name = item['values'][1]
        
        if messagebox.askyesno("Confirm Delete", f"Delete segment for {match_name}?"):
            self.db.delete_segment(segment_id)
            self.refresh_data()
    
    def show_context_menu(self, event):
        """Show context menu for treeview"""
        item = self.tree.identify_row(event.y)
        if item:
            self.tree.selection_set(item)
            
            context_menu = tk.Menu(self.root, tearoff=0)
            context_menu.add_command(label="Edit Family", command=self.edit_segment)
            context_menu.add_command(label="Delete Segment", command=self.delete_segment)
            context_menu.add_separator()
            context_menu.add_command(label="View Details", command=self.view_segment_details)
            
            try:
                context_menu.tk_popup(event.x_root, event.y_root)
            finally:
                context_menu.grab_release()
    
    def view_segment_details(self):
        """View detailed information about selected segment"""
        selection = self.tree.selection()
        if not selection:
            return
        
        item = self.tree.item(selection[0])
        segment_id = item['values'][0]
        
        # Get full segment data from database
        df = self.db.get_all_segments()
        segment = df[df['id'] == segment_id].iloc[0]
        
        details_window = tk.Toplevel(self.root)
        details_window.title(f"Segment Details - {segment['match_name']}")
        details_window.geometry("400x300")
        
        details_text = tk.Text(details_window, wrap=tk.WORD, padx=10, pady=10)
        details_text.pack(fill=tk.BOTH, expand=True)
        
        details_info = f"""Segment Details:

Match Name: {segment['match_name']}
Chromosome: {segment['chromosome']}
Start Location: {segment['start_location']:,} bp
End Location: {segment['end_location']:,} bp
Length: {(segment['end_location'] - segment['start_location'])/1000000:.2f} Mb
Centimorgans: {segment['centimorgans']} cM
Matching SNPs: {segment['matching_snps']:,}
Assigned Family: {segment['assigned_family'] or 'Not assigned'}
Source File: {segment['source_file'] or 'Unknown'}
Import Date: {segment['import_date']}
Notes: {segment['notes'] or 'None'}
"""
        
        details_text.insert('1.0', details_info)
        details_text.config(state=tk.DISABLED)
    
    def generate_chromosome_browser(self):
        """Generate chromosome browser visualization"""
        try:
            # Get data from database
            df = self.db.get_all_segments()
            
            if df.empty:
                messagebox.showwarning("No Data", "No segments found in database. Import some data first.")
                return
            
            self.viz_info.delete('1.0', tk.END)
            self.viz_info.insert(tk.END, f"Found {len(df)} segments in database\n")
            self.viz_info.insert(tk.END, f"Database columns: {list(df.columns)}\n\n")
            
            # Apply filters
            family_filter = self.viz_family_filter.get()
            if family_filter and family_filter != "All":
                df = df[df['assigned_family'] == family_filter]
                self.viz_info.insert(tk.END, f"Filtered to family '{family_filter}': {len(df)} segments\n")
            
            chr_selection = self.chr_selection.get()
            if chr_selection == "Autosomal Only":
                df = df[df['chromosome'].isin([str(i) for i in range(1, 23)])]
                self.viz_info.insert(tk.END, f"Filtered to autosomal chromosomes: {len(df)} segments\n")
            elif chr_selection != "All" and chr_selection in [str(i) for i in range(1, 23)] + ['X', 'Y']:
                df = df[df['chromosome'] == chr_selection]
                self.viz_info.insert(tk.END, f"Filtered to chromosome {chr_selection}: {len(df)} segments\n")
            
            if df.empty:
                messagebox.showwarning("No Data", "No segments match the selected filters.")
                return
            
            # Convert to format expected by DNAChromosomeBrowser
            # Create a copy and rename columns to match expected format
            browser_df = df.copy()
            browser_df = browser_df.rename(columns={
                'match_name': 'Match Name',
                'chromosome': 'Chromosome',
                'start_location': 'Start Location',
                'end_location': 'End Location',
                'centimorgans': 'Centimorgans',
                'matching_snps': 'Matching SNPs'
            })
            
            self.viz_info.insert(tk.END, f"Converted columns: {list(browser_df.columns)}\n")
            
            # Ensure all required columns exist
            required_columns = ['Match Name', 'Chromosome', 'Start Location', 'End Location', 'Centimorgans', 'Matching SNPs']
            missing_columns = [col for col in required_columns if col not in browser_df.columns]
            if missing_columns:
                error_msg = f"Missing required columns: {missing_columns}"
                self.viz_info.insert(tk.END, f"Error: {error_msg}\n")
                messagebox.showerror("Column Error", error_msg)
                return
            
            self.browser.data = browser_df
            
            # Generate visualization
            self.viz_info.insert(tk.END, "Generating chromosome browser visualization...\n")
            self.root.update()
            
            # Create the plot
            chromosomes = None
            if chr_selection != "All" and chr_selection != "Autosomal Only":
                chromosomes = [chr_selection]
            
            fig = self.browser.create_chromosome_browser(chromosomes=chromosomes)
            
            if fig:
                # Save to temporary HTML file and open in browser
                temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.html', delete=False)
                fig.write_html(temp_file.name)
                temp_file.close()
                
                webbrowser.open(f'file://{temp_file.name}')
                
                # Update info
                info_text = f"""Chromosome Browser Generated Successfully!

Data Summary:
- Total segments displayed: {len(browser_df)}
- Unique matches: {browser_df['Match Name'].nunique()}
- Chromosomes: {sorted(browser_df['Chromosome'].unique())}
- Families: {sorted(df['assigned_family'].dropna().unique())}
- Total cM: {browser_df['Centimorgans'].sum():.2f}

The visualization has opened in your web browser.
You can interact with the plot by:
- Hovering over segments for details
- Zooming in/out
- Panning across chromosomes
- Using the legend to hide/show families

"""
                self.viz_info.delete('1.0', tk.END)
                self.viz_info.insert('1.0', info_text)
            else:
                self.viz_info.insert(tk.END, "Error: Could not generate visualization.\n")
                
        except Exception as e:
            error_msg = f"Error generating visualization:\n{str(e)}"
            messagebox.showerror("Visualization Error", error_msg)
            self.viz_info.insert(tk.END, f"Error: {str(e)}\n")
            print(f"Debug - Full error: {e}")  # Console debugging
    
    def run_triangulation(self):
        """Run triangulation analysis with algorithm choice"""
        try:
            # Get parameters
            min_overlap_mb = float(self.min_overlap.get())
            min_group_size = int(self.min_group_size.get())
            min_cm = float(self.min_cm.get())
            max_cm = float(self.max_cm.get()) if self.max_cm.get() else float('inf')
            min_overlap_bp = int(min_overlap_mb * 1000000)
            algorithm = self.algorithm_choice.get()
            
            # Get data from database
            df = self.db.get_all_segments()
            
            if df.empty:
                messagebox.showwarning("No Data", "No segments found in database. Import some data first.")
                return
            
            # Show data size warning for large datasets
            if len(df) > 5000 and algorithm == "Original (Slow)":
                response = messagebox.askyesno(
                    "Large Dataset Warning", 
                    f"You have {len(df)} segments. The original algorithm may be very slow.\n\n"
                    "Recommended: Use 'Optimized (Fast)' algorithm for large datasets.\n\n"
                    "Continue with original algorithm?"
                )
                if not response:
                    return
            
            # Clear results display and show progress
            self.analysis_results.delete('1.0', tk.END)
            self.analysis_results.insert(tk.END, f"Starting triangulation analysis on {len(df)} segments...\n")
            self.analysis_results.insert(tk.END, f"Algorithm: {algorithm}\n")
            self.analysis_results.insert(tk.END, f"Parameters: {min_overlap_mb} Mb overlap, {min_group_size} min group size, {min_cm} min cM\n\n")
            
            # Apply cM filters to data
            if min_cm > 0 or max_cm < float('inf'):
                original_count = len(df)
                if min_cm > 0:
                    df = df[df['centimorgans'] >= min_cm]
                if max_cm < float('inf'):
                    df = df[df['centimorgans'] <= max_cm]
                filtered_count = len(df)
                
                filter_msg = f"Filtered by cM: {original_count} → {filtered_count} segments"
                if min_cm > 0 and max_cm < float('inf'):
                    filter_msg += f" (kept {min_cm}-{max_cm} cM range)"
                elif min_cm > 0:
                    filter_msg += f" (removed {original_count - filtered_count} segments < {min_cm} cM)"
                elif max_cm < float('inf'):
                    filter_msg += f" (removed {original_count - filtered_count} segments > {max_cm} cM)"
                
                self.analysis_results.insert(tk.END, filter_msg + "\n\n")
                
                if df.empty:
                    self.analysis_results.insert(tk.END, "No segments remain after cM filtering. Try adjusting the cM range.\n")
                    return
            
            # Start progress bar
            self.progress_bar.start()
            self.root.update()
            
            # Run analysis in separate thread to prevent GUI freezing
            def run_analysis():
                try:
                    if algorithm == "Simple (Debug)":
                        # Use simple algorithm with debugging
                        browser_df = df.rename(columns={
                            'match_name': 'Match Name',
                            'chromosome': 'Chromosome',
                            'start_location': 'Start Location',
                            'end_location': 'End Location',
                            'centimorgans': 'Centimorgans',
                            'matching_snps': 'Matching SNPs'
                        })
                        
                        groups = self.simple_engine.detect_triangulation_simple(
                            browser_df, min_overlap_bp, min_group_size, debug=True
                        )
                        
                        # Combine debug output with results
                        debug_text = "\n".join(self.simple_engine.debug_output)
                        report_text = self.simple_engine.create_triangulation_report_by_size(min_cm)
                        results_text = f"DEBUG OUTPUT:\n{debug_text}\n\n{'='*60}\n\n{report_text}"
                        
                    elif algorithm == "Optimized (Fast)":
                        # Convert to format expected by OptimizedTriangulationEngine
                        browser_df = df.rename(columns={
                            'match_name': 'Match Name',
                            'chromosome': 'Chromosome',
                            'start_location': 'Start Location',
                            'end_location': 'End Location',
                            'centimorgans': 'Centimorgans',
                            'matching_snps': 'Matching SNPs'
                        })
                        
                        # Run optimized analysis with debugging
                        groups = self.optimized_engine.detect_triangulation_groups_optimized(
                            browser_df, min_overlap_bp, min_group_size, debug=True
                        )
                        
                        # Generate report based on selected format
                        report_format = self.report_format.get()
                        if report_format == "Quick Summary":
                            results_text = self.optimized_engine.create_simple_group_list()
                        else:
                            results_text = self.optimized_engine.create_triangulation_report_by_size(min_cm)
                        
                    else:
                        # Use original algorithm
                        self.browser.data = df.rename(columns={
                            'match_name': 'Match Name',
                            'chromosome': 'Chromosome',
                            'start_location': 'Start Location',
                            'end_location': 'End Location',
                            'centimorgans': 'Centimorgans',
                            'matching_snps': 'Matching SNPs'
                        })
                        
                        groups = self.browser.detect_triangulation_groups(
                            min_overlap_bp=min_overlap_bp, 
                            min_group_size=min_group_size
                        )
                        
                        # Generate original report format
                        results_text = self.generate_original_report(groups)
                    
                    # Update GUI in main thread
                    self.root.after(0, self.update_triangulation_results, results_text)
                    
                except Exception as e:
                    error_msg = f"Error during analysis: {str(e)}"
                    self.root.after(0, self.show_triangulation_error, error_msg)
            
            # Start analysis thread
            analysis_thread = threading.Thread(target=run_analysis)
            analysis_thread.daemon = True
            analysis_thread.start()
            
        except ValueError as e:
            messagebox.showerror("Parameter Error", f"Invalid parameter values:\n{str(e)}")
        except Exception as e:
            messagebox.showerror("Analysis Error", f"Error starting triangulation analysis:\n{str(e)}")
    
    def update_triangulation_results(self, results_text):
        """Update triangulation results in main thread"""
        self.progress_bar.stop()
        self.analysis_results.delete('1.0', tk.END)
        self.analysis_results.insert('1.0', results_text)
    
    def show_triangulation_error(self, error_msg):
        """Show triangulation error in main thread"""
        self.progress_bar.stop()
        self.analysis_results.insert(tk.END, f"\nERROR: {error_msg}")
        messagebox.showerror("Analysis Error", error_msg)
    
    def generate_original_report(self, groups):
        """Generate report in original format for comparison"""
        if not groups:
            return "No triangulation groups found with the current parameters.\nTry reducing the minimum overlap or group size."
        
        results_text = f"TRIANGULATION ANALYSIS RESULTS (Original Algorithm)\n"
        results_text += f"{'='*60}\n\n"
        results_text += f"Found {len(groups)} potential triangulation groups:\n\n"
        
        for i, group in enumerate(groups, 1):
            results_text += f"GROUP {i}:\n"
            results_text += f"  Chromosome: {group[0]['Chromosome']}\n"
            
            # Calculate group boundaries
            min_start = min(member['Start Location'] for member in group)
            max_end = max(member['End Location'] for member in group)
            
            results_text += f"  Region: {min_start:,} - {max_end:,} bp\n"
            results_text += f"  Length: {(max_end - min_start)/1000000:.2f} Mb\n"
            results_text += f"  Members ({len(group)}):\n"
            
            for member in group:
                surname = self.browser.extract_surname(member['Match Name'])
                results_text += f"    • {member['Match Name']} ({surname})\n"
                results_text += f"      {member['Start Location']:,} - {member['End Location']:,} bp, {member['Centimorgans']} cM\n"
            
            # Check if all members are from same family
            surnames = set(self.browser.extract_surname(member['Match Name']) for member in group)
            if len(surnames) == 1:
                results_text += f"  ✓ All members from {list(surnames)[0]} family\n"
            else:
                results_text += f"  ⚠ Mixed families: {', '.join(surnames)}\n"
            
            results_text += "\n"
        
        return results_text
    
    def export_triangulation_results(self):
        """Export triangulation results to CSV"""
        if not hasattr(self.optimized_engine, 'triangulation_groups') or not self.optimized_engine.triangulation_groups:
            messagebox.showwarning("No Results", "No triangulation results to export. Run analysis first.")
            return
        
        filename = filedialog.asksaveasfilename(
            title="Save Triangulation Results",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        
        if filename:
            try:
                count = self.optimized_engine.export_optimized_results(filename)
                messagebox.showinfo("Export Complete", f"Exported {count} triangulation records to {filename}")
            except Exception as e:
                messagebox.showerror("Export Error", f"Error exporting results:\n{str(e)}")


def main():
    """Main function to run the GUI"""
    root = tk.Tk()
    app = DNAAnalysisGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()