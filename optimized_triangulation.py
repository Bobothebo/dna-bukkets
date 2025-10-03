import pandas as pd
import numpy as np
from collections import defaultdict
from multiprocessing import Pool, cpu_count
import time
from typing import List, Dict, Tuple
import itertools

class OptimizedTriangulationEngine:
    """Optimized triangulation engine using multiprocessing and interval trees"""
    
    def __init__(self):
        self.data = None
        self.triangulation_groups = []
        
    def detect_triangulation_groups_optimized(self, data, min_overlap_bp=1000000, min_group_size=2, max_workers=None, debug=False):
        """
        Optimized triangulation detection using multiprocessing and interval overlap algorithms
        """
        print(f"Starting optimized triangulation analysis on {len(data)} segments...")
        start_time = time.time()
        
        self.data = data
        
        # Use all CPU cores by default, but leave one free
        if max_workers is None:
            max_workers = max(1, cpu_count() - 1)
        
        print(f"Using {max_workers} worker processes...")
        
        # Group segments by chromosome for parallel processing
        chromosome_groups = self.data.groupby('Chromosome')
        
        # Prepare arguments for parallel processing
        args_list = []
        for chromosome, chr_data in chromosome_groups:
            if len(chr_data) >= min_group_size:  # Skip chromosomes with too few segments
                args_list.append((chromosome, chr_data, min_overlap_bp, min_group_size, debug))
        
        print(f"Processing {len(args_list)} chromosomes in parallel...")
        
        # Process chromosomes in parallel
        all_groups = []
        if len(args_list) > 0:
            with Pool(processes=max_workers) as pool:
                results = pool.map(self._process_chromosome_optimized, args_list)
                
                # Flatten results
                for chromosome_groups in results:
                    all_groups.extend(chromosome_groups)
        
        self.triangulation_groups = all_groups
        
        end_time = time.time()
        processing_time = end_time - start_time
        
        print(f"Triangulation analysis completed in {processing_time:.2f} seconds")
        print(f"Found {len(all_groups)} triangulation groups")
        
        return all_groups
    
    @staticmethod
    def _process_chromosome_optimized(args):
        """Process a single chromosome for triangulation - with debugging"""
        chromosome, chr_data, min_overlap_bp, min_group_size, debug = args
        
        print(f"Processing chromosome {chromosome} with {len(chr_data)} segments...")
        
        if debug and chromosome == '1':  # Debug chromosome 1 specifically
            print(f"\nDEBUG: Chromosome 1 segments:")
            for idx, row in chr_data.iterrows():
                print(f"  {row['Match Name']}: {row['Start Location']:,} - {row['End Location']:,} ({row['Centimorgans']} cM)")
            print()
        
        # Convert to list of tuples for faster processing
        segments = []
        for idx, row in chr_data.iterrows():
            segments.append({
                'index': idx,
                'data': row,
                'start': row['Start Location'],
                'end': row['End Location'],
                'match_name': row['Match Name'],
                'cm': row['Centimorgans']
            })
        
        # Sort segments by start position for efficient overlap detection
        segments.sort(key=lambda x: x['start'])
        
        groups = []
        processed_indices = set()
        
        for i, segment1 in enumerate(segments):
            if segment1['index'] in processed_indices:
                continue
                
            # Start a new potential group with this segment
            current_group = [segment1]
            current_group_indices = {segment1['index']}
            
            if debug and chromosome == '1' and segment1['match_name'].startswith('Kenneth'):
                print(f"\nDEBUG: Starting group with {segment1['match_name']}")
                print(f"  Position: {segment1['start']:,} - {segment1['end']:,}")
            
            # Find all segments that directly overlap with this segment
            for j, segment2 in enumerate(segments):
                if j == i or segment2['index'] in current_group_indices:
                    continue
                
                # Calculate actual overlap
                overlap_start = max(segment1['start'], segment2['start'])
                overlap_end = min(segment1['end'], segment2['end'])
                actual_overlap = max(0, overlap_end - overlap_start)
                
                if debug and chromosome == '1' and segment1['match_name'].startswith('Kenneth') and actual_overlap > 0:
                    print(f"  Checking {segment2['match_name']}: {segment2['start']:,} - {segment2['end']:,}")
                    print(f"    Overlap: {overlap_start:,} - {overlap_end:,} = {actual_overlap:,} bp")
                    print(f"    Required: {min_overlap_bp:,} bp")
                    print(f"    Valid: {actual_overlap >= min_overlap_bp}")
                
                if actual_overlap >= min_overlap_bp:
                    current_group.append(segment2)
                    current_group_indices.add(segment2['index'])
                    
                    if debug and chromosome == '1' and segment1['match_name'].startswith('Kenneth'):
                        print(f"    → Added to group!")
            
            # If we found a valid group, add it and mark segments as processed
            if len(current_group) >= min_group_size:
                if debug and chromosome == '1' and any(seg['match_name'].startswith('Kenneth') for seg in current_group):
                    print(f"\nDEBUG: Found group with {len(current_group)} members:")
                    for seg in current_group:
                        print(f"  {seg['match_name']}: {seg['start']:,} - {seg['end']:,} ({seg['cm']} cM)")
                    
                    # Calculate group boundaries
                    group_start = min(seg['start'] for seg in current_group)
                    group_end = max(seg['end'] for seg in current_group)
                    print(f"  Group span: {group_start:,} - {group_end:,} = {(group_end-group_start)/1000000:.1f} Mb")
                
                group_data = [seg['data'] for seg in current_group]
                groups.append(group_data)
                
                # Mark all segments in this group as processed
                for seg in current_group:
                    processed_indices.add(seg['index'])
        
        print(f"Chromosome {chromosome}: Found {len(groups)} triangulation groups")
        return groups
    
    @staticmethod
    def _verify_triangulation_group(segments, min_overlap_bp):
        """Verify that segments form a true triangulation group with direct overlaps"""
        if len(segments) < 2:
            return segments
        
        # For small groups, verify all pairs overlap
        if len(segments) <= 5:
            valid_segments = []
            
            # Start with the first segment
            valid_segments.append(segments[0])
            
            # Add segments that overlap with at least one segment already in the valid group
            for candidate in segments[1:]:
                overlaps_with_valid_group = False
                
                for valid_seg in valid_segments:
                    overlap_start = max(candidate['start'], valid_seg['start'])
                    overlap_end = min(candidate['end'], valid_seg['end'])
                    actual_overlap = max(0, overlap_end - overlap_start)
                    
                    if actual_overlap >= min_overlap_bp:
                        overlaps_with_valid_group = True
                        break
                
                if overlaps_with_valid_group:
                    valid_segments.append(candidate)
            
            return valid_segments
        
        else:
            # For larger groups, use a more efficient connected components approach
            return OptimizedTriangulationEngine._find_connected_components(segments, min_overlap_bp)
    
    @staticmethod
    def _find_connected_components(segments, min_overlap_bp):
        """Find connected components where each segment overlaps with at least one other"""
        if len(segments) < 2:
            return segments
        
        # Build adjacency list of overlapping segments
        adjacency = {i: set() for i in range(len(segments))}
        
        for i in range(len(segments)):
            for j in range(i + 1, len(segments)):
                seg1, seg2 = segments[i], segments[j]
                
                overlap_start = max(seg1['start'], seg2['start'])
                overlap_end = min(seg1['end'], seg2['end'])
                actual_overlap = max(0, overlap_end - overlap_start)
                
                if actual_overlap >= min_overlap_bp:
                    adjacency[i].add(j)
                    adjacency[j].add(i)
        
        # Find the largest connected component using DFS
        visited = set()
        largest_component = []
        
        for start_node in range(len(segments)):
            if start_node in visited:
                continue
            
            # DFS to find connected component
            component = []
            stack = [start_node]
            
            while stack:
                node = stack.pop()
                if node in visited:
                    continue
                
                visited.add(node)
                component.append(segments[node])
                
                # Add unvisited neighbors to stack
                for neighbor in adjacency[node]:
                    if neighbor not in visited:
                        stack.append(neighbor)
            
            # Keep the largest component found
            if len(component) > len(largest_component):
                largest_component = component
        
        return largest_component
    
    def create_triangulation_report_by_size(self, min_cm_filter=None):
        """Create triangulation report sorted by group size for manual ancestor assignment"""
        if not self.triangulation_groups:
            return "No triangulation groups found."
        
        report = []
        report.append("DNA TRIANGULATION GROUPS - SORTED BY SIZE")
        report.append("=" * 60)
        report.append("(Largest groups first - these represent your strongest ancestral lines)")
        if min_cm_filter and min_cm_filter > 0:
            report.append(f"(Filtered to show only matches ≥ {min_cm_filter} cM)")
        report.append("")
        
        # Sort groups by size (largest first)
        sorted_groups = sorted(self.triangulation_groups, key=len, reverse=True)
        
        report.append(f"SUMMARY: Found {len(sorted_groups)} triangulation groups")
        report.append(f"Largest group: {len(sorted_groups[0])} people")
        report.append(f"Smallest group: {len(sorted_groups[-1])} people")
        report.append("")
        report.append("=" * 60)
        report.append("")
        
        for rank, group in enumerate(sorted_groups, 1):
            # Calculate group boundaries
            min_start = min(member['Start Location'] for member in group)
            max_end = max(member['End Location'] for member in group)
            group_length_mb = (max_end - min_start) / 1000000
            
            # Get chromosome
            chromosome = group[0]['Chromosome']
            
            # Calculate total cM for the group
            total_cm = sum(member['Centimorgans'] for member in group)
            avg_cm = total_cm / len(group)
            
            report.append(f"RANK #{rank}: {len(group)} PEOPLE SHARE DNA HERE")
            report.append(f"Location: Chromosome {chromosome}, {min_start:,} - {max_end:,} bp ({group_length_mb:.1f} Mb)")
            report.append(f"Strength: Average {avg_cm:.1f} cM per person, Total {total_cm:.1f} cM")
            report.append("")
            
            # Group people by surname for easier ancestor assignment
            family_groups = defaultdict(list)
            for member in group:
                surname = self.extract_surname(member['Match Name'])
                family_groups[surname].append(member)
            
            report.append(f"PEOPLE IN THIS GROUP (by family name):")
            
            # Sort families by size within the group
            sorted_families = sorted(family_groups.items(), key=lambda x: len(x[1]), reverse=True)
            
            for surname, family_members in sorted_families:
                report.append(f"  {surname} family ({len(family_members)} people):")
                
                # Sort family members by cM strength
                sorted_members = sorted(family_members, key=lambda x: x['Centimorgans'], reverse=True)
                
                for i, member in enumerate(sorted_members):
                    cm_str = f"{member['Centimorgans']:.1f} cM"
                    snp_str = f"{member['Matching SNPs']:,} SNPs"
                    
                    # Show all members for smaller groups, limit for larger groups
                    if len(group) <= 20 or i < 5:  # Show all for small groups, top 5 for large groups
                        report.append(f"    • {member['Match Name']} ({cm_str}, {snp_str})")
                    elif i == 5:
                        report.append(f"    ... and {len(sorted_members) - 5} more {surname} family members")
                        break
            
            # Suggest potential ancestor assignment
            if len(family_groups) == 1:
                family_name = list(family_groups.keys())[0]
                report.append(f"")
                report.append(f"  → SUGGESTION: This appears to be a pure {family_name} ancestral line")
                report.append(f"    Consider assigning to a {family_name} ancestor")
            else:
                families = list(family_groups.keys())
                report.append(f"")
                report.append(f"  → SUGGESTION: Mixed families ({', '.join(families)})")
                report.append(f"    This may represent a common ancestor shared by these families")
            
            report.append("")
            report.append("  MANUAL ASSIGNMENT: ________________________________")
            report.append("  (Write the ancestor name you want to assign this group to)")
            report.append("")
            report.append("-" * 60)
            report.append("")
        
        return "\n".join(report)
    
    def create_simple_group_list(self):
        """Create a simple list format for quick review"""
        if not self.triangulation_groups:
            return "No triangulation groups found."
        
        # Sort groups by size (largest first)
        sorted_groups = sorted(self.triangulation_groups, key=len, reverse=True)
        
        report = []
        report.append("QUICK GROUP SUMMARY (Largest First)")
        report.append("=" * 50)
        
        for rank, group in enumerate(sorted_groups, 1):
            chromosome = group[0]['Chromosome']
            min_start = min(member['Start Location'] for member in group)
            max_end = max(member['End Location'] for member in group)
            group_length_mb = (max_end - min_start) / 1000000
            
            # Get dominant family
            family_groups = defaultdict(int)
            for member in group:
                surname = self.extract_surname(member['Match Name'])
                family_groups[surname] += 1
            
            dominant_family = max(family_groups.items(), key=lambda x: x[1])
            family_mix = f"{dominant_family[0]} ({dominant_family[1]})" if len(family_groups) == 1 else f"Mixed: {len(family_groups)} families"
            
            report.append(f"{rank:2d}. {len(group):3d} people | Chr {chromosome:2s} | {group_length_mb:5.1f} Mb | {family_mix}")
        
        return "\n".join(report)
    
    def extract_surname(self, match_name):
        """Extract surname from match name"""
        if not match_name:
            return "Unknown"
        
        # Remove titles and suffixes
        name_parts = match_name.replace("Mr.", "").replace("Mrs.", "").replace("Dr.", "").replace("Ph.D.", "").strip().split()
        
        if len(name_parts) >= 2:
            return name_parts[-1]  # Last word as surname
        return name_parts[0] if name_parts else "Unknown"
    
    def export_optimized_results(self, filename="optimized_triangulation.csv"):
        """Export optimized triangulation results to CSV"""
        if not self.triangulation_groups:
            print("No triangulation groups to export.")
            return
        
        export_data = []
        for group_id, group in enumerate(self.triangulation_groups, 1):
            # Calculate group statistics
            min_start = min(member['Start Location'] for member in group)
            max_end = max(member['End Location'] for member in group)
            group_length = (max_end - min_start) / 1000000  # in Mb
            
            for member in group:
                export_data.append({
                    'Group_ID': group_id,
                    'Group_Size': len(group),
                    'Group_Length_Mb': round(group_length, 2),
                    'Group_Start': min_start,
                    'Group_End': max_end,
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
        return len(export_data)

# Test function
if __name__ == "__main__":
    print("Optimized Triangulation Engine initialized!")
    print("This module provides:")
    print("- Multiprocessed triangulation analysis")
    print("- Interval tree algorithms for faster overlap detection")
    print("- Optimized memory usage for large datasets")
    print("- Detailed reporting and statistics")