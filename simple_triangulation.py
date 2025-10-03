import pandas as pd
import numpy as np
from collections import defaultdict

class SimpleTriangulationEngine:
    """Simple, debuggable triangulation engine with clear overlap logic"""
    
    def __init__(self):
        self.triangulation_groups = []
        
    def detect_triangulation_simple(self, data, min_overlap_bp=1000000, min_group_size=2, debug=False):
        """
        Simple triangulation detection with clear logic and debugging
        """
        debug_output = []
        debug_output.append(f"Starting SIMPLE triangulation analysis on {len(data)} segments...")
        
        # Group by chromosome
        chromosome_groups = data.groupby('Chromosome')
        all_groups = []
        
        for chromosome, chr_data in chromosome_groups:
            if len(chr_data) < min_group_size:
                continue
                
            debug_output.append(f"\nProcessing chromosome {chromosome} with {len(chr_data)} segments...")
            
            if debug and chromosome == '1':
                debug_output.append(f"\nDEBUG: All segments on chromosome 1:")
                for idx, row in chr_data.iterrows():
                    debug_output.append(f"  {row['Match Name']}: {row['Start Location']:,} - {row['End Location']:,} bp ({row['Centimorgans']} cM)")
                debug_output.append("")
            
            # Convert to simple list
            segments = []
            for idx, row in chr_data.iterrows():
                segments.append({
                    'name': row['Match Name'],
                    'start': row['Start Location'],
                    'end': row['End Location'],
                    'cm': row['Centimorgans'],
                    'data': row
                })
            
            # Find groups using simple pairwise comparison
            chr_groups, chr_debug = self._find_overlapping_groups_simple(segments, min_overlap_bp, min_group_size, debug and chromosome == '1')
            all_groups.extend(chr_groups)
            if debug and chromosome == '1':
                debug_output.extend(chr_debug)
        
        self.triangulation_groups = all_groups
        self.debug_output = debug_output
        debug_output.append(f"\nFound {len(all_groups)} total triangulation groups")
        return all_groups
    
    def _find_overlapping_groups_simple(self, segments, min_overlap_bp, min_group_size, debug=False):
        """Find overlapping groups using simple, clear logic"""
        groups = []
        used_segments = set()
        debug_output = []
        
        for i, segment1 in enumerate(segments):
            if i in used_segments:
                continue
            
            if debug:
                debug_output.append(f"\nDEBUG: Starting new group with {segment1['name']}")
                debug_output.append(f"  Position: {segment1['start']:,} - {segment1['end']:,} bp")
            
            # Start a new group
            current_group = [segment1]
            current_indices = {i}
            
            # Find all segments that overlap with segment1
            for j, segment2 in enumerate(segments):
                if j == i or j in used_segments:
                    continue
                
                # Calculate overlap
                overlap_start = max(segment1['start'], segment2['start'])
                overlap_end = min(segment1['end'], segment2['end'])
                actual_overlap = overlap_end - overlap_start
                
                if debug:
                    debug_output.append(f"  Checking {segment2['name']}: {segment2['start']:,} - {segment2['end']:,} bp")
                    debug_output.append(f"    Overlap calculation:")
                    debug_output.append(f"      max({segment1['start']:,}, {segment2['start']:,}) = {overlap_start:,}")
                    debug_output.append(f"      min({segment1['end']:,}, {segment2['end']:,}) = {overlap_end:,}")
                    debug_output.append(f"      Overlap length: {overlap_end:,} - {overlap_start:,} = {actual_overlap:,} bp")
                    debug_output.append(f"      Required minimum: {min_overlap_bp:,} bp")
                    debug_output.append(f"      Valid overlap: {actual_overlap >= min_overlap_bp}")
                
                if actual_overlap >= min_overlap_bp:
                    current_group.append(segment2)
                    current_indices.add(j)
                    if debug:
                        debug_output.append(f"    → ADDED to group!")
                else:
                    if debug:
                        debug_output.append(f"    → NOT added (insufficient overlap)")
            
            # If we found a valid group, save it
            if len(current_group) >= min_group_size:
                if debug:
                    debug_output.append(f"\nDEBUG: Group completed with {len(current_group)} members:")
                    for seg in current_group:
                        debug_output.append(f"  {seg['name']}: {seg['start']:,} - {seg['end']:,} bp ({seg['cm']} cM)")
                    
                    group_start = min(seg['start'] for seg in current_group)
                    group_end = max(seg['end'] for seg in current_group)
                    group_span = (group_end - group_start) / 1000000
                    debug_output.append(f"  Group span: {group_start:,} - {group_end:,} bp ({group_span:.1f} Mb)")
                
                # Convert back to original data format
                group_data = [seg['data'] for seg in current_group]
                groups.append(group_data)
                
                # Mark segments as used
                used_segments.update(current_indices)
            elif debug:
                debug_output.append(f"    Group too small ({len(current_group)} < {min_group_size}), discarded")
        
        return groups, debug_output
    
    def create_triangulation_report_by_size(self, min_cm_filter=None):
        """Create triangulation report sorted by group size"""
        if not self.triangulation_groups:
            return "No triangulation groups found."
        
        report = []
        report.append("SIMPLE DNA TRIANGULATION GROUPS - SORTED BY SIZE")
        report.append("=" * 60)
        report.append("(Largest groups first - using simple overlap detection)")
        if min_cm_filter and min_cm_filter > 0:
            report.append(f"(Filtered to show only matches ≥ {min_cm_filter} cM)")
        report.append("")
        
        # Sort groups by size (largest first)
        sorted_groups = sorted(self.triangulation_groups, key=len, reverse=True)
        
        report.append(f"SUMMARY: Found {len(sorted_groups)} triangulation groups")
        if sorted_groups:
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
                    
                    # Show details for each member including position
                    start_mb = member['Start Location'] / 1000000
                    end_mb = member['End Location'] / 1000000
                    pos_str = f"{start_mb:.1f}-{end_mb:.1f}Mb"
                    
                    if len(group) <= 20 or i < 5:
                        report.append(f"    • {member['Match Name']} ({cm_str}, {pos_str})")
                    elif i == 5:
                        report.append(f"    ... and {len(sorted_members) - 5} more {surname} family members")
                        break
            
            report.append("")
            report.append("-" * 60)
            report.append("")
        
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

# Test function
if __name__ == "__main__":
    print("Simple Triangulation Engine initialized!")
    print("This engine uses:")
    print("- Clear, step-by-step overlap detection")
    print("- Detailed debugging output")
    print("- Simple pairwise comparison logic")
    print("- No complex algorithms that could have bugs")