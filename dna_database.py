import sqlite3
import pandas as pd
import os
from datetime import datetime
from typing import List, Dict, Optional, Tuple

class DNADatabase:
    """Handles all database operations for DNA segment data"""
    
    def __init__(self, db_path="dna_analysis.db"):
        self.db_path = db_path
        self.init_database()
    
    def init_database(self):
        """Initialize database with required tables"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        # DNA segments table
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS dna_segments (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                match_name TEXT NOT NULL,
                chromosome TEXT NOT NULL,
                start_location INTEGER NOT NULL,
                end_location INTEGER NOT NULL,
                centimorgans REAL NOT NULL,
                matching_snps INTEGER NOT NULL,
                source_file TEXT,
                import_date TEXT NOT NULL,
                assigned_family TEXT,
                notes TEXT
            )
        ''')
        
        # Family/surname assignments table
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS family_assignments (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                surname TEXT NOT NULL UNIQUE,
                color_code TEXT,
                ancestral_line TEXT,
                notes TEXT,
                created_date TEXT NOT NULL
            )
        ''')
        
        # Triangulation groups table
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS triangulation_groups (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                group_name TEXT NOT NULL,
                chromosome TEXT NOT NULL,
                start_location INTEGER NOT NULL,
                end_location INTEGER NOT NULL,
                member_count INTEGER NOT NULL,
                confidence_score REAL,
                created_date TEXT NOT NULL,
                notes TEXT
            )
        ''')
        
        # Triangulation group members table
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS triangulation_members (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                group_id INTEGER NOT NULL,
                segment_id INTEGER NOT NULL,
                FOREIGN KEY (group_id) REFERENCES triangulation_groups (id),
                FOREIGN KEY (segment_id) REFERENCES dna_segments (id)
            )
        ''')
        
        # Import history table
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS import_history (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                filename TEXT NOT NULL,
                file_path TEXT NOT NULL,
                segments_imported INTEGER NOT NULL,
                import_date TEXT NOT NULL,
                file_hash TEXT
            )
        ''')
        
        conn.commit()
        conn.close()
        print(f"Database initialized: {self.db_path}")
    
    def import_csv_data(self, file_path: str, source_name: str = None) -> Tuple[int, int]:
        """Import DNA segment data from CSV file"""
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")
        
        # Read CSV
        df = pd.read_csv(file_path)
        
        # Validate required columns
        required_cols = ['Match Name', 'Chromosome', 'Start Location', 'End Location', 'Centimorgans', 'Matching SNPs']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")
        
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        import_date = datetime.now().isoformat()
        source_name = source_name or os.path.basename(file_path)
        
        # Check for duplicates
        imported_count = 0
        duplicate_count = 0
        
        for _, row in df.iterrows():
            # Check if segment already exists
            cursor.execute('''
                SELECT COUNT(*) FROM dna_segments 
                WHERE match_name = ? AND chromosome = ? AND start_location = ? AND end_location = ?
            ''', (row['Match Name'], str(row['Chromosome']), int(row['Start Location']), int(row['End Location'])))
            
            if cursor.fetchone()[0] == 0:
                # Extract surname for auto-assignment
                surname = self.extract_surname(row['Match Name'])
                
                cursor.execute('''
                    INSERT INTO dna_segments (
                        match_name, chromosome, start_location, end_location, 
                        centimorgans, matching_snps, source_file, import_date, assigned_family
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                ''', (
                    row['Match Name'], str(row['Chromosome']), int(row['Start Location']), 
                    int(row['End Location']), float(row['Centimorgans']), 
                    int(row['Matching SNPs']), source_name, import_date, surname
                ))
                imported_count += 1
            else:
                duplicate_count += 1
        
        # Record import history
        cursor.execute('''
            INSERT INTO import_history (filename, file_path, segments_imported, import_date)
            VALUES (?, ?, ?, ?)
        ''', (os.path.basename(file_path), file_path, imported_count, import_date))
        
        conn.commit()
        conn.close()
        
        print(f"Imported {imported_count} new segments, skipped {duplicate_count} duplicates")
        return imported_count, duplicate_count
    
    def extract_surname(self, match_name: str) -> str:
        """Extract surname from match name"""
        if not match_name:
            return "Unknown"
        
        # Remove titles and suffixes
        name_parts = match_name.replace("Mr.", "").replace("Mrs.", "").replace("Dr.", "").replace("Ph.D.", "").strip().split()
        
        if len(name_parts) >= 2:
            return name_parts[-1]  # Last word as surname
        return name_parts[0] if name_parts else "Unknown"
    
    def get_all_segments(self) -> pd.DataFrame:
        """Get all DNA segments as DataFrame"""
        conn = sqlite3.connect(self.db_path)
        df = pd.read_sql_query('''
            SELECT * FROM dna_segments ORDER BY chromosome, start_location
        ''', conn)
        conn.close()
        return df
    
    def get_segments_by_family(self, family_name: str) -> pd.DataFrame:
        """Get segments for specific family"""
        conn = sqlite3.connect(self.db_path)
        df = pd.read_sql_query('''
            SELECT * FROM dna_segments WHERE assigned_family = ? 
            ORDER BY chromosome, start_location
        ''', conn, params=[family_name])
        conn.close()
        return df
    
    def get_all_families(self) -> List[str]:
        """Get list of all assigned families"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        cursor.execute('SELECT DISTINCT assigned_family FROM dna_segments WHERE assigned_family IS NOT NULL')
        families = [row[0] for row in cursor.fetchall()]
        conn.close()
        return sorted(families)
    
    def update_segment_family(self, segment_id: int, family_name: str):
        """Update family assignment for a segment"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        cursor.execute('''
            UPDATE dna_segments SET assigned_family = ? WHERE id = ?
        ''', (family_name, segment_id))
        conn.commit()
        conn.close()
    
    def delete_segment(self, segment_id: int):
        """Delete a DNA segment"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        cursor.execute('DELETE FROM dna_segments WHERE id = ?', (segment_id,))
        conn.commit()
        conn.close()
    
    def get_database_stats(self) -> Dict:
        """Get database statistics"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        stats = {}
        
        # Total segments
        cursor.execute('SELECT COUNT(*) FROM dna_segments')
        stats['total_segments'] = cursor.fetchone()[0]
        
        # Unique matches
        cursor.execute('SELECT COUNT(DISTINCT match_name) FROM dna_segments')
        stats['unique_matches'] = cursor.fetchone()[0]
        
        # Chromosomes covered
        cursor.execute('SELECT COUNT(DISTINCT chromosome) FROM dna_segments')
        stats['chromosomes_covered'] = cursor.fetchone()[0]
        
        # Families assigned
        cursor.execute('SELECT COUNT(DISTINCT assigned_family) FROM dna_segments WHERE assigned_family IS NOT NULL')
        stats['families_assigned'] = cursor.fetchone()[0]
        
        # Import history
        cursor.execute('SELECT COUNT(*) FROM import_history')
        stats['import_sessions'] = cursor.fetchone()[0]
        
        # Total cM
        cursor.execute('SELECT SUM(centimorgans) FROM dna_segments')
        result = cursor.fetchone()[0]
        stats['total_centimorgans'] = result if result else 0
        
        conn.close()
        return stats
    
    def get_import_history(self) -> pd.DataFrame:
        """Get import history"""
        conn = sqlite3.connect(self.db_path)
        df = pd.read_sql_query('''
            SELECT filename, segments_imported, import_date 
            FROM import_history 
            ORDER BY import_date DESC
        ''', conn)
        conn.close()
        return df
    
    def search_segments(self, search_term: str, search_type: str = "match_name") -> pd.DataFrame:
        """Search segments by various criteria"""
        conn = sqlite3.connect(self.db_path)
        
        if search_type == "match_name":
            query = "SELECT * FROM dna_segments WHERE match_name LIKE ? ORDER BY chromosome, start_location"
            params = [f"%{search_term}%"]
        elif search_type == "family":
            query = "SELECT * FROM dna_segments WHERE assigned_family LIKE ? ORDER BY chromosome, start_location"
            params = [f"%{search_term}%"]
        elif search_type == "chromosome":
            query = "SELECT * FROM dna_segments WHERE chromosome = ? ORDER BY start_location"
            params = [search_term]
        else:
            # Default to match name
            query = "SELECT * FROM dna_segments WHERE match_name LIKE ? ORDER BY chromosome, start_location"
            params = [f"%{search_term}%"]
        
        df = pd.read_sql_query(query, conn, params=params)
        conn.close()
        return df

# Test the database
if __name__ == "__main__":
    db = DNADatabase()
    stats = db.get_database_stats()
    print("Database Statistics:")
    for key, value in stats.items():
        print(f"  {key}: {value}")