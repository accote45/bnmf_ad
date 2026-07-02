#!/usr/bin/env python3
"""
GWAS Extractor Script

Downloads GWAS summary statistics files from URLs listed in a Google Sheet.
Only downloads files where the "Downloaded" column is "No", and renames
them according to the "File Name" column after successful download.
"""

import csv
import subprocess
import sys
import os
import urllib.request
import urllib.parse
import tempfile
import gzip
from pathlib import Path


def download_google_sheet_csv(sheet_id, output_file="temp_sheet.csv"):
    """
    Download a Google Sheet as CSV.
    
    Args:
        sheet_id: The Google Sheet ID from the URL
        output_file: Path to save the CSV file
    """
    # Google Sheets CSV export URL format
    csv_url = f"https://docs.google.com/spreadsheets/d/{sheet_id}/export?format=csv&gid=0"
    
    print(f"Downloading Google Sheet as CSV...")
    try:
        urllib.request.urlretrieve(csv_url, output_file)
        print(f"Sheet downloaded to {output_file}")
        return output_file
    except Exception as e:
        print(f"Error downloading Google Sheet: {e}")
        sys.exit(1)


def parse_sheet_and_download(csv_file, download_dir):
    """
    Parse the CSV file and download files where Downloaded column is "No".
    
    Args:
        csv_file: Path to the CSV file
        download_dir: Directory to save downloaded files
    """
    # Create download directory if it doesn't exist
    Path(download_dir).mkdir(parents=True, exist_ok=True)
    
    # Column indices (0-indexed, accounting for header row)
    # Column I = index 8 (Download Link)
    # Column L = index 11 (Downloaded)
    # Column M = index 12 (File Name)
    
    COL_DOWNLOAD_LINK = 8  # Column I
    COL_DOWNLOADED = 11     # Column L
    COL_FILE_NAME = 12      # Column M
    
    downloaded_count = 0
    skipped_count = 0
    error_count = 0
    
    with open(csv_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        
        # Skip header rows (first 2 rows based on the sheet structure)
        next(reader)  # Skip "Market Scan of Cardiometabolic GWAS" row
        next(reader)  # Skip column header row
        
        for row_num, row in enumerate(reader, start=3):
            # Skip empty rows
            if not row or len(row) <= COL_FILE_NAME:
                continue
            
            download_link = row[COL_DOWNLOAD_LINK].strip() if len(row) > COL_DOWNLOAD_LINK else ""
            downloaded_status = row[COL_DOWNLOADED].strip() if len(row) > COL_DOWNLOADED else ""
            file_name = row[COL_FILE_NAME].strip() if len(row) > COL_FILE_NAME else ""
            
            # Check if download link is unavailable
            if download_link.upper() == "UNAVAILABLE":
                print("ftp link is unavailable, please download manually")
                sys.exit(1)
            
            # Skip if no download link
            if not download_link:
                continue
            
            # Only download if status is "No" or empty (not yet downloaded)
            if downloaded_status and downloaded_status.lower() not in ["no", ""]:
                print(f"Row {row_num}: Already downloaded, skipping...")
                skipped_count += 1
                continue
            
            # Skip if no file name specified
            if not file_name:
                print(f"Row {row_num}: No file name specified, skipping...")
                skipped_count += 1
                continue
            
            print(f"\nRow {row_num}: Downloading from {download_link}")
            print(f"  Target filename: {file_name}")
            
            # Create a temporary directory for wget to download with original filename
            temp_dir = tempfile.mkdtemp(dir=download_dir)
            original_cwd = os.getcwd()
            
            try:
                # Change to temp directory and let wget download with default filename
                os.chdir(temp_dir)
                
                # Run wget without -O to preserve original filename
                wget_cmd = [
                    "wget",
                    download_link
                ]
                
                result = subprocess.run(
                    wget_cmd,
                    capture_output=True,
                    text=True,
                    check=True
                )
                
                os.chdir(original_cwd)
                
                # Find the downloaded file (wget will use the filename from URL)
                downloaded_files = os.listdir(temp_dir)
                if not downloaded_files:
                    print(f"  ✗ Error: File was not downloaded")
                    error_count += 1
                    # Clean up temp directory
                    for file in os.listdir(temp_dir):
                        os.remove(os.path.join(temp_dir, file))
                    os.rmdir(temp_dir)
                    continue
                
                # Get the downloaded file path
                downloaded_filename = downloaded_files[0]
                temp_path = os.path.join(temp_dir, downloaded_filename)
                
                # Check if the downloaded file is a gzip file
                is_gzip_file = False
                try:
                    # Try to open as gzip to verify it's a valid gzip file
                    with gzip.open(temp_path, 'rb') as f:
                        f.read(1)  # Try to read to verify it's valid gzip
                    # If we get here, it's a valid gzip file
                    is_gzip_file = True
                except (gzip.BadGzipFile, OSError):
                    # Not a valid gzip file
                    is_gzip_file = False
                
                # Treat as .txt.gz if: valid gzip AND (target filename ends with .txt.gz OR downloaded filename contains "gz")
                should_rename = is_gzip_file and (file_name.endswith('.txt.gz') or 'gz' in downloaded_filename.lower())
                
                if should_rename:
                    # Rename to target filename only if it's .txt.gz
                    final_path = os.path.join(download_dir, file_name)
                    os.rename(temp_path, final_path)
                    os.rmdir(temp_dir)
                    print(f"  ✓ Successfully downloaded and renamed to {final_path}")
                else:
                    # Move file to download_dir with original filename (preserving extension)
                    final_path = os.path.join(download_dir, downloaded_filename)
                    os.rename(temp_path, final_path)
                    os.rmdir(temp_dir)
                    print(f"  ✓ Successfully downloaded to {final_path} (not .txt.gz, keeping original name)")
                
                downloaded_count += 1
                
            except subprocess.CalledProcessError as e:
                # Restore working directory
                os.chdir(original_cwd)
                print(f"  ✗ Error downloading: {e}")
                print(f"  stderr: {e.stderr}")
                error_count += 1
                # Clean up temp directory
                if os.path.exists(temp_dir):
                    for file in os.listdir(temp_dir):
                        os.remove(os.path.join(temp_dir, file))
                    os.rmdir(temp_dir)
            except Exception as e:
                # Restore working directory
                os.chdir(original_cwd)
                print(f"  ✗ Unexpected error: {e}")
                error_count += 1
                # Clean up temp directory
                if os.path.exists(temp_dir):
                    for file in os.listdir(temp_dir):
                        os.remove(os.path.join(temp_dir, file))
                    os.rmdir(temp_dir)
            except Exception as e:
                print(f"  ✗ Unexpected error: {e}")
                error_count += 1
    
    print(f"\n{'='*60}")
    print(f"Download Summary:")
    print(f"  Successfully downloaded: {downloaded_count}")
    print(f"  Skipped: {skipped_count}")
    print(f"  Errors: {error_count}")
    print(f"{'='*60}")


def main():
    """Main function."""
    # Google Sheet ID from the URL
    sheet_id = "10D4g2csnbHck4tZ5BR0HWbf6k8a5hZnNjFbTOsJJuac"
    
    # Default: project root / sumstats/
    project_root = Path(__file__).resolve().parent.parent.parent
    download_dir = project_root / "sumstats"
    if len(sys.argv) > 1:
        download_dir = Path(sys.argv[1])
    
    # Download the Google Sheet as CSV
    csv_file = download_google_sheet_csv(sheet_id)
    
    try:
        # Parse and download files
        parse_sheet_and_download(csv_file, download_dir)
    finally:
        # Clean up temporary CSV file
        if os.path.exists(csv_file):
            os.remove(csv_file)
            print(f"\nCleaned up temporary file: {csv_file}")


if __name__ == "__main__":
    main()
