#!/usr/bin/env python3
"""
MAF file indexing script with error handling and diagnostics
Indexes MAF files and reports any issues found
Handles gzipped MAF files automatically
"""

import argparse
import os
import sys
import gzip
import shutil
from Bio import AlignIO

def is_gzipped(filepath):
    """Check if a file is gzipped by reading magic bytes"""
    try:
        with open(filepath, 'rb') as f:
            magic = f.read(2)
            return magic == b'\x1f\x8b'  # Gzip magic number
    except:
        return False

def decompress_maf_file(gzipped_file, output_file=None):
    """Decompress a gzipped MAF file"""
    if output_file is None:
        # Decompress in place (remove .gz extension)
        output_file = gzipped_file.rstrip('.gz')
    
    try:
        with gzip.open(gzipped_file, 'rb') as f_in:
            with open(output_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        return output_file
    except Exception as e:
        raise Exception(f"Failed to decompress {gzipped_file}: {str(e)}")

def get_maf_file_path(maf_file):
    """Get the actual MAF file path, decompressing if necessary"""
    if not os.path.exists(maf_file):
        return None, "File does not exist"
    
    original_file = maf_file
    renamed = False
    
    # Check if file is gzipped (by magic bytes, not extension)
    if is_gzipped(maf_file):
        # If file doesn't end with .gz, it has wrong extension - rename it
        if not maf_file.endswith('.gz'):
            # Rename to add .gz extension
            renamed_file = maf_file + '.gz'
            try:
                os.rename(maf_file, renamed_file)
                maf_file = renamed_file  # Update to use renamed file
                renamed = True
            except Exception as e:
                return None, f"rename_error: Could not rename {original_file} to {renamed_file}: {str(e)}"
        
        # Now handle normal .gz file
        decompressed_path = maf_file[:-3]  # Remove .gz extension
        if os.path.exists(decompressed_path):
            # Check if decompressed is newer than compressed
            if os.path.getmtime(decompressed_path) >= os.path.getmtime(maf_file):
                status = "decompressed_exists"
                if renamed:
                    status = "renamed_and_decompressed_exists"
                return decompressed_path, status
        
        # Need to decompress
        try:
            decompressed_path = decompress_maf_file(maf_file)
            status = "decompressed"
            if renamed:
                status = "renamed_and_decompressed"
            return decompressed_path, status
        except Exception as e:
            return None, f"decompression_error: {str(e)}"
    else:
        return maf_file, "not_compressed"

def diagnose_maf_file(maf_file):
    """Diagnose issues with a MAF file"""
    issues = []
    info = {}
    
    # Check if file exists
    if not os.path.exists(maf_file):
        return {"status": "error", "error": "File does not exist", "issues": []}
    
    # Check if gzipped
    is_gz = is_gzipped(maf_file)
    info["is_gzipped"] = is_gz
    
    # Get file size
    file_size = os.path.getsize(maf_file)
    info["size_bytes"] = file_size
    info["size_mb"] = round(file_size / (1024 * 1024), 2)
    
    if file_size == 0:
        issues.append("File is empty (0 bytes)")
        return {"status": "error", "error": "Empty file", "issues": issues, "info": info}
    
    # If gzipped, try to get decompressed size estimate
    if is_gz:
        try:
            with gzip.open(maf_file, 'rb') as f:
                # Read a bit to verify it's valid gzip
                f.read(1000)
            issues.append("File is gzipped and will be decompressed before indexing")
        except Exception as e:
            issues.append(f"File appears to be gzipped but cannot be read: {str(e)}")
            return {"status": "error", "error": "Invalid gzip file", "issues": issues, "info": info}
    
    # Check if file is too small (likely incomplete)
    if file_size < 1024:  # Less than 1 KB
        issues.append("File is suspiciously small (< 1 KB), may be incomplete")
    
    return {"status": "ok" if len(issues) == 0 else "warning", "issues": issues, "info": info}

def index_maf_file(maf_file, idx_file, target_seqname):
    """Index a single MAF file with error handling"""
    diagnosis = diagnose_maf_file(maf_file)
    
    if diagnosis["status"] == "error":
        return {
            "status": "error",
            "chromosome": os.path.basename(maf_file).replace(".maf", "").replace(".gz", ""),
            "error": diagnosis.get("error", "Unknown error"),
            "issues": diagnosis.get("issues", []),
            "info": diagnosis.get("info", {})
        }
    
    # Get actual MAF file (decompress if needed, rename if incorrectly named)
    actual_maf_file, file_status = get_maf_file_path(maf_file)
    if actual_maf_file is None:
        return {
            "status": "error",
            "chromosome": os.path.basename(maf_file).replace(".maf", "").replace(".gz", ""),
            "error": "Could not prepare MAF file",
            "error_detail": file_status,
            "issues": diagnosis.get("issues", []) + [f"File preparation error: {file_status}"],
            "info": diagnosis.get("info", {})
        }
    
    # Add status info to issues
    if file_status == "decompressed":
        diagnosis["issues"].append("File was decompressed from gzip format")
    elif "renamed" in file_status:
        diagnosis["issues"].append("File was renamed to add .gz extension (was incorrectly named)")
        if "decompressed" in file_status:
            diagnosis["issues"].append("File was then decompressed")
    
    # Check if index already exists and is newer than MAF file
    if os.path.exists(idx_file):
        maf_mtime = os.path.getmtime(actual_maf_file)
        idx_mtime = os.path.getmtime(idx_file)
        if idx_mtime > maf_mtime:
            return {
                "status": "skipped",
                "chromosome": os.path.basename(maf_file).replace(".maf", "").replace(".gz", ""),
                "message": "Index already exists and is up to date",
                "info": diagnosis.get("info", {})
            }
    
    # Try to create index
    try:
        idx = AlignIO.MafIO.MafIndex(idx_file, actual_maf_file, target_seqname)
        return {
            "status": "success",
            "chromosome": os.path.basename(maf_file).replace(".maf", "").replace(".gz", ""),
            "message": "Index created successfully",
            "issues": diagnosis.get("issues", []),
            "info": diagnosis.get("info", {})
        }
    except UnicodeDecodeError as e:
        return {
            "status": "error",
            "chromosome": os.path.basename(maf_file).replace(".maf", "").replace(".gz", ""),
            "error": "UnicodeDecodeError",
            "error_detail": str(e),
            "issues": diagnosis.get("issues", []) + [f"Encoding error: {str(e)}"],
            "info": diagnosis.get("info", {}),
            "recommendation": "File may be corrupted or incomplete. Try re-downloading."
        }
    except Exception as e:
        return {
            "status": "error",
            "chromosome": os.path.basename(maf_file).replace(".maf", "").replace(".gz", ""),
            "error": type(e).__name__,
            "error_detail": str(e),
            "issues": diagnosis.get("issues", []) + [f"Indexing error: {str(e)}"],
            "info": diagnosis.get("info", {}),
            "recommendation": "Check file integrity and format."
        }

def main():
    parser = argparse.ArgumentParser(description='Index MAF files with diagnostics')
    parser.add_argument('--maf-dir', required=True, help='Directory containing MAF files')
    parser.add_argument('--chromosomes', nargs='+', help='Specific chromosomes to index (default: all found)')
    parser.add_argument('--target-seqname', default='hg38', help='Target sequence name for indexing')
    
    args = parser.parse_args()
    
    if not os.path.isdir(args.maf_dir):
        print("Error: MAF directory not found: %s" % args.maf_dir)
        sys.exit(1)
    
    # Find MAF files (including .maf.gz)
    if args.chromosomes:
        maf_files = []
        for chr in args.chromosomes:
            # Try both .maf and .maf.gz
            maf_path = os.path.join(args.maf_dir, "%s.maf" % chr)
            maf_gz_path = os.path.join(args.maf_dir, "%s.maf.gz" % chr)
            if os.path.exists(maf_gz_path):
                maf_files.append(maf_gz_path)
            elif os.path.exists(maf_path):
                maf_files.append(maf_path)
    else:
        all_files = os.listdir(args.maf_dir)
        maf_files = []
        for f in all_files:
            if f.endswith('.maf') or f.endswith('.maf.gz'):
                maf_files.append(os.path.join(args.maf_dir, f))
    
    if len(maf_files) == 0:
        print("Error: No MAF files found in %s" % args.maf_dir)
        sys.exit(1)
    
    maf_files.sort()
    
    results = {
        "success": [],
        "error": [],
        "skipped": [],
        "warning": []
    }
    
    print("Indexing MAF files...")
    print("=" * 60)
    
    for maf_file in maf_files:
        chr_name = os.path.basename(maf_file).replace(".maf", "").replace(".gz", "")
        idx_file = os.path.join(args.maf_dir, "%s.mafindex" % chr_name)
        target_seqname = "%s.%s" % (args.target_seqname, chr_name)
        
        print("\nProcessing %s..." % chr_name)
        result = index_maf_file(maf_file, idx_file, target_seqname)
        
        if result["status"] == "success":
            print("  ✓ Success")
            results["success"].append(result)
        elif result["status"] == "skipped":
            print("  ⊘ Skipped (index already exists)")
            results["skipped"].append(result)
        elif result["status"] == "error":
            print("  ✗ Error: %s" % result.get("error", "Unknown"))
            if "error_detail" in result:
                print("    Detail: %s" % result["error_detail"])
            if "recommendation" in result:
                print("    Recommendation: %s" % result["recommendation"])
            results["error"].append(result)
        
        if result.get("issues"):
            print("  Warnings/Info:")
            for issue in result["issues"]:
                print("    - %s" % issue)
        
        # Show rename status if file was renamed
        if "renamed" in str(result.get("issues", [])):
            print("  Note: File had incorrect extension and was renamed to .gz")
        
        if "info" in result and "size_mb" in result["info"]:
            print("  File size: %.2f MB" % result["info"]["size_mb"])
            if result["info"].get("is_gzipped"):
                print("  (File is gzipped)")
    
    # Summary
    print("\n" + "=" * 60)
    print("Indexing Summary:")
    print("  Success: %d" % len(results["success"]))
    print("  Skipped: %d" % len(results["skipped"]))
    print("  Errors: %d" % len(results["error"]))
    
    if len(results["error"]) > 0:
        print("\nFiles with errors:")
        for err in results["error"]:
            print("  - %s: %s" % (err["chromosome"], err.get("error", "Unknown")))
        sys.exit(1)
    
    if len(results["success"]) == 0 and len(results["skipped"]) == 0:
        print("\nNo files were successfully indexed!")
        sys.exit(1)
    
    print("\nAll MAF files indexed successfully!")
    sys.exit(0)

if __name__ == "__main__":
    main()
