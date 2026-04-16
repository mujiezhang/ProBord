#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
ProBord CLI: Command-line interface with subcommands.

Usage:
    probord run -hf <host_fasta> -vf <virus_info> -wd <output_dir> -db <blastn_db> [-cv <checkv_db>]
    probord prepare_db <genus_name> [bacteria|archaea] [-t <threads>] [-o <output_dir>]
"""

import argparse
import subprocess
import os
import sys
import glob
import shutil


__version__ = "1.0"

def setup_run_parser(subparsers):
    """Set up the 'run' subcommand parser."""
    run_parser = subparsers.add_parser('run', 
        help='Run the ProBord pipeline for provirus border delimitation',
        description='Run the ProBord pipeline for provirus border delimitation')
    
    required = run_parser.add_argument_group('Required arguments')
    required.add_argument('-hf','--host_fasta', type=str, required=True, 
        help='Host genome/contig file containing provirus (FASTA format)', metavar='<path>')
    required.add_argument('-vf', '--virus_information', type=str, required=True, 
        help='A tab-delimited file with columns: viral_name, host_contig, start, end', metavar='<path>')
    required.add_argument('-wd', '--working_path', type=str, required=True, 
        help='Path to the output directory', metavar='<path>')
    required.add_argument('-db', '--blastn_db', type=str, help='Path to the BLASTn database for attB detection', metavar='<path>')
    
    optional = run_parser.add_argument_group('Optional arguments')
    optional.add_argument('-cv', '--checkv_db', type=str, help='Path to the CheckV database', metavar='<path>')
    optional.add_argument('-s', '--score', type=int, default=20, help='Cutoff for attB score (default: 20)', metavar='<int>')
    optional.add_argument('-t','--threads', type=int, default=8, help='Number of threads to use (default: 8)', metavar='<int>')
    optional.add_argument('-k','--keep-temp', action='store_true', help='Keep temporary files after the run')


def setup_prepare_db_parser(subparsers):
    """Set up the 'prepare_db' subcommand parser."""
    db_parser = subparsers.add_parser('prepare_db',
        help='Download NCBI genomes and build BLAST database for a genus',
        description='Download NCBI genomes and build BLAST database for a genus')
    
    db_parser.add_argument('genus', type=str, help='Genus name (case-sensitive)')
    db_parser.add_argument('group', type=str, nargs='?', default='bacteria', 
        choices=['bacteria', 'archaea'], help='Domain: bacteria or archaea (default: bacteria)')
    db_parser.add_argument('-t', '--threads', type=int, default=4, 
        help='Number of download threads (default: 4)', metavar='<int>')
    db_parser.add_argument('-o', '--output', type=str, default=None,
        help='Custom output directory (default: use genus name)', metavar='<path>')


def cmd_run(args):
    """Execute the 'run' subcommand by invoking probord.main()."""
    # Reconstruct sys.argv for probord.py's argument parser
    run_argv = [
        '-hf', args.host_fasta,
        '-vf', args.virus_information,
        '-wd', args.working_path,
    ]
    if args.blastn_db:
        run_argv.extend(['-db', args.blastn_db])
    if args.checkv_db:
        run_argv.extend(['-cv', args.checkv_db])
    if args.score != 20:
        run_argv.extend(['-s', str(args.score)])
    if args.threads != 8:
        run_argv.extend(['-t', str(args.threads)])
    if args.keep_temp:
        run_argv.append('-k')
    
    # Override sys.argv so probord.py's argument parser picks up the arguments
    sys.argv = ['probord'] + run_argv
    
    from probord.probord import main
    main()


def cmd_prepare_db(args):
    """Execute the 'prepare_db' subcommand."""
    output_dir = args.output if args.output else args.genus
    
    print("=" * 54)
    print(f"Building BLAST database for: {args.genus} ({args.group})")
    print(f"Output directory: {output_dir}")
    print(f"Download threads: {args.threads}")
    print("=" * 54)
    
    os.makedirs(output_dir, exist_ok=True)
    original_dir = os.getcwd()
    
    try:
        os.chdir(output_dir)
        
        # Step 1: Download genomes
        print(f"[1/4] Downloading {args.genus} genomes from {args.group}...")
        cmd = [
            "ncbi-genome-download",
            "-F", "fasta",
            "--genera", args.genus,
            "-o", args.genus,
            "--flat-output",
            "-P",
            "-p", str(args.threads),
            args.group
        ]
        subprocess.run(cmd, check=True)
        
        # Check download result
        if not os.path.exists(args.genus) or not os.listdir(args.genus):
            print("Error: No genomes downloaded. Possible reasons:")
            print("  - Genus name misspelled (case-sensitive)")
            print(f"  - No genomes available for {args.genus} in {args.group}")
            print("  - NCBI API connection issue")
            sys.exit(1)
        
        # Step 2: Process genomes
        print("[2/4] Processing downloaded genomes...")
        for gz_file in glob.glob(os.path.join(args.genus, "*.gz")):
            subprocess.run(["gunzip", "-f", gz_file], check=False)
        
        fna_files = glob.glob(os.path.join(args.genus, "*.fna"))
        if not fna_files:
            print("Error: No FASTA files found after decompression.")
            sys.exit(1)
        
        combined_fna = f"{args.genus}.fna"
        with open(combined_fna, 'w') as outf:
            for fna in fna_files:
                with open(fna) as inf:
                    shutil.copyfileobj(inf, outf)
        
        # Step 3: Build BLAST database
        print("[3/4] Building BLAST database (this may take several minutes)...")
        cmd = [
            "makeblastdb",
            "-in", combined_fna,
            "-dbtype", "nucl",
            "-out", args.genus,
            "-title", f"{args.genus}_DB"
        ]
        subprocess.run(cmd, check=True)
        
        # Step 4: Clean up
        print("[4/4] Cleaning temporary files...")
        shutil.rmtree(args.genus, ignore_errors=True)
        os.remove(combined_fna)
        
        print("=" * 54)
        print(f"Successfully created BLAST database for {args.genus}")
        print(f"Database path: {os.path.join(output_dir, args.genus)}")
        print(f"Use this path as the -db argument: {os.path.join(output_dir, args.genus, args.genus)}")
        
    except subprocess.CalledProcessError as e:
        print(f"Error: Command failed: {e}")
        sys.exit(1)
    finally:
        os.chdir(original_dir)


def main():
    """Main entry point for ProBord CLI."""
    parser = argparse.ArgumentParser(
        prog='probord',
        description=f'ProBord: Provirus Border Delimiter v{__version__}',
    )
    parser.add_argument('-v', '--version', action='version', version=f'ProBord v{__version__}')
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    subparsers.required = True
    
    setup_run_parser(subparsers)
    setup_prepare_db_parser(subparsers)
    
    args = parser.parse_args()
    
    if args.command == 'run':
        cmd_run(args)
    elif args.command == 'prepare_db':
        cmd_prepare_db(args)


if __name__ == "__main__":
    main()