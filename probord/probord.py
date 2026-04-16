#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
ProBord: A tool for precise delimitation of provirus borders in host genomes.
"""

import os
import re
import sys
import time
import shutil
import argparse
import subprocess
import resource
from Bio import SeqIO
from collections import defaultdict

# --- Constants ---
__version__ = "1.0"
PROG_NAME = "ProBord"
LOG_FILE_NAME = "probord.log"
QUALIFIED_PROVIRUS_FNA = "qualified_provirus.fna"
EXTENDED_5KB_FNA = "extended_5kb_after_checkv.fna"
FINAL_CANDIDATE_ATTB_FNA = "final_merged_candidate_attB.fna"
ATT_PREDICTION_TSV = "att_prediction.tsv"
DEFAULT_SCORE_CUTOFF = 20
DEFAULT_THREADS = 8
DEFAULT_DISTANCE_CUT = 100
TOTAL_STEPS = 3 # Total number of main steps in the pipeline

# --- Argument Parsing ---
def setup_arg_parser():
    """Set up the command-line argument parser."""
    custom_usage = f"%(prog)s -hf <host_fasta> -vf <virus_info> -wd <output_dir> -db <blastn_db> [-cv <checkv_db>]"

    class CustomHelpFormatter(argparse.RawTextHelpFormatter):
        """Custom help formatter to add a title to the help message."""
        def format_help(self):
            # Avoid adding title to version message
            if any(arg in sys.argv for arg in ['-v', '--version']):
                return super().format_help()
            return f"{PROG_NAME}: Provirus Border Delimiter\n\n{super().format_help()}"

    parser = argparse.ArgumentParser(
        description='',
        usage=custom_usage,
        formatter_class=CustomHelpFormatter,
        add_help=False # Custom help handling
    )

    # Required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-hf','--host_fasta', type=str, required=True, help='Host genome/contig file containing provirus (FASTA format)', metavar='<path>')
    required.add_argument('-vf', '--virus_information', type=str, required=True, help='A tab-delimited file with columns: viral_name, host_contig, start, end', metavar='<path>')
    required.add_argument('-wd', '--working_path', type=str, required=True, help='Path to the output directory', metavar='<path>')
    required.add_argument('-db', '--blastn_db', type=str, help='Path to the BLASTn database for attB detection', metavar='<path>')
    # Optional arguments
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('-cv', '--checkv_db', type=str, help='Path to the CheckV database', metavar='<path>')
    #optional.add_argument('-db', '--blastn_db', type=str, help='Path to the BLASTn database for attB detection', metavar='<path>')
    optional.add_argument('-s', '--score', type=int, default=DEFAULT_SCORE_CUTOFF, help=f'Cutoff for attB score (default: {DEFAULT_SCORE_CUTOFF})', metavar='<int>')
    optional.add_argument('-t','--threads', type=int, default=DEFAULT_THREADS, help=f'Number of threads to use (default: {DEFAULT_THREADS})', metavar='<int>')
    optional.add_argument('-k','--keep-temp', action='store_true', help='Keep temporary files after the run')
    
    # Information arguments
    info = parser.add_argument_group('Information')
    info.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    info.add_argument('-v','--version', action='version', version=f'{PROG_NAME} v{__version__}')
    
    return parser.parse_args()

args = setup_arg_parser()

# --- Global Variables & Setup ---
working_path = args.working_path
host_fasta_path = args.host_fasta
blast_db_path = args.blastn_db
virus_info_path = args.virus_information

# The host_dict will be loaded inside the main function to ensure
# all functions are defined before heavy I/O operations.
log_file_path = os.path.join(working_path, LOG_FILE_NAME)

# --- Logging Setup ---
os.makedirs(working_path, exist_ok=True)

def log_message(message: str):
    """
    Prints a message to the console and writes it to the log file.
    
    Args:
        message (str): The message to log.
    """
    ts = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    line = f"[{ts}] {message}"
    print(line, flush=True)
    with open(log_file_path, 'a') as lf:
        lf.write(line + "\n")

def log_step(step_idx: int = None, total_steps: int = None, custom_msg: str = None):
    """
    A decorator to log the start of a pipeline step.

    Args:
        step_idx (int, optional): The index of the current step.
        total_steps (int, optional): The total number of steps.
        custom_msg (str, optional): A custom message to display.
    """
    def decorator(func):
        def wrapper(*args, **kwargs):
            prefix = f"[{step_idx}/{total_steps}] " if step_idx and total_steps else ""
            msg = custom_msg if custom_msg else f"Running {func.__name__}"
            log_message(f"{prefix}{msg}")
            return func(*args, **kwargs)
        return wrapper
    return decorator

def log_program_stats():
    """Logs the total runtime and peak memory usage of the program."""
    program_end_time = time.time()
    total_time = program_end_time - program_start_time
    
    # Get peak memory usage for both the main script and its children (e.g., BLASTn).
    # This provides a more accurate measure of the workflow's total memory footprint.
    # resource.getrusage reports in kilobytes on Linux and bytes on macOS.
    
    max_memory_mb = 0
    try:
        # Get peak memory of the main Python process
        self_mem_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        
        # Get the peak memory across all completed children processes
        children_mem_usage = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss

        if sys.platform == "darwin":
            # Convert bytes to MB for macOS
            self_mem_mb = self_mem_usage / 1024 / 1024
            children_mem_mb = children_mem_usage / 1024 / 1024
        else: # Primarily for Linux
            # Convert kilobytes to MB for Linux
            self_mem_mb = self_mem_usage / 1024
            children_mem_mb = children_mem_usage / 1024

        # The true peak memory is the max of the parent's peak and the overall children's peak.
        max_memory_mb = max(self_mem_mb, children_mem_mb)

    except Exception as e:
        log_message(f"Warning: Could not determine peak memory usage. Reason: {e}")

    # Convert MB to GB for reporting
    max_memory_gb = max_memory_mb / 1024

    # Format time for readability
    hours, remainder = divmod(total_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    
    if hours > 0:
        time_str = f"{int(hours)}h {int(minutes)}m {seconds:.2f}s"
    elif minutes > 0:
        time_str = f"{int(minutes)}m {seconds:.2f}s"
    else:
        time_str = f"{seconds:.2f}s"
    
    success_line = "✨🌟" * 6 + " SUCCESS " + "🌟✨" * 6
    completion_info = f"{PROG_NAME} completed - Total time: {time_str}, Max memory usage: {max_memory_gb:.2f} GB"
    
    print(f"\n{success_line}", flush=True)
    print(completion_info, flush=True)
    with open(log_file_path, 'a') as lf:
        lf.write(f"\n{success_line}\n")
        lf.write(completion_info + "\n")

# --- Core Functions ---

@log_step(1, TOTAL_STEPS, custom_msg="Preprocessing virus information and running CheckV")
def preprocessing_virus_info(host_dict, distance_cut=DEFAULT_DISTANCE_CUT):
    """
    Parses the initial virus information, extracts sequences, and runs CheckV
    to remove host contamination.

    Args:
        distance_cut (int): Minimum required distance from a provirus end to a contig end.

    Returns:
        tuple: A tuple containing:
            - dict: A dictionary of qualified proviruses.
            - dict: A dictionary containing the original input information.
    """
    preprocess_dir = os.path.join(working_path, "preprocess")
    if os.path.exists(preprocess_dir):
        log_message(f"Error: Preprocessing directory '{preprocess_dir}' already exists. Please remove or rename it.")
        sys.exit(1)
    os.makedirs(preprocess_dir)
    
    original_info = {}
    original_virus_fna = os.path.join(preprocess_dir, "original_virus.fna")
    
    # Read the input virus information and extract sequences
    with open(virus_info_path, 'r') as vi, open(original_virus_fna, "w") as f:
        for i in vi:
            # Skip header or empty lines
            if i.startswith('virus') or not i.strip():
                continue
            
            line = i.strip().split('\t')
            contig = line[1]
            if contig not in host_dict:
                log_message(f"Error: Contig '{contig}' from virus info file not found in host FASTA file. Exiting.")
                sys.exit(1)
            host_seq = host_dict[contig].seq
            start = int(line[-2])
            end = int(line[-1])
            
            # Sanitize virus name for use in file names
            v_name=re.sub(r'[^a-zA-Z0-9_-]', '--', line[0])+'__'+contig+'__'+str(start)+'__'+str(end)
            v_seq=host_seq[start-1:end]
            original_info[v_name] = line
            f.write(f">{v_name}\n{v_seq}\n")

    qualified = {}
    if args.checkv_db:
        # Run CheckV to trim proviruses
        checkv_dir = os.path.join(preprocess_dir, "checkv")
        cmd = [
            "checkv", "end_to_end",
            original_virus_fna,
            checkv_dir,
            "-d", args.checkv_db,
            "-t", str(args.threads),
            "--quiet",
            "--remove_tmp"
        ]
        #log_message(f"Running CheckV command: {' '.join(cmd)}")
        try:
            # Using subprocess.run for better error handling
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            log_message(f"Error running CheckV: {e.stderr}")
            sys.exit(1)

        provirus_fna = os.path.join(checkv_dir, "proviruses.fna")
        virus_fna = os.path.join(checkv_dir, "viruses.fna")
        
        # Process viruses classified as complete by CheckV
        if os.path.exists(virus_fna):
            with open(virus_fna, 'r') as vi:
                for record in SeqIO.parse(vi, "fasta"):
                    name = record.id
                    contig = name.split('__')[-3]
                    host_seq = host_dict[contig].seq
                    host_len = len(host_seq)
                    start = int(name.split('__')[-2])
                    end = int(name.split('__')[-1])
                    
                    gap_l = start - 1
                    gap_r = host_len - end
                    
                    # Ensure there is enough flanking sequence for analysis
                    if gap_l > distance_cut and gap_r > distance_cut:
                        qualified[name] = [start, end, gap_l, gap_r, 'prophage_no_trim', host_len, contig]
                    else:
                        log_message(f'Warning: Skipping {name} due to insufficient flanking host sequence.')

        # Process sequences identified as proviruses (potentially trimmed) by CheckV
        if os.path.exists(provirus_fna):
            with open(provirus_fna, 'r') as pro:
                 for record in SeqIO.parse(pro, "fasta"):
                    header = record.description
                    name = record.id
                    contig = name.split('__')[-3]
                    host_seq = host_dict[contig].seq
                    host_len = len(host_seq)

                    # Extract coordinates from CheckV header (e.g., "1-12345/12345")
                    coords_match = re.search(r'(\d+)-(\d+)/(\d+)', header)
                    if not coords_match:
                        log_message(f"Warning: Could not parse coordinates from CheckV provirus header: {header}")
                        continue
                    
                    s, e, _ = map(int, coords_match.groups())
                    
                    original_start = int(name.split('__')[-2])
                    
                    # CheckV coordinates are 1-based relative to the input sequence
                    start = original_start + s - 1
                    end = original_start + e - 1
                    
                    gap_l = start - 1
                    gap_r = host_len - end
                    
                    # Ensure there is enough flanking sequence after trimming
                    if gap_l > distance_cut and gap_r > distance_cut:
                        qualified[name] = [start, end, gap_l, gap_r, 'phage_trim', host_len, contig]
                    else:
                        log_message(f'Warning: Skipping {name} due to insufficient flanking host sequence after trimming.')

    else:
        # If CheckV is not used, qualify all input viruses that meet the distance cut
        log_message('CheckV database not provided, skipping CheckV analysis.')
        with open(original_virus_fna, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                name = record.id
                contig = name.split('__')[-3]
                host_seq = host_dict[contig].seq
                host_len = len(host_seq)
                start = int(name.split('__')[-2])
                end = int(name.split('__')[-1])

                gap_l = start - 1
                gap_r = host_len - end

                if gap_l > distance_cut and gap_r > distance_cut:
                    qualified[name] = [start, end, gap_l, gap_r, 'prophage_no_trim', host_len, contig]
                else:
                    log_message(f'Warning: Skipping {name} due to insufficient flanking host sequence.')

    #log_message(f"Found {len(qualified)} qualified proviruses for border analysis.")
    return qualified, original_info


def extend_qualified_prophage_fna(host_dict, extend_inner=20000, extend_outer=5000):
    """
    For each qualified prophage, extract the flanking regions from the host genome.
    These extended regions are used to search for attachment sites (att).

    Args:
        extend_inner (int): The distance to extend inward from the prophage boundary.
        extend_outer (int): The distance to extend outward from the prophage boundary.

    Returns:
        tuple: A tuple containing:
            - dict: A dictionary with extended sequence information.
            - dict: The original input information.
    """
    qualified, original_info = preprocessing_virus_info(host_dict, distance_cut=DEFAULT_DISTANCE_CUT)
    if not qualified:
        log_message("No qualified proviruses found. Exiting.")
        sys.exit(0)

    extend_virus = {}
    qualified_provirus_path = os.path.join(working_path, QUALIFIED_PROVIRUS_FNA)
    with open(qualified_provirus_path, 'w') as v_fna:
        for name, data in qualified.items():
            contig = data[-1]
            start = data[0]
            end = data[1]
            host_len = data[-2]
            host_seq = host_dict[contig].seq
            
            # Write the qualified provirus sequence itself
            fna = host_seq[start-1:end]
            header = f">{name} {' '.join(map(str, data))}"
            v_fna.write(f"{header}\n{fna}\n")
        
            # Define the boundaries for the extended regions
            mid = (start + end) // 2

            # Right flank of the left boundary
            start_r = min(start + extend_inner, mid)
            # Left flank of the left boundary
            start_l = max(start - extend_outer, 1)

            # Left flank of the right boundary
            end_l = max(end - extend_inner, mid)
            # Right flank of the right boundary
            end_r = min(end + extend_outer, host_len)

            # Create FASTA entries for the left and right extended regions
            extend_name_l = f"{name}|-|{start_l}-{start_r}"
            extend_l_fna = str(host_seq[start_l-1:start_r])
            extend_name_r = f"{name}|-|{end_l}-{end_r}"
            extend_r_fna = str(host_seq[end_l-1:end_r])
    
            extend_virus[name.replace('|','__')] = [extend_name_l, extend_l_fna,
                                    extend_name_r, extend_r_fna,
                                    contig]
    
    #log_message(f"Extended flanking regions for {len(extend_virus)} proviruses.")
    return extend_virus, original_info

def run_makeblastdb(virus_name):
    """Runs makeblastdb for a given virus sequence."""
    db_in = os.path.join(working_path, "extended_fna", f"{virus_name}_r.fna")
    db_out = os.path.join(working_path, "extended_fna", f"{virus_name}_r_db")
    cmd = [
        "makeblastdb",
        "-in", db_in,
        "-dbtype", "nucl",
        "-out", db_out
    ]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        log_message(f"Error running makeblastdb for {virus_name}: {e}")
        sys.exit(1)

def run_blastn(virus_name, threads=args.threads):
    """Runs blastn to find repeats between flanking regions."""
    db_path = os.path.join(working_path, "extended_fna", f"{virus_name}_r_db")
    query_path = os.path.join(working_path, "extended_fna", f"{virus_name}_l.fna")
    out_path = os.path.join(working_path, "extended_fna", f"{virus_name}_l_vs_r.out")
    cmd = [
        "blastn",
        "-db", db_path,
        "-query", query_path,
        "-out", out_path,
        "-evalue", "100",
        "-word_size", "6",
        "-perc_identity", "80",
        "-dust", "no",
        "-strand", "plus",
        "-outfmt", "6",
        "-num_threads", str(threads)
    ]
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        log_message(f"Error running blastn for {virus_name}: {e}")
        sys.exit(1)


def base_by_base(att_l, att_r, fna1, fna2, boundary1=0, boundary2=0, min_len=12, max_len=30):
    """
    Performs a base-by-base comparison to find short, identical sequences
    between two larger sequences (fna1 and fna2), avoiding regions already
    identified by BLAST (att_l, att_r).

    Args:
        att_l (list): List of start/end coordinates of BLAST hits in fna1.
        att_r (list): List of start/end coordinates of BLAST hits in fna2.
        fna1 (str): The first sequence (left flanking region).
        fna2 (str): The second sequence (right flanking region).
        boundary1 (int): The genomic start coordinate of fna1.
        boundary2 (int): The genomic start coordinate of fna2.
        min_len (int): The minimum length of identical sequence to find.
        max_len (int): The maximum length of identical sequence to find.

    Returns:
        tuple: Two lists containing the coordinates of newly found identical sites.
    """
    att_l2 = [i for i in att_l]; att_r2 = [i for i in att_r]
    attl_new, attr_new = [], []
    seq_index_fna2 = {}

    # Create an index of all possible substrings in fna2 for quick lookup
    for i in range(max_len, min_len - 1, -1):
        for m in range(len(fna2) - i + 1):
            seq = fna2[m:m+i]
            if seq in seq_index_fna2:
                seq_index_fna2[seq].append(m)
            else:
                seq_index_fna2[seq] = [m]

    # Search for matching substrings in fna1
    for i in range(max_len, min_len - 1, -1):
        for j in range(len(fna1) - i + 1):
            seq1 = fna1[j:j+i]
            
            # Skip sequences containing 'N'
            if 'N' in seq1:
                continue

            if seq1 in seq_index_fna2:
                for m in seq_index_fna2[seq1]:
                    # Check if this new match overlaps with an existing BLAST hit
                    is_overlapping = any(
                        att_l2[n][0] <= (j + boundary1) and (j + i - 1 + boundary1) <= att_l2[n][1] and \
                        att_r2[n][0] <= (m + boundary2) and (m + i - 1 + boundary2) <= att_r2[n][1] 
                        for n in range(len(att_l2))
                    )

                    if not is_overlapping:
                        att_l2.append([j + boundary1, j + i + boundary1 - 1])
                        att_r2.append([m + boundary2, m + i + boundary2 - 1])
                        attl_new.append([j + 1, j + i])
                        attr_new.append([m + 1, m + i])

    return attl_new, attr_new

def process_blast_result(virus_name):
    """Parses blastn output to extract potential att sites."""
    att_l, att_r = [], []
    blast_out_path = os.path.join(working_path, 'extended_fna', f'{virus_name}_l_vs_r.out')
    
    if not os.path.exists(blast_out_path):
        log_message(f"Warning: BLAST output not found for {virus_name}. Skipping.")
        return att_l, att_r
        
    with open(blast_out_path, 'r') as blast:
        for i in blast:
            parts = i.strip().split('\t')
            # Ensure the line has enough columns and alignment length is > 11
            if len(parts) >= 10 and int(parts[3]) > 11:
                att_l.append([int(parts[6]), int(parts[7])])
                att_r.append([int(parts[8]), int(parts[9])])
    return att_l, att_r

def merge_attB(att_l, att_r, method, boundary_l, boundary_r, host_seq, virus_name, extend_len=300, file_name="merged_candidate_attB.fna"):
    """
    Takes a list of paired attachment sites, calculates their genomic coordinates,
    extracts the surrounding sequences, and writes them to a FASTA file.

    Args:
        att_l (list): List of coordinates on the left side.
        att_r (list): List of coordinates on the right side.
        method (str): The method used to find these sites (e.g., 'blastn').
        boundary_l (int): Genomic start of the left search region.
        boundary_r (int): Genomic start of the right search region.
        host_seq (str): The full host contig sequence.
        virus_name (str): The name of the virus being processed.
        extend_len (int): Length of flanking sequence to include around the att site.
        file_name (str): The output FASTA file name.

    Returns:
        tuple: Two lists containing the absolute genomic coordinates of the att sites.
    """
    attl_loc = []
    attr_loc = []
    for i in range(len(att_l)):
        # Convert relative coordinates to absolute genomic coordinates
        start1 = boundary_l + att_l[i][0] - 1
        end1 = boundary_l + att_l[i][1] - 1
        
        start2 = boundary_r + att_r[i][0] - 1
        end2 = boundary_r + att_r[i][1] - 1

        attl_loc.append([start1, end1])
        attr_loc.append([start2, end2])

        # Extract the sequences
        attl = host_seq[start1-1:end1]
        attr = host_seq[start2-1:end2]
        att = attl # The attachment site sequence

        # Extract the flanking regions for the attB candidate
        attB_r = host_seq[end2 : end2 + extend_len]
        if start1 > extend_len:
            attB_l = host_seq[start1 - extend_len - 1 : start1 - 1]
        else:
            attB_l = host_seq[:start1 - 1]

        attB = attB_l + att + attB_r
        
        # Create a detailed FASTA header
        position = f"{start1}-{end1}-{start2}-{end2}-{len(attB_l)}-{len(attB_r)}-{len(att)}"
        output_path = os.path.join(working_path, file_name)
        with open(output_path, "a") as f:
            f.write(f">{virus_name}_{position} {attl}-{attr} {method}\n{attB}\n")
    return attl_loc, attr_loc

def get_length(host_dict):
    """Reads the extended FASTA file and returns sequence lengths and anchor points."""
    a, c = {}, {}
    extended_fna_path = os.path.join(working_path, EXTENDED_5KB_FNA)
    
    if not os.path.exists(extended_fna_path):
        log_message(f"Error: Required file not found: {extended_fna_path}")
        sys.exit(1)

    # This can be simplified using SeqIO.parse
    with open(extended_fna_path, "r") as f:
        for i in f:
            i = i.strip()
            if i.startswith('>'):
                parts = i.split()
                key = parts[0][1:]
                host_info = parts[1]
                host = '_'.join(host_info.split('_')[:-1])
                anchor = int(host_info.split('_')[-1].split('-')[0])
                a[key] = ""
                c[key] = [host, anchor]
            else:
                a[key] += i
    
    b = {key: len(seq) for key, seq in a.items()}
    return b, c

def calculate_richness(length):
    """
    Calculates a 'richness' score for each position in the extended sequences.
    Richness is defined by the number of overlapping self-BLAST hits at that position,
    which indicates repetitive regions.

    Args:
        length (dict): A dictionary mapping sequence IDs to their lengths.

    Returns:
        dict: A dictionary mapping sequence IDs to a list of richness scores.
    """
    richness = {i: [0] * length[i] for i in length.keys()}
    hot_blast_dir = os.path.join(working_path, "hot_blast")

    # This function processes all BLAST outputs together.
    # It might be more robust to process them one by one if they are independent.
    blast_results = {}
    for filename in os.listdir(hot_blast_dir):
        if not filename.endswith("_mix_outfmt.txt"):
            continue
        
        with open(os.path.join(hot_blast_dir, filename)) as f:
            for line in f:
                j = line.split("\t")
                if len(j) < 10: continue # Skip malformed lines
                query_id, hit_id = j[0], j[1]
                coords = j[6:10]
                
                if query_id not in blast_results:
                    blast_results[query_id] = {}
                if hit_id not in blast_results[query_id]:
                    blast_results[query_id][hit_id] = []
                blast_results[query_id][hit_id].append(coords)

    # Process blast results to calculate richness
    for query_id, hits in blast_results.items():
        min_len = length.get(query_id, 0) / 3
        if min_len == 0: continue # Skip if length is zero

        for hit_alignments in hits.values():
            if len(hit_alignments) <= 1:
                continue
                
            intervals = []
            for coords in hit_alignments:
                try:
                    s_bp = min(int(coords[2]), int(coords[3]))
                    e_bp = max(int(coords[2]), int(coords[3]))
                    qs1 = min(int(coords[0]), int(coords[1]))
                    qe1 = max(int(coords[0]), int(coords[1]))
                    intervals.append((s_bp, e_bp, qs1, qe1))
                except (ValueError, IndexError):
                    continue # Skip malformed lines

            intervals.sort()
            
            # Check for overlaps between intervals
            for i, (s1, e1, qs1, qe1) in enumerate(intervals[:-1]):
                for s2, e2, qs2, qe2 in intervals[i+1:]:
                    if s2 > e1:
                        break # Optimization since list is sorted
                    
                    overlap = set(range(s1, e1 + 1)) & set(range(s2, e2 + 1))
                    # Check for significant overlap and distance between query starts
                    if len(overlap) > 4 and (abs(qs2 - qe1) > min_len or abs(qs1 - qe2) > min_len):
                        # Ensure indices are within bounds before incrementing richness
                        if qe1 <= len(richness[query_id]) and qe2 <= len(richness[query_id]):
                           for k in range(qs1 - 1, qe1): richness[query_id][k] += 1
                           for k in range(qs2 - 1, qe2): richness[query_id][k] += 1

    return richness

def change_point(length):
    """
    Identifies change points in the richness score, which correspond to the
    boundaries of highly repetitive regions (potential att sites).

    Args:
        length (dict): A dictionary mapping sequence IDs to their lengths.

    Returns:
        dict: A dictionary mapping sequence IDs to the coordinates of the hot region.
    """
    richness = calculate_richness(length)
    hot_region = {}
    
    # Use a sliding window to find the point of maximum change in richness
    window = 50
    for seq_id, rich_values in richness.items():
        seq_len = len(rich_values)
        if seq_len < 2 * (window + 50): # Ensure sequence is long enough
            continue
        mid_point = seq_len // 2
        
        # Helper to calculate window means
        def get_window_means(start, end):
            means = []
            if end - start < window: return means
            for i in range(start, end - window):
                window_mean = sum(rich_values[i:i+window]) / window
                means.append((i + window, window_mean))
            return means
            
        # --- Find change point in the left half ---
        left_means = get_window_means(0, mid_point - 50)
        left, max_left_rate = 0, 0
        for pos, mean1 in left_means:
            # Check if the next window is within bounds
            if pos + window > len(rich_values): continue
            mean2 = sum(rich_values[pos:pos+window]) / window
            rate = mean1 - mean2 # We are looking for a drop-off
            if rate > max_left_rate:
                max_left_rate = rate
                left = pos
                
        # --- Find change point in the right half ---
        right_means = get_window_means(mid_point, seq_len - 50)
        right, max_right_rate = 0, 0
        for pos, mean1 in right_means:
            # Check if the next window is within bounds
            if pos + window > len(rich_values): continue
            mean2 = sum(rich_values[pos:pos+window]) / window  
            rate = mean2 - mean1 # We are looking for a drop-off (from virus to host)
            if rate > max_right_rate:
                max_right_rate = rate
                right = pos
                
        if left and right:
            # Define the "hot region" around the change points
            hot_region[seq_id] = [left - 30, left + 5, right - 5, right + 30]
            
    return hot_region

def find_long_att_sites(host_dict, extend_virus, virus_anchor):
    """
    Finds long attachment sites (>=12 bp) by blasting left and right flanking
    regions against each other and using base-by-base comparison.
    """
    long_att_pairs = {}
    final_candidate_path = os.path.join(working_path, FINAL_CANDIDATE_ATTB_FNA)
    extended_fna_dir = os.path.join(working_path, "extended_fna")

    for i in extend_virus.keys():
        #log_message(f"  - Processing {i} for long att sites...")
        long_att_pairs[i] = [[], []]
        fna1 = extend_virus[i][1]
        fna2 = extend_virus[i][3]
        name1 = extend_virus[i][0]
        name2 = extend_virus[i][2]
        boundary_l = int(name1.split('|-|')[-1].split('-')[0])
        boundary_r = int(name2.split('|-|')[-1].split('-')[0])
        host = virus_anchor[i][0]
        host_seq = str(host_dict[host].seq)

        # BLAST left flank vs. right flank
        run_makeblastdb(i)
        run_blastn(i, threads=args.threads)

        # Process BLAST results
        att_l, att_r = process_blast_result(i)
        # Find additional short repeats not found by BLAST
        attl_new, attr_new = base_by_base(att_l, att_r, fna1, fna2, boundary1=1, boundary2=1, min_len=12, max_len=30)

        # Merge and store results
        attloc1, attloc2 = merge_attB(att_l, att_r, "blastn", boundary_l, boundary_r, host_seq, i, extend_len=100, file_name=FINAL_CANDIDATE_ATTB_FNA)
        attloc3, attloc4 = merge_attB(attl_new, attr_new, "base_comparing", boundary_l, boundary_r, host_seq, i, extend_len=100, file_name=FINAL_CANDIDATE_ATTB_FNA)

        long_att_pairs[i][0].extend(attloc1 + attloc3)
        long_att_pairs[i][1].extend(attloc2 + attloc4)

        # Clean up intermediate blast files for this virus
        for ext in ['_r_db.nsq', '_r_db.nin', '_r_db.nhr', '_l.fna', '_r.fna', '_l_vs_r.out']:
            try:
                os.remove(os.path.join(extended_fna_dir, f"{i}{ext}"))
            except OSError:
                pass
    return long_att_pairs

def find_short_att_sites_in_hot_regions(host_dict, hot_region, long_att_pairs, virus_anchor):
    """Finds short attachment sites (5-11 bp) within pre-defined hot regions."""
    #log_message(f"Found {len(hot_region)} potential hot regions for short attB search.")
    for i in hot_region.keys():
        virus_name = i
        p = hot_region[i]
        host, anchor = virus_anchor[i]

        boundary_l = p[0] + anchor + 1
        boundary_r = p[2] + anchor + 1
        host_seq = str(host_dict[host].seq)
        fna1 = host_seq[p[0] + anchor : p[1] + anchor]
        fna2 = host_seq[p[2] + anchor : p[3] + anchor]

        # Use the base-by-base method with a smaller length range
        attl, attr = base_by_base(long_att_pairs[i][0], long_att_pairs[i][1], fna1, fna2, boundary1=boundary_l, boundary2=boundary_r, min_len=5, max_len=30)

        attl2, attr2 = [], []
        for j in range(len(attl)):
            # Filter for short repeats
            if (attl[j][1] - attl[j][0] + 1) < 12:
                attl2.append(attl[j])
                attr2.append(attr[j])

        if attl2:
             merge_attB(attl2, attr2, "base_comparing_short", boundary_l, boundary_r, host_seq, virus_name, extend_len=100, file_name=FINAL_CANDIDATE_ATTB_FNA)


def run_att_search_pipeline(host_dict):
    """
    The main workflow function that coordinates the search for att sites.
    It identifies long repeats, finds hot regions, searches for short repeats
    within those regions, and then blasts all candidates against a database.

    Returns:
        dict: The original input information dictionary.
    """
    extend_virus, original_info = extend_qualified_prophage_fna(host_dict, extend_inner=20000, extend_outer=5000)
    
    # Exit if no qualified proviruses were found
    if not extend_virus:
        log_message("No qualified proviruses found to process. Exiting.")
        return None, None # Return None to signal main to exit

    extended_fna_dir = os.path.join(working_path, "extended_fna")
    os.makedirs(extended_fna_dir, exist_ok=True)
    
    # Write the extended sequences to a single file for the hot-spot finding BLAST
    extended_5kb_path = os.path.join(working_path, EXTENDED_5KB_FNA)
    with open(extended_5kb_path, "w") as f:
        for i in extend_virus.keys():
            start = int(extend_virus[i][0].split('|-|')[-1].split('-')[0])
            end = int(extend_virus[i][2].split('|-|')[-1].split('-')[1])
            fna = str(host_dict[extend_virus[i][4]].seq[start-1:end])
            name = str(extend_virus[i][4])+'_'+str(start)+'-'+str(end)
            f.write('>' + i + ' ' + name + '\n' + fna + '\n' )

            # Also write individual files for the left vs. right flank BLAST
            l_fna_path = os.path.join(extended_fna_dir, f'{i}_l.fna')
            r_fna_path = os.path.join(extended_fna_dir, f'{i}_r.fna')
            with open(l_fna_path, 'w') as final_l, open(r_fna_path, 'w') as final_r:
                name_l = i + '_l' + ' ' + extend_virus[i][0]
                name_r = i + '_r' + ' ' + extend_virus[i][2]
                final_l.write('>'+name_l +'\n' + extend_virus[i][1] + '\n')
                final_r.write('>'+name_r +'\n' + extend_virus[i][3] + '\n')
    
    # --- Step 2: Identify Candidate attB Sites ---
    log_message(f"[2/{TOTAL_STEPS}] Identifying candidate attB sites")
    length, virus_anchor = get_length(host_dict)
    
    final_candidate_path = os.path.join(working_path, FINAL_CANDIDATE_ATTB_FNA)
    if os.path.exists(final_candidate_path):
        os.remove(final_candidate_path) # Clean up previous results before appending

    # --- Part 1: Find long, identical repeats (>=12 bp) ---
    long_att_pairs = find_long_att_sites(host_dict, extend_virus, virus_anchor)

    # --- Part 2: Find "hot regions" by blasting extended regions against each other ---
    hot_blast_dir = os.path.join(working_path, "hot_blast")
    os.makedirs(hot_blast_dir, exist_ok=True)

    hot_blast_out = os.path.join(hot_blast_dir, "hot_mix_outfmt.txt")
    cmd_hot_blast = [
        "blastn", "-db", blast_db_path,
        "-query", extended_5kb_path,
        "-out", hot_blast_out,
        "-max_target_seqs", "1000000",
        "-num_threads", str(args.threads),
        "-evalue", "1e-10",
        "-outfmt", "6"
    ]
    if blast_db_path:
        #log_message("Running BLAST to find hot regions...")
        subprocess.run(cmd_hot_blast, check=True)
    else:
        log_message("No blastn database provided, skipping hot region detection.")
    
    # --- Part 3: Find short, identical repeats (5-11 bp) within hot regions ---
    hot_region = change_point(length)
    find_short_att_sites_in_hot_regions(host_dict, hot_region, long_att_pairs, virus_anchor)

    # --- Step 3: BLAST all candidates against the DB and Score ---
    attb_blast_dir = os.path.join(working_path, "attB_blast")
    os.makedirs(attb_blast_dir, exist_ok=True)
    
    log_message(f"[3/{TOTAL_STEPS}] Comparing candidate attB sites against database and scoring")
    
    attb_blast_out = os.path.join(attb_blast_dir, "attB_mix_outfmt.txt")
    if blast_db_path and os.path.exists(final_candidate_path) and os.path.getsize(final_candidate_path) > 0:
        cmd_attb_blast = [
            "blastn", "-db", blast_db_path,
            "-query", final_candidate_path,
            "-out", attb_blast_out,
            "-max_target_seqs", "100000",
            "-num_threads", str(args.threads),
            "-qcov_hsp_perc", "70",
            "-evalue", "1e-10",
            "-outfmt", "'6 qseqid sseqid pident length qstart qend sstart send evalue bitscore mismatch gaps frames qseq sseq'"
        ]
        # Note: The final part of outfmt is tricky with subprocess list args. Using shell=True for this one command.
        subprocess.run(" ".join(cmd_attb_blast), check=True, shell=True)
    else:
        log_message("Skipping final attB BLAST - no database provided or no candidates found.")

    return original_info, True


def sort_top_hit(a, b, c, top=1):
    """Sorts hits by score, bitscore, and details."""
    a, b, c = zip(*sorted(zip(a, b, c),reverse=True))
    a = a[:top]; b = b[:top]; c = c[:top]
    return a, b, c

def count_difference(str1, str2):
    """Calculates the number of differences between two strings."""
    return sum(char1 != char2 for char1, char2 in zip(str1, str2))


def attB_scoring(blast_file, ratio_cut=0, score_cut=DEFAULT_SCORE_CUTOFF, ratio2_cut=0.5):
    """
    Scores potential attB sites based on the final BLAST results against a database.
    The score is based on arm length, symmetry, and mismatches.

    Args:
        blast_file (str): Path to the BLAST output file.
        ratio_cut (float): Cutoff for the ratio of short arm to long arm length.
        score_cut (int): Cutoff for the final attB score.
        ratio2_cut (float): Cutoff for the ratio of mismatches between arms.

    Returns:
        dict: A dictionary of high-scoring attB hits.
    """
    if not os.path.exists(blast_file):
        log_message(f"Warning: BLAST file not found for scoring: {blast_file}")
        return {}
    attB_hit = {}
    with open(blast_file, "r") as blast:
        for j in blast:
            j = j.strip('\n'); i = j.split('\t')
            if len(i) < 15: continue # Malformed line
            
            qseq = i[-2]; sseq = i[-1]
            qs = int(i[4]); qe = int(i[5])
            
            # Extract info from the query ID
            header_parts = i[0].split('_')[-1].split('-')
            attlen = int(header_parts[-1])
            llen = int(header_parts[-3])
            
            larm = llen - qs + 1
            rarm = qe - attlen - llen
            total_mis = int(i[10]); total_gap = int(i[11])
            identity = float(i[2])
            
           
            if larm >= 20 and rarm >= 20:
                if total_gap == 0:
                    # Simple case: no gaps
                    llarm = larm; rrarm = rarm
                    short_arm = min(larm, rarm)
                    long_arm = max(larm, rarm)
                    ratio = short_arm / long_arm if long_arm > 0 else 0
                    score = short_arm * ratio - total_mis

                else:
                    # Complex case: recalculate arm lengths considering gaps
                    flag1=0; flag2=0
                    
                    for base in qseq:
                        flag1 += 1
                        if base != '-':
                            flag2 += 1
                            if flag2 == larm:
                                llarm = flag1
                            if flag2 == (larm + attlen):
                                break
                    rrarm = len(qseq) - flag1 

                    short_arm = min(llarm, rrarm)
                    long_arm = max(llarm, rrarm)
                    ratio = short_arm / long_arm if long_arm > 0 else 0
                    score = short_arm * ratio - total_mis - total_gap
                    
                l_dif = count_difference(qseq[:llarm], sseq[:llarm])
                r_dif = count_difference(qseq[-rrarm:], sseq[-rrarm:])
                att1 = qseq[llarm:len(qseq)-rrarm]
                att2 = sseq[llarm:len(qseq)-rrarm]
                
                core_gap = att1.count('-') + att2.count('-')
                core_dif = count_difference(att1, att2)
                
                # Ratio of mismatches in arms - promotes symmetric dissimilarity
                if (total_mis + total_gap) > 10 and max(l_dif, r_dif) != 0:
                    ratio2 = min(l_dif, r_dif) / max(l_dif, r_dif)
                else:
                    ratio2 = 1
                
                # Apply filters
                if ratio >= ratio_cut and score >= score_cut and ratio2 > ratio2_cut and core_gap < 3:
                    # Different criteria for long vs. short att sites
                    if attlen > 11:
                        # Store hit if it passes filters
                        if i[0] not in attB_hit:
                            attB_hit[i[0]] = [[], [], []] # details, bitscores, scores
                        attB_hit[i[0]][0].append(j + '\t' + str(ratio) + '\t' + str(l_dif)+'|'+str(r_dif)+ '\t'  + str(ratio2) + '\t' + str(score)+'\t'+str(core_dif)+'|'+att1+'|'+att2)
                        attB_hit[i[0]][1].append(int(float(i[9]))) # bitscore
                        attB_hit[i[0]][2].append(score)
                            
                    elif attlen <= 11:
                        # Stricter criteria for short att sites (must be near-perfect match)
                        if total_gap == 0 and identity >= 99:
                            if i[0] not in attB_hit:
                                attB_hit[i[0]] = [[], [], []]
                            attB_hit[i[0]][0].append(j + '\t' + str(ratio) + '\t' + str(l_dif)+'|'+str(r_dif)+ '\t'  + str(ratio2) + '\t' + str(score)+'\t'+str(core_dif)+'|'+att1+'|'+att2)
                            attB_hit[i[0]][1].append(int(float(i[9])))
                            attB_hit[i[0]][2].append(score)
                                
    # for i in attB_hit.keys():
    #     a = attB_hit[i][2]; b = attB_hit[i][1]; c = attB_hit[i][0];
    #     attB_hit[i] = [c, b, a]
    return attB_hit

def process_and_write_results(original_info):
    """
    Processes the final scored attB hits, selects the best one for each provirus,
    and writes the final output table.
    
    Args:
        original_info (dict): Dictionary mapping temporary names to original input data.
    """
    host_dict = SeqIO.to_dict(SeqIO.parse(host_fasta_path, "fasta"))
    att_prediction_path = os.path.join(working_path, ATT_PREDICTION_TSV)
    with open(att_prediction_path, "w") as f:
        header = [
            "original_name", "host_contig", "original_start", "original_end",
            "temp_name", "attL_start", "attL_end", "attR_start", "attR_end",
            "att_length", "attL_sequence", "attR_sequence", "attB_score"
        ]
        f.write("\t".join(header) + "\n")
        
        attb_blast_dir = os.path.join(working_path, "attB_blast")
        attB_all = defaultdict(lambda: ([], [], [])) # details, bitscores, scores
        
        # Aggregate all BLAST results from the scoring step
        if os.path.exists(attb_blast_dir):
            for item in os.listdir(attb_blast_dir):
                if item.endswith("_mix_outfmt.txt"):    
                    blast_file_path = os.path.join(attb_blast_dir, item)
                    attB_hit = attB_scoring(
                        blast_file_path, 
                        score_cut=args.score, 
                        ratio_cut=0, 
                        ratio2_cut=0
                    )
                    for j, (details, bitscores, scores) in attB_hit.items():
                        attB_all[j][0].extend(details)
                        attB_all[j][1].extend(bitscores)
                        attB_all[j][2].extend(scores)

        # --- Refined Selection Logic ---
        
        # 1. For each candidate att site, find its single best hit based on the three-level sort.
        best_hit_for_each_candidate = {}
        for candidate_id, (details_list, bitscores, scores) in attB_all.items():
            if not scores: continue
            
            # Create a list of tuples for robust sorting
            sortable_hits = []
            for i in range(len(scores)):
                details_str = details_list[i]
                # Extract att_length from the header of the details string
                header = details_str.split('\t')[0]
                att_len = int(header.split('-')[-1])
                sortable_hits.append( (scores[i], bitscores[i], att_len, details_str) )
            
            # Sort using the three-level criteria: score -> bitscore -> att_length
            sortable_hits.sort(key=lambda x: (x[0], x[1], x[2]), reverse=True)
            
            # Store the details of the best hit for this candidate
            best_hit_for_each_candidate[candidate_id] = sortable_hits[0][3] # The full details string

        # 2. For each original provirus, select the best candidate among its potential att sites.
        final_selection = {}
        for candidate_id, details_string in best_hit_for_each_candidate.items():
            provirus_name = '_'.join(candidate_id.split('_')[:-1])
            
            # Extract metrics from the winning hit's details string
            score = float(details_string.split('\t')[-2])
            bitscore = float(details_string.split('\t')[9])
            att_length = int(details_string.split('\t')[0].split('-')[-1])
            
            current_best = final_selection.get(provirus_name)
            if not current_best or \
               score > current_best['score'] or \
               (score == current_best['score'] and bitscore > current_best['bitscore']) or \
               (score == current_best['score'] and bitscore == current_best['bitscore'] and att_length > current_best['att_length']):
                
                final_selection[provirus_name] = {
                    'score': score,
                    'bitscore': bitscore,
                    'att_length': att_length,
                    'details': details_string,
                    'temp_name': candidate_id
                }

        # 3. Write the final, unique result for each provirus to the file.
        for name, selection_data in sorted(final_selection.items()):
            best_hit_details = selection_data['details']
            temp_name = selection_data['temp_name']
            
            att_info = best_hit_details.split('\t')
            att_header_parts = att_info[0].split('_')[-1].split('-')
            attl_r_s_e = att_header_parts[:4]
            att_length = att_header_parts[-1]
            score_val = att_info[-2]

            # Extract att sequences directly from host genome using coordinates
            contig = original_info[name][1]
            host_seq = str(host_dict[contig].seq)
            attl_start = int(attl_r_s_e[0])
            attl_end = int(attl_r_s_e[1])
            attr_start = int(attl_r_s_e[2])
            attr_end = int(attl_r_s_e[3])
            attl_seq = host_seq[attl_start - 1:attl_end]
            attr_seq = host_seq[attr_start - 1:attr_end]
            
            if name in original_info:
                output_line = '\t'.join(original_info[name]) + '\t' + temp_name + '\t' + \
                              '\t'.join(attl_r_s_e) + '\t' + att_length + '\t' + \
                              attl_seq + '\t' + attr_seq + '\t' + score_val + '\n'
                f.write(output_line)
            else:
                log_message(f"Warning: Could not find original info for temp_name '{name}'. Skipping from final report.")


def main():
    """Main function to run the ProBord pipeline."""
    # Setup initial logging
    global program_start_time
    program_start_time = time.time()
    try:
        log_message("Command: " + ' '.join(sys.argv))
        log_message(f"{PROG_NAME} started. Version: {__version__}")
    except Exception:
        pass

    #log_message("Loading host genome...")
    try:
        host_dict = SeqIO.to_dict(SeqIO.parse(host_fasta_path, "fasta"))
    except FileNotFoundError:
        log_message(f"Error: Host FASTA file not found at '{host_fasta_path}'")
        sys.exit(1)
    #log_message("Host genome loaded.")

    # This function is the main driver of the pipeline
    original_info, _ = run_att_search_pipeline(host_dict)
    
    # If the pipeline returned nothing, it means no qualified proviruses were found.
    if original_info is None:
        log_program_stats() # Log stats and exit
        return

    # Process the results and write the final output file
    process_and_write_results(original_info)

    # Clean up temporary files unless requested not to
    if not args.keep_temp:
        #log_message("Cleaning up temporary files...")
        shutil.rmtree(os.path.join(working_path, "extended_fna"), ignore_errors=True)
        shutil.rmtree(os.path.join(working_path, "preprocess"), ignore_errors=True)
        shutil.rmtree(os.path.join(working_path, "hot_blast"), ignore_errors=True)
        #shutil.rmtree(os.path.join(working_path, "attB_blast"), ignore_errors=True)
        
        for f_to_del in [QUALIFIED_PROVIRUS_FNA, EXTENDED_5KB_FNA, FINAL_CANDIDATE_ATTB_FNA]:
            try:
                os.remove(os.path.join(working_path, f_to_del))
            except OSError:
                pass
    else:
        log_message("--keep-temp specified: skipping cleanup.")
    
    # Log program stats
    log_program_stats()


if __name__ == "__main__":
    main()
