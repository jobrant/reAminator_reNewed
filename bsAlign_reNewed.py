#!/usr/bin/env python3

import os
import warnings
import sqlite3
from typing import List, Dict, Generator, Optional, Iterator, Tuple, Any
from dataclasses import dataclass
import multiprocessing
from itertools import islice
from collections import Counter
from subprocess import run, PIPE
import subprocess
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
import shutil
from dataclasses import dataclass
import tempfile
import time
import sys
import argparse
import logging
import cProfile
import pstats
from functools import partial
import re

import deaminate

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


# Check that necessary executables are loaded and in the path
def check_executables():
    required = ['makeblastdb', 'blastn', 'convert2blastmask']
    missing = [cmd for cmd in required if shutil.which(cmd) is None]
    if missing:
        logger.error(f"Required executables not found: {', '.join(missing)}. Ensure BLAST+ module is loaded.")
        raise EnvironmentError(f"Required executables not found: {', '.join(missing)}. "
                               f"Ensure BLAST+ module is loaded.")
    

# Function to swap file extensions
def change_extension(file_path: str, new_extension: str) -> str:
    """
    Change the file extension of the given file path to the new extension.
    
    Parameters:
        file_path (str): The original file path.
        new_extension (str): The new file extension (without a dot).
    
    Returns:
        str: The file path with the new extension.
    """
    base = os.path.splitext(file_path)[0]
    return f"{base}.{new_extension}"


# How to handle ambigous IUPAC DNA nucleotide codes
# this was never implemented, but I try to implement in future updates
iupac_codes = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC',
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG',
    'N': 'ACGT'
}

def get_iupac_intersection(code1, code2):
    bases1 = iupac_codes.get(code1, '')
    bases2 = iupac_codes.get(code2, '')
    intersection = set(bases1) & set(bases2)
    return ''.join(sorted(intersection)) if intersection else 'N'

ambig_iupac = {i + j: get_iupac_intersection(i, j)
               for i in iupac_codes for j in iupac_codes}


# checks for duplicated headers in the sequence reads
def validate_read_ids(records: str) -> bool:
    """
    Check for duplicated headers in the sequence reads.
    
    Parameters:
        records (str): Path to the FASTA file.
    
    Returns:
        bool: True if no duplicates found, False otherwise.
    """
    try:
        with open(records, 'r') as source:
            ids = [record.id for record in SeqIO.parse(source, 'fasta')]
    
        duplicates = [id for id, count in Counter(ids).items() if count > 1]
    
        if duplicates:
            logger.warning(f"Found {len(duplicates)} duplicate read IDs.")
            logger.warning(f"First few duplicates: {duplicates[:5]}")
            return False
        return True
    
    except IOError as e:
        logger.error(f"Error reading file {records}: {e}")
        logger.error(f"Does file exist? {os.path.exists(records)}")
        logger.error(f"File permissions: {oct(os.stat(records).st_mode)[-3:]}" if os.path.exists(records) else "N/A")
        return False

# **** the new rewrite of conversion as a C extension is not working well. Adding back in the functions that use the new C implimentation of fastools ****
# rewrote process_strands and deaminate functions as C extensions; imported at begining
# def process_strand(source, strand):
#     return deaminate.deaminate(source, strand)


# # Remove any non-ATCG characters
# def clean_sequence(seq):
#     return re.sub(r'[^ATCG]', '', seq.upper())


# def parallel_deaminate(source, strands='ab'):
#     with multiprocessing.Pool() as pool:
#         func = partial(process_strand, source)
#         results = pool.map(func, strands)
    
#     # Separate 'a' and 'b' strands
#     a_seqs = []
#     b_seqs = []
#     for strand, result in zip(strands, results):
#         for id, seq in result:
#             cleaned_seq = clean_sequence(seq)
#             if cleaned_seq:  # Only add if sequence is not empty after cleaning
#                 if strand == 'a':
#                     a_seqs.append((id.strip(), cleaned_seq))
#                 elif strand == 'b':
#                     b_seqs.append((id.strip(), cleaned_seq))
#             else:
#                 logger.warning(f"Empty sequence after cleaning for ID: {id}")
    
#     # Combine 'a' strands first, then 'b' strands
#     all_seqs = a_seqs + b_seqs
    
#     output_file = change_extension(source, ''.join(strands))
#     try:
#         with open(output_file, 'w') as handle:
#             for id, seq in all_seqs:
#                 try:
#                     record = SeqRecord(seq=Seq(seq), id=id, description='')
#                     SeqIO.write(record, handle, 'fasta')
#                 except Exception as e:
#                     logger.error(f"Error creating SeqRecord for ID {id}: {e}")
#     except IOError as e:
#         logger.error(f"Error writing to file {output_file}: {e}")
#         raise
    
#     return output_file


# Converts sequences using new C implementation of fastools convert
def process_strand(args: Tuple[str, str]) -> List[SeqIO.SeqRecord]:
    source, strand = args
    
    # Determine conversion type based on strand
    conv_type = "GA" if strand == 'b' else ""
    
    # Run faconvert and capture output
    try:
        result = run(['faconvert', source, conv_type], stdout=PIPE, text=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running faconvert: {e}")
        return
    
    seqs = []
    for seq in SeqIO.parse(StringIO(result.stdout), 'fasta'):
        if strand == 'b':
            seq.seq = seq.seq.reverse_complement()
        seq.id += f".{strand}"
        seq.description = ''
        seqs.append(seq)
    
    return seqs

def deaminate(source: str, strand: str = 'ab') -> None:
    """
    Deaminate sequences using faconvert.
    
    Parameters:
        source (str): Path to the input FASTA file.
        strand (str): Strand(s) to process ('a', 'b', or 'ab').
    """
    with multiprocessing.Pool() as pool:
        results = pool.map(process_strand, [(source, s) for s in strand])
    
    seqs = [seq for result in results for seq in result]
    
    output_file = change_extension(source, strand)
    try:
        with open(output_file, 'w') as handle:
            SeqIO.write(seqs, handle, 'fasta')
    except IOError as e:
        logger.error(f"Error writing to file {output_file}: {e}")
        raise


def run_blast(query: str, subject: str, xml_output: str, mask: bool = False, num_threads: int = 8) -> None:
    """
    Run BLAST alignment with optional masking.
    Parameters:
        query: Path to query sequence file
        subject: Path to subject sequence file
        xml_output: Path for output XML file
        mask: Whether to use masking
        num_threads: Number of threads to use for BLAST
    """
    check_executables()

    def run_command(command: List[str], shell: bool = False) -> None:
        try:
            logger.info(f"Running command: {' '.join(command)}")
            result = subprocess.run(command, check=True, shell=shell, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            logger.info(f"Command output: {result.stdout}")
            if result.stderr:
                logger.warning(f"Command stderr: {result.stderr}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error running command: {' '.join(command)}")
            logger.error(f"Stdout: {e.stdout}")
            logger.error(f"Stderr: {e.stderr}")
            raise

    # Check if files exist
    if not os.path.exists(query):
        raise FileNotFoundError(f"Query file not found: {query}")
    if not os.path.exists(subject):
        raise FileNotFoundError(f"Subject file not found: {subject}")

    # Prepare subject database
    db_files = [f"{subject}.{ext}" for ext in ['nhr', 'nin', 'nsq']]
    db_exists = all(os.path.exists(f) for f in db_files)

    if not db_exists or mask:
        logger.info("Creating BLAST database...")
        if mask:
            run_command([
                'convert2blastmask', '-in', subject, 
                '-out', change_extension(subject, 'asnb'),
                '-masking_algorithm', 'repeat',
                '-masking_options', 'repeatmasker,default',
                '-outfmt', 'maskinfo_asn1_bin', 
                '-parse_seqids'
            ])
            run_command([
                'makeblastdb', '-dbtype', 'nucl', 
                '-parse_seqids', '-in', subject,
                '-mask_data', change_extension(subject, 'asnb')
            ])
        else:
            run_command(['makeblastdb', '-dbtype', 'nucl', '-parse_seqids', '-in', subject])
    else:
        logger.info("BLAST database already exists, skipping creation.")

    # Prepare BLAST command
    blast_command = [
        'blastn', '-task', 'megablast',
        '-dust', '50 64 1',
        '-num_threads', str(num_threads),
        '-outfmt', '5',
        '-query', query,
        '-db', subject,
        '-out', xml_output
    ]

    if mask:
        blast_command.extend(['-db_soft_mask', '40'])

    # Run BLAST
    run_command(blast_command, shell=False)

    # Check if output was created
    if not os.path.exists(xml_output):
        raise FileNotFoundError(f"BLAST did not create output file: {xml_output}")

# Usage
# run_blast('query.fa', 'subject.fa', 'output.xml', mask=True)


# Parse function for blastn XML output
@dataclass
class BlastHit:
    subject_accession: str
    query_id: str
    expect_value: float
    subject_start: int
    subject_end: int
    query_start: int
    query_end: int
    subject_sequence: str
    query_sequence: str

def parse_blast_xml(file_path: str, e_value_threshold: float = 0.1) -> Iterator[BlastHit]:
    """
    Parse BLAST XML output file and yield high-quality hits.
    
    Parameters:
        file_path: Path to the BLAST XML file
        e_value_threshold: E-value threshold relative to the best hit (default: 0.1)
    
    Yield: 
        BlastHit objects for qualifying alignments
    """
    try:
        with open(file_path, 'r') as xml_file:
            for record in NCBIXML.parse(xml_file):
                if not record.alignments:
                    continue

                # Get all HSPs and their E-values
                all_hsps = [
                    (alignment.hsps[0].expect, i, alignment)
                    for i, alignment in enumerate(record.alignments)
                    if alignment.hsps
                ]

                if not all_hsps:
                    continue

                # Sort by E-value and get the best E-value
                all_hsps.sort(key=lambda x: x[0])
                best_e_value = all_hsps[0][0]

                for e_value, _, alignment in all_hsps:
                    if e_value > best_e_value * (1 / e_value_threshold):
                        break

                    hsp = alignment.hsps[0]
                    yield BlastHit(
                        subject_accession=alignment.accession,
                        query_id=record.query.split()[0],
                        expect_value=hsp.expect,
                        subject_start=hsp.sbjct_start,
                        subject_end=hsp.sbjct_end,
                        query_start=hsp.query_start,
                        query_end=hsp.query_end,
                        subject_sequence=hsp.sbjct,
                        query_sequence=hsp.query
                    )
    except IOError as e:
        logger.error(f"Error reading BLAST XML file {file_path}: {e}")
        raise


def reaminate(ref_seq: str, query_seq: str, conv_ref: str, conv_query: str) -> Seq:
    """
    Reconstruct methylation status of query sequence based on alignment with reference.

    Parameters:
        ref_seq: Original reference sequence
        query_seq: Original query sequence
        conv_ref: Converted reference sequence (C->T or G->A)
        conv_query: Converted query sequence (C->T or G->A)
    
    Returns:
        Seq: Reconstructed query sequence with methylation information
    """
    def insert_gaps(seq: List[str], conv: str) -> None:
        for i, char in enumerate(conv):
            if char == '-':
                seq.insert(i, '-')

    def check_surrounding(seq: List[str], index: int, char: str) -> bool:
        start = max(0, index - 2)
        end = min(len(seq), index + 3)
        return char in seq[start:end]

    # Convert sequences to lists for easier manipulation
    ref_list = list(ref_seq)
    query_list = list(query_seq)

    # Align sequences based on converted sequence alignment
    insert_gaps(ref_list, conv_ref)
    insert_gaps(query_list, conv_query)

    assert len(ref_list) == len(query_list) == len(conv_ref) == len(conv_query), "Sequences must have equal length after alignment"

    final = query_list.copy()

    for i in range(len(ref_list)):
        if ref_list[i] == 'C' and ref_list[i] != query_list[i]:
            if query_list[i] == '-':  # Deletion in query
                if check_surrounding(query_list, i, 'C'):
                    final[i] = 'C'  # Methylated
                elif check_surrounding(query_list, i, 'T'):
                    final[i] = 'T'  # Unmethylated
                else:
                    final[i] = 'N'  # Ambiguous
            elif query_list[i] == 'T':
                if check_surrounding(query_list, i, 'C') and check_surrounding(ref_list, i, '-'):
                    final[i] = 'C'  # Likely methylated insertion
                else:
                    final[i] = 't'  # Unmethylated

    # Remove insertions in reference
    final = [base for ref_base, base in zip(ref_list, final) if ref_base != '-']

    return Seq(''.join(final))

'''reconstructed_seq = reaminate(reference_sequence, query_sequence, converted_reference, converted_query)
print(reconstructed_seq)'''


def put_data(seqs: str, refs: str, dest: str, mask: bool = False, strand: str = 'ab', update: bool = False, batch_size: int = 5000):
    """
    Process sequence alignments and store results in an SQLite database.

    Parameters:
        seqs: Path to query sequences file
        refs: Path to reference sequences file
        dest: Path for output database
        mask: Whether to use masking in BLAST
        strand: Strand(s) to process ('a', 'b', or 'ab')
        update: Whether to update an existing database
    """
    temp_index = tempfile.mktemp(dir=os.environ['HOME'])

    # Ensure unique identifiers
    seqs_valid = validate_read_ids(seqs)
    refs_valid = validate_read_ids(refs)
    
    if not seqs_valid or not refs_valid:
        logger.error("Duplicate read IDs found. Please check your input files.")
        return

    # Define file paths
    db_path = change_extension(dest, 'db')
    xml_path = change_extension(seqs, 'xml')
    query_path = change_extension(seqs, strand)
    subject_path = change_extension(refs, 'ab')

    # Initial setup
    if not update and os.path.exists(db_path):
        os.remove(db_path)

    # Run BLAST if necessary
    if not os.path.exists(xml_path):
        logger.info('Converting bases...')
        start_time = time.time()
        deaminate(seqs, strand)
        deaminate(refs)
        
        logger.info(f'Conversion completed in {time.time() - start_time:.2f} seconds')

        logger.info('Aligning...')
        start_time = time.time()
        run_blast(query_path, subject_path, xml_path, mask)
        logger.info(f'Alignment completed in {time.time() - start_time:.2f} seconds')
    else:
        logger.info('XML file already exists. Delete for fresh alignment.')

    logger.info('Filing data...')
    
    
    # Create sequence index
    try:
        seqs_index = SeqIO.index_db(temp_index, [seqs, refs], 'fasta')
    except Exception as e:
        logger.error(f"Error creating sequence index: {e}")
        return

    # Set up SQLite database
    try:
        conn = sqlite3.connect(db_path)
        setup_database(conn)
    except sqlite3.Error as e:
        logger.error(f"Error setting up SQLite database: {e}")
        return

    # Process BLAST results
    try:
        records_to_insert = []
        records_to_update = []
        parse_start = time.time()
        for hit in parse_blast_xml(xml_path):
            start_time = time.time()
            query_parts = hit.query_id.split('.')
            subject_parts = hit.subject_accession.split('.')
            query_strand, subject_strand = query_parts[-1], subject_parts[-1]
            query_id, subject_id = '.'.join(query_parts[:-1]), '.'.join(subject_parts[:-1])

            if hit.subject_start > hit.subject_end:
                continue  # Discard wrong strand alignments

            try:
                query_seq_obj = seqs_index[query_id]
                subject_seq_obj = seqs_index[subject_id]
            except KeyError as e:
                logger.warning(f"KeyError: {e}. Skipping this alignment.")
                continue

            # Handle reverse complement
            if query_strand == 'b':
                query_seq_obj = query_seq_obj.reverse_complement()
            if subject_strand == 'b':
                subject_seq_obj = subject_seq_obj.reverse_complement()
            
             # Reconstruct methylation status
            try:
                # reaminate_start = time.time()
                reconstructed_seq = reaminate(
                    str(subject_seq_obj.seq[hit.subject_start - 1:hit.subject_end]).upper(),
                    str(query_seq_obj.seq[hit.query_start - 1:hit.query_end]).upper(),
                    hit.subject_sequence,
                    hit.query_sequence
                )
                # logger.info(f"reAminate time: {time.time() - reaminate_start:.2f} seconds")
            except AssertionError as e:
                logger.error(f"Assertion error in reaminate: {e}. Skipping this alignment.")
                continue

            # Pad sequence
            padded_seq = '-' * (hit.subject_start - 1) + str(reconstructed_seq) + '-' * len(subject_seq_obj.seq[hit.subject_end:])

            # Restore original orientation
            if subject_strand == 'b':
                padded_seq = str(Seq(padded_seq).reverse_complement())

            # Prepare record
            query_id, experiment = (query_id.split('|') + [''])[:2]
            record = (query_id, experiment, subject_id, hit.expect_value, subject_strand, padded_seq.upper())
                        
            if update:
                records_to_update.append(record)
            else:
                records_to_insert.append(record)

            # Bulk insert or update when batch size is reached
            if len(records_to_insert) >= batch_size or len(records_to_update) >= batch_size:
               logger.info(f"Time spent processing before insertion: {time.time() - start_time:.2f} seconds")
               if update:
                   insert_start = time.time()
                   bulk_update_records(conn, records_to_update, batch_size=batch_size)
                   logger.info(f"Actual insertion time: {time.time() - insert_start:.2f} seconds")
                   records_to_update = []
               else:
                   insert_start = time.time()
                   bulk_insert_sequences(conn, records_to_insert, batch_size=batch_size)
                   records_to_insert = []
                   logger.info(f"Total time for this batch: {time.time() - insert_start:.2f} seconds")
                   logger.info(f"Parse blast xml time: {time.time() - parse_start:.2f} seconds")

        # Insert or update remaining records
        if update:
            bulk_update_records(conn, records_to_update, batch_size=batch_size)
        else:
            bulk_insert_sequences(conn, records_to_insert, batch_size=batch_size)

    except Exception as e:
        logger.error(f"Error processing BLAST results: {e}")
    finally:
        # Add references to database
        add_references_to_db(conn, seqs_index)

        logger.info(f'Data filing completed in {time.time() - parse_start:.2f} seconds')
        print_alignment_stats(conn)

        conn.commit()
        conn.close()
        os.remove(temp_index)

def setup_database(conn: sqlite3.Connection):
    """Set up the SQLite database schema"""
    try:
        cursor = conn.cursor()
        cursor.execute('''CREATE TABLE IF NOT EXISTS records 
                        (id INTEGER PRIMARY KEY NOT NULL, 
                        read TEXT, expt TEXT, locus TEXT, 
                        expect REAL, strand TEXT, sequence TEXT)''')
        cursor.execute('''CREATE TABLE IF NOT EXISTS loci 
                        (id INTEGER PRIMARY KEY NOT NULL, locus TEXT, sequence TEXT)''')
        cursor.execute('CREATE INDEX IF NOT EXISTS loc_idx ON loci (locus)')
        cursor.execute('CREATE INDEX IF NOT EXISTS els_idx ON records (expt, locus, strand)')
        cursor.execute('CREATE INDEX IF NOT EXISTS rls_idx ON records (read, locus, strand)')
        conn.commit()
    except sqlite3.Error as e:
        logger.error(f"Error setting up SQLite database schema: {e}")
        conn.rollback()
        raise


def bulk_insert_sequences(conn: sqlite3.Connection, records: List[Tuple], batch_size: int = 1000):
    """Insert multiple sequence records into the database."""
    try:
        cursor = conn.cursor()
        start_time = time.time()
        total_records = len(records)
        
        cursor.execute('BEGIN TRANSACTION')
        for i in range(0, total_records, batch_size):
            batch = records[i:i+batch_size]
            cursor.executemany('''
                INSERT INTO records (read, expt, locus, expect, strand, sequence)
                VALUES (?, ?, ?, ?, ?, ?)
            ''', batch)
        cursor.execute('COMMIT')
        
        end_time = time.time()
        duration = end_time - start_time
        logger.info(f"Successfully inserted {total_records} records in {duration:.2f} seconds")
        logger.info(f"Average insertion rate: {total_records / duration:.2f} records/second")
        logger.info(f"Batch size used: {batch_size}")
    except sqlite3.Error as e:
        logger.error(f"Error inserting records: {e}")
        conn.rollback()
        raise

def add_references_to_db(conn: sqlite3.Connection, seqs_index: Dict[str, SeqIO.SeqRecord]):
    """Add new reference sequences to the database."""
    try:
        cursor = conn.cursor()
        
        # Get existing loci
        cursor.execute('SELECT locus FROM loci')
        existing_loci = set(row[0] for row in cursor.fetchall())
        
        # Get new loci
        cursor.execute('SELECT DISTINCT locus FROM records')
        all_loci = set(row[0] for row in cursor.fetchall())
        new_loci = all_loci - existing_loci
        
        logger.info(f"Adding {len(new_loci)} new reference sequences to the database")
        
        new_references = []
        for locus in new_loci:
            if locus not in seqs_index:
                logger.warning(f"Locus {locus} not found in seqs_index, skipping")
                continue
            
            seq_record = seqs_index[locus]
            sequence_str = str(seq_record.seq)
            new_references.append((seq_record.id, sequence_str))
        
        cursor.executemany("INSERT INTO loci (locus, sequence) VALUES (?, ?)", new_references)
        conn.commit()
        logger.info("Successfully added reference sequences to the database")
        
    except sqlite3.Error as e:
        logger.error(f"Error adding references to database: {e}")
        conn.rollback()
        raise
    except KeyError as e:
        logger.error(f"KeyError when accessing seqs_index: {e}")
        conn.rollback()
        raise
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        conn.rollback()
        raise


# rewrote for adding records in bulk instead of one at a time. This should save a huge amount of time!
def bulk_update_records(conn: sqlite3.Connection, records: List[Tuple], batch_size: int = 1000):
    """Update multiple existing records in the database."""
    try:
        cursor = conn.cursor()
        start_time = time.time()
        total_records = len(records)
        
        cursor.execute('BEGIN TRANSACTION')
        for i in range(0, total_records, batch_size):
            batch = records[i:i+batch_size]
            cursor.executemany('''
                UPDATE records
                SET expt = ?, locus = ?, expect = ?, strand = ?, sequence = ?
                WHERE read = ?
            ''', batch)
        cursor.execute('COMMIT')
        
        end_time = time.time()
        duration = end_time - start_time
        logger.info(f"Successfully updated {total_records} records in {duration:.2f} seconds")
        logger.info(f"Average update rate: {total_records / duration:.2f} records/second")
    except sqlite3.Error as e:
        logger.error(f"Error updating records: {e}")
        conn.rollback()
        raise


def print_alignment_stats(conn: sqlite3.Connection):
    """Print alignment statistics."""
    cursor = conn.cursor()
    try:
        num_reads = len(list(cursor.execute('SELECT read FROM records GROUP BY read')))
        num_loci = len(list(cursor.execute('SELECT locus FROM records GROUP BY locus')))
        logger.info(f'Aligned {num_reads} sequences to {num_loci} loci')
    except sqlite3.Error as e:
        logger.error(f"Error retrieving alignment statistics: {e}")


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Align and analyze bisulfite-converted DNA sequences.')
    parser.add_argument('seqs', help='FASTA-format file of read sequences')
    parser.add_argument('refs', help='FASTA-format file of reference sequences')
    parser.add_argument('-Save', help='output file name')
    parser.add_argument('-mask', action='store_true', help='enable lowercase masking in references')
    parser.add_argument('-pair', type=str, help='2nd FASTA-format file for paired-end sequencing')
    parser.add_argument('-update', action='store_true', help='update existing database')
    parser.add_argument('-1', '--a1_or_b1', action='append_const', const='a',
                        dest='strand', help='convert BS-read G to A, i.e. seq. primer=a1 or b1')
    parser.add_argument('-2', '--a2_or_b2', action='append_const', const='b',
                        dest='strand', help='convert BS-read C to T, i.e. seq. primer=a2 or b2')
    parser.add_argument('--profile', action='store_true', help='Enabling Profiling')
    parser.add_argument('--batch-size', type=int, default=5000, help='Batch size for database operations (default: 5000)')
    return parser.parse_args()


def determine_strand(args: argparse.Namespace) -> str:
    """Determine the strand(s) to process based on command-line arguments."""
    if not args.strand:
        return 'ab'
    return ''.join(sorted(args.strand)).lower()


def generate_output_filename(args: argparse.Namespace) -> str:
    """Generate an output filename if not provided."""
    if args.Save:
        return args.Save
    
    seqs_base = os.path.splitext(os.path.basename(args.seqs))[0][:10]
    refs_base = os.path.splitext(os.path.basename(args.refs))[0][:10]
    return os.path.join(os.path.dirname(args.seqs), f"{seqs_base}.{refs_base}")


def print_arguments(args: Dict[str, Any]) -> None:
    """Print the parsed command-line arguments."""
    for arg, val in sorted(args.items()):
        logger.info(f"{arg}: {val}")
    logger.info("")


def main() -> None:
    """Main function to run the script."""
    args = parse_arguments()
    args.strand = determine_strand(args)
    args.Save = generate_output_filename(args)

    print_arguments(vars(args))

    # Info for testing to make sure things are working correctly
    # logger.info(f"Current working directory: {os.getcwd()}")
    # logger.info(f"Query file: {os.path.abspath(args.seqs)}")
    # logger.info(f"Reference file: {os.path.abspath(args.refs)}")
    
    try:
        logger.info(f"Current working directory: {os.getcwd()}")
        logger.info(f"Seqs file absolute path: {os.path.abspath(args.seqs)}")
        logger.info(f"Refs file absolute path: {os.path.abspath(args.refs)}")
        logger.info(f"Seqs file exists: {os.path.exists(args.seqs)}")
        logger.info(f"Refs file exists: {os.path.exists(args.refs)}")

        if args.pair:
            logger.info(f"Pair Seqs file absolute path: {os.path.abspath(args.pair)}")
            logger.info(f"Pair Seqs file exists: {os.path.exists(args.pair)}")
    except Exception as e:
        logger.error(f"Error checking file existence: {e}")

    try:
        start_time = time.time()
        put_data(args.seqs, args.refs, args.Save, args.mask, args.strand, args.update, args.batch_size)
        logger.info(f"R1 processing finished in {time.time() - start_time:.2f} seconds")
        if args.pair:
            pair_time = time.time()
            paired_strand = {'a': 'b', 'b': 'a', 'ab': 'ab'}[args.strand]
            put_data(args.pair, args.refs, args.Save, args.mask, paired_strand, True, args.batch_size)
            logger.info(f"R2 processing finished in {time.time() - pair_time:.2f} seconds")
        logger.info(f"Total time {time.time() - start_time:.2f} seconds")
    except Exception as e:
        logger.error(f"An error occurred: {e}")
        sys.exit(1)
    
    logger.info(f"bsAlign_reNewed.py successfully finished")

if __name__ == '__main__':
    # cProfile.run('main()', 'output.prof')
    main()
    # Print sorted stats
    # p = pstats.Stats('ouput.prof')
    # p.sort_stats('cumulative').print_stats(30)
