#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
Supermatrix_Pipeline-NM123.py
--------------------------------------------------
Finalized workflow (Version FINAL 2025-11 + CD-HIT Target PreFilter):

Overall Workflow:
    Step0: Protein FASTA files under `target_dir` → CD-HIT (Reduce redundancy) (c=0.98) 

    Step1: Ref-FASTA → QC → CD-HIT → Working_Dir/ref.fasta
           (Automated filtering + clustering to derive representative reference sequences)

    Step2: All sequences in CDHIT_Target/ vs ref.fasta
           → AutoBlastPExtractor → FastaMerged/{gene}.fasta
           (Genome-wide BLASTp screening and merging per gene family)

    Step3: FastaMerged/{gene}.fasta + Pfam_hmmserach-config
           → Pfam_hmmserach-config_Results/{gene}.fasta
           (Domain-level HMM validation and refinement)

    Step4: Pfam_hmmserach-config_Results/{gene}.fasta
           → MAFFT (auto) + TrimAl (auto) + FastTree + phylopypruner
           → MAFFT (localpair) + TrimAl (gt 0.2)
           → concatenated supermatrix (need iq-tree LG+C60+F+G) + partition + statistics
           (Two-stage alignment, trimming, phylogenetic pruning, concatenation)

All logs are centrally managed through PipelineLogger.
"""

# ==============================================================
#                     User Config
# ==============================================================

taget_dir = "Demo_data"
ref_seq_dir = "Ref-datasets/NM123/Ref-FASTA"
ref_hmm_dir = "Ref-datasets/NM123/Ref-HMM"
supermatrix_outdir = "Demo_data-NM123"
prefix_concatenation = "Demo_data-NM123"

BLAST_THREADS_PER_TASK = 3
BLAST_MAX_PARALLEL = 20

# CD-HIT target (c=0.98) output directory name
CDHIT_TARGET_DIRNAME = "Working_Dir/CDHIT_Target"
# ==============================================================
#                     Imports
# ==============================================================

import os
import sys
import time
import csv
import shutil
import subprocess
import statistics
from glob import glob
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
from Bio import SeqIO


# ==============================================================
#                   Pipeline Logger (Unified Logging)
# ==============================================================

class PipelineLogger:
    """
    Unified pipeline logger providing:
        - Pipeline-level timing
        - Formatted time conversion
        - Standardized progress display
        - Automatic writing to persistent log file

    Every pipeline step calls this logger to ensure consistent reporting.
    """

    def __init__(self, logfile_path):
        self.logfile_path = logfile_path
        self.pipeline_start = time.time()

        log_dir = os.path.dirname(os.path.abspath(self.logfile_path))
        if log_dir and not os.path.exists(log_dir):
            os.makedirs(log_dir, exist_ok=True)

    @staticmethod
    def fmt_time(sec):
        """Return human-readable H:M:S time."""
        sec = int(sec)
        h = sec // 3600
        m = (sec % 3600) // 60
        s = sec % 60
        return f"{h:02d}:{m:02d}:{s:02d}"

    def write(self, msg):
        """Write a message to both stdout and log file."""
        line = str(msg)
        print(line)
        with open(self.logfile_path, "a", encoding="utf-8") as f:
            f.write(line + "\n")

    def progress(self, script_name, step_label, current, total, step_start_time):
        """
        Display and record progress of a running step.

        Parameters:
            script_name: Name of subscript/module (e.g. Script1)
            step_label:  Description of current step
            current:     Completed tasks
            total:       Total tasks
            step_start_time: Timestamp when the step started
        """
        if total <= 0:
            total = 1

        step_elapsed = time.time() - step_start_time
        total_elapsed = time.time() - self.pipeline_start
        pct = current / total * 100.0

        line = (
            f"[{script_name}] [{step_label}]  "
            f"{current}/{total} ({pct:5.1f}%) | "
            f"Step elapsed {self.fmt_time(step_elapsed)} | "
            f"Total elapsed {self.fmt_time(total_elapsed)}"
        )
        self.write(line)
# ==============================================================
# Script1: AutoFilterCDHIT — Ref-FASTA → QC → CD-HIT → ref.fasta
# ==============================================================

class AutoFilterCDHIT:
    """
    Script1 performs reference FASTA preprocessing, including:

    1. QC on each gene FASTA:
        - Remove sequences containing ambiguous amino acid "X"
        - Remove sequences outside dynamic length thresholds
          (mean ± SD × multiplier)
        - Remove gaps

    2. CD-HIT clustering for each gene family to obtain
       representative (non-redundant) sequences.

    3. Merge all representative sequences into:
           Working_Dir/ref.fasta

    This file becomes the BLAST reference database in Step2.
    """

    def __init__(self,
                 working_dir,
                 threads=12,
                 identity=0.7,
                 wordsize=5,
                 sd_multiplier=2,
                 output_merged="ref.fasta",
                 logger=None):

        self.working_dir = os.path.abspath(working_dir)
        self.THREADS = threads
        self.CDHIT_IDENTITY = identity
        self.CDHIT_WORDSIZE = wordsize
        self.SD_MULTIPLIER = sd_multiplier
        self.OUTPUT_MERGED = output_merged
        self.logger = logger

        # Supported reference file extensions
        self.support_ext = ("*.fasta", "*.fa", "*.faa", "*.fna", "*.pep")


    # --------------------------------------------------------------
    # QC routine — remove X, strip gaps, dynamic length filtering
    # --------------------------------------------------------------

    def filter_fasta(self, fasta_file):
        """
        Perform quality control on a reference FASTA.

        Steps:
        - Remove sequences containing 'X'
        - Remove sequences with extreme lengths based on mean ± SD×multiplier
        - Remove alignment gaps (-)
        """
        records = list(SeqIO.parse(fasta_file, "fasta"))
        if not records:
            return None

        # Remove sequences containing X
        noX = [r for r in records if "X" not in str(r.seq).upper()]
        if not noX:
            return None

        # Length-based filtering
        lens = [len(str(r.seq).replace("-", "")) for r in noX]
        if len(lens) < 2:
            return None

        mean_len = statistics.mean(lens)
        sd_len = statistics.stdev(lens)
        low = mean_len - self.SD_MULTIPLIER * sd_len
        high = mean_len + self.SD_MULTIPLIER * sd_len

        filtered = []
        for r in noX:
            seq_ng = str(r.seq).replace("-", "")
            if low <= len(seq_ng) <= high:
                r.seq = type(r.seq)(seq_ng)
                filtered.append(r)

        if not filtered:
            return None

        out_file = os.path.splitext(fasta_file)[0] + "-filtered.fasta"
        SeqIO.write(filtered, out_file, "fasta")
        return out_file


    # --------------------------------------------------------------
    # CD-HIT wrapper
    # --------------------------------------------------------------

    def run_cdhit(self, fasta_file):
        """
        Run CD-HIT on a filtered FASTA file to extract representative sequences.
        """
        gene = os.path.splitext(os.path.basename(fasta_file))[0].replace("-filtered", "")
        out_rep = f"{gene}-filtered_representatives.fasta"

        cmd = [
            "cd-hit",
            "-T", str(self.THREADS),
            "-c", str(self.CDHIT_IDENTITY),
            "-n", str(self.CDHIT_WORDSIZE),
            "-i", fasta_file,
            "-o", out_rep
        ]

        subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True
        )

        return out_rep


    # --------------------------------------------------------------
    # Merge representative sequences into ref.fasta
    # --------------------------------------------------------------

    def merge_representatives(self, rep_files):
        """
        Combine all CD-HIT representative sequences from individual genes
        into a single reference database:
            Working_Dir/ref.fasta
        """
        out_path = os.path.join(self.working_dir, self.OUTPUT_MERGED)
        total = 0

        with open(out_path, "w") as out:
            for f in rep_files:
                gene = os.path.basename(f).split("-filtered_representatives.fasta")[0]
                for rec in SeqIO.parse(f, "fasta"):
                    rec.id = f"{gene}@{rec.id}"  # tag each sequence with gene name
                    rec.description = ""
                    SeqIO.write(rec, out, "fasta")
                    total += 1

        return out_path, total


    # --------------------------------------------------------------
    # run() — main entry for Script1
    # --------------------------------------------------------------

    def run(self):
        """
        Main routine:
        - Scan reference directory
        - QC each FASTA
        - Run CD-HIT
        - Merge representatives into ref.fasta
        """
        fasta_files = []
        for ext in self.support_ext:
            fasta_files.extend(glob(os.path.join(self.working_dir, ext)))

        fasta_files = sorted(fasta_files)
        total = len(fasta_files)

        if total == 0:
            if self.logger:
                self.logger.write(
                    "[Script1] Warning: No reference sequences found "
                    "(*.fasta/fa/fna/faa/pep)"
                )
            return

        if self.logger:
            self.logger.write(f"[Script1] Number of reference gene families: {total}")

        rep_files = []
        temp_files = []
        step_start = time.time()

        for idx, fa in enumerate(fasta_files, start=1):
            # QC
            filtered = self.filter_fasta(fa)
            if filtered:
                temp_files.append(filtered)

                # CD-HIT
                rep = self.run_cdhit(filtered)
                rep_files.append(rep)
                temp_files.append(rep)

                clstr = rep + ".clstr"
                if os.path.exists(clstr):
                    temp_files.append(clstr)

            # Progress display
            if self.logger:
                self.logger.progress(
                    "Script1",
                    "QC + CD-HIT",
                    idx,
                    total,
                    step_start
                )

        # Merge representative sequences
        if rep_files:
            out_path, total_rep = self.merge_representatives(rep_files)
            if self.logger:
                self.logger.write(
                    f"[Script1] Generated ref.fasta: {out_path} "
                    f"(Total representative sequences: {total_rep})"
                )
        else:
            if self.logger:
                self.logger.write(
                    "[Script1] Warning: No representative sequences obtained; ref.fasta not generated."
                )

        # Cleanup temporary files
        for f in temp_files:
            try:
                os.remove(f)
            except Exception:
                pass
# ==============================================================
# Script2: AutoBlastPExtractor — BLASTp vs ref.fasta
# ==============================================================

class AutoBlastPExtractor:
    """
    Step2:
    Perform BLASTp queries of all CD-HIT–filtered genome protein FASTA files
    against the representative reference database (ref.fasta) generated in Step1.

    Features:
        - Automatic BLAST database creation
        - Multi-threaded BLAST execution across genomes
        - Hit filtering by identity and coverage thresholds
        - Best-hit extraction
        - Gene-level grouping (multi-reference support)
        - Exporting all accepted hits into per-gene FASTA files
        - Merging of per-genome FASTA into gene-level FASTA (FastaMerged/)
        - Construction of Genome × Gene presence/absence and count matrix

    Output directory structure (inside results_dir):
        Fasta/
        FastaMerged/
        Csv/
        Csv-threshold/
        Csv-threshold-best/
        Csv-threshold-best-clean/
        BLASTp_Statistics.csv
        BLASTp_Statistics_Normalized.csv
    """

    def __init__(self,
                 ref_fasta,
                 ref_db_prefix,
                 query_dir,
                 evalue=1e-5,
                 identity=30,
                 coverage=70,
                 threads_per_blast=2,
                 parallel_jobs=8,
                 allow_multi_ref=True,
                 add_genome_prefix=True,
                 force_rebuild_db=False,
                 results_dir=None,
                 logger=None):

        # Absolute paths
        self.REF_FASTA = os.path.abspath(ref_fasta)
        self.REF_DB = os.path.abspath(ref_db_prefix)
        self.QUERY_DIR = os.path.abspath(query_dir)

        # Threshold parameters
        self.E_VALUE = evalue
        self.IDENTITY_THRESHOLD = identity
        self.COVERAGE_THRESHOLD = coverage
        self.THREADS_PER_BLAST = threads_per_blast

        # BLAST outfmt
        self.OUTFMT = (
            "10 qseqid qlen sseqid slen pident qcovs "
            "qstart qend sstart send evalue bitscore"
        )

        # Parallel BLAST settings
        self.PARALLEL_JOBS = parallel_jobs
        self.ALLOW_MULTI_REF = allow_multi_ref
        self.ADD_GENOMEPREFIX = add_genome_prefix
        self.FORCE_REBUILD_DB = force_rebuild_db

        self.logger = logger

        # Output directory settings
        if results_dir is None:
            self.RESULTS_DIR = os.path.join(
                os.path.dirname(self.REF_FASTA),
                "AutoBlastPExtractor_Results"
            )
        else:
            self.RESULTS_DIR = os.path.abspath(results_dir)

        self.SUBDIRS = {
            "fasta": os.path.join(self.RESULTS_DIR, "Fasta"),
            "merged": os.path.join(self.RESULTS_DIR, "FastaMerged"),
            "csv": os.path.join(self.RESULTS_DIR, "Csv"),
            "thr": os.path.join(self.RESULTS_DIR, "Csv-threshold"),
            "best": os.path.join(self.RESULTS_DIR, "Csv-threshold-best"),
            "clean": os.path.join(self.RESULTS_DIR, "Csv-threshold-best-clean"),
        }

        for p in self.SUBDIRS.values():
            os.makedirs(p, exist_ok=True)
        os.makedirs(self.RESULTS_DIR, exist_ok=True)

        # Track completed genomes
        self.DONE_FILE = os.path.join(self.RESULTS_DIR, "_done.txt")
        self.done_set = set()
        if os.path.exists(self.DONE_FILE):
            with open(self.DONE_FILE) as f:
                self.done_set = {line.strip() for line in f if line.strip()}


    # --------------------------------------------------------------
    # makeblastdb
    # --------------------------------------------------------------

    def make_blastdb_if_needed(self):
        """
        Create BLAST database from ref.fasta when:
        - Database does not exist, or
        - force_rebuild_db == True
        """
        need = [self.REF_DB + ext for ext in [".pin", ".psq", ".phr"]]
        if all(os.path.exists(f) for f in need) and not self.FORCE_REBUILD_DB:
            return

        cmd = [
            "makeblastdb",
            "-in", self.REF_FASTA,
            "-dbtype", "prot",
            "-out", self.REF_DB
        ]
        subprocess.run(cmd,
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL,
                       check=True)


    # --------------------------------------------------------------
    # sseqid cleanup helper
    # --------------------------------------------------------------

    @staticmethod
    def clean_sseqid(s):
        """
        Clean subject sequence IDs returned by BLASTp.

        Typically removes extra fields like:
            ref|WP_12345.1|
        → WP_12345.1
        """
        if "|" in s:
            return s.split("|")[1]
        return s.split()[0]

    @staticmethod
    def load_ref_ids(path):
        """
        Load all IDs from reference FASTA.
        """
        ids = []
        for rec in SeqIO.parse(path, "fasta"):
            ids.append(rec.id.split()[0])
        return ids, set(ids)


    # --------------------------------------------------------------
    # Export accepted hit sequences to fasta/
    # --------------------------------------------------------------

    def extract_hits_to_fasta(self, genome, faa, df_clean, output_dir, add_prefix=False):
        """
        Export accepted BLASTp hits from a genome into gene-wise FASTA files.

        If add_prefix=True:
            sequence IDs become:
                genome@sequence_id
        """
        seqs = SeqIO.to_dict(SeqIO.parse(faa, "fasta"))
        group_col = "gene" if (self.ALLOW_MULTI_REF and "gene" in df_clean.columns) else "sseqid"

        for gene, group in df_clean.groupby(group_col):
            out_fa = os.path.join(output_dir, f"{gene}_x_{genome}_vsref.fasta")

            with open(out_fa, "a") as fw:
                for qid in group["qseqid"]:
                    if qid in seqs:
                        rec = seqs[qid]
                        if add_prefix:
                            rec.id = f"{genome}@{rec.id}"
                            rec.description = ""
                        SeqIO.write(rec, fw, "fasta")


    # --------------------------------------------------------------
    # Process one genome FASTA through BLASTp
    # --------------------------------------------------------------

    def _process_one_genome(self, faa_file):
        genome = os.path.splitext(os.path.basename(faa_file))[0]

        if genome in self.done_set:
            return genome

        out_csv = os.path.join(self.SUBDIRS["csv"], f"{genome}_vs_ref.csv")
        thr_csv = os.path.join(self.SUBDIRS["thr"], f"{genome}_vs_ref-threshold.csv")
        best_csv = os.path.join(self.SUBDIRS["best"], f"{genome}_vs_ref-threshold-best.csv")
        clean_csv = os.path.join(self.SUBDIRS["clean"], f"{genome}_vs_ref-threshold-best-clean.csv")

        # ------------------------------------------------------------------
        # BLASTp
        # ------------------------------------------------------------------
        if not os.path.exists(out_csv):
            subprocess.run(
                [
                    "blastp",
                    "-query", faa_file,
                    "-db", self.REF_DB,
                    "-outfmt", self.OUTFMT,
                    "-evalue", str(self.E_VALUE),
                    "-num_threads", str(self.THREADS_PER_BLAST),
                    "-out", out_csv
                ],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=True
            )

        # Load BLAST result
        df = pd.read_csv(out_csv, header=None)
        df.columns = [
            "qseqid", "qlen", "sseqid", "slen", "pident", "qcovs",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"
        ]

        # ------------------------------------------------------------------
        # Threshold filtering
        # ------------------------------------------------------------------
        df_thr = df[(df["pident"] >= self.IDENTITY_THRESHOLD) &
                    (df["qcovs"] >= self.COVERAGE_THRESHOLD)]
        df_thr.to_csv(thr_csv, index=False)

        if df_thr.empty:
            df_thr.to_csv(best_csv, index=False)
            df_thr.to_csv(clean_csv, index=False)
            return genome

        # Select best hit per query (lowest evalue)
        best = (
            df_thr.sort_values("evalue")
                  .groupby("qseqid", as_index=False)
                  .first()
        )
        best.to_csv(best_csv, index=False)

        # ------------------------------------------------------------------
        # Reference ID filtering + gene name extraction
        # ------------------------------------------------------------------
        _, ref_set = self.load_ref_ids(self.REF_FASTA)
        best["sseqid"] = best["sseqid"].apply(self.clean_sseqid)
        clean = best[best["sseqid"].isin(ref_set)].copy()

        # Gene name extraction (from refseq format: "gene@id")
        if self.ALLOW_MULTI_REF:
            clean["gene"] = clean["sseqid"].apply(
                lambda x: x.split("@")[0] if "@" in x else x
            )
            clean = clean[["qseqid", "qlen", "gene"] +
                          [c for c in clean.columns
                           if c not in ("qseqid", "qlen", "gene")]]

        # Keep top 10 bitscore per gene
        if "gene" in clean.columns:
            clean = (
                clean.sort_values("bitscore", ascending=False)
                     .groupby("gene", as_index=False)
                     .head(10)
            )

        clean.to_csv(clean_csv, index=False)

        # Export sequences
        self.extract_hits_to_fasta(
            genome,
            faa_file,
            clean,
            self.SUBDIRS["fasta"],
            add_prefix=self.ADD_GENOMEPREFIX
        )

        with open(self.DONE_FILE, "a") as f:
            f.write(genome + "\n")

        return genome


    # --------------------------------------------------------------
    # Merge sequences gene-wise → FastaMerged/
    # --------------------------------------------------------------

    def _merge_fasta_by_gene(self):
        """
        Merge all gene-specific FASTA files from all genomes into:
            FastaMerged/{gene}.fasta
        """
        gene_fastas = glob(os.path.join(self.SUBDIRS["fasta"], "*.fasta"))
        gene_map = {}

        for fa in gene_fastas:
            fname = os.path.basename(fa)
            gene = fname.split("_x_")[0]
            gene_map.setdefault(gene, []).append(fa)

        total = len(gene_map)
        if total == 0:
            return

        step_start = time.time()
        idx = 0

        for gene, paths in gene_map.items():
            idx += 1
            merged = os.path.join(self.SUBDIRS["merged"], f"{gene}.fasta")

            with open(merged, "w") as out:
                for p in paths:
                    for rec in SeqIO.parse(p, "fasta"):
                        SeqIO.write(rec, out, "fasta")

            if self.logger:
                self.logger.progress(
                    "Script2",
                    "GeneMerge",
                    idx,
                    total,
                    step_start
                )


    # --------------------------------------------------------------
    # Build Genome × Gene matrix
    # --------------------------------------------------------------

    def _build_statistics_matrix(self):
        """
        Build Genome × Gene matrix from cleaned BLAST results.
        Two output matrices:
            BLASTp_Statistics.csv         — raw counts
            BLASTp_Statistics_Normalized.csv — presence/absence (0/1)
        """
        clean_files = glob(os.path.join(self.SUBDIRS["clean"],
                                        "*_vs_ref-threshold-best-clean.csv"))

        records = []
        for f in clean_files:
            genome = os.path.basename(f).replace("_vs_ref-threshold-best-clean.csv", "")
            df = pd.read_csv(f)
            if df.empty:
                continue

            col = "gene" if "gene" in df.columns else "sseqid"
            for g in df[col]:
                records.append((genome, g))

        if not records:
            return

        df_all = pd.DataFrame(records, columns=["Genome", "Gene"])

        matrix = (
            df_all.groupby(["Genome", "Gene"])
                  .size()
                  .unstack(fill_value=0)
                  .astype(int)
        )

        matrix_path = os.path.join(self.RESULTS_DIR, "BLASTp_Statistics.csv")
        matrix_norm_path = os.path.join(self.RESULTS_DIR, "BLASTp_Statistics_Normalized.csv")

        matrix.to_csv(matrix_path, encoding="utf-8-sig")
        (matrix > 0).astype(int).to_csv(matrix_norm_path, encoding="utf-8-sig")


    # --------------------------------------------------------------
    # run() — main entry point of Script2
    # --------------------------------------------------------------

    def run(self):
        """
        Main routine for Script2:
            - Create BLAST database
            - Locate target FASTA files
            - Run BLASTp in parallel
            - Build statistics matrix
            - Merge gene-wise FASTA
        """
        self.make_blastdb_if_needed()

        # Step2 scans only CDHIT_Target directory
        faa_files = []
        for ext in ("*.fasta", "*.fa", "*.faa", "*.fna"):
            faa_files.extend(glob(os.path.join(self.QUERY_DIR, "**", ext),
                                  recursive=True))

        faa_files = sorted(f for f in faa_files if os.path.isfile(f))

        total = len(faa_files)
        if total == 0:
            if self.logger:
                self.logger.write(
                    "[Script2] Warning: No target sequences found in CDHIT_Target (*.faa/fa/fasta/fna)"
                )
            return

        if self.logger:
            self.logger.write(
                f"[Script2] Total BLAST tasks: {total} genomes (CD-HIT filtered)"
            )

        step_start = time.time()
        completed = 0

        # Parallel BLAST execution
        with ThreadPoolExecutor(max_workers=self.PARALLEL_JOBS) as exe:
            future_map = {
                exe.submit(self._process_one_genome, faa): faa
                for faa in faa_files
            }

            for fut in as_completed(future_map):
                _ = fut.result()
                completed += 1

                if self.logger:
                    self.logger.progress(
                        "Script2",
                        "BLASTp",
                        completed,
                        total,
                        step_start
                    )

        # Build Genome × Gene matrix
        self._build_statistics_matrix()

        # Merge sequences into gene-level FASTA files
        self._merge_fasta_by_gene()
# ==============================================================
# Script3: PfamGenomeScreen — Pfam HMM-based refinement
# ==============================================================

class PfamGenomeScreen:
    """
    Script3: Domain-level HMM validation and refinement.

    Purpose:
    After BLAST-based extraction (Step2), each gene family may contain
    false positives or sequences that only partially match the target domain.
    This module performs a second round of filtering using Pfam/TIGRFAM
    HMMs (via hmmsearch), based on:

        - Trusted cutoff (--cut_tc)
        - Domain coverage filtering (match_length / hmm_length ≥ threshold)
        - e-value threshold

    A configuration file named:
        Pfam_hmmserach-config.txt
    located under ref_hmm_dir/
    defines which HMM(s) belong to each gene.

    Format (strict two-column):
        gene_fasta_name    hmm1,hmm2,hmm3
    Example:
        rplB.fasta          PF12345.hmm,PF67890.hmm

    Only sequences that pass ANY of the assigned HMMs
    (union mode) will be retained.
    """

    def __init__(self,
                 fasta_merged_dir,
                 ref_hmm_dir,
                 working_dir,
                 cpu=4,
                 evalue=1e-5,
                 domain_cov=0.5,
                 logger=None):

        # Directories
        self.fasta_merged_dir = os.path.abspath(fasta_merged_dir)
        self.ref_hmm_dir = os.path.abspath(ref_hmm_dir)
        self.working_dir = os.path.abspath(working_dir)

        # Parameters
        self.CPU = cpu
        self.EVALUE = evalue
        self.DOMAIN_COV = domain_cov
        self.logger = logger

        # Output directory
        self.result_dir = os.path.join(self.working_dir,
                                       "Pfam_hmmserach-config_Results")
        os.makedirs(self.result_dir, exist_ok=True)

        # Config file (two-column mapping)
        self.config_path = os.path.join(self.ref_hmm_dir,
                                        "Pfam_hmmserach-config.txt")

        # Mapping: gene → list of HMM paths
        self.mapping = self.load_config(self.config_path)


    # --------------------------------------------------------------
    # Load Pfam_hmmserach-config (two-column)
    # --------------------------------------------------------------

    def load_config(self, path):
        """
        Parse the Pfam_hmmserach-config.txt file.

        Expected format:
            <gene_fasta_name>   <comma-separated HMM names>

        Example:
            rplB.fasta  PF00203.hmm,PF03947.hmm
        """
        if not os.path.isfile(path):
            raise RuntimeError(f"[Script3] Config file not found: {path}")

        mapping = {}
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                if not line.strip():
                    continue

                parts = line.strip().split()
                if len(parts) != 2:
                    raise RuntimeError(
                        f"[Script3] Invalid config format (two columns required): {line}"
                    )

                gene_fa, hmm_field = parts
                gene = os.path.splitext(gene_fa)[0]

                hmm_names = [x.strip() for x in hmm_field.split(",") if x.strip()]
                if not hmm_names:
                    raise RuntimeError(f"[Script3] No HMMs specified for: {gene_fa}")

                hmm_paths = []
                for hmm_name in hmm_names:
                    hmm_path = os.path.join(self.ref_hmm_dir, hmm_name)
                    if not os.path.isfile(hmm_path):
                        raise RuntimeError(f"[Script3] HMM not found: {hmm_path}")
                    hmm_paths.append(hmm_path)

                mapping[gene] = hmm_paths

        return mapping


    # --------------------------------------------------------------
    # domtblout parser — extract passing hits
    # --------------------------------------------------------------

    def parse_domtblout(self, domtblout_path):
        """
        Parse HMMER domtblout and extract IDs that pass:
            - E-value threshold
            - Domain coverage ≥ DOMAIN_COV

        Domain coverage = (hmm_to - hmm_from + 1) / hmm_length.
        """
        hits = set()
        if not os.path.isfile(domtblout_path):
            return hits

        with open(domtblout_path, "r") as f:
            for line in f:
                if not line or line.startswith("#"):
                    continue

                parts = line.strip().split()
                if len(parts) < 23:
                    continue

                try:
                    # Indices based on HMMER domtblout format
                    evalue = float(parts[6])
                    hmm_len = int(parts[5])
                    hmm_from = int(parts[15])
                    hmm_to = int(parts[16])

                    match_len = hmm_to - hmm_from + 1
                    cov = match_len / hmm_len if hmm_len > 0 else 0.0

                    if evalue <= self.EVALUE and cov >= self.DOMAIN_COV:
                        hits.add(parts[0])  # target sequence ID
                except Exception:
                    continue

        return hits


    # --------------------------------------------------------------
    # Process one gene family (union mode across all assigned HMMs)
    # --------------------------------------------------------------

    def _process_one_gene(self, gene):
        """
        For each gene family:
            - Run hmmsearch for all assigned HMMs
            - Parse domtblout results
            - Combine hits (union mode)
            - Write filtered sequences to result_dir/{gene}.fasta
        """
        gene_fasta = os.path.join(self.fasta_merged_dir, f"{gene}.fasta")
        if not os.path.isfile(gene_fasta):
            if self.logger:
                self.logger.write(
                    f"[Script3] Warning: skipping {gene}, file not found: {gene}.fasta"
                )
            return

        hmm_files = self.mapping.get(gene, [])
        if not hmm_files:
            if self.logger:
                self.logger.write(f"[Script3] Warning: no HMMs configured for {gene}")
            return

        tmp_dir = os.path.join(self.working_dir, f"PfamTmp_{gene}")
        os.makedirs(tmp_dir, exist_ok=True)

        union_hits = set()

        # Run hmmsearch for each HMM
        for hmm_file in hmm_files:
            hmm_basename = os.path.basename(hmm_file)
            domtblout_path = os.path.join(
                tmp_dir, f"{gene}_{hmm_basename}.domtblout"
            )

            cmd = [
                "hmmsearch",
                "--cpu", str(self.CPU),
                "--domtblout", domtblout_path,
                "--cut_tc",
                hmm_file,
                gene_fasta
            ]

            subprocess.run(
                cmd,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=True
            )

            hits = self.parse_domtblout(domtblout_path)
            union_hits |= hits

        # Write filtered sequences
        out_fa = os.path.join(self.result_dir, f"{gene}.fasta")
        with open(out_fa, "w") as fw:
            for rec in SeqIO.parse(gene_fasta, "fasta"):
                if rec.id in union_hits:
                    SeqIO.write(rec, fw, "fasta")

        # Remove temporary files
        shutil.rmtree(tmp_dir, ignore_errors=True)


    # --------------------------------------------------------------
    # run() — main entry point of Script3
    # --------------------------------------------------------------

    def run(self):
        """
        Main routine:
            - Load config mapping
            - For each gene, run HMM-based refinement
            - Write results to result_dir
        """
        genes = sorted(self.mapping.keys())
        total_genes = len(genes)

        if total_genes == 0:
            if self.logger:
                self.logger.write("[Script3] Warning: no genes in config file")
            return

        if self.logger:
            self.logger.write(
                f"[Script3] Pfam HMM refinement (union mode): {total_genes} gene families"
            )

        step_start = time.time()

        for idx, gene in enumerate(genes, start=1):
            self._process_one_gene(gene)

            if self.logger:
                self.logger.progress(
                    "Script3",
                    "HMMER_Gene",
                    idx,
                    total_genes,
                    step_start
                )

        if self.logger:
            self.logger.write(
                f"[Script3] Completed. Results written to {self.result_dir}"
            )
# ==============================================================
# Script4: AutoMSAConcatenation
# Two-stage MSA → trimming → phylogeny → pruning → concatenation
# ==============================================================

class AutoMSAConcatenation:
    """
    Script4:

    This module performs two-stage multiple sequence alignment (MSA),
    trimming, tree reconstruction, phylogenetic pruning (phylopypruner),
    and supermatrix concatenation.

    Workflow overview:

    Primary stage:
        1. Remove terminal stop codons
        2. MAFFT alignment (global, fast or auto mode)
        3. TrimAl trimming
        4. FastTree inference
        5. Keep only valid gene families where aln + tree pair is complete
        6. phylopypruner filtering

    Secondary stage:
        1. MAFFT --localpair (more accurate realignment)
        2. TrimAl trimming
        3. Concatenation into supermatrix
        4. FastTree tree for the final concatenated alignment
        5. Generate:
               - Gene_Count_Matrix.csv
               - partition.txt
               - statistics.csv
               - {prefix}-concatenation_bygene.nex

    Notes:
        - IQ-TREE has been removed; only FastTree is used.
        - All intermediate and final results are stored under:
              fasta_dir/Auto_MSA_Tree_Results/
        - All concatenated outputs are copied to the supermatrix directory
          in Step4 of the main pipeline.
    """

    def __init__(self, fasta_dir, logger=None):

        self.fasta_merged_dir = os.path.abspath(fasta_dir)
        self.logger = logger

        # Parallelization configuration
        self.THREADS = 5
        self.MAX_PARALLEL_FASTA = 4
        self.MAX_PARALLEL_SECONDARY = 4

        # Primary-stage parameters
        self.MAFFT_MODE = "auto"
        self.TRIMAL_MODE = "auto"

        # Secondary-stage parameters
        self.Secondary_MAFFT_MODE = "localpair"
        self.Secondary_TRIMAL_MODE = "gt"

        # phylopypruner
        self.PHYLOPYP_RUN = True

        # FastTree parameters
        self.FASTTREE_ARGS = "-lg -gamma -boot 1000"

        # Output directory structure
        self.results_root = os.path.join(self.fasta_merged_dir,
                                         "Auto_MSA_Tree_Results")
        self.tree_dir = os.path.join(self.results_root, "Tree_files")
        self.inter_dir = os.path.join(self.results_root, "MSA_Tree_Intermediate")
        self.trimmed_dir = os.path.join(self.results_root, "Trimed_fasta_files")
        self.fasta_tree_dir = os.path.join(self.results_root, "fasta_tree")

        for d in [
            self.results_root,
            self.tree_dir,
            self.inter_dir,
            self.trimmed_dir,
            self.fasta_tree_dir
        ]:
            os.makedirs(d, exist_ok=True)


    # ==============================================================
    # Utility functions
    # ==============================================================

    @staticmethod
    def run_cmd(cmd):
        """
        Execute shell command silently (stdout/stderr suppressed).
        """
        subprocess.run(cmd, shell=True, check=True,
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)

    @staticmethod
    def strip_terminal_stop(in_fa, out_fa):
        """
        Remove trailing '*' (terminal stop codons) from amino acid sequences.

        Input FASTA → Output FASTA with all terminal '*' removed.
        """
        seq = ""
        hdr = None

        def flush(h, s, fout):
            if h is None:
                return
            s = s.rstrip("*")
            if not s:
                return
            fout.write(h + "\n")
            for i in range(0, len(s), 60):
                fout.write(s[i:i+60] + "\n")

        with open(in_fa) as fin, open(out_fa, "w") as fout:
            for line in fin:
                L = line.strip()
                if not L:
                    continue

                if L.startswith(">"):
                    flush(hdr, seq, fout)
                    hdr = L
                    seq = ""
                else:
                    seq += L

            flush(hdr, seq, fout)

        return out_fa


    @staticmethod
    def detect_fasttree():
        """
        Locate FastTree executable:
            FastTreeMP > fasttreemp > FastTree > fasttree

        Raises error if neither is found.
        """
        for name in ["FastTreeMP", "fasttreemp", "FastTree", "fasttree"]:
            if shutil.which(name):
                return name
        raise RuntimeError("[Script4] FastTree/FastTreeMP not found.")


    # ==============================================================
    # FastTree input requirements
    # ==============================================================

    def fasttree_input_valid(self, trimmed_fasta, gene):
        """
        Validate alignment before FastTree usage:

        Criteria:
            - Number of sequences ≥ 4
            - No sequence of length 0
            - All sequences must have identical aligned length
        """
        seqs = list(SeqIO.parse(trimmed_fasta, "fasta"))
        n = len(seqs)

        if n < 4:
            return False, f"Sequence count {n} < 4"

        lengths = [len(str(rec.seq)) for rec in seqs]
        if any(L == 0 for L in lengths):
            return False, "Contains empty sequences"

        if len(set(lengths)) != 1:
            return False, "Alignment lengths inconsistent"

        return True, ""


    # ==============================================================
    # Generic parallel execution engine
    # ==============================================================

    def run_stage_parallel(self, label, tasks, worker, max_workers, script_label):
        """
        Parallel execution helper for all primary and secondary stage tasks.

        Parameters:
            label:          Step label (e.g. "Primary MAFFT")
            tasks:          List of tasks (tuples)
            worker:         Function to process one task
            max_workers:    Number of parallel threads
            script_label:   Prefix label (e.g. "Script4")
        """
        total = len(tasks)
        if total == 0:
            if self.logger:
                self.logger.write(f"[{script_label}] [{label}] No tasks.")
            return []

        if self.logger:
            self.logger.write(f"[{script_label}] [{label}] Starting: {total} tasks")

        results = []
        completed = 0
        start = time.time()

        with ThreadPoolExecutor(max_workers=max_workers) as exe:
            fut_map = {exe.submit(worker, t): t for t in tasks}

            for fut in as_completed(fut_map):
                res = fut.result()
                results.append(res)
                completed += 1

                if self.logger:
                    self.logger.progress(
                        script_label,
                        label,
                        completed,
                        total,
                        start
                    )

        return results


    # ==============================================================
    # Primary stage: MAFFT
    # ==============================================================

    def worker_mafft_primary(self, task):
        """
        Worker function for primary MAFFT alignment.
        """
        gene, in_fa, gene_dir = task
        aln = os.path.join(gene_dir, f"{gene}.aln.fasta")

        cmd = (
            f"mafft --thread {self.THREADS} --anysymbol --reorder "
            f"{in_fa} > {aln}"
        )
        self.run_cmd(cmd)
        return gene


    # ==============================================================
    # Primary stage: TrimAl
    # ==============================================================

    def worker_trimal_primary(self, task):
        """
        Worker for primary TrimAl trimming.
        """
        gene, gene_dir = task
        aln = os.path.join(gene_dir, f"{gene}.aln.fasta")
        trimmed = os.path.join(gene_dir, f"{gene}.trimmed.fasta")

        if self.TRIMAL_MODE == "auto":
            cmd = f"trimal -automated1 -in {aln} -out {trimmed}"
        else:
            cmd = f"trimal -in {aln} -out {trimmed} -gt 0.2"

        self.run_cmd(cmd)
        return gene


    # ==============================================================
    # Primary stage: FastTree
    # ==============================================================

    def worker_tree_primary(self, task):
        """
        Worker for primary FastTree phylogeny:
            trimmed alignment → FastTree tree
        """
        gene, gene_dir = task
        trimmed = os.path.join(gene_dir, f"{gene}.trimmed.fasta")

        ok, reason = self.fasttree_input_valid(trimmed, gene)
        if not ok:
            if self.logger:
                self.logger.write(f"[Script4] Skipping {gene}: {reason}")
            return gene

        ft = self.detect_fasttree()
        treefile = os.path.join(gene_dir, f"{gene}.fasttree.tre")

        cmd = f"{ft} {self.FASTTREE_ARGS} < {trimmed} > {treefile}"
        self.run_cmd(cmd)
        return gene
    # ==============================================================
    # Secondary stage: MAFFT (localpair)
    # ==============================================================

    def worker_secondary_mafft(self, task):
        """
        Worker for secondary MAFFT alignment using --localpair and
        high-iteration refinement, typically on pruned alignments
        from phylopypruner.
        """
        gene, in_fa, out_fa = task
        cmd = (
            f"mafft --anysymbol --reorder --localpair --maxiterate 1000 "
            f"--thread {self.THREADS} {in_fa} > {out_fa}"
        )
        self.run_cmd(cmd)
        return gene


    # ==============================================================
    # Secondary stage: TrimAl
    # ==============================================================

    def worker_secondary_trimal(self, task):
        """
        Worker for secondary TrimAl trimming.
        """
        gene, in_fa, out_fa = task

        if self.Secondary_TRIMAL_MODE == "gt":
            cmd = f"trimal -in {in_fa} -out {out_fa} -gt 0.2"
        else:
            cmd = f"trimal -automated1 -in {in_fa} -out {out_fa}"

        self.run_cmd(cmd)
        return gene


    # ==============================================================
    # Primary stage controller:
    # MAFFT → TrimAl → FastTree → phylopypruner
    # ==============================================================

    def run_primary_stage(self):
        """
        Run the primary pipeline stage:

            1. Filter gene families with ≤ 1 sequence.
            2. Remove terminal stop codons.
            3. Primary MAFFT alignment.
            4. Primary TrimAl trimming.
            5. FastTree phylogeny.
            6. phylopypruner on (aln + tree) pairs.

        Returns:
            List of genes that successfully passed primary stage
            and entered phylopypruner.
        """
        fasta_files = [
            f for f in os.listdir(self.fasta_merged_dir)
            if f.lower().endswith((".fasta", ".fa", ".faa"))
        ]

        genes = []
        clean_fastas = {}

        # Remove terminal stop codons and prepare per-gene directories
        for f in fasta_files:
            gene = os.path.splitext(f)[0]
            in_fa = os.path.join(self.fasta_merged_dir, f)

            seq_count = sum(1 for _ in SeqIO.parse(in_fa, "fasta"))
            if seq_count <= 1:
                if self.logger:
                    self.logger.write(
                        f"[Script4] Skipping {gene}: sequence count = {seq_count}"
                    )
                continue

            out_fa = os.path.join(self.results_root, f"{gene}.nostop.fasta")
            self.strip_terminal_stop(in_fa, out_fa)

            gene_dir = os.path.join(self.inter_dir, gene)
            os.makedirs(gene_dir, exist_ok=True)

            genes.append(gene)
            clean_fastas[gene] = (out_fa, gene_dir)

        if not genes:
            return []

        # ---------------- Primary MAFFT ----------------
        mafft_tasks = [
            (g, clean_fastas[g][0], clean_fastas[g][1]) for g in genes
        ]
        self.run_stage_parallel(
            "Primary MAFFT",
            mafft_tasks,
            self.worker_mafft_primary,
            self.MAX_PARALLEL_FASTA,
            "Script4"
        )

        # ---------------- Primary TrimAl ----------------
        trimal_tasks = [(g, clean_fastas[g][1]) for g in genes]
        self.run_stage_parallel(
            "Primary TrimAl",
            trimal_tasks,
            self.worker_trimal_primary,
            self.MAX_PARALLEL_FASTA,
            "Script4"
        )

        # ---------------- Primary FastTree ----------------
        tree_tasks = [(g, clean_fastas[g][1]) for g in genes]
        self.run_stage_parallel(
            "Primary FastTree",
            tree_tasks,
            self.worker_tree_primary,
            self.MAX_PARALLEL_FASTA,
            "Script4"
        )

        # --------------------------------------------------
        # Critical logic:
        # Only gene families with both aln and tree files
        # proceed to phylopypruner and secondary stage.
        # --------------------------------------------------
        valid_pairs = []

        for g in genes:
            gene_dir = clean_fastas[g][1]

            aln = os.path.join(gene_dir, f"{g}.aln.fasta")
            treefile = os.path.join(gene_dir, f"{g}.fasttree.tre")

            if os.path.isfile(aln) and os.path.isfile(treefile):
                shutil.copy(aln, os.path.join(self.fasta_tree_dir, f"{g}.aln.fasta"))
                shutil.copy(treefile, os.path.join(self.fasta_tree_dir, f"{g}.tre"))
                shutil.copy(treefile, os.path.join(self.tree_dir, f"{g}.tre"))

                valid_pairs.append(g)
            else:
                if self.logger:
                    self.logger.write(
                        f"[Script4] Skipping {g}: missing aln/tree pair, "
                        f"not passed to phylopypruner."
                    )

        if self.logger:
            self.logger.write(
                f"[Script4] Primary valid genes with aln+tree pairs: {len(valid_pairs)}"
            )

        # ---------------- phylopypruner ----------------
        if self.PHYLOPYP_RUN and valid_pairs:
            abs_path = os.path.abspath(self.fasta_tree_dir)
            cmd = f"QT_QPA_PLATFORM=offscreen phylopypruner --dir {abs_path}"
            self.run_cmd(cmd)

        return valid_pairs


    # ==============================================================
    # Secondary stage:
    #   Load phylopypruner output → secondary MAFFT/TrimAl → concat
    # ==============================================================

    def run_secondary_stage(self):
        """
        Run secondary alignment and concatenation:

            - Read pruned alignments from phylopypruner
            - Secondary MAFFT (--localpair)
            - Secondary TrimAl
            - Build supermatrix
            - Build FastTree tree for supermatrix
            - Generate per-gene statistics, partition file, and Nexus
        """
        pruned_dir = os.path.join(
            self.fasta_tree_dir,
            "phylopypruner_output",
            "output_alignments"
        )

        if not os.path.isdir(pruned_dir):
            if self.logger:
                self.logger.write(
                    "[Script4] Warning: phylopypruner output not found."
                )
            return

        # Secondary directory structure
        sec_root = os.path.join(self.fasta_merged_dir, "Secondary_MAFFT_TRIMAL")
        mafft_dir = os.path.join(sec_root, "Mafft")
        trim_dir = os.path.join(sec_root, "Trim")
        conc_dir = os.path.join(sec_root, "Concatenation")

        for d in [sec_root, mafft_dir, trim_dir, conc_dir]:
            os.makedirs(d, exist_ok=True)

        # Collect pruned alignments
        pruned_files = [
            f for f in os.listdir(pruned_dir)
            if f.endswith(".aln_pruned.fasta")
        ]
        if not pruned_files:
            return

        gene_list = []
        sec_mafft_tasks = []

        for f in pruned_files:
            gene = f.replace(".aln_pruned.fasta", "")
            gene_list.append(gene)

            in_fa = os.path.join(pruned_dir, f)
            out_fa = os.path.join(mafft_dir, f"{gene}.m.fasta")

            sec_mafft_tasks.append((gene, in_fa, out_fa))

        # ---------------- Secondary MAFFT (localpair) ----------------
        self.run_stage_parallel(
            "Secondary MAFFT",
            sec_mafft_tasks,
            self.worker_secondary_mafft,
            self.MAX_PARALLEL_SECONDARY,
            "Script4"
        )

        # ---------------- Secondary TrimAl ----------------
        sec_trimal_tasks = []
        for gene in gene_list:
            in_fa = os.path.join(mafft_dir, f"{gene}.m.fasta")
            out_fa = os.path.join(trim_dir, f"{gene}.m_trim.fasta")
            sec_trimal_tasks.append((gene, in_fa, out_fa))

        self.run_stage_parallel(
            "Secondary TrimAl",
            sec_trimal_tasks,
            self.worker_secondary_trimal,
            self.MAX_PARALLEL_SECONDARY,
            "Script4"
        )

        # ==================================================
        # Build concatenated supermatrix
        # ==================================================
        gene_to_seq = {}
        genomes = set()

        # Parse each secondary-trimmed gene alignment
        for gene in gene_list:
            fa = os.path.join(trim_dir, f"{gene}.m_trim.fasta")
            if not os.path.isfile(fa):
                continue

            gene_to_seq[gene] = {}

            hdr = None
            seq = ""
            with open(fa) as f:
                for line in f:
                    L = line.strip()
                    if not L:
                        continue
                    if L.startswith(">"):
                        if hdr:
                            g = hdr[1:].split("@")[0]
                            genomes.add(g)
                            gene_to_seq[gene][g] = seq
                        hdr = L
                        seq = ""
                    else:
                        seq += L

            if hdr:
                g = hdr[1:].split("@")[0]
                genomes.add(g)
                gene_to_seq[gene][g] = seq

        genomes = sorted(genomes)
        if not genomes:
            return

        # Alignment length for each gene
        gene_len = {
            gene: (
                len(next(iter(gene_to_seq[gene].values())))
                if gene_to_seq.get(gene) else 0
            )
            for gene in gene_list
        }

        # --------------------------------------------------
        # Supermatrix output: use global prefix_concatenation
        # --------------------------------------------------
        concat_fa = os.path.join(
            conc_dir,
            f"{prefix_concatenation}-concatenation.fasta"
        )

        with open(concat_fa, "w") as fw:
            for g in genomes:
                fw.write(f">{g}\n")
                parts = []
                for gene in gene_list:
                    L = gene_len.get(gene, 0)
                    if L == 0:
                        continue
                    if g in gene_to_seq.get(gene, {}):
                        parts.append(gene_to_seq[gene][g])
                    else:
                        parts.append("?" * L)
                fw.write("".join(parts) + "\n")

        # ==================================================
        # FastTree for concatenated supermatrix (FastTree-only)
        # ==================================================
        ft = self.detect_fasttree()
        concat_tree = os.path.join(
            conc_dir,
            f"{prefix_concatenation}-concatenation.fasttree-ref.tre"
        )

        cmd = f"{ft} {self.FASTTREE_ARGS} < {concat_fa} > {concat_tree}"
        self.run_cmd(cmd)

        # ==================================================
        # Generate Gene_Count_Matrix, partition, statistics, and Nexus
        # ==================================================
        self._write_gene_count_matrix(conc_dir, genomes, gene_list, gene_to_seq)
        self._write_partition_and_statistics(
            conc_dir, genomes, gene_list, gene_len, gene_to_seq
        )
        self._write_nexus_by_gene(
            conc_dir, concat_fa, genomes, gene_list, gene_len
        )

        if self.logger:
            self.logger.write(
                "[Script4] Secondary stage completed. "
                "All supermatrix files have been generated."
            )


    # ==============================================================
    # Gene_Count_Matrix.csv
    # ==============================================================

    def _write_gene_count_matrix(self, conc_dir, genomes, gene_list, gene_to_seq):
        """
        Write presence/absence matrix (1/0) of genes per genome:
            Gene_Count_Matrix.csv
        """
        out_path = os.path.join(conc_dir, "Gene_Count_Matrix.csv")
        with open(out_path, "w", newline="") as fw:
            writer = csv.writer(fw)
            writer.writerow(["Genome"] + list(gene_list))

            for g in genomes:
                row = [g]
                for gene in gene_list:
                    present = 1 if g in gene_to_seq.get(gene, {}) else 0
                    row.append(present)
                writer.writerow(row)


    # ==============================================================
    # partition.txt + statistics.csv
    # ==============================================================

    def _write_partition_and_statistics(self,
                                        conc_dir,
                                        genomes,
                                        gene_list,
                                        gene_len,
                                        gene_to_seq):
        """
        Generate:
            - partition.txt  (per-gene alignment blocks)
            - statistics.csv (per-gene occupancy and coverage)
        """
        part_path = os.path.join(conc_dir, "partition.txt")
        start = 1

        with open(part_path, "w") as fw:
            for gene in gene_list:
                L = gene_len.get(gene, 0)
                if L <= 0:
                    continue
                end = start + L - 1
                fw.write(f"{gene} = {start}-{end}\n")
                start = end + 1

        stat_path = os.path.join(conc_dir, "statistics.csv")
        n_total = len(genomes)

        with open(stat_path, "w", newline="") as fw:
            writer = csv.writer(fw)

            writer.writerow([
                "Gene", "Alignment_Length",
                "Num_Genomes_Total",
                "Num_Genomes_Present",
                "Num_Genomes_Missing",
                "Frac_Genomes_Present",
                "Frac_Genomes_Missing"
            ])

            for gene in gene_list:
                L = gene_len.get(gene, 0)
                present = len(gene_to_seq.get(gene, {}))
                missing = n_total - present
                frac_present = present / n_total if n_total else 0
                frac_missing = missing / n_total if n_total else 0

                writer.writerow([
                    gene, L,
                    n_total, present, missing,
                    f"{frac_present:.4f}",
                    f"{frac_missing:.4f}"
                ])


    # ==============================================================
    # concatenation_bygene.nex (Nexus with gene charsets)
    # ==============================================================

    def _write_nexus_by_gene(self,
                             conc_dir,
                             concat_fa,
                             genomes,
                             gene_list,
                             gene_len):
        """
        Generate Nexus file with per-gene charsets:
            {prefix_concatenation}-concatenation_bygene.nex
        """
        nex_path = os.path.join(
            conc_dir,
            f"{prefix_concatenation}-concatenation_bygene.nex"
        )

        taxon_to_seq = {}
        hdr = None
        seq = ""

        with open(concat_fa) as f:
            for line in f:
                L = line.strip()
                if not L:
                    continue
                if L.startswith(">"):
                    if hdr:
                        taxon_to_seq[hdr[1:]] = seq
                    hdr = L
                    seq = ""
                else:
                    seq += L

            if hdr:
                taxon_to_seq[hdr[1:]] = seq

        taxa = genomes
        nchar = len(taxon_to_seq[taxa[0]])

        with open(nex_path, "w") as fw:
            fw.write("#NEXUS\n\n")
            fw.write("begin data;\n")
            fw.write(f" dimensions ntax={len(taxa)} nchar={nchar};\n")
            fw.write(" format datatype=protein gap=- missing=?;\n")
            fw.write(" matrix\n")

            for t in taxa:
                fw.write(f"{t}    {taxon_to_seq[t]}\n")

            fw.write(";\nend;\n\n")

            fw.write("begin sets;\n")
            start = 1
            for gene in gene_list:
                L = gene_len.get(gene, 0)
                if L > 0:
                    end = start + L - 1
                    fw.write(f" charset {gene} = {start}-{end};\n")
                    start = end + 1
            fw.write("end;\n")


    # ==============================================================
    # run() — primary + secondary stage controller
    # ==============================================================

    def run(self):
        """
        Entry point for Script4:

            1. Run primary stage:
                - nostop FASTA
                - primary MSA / TrimAl / FastTree
                - phylopypruner

            2. Run secondary stage:
                - secondary MSA / TrimAl
                - concatenation
                - gene matrix, partition, statistics, Nexus
        """
        genes = self.run_primary_stage()
        if genes:
            self.run_secondary_stage()
# ==============================================================
# Pipeline: SupermatrixPipeline — Main workflow controller
# ==============================================================

class SupermatrixPipeline:
    """
    Main pipeline controller orchestrating all processing steps:

        Step0: CD-HIT filtering of target genome protein FASTA files
        Step1: Reference FASTA → QC → CD-HIT representative sequences
        Step2: BLASTp screening and extraction (AutoBlastPExtractor)
        Step3: Domain-level refinement using Pfam/TIGRFAM HMMs
        Step4: Two-stage MSA, trimming, tree inference, pruning, and
               concatenated supermatrix generation

    After Step4:
        - Clean up temporary files and intermediate directories
        - Finalize and save log file into supermatrix output directory
    """

    def __init__(self):
        self.base_dir = os.path.abspath(os.getcwd())

        # User-defined directories
        self.taget_dir_name = taget_dir
        self.ref_seq_dir_name = ref_seq_dir
        self.ref_hmm_dir_name = ref_hmm_dir
        self.supermatrix_outdir_name = supermatrix_outdir

        # Supermatrix output directory
        self.supermatrix_outdir = os.path.join(
            self.base_dir,
            self.supermatrix_outdir_name
        )
        os.makedirs(self.supermatrix_outdir, exist_ok=True)

        # Working directory
        self.working_dir = os.path.join(
            self.supermatrix_outdir,
            "Working_Dir"
        )
        os.makedirs(self.working_dir, exist_ok=True)

        # CD-HIT target directory (Step0 output)
        self.cdhit_target_dir = os.path.join(
            self.supermatrix_outdir,
            CDHIT_TARGET_DIRNAME
        )
        os.makedirs(self.cdhit_target_dir, exist_ok=True)

        # Absolute paths for target and reference directories
        self.taget_dir = os.path.join(self.base_dir, self.taget_dir_name)
        self.ref_seq_dir = os.path.join(self.base_dir, self.ref_seq_dir_name)
        self.ref_hmm_dir = os.path.join(self.base_dir, self.ref_hmm_dir_name)

        # Logging
        self.pipeline_log_tmp = os.path.join(
            self.working_dir,
            "Supermatrix_Pipeline.log.tmp"
        )
        self.logger = PipelineLogger(self.pipeline_log_tmp)

        # Step1 output paths
        self.ref_fasta = os.path.join(self.working_dir, "ref.fasta")
        self.ref_db_prefix = os.path.join(self.working_dir, "ref")


    # --------------------------------------------------------------
    # Step0: CD-HIT (c = 0.98) for target genome protein sequences
    # --------------------------------------------------------------

    def prepare_target_cdhit(self):
        """
        Step0:
            Perform CD-HIT (c = 0.98) on all protein FASTA files
            under target_dir. Results are placed into:

                Working_Dir/CDHIT_Target/

            Original genome files remain unchanged.
        """
        self.logger.write(
            "[Pipeline] ===== Step0: CD-HIT (c=0.98) preprocessing of target genome sequences ====="
        )

        patterns = ["*.faa", "*.fa", "*.fasta", "*.fna"]
        files = []
        for p in patterns:
            files.extend(glob(os.path.join(self.taget_dir, "**", p),
                              recursive=True))

        files = sorted(f for f in files if os.path.isfile(f))
        if not files:
            self.logger.write(
                "[Pipeline] Warning: No genome sequences found in taget_dir (*.faa/fa/fasta/fna)"
            )
            return

        self.logger.write(f"[Pipeline] Step0: {len(files)} genome files to be CD-HIT filtered")
        step_start = time.time()

        for idx, fa in enumerate(files, start=1):
            genome = os.path.basename(fa)
            out_fa = os.path.join(self.cdhit_target_dir, genome)
            tmp_out = out_fa + ".cdhit98.tmp"

            cmd = [
                "cd-hit",
                "-i", fa,
                "-o", tmp_out,
                "-c", "0.98",
                "-n", "5",
                "-T", "5"
            ]

            subprocess.run(
                cmd,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=True
            )

            shutil.move(tmp_out, out_fa)

            clstr = tmp_out + ".clstr"
            if os.path.exists(clstr):
                os.remove(clstr)

            self.logger.progress(
                "Pipeline",
                "CD-HIT taget_dir",
                idx,
                len(files),
                step_start
            )

        self.logger.write(
            "[Pipeline] Step0 completed: CDHIT_Target contains all filtered sequences."
        )


    # --------------------------------------------------------------
    # Step1: Reference FASTA → Working_Dir → QC + CD-HIT
    # --------------------------------------------------------------

    def prepare_ref_cdhit_input(self):
        """
        Copy reference FASTA files into Working_Dir for Script1 processing.
        """
        patterns = ["*.fasta", "*.fa", "*.faa", "*.pep"]
        files = []

        for p in patterns:
            files.extend(glob(os.path.join(self.ref_seq_dir, p)))

        if not files:
            self.logger.write(
                "[Pipeline] Warning: No reference FASTA files found (*.fasta/fa/faa/pep)"
            )
            return

        self.logger.write(
            f"[Pipeline] Copying {len(files)} reference FASTA files → Working_Dir"
        )

        for f in files:
            shutil.copy2(
                f,
                os.path.join(self.working_dir, os.path.basename(f))
            )

    def run_step1(self):
        """
        Run AutoFilterCDHIT on reference FASTA files to generate:

            Working_Dir/ref.fasta
        """
        self.logger.write("\n[Pipeline] ===== Step1: AutoFilterCDHIT =====")
        self.prepare_ref_cdhit_input()

        s1 = AutoFilterCDHIT(
            working_dir=self.working_dir,
            threads=5,
            identity=0.7,
            wordsize=5,
            sd_multiplier=2,
            output_merged="ref.fasta",
            logger=self.logger
        )
        s1.run()


    # --------------------------------------------------------------
    # Step2: AutoBlastPExtractor
    # --------------------------------------------------------------

    def run_step2(self):
        """
        Run BLASTp screening of CDHIT_Target sequences against ref.fasta.
        Produces merged per-gene FASTA files and BLAST statistics.
        """
        self.logger.write("\n[Pipeline] ===== Step2: AutoBlastPExtractor =====")

        results_dir = os.path.join(
            self.working_dir,
            "AutoBlastPExtractor_Results"
        )

        s2 = AutoBlastPExtractor(
            ref_fasta=self.ref_fasta,
            ref_db_prefix=self.ref_db_prefix,
            query_dir=self.cdhit_target_dir,
            evalue=1e-5,
            identity=30,
            coverage=70,
            threads_per_blast=BLAST_THREADS_PER_TASK,
            parallel_jobs=BLAST_MAX_PARALLEL,
            allow_multi_ref=True,
            add_genome_prefix=True,
            force_rebuild_db=False,
            results_dir=results_dir,
            logger=self.logger,
        )
        s2.run()

        # Remove ref.fasta & BLAST database files
        if os.path.exists(self.ref_fasta):
            os.remove(self.ref_fasta)
        for ext in [".phr", ".pin", ".psq"]:
            f = self.ref_db_prefix + ext
            if os.path.exists(f):
                os.remove(f)


    # --------------------------------------------------------------
    # Step3: PfamGenomeScreen
    # --------------------------------------------------------------

    def run_step3(self):
        """
        Run HMM-based domain validation (PfamGenomeScreen).
        """
        self.logger.write("\n[Pipeline] ===== Step3: PfamGenomeScreen =====")

        config_path = os.path.join(
            self.ref_hmm_dir,
            "Pfam_hmmserach-config.txt"
        )
        if not os.path.isfile(config_path):
            raise RuntimeError(
                f"[Pipeline] Pfam_hmmserach-config.txt not found: {config_path}"
            )

        s3 = PfamGenomeScreen(
            fasta_merged_dir=os.path.join(
                self.working_dir,
                "AutoBlastPExtractor_Results",
                "FastaMerged"
            ),
            ref_hmm_dir=self.ref_hmm_dir,
            working_dir=self.working_dir,
            logger=self.logger
        )
        s3.run()


    # --------------------------------------------------------------
    # Step4: AutoMSAConcatenation
    # --------------------------------------------------------------

    def run_step4(self):
        """
        Run Script4 to perform two-stage alignment, trimming,
        tree inference, pruning, and supermatrix concatenation.
        """
        self.logger.write("\n[Pipeline] ===== Step4: AutoMSAConcatenation =====")

        fasta_dir = os.path.join(
            self.working_dir,
            "Pfam_hmmserach-config_Results"
        )
        if not os.path.isdir(fasta_dir):
            self.logger.write(f"[Pipeline] Warning: directory not found: {fasta_dir}")
            return

        s4 = AutoMSAConcatenation(
            fasta_dir=fasta_dir,
            logger=self.logger
        )
        s4.run()

        # Copy concatenation outputs to final supermatrix directory
        concat_dir = os.path.join(
            fasta_dir,
            "Secondary_MAFFT_TRIMAL",
            "Concatenation"
        )

        if os.path.isdir(concat_dir):
            for f in os.listdir(concat_dir):
                src = os.path.join(concat_dir, f)
                dst = os.path.join(self.supermatrix_outdir, f)
                shutil.copy2(src, dst)

            self.logger.write("[Pipeline] Concatenation results copied to supermatrix directory.")
        else:
            self.logger.write("[Pipeline] Warning: Concatenation directory not found.")


    # --------------------------------------------------------------
    # Final cleanup after Step4
    # --------------------------------------------------------------

    def finalize_log(self, total_seconds):
        """Write final timing info and copy log file to final directory."""
        h = total_seconds // 3600
        m = (total_seconds % 3600) // 60
        s = total_seconds % 60

        self.logger.write(f"[Pipeline] Total time: {h:02d}h:{m:02d}m:{s:02d}s")

        final_log = os.path.join(
            self.supermatrix_outdir,
            "Supermatrix_Pipeline.log"
        )
        shutil.copy2(self.pipeline_log_tmp, final_log)
        self.logger.write(f"[Pipeline] Log saved to: {final_log}")

    # --------------------------------------------------------------
    # All-in-one execution
    # --------------------------------------------------------------

    def run_all(self):
        """
        Execute the entire pipeline from Step0 to Step4,
        then perform cleanup and finalize logs.
        """
        global_start = time.time()

        # ---------------- Step0 ----------------
        try:
            self.prepare_target_cdhit()
        except Exception as e:
            self.logger.write(f"[Pipeline] Step0 error: {e}")
            raise

        # ---------------- Step1 ----------------
        try:
            self.run_step1()
        except Exception as e:
            self.logger.write(f"[Pipeline] Step1 error: {e}")
            raise

        # ---------------- Step2 ----------------
        try:
            self.run_step2()
        except Exception as e:
            self.logger.write(f"[Pipeline] Step2 error: {e}")
            raise

        # ---------------- Step3 ----------------
        try:
            self.run_step3()
        except Exception as e:
            self.logger.write(f"[Pipeline] Step3 error: {e}")
            raise

        # ---------------- Step4 ----------------
        try:
            self.run_step4()
        except Exception as e:
            self.logger.write(f"[Pipeline] Step4 error: {e}")
            raise

        # ==========================================================
        # Cleanup after Step4
        # ==========================================================
        try:
            self.logger.write("[Pipeline] ===== Starting cleanup tasks =====")

            # Remove all non-directory files in Working_Dir
            for item in os.listdir(self.working_dir):
                p = os.path.join(self.working_dir, item)
                if os.path.isfile(p):
                    try:
                        os.remove(p)
                    except Exception:
                        pass

            # Remove CDHIT_Target
            if os.path.isdir(self.cdhit_target_dir):
                shutil.rmtree(self.cdhit_target_dir, ignore_errors=True)

            # Remove Blast CSV directories
            s2_res = os.path.join(
                self.working_dir,
                "AutoBlastPExtractor_Results"
            )
            rm_dirs = ["Csv", "Csv-threshold", "Csv-threshold-best"]
            for d in rm_dirs:
                p = os.path.join(s2_res, d)
                if os.path.isdir(p):
                    shutil.rmtree(p, ignore_errors=True)

            # Remove *.nostop.fasta in Pfam results
            p3 = os.path.join(
                self.working_dir,
                "Pfam_hmmserach-config_Results"
            )
            for f in glob(os.path.join(p3, "*.nostop.fasta")):
                try:
                    os.remove(f)
                except Exception:
                    pass

            # Remove *.fasta in Pfam results (but not concatenated outputs)
            for f in glob(os.path.join(p3, "*.fasta")):
                try:
                    os.remove(f)
                except Exception:
                    pass

            self.logger.write("[Pipeline] Cleanup tasks completed.")

        except Exception as e:
            self.logger.write(f"[Pipeline] Cleanup error: {e}")

        # ==========================================================
        # Finalize logging
        # ==========================================================
        total_seconds = int(time.time() - global_start)
        self.finalize_log(total_seconds)
# ==============================================================
# main() — Entry point
# ==============================================================

def main():
    """
    Main entry for running the entire supermatrix pipeline.
    """
    pipeline = SupermatrixPipeline()
    pipeline.run_all()


if __name__ == "__main__":
    main()
