#!/usr/bin/env python3
# local_coverct_inputs.py
#
# Fast, streaming DP feature extraction for BAM/CRAM and VCF.
# - BAM/CRAM coverage via `samtools depth` (compiled C; no intermediates).
# - VCF variant DP/AD summaries via pysam (htslib).
# - Multi-process across files, per-process samtools threads.
# - Concurrency can be capped by --ram-gb *and* a per-job memory heuristic.

import argparse
import concurrent.futures as cf
import math
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

try:
    import numpy as np
except Exception as e:
    print("ERROR: numpy is required. Install with `pip install numpy`.", file=sys.stderr)
    raise

# Optional: psutil for smarter RAM gating. We keep it optional.
try:
    import psutil
except Exception:
    psutil = None

# ---------- Config ----------
MAX_DP = 10000  # depths >= MAX_DP are clamped into MAX_DP bin
DP_THRESHOLDS = [10, 20, 30, 50, 100, 200, 300, 500, 1000, 1500, 2000, 3000, 5000, 10000]
SPAN_BINS = [
    (0, 0, 0),       # DP == 0
    (1, 9, 1),       # DP 1-9
    (10, 49, 10),    # DP 10-49
    (50, 99, 50),    # DP 50-99
    (100, 999, 100), # DP 100-999
]
EPS = 1e-4

# ---------- Utilities ----------
def log2_safe(x: float) -> float:
    return math.log(x + EPS, 2)

def as_gb(val: float) -> float:
    return float(val) / (1024.0 ** 3)

def basename_no_ext(p: str) -> str:
    b = os.path.basename(p)
    return b[:-len(os.path.splitext(b)[1])] if "." in b else b

# ---------- BAM/CRAM: coverage summarization ----------
@dataclass
class BamSummary:
    sample_name: str
    file_path: str
    total_positions: int
    dp_mean_DP: float
    dp_sd_DP: float
    dp_chaos: Optional[float]
    fractions: Dict[str, float]      # fraction_below_* features
    cs: Dict[str, float]             # cs* features
    # "span" analogs
    span_span_count: int
    span_weighted_span_count: float
    span_log_weighted_span_count: float
    span_mean_DP: float
    span_sd_DP: float
    span_fraction_below_1: float
    span_fraction_below_10: float
    span_fraction_below_50: float
    span_fraction_below_100: float
    span_fraction_below_1000: float
    span_cs1: float
    span_cs10: float
    span_cs50: float
    span_cs100: float
    span_cs1000: float
    span_chaos: Optional[float]

def run_samtools_depth(
    bam: str, samtools: str, threads: int, region: Optional[str], include_all: bool
) -> Iterable[str]:
    """
    Stream 'samtools depth' lines: CHR POS DP
    - Use '-aa' to include zero-coverage positions across the entire reference contigs in the file header.
    - '-Q 0 -d 0' to avoid quality/depth caps (we clamp in Python).
    - '-@ threads' for parallel bgzip decompression where possible.
    """
    if shutil.which(samtools) is None:
        raise RuntimeError(f"`{samtools}` not found in PATH.")

    cmd = [samtools, "depth"]
    if include_all:
        cmd.append("-aa")  # absolutely all positions, including zeros
    else:
        cmd.append("-a")   # include zeros within covered regions

    cmd += ["-Q", "0", "-d", "0", "-@", str(max(1, threads))]

    if region:
        cmd += [bam, region]
    else:
        cmd += [bam]

    # Use text mode stream
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1)
    # stream stdout lines
    if p.stdout is None:
        raise RuntimeError("Failed to obtain stdout from samtools.")

    for line in p.stdout:
        yield line

    # ensure process finished
    p.stdout.close()
    ret = p.wait()
    if ret != 0:
        err = p.stderr.read() if p.stderr else ""
        raise RuntimeError(f"samtools depth failed [{ret}] on {bam}.\n{err}")

def summarize_depth_stream(lines: Iterable[str], sample_name: str, file_path: str) -> BamSummary:
    # histogram over 0..MAX_DP (clamped)
    hist = np.zeros(MAX_DP + 1, dtype=np.int64)
    total = 0
    sum_dp = 0
    sum_dp2 = 0

    span_count = 0
    prev_dp = None

    for ln in lines:
        if not ln or ln[0] == "#":
            continue
        # depth format: chrom \t pos \t depth
        parts = ln.split()
        if len(parts) < 3:
            continue
        try:
            dp = int(parts[2])
        except ValueError:
            continue
        if dp < 0:
            dp = 0
        if dp > MAX_DP:
            dp = MAX_DP

        # update run-length span count (count changes in DP)
        if prev_dp is None or dp != prev_dp:
            span_count += 1
            prev_dp = dp

        hist[dp] += 1
        total += 1
        sum_dp += dp
        sum_dp2 += dp * dp

    if total == 0:
        # Empty result; return zeros
        zero_frac = {f"fraction_below_{t}": 0.0 for t in DP_THRESHOLDS}
        zero_cs = {f"cs{t}": log2_safe(0.0) for t in DP_THRESHOLDS}
        return BamSummary(
            sample_name=sample_name,
            file_path=file_path,
            total_positions=0,
            dp_mean_DP=0.0,
            dp_sd_DP=0.0,
            dp_chaos=None,
            fractions=zero_frac,
            cs=zero_cs,
            span_span_count=0,
            span_weighted_span_count=0.0,
            span_log_weighted_span_count=0.0,
            span_mean_DP=0.0,
            span_sd_DP=0.0,
            span_fraction_below_1=0.0,
            span_fraction_below_10=0.0,
            span_fraction_below_50=0.0,
            span_fraction_below_100=0.0,
            span_fraction_below_1000=0.0,
            span_cs1=log2_safe(0.0),
            span_cs10=log2_safe(0.0),
            span_cs50=log2_safe(0.0),
            span_cs100=log2_safe(0.0),
            span_cs1000=log2_safe(0.0),
            span_chaos=None,
        )

    mean = float(sum_dp) / float(total)
    var = max(0.0, float(sum_dp2) / float(total) - mean * mean)
    sd = math.sqrt(var)
    chaos = (-mean / sd) if sd > 0 else None

    cdf = np.cumsum(hist, dtype=np.int64)  # length MAX_DP+1; cdf[k] = count of dp <= k

    # Fractions "< t" (to match your SQL COUNTIF(DP < t))
    fractions = {}
    cs_vals = {}
    for t in DP_THRESHOLDS:
        idx = min(MAX_DP, t - 1)  # "< t" == "<= t-1"
        below = int(cdf[idx]) if idx >= 0 else 0
        frac = below / float(total)
        fractions[f"fraction_below_{t}"] = frac
        cs_vals[f"cs{t}"] = log2_safe(frac)

    # Span-like aggregates on a few fixed boundaries (<= x)
    # <=0, <=9, <=49, <=99, <=999
    def frac_le(k: int) -> float:
        k = min(k, MAX_DP)
        return float(int(cdf[k])) / float(total)

    span_fraction_below_1 = frac_le(0)
    span_fraction_below_10 = frac_le(9)
    span_fraction_below_50 = frac_le(49)
    span_fraction_below_100 = frac_le(99)
    span_fraction_below_1000 = frac_le(999)

    span_cs1 = log2_safe(span_fraction_below_1)
    span_cs10 = log2_safe(span_fraction_below_10)
    span_cs50 = log2_safe(span_fraction_below_50)
    span_cs100 = log2_safe(span_fraction_below_100)
    span_cs1000 = log2_safe(span_fraction_below_1000)

    # Weighted span counts per your l_bnd definition
    def range_count(lo: int, hi: int) -> int:
        lo = max(0, lo)
        hi = min(MAX_DP, hi)
        if lo == 0:
            return int(cdf[hi])
        else:
            return int(cdf[hi] - cdf[lo - 1])

    weighted_span = 0.0
    log_weighted_span = 0.0
    for lo, hi, l_bnd in SPAN_BINS:
        cnt = range_count(lo, hi)
        weighted_span += l_bnd * cnt
        if l_bnd > 0:
            log_weighted_span += (1.0 / math.log10(l_bnd + 1.0)) * cnt
        # l_bnd == 0 contributes nothing to log-weighted span

    return BamSummary(
        sample_name=sample_name,
        file_path=file_path,
        total_positions=int(total),
        dp_mean_DP=mean,
        dp_sd_DP=sd,
        dp_chaos=chaos,
        fractions=fractions,
        cs=cs_vals,
        span_span_count=int(span_count),
        span_weighted_span_count=float(weighted_span),
        span_log_weighted_span_count=float(log_weighted_span),
        span_mean_DP=mean,
        span_sd_DP=sd,
        span_fraction_below_1=span_fraction_below_1,
        span_fraction_below_10=span_fraction_below_10,
        span_fraction_below_50=span_fraction_below_50,
        span_fraction_below_100=span_fraction_below_100,
        span_fraction_below_1000=span_fraction_below_1000,
        span_cs1=span_cs1,
        span_cs10=span_cs10,
        span_cs50=span_cs50,
        span_cs100=span_cs100,
        span_cs1000=span_cs1000,
        span_chaos=chaos,
    )

def process_bam_file(
    bam_path: str, samtools: str, threads: int, include_all: bool, region: Optional[str], sample_name_from_rg: bool
) -> Dict[str, str]:
    sample_name = basename_no_ext(bam_path)

    # Optional: take SM from RG header as sample name (first RG:SM encountered)
    if sample_name_from_rg:
        try:
            import pysam
            with pysam.AlignmentFile(bam_path, "rb") as af:
                for rg in af.header.get("RG", []):
                    sm = rg.get("SM")
                    if sm:
                        sample_name = sm
                        break
        except Exception:
            pass  # fall back to filename

    summary = summarize_depth_stream(
        run_samtools_depth(bam_path, samtools, threads, region, include_all),
        sample_name=sample_name,
        file_path=os.path.abspath(bam_path),
    )

    # Flatten to row dict compatible with CSV
    row: Dict[str, str] = {
        "sample_name": summary.sample_name,
        "file_path": summary.file_path,
        "total_positions": str(summary.total_positions),
        "dp_mean_DP": f"{summary.dp_mean_DP:.6g}",
        "dp_sd_DP": f"{summary.dp_sd_DP:.6g}",
        "dp_chaos": "" if summary.dp_chaos is None else f"{summary.dp_chaos:.6g}",
        "span_span_count": str(summary.span_span_count),
        "span_weighted_span_count": f"{summary.span_weighted_span_count:.6g}",
        "span_log_weighted_span_count": f"{summary.span_log_weighted_span_count:.6g}",
        "span_mean_DP": f"{summary.span_mean_DP:.6g}",
        "span_sd_DP": f"{summary.span_sd_DP:.6g}",
        "span_fraction_below_1": f"{summary.span_fraction_below_1:.6g}",
        "span_fraction_below_10": f"{summary.span_fraction_below_10:.6g}",
        "span_fraction_below_50": f"{summary.span_fraction_below_50:.6g}",
        "span_fraction_below_100": f"{summary.span_fraction_below_100:.6g}",
        "span_fraction_below_1000": f"{summary.span_fraction_below_1000:.6g}",
        "span_cs1": f"{summary.span_cs1:.6g}",
        "span_cs10": f"{summary.span_cs10:.6g}",
        "span_cs50": f"{summary.span_cs50:.6g}",
        "span_cs100": f"{summary.span_cs100:.6g}",
        "span_cs1000": f"{summary.span_cs1000:.6g}",
        "span_chaos": "" if summary.span_chaos is None else f"{summary.span_chaos:.6g}",
    }
    for t in DP_THRESHOLDS:
        row[f"fraction_below_{t}"] = f"{summary.fractions[f'fraction_below_{t}']:.6g}"
        row[f"cs{t}"] = f"{summary.cs[f'cs{t}']:.6g}"

    return row

# ---------- VCF: variant-site DP/AD summarization ----------
@dataclass
class VcfSummary:
    sample_name: str
    file_path: str
    total_variants: int
    major_variants: int
    minor_variants: int
    proportion_minor_variants: float
    mean_DP: float
    sd_DP: float
    chaos: Optional[float]
    fractions: Dict[str, float]
    cs: Dict[str, float]

def process_vcf_file(vcf_path: str) -> List[Dict[str, str]]:
    try:
        import pysam
    except Exception:
        print("ERROR: pysam is required for VCF mode. Install with `pip install pysam`.", file=sys.stderr)
        raise

    vf = pysam.VariantFile(vcf_path)
    sample_names = list(vf.header.samples)
    out_rows: List[Dict[str, str]] = []

    # If no samples, treat as single "unsampled" stream using INFO.DP only (limited)
    if not sample_names:
        sample_names = ["_no_sample_"]

    # Prepare per-sample accumulators
    acc = {}
    for sm in sample_names:
        acc[sm] = {
            "dp_values": [],  # variant-site DP
            "total": 0,
            "major": 0,
            "minor": 0,
        }

    for rec in vf:
        for sm in sample_names:
            # Extract counts
            if sm == "_no_sample_":
                dp = rec.info.get("DP")
                alt_count = None
            else:
                sdata = rec.samples[sm]
                dp = sdata.get("DP", rec.info.get("DP", None))
                ad = sdata.get("AD", None)  # array: [ref, alt1, alt2...]
                ao = sdata.get("AO", None)
                ro = sdata.get("RO", None)
                alt_count = None

                if ad is not None and len(ad) >= 2:
                    # use the alt allele with max count
                    try:
                        alt_count = max(int(x) for x in ad[1:] if x is not None)
                    except ValueError:
                        alt_count = None
                    if dp is None:
                        # DP might be missing; reconstruct if ref+alts available
                        try:
                            dp = int(ad[0]) + sum(int(x) for x in ad[1:] if x is not None)
                        except Exception:
                            pass
                elif ao is not None and ro is not None:
                    # FreeBayes-style AO (alts), RO (ref)
                    try:
                        if isinstance(ao, (list, tuple)):
                            alt_count = max(int(x) for x in ao if x is not None)
                            alt_sum = sum(int(x) for x in ao if x is not None)
                        else:
                            alt_count = int(ao)
                            alt_sum = int(ao)
                        ref_count = int(ro)
                        if dp is None:
                            dp = ref_count + alt_sum
                    except Exception:
                        alt_count = None

            # record
            if dp is None:
                continue
            try:
                dp = int(dp)
            except Exception:
                continue
            if dp <= 0:
                # no meaningful alt fraction when DP=0
                acc[sm]["dp_values"].append(0)
                acc[sm]["total"] += 1
                continue

            frac = None
            if alt_count is not None:
                try:
                    frac = float(alt_count) / float(dp) if dp > 0 else None
                except Exception:
                    frac = None

            acc[sm]["dp_values"].append(dp)
            acc[sm]["total"] += 1
            if frac is not None:
                if frac > 0.5:
                    acc[sm]["major"] += 1
                elif frac < 0.5:
                    acc[sm]["minor"] += 1
                # exactly 0.5 ignored (neither major nor minor)

    for sm in sample_names:
        dps = np.array(acc[sm]["dp_values"], dtype=np.int64)
        total = int(acc[sm]["total"])
        major = int(acc[sm]["major"])
        minor = int(acc[sm]["minor"])

        if total > 0:
            mean = float(dps.mean())
            sd = float(dps.std(ddof=0))
            chaos = (-mean / sd) if sd > 0 else None
        else:
            mean = 0.0
            sd = 0.0
            chaos = None

        # fractions "< t" over variant DP values
        fractions = {}
        cs_vals = {}
        if total > 0:
            # histogram/cdf up to MAX_DP clamp like BAM
            hist = np.bincount(np.clip(dps, 0, MAX_DP), minlength=MAX_DP+1)
            cdf = np.cumsum(hist)
            for t in DP_THRESHOLDS:
                idx = min(MAX_DP, t - 1)
                below = int(cdf[idx]) if idx >= 0 else 0
                frac = below / float(total)
                fractions[f"fraction_below_{t}"] = frac
                cs_vals[f"cs{t}"] = log2_safe(frac)
        else:
            for t in DP_THRESHOLDS:
                fractions[f"fraction_below_{t}"] = 0.0
                cs_vals[f"cs{t}"] = log2_safe(0.0)

        prop_minor = (minor / total) if total > 0 else 0.0

        row = {
            "sample_name": f"{basename_no_ext(vcf_path)}:{sm}",
            "file_path": os.path.abspath(vcf_path),
            "total_variants": str(total),
            "major_variants": str(major),
            "minor_variants": str(minor),
            "proportion_minor_variants": f"{prop_minor:.6g}",
            "mean_DP": f"{mean:.6g}",
            "sd_DP": f"{sd:.6g}",
            "chaos": "" if chaos is None else f"{chaos:.6g}",
        }
        for t in DP_THRESHOLDS:
            row[f"fraction_below_{t}"] = f"{fractions[f'fraction_below_{t}']:.6g}"
            row[f"cs{t}"] = f"{cs_vals[f'cs{t}']:.6g}"

        out_rows.append(row)

    return out_rows

# ---------- Orchestration ----------
def compute_worker_count(max_procs: int, ram_gb: Optional[float], ram_per_job_mb: int) -> int:
    if max_procs is None or max_procs <= 0:
        max_procs = os.cpu_count() or 1
    if ram_gb is None:
        return max(1, int(max_procs))
    jobs_by_ram = max(1, int((ram_gb * 1024.0) // float(ram_per_job_mb)))
    # Optional: cap by current free RAM if psutil is available
    if psutil is not None:
        avail_gb = as_gb(psutil.virtual_memory().available)
        jobs_by_avail = max(1, int((avail_gb * 1024.0) // float(ram_per_job_mb)))
        jobs_by_ram = min(jobs_by_ram, jobs_by_avail)
    return max(1, min(max_procs, jobs_by_ram))

def detect_mode_from_ext(path: str) -> str:
    ext = os.path.splitext(path)[1].lower()
    if ext in (".bam", ".cram", ".sam"):
        return "bam"
    if ext in (".vcf", ".vcf.gz", ".bcf"):
        return "vcf"
    return "auto"

def main():
    ap = argparse.ArgumentParser(description="Fast DP features from BAM/CRAM (coverage) and VCF (variant DP).")
    ap.add_argument("--inputs", nargs="+", required=True, help="List of BAM/CRAM/VCF files (mix allowed).")
    ap.add_argument("--mode", choices=["auto", "bam", "vcf"], default="auto",
                    help="Force mode or auto-detect by extension.")
    ap.add_argument("--output-csv", required=True, help="Path to output CSV.")
    ap.add_argument("--max-procs", type=int, default=0, help="Max worker processes (default: all cores).")
    ap.add_argument("--ram-gb", type=float, default=None, help="Cap total RAM for workers (GB).")
    ap.add_argument("--ram-per-job-mb", type=int, default=256, help="Heuristic RAM per job (MB).")
    ap.add_argument("--samtools", default="samtools", help="Path to samtools.")
    ap.add_argument("--samtools-threads", type=int, default=2, help="Threads passed to samtools -@.")
    ap.add_argument("--include-all-positions", action="store_true",
                    help="Use 'samtools depth -aa' to include all reference positions (recommended for viral genomes).")
    ap.add_argument("--region", default=None, help="Optional region (e.g., 'MN908947.3' or 'chr1:1-100000').")
    ap.add_argument("--sample-name-from-rg", action="store_true",
                    help="For BAM/CRAM: use RG:SM as sample_name if present.")

    args = ap.parse_args()

    # Prepare CSV
    out_path = args.output_csv
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    out = open(out_path, "w", buffering=1)
    wrote_header = False

    def write_row(row: Dict[str, str]):
        nonlocal wrote_header
        if not wrote_header:
            out.write(",".join(row.keys()) + "\n")
            wrote_header = True
        out.write(",".join(str(row[k]) for k in row.keys()) + "\n")

    # Determine concurrency
    workers = compute_worker_count(args.max_procs, args.ram_gb, args.ram_per_job_mb)

    # Work queue
    tasks: List[Tuple[str, str]] = []
    for p in args.inputs:
        if not os.path.exists(p):
            print(f"WARNING: input not found: {p}", file=sys.stderr)
            continue
        m = args.mode if args.mode != "auto" else detect_mode_from_ext(p)
        if m == "auto":
            # try BAM first
            m = "bam"
        tasks.append((p, m))

    # Process
    with cf.ProcessPoolExecutor(max_workers=workers) as ex:
        futs = []
        for path, mode in tasks:
            if mode == "bam":
                futs.append(
                    ex.submit(
                        process_bam_file, path, args.samtools, args.samtools_threads,
                        args.include_all_positions, args.region, args.sample_name_from_rg
                    )
                )
            elif mode == "vcf":
                futs.append(ex.submit(process_vcf_file, path))
            else:
                # Default to BAM behavior
                futs.append(
                    ex.submit(
                        process_bam_file, path, args.samtools, args.samtools_threads,
                        args.include_all_positions, args.region, args.sample_name_from_rg
                    )
                )

        # Gather as they complete
        for f in cf.as_completed(futs):
            try:
                res = f.result()
            except Exception as e:
                print(f"ERROR: {e}", file=sys.stderr)
                continue

            if isinstance(res, dict):
                write_row(res)
            elif isinstance(res, list):
                for row in res:
                    write_row(row)
            else:
                print("WARNING: Unknown result type.", file=sys.stderr)

    out.close()
    print(f"Done. Wrote {out_path}")

if __name__ == "__main__":
    main()
