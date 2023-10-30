use needletail::*;

use clap::{self, Parser};

use std::path::{Path, PathBuf};
use std::time::Instant;

use simple_saca::*;

fn main() {
    let start = Instant::now();
    let start_mem = max_mem_usage_mb();
    let args = Args::parse();
    eprintln!("{args:?}");

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    let seq = read_fasta(&args.fasta);
    eprintln!("Sequence length: {}", seq.len());
    let seq = seq
        .into_iter()
        .filter(|&b| b.to_ascii_uppercase() != b'N')
        .collect::<Vec<_>>();
    eprintln!("Sequence length (no Ns): {}", seq.len());

    let start_saca = Instant::now();

    if args.divsufsort {
        use libdivsufsort_rs::divsufsort64;
        let suffix_array = divsufsort64(&seq).unwrap();
        eprintln!("Suffix array length: {}", suffix_array.len());
    } else {
        let suffix_array = match args.ctx {
            125 => SuffixArray::<5>::new_packed::<125>(&seq, args.k, args.bucket_threads),
            250 => SuffixArray::<5>::new_packed::<250>(&seq, args.k, args.bucket_threads),
            500 => SuffixArray::<5>::new_packed::<500>(&seq, args.k, args.bucket_threads),
            1000 => SuffixArray::<5>::new_packed::<1000>(&seq, args.k, args.bucket_threads),
            _ => panic!("Context length of {} is not supported!", args.ctx),
        };
        eprintln!("Suffix array length: {}", suffix_array.idxs().len());
    }

    let elapsed_saca = start_saca.elapsed().as_secs_f64();
    eprintln!("Suffix array construction run time (s): {elapsed_saca}");

    let elapsed = start.elapsed().as_secs_f64();
    eprintln!("Total run time (s): {elapsed}");
    let mem = (max_mem_usage_mb() - start_mem).max(0.0);
    eprintln!("Peak memory usage (MB): {mem}");
}

fn read_fasta(path: &Path) -> Vec<u8> {
    let mut r = parse_fastx_file(path).unwrap();
    let mut seq = Vec::new();

    while let Some(record) = r.next() {
        let record = record.unwrap();
        seq.extend_from_slice(&record.seq());
    }

    seq
}

fn max_mem_usage_mb() -> f64 {
    let rusage = unsafe {
        let mut rusage = std::mem::MaybeUninit::uninit();
        libc::getrusage(libc::RUSAGE_SELF, rusage.as_mut_ptr());
        rusage.assume_init()
    };
    let maxrss = rusage.ru_maxrss as f64;
    if cfg!(target_os = "macos") {
        maxrss / 1024.0 / 1024.0
    } else {
        maxrss / 1024.0
    }
}

#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    /// Input fasta file.
    fasta: PathBuf,
    /// Number of threads to use.
    #[arg(short, long, default_value_t = 32)]
    threads: usize,
    /// Number of threads to use for bucketing.
    #[arg(short, long, default_value_t = 32)]
    bucket_threads: usize,
    /// Number of base pairs to use for bucketing.
    #[arg(short, long, default_value_t = 10)]
    k: usize,
    /// Number of base pairs of context to use for sorting.
    ///
    /// Supported multiples of 125: 125, 250, 500, 1000
    #[arg(short, long, default_value_t = 250)]
    ctx: usize,
    /// Run 64-bit libdivsufsort instead.
    #[arg(long)]
    divsufsort: bool,
}
