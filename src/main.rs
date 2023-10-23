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

    let (_, seq) = read_fasta(&args.fasta);
    eprintln!("Sequence length: {}", seq.len());
    let seq = seq
        .into_iter()
        .filter(|&b| b.to_ascii_uppercase() != b'N')
        .collect::<Vec<_>>();
    eprintln!("Sequence length (no Ns): {}", seq.len());

    let start_saca = Instant::now();
    let suffix_array = SuffixArray::<5>::new_packed::<248>(&seq, 9);
    eprintln!("Suffix array length: {}", suffix_array.idxs().len());
    let elapsed_saca = start_saca.elapsed().as_secs_f64();
    eprintln!("Suffix array construction run time (s): {elapsed_saca}");

    let elapsed = start.elapsed().as_secs_f64();
    eprintln!("Total run time (s): {elapsed}");
    let mem = (max_mem_usage_mb() - start_mem).max(0.0);
    eprintln!("Peak memory usage (MB): {mem}");
}

fn read_fasta(path: &Path) -> (String, Vec<u8>) {
    let mut r = parse_fastx_file(path).unwrap();

    if let Some(record) = r.next() {
        let record = record.unwrap();
        let id = std::str::from_utf8(record.id()).unwrap().to_owned();
        let seq = record.seq().into_owned();
        (id, seq)
    } else {
        panic!("Empty file!")
    }
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
    #[arg(short, long, default_value_t = 16)]
    threads: usize,
}