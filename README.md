# simple-saca
Hardware go brrr bounded context suffix array construction algorithm.

## Benchmark
Results for T2T-CHM13+Y (complete human genome sequence, 3.1 billion base pairs) on AWS c7i.8xlarge (32 vCPUs, 64GB memory):
* **32 threads, 248bp context: 61s, 18.7GB**
* 16 threads, 248bp context: 83s, 18.6GB
* 8 threads, 248bp context: 159s, 18.6GB
* **32 threads, 992bp context: 64s, 18.7GB**
* libdivsufsort (64-bit): 418s, 26.7GB

The peak memory usage is barely larger than the space needed to simply store the genome and the suffix array.

Note that 10-mer buckets are used, and the largest bucket has around 3 million suffixes. This is the largest bucket
that that will be sorted with comparison-based sorting.

## Why?
I was curious how fast bounded context suffix array sorting can be optimized by using a naive sorting approach,
but making full use of CPU parallelism.

Since this algorithm uses bucketing and comparison-based sorting, it should be easily adaptable to non-standard
suffix array use cases (eg., spaced seeds) by using different functions for extracting a "kmer" for bucketing and
comparing two suffixes.

## Algorithm overview
Note that this algorithm does not output fully sorted suffixes.
Each suffix is sorted up to a fixed number of base pairs (the "context").
This is useful for many genomics seeding methods.

1. Convert genome to 2-bit representation, in reverse to make later lexicographic comparisons easier.
2. Get the first `k` base pairs (kmer) for each suffix and count them separately in parallel.
The idea is to bucket suffixes by its kmer.
3. Prefix sum the counts to get bucket boundaries for each kmer.
4. Place suffixes into kmer buckets in parallel.
5. Sort each bucket separately in parallel by using the comparison-based sort from Rust's standard library.
Suffixes are compared lexicographically up to the bounded context length by using AVX2 SIMD.

Suffix indexes and kmer counts are stored using 40-bit integers to save space.

## Run
1. Clone this repo and install Rust.
2. Make sure you are running this on x86 CPUs supporting AVX2.
3. `cargo run --release -- genome.fasta.gz`

Use the `third_party` feature flag to benchmark other suffix array construction libraries:
```
cargo run --release --features third_party -- genome.fasta.gz
```

The binary only benchmarks the suffix array construction algorithm, but you can use this crate as a library.

Use `--help` to see all options. You can adjust the number of threads
and bounded context length.
