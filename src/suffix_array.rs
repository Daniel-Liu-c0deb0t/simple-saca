use rayon::prelude::*;

use std::cmp::Ordering;
use std::time::Instant;

#[cfg(target_arch = "x86")]
use std::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

use crate::compact_vec::*;

pub struct SuffixArray<const BYTES: usize> {
    idxs: CompactVec<BYTES>,
    k: usize,
    ctx: usize,
}

impl<const BYTES: usize> SuffixArray<BYTES> {
    pub fn new_packed<const CTX: usize>(bytes: &[u8], k: usize) -> Self {
        let idxs = unsafe { Self::sort_packed::<CTX>(bytes, k) };

        Self { idxs, k, ctx: CTX }
    }

    #[target_feature(enable = "avx2")]
    unsafe fn sort_packed<const CTX: usize>(bytes: &[u8], k: usize) -> CompactVec<BYTES> {
        let start = Instant::now();
        let k = k * 2;
        let packed = RevPacked::new(bytes);
        let bytes_no_ctx = &bytes[..bytes.len() - CTX];
        let mask = (1u32 << k) - 1;

        let mut counts = CompactVec::<BYTES>::new(1 << k);

        for i in 0..bytes_no_ctx.len() {
            let kmer = packed.load_24(i) & mask;
            let count = (*counts.as_ptr().add(kmer as usize)).get_usize();
            (*counts.as_mut_ptr().add(kmer as usize)).set_usize(count + 1);
        }

        let mut sum = 0;
        let mut max = 0;

        for c in counts.iter_mut() {
            let curr = c.get_usize();
            c.set_usize(sum);
            sum += curr;
            max = max.max(curr);
        }

        let mut sorted = CompactVec::<BYTES>::new(bytes_no_ctx.len());

        for i in 0..bytes_no_ctx.len() {
            let kmer = packed.load_24(i) & mask;
            let idx = (*counts.as_ptr().add(kmer as usize)).get_usize();

            (*sorted.as_mut_ptr().add(idx)).set_usize(i);
            (*counts.as_mut_ptr().add(kmer as usize)).set_usize(idx + 1);
        }

        let elapsed1 = start.elapsed().as_secs_f64();

        let start = Instant::now();
        let sorted_ptr = MutPtr(sorted.as_mut_ptr());

        (0..(1 << k)).into_par_iter().for_each(|i| {
            let start = if i == 0 {
                0
            } else {
                (*counts.as_ptr().add(i - 1)).get_usize()
            };
            let end = (*counts.as_ptr().add(i)).get_usize();
            let ptr = sorted_ptr;
            let slice = unsafe { std::slice::from_raw_parts_mut(ptr.0.add(start), end - start) };

            slice.sort_by(|a_idx, b_idx| unsafe {
                simd_cmp_packed::<CTX>(&packed, a_idx.get_usize(), b_idx.get_usize())
            });
        });

        let elapsed2 = start.elapsed().as_secs_f64();
        eprintln!("Bucket sort run time (s): {elapsed1}");
        eprintln!("Parallel sort run time (s): {elapsed2}");
        eprintln!("Largest bucket / total: {max} / {sum}");

        sorted
    }

    pub fn new_bytes<const CTX: usize>(bytes: &[u8]) -> Self {
        let idxs = unsafe { Self::sort_bytes::<CTX>(bytes) };

        Self {
            idxs,
            k: 0,
            ctx: CTX,
        }
    }

    #[target_feature(enable = "avx2")]
    unsafe fn sort_bytes<const CTX: usize>(bytes: &[u8]) -> CompactVec<BYTES> {
        let bytes_no_ctx = &bytes[..bytes.len() - CTX];

        let mut sorted = CompactVec::<BYTES>::new(bytes_no_ctx.len());

        for i in 0..bytes_no_ctx.len() {
            sorted[i].set_usize(i);
        }

        sorted.par_sort_by(|a_idx, b_idx| unsafe {
            simd_cmp_bytes::<CTX>(bytes, a_idx.get_usize(), b_idx.get_usize())
        });

        sorted
    }

    pub fn new<const CTX: usize>(seeds: &[u16], k: usize) -> Self {
        assert!(k <= 16);

        let idxs = unsafe { Self::sort::<CTX>(seeds, k) };

        Self { idxs, k, ctx: CTX }
    }

    #[target_feature(enable = "avx2")]
    unsafe fn sort<const CTX: usize>(seeds: &[u16], k: usize) -> CompactVec<BYTES> {
        let seeds_no_ctx = &seeds[..seeds.len() - CTX];

        let mut counts = CompactVec::<BYTES>::new(1 << k);

        for &s in seeds_no_ctx {
            let count = counts[s as usize].get_usize();
            counts[s as usize].set_usize(count + 1);
        }

        let mut seed_to_idx = CompactVec::<BYTES>::new((1 << k) + 1);
        let mut idx = 0;
        let mut sum = 0;

        seed_to_idx[idx].set_usize(sum);
        idx += 1;
        for c in counts.iter() {
            sum += c.get_usize();
            seed_to_idx[idx].set_usize(sum);
            idx += 1;
        }

        let mut sorted = CompactVec::<BYTES>::new(seeds_no_ctx.len());

        for (i, &s) in seeds_no_ctx.iter().enumerate() {
            let end = seed_to_idx[(s as usize) + 1].get_usize();
            let count = counts[s as usize].get_usize();

            sorted[end - count].set_usize(i);
            counts[s as usize].set_usize(count - 1);
        }

        drop(counts);

        let sorted_ptr = MutPtr(sorted.as_mut_ptr());

        (0..(1 << k)).into_par_iter().for_each(|i| {
            let start = seed_to_idx[i].get_usize();
            let end = seed_to_idx[i + 1].get_usize();
            let ptr = sorted_ptr;
            let slice = unsafe { std::slice::from_raw_parts_mut(ptr.0.add(start), end - start) };
            slice.sort_by(|a_idx, b_idx| unsafe {
                simd_cmp::<CTX>(seeds, a_idx.get_usize() + 1, b_idx.get_usize() + 1)
            });
        });

        sorted
    }

    pub fn idxs(&self) -> &[Int<BYTES>] {
        &self.idxs
    }

    pub fn k(&self) -> usize {
        self.k
    }

    pub fn ctx(&self) -> usize {
        self.ctx
    }
}

struct RevPacked {
    data: Vec<u8>,
    len: usize,
}

static LUT: [u8; 128] = {
    let mut l = [0u8; 128];
    l[b'A' as usize] = 0b00;
    l[b'C' as usize] = 0b01;
    l[b'G' as usize] = 0b10;
    l[b'T' as usize] = 0b11;
    l[b'a' as usize] = 0b00;
    l[b'c' as usize] = 0b01;
    l[b'g' as usize] = 0b10;
    l[b't' as usize] = 0b11;
    l
};

impl RevPacked {
    pub fn new(bytes: &[u8]) -> Self {
        let padded_len = bytes.len() + 4;
        let len = (padded_len + 3) / 4;
        let mut data = vec![0u8; len];

        for (i, &b) in bytes.iter().enumerate() {
            let i = padded_len - i - 1;
            unsafe {
                *data.as_mut_ptr().add(i / 4) |= *LUT.as_ptr().add(b as usize) << ((i % 4) * 2);
            }
        }

        Self {
            data,
            len: padded_len,
        }
    }

    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn load_248(&self, idx: usize) -> __m256i {
        let idx = self.len - idx - 128;
        let i = idx / 4;
        let j = idx % 4;
        let val = _mm256_loadu_si256(self.data.as_ptr().add(i) as _);

        // shift left by bits
        let left_shift = _mm256_set1_epi64x((j * 2) as _);
        let hi = _mm256_sllv_epi64(val, left_shift);
        let right_shift = _mm256_set1_epi64x(((32 - j) * 2) as _);
        let lo = _mm256_srlv_epi64(_mm256_permute4x64_epi64(val, 0b10_01_00_11), right_shift);

        let mask = _mm256_set_epi8(
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, 0,
        );

        _mm256_and_si256(_mm256_or_si256(hi, lo), mask)
    }

    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn load_24(&self, idx: usize) -> u32 {
        let idx = self.len - idx - 16;
        let i = idx / 4;
        let j = idx % 4;
        let val = std::ptr::read_unaligned(self.data.as_ptr().add(i) as *const u32);
        (val >> ((4 - j) * 2)) & 0x00FFFFFFu32
    }
}

#[inline]
#[target_feature(enable = "avx2")]
unsafe fn simd_cmp_packed<const CTX: usize>(
    packed: &RevPacked,
    a_idx: usize,
    b_idx: usize,
) -> Ordering {
    const L: usize = 128 - 4;
    let mut a_i = a_idx;
    let mut b_i = b_idx;

    for _ in 0..(CTX / L) {
        let a = packed.load_248(a_i);
        let b = packed.load_248(b_i);

        let eq = _mm256_cmpeq_epi8(a, b);
        let neq_mask = !(_mm256_movemask_epi8(eq) as u32);

        if neq_mask != 0 {
            let msb_mask = 1u32 << (31 - neq_mask.leading_zeros());
            let gt = _mm256_max_epu8(a, b);
            let gt = _mm256_cmpeq_epi8(gt, a);
            let gt_mask = _mm256_movemask_epi8(gt) as u32;

            if (msb_mask & gt_mask) > 0 {
                return Ordering::Greater;
            } else {
                return Ordering::Less;
            }
        }

        a_i += L;
        b_i += L;
    }

    a_i.cmp(&b_i)
}

#[inline]
#[target_feature(enable = "avx2")]
unsafe fn simd_cmp_bytes<const CTX: usize>(bytes: &[u8], a_idx: usize, b_idx: usize) -> Ordering {
    const L: usize = 32;
    let ptr = bytes.as_ptr();
    let mut a_i = a_idx;
    let mut b_i = b_idx;

    for _ in 0..(CTX / L) {
        let a = _mm256_loadu_si256(ptr.add(a_i) as _);
        let b = _mm256_loadu_si256(ptr.add(b_i) as _);

        let eq = _mm256_cmpeq_epi8(a, b);
        let neq_mask = !(_mm256_movemask_epi8(eq) as u32);

        if neq_mask != 0 {
            // TODO: try branchless
            let lsb_mask = neq_mask & neq_mask.wrapping_neg();
            let gt = _mm256_max_epu8(a, b);
            let gt = _mm256_cmpeq_epi8(gt, a);
            let gt_mask = _mm256_movemask_epi8(gt) as u32;

            if (lsb_mask & gt_mask) > 0 {
                return Ordering::Greater;
            } else {
                return Ordering::Less;
            }
        }

        a_i += L;
        b_i += L;
    }

    Ordering::Equal
}

#[inline]
#[target_feature(enable = "avx2")]
unsafe fn simd_cmp<const CTX: usize>(seeds: &[u16], a_idx: usize, b_idx: usize) -> Ordering {
    const L: usize = 16;
    let ptr = seeds.as_ptr();
    let mut a_i = a_idx;
    let mut b_i = b_idx;

    for _ in 0..(CTX / L) {
        let a = _mm256_loadu_si256(ptr.add(a_i) as _);
        let b = _mm256_loadu_si256(ptr.add(b_i) as _);

        let eq = _mm256_cmpeq_epi16(a, b);
        let neq_mask = !(_mm256_movemask_epi8(eq) as u32);

        if neq_mask != 0 {
            let lsb_mask = neq_mask & neq_mask.wrapping_neg();
            let gt = _mm256_max_epu16(a, b);
            let gt = _mm256_cmpeq_epi16(gt, a);
            let gt_mask = _mm256_movemask_epi8(gt) as u32;

            if (lsb_mask & gt_mask) > 0 {
                return Ordering::Greater;
            } else {
                return Ordering::Less;
            }
        }

        a_i += L;
        b_i += L;
    }

    Ordering::Equal
}

#[derive(Copy, Clone)]
struct MutPtr<const BYTES: usize>(*mut Int<BYTES>);
unsafe impl<const BYTES: usize> std::marker::Send for MutPtr<BYTES> {}
unsafe impl<const BYTES: usize> std::marker::Sync for MutPtr<BYTES> {}