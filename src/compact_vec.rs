#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(transparent)]
pub struct Int<const BYTES: usize>([u8; BYTES]);

impl<const BYTES: usize> Int<BYTES> {
    #[inline(always)]
    pub fn get_usize(&self) -> usize {
        let mut res = 0u64;
        unsafe {
            std::ptr::copy_nonoverlapping(self.0.as_ptr(), &mut res as *mut _ as _, BYTES);
        }
        res as usize
    }

    #[inline(always)]
    pub fn set_usize(&mut self, val: usize) {
        unsafe {
            std::ptr::copy_nonoverlapping(&val as *const _ as _, self.0.as_mut_ptr(), BYTES);
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CompactVec<const BYTES: usize> {
    data: Vec<Int<BYTES>>,
}

impl<const BYTES: usize> CompactVec<BYTES> {
    pub fn new(len: usize) -> Self {
        assert!(BYTES <= 8);

        Self {
            data: vec![Int([0u8; BYTES]); len + if len > 0 { 8 } else { 0 }],
        }
    }
}

impl<const BYTES: usize> std::ops::Deref for CompactVec<BYTES> {
    type Target = [Int<BYTES>];

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl<const BYTES: usize> std::ops::DerefMut for CompactVec<BYTES> {
    fn deref_mut(&mut self) -> &mut [Int<BYTES>] {
        &mut self.data
    }
}
