use prng_mt::mt19937::MT19937;
pub mod init;
pub mod updates;
pub mod measure;
pub mod random;

/* ==================================================================
    These two macros are to drop the safe checks in vectors
================================================================== */
#[macro_export]
macro_rules! unsafe_get {
    ($vec:expr, $idx:expr) => {
        $vec.get_unchecked($idx).clone()
    };

    ($vec:expr, $row:expr, $col:expr) => {
        $vec.get_unchecked($row).get_unchecked($col).clone()
    };
}

#[macro_export]
macro_rules! unsafe_idx {
    ($vec:expr, $idx:expr) => {
        *$vec.get_unchecked_mut($idx)
    };

    ($vec:expr, $row:expr, $col:expr) => {
        *$vec.get_unchecked_mut($row).get_unchecked_mut($col)
    };
}

pub struct HeisenbergModel {
    // ----------------------------------------------------------------
    //  Basic params
    // ----------------------------------------------------------------
    pub l: usize,               // length of the 1D ring
    pub beta: f64,              // the inverse temperature

    n: usize,                   // number of non-identity operators
    pub m: usize,                   // truncation order of the series

    // --------------------------------------------------------
    //  Lattice
    // --------------------------------------------------------
    pub num_sites: usize,             // number of sites
    num_bonds: usize,                 // number of bonds
    b_sites: Vec<Vec<usize>>,         // record the linking of the lattice

    // --------------------------------------------------------
    //  Two frequently-used factors
    // --------------------------------------------------------
    prob_add_factor: f64,
    prob_remove_factor: f64,

    // -------------------------------------------------
    //  PRNG with MT19937
    // -------------------------------------------------
    rng: MT19937,

    // -----------------------------------------------
    //  Data structures for configuration updates
    // -----------------------------------------------
    spins: Vec<i32>,
    op_string: Vec<i32>,
    v_first: Vec<i32>,
    v_last: Vec<i32>,
    vertex_list: Vec<i32>,

    // ----------------------------------------
    //  Measurements related
    // ----------------------------------------
    pub energy: f64,
    energy2: f64,              
    pub heat_capacity: f64,
    pub zz_correlation: Vec<f64>,   
}