use prng_mt::MT19937;
use crate::hm::HeisenbergModel;

impl HeisenbergModel {
    pub fn new(para_l: usize, para_beta: f64, para_seed: u32) -> Self {
        Self {
            // ----------------------------------------------------------------
            //  Basic params
            // ----------------------------------------------------------------
            l: para_l,          
            beta: para_beta,
            n: 0,   
            m: 10,

            // --------------------------------------------------------
            //  Lattice (PBC)
            // --------------------------------------------------------
            num_sites: para_l,
            num_bonds: para_l,            
            b_sites: (0..para_l).map(|i| vec![i, (i + 1) % para_l]).collect(),

            // --------------------------------------------------------
            //  Two frequently-used factors
            // --------------------------------------------------------
            prob_add_factor: 0.0,       // initialized in "init()"
            prob_remove_factor: 0.0,    // initialized in "init()"

            // ----------------------------------------
            //  Random number generator
            // ----------------------------------------
            rng: MT19937::new(para_seed),

            // ----------------------------------------
            //  Data structures for configuration
            // ----------------------------------------
            spins: Vec::new(),
            op_string: Vec::new(),
            v_first: Vec::new(),
            v_last: Vec::new(),
            vertex_list: Vec::new(),

            // ----------------------------------------
            //  For measurements
            // ----------------------------------------
            energy: 0.0,
            energy2: 0.0,
            heat_capacity: 0.0,
            zz_correlation: vec![0.0; para_l],
        }
    }

    pub fn initialization(&mut self){
        assert!(self.l % 2 == 0);
        self.prob_add_factor = 0.5 * self.beta * self.num_bonds as f64;
        self.prob_remove_factor = 1.0 / self.prob_add_factor;

        // Initialize other data structures
        self.spins = (0..self.num_sites).map(|_| if self.rand_prob() > 0.5 { 1 } else { -1 }).collect();
        self.op_string = vec![-1; self.m];
        self.v_first = vec![-1; self.num_sites];
        self.v_last = vec![-1; self.num_sites];
        self.vertex_list = vec![-1; 4 * self.m];
    }
}