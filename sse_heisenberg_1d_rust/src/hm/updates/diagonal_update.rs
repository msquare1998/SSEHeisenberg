use crate::hm::HeisenbergModel;
use crate::{unsafe_get, unsafe_idx};

/* ==================================================================
    I consider:
        op_string[p] := 2 * b[p] + a[p]
    where
        for the identity operator, op_string[p] = -1, and a[p] = -1;
        for the diagonal operator, a[p] = 0;
        for the off-diagonal operator, a[p] = 1;
================================================================== */

impl HeisenbergModel {
    pub fn diag_update(&mut self) {
        let mut new_bond;
        let mut the_bond;
        let mut op;

        unsafe {
            for p in 0..self.m {
                op = unsafe_get!(self.op_string, p);

                // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                //  Get identity ==> consider adding a diag-operator
                // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                if op < 0 {
                    new_bond = self.rand_bond();

                    if unsafe_get!(self.spins, unsafe_get!(self.b_sites, new_bond, 0)) != unsafe_get!(self.spins, unsafe_get!(self.b_sites, new_bond, 1)) {
                        if (self.prob_add_factor >= (self.m - self.n) as f64)
                            ||
                            (self.prob_add_factor >= self.rand_prob() * (self.m - self.n) as f64) {
                            unsafe_idx!(self.op_string, p) = 2 * new_bond as i32;
                            self.n += 1;
                        }
                    }
                }

                // :::::::::::::::::::::::::::::::::::::::::::::::::::
                //  Get diag ==> consider removing it
                // :::::::::::::::::::::::::::::::::::::::::::::::::::
                else if op % 2 == 0 {
                    if (self.prob_remove_factor * (self.m - self.n + 1) as f64 >= 1.0)
                        ||
                        (self.rand_prob() <= self.prob_remove_factor * (self.m - self.n + 1) as f64) {
                        unsafe_idx!(self.op_string, p) = -1;
                        self.n -= 1;
                    }
                }

                // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                //  Meet off-diag ==> propagate the spins only
                // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                else {
                    the_bond = (op / 2) as usize;
                    unsafe_idx!(self.spins, unsafe_get!(self.b_sites, the_bond, 0)) *= -1;
                    unsafe_idx!(self.spins, unsafe_get!(self.b_sites, the_bond, 1)) *= -1;
                }
            }
        }
    }
}