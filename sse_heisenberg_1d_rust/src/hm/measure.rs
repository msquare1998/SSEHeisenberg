use crate::hm::HeisenbergModel;
use crate::{unsafe_get, unsafe_idx};

impl HeisenbergModel {
    pub fn get_zz_correlation(&self, s_i: usize, s_j: usize) -> i32 {
        self.spins[s_i] * self.spins[s_j]
    }
}

impl HeisenbergModel {
    pub fn init_measure(&mut self) {
        self.energy = 0.0;
        self.energy2 = 0.0;

        unsafe {
            for i in 0..self.num_sites {
                unsafe_idx!(self.zz_correlation, i) = 0.0;
            }
        }
    }

    pub fn measure(&mut self) {
        unsafe {
            for s in 0..self.num_sites {
                unsafe_idx!(self.zz_correlation, s) += self.get_zz_correlation(0, s) as f64;
            }
        }

        self.energy += self.n as f64;
        self.energy2 += (self.n * self.n) as f64;
    }

    pub fn statisticize(&mut self, mc_steps: usize) {
        unsafe {
            for s in 0..self.num_sites {
                unsafe_idx!(self.zz_correlation, s) = unsafe_get!(self.zz_correlation, s) / mc_steps as f64;
            }
        }

        self.energy /= mc_steps as f64;
        self.energy2 /= mc_steps as f64;

        // C = <n^2> - <n>^2 - <n>
        self.heat_capacity = self.energy2 - self.energy * self.energy - self.energy;
        // E = -<n> / beta
        self.energy /= -self.beta;
    }
}