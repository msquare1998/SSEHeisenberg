use crate::hm::HeisenbergModel;

impl HeisenbergModel {
    pub fn rand_prob(&mut self) -> f64 {
        (self.rng.next() % u32::MAX) as f64 / (u32::MAX as f64)
    }

    pub fn rand_bond(&mut self) -> usize {
        self.rng.next() as usize % self.num_bonds
    }
}