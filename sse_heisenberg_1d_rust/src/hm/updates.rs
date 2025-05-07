use crate::hm::HeisenbergModel;
pub mod diagonal_update;
pub mod loop_update;

impl HeisenbergModel {
    pub fn mc_update(&mut self) {
        self.diag_update();
        self.loop_update();
    }

    pub fn adjust_m(&mut self) {
        let new_m = self.n + self.n / 3;

        if self.m < new_m {
            self.op_string.extend(vec![-1; new_m - self.m]);
            self.m = new_m;
            self.vertex_list = vec![-1; 4 * new_m];
        }
    }
}