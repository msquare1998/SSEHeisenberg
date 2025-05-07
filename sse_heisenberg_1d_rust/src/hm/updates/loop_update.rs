use crate::hm::HeisenbergModel;
use crate::{unsafe_get, unsafe_idx};

impl HeisenbergModel {
    pub fn make_vertex_lists(&mut self) {
        let mut op;
        let mut s0;
        let mut s1;
        let mut b_p;
        let mut v_leg0;
        let mut s0_v_last;
        let mut s1_v_last;

        // ::::::::::::::::::::::::::::::::::::::::::::::::::::
        //  Reset "vertex_list" and "v_first", "v_last"
        // ::::::::::::::::::::::::::::::::::::::::::::::::::::
        unsafe {
            for v in 0..4 * self.m {
                unsafe_idx!(self.vertex_list, v) = -1;
            }

            for s in 0..self.num_sites {
                unsafe_idx!(self.v_first, s) = -1;
                unsafe_idx!(self.v_last, s) = -1;
            }

        }

        // ::::::::::::::::::::::::::::::::::::
        //  Linking vertices sequentially
        // ::::::::::::::::::::::::::::::::::::
        unsafe {
            for p in 0..self.m {
                op = unsafe_get!(self.op_string, p);

                if op != -1 {
                    b_p = (op / 2) as usize;
                    s0 = unsafe_get!(self.b_sites, b_p, 0);
                    s1 = unsafe_get!(self.b_sites, b_p, 1);
                    v_leg0 = 4 * p;
                    s0_v_last = unsafe_get!(self.v_last, s0);
                    s1_v_last = unsafe_get!(self.v_last, s1);

                    // --------------------------------------------------------------------
                    if s0_v_last > -1 {
                        // which means this is not the 1st time we approach this site
                        //      then in the horizontal direction, we can connect "v_leg0" to a former leg
                        //      whose label is saved in vLast
                        unsafe_idx!(self.vertex_list, s0_v_last as usize) = v_leg0 as i32;
                        unsafe_idx!(self.vertex_list, v_leg0) = s0_v_last;
                    }

                    else  {
                        // which means this is the first time we approach this site
                        // then "v_leg0" becomes the first leg of this site
                        unsafe_idx!(self.v_first, s0) = v_leg0 as i32;
                    }

                    // of course, we have to update "v_Last" of this site for this "p"
                    unsafe_idx!(self.v_last, s0) = (v_leg0 + 2) as i32;

                    // --------------------------------------------------------------------
                    if s1_v_last > -1 {
                        unsafe_idx!(self.vertex_list, s1_v_last as usize) = (v_leg0 + 1) as i32;
                        unsafe_idx!(self.vertex_list, v_leg0 + 1) = s1_v_last;
                    }

                    else {
                        unsafe_idx!(self.v_first, s1) = (v_leg0 + 1) as i32;
                    }

                    unsafe_idx!(self.v_last, s1) = (v_leg0 + 3) as i32;
                }
            }
        }

        // ::::::::::::::::::::::::::::::::::::
        //  Make the PBC correction
        // ::::::::::::::::::::::::::::::::::::
        let mut s_v_first;
        let mut s_v_last;

        unsafe {
            for s in 0..self.num_sites {
                s_v_first = unsafe_get!(self.v_first, s);
                
                if s_v_first != -1 {
                    s_v_last = unsafe_get!(self.v_last, s);
                    unsafe_idx!(self.vertex_list, s_v_first as usize) = s_v_last;
                    unsafe_idx!(self.vertex_list, s_v_last as usize) = s_v_first;
                }
            }
        }
    }

    pub fn loop_update(&mut self) {
        self.make_vertex_lists();

        let mut v_head;
        let mut v_tail;

        // ::::::::::::::::::::::::::::::::::
        //  Update operators
        // ::::::::::::::::::::::::::::::::::
        unsafe {
            for v in (0..4 * self.m).step_by(2) {
                if unsafe_get!(self.vertex_list, v) < 0 {
                    continue;
                }

                v_head = v;

                if self.rand_prob() < 0.5 {
                    loop {
                        unsafe_idx!(self.op_string, v_head / 4) ^= 1;                   // diagonal (0) <---> off-diagonal (1)
                        unsafe_idx!(self.vertex_list, v_head) = -2;                     // mark this leg as "traversed" with "-2"
                        v_tail = v_head ^ 1;                                            // move to its neighbor
                        v_head = unsafe_get!(self.vertex_list, v_tail) as usize;        // go to another spin in another operator
                        unsafe_idx!(self.vertex_list, v_tail) = -2;                     // mark this leg as "traversed"

                        if v_head == v {
                            break;
                        }
                    }
                }

                else {
                    loop {
                        unsafe_idx!(self.vertex_list, v_head) = -1;
                        v_tail = v_head ^ 1;
                        v_head = unsafe_get!(self.vertex_list, v_tail) as usize;
                        unsafe_idx!(self.vertex_list, v_tail) = -1;
                        if v_head == v {
                            break;
                        }
                    }
                }
            }
        }

        // ::::::::::::::::::::::::::::::::::
        //  Update spins
        // ::::::::::::::::::::::::::::::::::
        unsafe {
            for i in 0..self.num_sites {
                // Flip the spins whose corresponding operator has been changed
                if unsafe_get!(self.v_first, i) != -1 {
                    if unsafe_idx!(self.vertex_list, unsafe_get!(self.v_first, i) as usize) == -2 {
                        unsafe_idx!(self.spins, i) *= -1;
                    }
                }

                // Flip those isolated spins with prob 50%
                else {
                    if self.rand_prob() < 0.5 {
                        unsafe_idx!(self.spins, i) *= -1;
                    }
                }
            }
        }
    }
}