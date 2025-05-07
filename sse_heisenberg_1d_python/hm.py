import random
import numpy as np

class Heisenberg:
    def __init__(self, para_l, para_beta, para_seed):
        # ------------------------------------------------
        #   Basic params
        # ------------------------------------------------
        random.seed(para_seed)      # read the seed
        assert para_l % 2 == 0      # bipartite condition
        self.l = para_l             # length of the 1D ring
        self.beta = para_beta       # the inverse temperature
        self.n = 0                  # number of non-identity operators
        self.m = 10                 # truncation order of the series

        # ------------------------------------------------
        #   Lattice (PBC)
        # ------------------------------------------------
        self.num_sites = para_l
        self.num_bonds = para_l
        self.b_sites = [[i, (i + 1) % para_l] for i in range(para_l)]

        # ------------------------------------------------
        #   Two frequently-used factors
        # ------------------------------------------------
        self.prob_add_factor = 0.5 * self.beta * self.num_bonds
        self.prob_remove_factor = 1.0 / self.prob_add_factor

        # ------------------------------------------------
        #   Data structures for configuration
        # ------------------------------------------------
        self.spins = [1 if self.rand_prob() > 0.5 else -1 for _ in range(self.l)]
        self.op_string = np.full(self.m, -1)
        self.v_first = np.full(self.num_sites, -1)
        self.v_last = np.full(self.num_sites, -1)
        self.vertex_list = np.full(4 * self.m, -1)

        # ------------------------------------------------
        #   Measurements related
        # ------------------------------------------------
        self.energy = 0.0
        self.zz_correlation = [0.0 for _ in range(self.l)]

    def mc_update(self):
        self.diagonal_update()
        self.loop_update()

    def diagonal_update(self):
        # ---------------------------------------------------------------------------
        #   I consider:
        #         op_string[p] := 2 * b[p] + a[p]
        #   where
        #         for the identity operator, op_string[p] = -1, and a[p] = -1;
        #         for the diagonal operator, a[p] = 0;
        #         for the off-diagonal operator, a[p] = 1;
        # ---------------------------------------------------------------------------
        for t in range(self.m):
            op = self.op_string[t]

            # Get identity ==> consider adding a diag-operator
            if op < 0:
                new_bond = self.rand_bond()
                if self.spins[self.b_sites[new_bond][0]] != self.spins[self.b_sites[new_bond][1]]:
                    if self.prob_add_factor >= (self.m - self.n) or self.prob_add_factor >= self.rand_prob() * (self.m - self.n):
                        self.op_string[t] = 2 * new_bond
                        self.n += 1

            # Get diag ==> consider removing it
            elif op % 2 == 0:
                if self.prob_remove_factor * (self.m - self.n + 1) >= 1.0 or self.rand_prob() <= self.prob_remove_factor * (self.m - self.n + 1):
                    self.op_string[t] = -1
                    self.n -= 1

            # Meet off-diag ==> propagate the spins only
            else:
                the_bond = int(op / 2)
                self.spins[self.b_sites[the_bond][0]] *= -1
                self.spins[self.b_sites[the_bond][1]] *= -1

    def adjust_m(self):
        new_m = self.n + int(self.n / 3)

        if self.m < new_m:
            self.op_string = np.append(self.op_string, [-1] * (new_m - self.m))
            self.m = new_m
            self.vertex_list = np.full(4 * self.m, -1)

    def make_vertex_list(self):
        # Reset "vertex_list" and "v_first", "v_last"
        for v in range(4 * self.m):
            self.vertex_list[v] = -1
        for s in range(self.num_sites):
            self.v_first[s] = -1
            self.v_last[s] = -1

        # Linking vertices sequentially
        for t in range(self.m):
            op = self.op_string[t]
            if op != -1:
                b_p = int(op / 2)
                s0 = self.b_sites[b_p][0]
                s1 = self.b_sites[b_p][1]
                v_leg0 = 4 * t

                s0_v_last = self.v_last[s0]
                if s0_v_last > -1:
                    # which means this is not the 1st time we approach this site
                    #   then in the time direction, we can connect "v_leg0" to a former leg
                    self.vertex_list[s0_v_last] = v_leg0
                    self.vertex_list[v_leg0] = s0_v_last
                else:
                    self.v_first[s0] = v_leg0
                self.v_last[s0] = v_leg0 + 2

                s1_v_last = self.v_last[s1]
                if s1_v_last > -1:
                    self.vertex_list[s1_v_last] = v_leg0 + 1
                    self.vertex_list[v_leg0 + 1] = s1_v_last
                else:
                    self.v_first[s1] = v_leg0 + 1
                self.v_last[s1] = v_leg0 + 3

        # Make the PBC correction
        for s in range(self.num_sites):
            s_v_first = self.v_first[s]

            if s_v_first != -1:
                s_v_last = self.v_last[s]
                self.vertex_list[s_v_last] = s_v_first
                self.vertex_list[s_v_first] = s_v_last

    def loop_update(self):
        self.make_vertex_list()

        # ---------------------------------------
        #   Update operators
        # ---------------------------------------
        for v in range(0, 4 * self.m, 2):
            if self.vertex_list[v] < 0:
                continue        # skip null legs and the legs that have been processed

            v0 = v
            if self.rand_prob() <= 0.5:
                while True:
                    self.op_string[int(v0 / 4)] ^= 1
                    self.vertex_list[v0] = -2
                    v1 = v0 ^ 1
                    v0 = self.vertex_list[v1]
                    self.vertex_list[v1] = -2

                    if v0 == v:
                        break

            else:
                while True:
                    self.vertex_list[v0] = -1
                    v1 = v0 ^ 1
                    v0 = self.vertex_list[v1]
                    self.vertex_list[v1] = -1
                    if v0 == v:
                        break

        # ---------------------------------------
        #   Update spins
        # ---------------------------------------
        for s in range(self.num_sites):
            if self.v_first[s] != -1:
                if self.vertex_list[self.v_first[s]] == -2:
                    self.spins[s] *= -1

            else:
                if self.rand_prob() <= 0.5:
                    self.spins[s] *= -1

    def rand_prob(self):
        return random.random()

    def rand_bond(self):
        return random.randrange(self.num_bonds)

    def get_zz_correlation(self, s_i, s_j):
        return self.spins[s_i] * self.spins[s_j]

    def init_measure(self):
        self.energy = 0.0
        for s in range(self.num_sites):
            self.zz_correlation[s] = 0.0

    def measure(self):
        self.energy += float(self.n)
        for s in range(self.num_sites):
            self.zz_correlation[s] += self.get_zz_correlation(0, s)

    def statisticize(self, num_samples):
        self.energy /= float(-self.beta * num_samples)
        for s in range(self.num_sites):
            self.zz_correlation[s] /= float(num_samples)