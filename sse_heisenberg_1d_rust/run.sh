l=4
beta=4
num_thm=20000
num_stat=10000
num_bins=5
cargo build --release
./target/release/sse_heisenberg_1d_rust $l $beta $num_thm $num_stat $num_bins