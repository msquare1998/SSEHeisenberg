# Stochastic series expansion quantum Monte Carlo algorithm for S=1/2 AFM Heisenberg model
This project reproduces the same algorithm as that in https://physics.bu.edu/~sandvik/programs/ssebasic/ssebasic.html

## About Python
The Python version is just for beginners to learn the algorithm. Practically, the Python program is much slower in Monte Carlo tasks. 

## Unsafe macro in Rust
You can replace all the unsafe blocks in my code with the corresponding safe versions. I had some tests on my MacBook Pro (M1), and found that it would lead to around 5% performance loss if you do this. I suggest the safe version, which is reassuring (otherwise why don't we use C/C++ instead?), though I upload the unsafe version here for those who need ultimate performance.
