# AVX2-BDGL-bucketer

## An Implementation of Locality-Sensitive Hashing over the sphere à-la Becker-Ducas-Gama-Laarhoven, abusing the AVX2 instruction set

Contributors:
- Léo Ducas
- Wessel van Woerden
- Joe Rowell

### Disclaimer

This is not pretty code. 

### Description

This class implements a Locality-Sensitive-Hashing function over the sphere, roughly following the principle of the Becker-Ducas-Gama-Laarhoven algorithm.

**New directions in nearest neighbor searching with applications to lattice sieving**  
*Anja Becker and Léo Ducas and Nicolas Gama and Thijs Laarhoven*
https://eprint.iacr.org/2015/1128

Numerous liberties have been taken from the original algorithm:

- The function returns a fixed number M of buckets, rather than the set of buckets given by a fixed angle threshold. These are not necessarily exactly the best M buckets due to vectorization constraints.
- The subbucket directions of each block are not uniformly random; they are ternary vectors of weight 16 or 32, and are related by some underlying Hadamard structure over subsets of 16 or 32 buckets.

These adaptation permits severe AVX2 optimizations, discussed in more detailed in

**Advanced lattice Sieving on GPU, with Tensor Cores**  
*Leo Ducas and Marc Stevens and Wessel van Woerden*
(To appear)

### Interface

The main class is `ProductLSH`

Initialization:
``` C++
lsh = ProductLSH(n, b, C, M, seed);
```
- `n` is the dimension of the space
- `b` is the number of ``blocks''
- `C` is the *desired* number of buckets (i.e. possible value for the hash function). Due to block constraint C = 2^(b-1) * C_1 * C_2 * ... * C_b, the actual number of buckets may be somewhat lower. Actual value accessible via `lsh.codesize'
- `M` is the number of desired output hashvalues.
- `seed` for controlled pseudo-randomness.

The class provides a single function:
``` C++
int32_t res[M]; 
lsh.hash(v, res);
```
which hash the vector v, and store the M outputs in res.


Hashing time is essentially proportional to C^1/b. Quality of the bucket decrease as b increases; i.e., for the same C, vectors in the same bucket may be further appart for larger b.


Remarks:
- The hash function is symmetric: `hash(v) = hash(-v)`
- Only implemented for parameters such that n/b is in the range [17,128].

### Testing/Benchmarking

Compile and run `test_lsh` for some statistics with various parameters and dimension. 
Statistic interpretation not included.

