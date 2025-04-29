# operator_cg
This repository contains MATLAB codes implementing a Conjugate Gradient (CG) method applied to an ordinary differential equation:

-(k(x) u'(x))' + c(x) u(x) = f(x) on Ω = (-1,1),  
f ∈ L²(Ω), k, c ∈ L^∞(Ω), k(x) > 0, c(x) ≥ 0

The code is built using the [Chebfun](https://www.chebfun.org/) toolbox.
To use these codes, you will need to have MATLAB and the Chebfun toolbox installed.

The repository contains the following functions:
- cgh1 - CG method, when H₀¹(Ω) inner product is considered 
- cgh1_ab Generalization of cgh1 to the case Ω = (a, b)
- cgl2 CG method, when L² inner product is considered
- myprojection - used in cgl2, computes projection onto polynomial space spanned by shifted Chebyshev polynomials
- pcg2 – Modified version of the Chebfun `pcg` function that saves approximations of the solution at each iteration  
  - Based on: [Chebfun PCG source](https://github.com/chebfun/chebfun/blob/master/%40chebop/pcg.m)  
  - License: BSD 2-Clause (see header in `pcg2.m` for details)

The repository also contains four examples that were used for experiments in the diploma thesis.
