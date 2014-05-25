novo_muta
=========

@author Melissa Ip

An implementation of probabilistic methods for detecting de novo mutations from nuclear family genome data.

See Cartwright et al.: Family-Based Method for Capturing De Novo Mutations:
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3728889/

This project uses C++11, Eigen library, and GNU Scientific Library. It compiles using Clang.

Download Eigen and put the source code in the same directory as the novo_muta repository.
http://eigen.tuxfamily.org/

Download GSL.
http://www.gnu.org/software/gsl/

There are currently two versions of the trio model. Master is the trio model that uses customized Dirichlet-multinomial approximations. The infinite sites model branch is the trio model that uses simpler multinomial approximations and an infinite sites model. The simulation program is based on the trio model that uses the Dirichlet-multinomial.