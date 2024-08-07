# LCP for Succinct k-Spectra
This is the code for the paper [Longest Common Prefix Arrays for Succinct k-Spectra](https://link.springer.com/chapter/10.1007/978-3-031-43980-3_1) by J. N. Alanko, E. Biagi,  S. J. Puglisi, SPIRE 2023. 

# Building

```
git submodule update --init --recursive
cd SBWT/build

# Here we are compiling with gcc-10 and g++-10. Replace these with whatever
# version of g++ you have. Clang will not work.
cmake .. -DCMAKE_C_COMPILER=$(which gcc-10) -DCMAKE_CXX_COMPILER=$(which g++-10)
make -j4

cd ../..
make benchmark --always-make CXX=g++-10
```

# Running

The code takes a plain-matrix sbwt file as input. You can generate one by running:

```
./SBWT/build/bin/sbwt build -i SBWT/example_data/coli3.fna -o index.sbwt -k 30
```

Then, you can run the benchmark with:

```
./benchmark index.sbwt <variant> out.txt
```

The possible variants are: naive, basic, superalphabet-2, superalphabet-4, linear, basic-parallel-(#threads).

By adding a fourth argument it is possible to save the LCS to a file. The following command saves the LCS to the file index-lcs.sdsl

```
./benchmark index.sbwt <variant> out.txt index-lcs
```

