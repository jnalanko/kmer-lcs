# LCP for Succinct k-Spectra
This is the code for the paper [Longest Common Prefix Arrays for Succinct k-Spectra](https://arxiv.org/abs/2306.04850) by J. N. Alanko, E. Biagi,  S. J. Puglisi, SPIRE 2023. 

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
./benchmark index.sbwt linear out.txt
```

where linear is one of: naive, basic, superalphabet-2, superalphabet-4, linear.
