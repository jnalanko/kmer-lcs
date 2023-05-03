# Building

```
git submodule update --init --recursive
cd SBWT/build
cmake ..
make -j4
cd ../..
make main --always-make
```

# Running

The code takes a plain-matrix sbwt file as input. You can generate one by running:

```
./SBWT/build/bin/sbwt build -i SBWT/example_data/coli3.fna -o index.sbwt -k 30
```

Then, you can run the LCS construction algorithms with:

```
./main index.sbwt
```