.PHONY: main reduce_sbwt_order

benchmark:
	g++ src/benchmark.cpp SBWT/build/libsbwt_static.a SBWT/build/external/sdsl-lite/build/lib/libsdsl.a -std=c++20 -I ./SBWT/sdsl-lite/include/ -O3 -I include -I ./SBWT/include -I ./SBWT/include/sbwt -I SBWT/build/external/sdsl-lite/build/external/libdivsufsort/include/ -g -o benchmark -Wno-deprecated-declarations -DNDEBUG	

reduce_sbwt_order:
	g++ src/reduce_sbwt_order.cpp SBWT/build/libsbwt_static.a SBWT/build/external/sdsl-lite/build/lib/libsdsl.a -std=c++20 -I ./SBWT/sdsl-lite/include/ -O3 -I include -I ./SBWT/include -I ./SBWT/include/sbwt -I SBWT/build/external/sdsl-lite/build/external/libdivsufsort/include/ -g -o reduce_sbwt_order -Wno-deprecated-declarations -DNDEBUG	