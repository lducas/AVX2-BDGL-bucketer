test_lsh: fht_lsh.cpp fht_lsh.h test_lsh.cpp
	g++ -Wall -Wextra -march=native -mtune=native -ftree-vectorize -funroll-loops -O3  -std=c++11 -c fht_lsh.cpp -o fht_lsh.o
	g++ -Wall -Wextra -march=native -mtune=native -ftree-vectorize -funroll-loops -O3  -std=c++11 -o test_lsh test_lsh.cpp fht_lsh.o
