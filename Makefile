CXXFLAGS = -O2 -march=native -Wall -Wextra -std=gnu++26 -fopenmp

all: main

main: main.cc
	$(CXX) $(CXXFLAGS) -o $@ $< -MMD

render: main
	./main > result.ppm

clean:
	rm -f result.ppm main main.d

-include main.d

.PHONY: all render clean
