CXX=clang++

all: a.out

a.out: main.cpp
	$(CXX) main.cpp

clean:
	rm -f a.out

re: clean all
