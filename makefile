all:
	g++ BfastConvert.cpp -lm -o bfastconvert
	g++ bowtie2convert.cpp -lm -o bowtie2convert
	g++ pAlign.cpp -lm -lpthread -o align
	g++ cgal.cpp -lm -o cgal
