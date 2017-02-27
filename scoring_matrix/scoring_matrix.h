#include <sstream>
#include <fstream>
#include <iostream>
#define ALPHABET_SIZE 26

class ScoringMatrix {
	int* name_to_num;
	int** matrix;
	int input_size;
	void read_matrix(std::string filename);
public:
	ScoringMatrix(std::string matrix_name);	
	ScoringMatrix() : matrix(NULL), name_to_num(NULL) { }
	~ScoringMatrix();
	int get_score(char a, char b);
};

