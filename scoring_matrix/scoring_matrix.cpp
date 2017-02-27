#include "scoring_matrix.h"

ScoringMatrix::ScoringMatrix(std::string matrix_name) : matrix(NULL), name_to_num(NULL) {
	std::string matrix_path = "/home/armen.abnousi/libs/seq_align/" + matrix_name;
	read_matrix(matrix_path);
}

ScoringMatrix::~ScoringMatrix() {
	if (matrix) {
		for (int i = 0; i < input_size; i++) {
			delete[] matrix[i];
		}
		delete[] matrix;
	}
	if (name_to_num) delete[] name_to_num;
}

void ScoringMatrix::read_matrix(std::string filename) {
	std::ifstream infile(filename);
	std::string line;
//	std::cout << "herenow" <<std::endl;
	infile >> input_size;
//	std::cout << input_size << std::endl;
	std::string num_to_name;
	std::getline(infile, num_to_name);
	std::getline(infile, num_to_name);
//	std::cout << "these are the alphabet:" << num_to_name << std::endl;
	name_to_num = new int[ALPHABET_SIZE];
	for (int i = 0; i < input_size; i++) {
		name_to_num[num_to_name[i] - 'A'] = i;
//		std::cout << "setting " << num_to_name[i] - 'A' << " to " << i << std::endl;
	}
	matrix = new int*[input_size];
	for (int i = 0; i < input_size; i++) {
		matrix[i] = new int[input_size];
	}
	int row_num = 0;
	while(std::getline(infile, line)) {
		//std::cout << "this is the line: " << line << std::endl;
		std::istringstream iss(line);
		for (int j = 0; j < input_size; j++) {
			iss >> matrix[row_num][j];
		}
		row_num++;
	}
//	std::cout << "matrix read" << std::endl;
}

int ScoringMatrix::get_score(char a, char b) {
//	std::cout << "checking " << a << " " << b << " " << a - 'A' << " " << b - 'A' << std::endl;
//	std::cout << name_to_num[a- 'A'] << " " << name_to_num[b - 'A'] << std::endl;
//	std::cout << "result is " << matrix[name_to_num[a - 'A']][name_to_num[b - 'A']] << std::endl;
	return matrix[name_to_num[a - 'A']][name_to_num[b - 'A']];
}
