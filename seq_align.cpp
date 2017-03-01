#include "seq_align.h"

SequenceAlignment::SequenceAlignment(std::string matrix, int gap_o, int gap) : 
					scoring_scheme(matrix), gap_opening_penalty(gap_o), gap_penalty(gap) {
		scores = new ScoringMatrix(scoring_scheme);
		seq1_start = 0;
		seq2_start = 0;
		seq1_end = 0;
		seq2_end = 0;
		exact_match_count = 0;
		//align(seq1, seq2, type);
}

void SequenceAlignment::align(std::string seq1, std::string seq2, std::string type) {
	this -> seq1 = seq1;
	this -> seq2 = seq2;
	this -> type = type;
	seq1_start = 0;
	seq2_start = 0;
	seq1_end = seq1.length() - 1;
	seq2_end = seq2.length() - 1;
	int nrow = seq1.length() + 1;
	int ncol = seq2.length() + 1;
	int** matrix = new int*[nrow];
	int** s1_gap = new int*[nrow];
	int** s2_gap = new int*[nrow];
	for (int i = 0; i < nrow; i++) {
		matrix[i] = new int[ncol];
		s1_gap[i] = new int[ncol];
		s2_gap[i] = new int[ncol];
	}
	int max_i, max_j;
	fill_matrix(seq1, seq2, matrix, s1_gap, s2_gap, type, nrow, ncol, &max_i, &max_j);
//	std::cout << "matrix filled" << std::endl;
	score = get_max(matrix[max_i][max_j], s1_gap[max_i][max_j], s2_gap[max_i][max_j]);
//	print_matrix(matrix, nrow, ncol);
//	print_matrix(s1_gap, nrow, ncol);
//	print_matrix(s2_gap, nrow, ncol);

	backtrack_matrix(seq1, seq2, matrix, s1_gap, s2_gap, type, max_i, max_j);
//	std::cout << "score = " << score << " " << matrix[seq1_end][seq2_end] << " " << seq1[seq1_end - 1] << " " << seq2[seq2_end - 1] << " " << seq2_end - seq2_start + 1 << " " << seq1_end - seq1_start + 1 << std::endl;
	int max_possible_score = seq1.length() > seq2.length() ? compute_max_score(seq1) : compute_max_score(seq2);
	score_percent = (double) score / max_possible_score;
//	std::cout << seq1_start << " " << seq1_end << " " << seq2_start << " " << seq2_end << std::endl;
	for (int i = 0; i < nrow; i++) {
		delete[] matrix[i];
		delete[] s1_gap[i];
		delete[] s2_gap[i];
	}
	delete[] matrix;
	delete[] s1_gap;
	delete[] s2_gap;
}

void SequenceAlignment::print_alignment() {}

std::string SequenceAlignment::get_alignment(int* s1_start, int* s1_end, int* s2_start, int* s2_end, int* score, double* score_percent, 
						int* identity_count) {
	*s1_start = this -> seq1_start;
	*s1_end = this -> seq1_end;
	*s2_start = this -> seq2_start;
	*s2_end = this -> seq2_end;
	*score = this -> score;
	*score_percent = this -> score_percent;
	*identity_count = exact_match_count;
	return alignment;
}

void SequenceAlignment::fill_matrix(std::string seq1, std::string seq2, int** matrix, int** s1_gaps, int** s2_gaps, 
					std::string type, int nrow, int ncol, int* max_i, int* max_j) { 
//	std::cout << "filling the matrix"  << std::endl;
	int min_possible_int = INT_MIN + 10000;
	for (int i = 0; i < nrow; i++) {
		matrix[i][0] = min_possible_int;
		
		//s1_gaps[i][0] = min_possible_int;
		//s2_gaps[i][0] = min_possible_int;
		s2_gaps[i][0] = min_possible_int;
		s1_gaps[i][0] = type == "local" ? 0 : gap_opening_penalty + (i * gap_penalty);
	}
	for (int j = 0; j < ncol; j++) {
		matrix[0][j] = min_possible_int;
		//s1_gaps[0][j] = min_possible_int;
		//s2_gaps[0][j] = min_possible_int;
		s1_gaps[0][j] = min_possible_int;
		s2_gaps[0][j] = type == "local" ? 0 : gap_opening_penalty + (j * gap_penalty);
	}
	matrix[0][0] = 0;
	s1_gaps[0][0] = min_possible_int;
	s2_gaps[0][0] = min_possible_int;
	char s1_char, s2_char;
//	std::cout << "all initialized" << std::endl;
	*max_i = nrow - 1;
	*max_j = ncol - 1;
	int max_score = 0;
	for (int i = 1; i < nrow; i++) {
		s1_char = seq1[i - 1];
		for (int j = 1; j < ncol; j++) {
			s2_char = seq2[j - 1];
			//std::cout << "checking i= " << i << " j= " << j << std::endl;
			//std::cout << "that gives us " << scores -> get_score(s1_char, s2_char) << std::endl;
			//std::cout << "getting score for " << s1_char << " and " << s2_char << std::endl;
			//std::cout << "that's : " << scores -> get_score(s1_char, s2_char) << std::endl;
			matrix[i][j] = scores -> get_score(s1_char, s2_char) + 
				get_max(matrix[i-1][j-1], s1_gaps[i-1][j-1], s2_gaps[i-1][j-1]);
			if (type == "local")
				if (matrix[i][j] < 0)
					matrix[i][j] = 0;	
			//std::cout << "computing max of " << gap_opening_penalty + gap_penalty + matrix[i-1][j] << " " << gap_penalty + s1_gaps[i-1][j] << " as: " << get_max(gap_opening_penalty + gap_penalty + matrix[i-1][j], gap_penalty + s1_gaps[i-1][j]) << std::endl;
			s1_gaps[i][j] = get_max(gap_opening_penalty + gap_penalty + matrix[i][j-1], 
						gap_penalty + s1_gaps[i][j-1],
						gap_opening_penalty + gap_penalty + s2_gaps[i][j-1]);
			s2_gaps[i][j] = get_max(gap_opening_penalty + gap_penalty + matrix[i-1][j],
						gap_opening_penalty + gap_penalty + s1_gaps[i-1][j],
						gap_penalty + s2_gaps[i-1][j]);
			int step_max = get_max(matrix[i][j], s1_gaps[i][j], s2_gaps[i][j]);
			if (type == "local" && max_score < step_max) {
				max_score = step_max;
				*max_i = i;
				*max_j = j;
			}
		}
	}
}

void SequenceAlignment::backtrack_matrix(std::string seq1, std::string seq2, int** matrix, int** s1_gaps, int** s2_gaps, 
						std::string type, int max_i, int max_j) {
	if (type == "local") {
		seq1_end = max_i;
		seq2_end = max_j;
	}
	int i = max_i;
	int j = max_j;
	alignment = "";
	exact_match_count = 0;
	int step_score = get_max(s1_gaps[i][j], s2_gaps[i][j], matrix[i][j]);
	while ( (type == "local" && step_score > 0) || (type == "global" && i != 0 && j != 0) ) {
//		std::cout << "score: " << step_score << " at:" << i << ", " << j << std::endl;
		if (step_score == get_max(gap_opening_penalty + gap_penalty + matrix[i-1][j], 
					gap_opening_penalty + gap_penalty + s1_gaps[i-1][j], 
					gap_penalty + s2_gaps[i-1][j])) {
			alignment += '2';
			if (step_score == gap_opening_penalty + gap_penalty + matrix[i-1][j]) step_score = matrix[i-1][j];
			else if (step_score == gap_opening_penalty + gap_penalty + s1_gaps[i-1][j]) step_score = s1_gaps[i-1][j];
			else step_score = s2_gaps[i-1][j];
			i--;
		} else if (step_score == get_max(gap_opening_penalty + gap_penalty + matrix[i][j-1],
						gap_penalty + s1_gaps[i][j-1],
						gap_opening_penalty + gap_penalty + s2_gaps[i][j-1])) {
			if (step_score == gap_opening_penalty + gap_penalty + matrix[i][j-1]) step_score = matrix[i][j-1];
			else if (step_score == gap_penalty + s1_gaps[i][j-1]) step_score = s1_gaps[i][j-1];
			else step_score = s2_gaps[i][j-1];
			alignment += '1';
			j--;
		} else if (step_score == scores -> get_score(seq1[i-1], seq2[j-1]) + 
					get_max(matrix[i-1][j-1], s1_gaps[i-1][j-1], s2_gaps[i-1][j-1])) {
			if (seq1[i-1] == seq2[j-1]) {
				alignment += '|';
				exact_match_count++;
			} else {
				alignment += '-';
			}
			int match_point = scores -> get_score(seq1[i-1], seq2[j-1]);
			if (step_score == match_point + matrix[i-1][j-1]) step_score = matrix[i-1][j-1];
			else if (step_score == match_point + s1_gaps[i-1][j-1]) step_score = s1_gaps[i-1][j-1];
			else step_score = s2_gaps[i-1][j-1];
			i--;
			j--;
		}
		//step_score = get_max(s1_gaps[i][j], s2_gaps[i][j], matrix[i][j]);
	}
	seq1_start = i;
	seq2_start = j;
	int length = alignment.length() - 1;
//	std::cout << "reverse alignment: " << alignment << std::endl;
	for (int i = 0; i < length - i; i++) {
		char temp = alignment[i];
		alignment[i] = alignment[length - i];
		alignment[length - i] = temp;
	}
//	std::cout << "alignmnet: " << alignment << std::endl;
}

int SequenceAlignment::get_max(int a, int b, int c) {
	int temp = a > b ? a : b;
	return get_max(temp, c);
}

int SequenceAlignment::get_max(int a, int b) {
	return a > b ? a : b;
}

void SequenceAlignment::print_matrix(int** matrix, int nrow, int ncol) {
	std::cout << std::endl;
	for (int i = 0; i < nrow; i++) {
		for (int j = 0; j < ncol; j++) {
			std::cout << '\t' << matrix[i][j];
		}
		std::cout << std::endl;
	}
}

int SequenceAlignment::compute_max_score(std::string seq) {
	int result = 0;
	for (int i = seq.length() - 1 ; i >= 0; i--) {
		result += scores -> get_score(seq[i], seq[i]);
	}
	return result;
}
