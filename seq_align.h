#pragma once
#include <limits.h>
#include <iostream>
#include "scoring_matrix/scoring_matrix.h"

class SequenceAlignment {
	private:
		int score;
		double score_percent;
		int exact_match_count;
		std::string type;
		std::string scoring_scheme;
		ScoringMatrix* scores;
		std::string alignment;
		int gap_opening_penalty;
		int gap_penalty;
		int seq1_start, seq1_end, seq2_start, seq2_end;
		std::string seq1, seq2;
	public:
		SequenceAlignment (std::string matrix, int gap_o, int gap);
		void align(std::string seq1, std::string seq2, std::string type);
		void print_alignment();
		std::string get_alignment(int* s1_start, int* s1_end, int* s2_start, int* s2_end, int* score, double* score_percent, 
						int* identity_count);
		void fill_matrix(std::string seq1, std::string seq2, int** matrix, int** s1_gaps, int** s2_gaps, 
					std::string type, int nrow, int ncol, int* max_i, int* max_j);
		void backtrack_matrix(std::string seq1, std::string seq2, int** matrix, int** s1_gaps, int** s2_gaps, 
					std::string type, int max_i, int max_j);
		int get_max(int a, int b, int c);
		int get_max(int a, int b);
		void print_matrix(int** matrix, int nrow, int ncol);
		int compute_max_score(std::string seq);
};
