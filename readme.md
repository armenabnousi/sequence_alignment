#DNA/protein sequence alignment</br>
##Smith-Waterman local and Needleman-Wunsch global algorithms with affine gap penalty</br>
</br>
This library provides a function for running Smith-Waterman or Needleman-Wunsch algorihtms for local/global sequence alginment. The scoring atrix should be provided as a text file as described below.</br>
</br>
**Usage:**</br>
To use this library you should declare an instance of `SeqAlign` object providing the name and path to the scoring matrix, a gap opening penalty and a gap extension penalty. Then you can perform the alignment by calling the `align` function, passing in the sequences to be aligned and the desired alignment type (local/global). This function does not return anything. Instead you can use the `print_alignment` function to print the alinged sequences or use the `get_alignment` function to receive the alignment score and the alignment start and end positions for the local alignments [start, end). Example:
```c++
SeqAlign seq_align("scoring_matrix/blosum62", -5, -1);
std::string seq1 = "NYLKDCIRVQTDLAKSHEYQ";
std::string seq2 = "GSYSAFGERDGNYLKDGIGNTWL";
seq_align.align(seq1, seq2, "local");
seq_align.print_alignment();
int seq1_start, seq2_start, seq1_end, seq2_end, alignment_score;
double score_percentage;
std::string alignment = seq_align.get_alignment(&seq1_start, &seq1_end, &seq2_start, @seq2_end, &alignment_score, &score_percentage);
```
To compute the score percentage, first we computed the maximum possible alignment score by picking the sequence with shorter length and performing a global alignment by itself. Then the score percentage is computed by dividing the original alignment score by the maximum possible score.</br>
</br>
**Scoring Matrix:**</br>
A file containing the scoring matrix should be provided. An example of such a file is provided in 'scoring_matrix/blosum62'. This file should contain the number of letters in the alphabet on the first line. The second line contains the letters of the alphabet in the same order that they appear in the rows/columns of the matrix. Note that these letters should be printed consecutively with no space in between them.</br> Starting from the third line the scores should be included in the file, intuitively, separated by spaces for representing columns and new lines for rows.</br>

