#DNA/protein sequence alignment</br>
##Smith-Waterman local and Needleman-Wunsch global alignment algorithms with affine gap penalty</br>
</br>
This library provides implementations for Smith-Waterman or Needleman-Wunsch algorihtms for local/global sequence alginment. The scoring matrix should be provided as a text file as described below.</br>
</br>
**Usage:**</br>
To use this library you should declare an instance of `SeqAlign` object providing the name and path to the scoring matrix, a gap opening penalty and a gap extension penalty. Then you can perform the alignment by calling the `align` function, passing in the sequences to be aligned and the desired alignment type (local/global). This function does not return anything, rather it only performs the alignment and the result can be fetched back through `print_alignment` and `get_alignment` functions. Example:
```c++
SeqAlign seq_align("scoring_matrix/blosum62", -5, -1);
std::string seq1 = "NYLKDCIRVQTDLAKSHEYQ";
std::string seq2 = "GSYSAFGERDGNYLKDGIGNTWL";
seq_align.align(seq1, seq2, "local");
seq_align.print_alignment();
int seq1_start, seq2_start, seq1_end, seq2_end, alignment_score;
double score_percentage;
std::string alignment = seq_align.get_alignment(&seq1_start, &seq1_end, &seq2_start, &seq2_end, &alignment_score, &score_percentage);
```
The `get_alignment` function returns a std::string of the alignment where '|' represents exact match, '-' represents a mismatch, '1' represents a gap in sequence1 (first sequence passed to `align` function), and '2' represents a gap in sequence2 (second sequence passed to `align` function). The parameters of this function are all used as output arguments, that will return the alignment score, start and end positions of the alignment on the two sequences, and score percentage and identity count as explained below.</br>
The `print_alignment` function only prints the sequence of matches, mimatches and gaps ('|', '-', '1', and '2' as explained above for `get_alignment` function.</br>
**score percentage**</br>
To compute the score percentage, first we computed the maximum possible alignment score by picking the sequence with shorter length and performing a global alignment by itself. Then the score percentage is computed by dividing the original alignment score by the maximum possible score.</br>
**identity count**</br>
This is the number of exact matches in the alignment found between the two given sequences.</br>
</br>
**Scoring Matrix:**</br>
A file containing the scoring matrix should be provided. An example of such a file is provided in 'scoring_matrix/blosum62'. This file should contain the number of letters in the alphabet on the first line. The second line contains the letters of the alphabet in the same order that they appear in the rows/columns of the matrix. Note that these letters should be printed consecutively with no space in between them.</br> Starting from the third line the scores should be included in the file, intuitively, separated by spaces for representing columns and new lines for rows.</br>

