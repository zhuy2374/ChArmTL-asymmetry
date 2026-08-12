#
# @Author  : zhuy
# @Date    : 2026/1/10 18:07
# @Desc    : Sequence similarity calculation
#


class SimilarytiCalculator:
    """
    Sequence similarity calculator.
    """

    @staticmethod
    def calculate_similarity(seq1, seq2, k):
        """
            Calculate the Levenshtein similarity of two sequence strings over the overlapping region.
            The similarity is defined as 1 - (edit_distance / overlap_length), where
            overlap_length = min(len(seq1), len(seq2)).
            If overlap_length == 0, return 0.0 to preserve the original behavior.

            Parameters:
                seq1, seq2: Sequence strings to compare.

            Returns:
                float: Similarity. 1.0 means identical, and 0.0 means fully different or no overlap.

            Examples:
                >>> calculate_similarityV2("ACGT", "ACGG")
                0.75   # edit distance=1, overlap=4 -> 1 - 1/4 = 0.75
                >>> calculate_similarityV2("AAAA", "TTTT")
                0.0    # edit distance=4, overlap=4 -> 1 - 4/4 = 0.0
                >>> calculate_similarityV2("", "ACGT")
                0.0    # overlap = 0 -> returns 0.0
            """
        # Use the overlapping region, consistent with the original behavior.
        overlap = min(len(seq1), len(seq2))
        if overlap == 0:
            return 0.0

        a = seq1[:overlap]
        b = seq2[:overlap]

        # Fast equality check.
        if a == b:
            return 1.0

        # Keep b as the shorter sequence to reduce memory usage.
        if len(b) > len(a):
            a, b = b, a  # Ensure len(a) >= len(b).

        # Calculate Levenshtein edit distance with rolling arrays.
        prev_row = list(range(len(b) + 1))
        for i, ca in enumerate(a, start=1):
            cur_row = [i] + [0] * len(b)
            for j, cb in enumerate(b, start=1):
                cost = 0 if ca == cb else 1
                # Delete: prev_row[j] + 1
                # Insert: cur_row[j-1] + 1
                # Replace: prev_row[j-1] + cost
                cur_row[j] = min(prev_row[j] + 1, cur_row[j - 1] + 1, prev_row[j - 1] + cost)
            prev_row = cur_row

        edit_distance = prev_row[-1]
        relative_distance = edit_distance / overlap  # Same normalization as the original code.
        similarity = 1.0 - relative_distance

        # Clamp numeric bounds.
        if similarity < 0.0:
            similarity = 0.0
        elif similarity > 1.0:
            similarity = 1.0

        return round(float(similarity), k)
