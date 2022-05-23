# DFR data for near-codewords

This folder contains data on decoding failure rates for near-codewords for `r = 659`. The JSON files contain data in the following format:

* `"set"`: which of the three defined classes of near-codeword is used. This is either `C`, `N`, or `2N`.
* `"bits"`: The security parameter, denoted `λ` in the paper.
* `"r"`: The circulant block length.
* `"d"`: The column weight, always set to 15.
* `"t"`: The error vector weight, always set to 18.
* `"iterations"`: The maximum number of iterations of the BGF decoder.
* `"weak_key_threshold"`: The weak key threshold parameter, denoted `T` in the paper and set to 3.
* `"keys"`: The number of random error vectors for which decoding was attempted for each value of `l`.
* `"dfrs"`: A two-dimensional array encoding the DFRs for various values of `l` (the first index) and `delta` (the second index). A triple of numbers is recorded for each `l` and `delta`: The first number is the number of decoding failures where the syndromes were equal, the second is the total number of decoding failures, the third is the number of decoding attempts.
* `"decoding_failures"`: A list of the decoding failures. Each decoding failure is encoded as an associative array with entries `"l"` and `"delta"` being the corresponding parameters, `"H0"` and `"H1"` the supports of the first columns of the two circulant blocks, `"e_in"` the support of the error vector, and `"e_out"` the support of the output of the decoder. The number of decoding failures recorded in this way was capped at 1000.
* `"syndrome_matches"`: A list of any decoding failures where the decoded vector has the same syndrome as the error vector. (This was not observed in this dataset and thus is always empty, except for class `C`, where this data was not recorded as such vectors are common.)
* `"bad_vectors"`: A list of any vectors where a runtime error occurred during decoding. (This did not occur and thus is always empty.)
* `"runtime"`: The approximate total runtime in seconds used to compute the data.
