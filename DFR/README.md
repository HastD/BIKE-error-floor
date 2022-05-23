# DFR data for random error vectors

This folder contains data on decoding failure rates for randomly generated key-error vector pairs (conditioned on the key not being a weak key). The JSON files contain data in the following format:

* `"bits"`: The security parameter, denoted `λ` in the paper.
* `"r"`: The circulant block length.
* `"d"`: The column weight, always set to 15.
* `"t"`: The error vector weight, always set to 18.
* `"iterations"`: The maximum number of iterations of the BGF decoder.
* `"weak_key_threshold"`: The weak key threshold parameter, denoted `T` in the paper. This is set to `3` except in the file ending in `-T4.json`.
* `"keys"`: The number of random error vectors for which decoding was attempted.
* `"dfr"`: A triple encoding the DFR. The first number is the number of decoding failures where the syndromes were equal, the second is the total number of decoding failures, the third is the number of decoding attempts.
* `"decoding_failures"`: A list of the decoding failures. Each decoding failure is encoded as an associative array with entries `"H0"` and `"H1"` being the supports of the first columns of the two circulant blocks, `"e_in"` the support of the error vector, and `"e_out"` the support of the output of the decoder.
* `"syndrome_matches"`: A list of any decoding failures where the decoded vector has the same syndrome as the error vector. (This was not observed in this dataset and thus is always empty.)
* `"bad_vectors"`: A list of any vectors where a runtime error occurred during decoding. (This did not occur and thus is always empty.)
* `"runtime"`: The approximate total runtime in seconds used to compute the data.
