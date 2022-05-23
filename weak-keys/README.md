# Weak key data

This folder contains data on the frequency and DFR of weak keys in certain classes. The file `weak-key-frequencies.json` lists the observed frequency of weak keys for various values of `r` (the first index) and `T` (the second index). Weak key frequencies are stored as an ordered pair where the first entry is the number of weak keys observed and the second entry is the number of randomly generated keys tested.

The remaining files record data on the DFR for certain classes of weak keys, encoded in the following format:

* `"bits"`: The security parameter, denoted `λ` in the paper.
* `"r"`: The circulant block length.
* `"d"`: The column weight, always set to 15.
* `"t"`: The error vector weight, always set to 18.
* `"iterations"`: The maximum number of iterations of the BGF decoder.
* `"weak_key_threshold"`: The weak key threshold parameter, denoted `T` in the paper.
* `"weak_key_type"`: This is equal to 0, 1, 2, or 3. A value of 1, 2, or 3 indicates the weak keys are of type I, II, or III, respectively, as described in Vasseur's thesis. A value of 0 indicates that keys were generated at random until a weak key (of any of the three types) was found.
* `"keys"`: The number of random error vectors for which decoding was attempted.
* `"dfr"`: A triple encoding the DFR. The first number is the number of decoding failures where the syndromes were equal, the second is the total number of decoding failures, the third is the number of decoding attempts.
* `"decoding_failures"`: A list of the decoding failures. Each decoding failure is encoded as an associative array with entries `"H0"` and `"H1"` being the supports of the first columns of the two circulant blocks, `"e_in"` the support of the error vector, and `"e_out"` the support of the output of the decoder. The list is capped at 1000 entries.
* `"syndrome_matches"`: A list of any decoding failures where the decoded vector has the same syndrome as the error vector. (This was not observed in this dataset and thus is always empty.)
* `"bad_vectors"`: A list of any vectors where a runtime error occurred during decoding. (This did not occur and thus is always empty.)
* `"runtime"`: The approximate total runtime in seconds used to compute the data.
