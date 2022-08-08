# Decoder and data used in the paper "A Study of Error Floor Behavior in QC-MDPC Codes"

This repository contains raw data for our analysis of error floor behavior of QC-MDPC codes at the 20-bit security level. The files are organized as follows:

* `DFR`: data on decoding failure rates for random key-error vector pairs for various values of `r` (the circulant block size).
* `NCW-DFR-###`: data on decoding failure rates for random keys and error vectors chosen from certain classes of near-codeword. The number corresponds to the value of `r`.
* `weak-keys`: data on the frequency of certain classes of weak key and their effect on the DFR.
* `decoder.sage`: an implementation of the BIKE decoder. Parts of this code were written by Paolo Santini.

