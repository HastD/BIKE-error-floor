import itertools
import random

def bgf_decoder(r, d, t, s, H0tr, H1tr, Nb_iter, tau, thres="AFFINE", SC=1, threshold_cache=None):
    """
    This is a black-grey-flip decoder as specified in the shades of grey paper.
    In particular, note that there is only one black masked step and one grey masked step.

    Inputs:
    - r = Circulant block size (integer)
    - d = Column weight of parity check matrix (integer)
    - t = weight of original error vector (integer)
    - s = Syndrome to decode (vector over GF(2))
    - H0tr = Support of first column of first circulant block (list of integers)
    - H1tr = Support of first column of second circulant block (list of integers)
    - Nb_iter = Maximum number of iterations to run the BGF algorithm (integer)
    - tau = Amount by which the threshold is lowered for the grey masked flip (The paper and BIKE use tau = 3) (integer)
    - (optional) thres = Threshold function to call (string)
                    - (default) thres == "AFFINE": Predetermined BIKE threshold function for high security levels
                    - thres == "COUNTER": Threshold is the maximum counter
                    - thres == "EXACT": Threshold is computed using the function in the original BIKE specification
    - (optional) SC = BIKE security class. Default is 1. (integer)

    Output: the result e of running the decoder. (vector over GF(2))
    """
    iter_index = 0
    e = vector(GF(2), 2*r) # The error vector starts out as all zeroes
    s_copy = copy(s) # Create a copy of s so that we don't modify the original syndrome
    ws = s.hamming_weight() # Weight of s

    while ws > 0 and iter_index < Nb_iter:
        # Compute threshold
        if thres.upper() == "AFFINE":
            T = affine_threshold(ws, SC)
        elif thres.upper() == "COUNTER":
            T = maximum_counter(r, s_copy, H0tr, H1tr)
        elif thres.upper() == "EXACT":
            if threshold_cache is not None:
                if (r, d, ws, t) in threshold_cache:
                    T = threshold_cache[(r, d, ws, t)]
                else:
                    T = exact_threshold_ineq(r, d, ws, t)
                    threshold_cache[(r, d, ws, t)] = T
            else:
                T = exact_threshold(r, d, ws, t)
        else:
            raise ValueError("Invalid threshold function")

        # Bitflip decoder step
        (e, s_copy, black0, black1, gray0, gray1) = BFiter(r, s_copy, H0tr, H1tr, e, tau, T)

        if iter_index == 0: # Only on the first iteration, do black & grey masked flips
            T = (d+1)/2 # Shades of grey paper uses this different T for masked flips. Is this correct?
            (e, s_copy) = BFmaskedIter(r, s_copy, H0tr, H1tr, e, tau, T, black0, black1)
            (e, s_copy) = BFmaskedIter(r, s_copy, H0tr, H1tr, e, tau, T, gray0, gray1)

        # Update iter index and threshold weight
        ws = s_copy.hamming_weight()
        iter_index += 1
    return e

def bgf_decoder_failure_test(r, s, H0tr, H1tr, e_supp, e_out):
    """
    In this test, we compare the output e_out of the above decoder with the original error vector with support e_supp.
    The output is referred to as eprime since it may not equal the original error vector.

    Outputs a list with the following elements:
    - True if decoding succeeded, False if decoding failed (boolean)
    - True if decoding yields vector with syndrome s, False otherwise
    - A string that explains the decoding failure (string)
    """
    e_vec = vector(GF(2), 2*r)
    for i in e_supp:
        e_vec[i] = 1 # vector form of the original error vector (since it's given as support)
    if syndrome(r, H0tr, H1tr, e_out.support()) == s:
        if e_vec == e_out:
            msg = "Decoding was successful."
            return True, True, msg
        else:
            msg = "e_out has syndrome s but e_out differs from original e."
            return False, True, msg # this appears to be really rare
    else:
        msg = "Decoding failed because eprime does not have syndrome s."
        return False, False, msg

################################################

def syndrome(r, H0tr, H1tr, e_supp):
    """
    Computes the syndrome associated to an error vector.

    Inputs:
    - r = Circulant block size (integer)
    - H0tr = Support of first column of first circulant block (list of integers)
    - H1tr = Support of first column of second circulant block (list of integers)
    - e_supp = Support of the error vector (list of integers)

    Output: Syndrome of e. (vector over GF(2))
    """
    s = vector(GF(2), r);
    for i in e_supp:
        if i < r:
            for j in H0tr:
                s[(i + j) % r] += 1
        else:
            for j in H1tr:
                s[(i + j) % r] += 1
    return s

##############################################

#### THRESHOLD FUNCTIONS ####

def affine_threshold(ws, SC):
    """
    Hard-coded BIKE affine threshold function.

    Inputs:
    - ws = Syndrome weight (integer)
    - SC = Security level (integer)

    Output: Threshold (integer)
    """
    if SC == 1:
        val1 = floor(0.0069722*ws + 13.53)
        val2 = 36
    elif SC == 3:
        val1 = floor(0.005265*ws + 15.2588)
        val2 = 52
    elif SC == 5:
        val1 = floor(0.00402312*ws + 17.8785)
        val2 = 69
    else:
        raise ValueError("Invalid security level")
    thresh = max(val1, val2)
    return thresh

def maximum_counter(r, s, H0tr, H1tr):
    """
    Gives the maximum counter of a syndrome; that is, the largest Schur product of a column of the parity check matrix
    and the syndrome.

    Inputs:
    - r = Circulant block size (integer)
    - s = Syndrome (vector over GF(2))
    - H0tr = Support of first column of first circulant block (list of integers)
    - H1tr = Support of first column of second circulant block (list of integers)

    Output: Threshold (integer)
    """
    max_upc = 0
    for i in range(r):
        upc0 = 0 #After the next for loop, upc is the Hamming weight of of the schur product of column i of H0 and s (i.e. the number of ones the vectors share) (upc = unsatisfied parity check)
        upc1 = 0
        for j in H0tr:
            # (i + j) % r represents moving to the next column of circulant block 0
            upc0 += s[(i + j) % r]
        for j in H1tr:
            upc1 += s[(i + j) % r]
        #print(str(upc0)+ " " + str(upc1))
        max_upc = max(max_upc, upc0, upc1)
    #print(str(max_upc))
    return max_upc

def exact_threshold(r, d, ws, t):
    """
    Computes the exact threshold as given in the 2017 BIKE specification.

    Inputs:
    - r = Circulant block size (integer)
    - d = Column weight of parity check matrix (integer)
    - ws = Syndrome weight (integer)
    - t = Weight of the original error vector (integer)

    Output: Threshold (integer)
    """
    n = 2 * r
    w = 2 * d
    n_minus_w = n - w
    l_range = range(3, min(t, w), 2)
    X = RR(r) * sum(RR(l - 1) * binomial(w, l) * binomial(n_minus_w, t - l) for l in l_range) / binomial(n, t)
    pi1 = (ws + X) / (t * d)
    pi0 = (w * ws - X) / ((n - t) * d)

    c = log((1 - pi0) / (1 - pi1))
    T_num = log(RR((n - t) / t)) + d * c
    T_den = log(pi1 / pi0) + c
    T = ceil(T_num / T_den)
    return T

def exact_threshold_ineq(r, d, ws, t):
    """
    Same as above but uses inequality instead of equation (1)
    See BIKE 2017 page 18 : https://bikesuite.org/files/BIKE.2017.11.30.pdf

    Comments:
    - No complex logs (cf exact_threshold)
    - Small discrepancy with precomputed_data.sage (double check)
    - exact_threshold is based on Julia Chaulet's thesis (in French)
    """
    n = 2 * r
    w = 2 * d
    n_minus_w = n - w
    n_minus_t = n - t
    l_range = range(3, min(t, w), 2)
    X = r * sum((l - 1) * binomial(w, l) * binomial(n_minus_w, t - l) for l in l_range) / binomial(n, t)
    pi1 = (ws + X) / (t * d)
    pi0 = (w * ws - X) / ((n - t) * d)

    T = 0
    while t * binomial(d, T) * pi1^T * (1 - pi1)^(d - T) < n_minus_t * binomial(d, T) * pi0^T * (1 - pi0)^(d - T):
        T += 1
    return T





#########################################################################

#### BITFLIP FUNCTIONS ####

def BFiter(r, s, H0tr, H1tr, e, tau, T):
    """
    This function performs one step of bitflip decoding.

    Inputs:
    - r = Circulant block size (integer)
    - s = Syndrome to decode (vector over GF(2))
    - e = Error vector (vector over GF(2))
    - H0tr = Support of first column of first circulant block (list of integers)
    - H1tr = Support of first column of second circulant block (list of integers)
    - tau = Amount by which the threshold is lowered for the grey masked flip (integer)
    - T = Threshold function to call (string)
                    - thres == "AFFINE": Predetermined BIKE threshold function for high security levels
                    - thres == "COUNTER": Threshold is the maximum counter
                    - thres == "EXACT": Threshold is computed as determined in the original BIKE specification
    - (optional) SC = BIKE security class (integer)

    Output: A list containing:
    - e = Updated error vector (vector over GF(2))
    - s = Updated syndrome (vector over GF(2))
    - black0 = Black list for first circulant block (list of integers)
    - black1 = Black list for second circulant block (list of integers)
    - grey0 = Grey list for first circulant block (list of integers)
    - grey1 = Grey list for second circulant block (list of integers)
    """
    black0 = []
    black1 = []
    gray0 = []
    gray1 = []

    #Flips for H0 (first circulant block)
    copy_s = s.change_ring(ZZ)
    for i in range(r):
        upc = 0 #After the next for loop, upc is the Hamming weight of of the schur product of column i of H0 and s (i.e. the number of ones the vectors share) (upc = unsatisfied parity check)
        for j in H0tr:
            new_pos = (i+j) % r #(i+j)%r represents moving to the next column of circulant block 0
            upc += copy_s[new_pos]
        if upc >= T: # T is the black error threshold
            e[i] += 1 # flip the bit in e
            black0.append(i) #record the flipped position in black0
            """
            The for loop below does the "s = H(c^T + e^T)" step.
            Since e starts as the zero vector and c never changes throughout the decoding process, the result of H(c^T+e^T)
            amounts to adding column i of H to the current s every time the ith bit of e is flipped. The most important thing
            to remember about this (since I've forgotten it like 800 times now) is that, when starting this routine,
            s is H(c^T+e^T) with the starting e that is passed to the function, NOT e = (0, 0, ... , 0).
            """
            for j in H0tr:
                new_pos = (i+j) % r
                s[new_pos] += 1
        else:
            if upc >= (T - tau): #T-tau is the grey error threshold
                gray0.append(i)

    # Flips for H1 (second circulant block). Works the same as previous part except it uses H1 instead of H0
    for i in range(r):
        upc = 0
        for j in H1tr:
            new_pos = (i+j) % r
            upc += copy_s[new_pos]
        if upc >= T:
            e[r+i] += 1 #the extra r here puts us over in the second block
            black1.append(i)
            for j in H1tr:
                new_pos = (i+j) % r
                s[new_pos] += 1
        else:
            if upc >= T - tau:
                gray1.append(i)

    return e, s, black0, black1, gray0, gray1


def BFmaskedIter(r, s, H0tr, H1tr, e, tau, T, set0, set1):
    """
    This function performs one step of masked bitflip decoding.

    Inputs:
    - r = Circulant block size (integer)
    - s = Syndrome to decode (vector over GF(2))
    - H0tr = Support of first column of first circulant block (list of integers)
    - H1tr = Support of first column of second circulant block (list of integers)
    - e = Error vector (vector over GF(2))
    - tau = Amount by which the threshold is lowered for the grey masked flip (integer)
    - T = Threshold (integer)
    - set0 = Masked locations for first circulant block (list of integers)
    - set1 = Masked locations for second circulant block (list of integers)

    Output: A list containing:
    - e = Updated error vector (vector over GF(2))
    - s = Updated syndrome (vector over GF(2))
    """
    copy_s = s.change_ring(ZZ)
    #Flips for H0 (first circulant block)
    for i in set0:
        upc = 0
        for j in H0tr:
            new_pos = (i+j) % r
            upc += copy_s[new_pos]
        if upc >= T:
            e[i] += 1
            for j in H0tr:
                new_pos = (i+j) % r
                s[new_pos] += 1

    #Flips for H1 (second circulant block)
    for i in set1:
        upc = 0
        for j in H1tr:
            new_pos = (i+j) % r
            upc += copy_s[new_pos]
        if upc >= T:
            e[r+i] += 1
            for j in H1tr:
                new_pos = (i+j) % r
                s[new_pos] += 1
    return e, s

######################################################################

def sample_C(r, H0tr, H1tr):
    return H1tr + [r + j for j in H0tr]

def sample_N(r, H0tr, H1tr):
    if random.randrange(2):
        return H0tr
    else:
        return [r + j for j in H1tr]

def sample_2N(r, l, H0tr, H1tr):
    shift = random.randrange(r)
    shared = -1
    while shared < l: # since the weight of sum is variable, we reselect until we're sure we can sample l values
        supp1 = sample_N(r, H0tr, H1tr)
        supp2 = sample_N(r, H0tr, H1tr)
         # Adding 2 vectors together, given their supports, is done with XOR (i.e. symmetric difference)
        sumN = list(set(supp1).symmetric_difference(blockwise_shift(r, supp2, shift)))
        shared = len(sumN)
    sumN.sort()
    return sumN # Also return shared, because: \delta = shared + t - 2ell. need to alter below in this case

def element_of_AtlS(r, t, l, H0tr, H1tr, S="N"):
    """
    This function gives a random element of the set A_{t,l}(S).
    It's the same as Algorithm 16.1 from the thesis.

    Inputs:
    - r = Circulant block size (integer)
    - t = Error vector size (integer)
    - l = Number of ones that the output should share with an element of S (integer)
          (Note: l can't exceed t, and for N can't exceed d.)
    - H0tr = Support of first column of first circulant block (list of integers)
    - H1tr = Support of first column of second circulant block (list of integers)
    - (optional) S = Set of error-prone vectors to consider. Default is "N".
                    - S == "N": (d,d) near-codewords
                    - S == "2N": Sum of two (d,d) near-codewords
                    - S == "C": Low-weight codewords

    Output: The support of an element of A_{t,l}(S) (list of integers)
    """
    if S.upper() == "C": # output has weight w = 2d
        c = sample_C(r, H0tr, H1tr)
    elif S.upper() == "N": # output has weight d
        c = sample_N(r, H0tr, H1tr)
    elif S.upper() == "2N": # output has weight a bit smaller than w
        c = sample_2N(r, l, H0tr, H1tr)
    else:
        raise ValueError("Invalid name for S")
    shift = random.randrange(r)
    list1 = random.sample(c, l)
    list2 = random.sample([x for x in range(2*r) if x not in c], t - l) #2*r is corrected from thesis where r is displayed
    samples = list1 + list2
    return (blockwise_shift(r, samples, shift),len(c))

def blockwise_shift(r, v_supp, shift):
    """
    Performs a blockwise circular shift of a vector with 2*r entries by s.
    Note that the input is the support of v, not v itself.

    Inputs:
    - r = Circulant block size (integer)
    - v_supp = Support of v (list of integers)
    - shift = Amount by which to shift (integer)

    Output: The support of the vector v circular shifted by shift (list of integers)
    """
    v = []
    for i in v_supp:
        if i < r:
            v.append((i + shift) % r)
        else:
            v.append((i + shift) % r + r)
    v.sort()
    return v

def sparse_vector(r, d):
    """
    Generates the support of a random vector with length r and weight d.

    Inputs:
    - r = Vector length (integer)
    - d = Vector weight (integer)

    Output: The support of a random vector with length r and weight d (list of integers)
    """
    supp = random.sample(range(r), d)
    supp.sort()
    return supp

def dist_supp(r, i, j):
    """
    Distance between two entries of support of column vector (from Ch 15)
    """
    return min((j - i) % r, (i - j) % r)

def is_weak_key(r, d, T, H0supp, H1supp):
    """
    Rejection algorithm for weak keys of Type I, II, III (Algorithm 15.3 p155 of thesis)

    Input:
    - r = size of circulant block
    - d = column weight
    - T = threshold used for rejection algorithm (different from the threshold in BGF)
    - H0supp, H1supp = support of first columns of circulant blocks

    Output:
    - Boolean output (True = weak key)
    """
    lst = (H0supp, H1supp)
    for i in range(2):
        Hi_supp = lst[i]
        S = [0]*int(r/2 + 1)
        # Iterate over ordered pairs of distinct elements of Hi_supp
        for j_k, j_l in itertools.combinations(Hi_supp, 2):
            # This is equivalent to `delta = dist_supp(r, i, j)`,
            # taking advantage of the fact that 0 <= i < j < r.
            diff = j_l - j_k
            delta = min(diff, r - diff)
            S[delta] += 1
            if S[delta] >= T:
                # Type I or II weak key
                return True
    H0s = set(H0supp)
    # Iterate over k, l in range(d)
    for k, l in itertools.product(range(d), repeat=2):
        shift = H0supp[k] - H1supp[l]
        H1s = {(i + shift) % r for i in H1supp}
        if len(H0s & H1s) >= T:
            # Type III weak key
            return True
    return False

def col_to_row(r, h):
    """
    Take a column generator (defined by weights) of a circulant matrix and returns a row generator.

    Input:
    - r = size of circulant block (integer)
    - h =  support of first column of a circulant matrix (list of integers)

    Output:
    - support of first row of a circulant matrix (list of integers)
    """
    l = [(-i) % r for i in h]
    l.sort()
    return l

def non_weak_key(r, d, T):
    """Generate a random PCH, but make sure it's not weak."""
    weak = True
    while weak:
        H0tr_supp = sparse_vector(r, d)
        H1tr_supp = sparse_vector(r, d)
        weak = is_weak_key(r, d, T, H0tr_supp, H1tr_supp)
    return H0tr_supp, H1tr_supp

def weak_key(r, d, T, weak_type=None):
    """Generate a random PCH, but make sure it is weak."""
    if weak_type is None or weak_type == 0:
        weak = False
        while not weak:
            H0tr_supp = sparse_vector(r, d)
            H1tr_supp = sparse_vector(r, d)
            weak = is_weak_key(r, d, T, H0tr_supp, H1tr_supp)
        return H0tr_supp, H1tr_supp
    elif weak_type == 1:
        return weak_key_type1(r, d, T)
    elif weak_type == 2:
        return weak_key_type2(r, d, T)
    elif weak_type == 3:
        return weak_key_type3(r, d, T)
    else:
        raise ValueError("Invalid value of weak_type parameter.")

def weak_key_type1(r, d, T):
    """Generate a random weak key of type I (Algorithm 15.1 in Vasseur's thesis)."""
    f = T + 1
    p = list(range(f)) + random.sample(range(f, r), d - f)
    delta = random.randrange(1, floor(r/2) + 1)
    l = random.randrange(r)
    h_weak = [(delta * (l + pk)) % r for pk in p]
    h_weak.sort()
    h_random = sparse_vector(r, d)
    if random.randrange(2) == 0:
        return h_weak, h_random
    else:
        return h_random, h_weak

def weak_key_type2(r, d, T):
    """Generate a random weak key of type II (Algorithm 15.2 in Vasseur's thesis)."""
    m = T + 1
    s = d - m
    a = [0] + random.sample(range(1, d), s - 1) + [d]
    a.sort()
    b = [0] + random.sample(range(1, r - d), s - 1) + [r - d]
    b.sort()
    o = [a[i+1] - a[i] for i in range(s)]
    z = [b[i+1] - b[i] for i in range(s)]
    delta = random.randrange(1, floor(r/2) + 1)
    l = random.randrange(z[1] + o[1])
    h_weak = []
    i = -l
    for j in range(s):
        i += z[j]
        for k in range(o[j]):
            h_weak.append((delta * (i + k)) % r)
        i += o[j]
    h_weak.sort()
    h_random = sparse_vector(r, d)
    if random.randrange(2) == 0:
        return h_weak, h_random
    else:
        return h_random, h_weak

def weak_key_type3(r, d, T):
    """Generate a random weak key of type III."""
    h_random = sparse_vector(r, d)
    l = random.randrange(r)
    overlap = [(j + l) % r for j in random.sample(h_random, T)]
    h_weak = overlap + random.sample(set(range(r)) - set(overlap), d - T)
    h_weak.sort()
    if random.randrange(2) == 0:
        return h_weak, h_random
    else:
        return h_random, h_weak

