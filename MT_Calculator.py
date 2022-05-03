##
# imports
import math
import warnings

from Bio import SeqUtils, Seq
from Bio import BiopythonWarning

# from Bio import BiopythonDeprecationWarning

##
# Thermodynamic lookup tables
# SantaLucia & Hicks (2004), Annu. Rev. Biophys. Biomol. Struct 33: 415-440
nn_table = {
    "init": (0.2, -5.7), "init_A/T": (2.2, 6.9), "init_G/C": (0, 0),
    "init_oneG/C": (0, 0), "init_allA/T": (0, 0), "init_5T/A": (0, 0),
    "sym": (0, -1.4),
    "AA/TT": (-7.6, -21.3), "AT/TA": (-7.2, -20.4), "TA/AT": (-7.2, -20.4),
    "CA/GT": (-8.5, -22.7), "GT/CA": (-8.4, -22.4), "CT/GA": (-7.8, -21.0),
    "GA/CT": (-8.2, -22.2), "CG/GC": (-10.6, -27.2), "GC/CG": (-9.8, -24.4),
    "GG/CC": (-8.0, -19.0)}

# Internal mismatch and inosine table
imm_table = {
    "AG/TT": (1.0, 0.9), "AT/TG": (-2.5, -8.3), "CG/GT": (-4.1, -11.7),
    "CT/GG": (-2.8, -8.0), "GG/CT": (3.3, 10.4), "GG/TT": (5.8, 16.3),
    "GT/CG": (-4.4, -12.3), "GT/TG": (4.1, 9.5), "TG/AT": (-0.1, -1.7),
    "TG/GT": (-1.4, -6.2), "TT/AG": (-1.3, -5.3), "AA/TG": (-0.6, -2.3),
    "AG/TA": (-0.7, -2.3), "CA/GG": (-0.7, -2.3), "CG/GA": (-4.0, -13.2),
    "GA/CG": (-0.6, -1.0), "GG/CA": (0.5, 3.2), "TA/AG": (0.7, 0.7),
    "TG/AA": (3.0, 7.4),
    "AC/TT": (0.7, 0.2), "AT/TC": (-1.2, -6.2), "CC/GT": (-0.8, -4.5),
    "CT/GC": (-1.5, -6.1), "GC/CT": (2.3, 5.4), "GT/CC": (5.2, 13.5),
    "TC/AT": (1.2, 0.7), "TT/AC": (1.0, 0.7),
    "AA/TC": (2.3, 4.6), "AC/TA": (5.3, 14.6), "CA/GC": (1.9, 3.7),
    "CC/GA": (0.6, -0.6), "GA/CC": (5.2, 14.2), "GC/CA": (-0.7, -3.8),
    "TA/AC": (3.4, 8.0), "TC/AA": (7.6, 20.2),
    "AA/TA": (1.2, 1.7), "CA/GA": (-0.9, -4.2), "GA/CA": (-2.9, -9.8),
    "TA/AA": (4.7, 12.9), "AC/TC": (0.0, -4.4), "CC/GC": (-1.5, -7.2),
    "GC/CC": (3.6, 8.9), "TC/AC": (6.1, 16.4), "AG/TG": (-3.1, -9.5),
    "CG/GG": (-4.9, -15.3), "GG/CG": (-6.0, -15.8), "TG/AG": (1.6, 3.6),
    "AT/TT": (-2.7, -10.8), "CT/GT": (-5.0, -15.8), "GT/CT": (-2.2, -8.4),
    "TT/AT": (0.2, -1.5),
    "AI/TC": (-8.9, -25.5), "TI/AC": (-5.9, -17.4), "AC/TI": (-8.8, -25.4),
    "TC/AI": (-4.9, -13.9), "CI/GC": (-5.4, -13.7), "GI/CC": (-6.8, -19.1),
    "CC/GI": (-8.3, -23.8), "GC/CI": (-5.0, -12.6),
    "AI/TA": (-8.3, -25.0), "TI/AA": (-3.4, -11.2), "AA/TI": (-0.7, -2.6),
    "TA/AI": (-1.3, -4.6), "CI/GA": (2.6, 8.9), "GI/CA": (-7.8, -21.1),
    "CA/GI": (-7.0, -20.0), "GA/CI": (-7.6, -20.2),
    "AI/TT": (0.49, -0.7), "TI/AT": (-6.5, -22.0), "AT/TI": (-5.6, -18.7),
    "TT/AI": (-0.8, -4.3), "CI/GT": (-1.0, -2.4), "GI/CT": (-3.5, -10.6),
    "CT/GI": (0.1, -1.0), "GT/CI": (-4.3, -12.1),
    "AI/TG": (-4.9, -15.8), "TI/AG": (-1.9, -8.5), "AG/TI": (0.1, -1.8),
    "TG/AI": (1.0, 1.0), "CI/GG": (7.1, 21.3), "GI/CG": (-1.1, -3.2),
    "CG/GI": (5.8, 16.9), "GG/CI": (-7.6, -22.0),
    "AI/TI": (-3.3, -11.9), "TI/AI": (0.1, -2.3), "CI/GI": (1.3, 3.0),
    "GI/CI": (-0.5, -1.3)}

# Terminal Mismatch
tmm_table = {
    "AA/TA": (-3.1, -7.8), "TA/AA": (-2.5, -6.3), "CA/GA": (-4.3, -10.7),
    "GA/CA": (-8.0, -22.5),
    "AC/TC": (-0.1, 0.5), "TC/AC": (-0.7, -1.3), "CC/GC": (-2.1, -5.1),
    "GC/CC": (-3.9, -10.6),
    "AG/TG": (-1.1, -2.1), "TG/AG": (-1.1, -2.7), "CG/GG": (-3.8, -9.5),
    "GG/CG": (-0.7, -19.2),
    "AT/TT": (-2.4, -6.5), "TT/AT": (-3.2, -8.9), "CT/GT": (-6.1, -16.9),
    "GT/CT": (-7.4, -21.2),
    "AA/TC": (-1.6, -4.0), "AC/TA": (-1.8, -3.8), "CA/GC": (-2.6, -5.9),
    "CC/GA": (-2.7, -6.0), "GA/CC": (-5.0, -13.8), "GC/CA": (-3.2, -7.1),
    "TA/AC": (-2.3, -5.9), "TC/AA": (-2.7, -7.0),
    "AC/TT": (-0.9, -1.7), "AT/TC": (-2.3, -6.3), "CC/GT": (-3.2, -8.0),
    "CT/GC": (-3.9, -10.6), "GC/CT": (-4.9, -13.5), "GT/CC": (-3.0, -7.8),
    "TC/AT": (-2.5, -6.3), "TT/AC": (-0.7, -1.2),
    "AA/TG": (-1.9, -4.4), "AG/TA": (-2.5, -5.9), "CA/GG": (-3.9, -9.6),
    "CG/GA": (-6.0, -15.5), "GA/CG": (-4.3, -11.1), "GG/CA": (-4.6, -11.4),
    "TA/AG": (-2.0, -4.7), "TG/AA": (-2.4, -5.8),
    "AG/TT": (-3.2, -8.7), "AT/TG": (-3.5, -9.4), "CG/GT": (-3.8, -9.0),
    "CT/GG": (-6.6, -18.7), "GG/CT": (-5.7, -15.9), "GT/CG": (-5.9, -16.1),
    "TG/AT": (-3.9, -10.5), "TT/AG": (-3.6, -9.8)}

# # Dangling ends table
# # Bommarito et al. (2000), Nucl Acids Res 28: 1929-1934
de_table = {
    "AA/.T": (0.2, 2.3), "AC/.G": (-6.3, -17.1), "AG/.C": (-3.7, -10.0),
    "AT/.A": (-2.9, -7.6), "CA/.T": (0.6, 3.3), "CC/.G": (-4.4, -12.6),
    "CG/.C": (-4.0, -11.9), "CT/.A": (-4.1, -13.0), "GA/.T": (-1.1, -1.6),
    "GC/.G": (-5.1, -14.0), "GG/.C": (-3.9, -10.9), "GT/.A": (-4.2, -15.0),
    "TA/.T": (-6.9, -20.0), "TC/.G": (-4.0, -10.9), "TG/.C": (-4.9, -13.8),
    "TT/.A": (-0.2, -0.5),
    ".A/AT": (-0.7, -0.8), ".C/AG": (-2.1, -3.9), ".G/AC": (-5.9, -16.5),
    ".T/AA": (-0.5, -1.1), ".A/CT": (4.4, 14.9), ".C/CG": (-0.2, -0.1),
    ".G/CC": (-2.6, -7.4), ".T/CA": (4.7, 14.2), ".A/GT": (-1.6, -3.6),
    ".C/GG": (-3.9, -11.2), ".G/GC": (-3.2, -10.4), ".T/GA": (-4.1, -13.1),
    ".A/TT": (2.9, 10.4), ".C/TG": (-4.4, -13.1), ".G/TC": (-5.2, -15.0),
    ".T/TA": (-3.8, -12.6)}


#

##
# a function for checking the input
def _check(seq):
    baseset = ("A", "C", "G", "T", "I")
    seq = "".join([base for base in seq if base in baseset])
    return seq


##
def _key_error(neighbors, strict):
    """Throw an error or a warning if there is no data for the neighbors (PRIVATE)."""
    # We haven't found the key in the tables
    if strict:
        raise ValueError("no thermodynamic data for neighbors %r available" % neighbors)
    else:
        warnings.warn(
            "no themodynamic data for neighbors %r available. "
            "Calculation will be wrong" % neighbors,
            BiopythonWarning,
        )


##
# Salt correction: If any of K, Tris, Mg and dNTPS is non-zero, a 'sodium-equivalent' concentration is calculated
# according to von Ahsen et al. (2001, Clin Chem 47: 1956-1961): [Na_eq] = [Na+] + [K+] + [Tris]/2 + 120*([Mg2+] - [
# dNTPs])^0.5
# Correction for deltaS: 0.368 x (N-1) x ln[Na+]
# (SantaLucia (1998), Proc Natl Acad Sci USA 95: 1460-1465)
def salt_correction(Na=0, K=0, Tris=0, Mg=0, dNTPs=0, seq=None):
    if seq:
        seq = str(seq)
    corr = 0
    Mon = Na + K + Tris / 2.0  # Note: all these values are millimolar
    mg = Mg * 1e-3  # Lowercase ions (mg, mon, dntps) are molar
    # Na equivalent according to von Ahsen et al. (2001):
    if sum((K, Mg, Tris, dNTPs)) > 0 and dNTPs < Mg:
        # dNTPs bind Mg2+ strongly. If [dNTPs] is larger or equal than
        # [Mg2+], free Mg2+ is considered not to be relevant.
        Mon += 120 * math.sqrt(Mg - dNTPs)
    mon = Mon * 1e-3
    # Note: math.log = ln(), math.log10 = log()
    if not mon:
        raise ValueError(
            "Total ion concentration of zero is not allowed in this method."
        )
    corr = 0.368 * (len(seq) - 1) * math.log(mon)
    return corr


##
# Calculate the Tm using NN method. DNA concentration in [nM], ions in [mM].
def Tm_NN(seq, check=True, strict=True, c_seq=None,
          shift=0, dnac1=10000, dnac2=10000, Na=0, K=50, Tris=1,
          Mg=1,
          dNTPs=0,
          ):
    seq = str(seq)
    if not c_seq:
        c_seq = Seq.Seq(seq).complement()
    c_seq = str(c_seq)
    if check:
        seq = _check(seq)
        c_seq = _check(c_seq)
    tmp_seq = seq
    tmp_cseq = c_seq
    delta_h = 0
    delta_s = 0
    d_h = 0  # Names for indexes
    d_s = 1  # 0 and 1

    # Dangling ends?
    if shift or len(seq) != len(c_seq):
        # Align both sequences using the shift parameter
        if shift > 0:
            tmp_seq = "." * shift + seq
        if shift < 0:
            tmp_cseq = "." * abs(shift) + c_seq
        if len(tmp_cseq) > len(tmp_seq):
            tmp_seq += (len(tmp_cseq) - len(tmp_seq)) * "."
        if len(tmp_cseq) < len(tmp_seq):
            tmp_cseq += (len(tmp_seq) - len(tmp_cseq)) * "."
        # Remove 'over-dangling' ends
        while tmp_seq.startswith("..") or tmp_cseq.startswith(".."):
            tmp_seq = tmp_seq[1:]
            tmp_cseq = tmp_cseq[1:]
        while tmp_seq.endswith("..") or tmp_cseq.endswith(".."):
            tmp_seq = tmp_seq[:-1]
            tmp_cseq = tmp_cseq[:-1]
        # Now for the dangling ends
        if tmp_seq.startswith(".") or tmp_cseq.startswith("."):
            left_de = tmp_seq[:2] + "/" + tmp_cseq[:2]
            try:
                delta_h += de_table[left_de][d_h]
                delta_s += de_table[left_de][d_s]
            except KeyError:
                _key_error(left_de, strict)
            tmp_seq = tmp_seq[1:]
            tmp_cseq = tmp_cseq[1:]
        if tmp_seq.endswith(".") or tmp_cseq.endswith("."):
            right_de = tmp_cseq[-2:][::-1] + "/" + tmp_seq[-2:][::-1]
            try:
                delta_h += de_table[right_de][d_h]
                delta_s += de_table[right_de][d_s]
            except KeyError:
                _key_error(right_de, strict)
            tmp_seq = tmp_seq[:-1]
            tmp_cseq = tmp_cseq[:-1]

    # Now for terminal mismatches
    left_tmm = tmp_cseq[:2][::-1] + "/" + tmp_seq[:2][::-1]
    if left_tmm in tmm_table:
        delta_h += tmm_table[left_tmm][d_h]
        delta_s += tmm_table[left_tmm][d_s]
        tmp_seq = tmp_seq[1:]
        tmp_cseq = tmp_cseq[1:]
    right_tmm = tmp_seq[-2:] + "/" + tmp_cseq[-2:]
    if right_tmm in tmm_table:
        delta_h += tmm_table[right_tmm][d_h]
        delta_s += tmm_table[right_tmm][d_s]
        tmp_seq = tmp_seq[:-1]
        tmp_cseq = tmp_cseq[:-1]

    # Now everything 'unusual' at the ends is handled and removed and we can
    # look at the initiation.
    # One or several of the following initiation types may apply:

    # Type: General initiation value
    delta_h += nn_table["init"][d_h]
    delta_s += nn_table["init"][d_s]

    # Type: Duplex with no (allA/T) or at least one (oneG/C) GC pair
    if SeqUtils.GC(seq) == 0:
        delta_h += nn_table["init_allA/T"][d_h]
        delta_s += nn_table["init_allA/T"][d_s]
    else:
        delta_h += nn_table["init_oneG/C"][d_h]
        delta_s += nn_table["init_oneG/C"][d_s]

    # Type: Penalty if 5' end is T
    if seq.startswith("T"):
        delta_h += nn_table["init_5T/A"][d_h]
        delta_s += nn_table["init_5T/A"][d_s]
    if seq.endswith("A"):
        delta_h += nn_table["init_5T/A"][d_h]
        delta_s += nn_table["init_5T/A"][d_s]

    # Type: Different values for G/C or A/T terminal basepairs
    ends = seq[0] + seq[-1]
    AT = ends.count("A") + ends.count("T")
    GC = ends.count("G") + ends.count("C")
    delta_h += nn_table["init_A/T"][d_h] * AT
    delta_s += nn_table["init_A/T"][d_s] * AT
    delta_h += nn_table["init_G/C"][d_h] * GC
    delta_s += nn_table["init_G/C"][d_s] * GC

    # Finally, the 'zipping'
    for basenumber in range(len(tmp_seq) - 1):
        neighbors = (
                tmp_seq[basenumber: basenumber + 2]
                + "/"
                + tmp_cseq[basenumber: basenumber + 2]
        )
        if neighbors in imm_table:
            delta_h += imm_table[neighbors][d_h]
            delta_s += imm_table[neighbors][d_s]
        elif neighbors[::-1] in imm_table:
            delta_h += imm_table[neighbors[::-1]][d_h]
            delta_s += imm_table[neighbors[::-1]][d_s]
        elif neighbors in nn_table:
            delta_h += nn_table[neighbors][d_h]
            delta_s += nn_table[neighbors][d_s]
        elif neighbors[::-1] in nn_table:
            delta_h += nn_table[neighbors[::-1]][d_h]
            delta_s += nn_table[neighbors[::-1]][d_s]
        else:
            # We haven't found the key...
            _key_error(neighbors, strict)

    k = (dnac1 - (dnac2 / 2.0)) * 1e-9
    R = 1.987  # universal gas constant in Cal/degrees C*Mol
    corr = salt_correction(Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, seq=seq)
    delta_s += corr
    melting_temp = (1000 * delta_h) / (delta_s + (R * (math.log(k)))) - 273.15

    return melting_temp


##
print(Tm_NN(seq='ACACCTTAATCACCGCTTCA', c_seq='TGTGGAATTAGTGGCGAAG'))
