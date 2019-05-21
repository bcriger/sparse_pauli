from functools import reduce #PY3 compatibility
from itertools import chain, combinations, product
from operator import add, xor

__all__ = [
            "Pauli", "I", "X", "Y", "Z", "local_group",
            "generated_group", "char", "str_pauli"
        ]

PHASES_README = """
For multiplying Paulis by complex numbers, we have a dict that
basically gives us the complex log.
"""
PHASES = {1 : 0, 1j : 1, -1 : 2, -1j : 3}

def char(elem, pauli):
    """
    Spits out an X/Y/Z depending on what sets elem is a member of.
    """
    
    if elem in pauli.x_set:
        return 'Y' if elem in pauli.z_set else 'X'
    
    return 'Z' if elem in pauli.z_set else 'I'


class Pauli(object):
    """
    Stores the X and Z support of a Pauli as two sets, with the phase 
    as an integer between 0 and 3.

    Supports multiplication, commutation, Hadamard and CNOT, and a few
    other operations.

    Implied internal representation is
    (-i)^ph * (X_{x_set}) * (Z_{z_set}).
    """
    def __init__(self, x_set={}, z_set={}, ph=0):
        self.x_set = set(x_set)
        self.z_set = set(z_set)
        self.ph = ph % 4 

    # basic info/operations
    def copy(self):
        return Pauli(self.x_set, self.z_set)

    def support(self):
        """
        The set of all qubits on which the Pauli is non-trivial.
        """
        bit_set = self.x_set.union(self.z_set)
        try: 
            return list(sorted(bit_set))
        except TypeError:
            return list(bit_set)

    def weight(self):
        return len(self.support())

    def xz_pair(self):
        """
        For convenience, returns an X Pauli and a Z Pauli whose 
        product is `self`.
        """
        return Pauli(self.x_set, {}), Pauli({}, self.z_set)

    def str_sprt_pair(self):
        """
        Returns a tuple (str, sprt), where str is a string, and sprt
        is a sorted (if possible) tuple of qubit labels on which the
        Pauli exists.
        """
        sprt = tuple(self.support())
        string = reduce(add, [char(elem, self) for elem in sprt], '')
        return string, sprt
    
    # printing/comparison
    def __repr__(self):
        total_phase = (self.ph - len(self.x_set & self.z_set)) % 4 
        string = {0 : '', 1 : 'i', 2 : '-', 3 : '-i'}[total_phase]
        
        sprt = self.support()

        if len(sprt) == 0:
            return 'I'
        
        for elem in sprt:
            ch = char(elem, self)
            string += '{}[{}] '.format(ch, elem)
        
        #strip trailing space
        return string[:-1]


    def __eq__(self, othr, sign=False):
        strings_eq = (self.x_set == othr.x_set) & (self.z_set == othr.z_set)
        if sign:
            return strings_eq & (self.ph == othr.ph)
        else:
            return strings_eq
    
    def __ne__(self, othr, sign=False):
        strings_neq = (self.x_set != othr.x_set) | (self.z_set != othr.z_set)
        if strings_neq or not(sign):
            return strings_neq
        else:
            return strings_neq | (self.ph != othr.ph)

    def __hash__(self):
        return hash((tuple(self.x_set), tuple(self.z_set)))
    
    #actual math
    def com(self, othr):
        """
        This method suffers from the usual ambiguity,
        com(self, othr) == 1 means they ANTIcommute.
        """
        return (len(self.x_set & othr.z_set) +
                    len(self.z_set & othr.x_set)) % 2

    def __mul__(self, othr):
        if type(othr) == Pauli:
            com_term = len(self.z_set & othr.x_set) % 2
            new_ph = self.ph + othr.ph + 2 * com_term
            return Pauli(self.x_set ^ othr.x_set,
                            self.z_set ^ othr.z_set, new_ph)
        else:
            return self.__rmul__(othr)

    def __rmul__(self, othr):
        # try number
        try:
            new_ph = self.ph + PHASES[othr]
            return Pauli(self.x_set, self.z_set, new_ph)
        except KeyError:
            raise ValueError("Paulis can only be multiplied by "
                "{}, {} entered.".format(list(PHASES.keys()), othr))

    def __neg__(self):
        return -1 * self

    def __call__(self,othr):
        return -othr if self.com(othr) else othr

    def cnot(self, ctrl_targs):
        """
        acts a cnot on pairs of qubits given by the set of tuples 
        ctrl_targs.
        """
        for ctrl, targ in ctrl_targs:
            if ctrl in self.x_set:
                self.x_set ^= set([targ])
            if targ in self.z_set:
                self.z_set ^= set([ctrl])
        pass
    
    def cz(self, prs):
        """
        acts a cz on pairs of qubits given by the set of tuples 
        prs.
        """
        for q_0, q_1 in prs:
            if q_0 in self.x_set:
                self.z_set ^= set([q_1])
            if q_1 in self.x_set:
                self.z_set ^= set([q_0])
        pass
    
    def h(self, qs):
        """
        acts a Hadamard on each bit in qs.
        """
        switches = (self.x_set ^ self.z_set) & set(qs)
        self.x_set ^= switches
        self.z_set ^= switches
        self.ph += 2 * len(self.x_set & self.z_set & set(qs))
        self.ph %= 4
        pass

    def p(self, qs):
        switches = self.x_set & set(qs)
        self.z_set ^= switches
        self.ph = (self.ph + len(switches)) % 4
        pass

    def prep(self, qs):
        """
        Flawlessly prepares qubits in the list `qs` in the +1 
        eigenstate a basis which is not tracked by this object.
        The sense in which 'preparation' is meant here is in the sense
        "providing an error free qubit at a specific location". 
        To do arbitrary GK preparations, I'd have to simulate what 
        happens to a tableau, which is outside the scope of this 
        library.
        """
        self.x_set -= set(qs)
        self.z_set -= set(qs)
        pass

    def meas(self, qs, basis='Z'):
        """
        Flawlessly 'measures' qubits in the given basis. I'm not 
        modeling stabiliser states with these Paulis, I'm modeling 
        errors. This means that the measurement result is 1 iff the 
        Pauli being measured anticommutes with the Pauli I'm measuring.  
        """
        _basis_check(basis)
    
        anticom_set = self.x_set if basis == 'Z' else self.z_set
    
        return {0 :  set(qs) - anticom_set, 1 : set(qs) & anticom_set}
    

#---------------------------public functions--------------------------#
I = Pauli()

X = lambda sett: Pauli(x_set=set(sett))

Y = lambda sett: Pauli(x_set=set(sett), z_set=set(sett), ph=len(sett))

Z = lambda sett: Pauli(z_set=set(sett))

def local_group(support):
    """
    Generator iterating over the entire Pauli group on bits in the 
    iterable `support`.
    """
    for set_pr in product(_powerset(support), repeat=2):
        yield Pauli(*set_pr)

def generated_group(x_sets, z_sets):
    """
    To generate a restricted group directly, rather than filtering
    a local group, we can start from lists of sets of acceptable Xs and
    Zs.
    """
    for x_lst, z_lst in product(_powerset(x_sets), _powerset(z_sets)):
        x_set = reduce(xor, x_lst, set())
        z_set = reduce(xor, z_lst, set())
        yield Pauli(x_set, z_set)

def str_pauli(string, support=None, error_check=True):
    
    string = string.upper()
    
    if error_check:
        if not all([ltr in 'IXYZ' for ltr in string]):
            raise InputError("All letters in input string must be "
                "I, X, Y, or Z. {} entered.".format(string) + 
                "Note: Phases not supported.")

    if support is None:
        support = range(len(string))

    new_pauli = I

    for idx, ltr in enumerate(string):
    
        if ltr in 'XY':
            new_pauli.x_set |= {support[idx]}
        
        if ltr in 'YZ': # *not* an elif
            new_pauli.z_set |= {support[idx]}
    
    return new_pauli

#---------------------------------------------------------------------#

#-------------------------private functions---------------------------#

def _basis_check(basis):
    if basis not in ['X', 'Z']:
        raise ValueError("basis must be 'X' or 'Z', "
                            "{} entered.".format(basis))

def _powerset(iterable):
    """
    From the itertools recipe page:
    powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)
    """
    s = list(iterable)
    return chain.from_iterable(combinations(s, r)
                                for r in range(len(s)+1))

#---------------------------------------------------------------------#
