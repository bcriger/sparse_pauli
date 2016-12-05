__all__ = ["Pauli"]

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

    #printing/comparison
    def __repr__(self):
        total_phase = (self.ph - len(self.x_set & self.z_set)) % 4 
        string = {0 : '', 1 : 'i', 2 : '-', 3 : '-i'}[total_phase]
        
        try:
            support = sorted(list(self.x_set | self.z_set))
        except TypeError:
            support = list(self.x_set | self.z_set)

        if len(support) == 0:
            return 'I'
        
        for elem in support:
            if elem in self.x_set:
                char = 'Y' if elem in self.z_set else 'X'
            else:
                char = 'Z'

            string += '{}[{}] '.format(char, elem)
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
    
    #actual math
    def com(self, othr):
        """
        This method suffers from the usual ambiguity,
        com(self, othr) == 1 means they ANTIcommute.
        """
        return (len(self.x_set & othr.z_set) +
                    len(self.z_set & othr.x_set)) % 2

    def __mul__(self, othr):
        new_ph = self.ph + othr.ph + 2 * self.com(othr)
        return Pauli(self.x_set ^ othr.x_set,
                        self.z_set ^ othr.z_set, new_ph)


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
        switches = self.x_set & qs
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
        basis_check(basis)
    
        anticom_set = self.x_set if basis == 'Z' else self.z_set
    
        return {0 :  set(qs) - anticom_set, 1 : set(qs) & anticom_set}
    
    def copy(self):
        return Pauli(self.x_set, self.z_set)

    def support(self):
        """
        The set of all qubits on which the Pauli is non-trivial.
        """
        return self.x_set.union(self.z_set)

    def weight(self):
        return len(self.support())

    def xz_pair(self):
        """
        For convenience, returns an X Pauli and a Z Pauli whose 
        product is `self`.
        """
        return Pauli(self.x_set, {}), Pauli({}, self.z_set)


    
#---------------------------------------------------------------------#
def basis_check(basis):
    if basis not in ['X', 'Z']:
        raise ValueError("basis must be 'X' or 'Z', "
                            "{} entered.".format(basis))