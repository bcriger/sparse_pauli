__all__ = ["Pauli"]

class Pauli(object):
    """
    Stores the X and Z support of a (phaseless) Pauli as two sets.

    Supports multiplication, commutation, Hadamard and CNOT, but 
    nothing else.
    """
    def __init__(self, x_set={}, z_set={}):
        self.x_set = set(x_set)
        self.z_set = set(z_set)

    #printing/comparison
    def __repr__(self):
        string = ''
        
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

    def __eq__(self, othr):
        return (self.x_set == othr.x_set) & (self.z_set == othr.z_set)
    
    def __ne__(self, othr):
        return (self.x_set != othr.x_set) | (self.z_set != othr.z_set)
    
    #actual math
    def __mul__(self, other):
        return Pauli(self.x_set ^ other.x_set,
                        self.z_set ^ other.z_set)

    def com(self, other):
        """
        This method suffers from the usual ambiguity,
        com(self, other) == 1 means they ANTIcommute.
        """
        return (len(self.x_set & other.z_set) +
                    len(self.z_set & other.x_set)) % 2

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
        pass

    def prep(self, qs, basis='Z'):
        """
        Flawlessly prepares qubits in the list `qs` in the +1 
        eigenstate of the selected basis, which is Z by default, and
        must be either Z or X.
        """
        basis_check(basis)
        
        if basis == 'X':
            self.x_set |= set(qs)
        elif basis == 'Z':
            self.z_set |= set(qs) 

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
        return self.x_set.union(self.z_set)

    def I():
        return Pauli({}, {})

#---------------------------------------------------------------------#
def basis_check(basis):
    if basis not in ['X', 'Z']:
        raise ValueError("basis must be 'X' or 'Z', "
                            "{} entered.".format(basis))