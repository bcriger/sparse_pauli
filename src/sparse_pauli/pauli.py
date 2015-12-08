__all__ = ["Pauli"]

class Pauli(object):
    """
    Stores the X and Z support of a (phaseless) Pauli as two sets.

    Supports multiplication, commutation, Hadamard and CNOT, but 
    nothing else.
    """
    def __init__(self, x_set, z_set):
        self.x_set = set(x_set)
        self.z_set = set(z_set)

    #printing/comparison
    def __str__(self):
        string = ''
        
        support = sorted(list(self.x_set | self.z_set))
        
        if len(support) == 0:
            return 'I'
        
        for elem in support:
            if elem in self.x_set:
                char = 'Y' if elem in self.z_set else 'X'
            else:
                char = 'Z'

            string += '{}[{}] '.format(char, elem)
        return string

    def __eq__(self, other):
        if self.x_set == other.x_set:
            return self.z_set == other.z_set
        else:
            return False
    
    #actual math
    def __mul__(self, other):
        return Pauli(self.x_set ^ other.x_set,
                        self.z_set ^ other.z_set)

    def com(self, other):
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
    
    def h(self, qs):
        """
        acts a Hadamard on each bit in qs.
        """
        switches = (self.x_set ^ self.z_set) & set(qs)
        self.x_set ^= switches
        self.z_set ^= switches
        pass
