import sparse_pauli as sp

x = sp.Pauli({0}, {})
z = sp.Pauli({}, {0})
i = sp.Pauli({}, {})
y = sp.Pauli({0}, {0})
ix = sp.Pauli({1}, {})
iz = sp.Pauli({}, {1})

def com_check(pauli_1, pauli_2, desired_output):
    return pauli_1.com(pauli_2) == desired_output

def test_com_zz():
    assert com_check(z, z, 0)

def test_com_xx():
    assert com_check(x, x, 0)

def test_com_xz():
    assert com_check(x, z, 1)

def test_com_zx():
    assert com_check(z, x, 1)

def mul_test_yy():
    assert x*z*ix*iz == sp.Pauli({0, 1}, {0, 1})

def h_test():
    big_pauli = sp.Pauli(set(range(4)), {})
    bits = {0, 2}
    big_pauli.h(bits)
    assert big_pauli == sp.Pauli({1, 3}, {0, 2})

def cnot_test():
    big_pauli = sp.Pauli(set(range(6)), set(range(6)))
    big_pauli.cnot({(0, 1),(2, 3),(4, 5)})
    assert big_pauli == sp.Pauli({0, 2, 4}, {1, 3, 5})

def cz_test():
    big_pauli = sp.Pauli(set(range(6)), {})
    big_pauli.cz({(0, 1),(2, 3),(4, 5)})
    assert big_pauli == sp.Pauli(set(range(6)), set(range(6)))
