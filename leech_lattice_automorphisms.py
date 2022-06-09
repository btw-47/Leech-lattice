import cypari2
pari = cypari2.Pari()

from sage.arith.functions import lcm
from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix
from sage.misc.functional import denominator
from sage.modules.free_module_element import vector
from sage.quadratic_forms.quadratic_form import QuadraticForm
from sage.rings.integer_ring import ZZ
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.rational_field import QQ


_L = matrix([[8]+[0]*23, [4,4]+[0]*22, [4,0,4]+[0]*21, [4,0,0,4]+[0]*20, [4,0,0,0,4]+[0]*19, [4,0,0,0,0,4]+[0]*18, [4,0,0,0,0,0,4]+[0]*17, [2]*8+[0]*16, [4,0,0,0,0,0,0,0,4]+[0]*15, [4,0,0,0,0,0,0,0,0,4]+[0]*14, [4,0,0,0,0,0,0,0,0,0,4]+[0]*13, [2,2,2,2,0,0,0,0,2,2,2,2]+[0]*12, [4,0,0,0,0,0,0,0,0,0,0,0,4]+[0]*11, [2,2,0,0,2,2,0,0,2,2,0,0,2,2]+[0]*10, [2,0,2,0,2,0,2,0,2,0,2,0,2,0,2]+[0]*9, [2,0,0,2,2,0,0,2,2,0,0,2,2,0,0,2]+[0]*8, [4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4]+[0]*7, [2,0,2,0,2,0,0,2,2,2,0,0,0,0,0,0,2,2]+[0]*6, [2,0,0,2,2,2,0,0,2,0,2,0,0,0,0,0,2,0,2]+[0]*5, [2,2,0,0,2,0,2,0,2,0,0,2,0,0,0,0,2,0,0,2,0,0,0,0],[0,2,2,2,2,0,0,0,2,0,0,0,2,0,0,0,2,0,0,0,2,0,0,0],[0,0,0,0,0,0,0,0,2,2,0,0,2,2,0,0,2,2,0,0,2,2,0,0],[0,0,0,0,0,0,0,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0],[-3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]])
_Lt = _L.transpose()
_Leech = _L*_Lt / 8

def Leech():
    return _Leech


class LeechLatticeAutomorphism:

    def __init__(self, g):
        self.__matrix = matrix(QQ, g)

    def __repr__(self):
        return 'Leech lattice automorphism\n%s'%self.matrix()

    def matrix(self):
        return self.__matrix

    ### Arithmetic operations

    def __invert__(self):
        return LeechLatticeAutomorphism(self.matrix().inverse())

    def __mul__(self, v):
        if isinstance(v, LeechLatticeAutomorphism):
            return LeechLatticeAutomorphism(self.matrix() * v.matrix())
        return self.matrix() * v

    def __neg__(self):
        return LeechLatticeAutomorphism(-self.matrix())

    def __pow__(self, n):
        if n == 1:
            return self
        elif n == 0:
            return LeechLatticeAutomorphism(identity_matrix(QQ, 24))
        elif n > 1:
            n_half = n // 2
            return self.__pow__(n_half) * self.__pow__(n - n_half)
        elif n < 0:
            return (~self).__pow__(-n)

    def __rmul__(self, v):
        if isinstance(v, LeechLatticeAutomorphism):
            return LeechLatticeAutomorphism(v.matrix() * self.matrix())
        return v * self.matrix()

    def transpose(self):
        return self.matrix().transpose()

    ### Basic properties

    def cycle_shape(self):
        try:
            return self.__cycle_shape
        except AttributeError:
            pass
        g = self.matrix()
        f = g.charpoly()
        i = 1
        s = ''
        d = {}
        R, t = PowerSeriesRing(QQ, 't', 100).objgen()
        f = (R([round(x) for x in f.list()])).log()
        while f and i <= 100:
            h = (1 - t**i).log()
            n = f[i]
            if n:
                s += '*%s^%s'%(i, -n)
            f += n*h
            d[i] = -n
            i += 1
        s = s[1:]
        self.__cycle_shape_dict = d
        self.__cycle_shape = s
        return s

    def cycle_shape_dict(self):
        try:
            return self.__cycle_shape_dict
        except AttributeError:
            _ = self.cycle_shape()
            return self.__cycle_shape_dict

    def level(self):
        try:
            return self.__level
        except AttributeError:
            pass
        d = self.cycle_shape_dict()
        h = max(c for c, u in d.items() if u)
        l = sum(u / c for c, u in d.items())
        n = h * denominator(h * l / 24)
        self.__level = n
        return n

    def order(self):
        try:
            return self.__order
        except AttributeError:
            pass
        d = self.cycle_shape_dict()
        n = lcm([x for x in d.keys() if d[x]])
        self.__order = n
        return n

    ### fixed points

    def fixed_lattice(self):
        V = (self.matrix().transpose() - identity_matrix(24)).integer_kernel()
        return LeechSublattice(matrix(ZZ, V.basis_matrix()))

    def Lambda_plus(self, d):
        n = self.order()
        b = (self ** d).fixed_lattice()
        if n % 2 or d % 2:
            return b, []
        b = b.basis_matrix()
        y = self ** (d // 2)
        a = matrix(ZZ, b * Leech() *  y.matrix() * b.transpose())
        I = [i for i, x in enumerate(a.diagonal()) if x % 2]
        c = identity_matrix(a.nrows())
        for i in range(len(I) - 1):
            c[I[i], I[i+1]] = 1
        if len(I) > 1:
            c[I[-1], I[-2]] = -1
        elif len(I) == 1:
            c[I[0], I[0]] = 2
        return LeechSublattice(c * b), I

class LeechSublattice:
    def __init__(self, basis_matrix):
        self.__basis_matrix = basis_matrix

    def __repr__(self):
        return 'Sublattice of the Leech lattice spanned by\n'+'\n'.join(list(map(str, self.basis_matrix().rows())))

    def basis_matrix(self):
        return self.__basis_matrix

    def gram_matrix(self):
        try:
            return self.__gram_matrix
        except AttributeError:
            x = self.basis_matrix()
            s = x * Leech() * x.transpose()
            self.__gram_matrix = s
            return s

    def _quadratic_form(self):
        return QuadraticForm(ZZ, self.gram_matrix())

    def _span(self):
        return span(self.basis_matrix())

    def intersection(self, other):
        return LeechSublattice( self._span().intersection(other._span()).basis_matrix() )


    ### attributes

    def discriminant(self):
        try:
            return self.__discriminant
        except AttributeError:
            d = self.gram_matrix().determinant()
            self.__discriminant = d
            return d

    def genus(self):
        return self._quadratic_form().global_genus_symbol()

    def level(self):
        try:
            return self.__level
        except AttributeError:
            l = self._quadratic_form().level()
            self.__level = l
            return l

    def rank(self):
        return self.basis_matrix().nrows()

    def short_vectors_in_dual(self, N):
        r"""
        Compute vectors in self's dual lattice of norm at most N as a list.
        """
        s = self.gram_matrix()
        s_inv = s.inverse()
        _, _, vs_list = pari(s_inv).qfminim(N + N, flag=2)
        vs_list = vs_list.sage().columns()
        vs_list = vs_list + [-x for x in vs_list] + [vector([0] * s.nrows())]
        return [s_inv * v for v in vs_list]


    ### vectors

    def __contains__(self, v):
        try:
            _ = self.basis_matrix().solve_left(v)
            return True
        except ValueError:
            return False

    def orthogonal_projection(self, v):
        try:
            p = self.__projection_matrix
        except AttributeError:
            s = self.gram_matrix()
            b = self.basis_matrix()
            p = b.transpose() * s.inverse() * b * Leech()
            self.__projection_matrix = p
        return p * v
