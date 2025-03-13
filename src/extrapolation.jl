using PhysicalConstants.CODATA2018: m_e, ħ

export MVL, MVLK, APL, APLK, P3, K3, G3S, G3M


# Z atomic number -- protons
x(V,V0) = (V / V0)^(1/3)
y(V,V0) = 1 - x(V,V0)

# Vinet
MVL(V,V0,K0,c) = 3*K0*x(V,V0)^-2 * y(V,V0) * exp(c * y(V,V0))
MVLK(V,V0,K0,c) = K0*x(V,V0)^-2*exp(c*y(V,V0))*(2 - x(V,V0) + c*x(V,V0)*y(V,V0))

# Holzapfel 98
const a_fg = ((3*π^2)^(2/3)/5)* ħ.val^2 / m_e.val * ((1e10)^5 / 1e9) # GPa * Å⁵
c_0(V0,K0,Z) = -log(3*K0 / (a_fg*(Z/V0)^(5/3)))
c_2(c0, K1) = 3/2*(K1 - 3) - c0
APL(V,V0,K0,c0,c2) = 3*K0*x(V,V0)^-5 * y(V,V0) * exp(c0*y(V,V0)) * (1 + x(V,V0)*c2*y(V,V0))
APLK(V,V0,K0,c0,c2) = K0/V*x(V,V0)^-2*exp(c0*y(V,V0))*(c2*V*(2+c0*(-2+x(V,V0))) + 5*V0 + V0*x(V,V0)*(-4+c0-c0*x(V,V0)) + c2*x(V,V0)*V0*(4-6*x(V,V0)+c0*x(V,V0)))

# Melinger-Cohen & Jeanloz '19 BM3
a(K0) = 9*K0
b(n, K0, K1) = -27*K0 * (2 + n - K1)
f(n, V, V0) = 1/n * ((V/V0)^(-n/3) - 1)
C(n, V, V0) = (1 + n*f(n,V,V0))^((3+n)/n)

# Pressure
P3(V, V0, n, K0, K1) = 1/3*f(n,V,V0)*C(n,V,V0) * (a(K0) + 1/2*b(n,K0,K1)*f(n,V,V0))
# Bulk Modulus (pressure dervative)
K3(V, V0, n, K0, K1) = 1/18*C(n,V,V0) * (2*a(K0) + (2*b(n,K0,K1) + a(K0)*(6 + 4*n))*f(n,V,V0) + 3*b(n,K0,K1)*(1 + n)*f(n,V,V0)^2)
# Shear Modulus
G3S(V, V0, K0, K1, G0, G1) = C(2,V,V0) * (G0 + (3*K0*G1 - 5*G0)*f(2,V,V0) + (K0*(9/2*K1 - 24 + 6*G1) - 14*G0)*f(2,V,V0)^2)
G3M(V, V0, K0, K1, G0, G1) = C(-2,V,V0) * (G0 + (3*K0*G1 - G0)*f(-2,V,V0) - (2*G0 + K0*(6 + 6*G1 + 9/2*K1))*f(-2,V,V0)^2)
