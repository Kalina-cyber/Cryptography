import random
import time

# Параметри кривої BN254
p = 0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47
q = 0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001
b = 3

# A1: Обчислення корисного ендоморфізму через ітерацію кубічного кореня
omega = pow(2, (p - 1) // 3, p)
while omega == 1:
    omega += 1
    omega = pow(omega, (p - 1) // 3, p)

# A2: Фі(P) = [lambda]P. Подаемо просто скалярне добуток з P до ф(P) і знаходимо лямбду

def point_add(P, Q):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and y1 != y2: return None
    if P == Q:
        l = (3 * x1 * x1) * pow(2 * y1, -1, p) % p
    else:
        l = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (l * l - x1 - x2) % p
    y3 = (l * (x1 - x3) - y1) % p
    return (x3, y3)

def point_double(P):
    return point_add(P, P)

def scalar_mult(k, P):
    return double_and_add(k, P)

def phi_endomorphism(P, omega):
    if P is None:
        return None
    x, y = P
    return (omega * x % p, y)

# Приклад точки P на кривій
x = 1
rhs = (x**3 + b) % p
while pow(rhs, (p - 1) // 2, p) != 1:
    x += 1
    rhs = (x**3 + b) % p

y = pow(rhs, (p + 1) // 4, p)
P = (x, y)

phi_P = phi_endomorphism(P, omega)

# Знаходимо lambda із рівності скалярів
lambda_glv = None
for lam in range(2, 1000):
    if scalar_mult(lam, P) == phi_P:
        lambda_glv = lam
        break

# A3: Розклад скаляра за допомогою ЕРС і методу з книги

def extended_euclid(a, b):
    s0, s1 = 1, 0
    t0, t1 = 0, 1
    r0, r1 = a, b
    seq = [(s0, t0, r0), (s1, t1, r1)]
    while r1:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
        t0, t1 = t1, t0 - q * t1
        seq.append((s1, t1, r1))
    return seq

def decompose_scalar(k, lam, n):
    seq = extended_euclid(n, lam)
    sqrt_n = int(n**0.5)
    l = max(i for i, (_, _, r) in enumerate(seq) if abs(r) >= sqrt_n)

    r_l, t_l = seq[l][2], seq[l][1]
    r_lp1, t_lp1 = seq[l+1][2], seq[l+1][1]
    r_lp2, t_lp2 = seq[l+2][2], seq[l+2][1]

    a1, b1 = r_lp1, -t_lp1
    if r_l**2 + t_l**2 <= r_lp2**2 + t_lp2**2:
        a2, b2 = r_l, -t_l
    else:
        a2, b2 = r_lp2, -t_lp2

    c1 = round(b2 * k / n)
    c2 = round(-b1 * k / n)

    k1 = (k - c1 * a1 - c2 * a2) % n
    k2 = (-c1 * b1 - c2 * b2) % n
    return k1, k2

# A4: Шамірова трюка для [k1]P + [k2]Q

def simultaneous_mult(k1, k2, P, Q):
    table = {
        (0, 0): None,
        (0, 1): Q,
        (1, 0): P,
        (1, 1): point_add(P, Q),
    }
    R = None
    for i in reversed(range(max(k1.bit_length(), k2.bit_length()))):
        R = point_double(R)
        b1 = (k1 >> i) & 1
        b2 = (k2 >> i) & 1
        if (b1, b2) != (0, 0):
            R = point_add(R, table[(b1, b2)])
    return R

# A5: GLV множення

def glv_multiply(P, alpha, q, lambda_glv, omega):
    phi_P = phi_endomorphism(P, omega)
    alpha1, alpha2 = decompose_scalar(alpha, lambda_glv, q)
    return simultaneous_mult(alpha1, alpha2, P, phi_P)

def double_and_add(k, P):
    R = None
    Q = P
    while k > 0:
        if k & 1:
            R = point_add(R, Q)
        Q = point_double(Q)
        k >>= 1
    return R

# Тестування та порівняння швидкості
alpha = random.randint(1, q - 1)

start_glv = time.time()
glv_result = glv_multiply(P, alpha, q, lambda_glv, omega)
end_glv = time.time()

start_dna = time.time()
dna_result = double_and_add(alpha, P)
end_dna = time.time()

print(f"GLV result correct: {glv_result == dna_result}")
print(f"GLV time: {end_glv - start_glv:.6f} seconds")
print(f"Double-and-add time: {end_dna - start_dna:.6f} seconds")
