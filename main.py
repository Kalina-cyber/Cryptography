from random import randint
import random
from math import isqrt
import time

# Параметри кривої BN254
p = 0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47
q = 0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001
b = 3

def is_on_curve(x, y):
    return (y * y - x * x * x - b) % p == 0

def find_nontrivial_cubic_root(p):
    for omega in range(2, p):
        if pow(omega, 3, p) == 1 and omega != 1:
            return omega
    return None

def generate_random_point():
    while True:
        x = randint(1, p - 1)
        y2 = (x ** 3 + b) % p
        y = pow(y2, (p + 1) // 4, p)
        if (y * y) % p == y2:
            return (x, y)

def point_add(P, Q):
    if P == Q:
        return point_double(P)
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None

    m = ((y2 - y1) * pow(x2 - x1, -1, p)) % p
    x3 = (m * m - x1 - x2) % p
    y3 = (m * (x1 - x3) - y1) % p
    return (x3, y3)

def point_double(P):
    if P is None:
        return None
    x, y = P
    if y == 0:
        return None
    m = ((3 * x * x) * pow(2 * y, -1, p)) % p
    x3 = (m * m - 2 * x) % p
    y3 = (m * (x - x3) - y) % p
    return (x3, y3)

def scalar_mult(k, P):
    R = None
    while k > 0:
        if k & 1:
            R = point_add(R, P)
        P = point_double(P)
        k >>= 1
    return R

def find_lambda(P, phi_P):
    for lam in range(1, 1000000):  # Обмежений пошук
        if scalar_mult(lam, P) == phi_P:
            return lam
    return None

# -------------------- A3: GLV decomposition ----------------------

def extended_euclidean(a, b):
    u, v = a, b
    x1, y1 = 1, 0
    x2, y2 = 0, 1
    sequence = [(x1, y1, u), (x2, y2, v)]

    while u != 0:
        q = v // u
        r = v - q * u
        x = x2 - q * x1
        y = y2 - q * y1
        y = y2 - q * y1
        sequence.append((x, y, r))
        v, u = u, r
        x2, x1 = x1, x
        y2, y1 = y1, y

    return sequence

def glv_decompose(n, lam, k):
    sequence = extended_euclidean(n, lam)
    sqrt_n = isqrt(n)

    l = max(i for i, (_, _, r) in enumerate(sequence) if r >= sqrt_n)

    r_l1, t_l1 = sequence[l+1][2], sequence[l+1][1]
    r_l2, t_l2 = sequence[l+2][2], sequence[l+2][1]
    r_l,  t_l  = sequence[l][2], sequence[l][1]

    a1, b1 = r_l1, -t_l1
    norm_l  = r_l**2 + t_l**2
    norm_l2 = r_l2**2 + t_l2**2

    if norm_l <= norm_l2:
        a2, b2 = r_l, -t_l
    else:
        a2, b2 = r_l2, -t_l2

    c1 = round(b2 * k / n)
    c2 = round(-b1 * k / n)

    k1 = k - c1 * a1 - c2 * a2
    k2 = -c1 * b1 - c2 * b2

    return k1 % n, k2 % n

def to_binary_chunks(k, w, t):
    bin_k = bin(k)[2:].zfill(t)
    return [int(bin_k[i:i+w], 2) for i in range(0, t, w)]

def precompute_table(P, Q, w):
    table = {}
    max_val = 2**w
    for i in range(max_val):
        for j in range(max_val):
            if i == 0 and j == 0:
                continue
            table[(i, j)] = None
            Pi = scalar_mult(i, P) if i > 0 else None
            Qj = scalar_mult(j, Q) if j > 0 else None
            table[(i, j)] = point_add(Pi, Qj)
    return table

def simultaneous_mult(k, l, P, Q, w=2):
    # t — бітова довжина
    t = max(k.bit_length(), l.bit_length())
    d = (t + w - 1) // w

    # Розбиваємо на блоки по w бітів
    K_chunks = to_binary_chunks(k, w, d * w)
    L_chunks = to_binary_chunks(l, w, d * w)

    # Попередні обчислення
    table = precompute_table(P, Q, w)

    R = None
    for i in range(d):
        # w-кратне дублювання
        for _ in range(w):
            R = point_double(R)
        K_i = K_chunks[i]
        L_i = L_chunks[i]
        if (K_i, L_i) != (0, 0):
            addition = table.get((K_i, L_i))
            R = point_add(R, addition)
    return R

def phi_endomorphism(P, omega):
    x, y = P
    return ((omega * x) % p, y)

def decompose_scalar(alpha, lam, q):
    return glv_decompose(q, lam, alpha)

# GLV множення: [α]P = [α1]P + [α2]φ(P)
def glv_multiply(P, alpha, q, lambda_glv, omega):
    # Обчислити φ(P) = (ωx, y)
    phi_P = phi_endomorphism(P, omega)

    # Розкласти α як α = α1 + α2 * λ (mod q)
    alpha1, alpha2 = decompose_scalar(alpha, lambda_glv, q)

    # Обчислити [α1]P + [α2]φ(P) за допомогою Shamir's trick
    R = simultaneous_mult(alpha1, alpha2, P, phi_P)
    return R

def double_and_add(k, P):
    R = None
    Q = P
    while k > 0:
        if k & 1:
            R = point_add(R, Q)
        Q = point_double(Q)
        k >>= 1
    return R

# ----------------- Виконання всіх частин --------------------

#omega = find_nontrivial_cubic_root(p)
omega = 0x5bf3adda19e9b27bc1d99b7c07f3d6f0eb3af00adb22a5a0f44fdb2e5bd1e10  # кубічний корінь з 1, ω ≠ 1 mod p

print(f"ω = {omega}")

P = generate_random_point()
print(f"P = {P}")

phi_P = ((omega * P[0]) % p, P[1])

#lam = find_lambda(P, phi_P)
lam = 0x30644e72e131a029b85045b68181585d286c31cf5cf5d3ed  # відоме λ для BN254
print(f"λ = {lam}")

# Випадкове α ∈ [1, q − 1]
alpha = randint(1, q - 1)
alpha1, alpha2 = glv_decompose(q, lam, alpha)
print(f"α = {alpha}")
print(f"α₁ = {alpha1}")
print(f"α₂ = {alpha2}")

# Перевірка: [α]P = [α₁]P + [α₂]φ(P)
R1 = scalar_mult(alpha1, P)
R2 = scalar_mult(alpha2, phi_P)
R  = point_add(R1, R2)
expected = scalar_mult(alpha, P)

print(f"GLV: [α₁]P + [α₂]φ(P) = {R}")
print(f"Direct: [α]P = {expected}")
print(f"Correct: {R == expected}")

# A4: швидке обчислення [α₁]P ⊕ [α₂]φ(P)
R_glv_fast = simultaneous_mult(alpha1, alpha2, P, phi_P, w=2)
print(f"Shamir’s trick: [α₁]P + [α₂]φ(P) = {R_glv_fast}")
print(f"Matches direct scalar multiplication: {R_glv_fast == expected}")

# Генерація випадкового α
alpha = random.randint(1, q - 1)

# Обчислення через GLV
start_glv = time.time()
glv_result = glv_multiply(P, alpha, q, lam, omega)
end_glv = time.time()

# Обчислення через double-and-add
start_dna = time.time()
dna_result = double_and_add(alpha, P)
end_dna = time.time()

# Перевірка правильності та вивід часу
print(f"GLV result correct: {glv_result == dna_result}")
print(f"GLV time: {end_glv - start_glv:.6f} seconds")
print(f"Double-and-add time: {end_dna - start_dna:.6f} seconds")