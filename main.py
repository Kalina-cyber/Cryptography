from time import time
from itertools import product

# Параметри кривої BN254
p = 0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47
q = 0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001
b = 3

"""
A1: Розглянемо ендоморфiзм φ : E(F_p) → E(F_p), що визначений як φ(x, y) = (ωx, y) 
для певної константи ω ∈ F_p. Напишiть програму, що за простим числом p, обчислює 
нетривiальне значення ω (ω≠1) таке, що φ дiйсно є ендоморфiзмом. Такий ендоморфiзмом
 в лiтературi називають корисним/GLV ендоморфiзмом.
"""


# A1: Пошук примітивного кубічного кореня з 1 (ω ≠ 1)
def find_omega():
    for x in range(2, p):
        if pow(x, 3, p) == 1 and x != 1:
            return x
    return None


# Операції на еліптичній кривій
def point_add(P, Q):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if P == Q:
        m = (3 * x1 * x1) * pow(2 * y1, -1, p) % p
    else:
        m = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (m * m - x1 - x2) % p
    y3 = (m * (x1 - x3) - y1) % p
    return (x3, y3)


def point_double(P):
    return point_add(P, P) if P is not None else None


def scalar_multiplication(k, P):
    result = None
    for bit in bin(k)[2:]:
        result = point_double(result)
        if bit == '1':
            result = point_add(result, P)
    return result


"""
A2: Оскiльки група E(F_p) є циклiчною, то справедливо φ(P) = [λ]P для деякої константи 
λ ∈ F_p i для кожного P ∈ E(F_p). Визначте λ для вашого ω.
"""


# A2: Пошук λ такого, що φ(P) = [λ]P
def find_lambda(omega, P, max_attempts=1 << 20):
    phiP = (omega * P[0] % p, P[1])
    for lam in range(1, max_attempts):
        if scalar_multiplication(lam, P) == phiP:
            return lam
    return None


"""
A3: Отже, нехай нам потрiбно обчислити [α]P, де α ∈ F_q. В такому разi головна iдея GLV 
полягає в розкладаннi α наступним чином: α = α_1 + α_2*λ (mod q), де бiтовий розмiр α_1 
та α_2 наближено вдвiчi менший нiж бiтовий розмiр α. Далi, достатньо обрахувати: [α]P = 
[α_1 + α_2*λ]P = [α_1]P ⊕ [α_2][λ]P = [α_1]P ⊕ [α_2]φ(P) Напишiть програму, що за допомогою
 α знаходить α_1, α_2 так, що |α_1|, |α_2| ≈ √q.

Алгоритм 3.47 з книги: Balanced length-two representation of a multiplier:
INPUT: Integers n, λ, k ∈ [0,n −1].
OUTPUT: Integers k_1, k_2 such that k = k_1 +k_2*λ mod n and |k_1|,|k_2| ≈ √n.
1. Run the extended Euclidean algorithm (Algorithm 2.19) with inputs n and λ. The algorithm 
produces a sequence of equations (s_i)*n +(t_i)*λ = r_i where s_0 = 1, t_0 = 0, r_0 = n, 
s_1 = 0, t_1 = 1, r_1 = λ, and the remainders ri and are non-negative and strictly decreasing. 
Let l be the greatest index for which r_l ≥ √n.
2. Set (a_1,b_1)←(r_(l+1),−t_(l+1)).
3. If (((r_l)^2) +((t_l)^2)) ≤ (((r_(l+2))^2) +((t_(l+2))^2)) then set (a_2,b_2)←(r_l,−t_l);
Else set (a_2,b_2)←(r_(l+2),−t_(l+2)).
4. Compute c_1 = [(b_2)k/n] and c_2 = [-(b_1)k/n].
5. Compute k_1 = k −c_1*a_1 −c_2*a_2 and k_2 = −c_1*b_1 −c_2*b_2.
6. Return(k_1, k_2).
"""


# A3: Розклад α = α1 + α2·λ (mod q)
def balanced_decomposition(k, lam, n):
    s0, t0, r0 = 1, 0, n
    s1, t1, r1 = 0, 1, lam
    seq = [(s0, t0, r0), (s1, t1, r1)]

    while r1 != 0:
        q_ = r0 // r1
        r2 = r0 - q_ * r1
        s2 = s0 - q_ * s1
        t2 = t0 - q_ * t1
        seq.append((s2, t2, r2))
        s0, t0, r0 = s1, t1, r1
        s1, t1, r1 = s2, t2, r2

    sqrt_n = int(n ** 0.5)
    l = max(i for i, (_, _, r) in enumerate(seq) if r >= sqrt_n)
    r_l, t_l = seq[l][2], seq[l][1]
    r_lp1, t_lp1 = seq[l + 1][2], seq[l + 1][1]

    if l + 2 < len(seq):
        r_lp2, t_lp2 = seq[l + 2][2], seq[l + 2][1]
        if r_l ** 2 + t_l ** 2 <= r_lp2 ** 2 + t_lp2 ** 2:
            a2, b2 = r_l, -t_l
        else:
            a2, b2 = r_lp2, -t_lp2
    else:
        a2, b2 = r_l, -t_l

    a1, b1 = r_lp1, -t_lp1
    c1 = (b2 * k) // n
    c2 = -(b1 * k) // n
    k1 = k - c1 * a1 - c2 * a2
    k2 = -c1 * b1 - c2 * b2
    return k1 % q, k2 % q


"""
A4: Обчислити лiнiйну комбiнацiю [α_1]P1 ⊕ [α_2]P_2 можна значно швидше, нiж спочатку окремо 
обчислити [α_1]P_1 та [α_2]P_2, а потiм їх додати. Саме тому вираз [α_1]P ⊕ [α_2]φ(P) можна 
обчислити вельми швидко. Напишiть програму, що реалiзує обчислення виразу [α_1]P_1 ⊕ [α_2]P_2 
за заданими P_1, P_2 ∈ E(F_p), α_1, α_2 ∈ F_q.

Алгоритм 3.48 з книги: Simultaneous multiple point multiplication: INPUT: Window width 
w, k = (k_(t−1),..., k_0)_2, l = (l_(t−1),...,l_0)_2, P, Q ∈ E(F_q ).
OUTPUT: kP +lQ.
1. Compute iP + jQ for all i, j ∈ [0,(2^w) −1].
2. Write k = (K^(d−1),..., K^1, K^0) and l = (L^(d−1),..., L^1, L^0) where each K^i , 
L^i is a bitstring of length w, and d = [t/w].
3. R←∞.
4. For i from d −1 downto 0 do
4.1 R←(2^w) R.
4.2 R← R +((K^i)P + (L^i)Q).
5. Return(R).
"""


# A4: Обчислення [α1]P + [α2]Q одночасно (Shamir’s trick)
def shamir_trick(alpha1, P, alpha2, Q, w=4):
    table = {}
    for i, j in product(range(2 ** w), repeat=2):
        Pi = scalar_multiplication(i, P) if i else None
        Qj = scalar_multiplication(j, Q) if j else None
        table[(i, j)] = point_add(Pi, Qj)

    def split_into_windows(k, w):
        bits = bin(k)[2:].zfill(((len(bin(k)) - 2 + w - 1) // w) * w)
        return [int(bits[i:i + w], 2) for i in range(0, len(bits), w)]

    a1_windows = split_into_windows(alpha1, w)
    a2_windows = split_into_windows(alpha2, w)

    R = None
    for i in range(len(a1_windows)):
        if R is not None:
            for _ in range(w):
                R = point_double(R)
        T = table[(a1_windows[i], a2_windows[i])]
        R = point_add(R, T)
    return R


"""
A5: Складiть усi пункти до цього моменту в один алгоритм, що реалiзує GLV множення. Ваша 
програма повинна приймати на вхiд P ∈ E(F_p), α ∈ F_q та повертати [α]P. Покажiть, що ваш 
алгоритм дiйсно швидший за метод подвоєнь i додавань. Для цього, реалiзуйте обидва алгоритми 
та порiвняйте їх швидкiсть.
"""


# A5: GLV множення
def glv_multiplication(alpha, P, w=4):
    omega = find_omega()
    lam = find_lambda(omega, P)
    #omega = 0x5bf3adda19e9b27bc1d99b7c07f3d6f0eb3af00adb22a5a0f44fdb2e5bd1e10
    #lam = 0x30644e72e131a029b85045b68181585d286c31cf5cf5d3ed
    if lam is None:
        raise ValueError("λ не знайдено")
    phiP = (omega * P[0] % p, P[1])
    a1, a2 = balanced_decomposition(alpha, lam, q)
    return shamir_trick(a1, P, a2, phiP, w)


P = (1, 2)
alpha = 0x123456789abcdef123456789abcdef123456789abcdef

start1 = time()
res1 = scalar_multiplication(alpha, P)
end1 = time()

start2 = time()
res2 = glv_multiplication(alpha, P)
end2 = time()

print("Метод подвоєнь і додавань:", res1, f"{end1 - start1:.6f} секунд")
print("Метод GLV:", res2, f"{end2 - start2:.6f} секунд")
print("Результати збігаються:", res1 == res2)