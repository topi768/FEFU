import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from sympy.codegen.rewriting import optimize
from scipy import optimize

# ------------------------
# Настройки и функция
# ------------------------
x_sym = sp.symbols('x')

a = 0.4
b = 0.9
x_star = 0.52
count_points = 11

# f_expr можно менять в одном месте
f_expr = x_sym**2 + sp.log(x_sym)
# f_expr = x_sym**2 - sp.log(x_sym + 2, 10)
# f_expr = x_sym**2 - sp.log(x_sym, 10)
# f_expr = x_sym**3 - sp.sin(x_sym)

# ------------------------
# Часть 1: определяем f
# ------------------------
f = sp.lambdify(x_sym, f_expr, modules=['numpy'])

x_points = np.linspace(a, b, count_points)
y_points = f(x_points)

x_fit = np.linspace(a, b, 500)
y_fit = f(x_fit)

# ------------------------
# Часть 2: Лагранж 1-го порядка
# ------------------------
if x_star <= x_points[0]:
    i = 0
elif x_star >= x_points[-1]:
    i = len(x_points) - 2
else:
    i = np.searchsorted(x_points, x_star) - 1

x_im1 = sp.Rational(str(x_points[i-1]))
x_i = sp.Rational(str(x_points[i]))
x_ip1 = sp.Rational(str(x_points[i+1]))

y_i_sym = sp.simplify(f_expr.subs(x_sym, x_i))
y_i1_sym = sp.simplify(f_expr.subs(x_sym, x_ip1))


x_star_sym = sp.Rational(str(x_star))
x_star_num = float(x_star_sym)

L1 = sp.simplify(
    y_i_sym * (x_star_sym - x_ip1) / (x_i - x_ip1) +
    y_i1_sym * (x_star_sym - x_i) / (x_ip1 - x_i)
)



# print(f"L1(x*): {L1} (символьный)")
# print(f"L1(x*) (численно, 20 знаков): {sp.N(L1, 20)}")


# ------------------------
# Часть 3: Производная
# ------------------------
f_prime_expr = sp.diff(f_expr, x_sym)
f_prime_expr2 = sp.diff(f_prime_expr, x_sym)

min_val = sp.minimum(f_prime_expr2, x_sym, sp.Interval(x_i , x_ip1))
max_val = sp.maximum(f_prime_expr2, x_sym, sp.Interval(x_i, x_ip1))



omega_2_expr = (x_sym - x_i) * (x_sym - x_ip1)
omega_2_val = sp.simplify(omega_2_expr.subs(x_sym, x_star_sym))
# print(f'omega_2(x_star): {sp.simplify(omega_2_expr.subs(x_sym, x_star_sym))}')

R_1min_val = sp.simplify(min_val * omega_2_val / 2)
R_1max_val = sp.simplify(max_val * omega_2_val / 2)

if R_1min_val > R_1max_val:
    boofer = R_1min_val
    R_1min_val = R_1max_val
    R_1max_val = boofer


# ------------------------
# Часть 4: R1
# ------------------------

R1_expr = sp.simplify(f_expr - L1 )
R1_star_val = sp.simplify(R1_expr.subs(x_sym, x_star_sym))

print(f'R1_star: {R1_star_val}')

min_val_R_1_num = sp.N(R_1min_val, 30)
max_val_R_1_num = sp.N(R_1max_val, 30)
R_1_star_num = sp.N(R1_star_val, 30)

print(f'R_1_star: {R_1_star_num}')

if min_val_R_1_num  <R_1_star_num < max_val_R_1_num:
    print("YES")
else:
    print("NO")


counted_x = sp.simplify(L1.subs(x_sym, x_star_sym))
counted_x = sp.N(counted_x, 30)

factual_x = sp.simplify(f_expr.subs(x_sym, x_star_sym))
factual_x = sp.N(factual_x, 30)
print(f'counted_x: { counted_x} f(x): {factual_x}')


# ------------------------
# Часть 5: L2
# ------------------------


L2 = (
            sp.simplify(f_expr.subs(x_sym, x_im1)) * (x_star_sym - x_i) * (x_star_sym - x_ip1) / ((x_im1 - x_i) * (x_im1 - x_ip1))
            + sp.simplify(f_expr.subs(x_sym, x_i)) * (x_star_sym - x_im1) * (x_star_sym - x_ip1) / ((x_i - x_im1) * (x_i - x_ip1))
            + sp.simplify(f_expr.subs(x_sym, x_ip1)) * (x_star_sym - x_im1) * (x_star_sym - x_i) / ((x_ip1 - x_im1) * (x_ip1 - x_i))
          )

L2_num = sp.N(L2, 30)

# ------------------------
# Часть 6: R2
# ------------------------

omega_3_expr = (x_sym - x_im1) * (x_sym - x_i) * (x_sym-x_ip1)
omega_3_val = sp.simplify(omega_3_expr.subs(x_sym, x_star_sym))

f_prime_expr3 = sp.diff(f_prime_expr2, x_sym)

min_val = sp.minimum(f_prime_expr3, x_sym, sp.Interval(x_im1 , x_ip1))
max_val = sp.maximum(f_prime_expr3, x_sym, sp.Interval(x_im1, x_ip1))

R_2min_val = sp.simplify(min_val * omega_3_val / 6)
R_2max_val = sp.simplify(max_val * omega_3_val / 6)

if R_2min_val > R_2max_val:
    boofer = R_2min_val
    R_2min_val = R_2max_val
    R_2max_val = boofer

# ------------------------
# Часть 7: проверка
# ------------------------

R2_expr = sp.simplify(f_expr - L2 )
R2_star_val = sp.simplify(R2_expr.subs(x_sym, x_star_sym))


min_val_R_2_num = sp.N(R_2min_val, 30)
max_val_R_2_num = sp.N(R_2max_val, 30)
R_2_star_num = sp.N(R2_star_val, 30)

# print(f'R_1_star: {R_2_star_num}')
# print(f"min_val_R_2_num: {min_val_R_2_num}")
# print(f"max_val_R_2_num: {max_val_R_2_num}")


if min_val_R_2_num < R_2_star_num < max_val_R_2_num:
    print("YES")
else:
    print("NO")



counted_x = sp.simplify(L2.subs(x_sym, x_star_sym))
counted_x = sp.N(counted_x, 30)

factual_x = sp.simplify(f_expr.subs(x_sym, x_star_sym))
factual_x = sp.N(factual_x, 30)
print(f'counted_x: { counted_x} f(x): {factual_x}')

# ------------------------
# Часть 8: Таблица разделённых разностей и многочлены Ньютона
# ------------------------

y_im1_sym = sp.simplify(f_expr.subs(x_sym, x_im1))
y_i_sym   = sp.simplify(f_expr.subs(x_sym, x_i))
y_ip1_sym = sp.simplify(f_expr.subs(x_sym, x_ip1))

# 0-й порядок (таблица)
y0 = y_im1_sym
y1 = y_i_sym
y2 = y_ip1_sym

# 1-я разделённая разность
f01 = sp.simplify((y1 - y0) / (x_i - x_im1))   # f[x0,x1]
f12 = sp.simplify((y2 - y1) / (x_ip1 - x_i))   # f[x1,x2]



# 2-я разделённая разность
f012 = sp.simplify((f12 - f01) / (x_ip1 - x_im1))  # f[x0,x1,x2]

# Печать "таблицы"
print("\nТаблица разделённых разностей для узлов:")
print(f"x_{i-1} = {x_im1}, f = {y0}")
print(f"x_{i}   = {x_i},   f = {y1}")
print(f"x_{i+1} = {x_ip1}, f = {y2}")
print(f"f[x0,x1] = {f01}")
print(f"f[x1,x2] = {f12}")
print(f"f[x0,x1,x2] = {f012}")


# --- Многочлены Ньютона ---
x = x_sym

# линейный (Newton first-order) на базе x0,x1 с коэффициентом f01
N1 = sp.simplify(y0 + f01 * (x - x_im1))

# квадратичный Newton (N2) на x0,x1,x2
N2 = sp.simplify(y0 + f01 * (x - x_im1) + f012 * (x - x_im1) * (x - x_i))

# значения в x_star
N1_xstar_sym = sp.simplify(N1.subs(x_sym, x_star_sym))
N2_xstar_sym = sp.simplify(N2.subs(x_sym, x_star_sym))

N1_xstar_num = sp.N(N1_xstar_sym, 30)
N2_xstar_num = sp.N(N2_xstar_sym, 30)

print("\nМногочлены Ньютона:")
print(f"N1(x) (симв.): {N1}")
print(f"N2(x) (симв.): {N2}")
print(f"N1(x*): {N1_xstar_sym}  (числ.: {N1_xstar_num})")
print(f"N2(x*): {N2_xstar_sym}  (числ.: {N2_xstar_num})")

# Сравнение с Lagrange (L1, L2) — L1 и L2 у тебя есть
# 1) символьная разность (должна быть 0)
diff_N1_L1 = sp.simplify(N1.subs(x_sym, x_star_sym) - L1)   # если L1 — значение в x*, использовать L1 directly
# Если L1 у тебя — символьное значение именно в x_star, то оставь как есть.
# Чтобы быть надежным, рассматриваем L1_sym и L2_sym как симв. формы. Если L1/ L2 — численные значения, сравниваем численно:
try:
    L1_sym = L1  # у тебя L1 — уже символьное значение в x*
    L2_sym = L2  # L2 — символьный (значение в x*)
except Exception:
    L1_sym = sp.simplify(y0 + f01 * (x_star_sym - x_im1))  # fallback
    L2_sym = sp.simplify(N2_xstar_sym)

# Для надёжности: сравним численно N1(x*), N2(x*) с L1 и L2
N1_vs_L1_diff = float(sp.N(N1_xstar_sym - L1, 30))
N2_vs_L2_diff = float(sp.N(N2_xstar_sym - L2, 30))

print("\nСравнение Newton <-> Lagrange в точке x*:")
print(f"N1(x*) - L1(x*): {sp.simplify(N1_xstar_sym - L1)}  (числ: {N1_vs_L1_diff})")
print(f"N2(x*) - L2(x*): {sp.simplify(N2_xstar_sym - L2)}  (числ: {N2_vs_L2_diff})")

# Если хочешь, проверим символьное тождество на полиномах (всей формуле, а не только в точке):
sym_diff_N2_L2_full = sp.simplify(N2 - sp.simplify(
    # build Lagrange polynomial symbolically from nodes for 'x' to compare
    y0 * (x - x_i) * (x - x_ip1) / ((x_im1 - x_i) * (x_im1 - x_ip1)) +
    y1 * (x - x_im1) * (x - x_ip1) / ((x_i - x_im1) * (x_i - x_ip1)) +
    y2 * (x - x_im1) * (x - x_i) / ((x_ip1 - x_im1) * (x_ip1 - x_i))
))
print("\nСимвольная разность N2 - Lagrange_polynomial (должно быть 0):")
print(sym_diff_N2_L2_full)

# краткий итог в читабельном виде
print("\nИтоги (численно, 30 значащих):")
print(f"N1(x*): {N1_xstar_num}")
print(f"L1(x*): {sp.N(L1,30)}")
print(f"разность N1-L1: {N1_vs_L1_diff}")

print(f"N2(x*): {N2_xstar_num}")
print(f"L2(x*): {sp.N(L2,30)}")
print(f"разность N2-L2: {N2_vs_L2_diff}")
# ------------------------
# Графики
# ------------------------
plt.plot(x_fit, y_fit, color="red", label="f(x)")
plt.plot(x_points, y_points, "o", color="blue", label="табличные точки")
plt.axvline(x_star_num, color="gray", linestyle="--", alpha=0.5)
plt.xlabel('Ось х')
plt.ylabel('Ось y')
plt.title('График функции и интерполяции')
plt.legend()
plt.grid()
plt.show()
