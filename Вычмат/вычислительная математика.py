import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from sympy.codegen.rewriting import optimize
from scipy import optimize
from table import table

i = 0

x_sym = sp.symbols('x')

a = 0.4
b = 0.9
x_star = 0.52
count_points = 11

# f_expr можно менять в одном месте
# f_expr = x_sym**2 - sp.log(x_sym + 2, 10)
# f_expr = x_sym**2 - sp.log(x_sym, 10)
f_expr = x_sym**3 - sp.sin(x_sym)

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

print(f'R_2_star_num: {R_2_star_num}')
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

y0 = y_im1_sym
y1 = y_i_sym
y2 = y_ip1_sym

f01 = sp.simplify((y1 - y0) / (x_i - x_im1))   # f[x0,x1]
f12 = sp.simplify((y2 - y1) / (x_ip1 - x_i))   # f[x1,x2]



f012 = sp.simplify((f12 - f01) / (x_ip1 - x_im1))  # f[x0,x1,x2]
f012 = sp.N(f012, 30)
f12 = sp.N(f12, 30)
x_im1 = sp.N(x_im1, 30)
x_ip1 = sp.N(x_ip1, 30)
f01 = sp.N(f01, 30)
y0 = sp.N(y0, 30)
y1 = sp.N(y1, 30)

print("Таблица разделённых разностей для узлов:")
print(f"x_{i-1} = {x_im1}, f = {y0}")
print(f"x_{i} = {x_i}, f = {y1}")
print(f"x_{i+1} = {x_ip1}, f = {y2}")
print(f"f[x0,x1] = {f01}")
print(f"f[x1,x2] = {f12}")
print(f"f[x0,x1,x2] = {f012}")


# Многочлены Ньютона
x = x_sym

N1 = sp.simplify(y0 + f01 * (x - x_im1))
N2 = sp.simplify(y0 + f01 * (x - x_im1) + f012 * (x - x_im1) * (x - x_i))

N1_xstar_sym = sp.simplify(N1.subs(x_sym, x_star_sym))
N2_xstar_sym = sp.simplify(N2.subs(x_sym, x_star_sym))

N1_xstar_num = sp.N(N1_xstar_sym, 30)
N2_xstar_num = sp.N(N2_xstar_sym, 30)

print("Многочлены Ньютона:")
# print(f"N1(x) (симв.): {N1}")
# print(f"N2(x) (симв.): {N2}")
print(f"N1(x*): {N1_xstar_sym} ")
print(f"N2(x*): {N2_xstar_sym} ")



N1_vs_L1_diff = float(sp.N(N1_xstar_sym - L1, 30))
N2_vs_L2_diff = float(sp.N(N2_xstar_sym - L2, 30))

print("Сравнение Ньютона и Лангранджа в точке x*:")
print(f"N1(x*) - L1(x*):  {N1_vs_L1_diff})")
print(f"N2(x*) - L2(x*): {N2_vs_L2_diff})")



# # краткий итог в читабельном виде
# print("\nИтоги (численно, 30 значащих):")
# print(f"N1(x*): {N1_xstar_num}")
# print(f"L1(x*): {sp.N(L1,30)}")
# print(f"разность N1-L1: {N1_vs_L1_diff}")
#
# print(f"N2(x*): {N2_xstar_num}")
# print(f"L2(x*): {sp.N(L2,30)}")
# print(f"разность N2-L2: {N2_vs_L2_diff}")
#








# ------------------------
# Графики
# ------------------------
L1_poly = sp.simplify(
    y_i_sym * (x_sym - x_ip1) / (x_i - x_ip1) +
    y_i1_sym * (x_sym - x_i) / (x_ip1 - x_i)
)

L2_poly = sp.simplify(
    y_im1_sym * (x_sym - x_i) * (x_sym - x_ip1) / ((x_im1 - x_i) * (x_im1 - x_ip1)) +
    y_i_sym   * (x_sym - x_im1) * (x_sym - x_ip1) / ((x_i - x_im1) * (x_i - x_ip1)) +
    y_ip1_sym * (x_sym - x_im1) * (x_sym - x_i) / ((x_ip1 - x_im1) * (x_ip1 - x_i))
)

L1_plot = sp.lambdify(x_sym, L1_poly, 'numpy')
L2_plot = sp.lambdify(x_sym, L2_poly, 'numpy')

N1_plot = sp.lambdify(x_sym, N1, 'numpy')
N2_plot = sp.lambdify(x_sym, N2, 'numpy')

# 1. Функция
plt.figure()
plt.plot(x_fit, y_fit, color="red", label="f(x)")
plt.plot(x_points, y_points, "o", color="blue", label="табличные точки")
plt.axvline(x_star_num, color="gray", linestyle="--", alpha=0.5)
plt.title('Функция и табличные точки')
plt.legend()
plt.grid()



# # 2. Лагранж L1
# plt.figure()
# plt.plot(x_fit, L1_plot(x_fit), label="Лагранж L1", linestyle="--")
# plt.title('Интерполяция Лагранжа L1')
# plt.legend()
# plt.grid()


# 3. Лагранж L2
plt.figure()
plt.plot(x_fit, L2_plot(x_fit), label="Лагранж L2", linestyle="-.")
plt.title('Интерполяция Лагранжа L2')
plt.legend()
plt.grid()


# # 4. Ньютон N1
# plt.figure()
# plt.plot(x_fit, N1_plot(x_fit), label="Ньютон N1", linestyle="--", alpha=0.7)
# plt.title('Интерполяция Ньютона N1')
# plt.legend()
# plt.grid()


# 5. Ньютон N2
plt.figure()
plt.plot(x_fit, N2_plot(x_fit), label="Ньютон N2", linestyle="-.", alpha=0.7)
plt.title('Интерполяция Ньютона N2')
plt.legend()
plt.grid()

#6
#
# x_big = np.linspace(0, 3, 500)  # от 0.1, чтобы избежать log(0)
# plt.figure()
# plt.plot(x_big, f(x_big), color="red", label="f(x) = x^2 + ln(x)")
# plt.title("Функция f(x) в большом масштабе")
# plt.legend()
# plt.grid()

# показываем всё сразу
plt.show()