import sympy as sp
x_sym = sp.symbols('x')

table = [
    {
        "f_expr": x_sym**2 + sp.log(x_sym),
        "a": 0.4,
        "b": 0.9,
        "x_star": 0.52,
        "x_pow2": 0.42,
        "x_pow3": 0.87,
        "x_pow4": 0.67
    },
    {
        "f_expr": x_sym**2 - sp.log(x_sym + 2, 10),
        "a": 0.5,
        "b": 1,
        "x_star": 0.53,
        "x_pow2": 0.52,
        "x_pow3": 0.97,
        "x_pow4": 0.73
    },
]
