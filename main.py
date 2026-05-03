import matplotlib.pyplot as plt

N = [100, 400, 1600, 6400, 25600, 102400, 409600, 1638400]

iter_time = [
    2.264977e-03,
    6.947994e-03,
    2.015615e-02,
    1.313379e-01,
    1.036980e+00,
    6.190235e+00,
    5.407916e+01,
    4.415511e+02
]

total_time = [
    6.991141e-03,
    2.853208e-02,
    5.333059e-02,
    2.644674e-01,
    1.643992e+00,
    8.633004e+00,
    6.388752e+01,
	5.166569e+02
]

plt.figure()
plt.loglog(N, iter_time, marker='o')
plt.xlabel("Размер системы N")
plt.ylabel("Число итераций")
plt.grid(True)
plt.title("Число итераций от размера системы")
plt.show()

plt.figure()
plt.loglog(N, total_time, marker='o')
plt.xlabel("Размер системы N")
plt.ylabel("Время построения решения, сек")
plt.grid(True)
plt.title("Время построения решения от размера системы")
plt.show()


import matplotlib.pyplot as plt

h = [
    0.0909091,
    4.761905e-02,
    2.439024e-02,
    1.234568e-02,
    6.211180e-03,
    3.115265e-03,
    1.560062e-03,
]

c_norm_error = [
    3.770941e-03,
    1.070875e-03,
    2.823413e-04,
    7.240855e-05,
    1.833261e-05,
    4.612189e-06,
    1.157134e-06,
]

l2_norm_error = [
    1.632250e-03,
    4.479220e-04,
    1.175001e-04,
    3.010391e-05,
    7.619714e-06,
    1.916910e-06,
    4.809401e-07,
]

# Эталонная линия O(h^2), масштабируем под первую точку C-нормы
oh2_c = [
    c_norm_error[0] * (hi / h[0])**2
    for hi in h
]

# Эталонная линия O(h^2), масштабируем под первую точку L2-нормы
oh2_l2 = [
    l2_norm_error[0] * (hi / h[0])**2
    for hi in h
]

plt.figure()
plt.loglog(h, c_norm_error, marker='o', label='C-норма')
plt.loglog(h, l2_norm_error, marker='o', label='L2-норма')
plt.loglog(h, oh2_c, '--', label='O(h^2) для C-нормы')
plt.loglog(h, oh2_l2, '--', label='O(h^2) для L2-нормы')

plt.xlabel('Шаг сетки h')
plt.ylabel('Погрешность')
plt.title('Погрешность решения от шага сетки')
plt.grid(True, which='both')
plt.legend()
plt.show()
