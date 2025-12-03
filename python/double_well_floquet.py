import numpy as np
import matplotlib.pyplot as plt


def initial_wavefunction(y, sigma=0.3):
    """
    Начальное состояние: гауссиан в левой яме (около y = -1).
    Можно поменять на правую яму или суперпозицию.
    """
    psi0 = np.exp(-0.5 * ((y + 1.0) / sigma) ** 2)
    # Нормировка
    dy = y[1] - y[0]
    norm = np.sqrt(np.trapezoid(np.abs(psi0) ** 2, dx=dy))
    return psi0 / norm


def build_grids(N=1024, y_max=4.0):
    """
    Строим периодический сеточный интервал [-y_max, y_max) с N точками
    и соответствующий спектр k для FFT.
    """
    # Периодическая сетка без повторения крайней точки
    y = np.linspace(-y_max, y_max, N, endpoint=False)
    dy = y[1] - y[0]
    # Волновые числа для FFT
    k = 2.0 * np.pi * np.fft.fftfreq(N, d=dy)
    return y, k, dy


def potential(y, tau, mu):
    """
    Безразмерный потенциал:
        V(y, tau) = (y^2 - 1)^2 + 4 mu cos(tau) y
    """
    return (y**2 - 1.0) ** 2 + 4.0 * mu * np.cos(tau) * y


def evolve_split_operator(
    epsilon=0.2,
    mu=0.3,
    N=1024,
    y_max=4.0,
    n_periods=10,
    steps_per_period=200,
):
    """
    Эволюция волновой функции по времени с помощью схемы split-operator (Strang splitting):
        U(Δτ) ≈ e^{-i Δτ V/2ε} e^{-i Δτ T/ε} e^{-i Δτ V/2ε}
    где
        T = -ε^2 d^2/dy^2
    в импульсном пространстве даёт множитель exp(-i ε k^2 Δτ).

    Возвращает:
        taus  - массив значений τ
        P_L   - вероятность в левой половине (y<0)
        P_R   - вероятность в правой половине (y>0)
    """
    # Сетки
    y, k, dy = build_grids(N=N, y_max=y_max)

    # Временной шаг и массив τ
    T_period = 2.0 * np.pi  # период по τ: 0..2π
    dtau = T_period / steps_per_period
    n_steps = n_periods * steps_per_period
    taus = np.linspace(0.0, n_periods * T_period, n_steps + 1)

    # Фиксированная кинетическая часть в k-пространстве:
    # из уравнения i ε ∂_τ ψ = (-ε^2 d^2/dy^2 + V) ψ
    # в k-пространстве T = ε^2 k^2 => фактор exp(-i ε k^2 Δτ)
    kinetic_phase = np.exp(-1j * epsilon * (k**2) * dtau)

    # Начальное состояние
    psi = initial_wavefunction(y)  # комплексный массив
    # для аккуратности приведем тип
    psi = psi.astype(np.complex128)

    # Массивы для вероятностей
    P_L = np.zeros_like(taus)
    P_R = np.zeros_like(taus)

    # Индексы для левой и правой половины
    left_mask = y < 0.0
    right_mask = y > 0.0

    # Функция для вычисления вероятностей
    def compute_probabilities(psi_state):
        prob_density = np.abs(psi_state) ** 2
        pL = np.trapezoid(prob_density[left_mask], dx=dy)
        pR = np.trapezoid(prob_density[right_mask], dx=dy)
        return pL, pR

    # Начальная вероятность
    P_L[0], P_R[0] = compute_probabilities(psi)

    # Основной цикл по времени
    for n in range(n_steps):
        tau_mid = taus[n] + 0.5 * dtau

        # Потенциал на середине шага
        V_mid = potential(y, tau_mid, mu)

        # Полушаг по потенциалу в координатном пространстве
        phase_V_half = np.exp(-0.5j * dtau * V_mid / epsilon)
        psi *= phase_V_half

        # Полный шаг по кинетике в k-пространстве
        psi_k = np.fft.fft(psi)
        psi_k *= kinetic_phase
        psi = np.fft.ifft(psi_k)

        # Ещё полушаг по потенциалу
        psi *= phase_V_half

        # Нормировка (накапливаемые численные ошибки можно подрезать)
        norm = np.sqrt(np.trapezoid(np.abs(psi) ** 2, dx=dy))
        psi /= norm

        # Сохраняем вероятности в левой и правой половине
        P_L[n + 1], P_R[n + 1] = compute_probabilities(psi)

    return taus, P_L, P_R, y, psi


def main():
    # Параметры можно подправлять
    epsilon = 0.2
    mu = 0.3
    N = 1024
    y_max = 4.0
    n_periods = 10
    steps_per_period = 200

    taus, P_L, P_R, y, psi_final = evolve_split_operator(
        epsilon=epsilon,
        mu=mu,
        N=N,
        y_max=y_max,
        n_periods=n_periods,
        steps_per_period=steps_per_period,
    )

    print(f"Simulation finished: epsilon={epsilon}, mu={mu}")
    print(f"Total time in tau units: {taus[-1]:.3f}")
    print(f"Final probabilities: P_L={P_L[-1]:.4f}, P_R={P_R[-1]:.4f}")

    # График вероятностей во времени
    plt.figure(figsize=(8, 5))
    plt.plot(taus, P_L, label="P_L (y<0)")
    plt.plot(taus, P_R, label="P_R (y>0)")
    plt.xlabel(r"$\tau$")
    plt.ylabel("Probability")
    plt.title("Floquet-driven double-well: left/right well occupation")
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()