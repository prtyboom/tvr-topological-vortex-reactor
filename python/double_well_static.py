import numpy as np
import matplotlib.pyplot as plt

try:
    from scipy.linalg import eigh
    have_scipy = True
except ImportError:
    have_scipy = False


def build_double_well_hamiltonian(
    N=400, y_max=3.0, epsilon=0.2
):
    """
    Строит дискретизированный гамильтониан для безразмерного
    оператора:
        H = -epsilon^2 d^2/dy^2 + (y^2 - 1)**2
    на отрезке [-y_max, y_max] с N точками.
    """
    y = np.linspace(-y_max, y_max, N)
    dy = y[1] - y[0]

    # Вторая производная (1D Лаплас, центральная разность O(dy^2))
    diag = np.full(N, -2.0)
    off = np.ones(N - 1)
    lap = (np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)) / (dy**2)

    # Потенциал двойной ямы (y^2 - 1)^2
    V = (y**2 - 1.0)**2
    H = - (epsilon**2) * lap + np.diag(V)
    return H, y


def compute_spectrum(H, k=6):
    """
    Вычисляет k наименьших собственных значений H.
    Если есть scipy.linalg.eigh — используем её,
    иначе fallback на numpy.linalg.eigh (медленнее, но работает).
    """
    if have_scipy:
        vals, vecs = eigh(H)
    else:
        vals, vecs = np.linalg.eigh(H)
    return vals[:k], vecs[:, :k]


def main():
    H, y = build_double_well_hamiltonian(N=400, y_max=3.0, epsilon=0.2)
    E, psi = compute_spectrum(H, k=4)

    print("Lowest energy levels (dimensionless):")
    for i, Ei in enumerate(E):
        print(f"  E[{i}] = {Ei:.6f}")

    # Потенциал и первые 2 волновые функции (масштабированные)
    V = (y**2 - 1.0)**2
    plt.figure(figsize=(8, 5))
    plt.plot(y, V, "k-", label="V(y) = (y^2 - 1)^2")

    scale = 0.5  # масштаб, чтобы волновые функции было видно на фоне потенциала
    for i in range(min(2, psi.shape[1])):
        psi_i = psi[:, i]
        # нормируем волновую функцию
        norm = np.trapz(np.abs(psi_i)**2, y)
        psi_i = psi_i / np.sqrt(norm)
        plt.plot(y, scale * psi_i.real + E[i], label=f"psi_{i}(y) + E[{i}]")

    plt.xlabel("y")
    plt.ylabel("Energy / Potential")
    plt.title("Double-well potential and lowest eigenstates (static case)")
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()