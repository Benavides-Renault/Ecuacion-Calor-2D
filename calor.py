"""Difusión térmica en una placa rectangular mediante diferencias finitas."""

import argparse
import math
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

ALFA = 1.0
LX = 1.0
LY = 1.0


def inicializar(nx, ny, condicion, frontera):
    """Genera una temperatura inicial: gauss, doble o seno."""
    x = np.linspace(0, LX, nx)
    y = np.linspace(0, LY, ny)
    X, Y = np.meshgrid(x, y)

    if condicion == "gauss":
        u = 100 * np.exp(-((X - 0.5)**2 + (Y - 0.5)**2) / 0.02)
    elif condicion == "doble":
        u = 100 * np.exp(-((X - 0.3)**2 + (Y - 0.4)**2) / 0.02)
        u += 70 * np.exp(-((X - 0.7)**2 + (Y - 0.6)**2) / 0.02)
    elif condicion == "seno":
        u = 100 * np.sin(np.pi * X / LX) * np.sin(np.pi * Y / LY)
    else:
        raise ValueError("Condición inicial: gauss, doble o seno")

    if frontera == "fija":
        u[[0, -1], :] = 0
        u[:, [0, -1]] = 0
    return u


def paso(u, rx, ry, frontera, salida=None):
    """Avanza un paso temporal; permite reutilizar la matriz de salida."""
    if salida is None:
        salida = np.empty_like(u)

    if frontera == "fija":
        centro = u[1:-1, 1:-1]
        salida[1:-1, 1:-1] = (
            centro
            + rx * (u[1:-1, 2:] - 2 * centro + u[1:-1, :-2])
            + ry * (u[2:, 1:-1] - 2 * centro + u[:-2, 1:-1])
        )
        salida[[0, -1], :] = 0
        salida[:, [0, -1]] = 0
    elif frontera == "aislada":
        p = np.pad(u, 1, mode="edge")
        salida[:] = (
            u
            + rx * (p[1:-1, 2:] - 2 * u + p[1:-1, :-2])
            + ry * (p[2:, 1:-1] - 2 * u + p[:-2, 1:-1])
        )
    else:
        raise ValueError("Frontera: fija o aislada")
    return salida


def simular(condicion="gauss", frontera="fija", nx=81, ny=81, tiempo=0.05):
    """Devuelve temperaturas en cuatro instantes y parámetros de integración."""
    if frontera not in ("fija", "aislada"):
        raise ValueError("Frontera: fija o aislada")
    if (nx < 3 or ny < 3 or not math.isfinite(tiempo) or tiempo <= 0
            or not math.isfinite(ALFA) or ALFA <= 0 or LX <= 0 or LY <= 0):
        raise ValueError("Malla, tiempo y difusividad deben ser positivos")

    dx = LX / (nx - 1)
    dy = LY / (ny - 1)
    dt_max = 1 / (2 * ALFA * (1/dx**2 + 1/dy**2))
    pasos = math.ceil(tiempo / (0.9 * dt_max))
    dt = tiempo / pasos
    rx = ALFA * dt / dx**2
    ry = ALFA * dt / dy**2

    u = inicializar(nx, ny, condicion, frontera)
    auxiliar = np.empty_like(u)
    instantes = {0, pasos//3, 2*pasos//3, pasos}
    cuadros = [(0.0, u.copy())]
    for n in range(1, pasos + 1):
        paso(u, rx, ry, frontera, salida=auxiliar)
        u, auxiliar = auxiliar, u
        if n in instantes:
            cuadros.append((n * dt, u.copy()))
    return cuadros, pasos, dt


def graficar(cuadros, destino="docs/calor_2d.png"):
    """Guarda mapas 2D de la evolución térmica."""
    fig, ejes = plt.subplots(1, len(cuadros), figsize=(12, 3),
                            constrained_layout=True)
    for eje, (t, u) in zip(np.atleast_1d(ejes), cuadros):
        imagen = eje.imshow(u, origin="lower", extent=[0, LX, 0, LY],
                            cmap="inferno", vmin=0, vmax=100)
        eje.set_title(f"t = {t:.4f}")
        eje.set_xlabel("x")
        eje.set_ylabel("y")
    fig.colorbar(imagen, ax=ejes, label="Temperatura")
    Path(destino).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destino, dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Ecuación de calor 2D")
    parser.add_argument("condicion", nargs="?", default="gauss")
    parser.add_argument("frontera", nargs="?", default="fija")
    parser.add_argument("nx", nargs="?", type=int, default=81)
    parser.add_argument("ny", nargs="?", type=int, default=81)
    parser.add_argument("tiempo", nargs="?", type=float, default=0.05)
    args = parser.parse_args()

    cuadros, pasos, dt = simular(args.condicion, args.frontera,
                                args.nx, args.ny, args.tiempo)
    graficar(cuadros)
    final = cuadros[-1][1]
    print(f"Pasos: {pasos}; dt: {dt:.8g}")
    print(f"Temperatura media: {final.mean():.10f}")
    print(f"Temperatura central: {final[args.ny//2, args.nx//2]:.10f}")
    print("Gráfico: docs/calor_2d.png")
