#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

// Un hilo ejecuta el cálculo secuencial; varios hilos utilizan OpenMP.
constexpr double ALFA = 1.0;
constexpr double LX = 1.0;
constexpr double LY = 1.0;
constexpr double PI = 3.14159265358979323846;

int main(int argc, char* argv[]) {
    try {
        const std::string condicion = argc > 1 ? argv[1] : "gauss";
        const std::string frontera = argc > 2 ? argv[2] : "fija";
        const int nx = argc > 3 ? std::stoi(argv[3]) : 81;
        const int ny = argc > 4 ? std::stoi(argv[4]) : 81;
        const double tiempo = argc > 5 ? std::stod(argv[5]) : 0.05;
        const int hilos = argc > 6 ? std::stoi(argv[6]) : 1;

        if (condicion != "gauss" && condicion != "doble" && condicion != "seno")
            throw std::invalid_argument("Condición inicial desconocida");
        if (frontera != "fija" && frontera != "aislada")
            throw std::invalid_argument("Frontera desconocida");
        if (nx < 3 || ny < 3 || !std::isfinite(tiempo) || tiempo <= 0 || hilos < 1)
            throw std::invalid_argument("Parámetros fuera de rango");

        const double dx = LX / (nx - 1);
        const double dy = LY / (ny - 1);
        const double dt_max = 1 / (2 * ALFA * (1/(dx*dx) + 1/(dy*dy)));
        const double pasos_reales = std::ceil(tiempo/(0.9*dt_max));
        if (pasos_reales > std::numeric_limits<int>::max())
            throw std::invalid_argument("Demasiados pasos temporales");
        const int pasos = static_cast<int>(pasos_reales);
        const double dt = tiempo / pasos;
        const double rx = ALFA * dt / (dx*dx);
        const double ry = ALFA * dt / (dy*dy);
        const bool fija = (frontera == "fija");

        std::vector<double> u(static_cast<size_t>(nx)*ny, 0.0);
        std::vector<double> v(u.size(), 0.0);
        for (int j = 0; j < ny; ++j) {
            const double y = j*dy;
            for (int i = 0; i < nx; ++i) {
                const double x = i*dx;
                const int k = j*nx + i;
                if (condicion == "seno") {
                    u[k] = 100*std::sin(PI*x/LX)*std::sin(PI*y/LY);
                } else {
                    const double a = condicion == "doble" ? 0.3 : 0.5;
                    const double b = condicion == "doble" ? 0.4 : 0.5;
                    u[k] = 100*std::exp(-((x-a)*(x-a)+(y-b)*(y-b))/0.02);
                    if (condicion == "doble")
                        u[k] += 70*std::exp(-((x-0.7)*(x-0.7)+(y-0.6)*(y-0.6))/0.02);
                }
                if (fija && (i == 0 || j == 0 || i == nx-1 || j == ny-1))
                    u[k] = 0;
            }
        }

        // Cada fila se escribe en v utilizando solamente la malla anterior u.
        auto actualizar_fila = [&](int j) {
            const int abajo = j == 0 ? j : j-1;
            const int arriba = j == ny-1 ? j : j+1;
            for (int i = 0; i < nx; ++i) {
                const int k = j*nx+i;
                if (fija && (i == 0 || i == nx-1 || j == 0 || j == ny-1)) {
                    v[k] = 0;
                    continue;
                }
                const int izq = i == 0 ? i : i-1;
                const int der = i == nx-1 ? i : i+1;
                v[k] = u[k] + rx*(u[j*nx+izq] - 2*u[k] + u[j*nx+der])
                              + ry*(u[abajo*nx+i] - 2*u[k] + u[arriba*nx+i]);
            }
        };

        const auto inicio = std::chrono::steady_clock::now();
        if (hilos == 1) {
            for (int n = 0; n < pasos; ++n) {
                for (int j = 0; j < ny; ++j)
                    actualizar_fila(j);
                u.swap(v);
            }
        } else {
            // El equipo de hilos se crea una sola vez.
            #pragma omp parallel num_threads(hilos)
            {
                for (int n = 0; n < pasos; ++n) {
                    #pragma omp for schedule(static)
                    for (int j = 0; j < ny; ++j)
                        actualizar_fila(j);

                    // Las barreras evitan leer una malla mientras se escribe.
                    #pragma omp single
                    u.swap(v);
                }
            }
        }
        const double segundos = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - inicio).count();
        const double promedio = std::accumulate(u.begin(), u.end(), 0.0)/u.size();
        const double centro = u[(ny/2)*nx + nx/2];

        std::cout << std::setprecision(12)
                  << "Hilos: " << hilos << "; pasos: " << pasos
                  << "; dt: " << dt << "; tiempo: " << segundos << " s\n"
                  << "Temperatura media: " << promedio
                  << "; centro: " << centro << '\n';
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
