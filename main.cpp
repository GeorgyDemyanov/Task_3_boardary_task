#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <fstream>


std::vector<double> multiple(std::vector<std::vector<double>> a, std::vector<double> b) {
    std::vector<double> res(a[0].size());
    for (int i = 0; i < a[0].size(); i++) {
        for (int j = 0; j < b.size(); j++) {
            res[i] += a[j][i] * b[j];
        }
    }
    return res;
}

std::vector<double> operator*(double a, std::vector<double> b) {
    for (double &i: b) {
        i *= a;
    }
    return b;
}

std::vector<double> operator+(std::vector<double> a, std::vector<double> b) {
    for (int i = 0; i < b.size(); i++) {
        b[i] += a[i];
    }
    return b;
}

std::ostream &operator<<(std::ostream &os, const std::vector<double> &a) {
    for (double i: a) {
        os << i << std::endl;
    }
    return os;
}

double difference(std::vector<std::vector<double>> a, std::vector<std::vector<double>> b) {
    double max = 0;
    for (int i = 0; i < a.size(); ++i) {
        for (int j = 0; j < a.size(); ++j) {
            if (max < std::abs(a[i][j] - b[i][j])) {
                max = abs(a[i][j] - b[i][j]);
            }
        }
    }
    return max;
}

std::vector<std::vector<double>>
k(std::vector<double> sol, double time, double h, std::vector<std::vector<double>> a, std::vector<double> c,
  std::function<std::vector<double>(double, std::vector<double>)> func) {
    std::vector<std::vector<double>> k(4);
    for (int i = 0; i < 4; ++i) {
        k[i] = std::vector<double>(2);
    }
    std::vector<double> ki(4);
    double diff = 1;
    while (diff > 1e-6) {
        std::vector<std::vector<double>> k_prev = k;
        for (int i = 0; i < k.size(); i++) {
            ki = func(time + c[i] * h, sol + h * multiple(k, a[i]));
            k[i] = ki;
        }
        diff = difference(k, k_prev);
    }
    return k;
}

std::vector<double>
step(double time, double h, std::vector<double> state, std::vector<std::vector<double>> a, std::vector<double> b,
     std::vector<double> c, std::function<std::vector<double>(double, std::vector<double>)> func) {
    return state + h * multiple(k(state, time, h, a, c, func), b);
}

std::vector<std::pair<double, double>> solution(double start, double end, double h, std::vector<double> state,
                                                std::function<std::vector<double>(double, std::vector<double>)> func) {
    double current_time = start;
    std::vector<std::pair<double, double>> result;
    std::vector<std::vector<double>> a = {{0,   0,   0, 0},
                                          {0.5, 0,   0, 0},
                                          {0,   0.5, 0, 0},
                                          {0,   0,   1, 0}};
    std::vector<double> b = {1. / 6, 1. / 3, 1. / 3, 1. / 6};
    std::vector<double> c = {0, 0.5, 0.5, 1};
    result.push_back({current_time, state[0]});
    while (current_time < end) {
        state = step(current_time, h, state, a, b, c, func);
        current_time += h;
        result.push_back({current_time, state[0]});

    }
    return result;
}

std::vector<double> f(double time, std::vector<double> sol) {
    return {sol[1], 2 - 6 * time + 2 * pow(time, 3) + (time * time - 3) * exp(time) * sin(time) * (1 + cos(time)) +
                    cos(time) * (exp(time) + (time * time - 1) + pow(time, 4) - 3 * time * time) -
                    (time * time - 3) * sol[1] + (time * time - 3) * cos(time) * sol[0]};
}

int main() {
    double a1 = M_PI * M_PI - 10;
    double a2 = M_PI * M_PI + 1;
    double m = (a1 + a2) / 2;
    double Rboard = M_PI*M_PI;

    std::vector<std::pair<double, double>> sol = solution(0, 3.1415, 0.001, {0, m}, f);
    while (std::fabs(sol[sol.size() - 1].second - Rboard) > 1e-6) {
        if (sol[sol.size() - 1].second > Rboard) {
            a2 = m;
            m = (a1 + a2) / 2;
        } else {
            if (sol[sol.size() - 1].second < Rboard) {
                a1 = m;
                m = (a1 + a2) / 2;
            }
        }

        sol = solution(0, 3.1415, 0.001, {0, m}, f);
    }
    std::cout<<m<<std::endl;
    std::ofstream fout("C:\\Users\\Demia\\CLionProjects\\Task 3 boardary-task\\data.txt");
    for (int i = 0; i < sol.size(); ++i) {
        std::cout << sol[i].first << " " << sol[i].second << std::endl;
        fout << sol[i].first << " " << sol[i].second << std::endl;
    }
    fout.close();
    return 0;
}
