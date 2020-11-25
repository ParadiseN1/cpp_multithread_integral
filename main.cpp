#include <iostream>
#include <cmath>
#include <thread>
#include <vector>
#include <mutex>

inline std::chrono::steady_clock::time_point get_current_time_fenced() {
    assert(std::chrono::steady_clock::is_steady &&
                   "Timer should be steady (monotonic).");
    std::atomic_thread_fence(std::memory_order_seq_cst);
    auto res_time = std::chrono::steady_clock::now();
    std::atomic_thread_fence(std::memory_order_seq_cst);
    return res_time;
}

template<class D>
inline long long to_s(const D& d)
{
    return std::chrono::duration_cast<std::chrono::microseconds>(d).count() / 1000000;
}

double de_zong(double x1, double x2){
    double res = 0.002;
    double sum = 0;
    double den = 0;
    for (int i = -2; i <= 2; ++i){
        for (int j = -2; j <= 2; ++j){
            den = 5 * (i + 2) + j + 3 + pow(x1 - 16 * j, 6) + pow(x2 - 16 * i, 6);
            sum += 1 / den;
        }
    }
    return pow(res + sum, -1);
}
double integral(double(*f)(double x1, double x2),double a, double b,double a1, double b1, double dx) {
    double area = 0;
    for (double x = a; x < b; x+=dx){
        for (double y = a1; y < b1; y+=dx){
            area += f(x, y) * dx * dx;
        }
    }
    return area;
}

double tintegral(double(*f)(double x1, double x2),double a, double b,double a1, double b1, double dx, double *res) {
    double area = 0;
    std::mutex mutex;
    for (double x = a; x < b; x+=dx){
        for (double y = a1; y < b1; y+=dx){
            area += f(x, y) * dx * dx;
        }
    }
    mutex.lock();
    *res += area;
    mutex.unlock();
    return area;
}

double thread_integral(double(*f)(double x1, double x2),double a, double b, double a1, double b1, double dx, int t){
    std::thread threads[t];
    double res = 0;
    for (int i = 0; i < t; ++i){
        threads[i] = std::thread(tintegral, f, a + (((b-a) / t) * i), b - (((b-a) / t) * (t - i - 1)), \
                                          a1, b1, dx, &res);
    }
    for (int i = 0; i < t; ++i){
        threads[i].join();
    }
    return (res);
}

int main()
{
    double dx = 0.1;
    std::cout.precision(9);
    auto  before = get_current_time_fenced();
    for (int i = 1; i < 4; ++i){

        std::cout << "it" << i << ": " << thread_integral(de_zong, -50, 50,-50, 50, dx/(2*i), 4) << "\n";
    }
    auto time_to_calculate = get_current_time_fenced() - before;
    std::cout << to_s(time_to_calculate) << "s\n";
    before = get_current_time_fenced();
    for (int i = 1; i < 4; ++i){
        std::cout << "it" << i << ": " << integral(de_zong, -50, 50,-50, 50, dx/(2*i)) << "\n";
    }
    time_to_calculate = get_current_time_fenced() - before;
    std::cout << to_s(time_to_calculate) << "s\n";
    return 0;
}
