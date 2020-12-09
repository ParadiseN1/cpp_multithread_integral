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

double two_sum (double &t, double a, double b) {
    double s = a+b;
    double bs = s-a;
    double as = s-bs;
    t = (b-bs) + (a-as);
    return s;
}
double fast_two_sum (double &t, double a, double b) {
    double s = a+b;
    t = b-(s-a);
    return s;
}

//double sum_kahan (const dvect_t &X) {
//    double s=0.0, c=0.0;
//    for (auto x: X) {
//        double y = x + c;
//        s = fast_two_sum (c, s, y);
//    }
//    return s;
//}




double de_zong(double x1, double x2){
    double res = 0.002;
    double sum = 0.0;
    double den = 0.0;
    //double e;
    double c = 0.0;
    for (int i = -2; i <= 2; ++i){
        for (int j = -2; j <= 2; ++j){
            den = 5 * (i + 2) + j + 3 + pow(x1 - 16 * j, 6) + pow(x2 - 16 * i, 6);
            double y = (1 / den) + c;
            sum = fast_two_sum (c, sum, y);

            //sum += two_sum(e, sum, 1 / den);
            //c+=e;

            //sum += 1/den;
        }
    }
    //sum += e;

    return pow(res + sum, -1);
}
double integral(double(*f)(double x1, double x2),double a, double b,double a1, double b1, double dx, double prev = 0) {
    double area = 0;
    for (double x = a; x < b; x+=dx){
        for (double y = a1; y < b1; y+=dx){
            area += f(x, y) * dx * dx;
        }
    }
    return area;
}

//double __fastcall sum_rump (const dvect_t &X) {
//    double s=0.0, c=0.0, e;
//    for (double x: X) {
//        s = two_sum (e, s, x);
//        c += e;
//    }
//    return s+c;
//}
double tintegral(double(*f)(double x1, double x2),double a, double b,double a1, double b1, double dx, double *res, double prev=0) {
    double area = 0;
    std::mutex mutex;
    double c=0.0, e;
    if (prev == 0)
        for (double x = a; x < b; x+=dx){
            for (double y = a1; y < b1; y+=dx){
               area = two_sum(e, area, f(x, y) * dx * dx);
               c += e;

            }
        }
    else
        for (double x = a+dx; x < b; x+=2*dx){
            for (double y = a1; y < b1; y+=dx){
                area = two_sum(e, area, f(x, y) * dx * dx);
                c += e;
            }
        }
    area += c;
    mutex.lock();
    *res += area;
    mutex.unlock();
    return area;
}

double thread_integral(double(*f)(double x1, double x2),double a, double b, double a1, double b1, double dx, int t, double prev =0){
    std::thread threads[t];
    double res = 0;
//    if (prev == 0)
    for (int i = 0; i < t; ++i){
        threads[i] = std::thread(tintegral, f, a + (((b-a) / t) * i), b - (((b-a) / t) * (t - i - 1)), \
                                          a1, b1, dx, &res, prev);
    }
    for (int i = 0; i < t; ++i){
        threads[i].join();
    }
    return (res + prev/2);
}

int main()
{
    int lower = -50;
    int upper = 50;
    double dx = abs(lower - upper) / 500.0;
    std::cout.precision(9);
    auto  before = get_current_time_fenced();
    int i = 1;
    double val = 0;
    double prev;
    while (i){
        prev = val;
        val = thread_integral(de_zong, lower, upper,lower, upper, dx, 4, prev);
        std::cout << i << " itetation" << ": " << val << ", dx:" << dx << ", error:" << abs(prev - val) / prev << "\n";
        if (abs(prev - val) / prev < 0.0001)
            break;
        dx /= 2.0;
        ++i;
    }
    auto time_to_calculate = get_current_time_fenced() - before;
    std::cout << to_s(time_to_calculate) << "s\n";
    return 0;
}
