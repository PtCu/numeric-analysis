#include <iostream>
#include "math.hpp"
#include <fstream>
#include <time.h>
#include <sys/time.h>

void solve(const Matrix<double, Dynamic, Dynamic> &m, const std::string &name)
{
    //幂法
    std::ofstream f("../ans.txt", std::ios::app); //输出文件
    double pld;
    std::vector<double> env;
    struct timeval tv;
    long long start, end;
    gettimeofday(&tv, NULL);
    start = tv.tv_sec * 1000 + tv.tv_usec / 1000;
    power_eng(m, pld, env);
    gettimeofday(&tv, NULL);
    end = tv.tv_sec * 1000 + tv.tv_usec / 1000;
    f << "幂法" << std::endl;
    f << name << "的最大特征值   " << pld << std::endl;
    f << name << "的对应的特征向量" << std::endl;
    VectorXd x;
    x.resize(env.size());
    for (size_t i = 0; i < env.size(); ++i)
    {
        f << env[i] << " ";
        x(i) = env[i];
    }
    f << "运行时间" << end - start << "ms" << std::endl;
    f << "误差" << (m * x - pld * x).squaredNorm() / x.squaredNorm() << std::endl;

    f << std::endl;

    //jacobi
    env.clear();
    gettimeofday(&tv, NULL);
    start = tv.tv_sec * 1000 + tv.tv_usec / 1000;
    jacobi_eng(m, env);
    gettimeofday(&tv, NULL);
    end = tv.tv_sec * 1000 + tv.tv_usec / 1000;
    f << "jacobi" << std::endl;
    f << name << "的特征值" << std::endl;

    sort(env.begin(), env.end());
    for (size_t i = 0; i < env.size(); ++i)
    {
        f << env[i] << " ";
    }
    f << "运行时间" << end - start << "ms" << std::endl;
    f << std::endl;

    //移动qr
    std::vector<std::complex<double>> env1;
    gettimeofday(&tv, NULL);
    start = tv.tv_sec * 1000 + tv.tv_usec / 1000;
    shifted_qr(m, env1);
    gettimeofday(&tv, NULL);
    end = tv.tv_sec * 1000 + tv.tv_usec / 1000;
    f << "移动QR" << std::endl;
    f << name << "的特征值" << std::endl;
    sort(env1.begin(), env1.end(), [](const std::complex<double> &a, const std::complex<double> &b) { return a.real() < b.real(); });
    for (size_t i = 0; i < env1.size(); ++i)
    {
        f << env1[i] << " ";
    }
    f << "运行时间" << end - start << "ms" << std::endl;
    f << std::endl;
    env1.clear();

}

int main()
{

    Matrix<double, 8, 8> A;
    A << 611, 196, -192, 407, -8, -52, -49, 29,
        196, 899, 113, -192, -71, -43, -8, -44,
        -192, 113, 899, 196, 61, 49, 8, 52,
        407, -192, 196, 611, 8, 44, 59, -23,
        -8, -71, 61, 8, 411, -599, 208, 208,
        -52, -43, 49, 44, -599, 411, 208, 208,
        -49, -8, 8, 59, 208, 208, 99, -911,
        29, -44, 52, -23, 208, 208, -911, 99;

    Matrix<double, 10, 10> B;

    B << 1 / 1.0, 1 / 2.0, 1 / 3.0, 1 / 4.0, 1 / 5.0, 1 / 6.0, 1 / 7.0, 1 / 8.0, 1 / 9.0, 1 / 10.0,
        1 / 2.0, 1 / 3.0, 1 / 4.0, 1 / 5.0, 1 / 6.0, 1 / 7.0, 1 / 8.0, 1 / 9.0, 1 / 10.0, 1 / 11.0,
        1 / 3.0, 1 / 4.0, 1 / 5.0, 1 / 6.0, 1 / 7.0, 1 / 8.0, 1 / 9.0, 1 / 10.0, 1 / 11.0, 1 / 12.0,
        1 / 4.0, 1 / 5.0, 1 / 6.0, 1 / 7.0, 1 / 8.0, 1 / 9.0, 1 / 10.0, 1 / 11.0, 1 / 12.0, 1 / 13.0,
        1 / 5.0, 1 / 6.0, 1 / 7.0, 1 / 8.0, 1 / 9.0, 1 / 10.0, 1 / 11.0, 1 / 12.0, 1 / 13.0, 1 / 14.0,
        1 / 6.0, 1 / 7.0, 1 / 8.0, 1 / 9.0, 1 / 10.0, 1 / 11.0, 1 / 12.0, 1 / 13.0, 1 / 14.0, 1 / 15.0,
        1 / 7.0, 1 / 8.0, 1 / 9.0, 1 / 10.0, 1 / 11.0, 1 / 12.0, 1 / 13.0, 1 / 14.0, 1 / 15.0, 1 / 16.0,
        1 / 8.0, 1 / 9.0, 1 / 10.0, 1 / 11.0, 1 / 12.0, 1 / 13.0, 1 / 14.0, 1 / 15.0, 1 / 16.0, 1 / 17.0,
        1 / 9.0, 1 / 10.0, 1 / 11.0, 1 / 12.0, 1 / 13.0, 1 / 14.0, 1 / 15.0, 1 / 16.0, 1 / 17.0, 1 / 18.0,
        1 / 10.0, 1 / 11.0, 1 / 12.0, 1 / 13.0, 1 / 14.0, 1 / 15.0, 1 / 16.0, 1 / 17.0, 1 / 18.0, 1 / 19.0;

    Matrix<double, 12, 12> C;
    C << 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1,
        11, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1,
        10, 10, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1,
        9, 9, 9, 9, 8, 7, 6, 5, 4, 3, 2, 1,
        8, 8, 8, 8, 8, 7, 6, 5, 4, 3, 2, 1,
        7, 7, 7, 7, 7, 7, 6, 5, 4, 3, 2, 1,
        6, 6, 6, 6, 6, 6, 6, 5, 4, 3, 2, 1,
        5, 5, 5, 5, 5, 5, 5, 5, 4, 3, 2, 1,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 2, 1,
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 1,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;

    Matrix<double, 20, 20> D;
    for (int i = 0; i < 20; i++)
    {
        for (int j = 0; j < 20; j++)
        {
            D(i, j) = sqrt(2.0f / 21.0f) * sin((double)(i + 1) * (j + 1) * acos(-1.0) / 21);
        }
    }

    Matrix<double, 50, 50> E;
    for (int i = 0; i < 50; i++)
    {
        for (int j = 0; j < 50; j++)
        {
            if (j < i)
                E(i, j) = -1;

            else if (i == j || j == 49)
                E(i, j) = 1;
            else
                E(i, j) = 0;
        }
    }
    Matrix<double, 10, 10> F;
    F << -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0,

    solve(A, "A");
    solve(B, "B");
    solve(C, "C");
    solve(D, "D");
    solve(E, "E");

    std::ofstream f("../ans.txt", std::ios::app);
    std::vector<std::complex<double>> env1;
    shifted_qr(F,env1);
    sort(env1.begin(), env1.end(), [](const std::complex<double> &a, const std::complex<double> &b) { return a.real() < b.real(); });
    for (size_t i = 0; i < env1.size(); ++i)
    {
        f << env1[i] << " ";
    }
    return 0;
}
