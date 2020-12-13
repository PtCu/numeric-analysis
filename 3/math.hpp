#include <eigen3/Eigen/Eigen>
#include <map>
#include <vector>
#include <fstream>

using namespace Eigen;

const int iter = 1000;
const double tol = 1e-10;
void power_eng(const MatrixXd &a, double &pld, std::vector<double> &env)
{
    MatrixXd A = MatrixXd(a);
    VectorXd u = VectorXd::Random(a.rows());
    VectorXd v = VectorXd(u);
    for (size_t j = 1; j <= iter; ++j)
    {
        v = u / u.lpNorm<Infinity>(); //利用无穷范数进行归一化
        u = A * v;
    }
    for (size_t i = 0; i < v.size(); ++i)
    {
        env.push_back(v(i));
    }
    //env = v;
    pld = u.lpNorm<Infinity>();
    /*
      MatrixXd A = MatrixXd(a);
    VectorXd u = VectorXd::Random(a.rows());
    for (size_t j = 1; j <= iter; ++j)
    {
        u = A * u;
        u = u / u.squaredNorm();
       
    }
    for (size_t i = 0; i < u.size(); ++i)
    {
        env.push_back(u(i));
    }
    //env = v;
    auto x = u.transpose()*A*u *(u.transpose()*u).inverse(); //瑞利商
    pld = x(0);
    */
}

void inv_power_eng(const MatrixXd &a, std::vector<double> &env, double &pld)
{
    MatrixXd A = MatrixXd(a);
    VectorXd u = VectorXd::Random(a.rows());
    VectorXd v = VectorXd(u);
    for (size_t j = 1; j <= iter; ++j)
    {
        v = u / u.lpNorm<Infinity>(); 
        u = A.inverse() * v;
    }
    for (size_t i = 0; i < v.size(); ++i)
    {
        env.push_back(v(i));
    }
    pld = u.lpNorm<Infinity>();
}

//a：原矩阵
//ev: 存放特征值
bool jacobi_eng(const MatrixXd &a, std::vector<double> &ev)
{
    MatrixXd A = MatrixXd(a);
    //pdbVec = MatrixXd::Identity(a.rows(),a.rows()); //存放特征向量
    ev.resize(a.cols());
    for (size_t num = 0; num < iter; ++num)
    {
        //在非对角线上找到绝对值最大的元素
        MatrixXd::Index maxRow = 0, maxCol = 0;
        double _max = 0;
        for (size_t i = 0; i < A.rows(); ++i)
        {
            for (size_t j = 0; j < A.rows(); ++j)
            {
                if (i != j && abs(A(i, j)) > _max)
                {
                    _max = abs(A(i, j));
                    maxRow = i;
                    maxCol = j;
                }
            }
        }
        if (_max < tol)
            break;
        // double max = A.maxCoeff(&maxRow, &maxCol);
        double App = A(maxRow, maxRow);
        double Apq = A(maxRow, maxCol);
        double Aqq = A(maxCol, maxCol);

        //计算旋转角度
        double Angle = 0.5 * atan2(-2 * Apq, Aqq - App);
        double SinTheta = sin(Angle);
        double CosTheta = cos(Angle);
        double Sin2Theta = sin(2 * Angle);
        double Cos2Theta = cos(2 * Angle);

        A(maxRow, maxRow) = App * CosTheta * CosTheta + Aqq * SinTheta * SinTheta + 2 * Apq * CosTheta * SinTheta;
        A(maxCol, maxCol) = App * SinTheta * SinTheta + Aqq * CosTheta * CosTheta - 2 * Apq * CosTheta * SinTheta;
        A(maxRow, maxCol) = 0.5 * (Aqq - App) * Sin2Theta + Apq * Cos2Theta;
        A(maxCol, maxRow) = A(maxRow, maxCol);

        //eigen下标从0开始
        for (size_t i = 0; i < A.cols(); ++i)
        {
            if ((i != maxRow) && (i != maxCol))
            {
                double t1 = A(i, maxRow);
                double t2 = A(i, maxCol);
                A(i, maxRow) = CosTheta * t1 + SinTheta * t2;
                A(i, maxCol) = -SinTheta * t1 + CosTheta * t2;
            }
        }

        for (size_t i = 0; i < A.cols(); ++i)
        {
            if ((i != maxRow) && (i != maxCol))
            {
                double t1 = A(maxRow, i);
                double t2 = A(maxCol, i);
                A(maxRow, i) = CosTheta * t1 + SinTheta * t2;
                A(maxCol, i) = -SinTheta * t1 + CosTheta * t2;
            }
        }
        // //计算特征向量
        // for (size_t i = 0; i < A.cols(); ++i)
        // {
        //     pdbVec(i, maxRow) = pdbVec(i, maxCol) * SinTheta + pdbVec(i, maxRow) * CosTheta;
        //     pdbVec(i, maxCol) = -pdbVec(i, maxCol) * SinTheta + pdbVec(i, maxRow) * CosTheta;
        // }
    }

    for (size_t i = 0; i < A.cols(); ++i)
    {
        ev[i] = A(i, i);
    }
}

//无移动qr
void unshifted_qr(const MatrixXd &a, std::vector<double> &ev)
{
    ev.resize(a.cols());
    MatrixXd Q = MatrixXd::Identity(a.cols(), a.cols());
    //Qbar为特征向量
    //MatrixXd Qbar = MatrixXd(Q);
    MatrixXd R = a;
    HouseholderQR<MatrixXd> qr;
    for (size_t i = 0; i < iter; ++i)
    {
        //QR分解
        qr.compute(R * Q);
        R = qr.matrixQR().triangularView<Eigen::Upper>();
        Q = qr.householderQ();
        //迭代
        // Qbar = Qbar * Q;
    }
    //特征值位于对角线上
    MatrixXd lam = R * Q;

    for (size_t i = 0; i < a.cols(); ++i)
    {
        ev[i] = lam(i, i);
    }
}

//高斯-海森伯格
Matrix<double, Dynamic, Dynamic> gauss_hessen(const Matrix<double, Dynamic, Dynamic> &a)
{
    std::ofstream f("../hs_record.txt", std::ios::app);
    int cols = a.cols();
    int rows = a.rows();
    Matrix<double, Dynamic, Dynamic> A;
    A.resize(rows, cols);
    A = a;
    f << std::endl
      << "初始A" << std::endl
      << A << std::endl;
    //从第二列开始，到倒数第二列结束
    for (size_t i = 1; i < A.cols() - 1; ++i)
    {
        MatrixXd G = MatrixXd::Identity(A.cols(), A.cols());
        MatrixXd G_ = MatrixXd::Identity(A.cols(), A.cols());
        // bool pass = true;

        // for (size_t k = i; k < A.cols(); ++k)
        // {
        //     if (abs(A(k, i - 1)) > tol)
        //     {
        //         //交换row(k)和row(i)
        //         if (k != i)
        //         {
        //             A.row(k).swap(A.row(i));
        //             A.col(k).swap(A.col(i));
        //         }

        //         pass = false;
        //         break;
        //     }
        // }
        //下面全为0
        // if (pass)
        //     continue;
        //保证A(i,i-1)不为0
        //构造矩阵G
        double some = 1;
        for (size_t k = i; k < A.rows(); ++k)
        {
            if (A(k, i - 1) != 0)
            {
                some = A(k, i - 1);
                break;
            }
        }
        for (size_t j = i + 1; j < A.rows(); ++j)
        {
            //计算a, b, c
            G(j, i) = -A(j, i - 1) / some;
            G_(j, i) = -G(j, i);
        }
        A = G * A.*G_;
    }
    return A;
}
//平移QR法
//计算方阵的实数和复数特征值
void shifted_qr(const MatrixXd &a, std::vector<std::complex<double>> &ev)
{
    ev.resize(a.cols());

    //最大迭代次数
    int kounttol = 1000;
    //行数
    int m = a.rows();
    int n = m;
    MatrixXcd A = gauss_hessen(a);
    HouseholderQR<MatrixXcd> qr;
    while (n > 1)
    {
        int kount = 0;
        double _max = 0;
        //第n行最大的元素
        for (size_t i = 0; i <= n - 2; ++i)
        {
            if (abs(A(n - 1, i)) > _max)
            {
                _max = abs(A(n - 1, i));
            }
        }
        while ((kount < kounttol) && _max > tol)
        {
            //记录QR的个数
            kount++;
            auto mu = A(n - 1, n - 1);
            //QR分解
            qr.compute(A - mu * MatrixXcd::Identity(n, n));
            MatrixXcd R = qr.matrixQR().triangularView<Eigen::Upper>();
            MatrixXcd Q = qr.householderQ();
            A = R * Q + mu * MatrixXcd::Identity(n, n);
            double _max = 0;
            //第n行最大的元素
            for (size_t i = 0; i <= n - 2; ++i)
            {
                if (abs(A(n - 1, i)) > _max)
                {
                    _max = abs(A(n - 1, i));
                }
            }
        }
        //具有分离的1*1的块
        if (kount < kounttol)
        {
            ev[n - 1] = (A(n - 1, n - 1));
            --n;
            A = A.block(0, 0, n, n).eval();
        }
        //存在复数根
        else
        {
            std::complex<double> disc = (A(n - 2, n - 2) - A(n - 1, n - 1)) * (A(n - 2, n - 2) - A(n - 1, n - 1)) + 4.0 * (A(n - 1, n - 2) * A(n - 2, n - 1));
            ev[n - 1] = 0.5 * (A(n - 2, n - 2) + A(n - 1, n - 1) + sqrt(disc));
            ev[n - 2] = 0.5 * (A(n - 2, n - 2) + A(n - 1, n - 1) - sqrt(disc));
            n -= 2;
            A = A.block(0, 0, n, n).eval();
        }
    }
    if (n > 0)
        ev[0] = A(0, 0);
}
